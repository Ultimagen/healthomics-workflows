version 1.0

import "structs.wdl" as Structs

task DNNCramToTensors {
  input {
    File input_cram
    File input_cram_index
    File featuremap_parquet
    String label  # "positive" or "negative"
    References references
    DeepSRSNVParams deep_srsnv_params
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = deep_srsnv_params.num_tensorize_workers * 4 + 8
    Int cpus = deep_srsnv_params.num_tensorize_workers + 2
  }

  Float input_size = size(input_cram, "GiB") + size(featuremap_parquet, "GiB")
  Int disk_size = ceil(input_size * 3 + 50)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    cram_to_tensors \
      --cram ~{input_cram} \
      --parquet ~{featuremap_parquet} \
      --label ~{label} \
      --output tensor_cache \
      --reference ~{references.ref_fasta} \
      --tensor-length ~{deep_srsnv_params.tensor_length} \
      --num-workers ~{deep_srsnv_params.num_tensorize_workers} \
      --shard-size ~{deep_srsnv_params.shard_size} \
      --channel-config ~{deep_srsnv_params.channel_registry} \
      --vocab-config ~{deep_srsnv_params.vocab_config}

    tar -chf tensor_cache.tar tensor_cache/
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File tensor_cache_tar = "tensor_cache.tar"
    File monitoring_log = "monitoring.log"
  }
}

task DNNCombineSplits {
  input {
    File positive_tensor_cache_tar
    File negative_tensor_cache_tar
    File training_interval_list
    DeepSRSNVParams deep_srsnv_params
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 128
    Int cpus = 16
  }

  Float input_size = size(positive_tensor_cache_tar, "GiB") + size(negative_tensor_cache_tar, "GiB")
  # combine_splits writes ~fold_data output, then re-tars each fold (duplicating it on disk),
  # so the peak footprint is roughly input + fold_data + fold_data-tar. Size generously.
  Int disk_size = ceil(input_size * 10 + 100)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    mkdir -p pos_cache neg_cache
    tar -xf ~{positive_tensor_cache_tar} -C pos_cache --strip-components=1
    tar -xf ~{negative_tensor_cache_tar} -C neg_cache --strip-components=1

    combine_splits \
      --positive pos_cache \
      --negative neg_cache \
      --training-regions ~{training_interval_list} \
      --k-folds ~{deep_srsnv_params.num_folds} \
      --holdout-chromosomes ~{deep_srsnv_params.holdout_chromosomes} \
      --random-seed ~{deep_srsnv_params.random_seed} \
      --output fold_data

    # Tar each fold separately and collect paths
    for i in $(seq 0 $((~{deep_srsnv_params.num_folds} - 1))); do
      tar -chf "fold_${i}.tar" -C fold_data "fold_${i}"
    done

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    Array[File] fold_tars = glob("fold_*.tar")
    File split_manifest = "fold_data/split_manifest.json"
    File monitoring_log = "monitoring.log"
  }
}

task DNNTrainFold {
  input {
    File fold_tar
    Int fold_idx
    DeepSRSNVParams deep_srsnv_params
    File? pretrained_checkpoint
    File stats_file
    Float mean_coverage
    File training_interval_list
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    # Training is GPU-bound: peak CPU across 2655 real folds was 2.2 vCPU (~6% of the
    # old 16). 4 vCPU leaves ample headroom and lets small folds land on g5.2xlarge/xlarge.
    Int cpus = 4
    Int? override_memory_gb   # optional: bypass the fold-size-based default (outlier folds)
  }

  Float input_size = size(fold_tar, "GiB")
  # Default calibrated to real peak memory across 2655 training folds (dev/prod/customer,
  # 4 mo): peak_mem ~= 3.25 + 1.53*fold_tar (r=0.986), so the request tracks the peak trend
  # with a ~3 GiB safety margin instead of the old 2x over-provision. Keeps every observed
  # fold <= 64 GiB (max request 57), so folds land on g5.4xlarge instead of g5.8xlarge.
  # Overridable via override_memory_gb for folds larger than anything seen.
  Int memory_gb = select_first([override_memory_gb, ceil(input_size * 1.53 + 9)])
  Int disk_size = ceil(input_size * 4 + 50)
  Boolean deterministic = select_first([deep_srsnv_params.deterministic, false])
  # Deterministic training runs on a single device (command passes --devices 1). Pin the requested
  # GPU count to 1 too, so the runtime doesn't over-allocate GPUs (extra cost / scheduling) that the
  # single-device training would never use.
  Int gpu_count = if deterministic then 1 else select_first([deep_srsnv_params.training_gpu_count, deep_srsnv_params.gpu_count])
  String gpu_type = select_first([deep_srsnv_params.gpu_type, "nvidia-tesla-t4"])
  String fold_basename = "~{base_file_name}_fold_~{fold_idx}"

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Deterministic cuBLAS GEMMs require this to be set before CUDA initializes.
    ~{if deterministic then "export CUBLAS_WORKSPACE_CONFIG=:4096:8" else "true"}

    mkdir -p fold_dir
    tar -xf ~{fold_tar} -C fold_dir --strip-components=1

    deep_srsnv_training \
      --fold-dir fold_dir \
      ~{if defined(pretrained_checkpoint) then "--pretrained-checkpoint " + pretrained_checkpoint else ""} \
      --training-regions ~{training_interval_list} \
      --stats-file ~{stats_file} \
      --mean-coverage ~{mean_coverage} \
      --epochs ~{deep_srsnv_params.epochs} \
      --patience ~{deep_srsnv_params.patience} \
      --batch-size ~{deep_srsnv_params.batch_size} \
      --lr-scheduler ~{deep_srsnv_params.lr_scheduler} \
      --learning-rate ~{deep_srsnv_params.learning_rate} \
      --random-seed ~{deep_srsnv_params.random_seed} \
      ~{if deterministic then "--deterministic" else ""} \
      ~{if deep_srsnv_params.use_amp then "--use-amp" else ""} \
      --use-tf32 \
      ~{if deterministic then "--devices 1" else "--devices auto"} \
      --channel-config ~{deep_srsnv_params.channel_registry} \
      --vocab-config ~{deep_srsnv_params.vocab_config} \
      --basename ~{fold_basename} \
      --output .

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
    gpuType: "nvidia-tesla-t4"
    gpuCount: gpu_count
    acceleratorType: gpu_type  #!UnknownRuntimeKey
    acceleratorCount: gpu_count  #!UnknownRuntimeKey
  }

  output {
    File dnn_metadata = "~{fold_basename}.srsnv_dnn_metadata.json"
    File dnn_checkpoint = glob("~{fold_basename}*.ckpt")[0]
    File dnn_onnx = "~{fold_basename}.dnn_model.onnx"
    File dnn_engine = "~{fold_basename}.dnn_model.engine"
    # TensorRT timing cache produced during engine serialization. Carried to inference_only so a
    # rebuilt engine reuses the same tactic selection (bit-identical on the same GPU + TRT version).
    File dnn_trt_timing_cache = "~{fold_basename}.dnn_model.engine.timing.cache"
    File dnn_featuremap_df = "~{fold_basename}.featuremap_df.parquet"
    File monitoring_log = "monitoring.log"
  }
}

task DNNRecalibrateFolds {
  input {
    Array[File] fold_parquets
    Array[File] fold_metadata
    File stats_file
    File training_interval_list
    Float mean_coverage
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Array[File] featuremap_parquets = []
    Int memory_gb = 32
    Int cpus = 8
  }

  Float fold_size = size(fold_parquets, "GiB")
  Float fm_size = size(featuremap_parquets, "GiB")
  Int disk_size = ceil((fold_size + fm_size) * 3 + 20)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    recalibrate_dnn_folds \
      --fold-parquets ~{sep=" " fold_parquets} \
      --fold-metadata ~{sep=" " fold_metadata} \
      --stats-file ~{stats_file} \
      --training-regions ~{training_interval_list} \
      --mean-coverage ~{mean_coverage} \
      --output-dir . \
      --basename ~{base_file_name} \
      ~{if length(featuremap_parquets) > 0 then "--featuremap-parquets" else ""} ~{sep=" " featuremap_parquets}

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File combined_featuremap_df = "~{base_file_name}.featuremap_df.parquet"
    File shared_lut_metadata = "~{base_file_name}.shared_lut_metadata.json"
    Array[File] updated_fold_metadata = glob("~{base_file_name}_fold_*.srsnv_dnn_metadata.json")
    File monitoring_log = "monitoring.log"
  }
}

task DNNVcfToParquet {
  input {
    File featuremap_vcf
    File featuremap_vcf_index
    File? inference_filters
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 32
    Int cpus = 4
  }

  Float input_size = size(featuremap_vcf, "GiB")
  Int disk_size = ceil(input_size * 3 + 20)
  String out_parquet = "~{base_file_name}.featuremap.parquet"

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    featuremap_to_dataframe \
      --input ~{featuremap_vcf} \
      --output ~{out_parquet} \
      --drop-format AD \
      ~{if defined(inference_filters) then "--read-filters-json " + select_first([inference_filters]) + " --read-filter-json-key filters_inference" else ""} \
      --verbose

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File featuremap_parquet = "~{out_parquet}"
    File monitoring_log = "monitoring.log"
  }
}

task DNNCramToTensorsInference {
  parameter_meta {
    num_folds: {
      help: "Override number of folds for region assignment (e.g., from length of fold_metadata array)",
      type: "Int?",
      category: "optional"
    }
  }
  input {
    File input_cram
    File input_cram_index
    File featuremap_parquet
    File training_interval_list
    References references
    DeepSRSNVParams deep_srsnv_params
    Int fold_idx
    Int? num_folds
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Int cpus = deep_srsnv_params.num_tensorize_workers + 2
  }

  Float parquet_size = size(featuremap_parquet, "GiB")
  Int memory_gb = ceil(parquet_size * 4.5 + 30)
  Float input_size = size(input_cram, "GiB") + parquet_size
  Int disk_size = ceil(input_size * 3 + 50)
  String holdout = select_first([deep_srsnv_params.holdout_chromosomes, "chr21"])
  Int effective_num_folds = select_first([num_folds, deep_srsnv_params.num_folds])
  String tensor_cache_tar_name = "tensor_cache.tar"

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    cram_to_tensors \
      --cram ~{input_cram} \
      --parquet ~{featuremap_parquet} \
      --label inference \
      --fold-idx ~{fold_idx} \
      --training-regions ~{training_interval_list} \
      --num-folds ~{effective_num_folds} \
      --random-seed ~{deep_srsnv_params.random_seed} \
      --holdout-chromosomes ~{holdout} \
      --reference ~{references.ref_fasta} \
      --output tensor_cache \
      --tensor-length ~{deep_srsnv_params.tensor_length} \
      --num-workers ~{deep_srsnv_params.num_tensorize_workers} \
      --shard-size ~{deep_srsnv_params.shard_size} \
      --compress \
      --channel-config ~{deep_srsnv_params.channel_registry} \
      --vocab-config ~{deep_srsnv_params.vocab_config}

    ls -ltr tensor_cache/

    # Singularity binds each input file as a separate argument, so a fold's ~1000 shards overflow
    # the exec argument limit in DNNFoldInference (index.json excluded: not in its tensor_cache/).
    if [ -n "${SINGULARITY_CONTAINER:-}${APPTAINER_CONTAINER:-}" ]; then
      tar -chf ~{tensor_cache_tar_name} --exclude index.json tensor_cache/
      rm -f tensor_cache/shard_*.pt.gz
    fi
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    Array[File] tensor_shards = glob("tensor_cache/shard_*.pt.gz")
    Array[File] tensor_cache_tar = glob(tensor_cache_tar_name)
    File tensor_cache_index = "tensor_cache/index.json"
    File monitoring_log = "monitoring.log"
  }
}

task DNNFoldInference {
  input {
    Array[File] tensor_shards = []
    # Set instead of tensor_shards when the producing task archived them (singularity backend).
    Array[File] tensor_cache_tar = []
    File fold_metadata
    File fold_checkpoint
    File fold_onnx_model
    File fold_engine
    # Optional TensorRT timing cache from the training run. When provided and the engine is rebuilt
    # (inference_only), it is passed to dnn_build_trt_engine so the rebuilt engine reuses the training
    # tactic selection -> bit-identical on the same GPU + TRT version.
    File? fold_timing_cache
    Int fold_idx
    DeepSRSNVParams deep_srsnv_params
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    # Rebuild the TensorRT engine from ONNX in this runtime before inference. Set true only for
    # inference_only (a prebuilt engine from another GPU/TRT env fails to deserialize here). In
    # full mode the engine was built by DNNTrainFold in this same pipeline, so no rebuild is needed.
    Boolean rebuild_engine_from_onnx = false
    Int memory_gb = 16
    Int cpus = 4
  }

  Float input_size = size(tensor_shards, "GiB") + size(tensor_cache_tar, "GiB")
  Int disk_size = ceil(input_size + size(tensor_cache_tar, "GiB") + 20)
  String out_parquet = "~{base_file_name}.fold_~{fold_idx}.predictions.parquet"
  Int gpu_count = 1
  String gpu_type = select_first([deep_srsnv_params.gpu_type, "nvidia-tesla-t4"])

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # GPU monitoring
    (while true; do
      nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.total \
        --format=csv,noheader 2>/dev/null | while IFS=, read -r idx util used total; do
        echo "MONITORING, [$(date)], GPU:$idx, %GPU_Util:$util, GPU_Mem:${used}/${total} MiB"
      done
      sleep 10
    done) >> monitoring.log 2>&1 &

    # Link tensor shards into a single directory for inference.
    # Use ln -sf to handle potential duplicates from HealthOmics caching.
    mkdir -p tensor_cache
    tar_list=~{write_lines(tensor_cache_tar)}
    if [ -s "$tar_list" ]; then
      tar -xf "$(head -n 1 "$tar_list")" -C tensor_cache --strip-components=1
    else
      cat ~{write_lines(tensor_shards)} | while read -r shard; do
        ln -sf "$shard" tensor_cache/
      done
    fi
    echo "Linked $(ls tensor_cache/ | wc -l) shards into tensor_cache/"

    # Copy model files to working directory so metadata can find them by relative path
    cp ~{fold_checkpoint} .
    cp ~{fold_onnx_model} .
    cp ~{fold_engine} .
    ~{if defined(fold_timing_cache) then "cp " + fold_timing_cache + " ." else "true"}

    # inference_only: the provided .engine was serialized in a different GPU/TensorRT runtime
    # and will fail to deserialize here. Rebuild it from this fold's ONNX in the current runtime,
    # overwriting the copied engine at the bare filename the metadata's trt_engine_path points to
    # (so dnn_fold_inference_from_cache resolves it with no metadata change). Each k-fold instance
    # rebuilds its own engine from its own ONNX. Skipped in full mode (engine already built in-env).
    # A provided timing cache is passed through so the rebuilt engine reuses the training tactic
    # selection (bit-identical on the same GPU + TRT version).
    if ~{true="true" false="false" rebuild_engine_from_onnx}; then
      echo "Rebuilding TRT engine from ONNX for cross-environment inference..."
      dnn_build_trt_engine \
        --onnx "$(basename ~{fold_onnx_model})" \
        --fold-metadata ~{fold_metadata} \
        --tensor-length ~{deep_srsnv_params.tensor_length} \
        ~{if defined(fold_timing_cache) then "--timing-cache \"$(basename " + fold_timing_cache + ")\"" else ""}
    fi

    dnn_fold_inference_from_cache \
      --tensor-cache tensor_cache/ \
      --fold-metadata ~{fold_metadata} \
      --output ~{out_parquet} \
      --batch-size ~{deep_srsnv_params.batch_size}

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
    gpuType: "nvidia-tesla-t4"
    gpuCount: gpu_count
    acceleratorType: gpu_type  #!UnknownRuntimeKey
    acceleratorCount: gpu_count  #!UnknownRuntimeKey
  }

  output {
    File fold_predictions_parquet = "~{out_parquet}"
    File monitoring_log = "monitoring.log"
  }
}

task DNNMergeAndAnnotate {
  input {
    File featuremap_vcf
    File featuremap_vcf_index
    Array[File] fold_predictions
    File fold_metadata_0
    DeepSRSNVParams deep_srsnv_params
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Int cpus = 4
  }

  Float input_size = size(featuremap_vcf, "GiB")
  # +50 base (was +20): peak memory is driven by MAX_PARALLEL_CHUNKS concurrent per-chunk
  # RN-explode joins, which scale with pileup depth, not VCF size — a small but deep targeted
  Int memory_gb = ceil(input_size * 1 + 50)
  Int disk_size = ceil(input_size * 5 + 20)
  String out_vcf = "~{base_file_name}.dnn_annotated.featuremap.vcf.gz"
  Float low_qual_threshold = select_first([deep_srsnv_params.low_qual_threshold, 40.0])

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    dnn_merge_and_annotate \
      --featuremap-vcf ~{featuremap_vcf} \
      --fold-predictions ~{sep=" " fold_predictions} \
      --fold-metadata ~{fold_metadata_0} \
      --output ~{out_vcf} \
      --low-qual-threshold ~{low_qual_threshold} \
      --chunk-size ~{select_first([deep_srsnv_params.dnn_merge_chunk_size, 2500000])} \
      --max-parallel-chunks ~{select_first([deep_srsnv_params.dnn_merge_max_parallel_chunks, 4])}

    bcftools index -f -t ~{out_vcf}

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File dnn_featuremap_vcf = "~{out_vcf}"
    File dnn_featuremap_vcf_index = "~{out_vcf}.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task DNNPrepareReport {
  input {
    File positive_featuremap_df
    File negative_featuremap_df
    File dnn_combined_featuremap_df
    File training_metadata
    File dnn_fold_0_metadata
    Array[File] dnn_fold_metadata
    Array[String] features
    String base_file_name
    String pipeline_version
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 16
    Int cpus = 4
  }

  Float input_size = size(positive_featuremap_df, "GiB") + size(negative_featuremap_df, "GiB") + size(dnn_combined_featuremap_df, "GiB")
  Int disk_size = ceil(input_size * 3 + 20)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    prepare_dnn_report \
      --training-parquet ~{positive_featuremap_df} \
      --negative-parquet ~{negative_featuremap_df} \
      --dnn-parquet ~{dnn_combined_featuremap_df} \
      --training-metadata ~{training_metadata} \
      --dnn-metadata ~{dnn_fold_0_metadata} \
      --dnn-fold-metadata ~{sep=" " dnn_fold_metadata} \
      --features ~{sep=":" features} \
      --pipeline-version ~{pipeline_version} \
      --docker-image ~{docker} \
      --output-dir . \
      --basename ~{base_file_name}

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File report_featuremap_df = "~{base_file_name}.featuremap_df.parquet"
    File report_metadata = "~{base_file_name}.srsnv_metadata.json"
    File monitoring_log = "monitoring.log"
  }
}

task DNNReport {
  input {
    File report_featuremap_df
    File report_metadata
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 16
    Int cpus = 4
  }

  Float input_size = size(report_featuremap_df, "GiB")
  Int disk_size = ceil(input_size + 10)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    srsnv_report \
      --featuremap-df ~{report_featuremap_df} \
      --srsnv-metadata ~{report_metadata} \
      --report-path . \
      --basename ~{base_file_name} \

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File dnn_report_html = "~{base_file_name}.report.html"
    File dnn_application_qc_h5 = "~{base_file_name}.single_read_snv.applicationQC.h5"
    File monitoring_log = "monitoring.log"
  }
}
