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
    Int memory_gb = deep_srsnv_params.tensorize_workers * 4 + 8
    Int cpus = deep_srsnv_params.tensorize_workers + 2
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
      --num-workers ~{deep_srsnv_params.tensorize_workers} \
      --tensorize-output-rows ~{deep_srsnv_params.tensorize_output_rows} \
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

# Shard-parallel tensorization: tensorize the reads in a single genomic interval (one snvfind shard),
# reading from the feature map VCF (VCF-direct) + CRAM with the Rust tensorizer. Designed for small
# machines (~1-2 CPU / 4-6 GiB) running in parallel; outputs a per-interval tensor cache directory tar.
task DNNCramToTensorsSharded {
  input {
    File input_cram
    File input_cram_index
    File featuremap_vcf
    File featuremap_vcf_index
    File interval_bed
    String label  # "positive" | "negative" | "inference"
    References references
    DeepSRSNVParams deep_srsnv_params
    # For training/positive: restrict to reads selected by the (filter+downsample) parquet.
    File? selection_parquet
    # For inference: reuse featuremap_to_dataframe's filter via this read-filters JSON.
    File? inference_filters_json
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = select_first([deep_srsnv_params.tensorize_task_memory_gb, 4])
    Int cpus = select_first([deep_srsnv_params.tensorize_task_cpus, 2])
  }

  Float input_size = size(input_cram, "GiB") + size(featuremap_vcf, "GiB")
  Int disk_size = ceil(input_size + 20)
  # This task fans out to hundreds of tiny preemptible VMs at once; at the global default of 1 preemptible
  # retry, mass GCP preemption terminal-fails shards and aborts the whole scatter. Floor at 3 retries (a
  # non-preemptible final attempt kicks in once preemptible attempts are exhausted). Honor a higher caller value.
  Int preemptible_effective = if preemptible_tries > 3 then preemptible_tries else 3

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    cram_to_tensors \
      --cram ~{input_cram} \
      --featuremap-vcf ~{featuremap_vcf} \
      --label ~{label} \
      --output tensor_cache \
      --reference ~{references.ref_fasta} \
      --interval-bed ~{interval_bed} \
      --tensorizer ~{select_first([deep_srsnv_params.tensorizer, "rust"])} \
      --tensor-length ~{deep_srsnv_params.tensor_length} \
      --tensorize-output-rows ~{deep_srsnv_params.tensorize_output_rows} \
      --num-workers ~{cpus} \
      --channel-config ~{deep_srsnv_params.channel_registry} \
      --vocab-config ~{deep_srsnv_params.vocab_config} \
      ~{if defined(selection_parquet) then "--selection-parquet " + select_first([selection_parquet]) else ""} \
      ~{if defined(inference_filters_json) then "--inference-filters-json " + select_first([inference_filters_json]) else ""}

    tar -chf tensor_cache_shard.tar tensor_cache/
  >>>

  runtime {
    preemptible: preemptible_effective
    # maxRetries covers NON-preemption transient failures (e.g. GCS localization/read IOExceptions) that
    # `preemptible` does not; a single flaky shard should not fail the whole scatter.
    maxRetries: 2
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File tensor_cache_shard_tar = "tensor_cache_shard.tar"
    File monitoring_log = "monitoring.log"
  }
}

# Gather per-interval tensor-cache shard tars into one deterministically-ordered tensor cache tar,
# so DNNCombineSplits (unchanged) consumes a cache bit-identical to a single non-sharded run.
task DNNConcatTensorShards {
  input {
    Array[File] tensor_cache_shard_tars
    String label
    DeepSRSNVParams deep_srsnv_params
    String docker
    Int preemptible_tries
    File monitoring_script
    # All shard tars are first extracted to disk (disk holds every shard cache); the concat step itself
    # then streams — appending in interval order with only ~one input + one output shard in MEMORY at a
    # time — so MEMORY is small and roughly constant while DISK must fit all shards (see disk_size below).
    Int memory_gb = 8
    Int cpus = 2
  }

  Float input_size = size(tensor_cache_shard_tars, "GiB")
  Int disk_size = ceil(input_size * 4 + 20)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    echo "Concatenating ~{label} tensor shards"

    # Extract each per-interval shard cache into its own directory.
    mkdir -p shards
    i=0
    for tar_file in ~{sep=" " tensor_cache_shard_tars}; do
      dest="shards/shard_dir_${i}"
      mkdir -p "$dest"
      tar -xf "$tar_file" -C "$dest"
      i=$((i+1))
    done

    # Concatenate with deterministic (CHROM,POS,RN) ordering.
    concat_tensor_shards \
      $(for d in shards/shard_dir_*/tensor_cache; do echo "--tensor-cache-dir $d"; done) \
      --output tensor_cache \
      --tensorize-output-rows ~{deep_srsnv_params.tensorize_output_rows}

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

task SelectFoldBeds {
  # Return the snvfind shard BEDs whose chromosome belongs to a given CV fold (fold 0 also owns holdout),
  # using the same split manifest as cram_to_tensors --fold-idx. Enables shard-parallel inference: the
  # caller scatters DNNCramToTensorsInference over this fold's BEDs. Tiny/cheap.
  input {
    Array[File] shard_beds
    File training_interval_list
    Int fold_idx
    Int num_folds
    Int random_seed
    String holdout_chromosomes
    String docker
    Int preemptible_tries
    File monitoring_script
  }

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # --copy-dir materializes the selected BEDs as fresh files under fold_beds/ (the task then globs those
    # real, task-produced files). We must NOT output the selected input BEDs' paths directly: on the
    # Cromwell/GCP backend that makes Cromwell try to delocalize an input's localized path, building a bogus
    # doubled gs:// URL ("matched no objects") so the downstream tensorize task can't localize its fold_bed
    # and exits 1. (Omics/miniwdl tolerate the path-list form, which is why this only failed on Cromwell.)
    mkdir -p fold_beds
    select_fold_beds \
      ~{sep=" " prefix("--bed ", shard_beds)} \
      --training-regions ~{training_interval_list} \
      --fold-idx ~{fold_idx} \
      --num-folds ~{num_folds} \
      --random-seed ~{random_seed} \
      --holdout-chromosomes ~{holdout_chromosomes} \
      --copy-dir fold_beds \
      --output fold_beds.json
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
  }

  output {
    # Glob the COPIED beds (real task-produced files) so Cromwell/GCP delocalizes them correctly.
    # (Do not read_json the selected input paths — see the --copy-dir note in the command.)
    Array[File] fold_beds = glob("fold_beds/*.bed")
    File monitoring_log = "monitoring.log"
  }
}

# Launder an optional File that originates from a conditional-scoped sub-workflow output into a
# REQUIRED task output. Cromwell cannot resolve a reference to an OPTIONAL call output (File?) from a
# deeper/sibling scope (e.g. inside a scatter) — it fails with "required input lookup failed". A REQUIRED
# output (File) is auto-exposed as a clean File? across scope boundaries (exactly why FeatureMapPrep's
# required `featuremap` output works but the optional `inference_filters` did not). So this task ALWAYS
# writes the output file (empty when the input is absent) and declares it `File out` (required). The
# caller passes `out` (always a valid File); DNNCramToTensorsInference skips --inference-filters-json when
# the file is EMPTY, preserving the score-everything semantics that a missing filter previously signalled.
task PassThroughOptionalFile {
  input {
    File? in_file
    String docker
    Int preemptible_tries
    File monitoring_script
  }
  Boolean present = defined(in_file)
  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    if ~{if present then "true" else "false"}; then
      cp ~{select_first([in_file, "/dev/null"])} passed_file.out
    else
      : > passed_file.out   # always create; EMPTY file means "no inference filters" (score everything)
    fi
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
  }
  output {
    # REQUIRED output (always exists). Empty file == no filters. Required so Cromwell resolves the
    # reference from inside the nested tensorize scatter.
    File out = "passed_file.out"
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
    # VCF-direct inference: read the feature map VCF directly (no DNNVcfToParquet), applying the same
    # inference read-filter by reusing featuremap_to_dataframe's filter code, then restrict to the fold.
    File featuremap_vcf
    File featuremap_vcf_index
    File? inference_filters_json
    # Optional shard restriction: when set, tensorize only this genomic interval (shard-parallel inference).
    # The fold filter (--fold-idx) still applies on top, so a bed outside the fold's chromosomes yields nothing.
    File? interval_bed
    File training_interval_list
    References references
    DeepSRSNVParams deep_srsnv_params
    Int fold_idx
    Int? num_folds
    String base_file_name
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = select_first([deep_srsnv_params.tensorize_task_memory_gb, 4])
    Int cpus = select_first([deep_srsnv_params.tensorize_task_cpus, 2])
  }

  Float input_size = size(input_cram, "GiB") + size(featuremap_vcf, "GiB")
  Int disk_size = ceil(input_size * 2 + 50)
  String holdout = select_first([deep_srsnv_params.holdout_chromosomes, "chr21"])
  Int effective_num_folds = select_first([num_folds, deep_srsnv_params.num_folds])
  # Fans out to hundreds of tiny preemptible VMs; floor preemptible retries at 3 so mass GCP preemption
  # doesn't terminal-fail shards and abort the scatter (a non-preemptible final attempt kicks in after).
  Int preemptible_effective = if preemptible_tries > 3 then preemptible_tries else 3

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Pass --inference-filters-json only when a NON-EMPTY filter file is provided. The upstream
    # PassThroughOptionalFile always emits a file (empty == "no filters"), so gate on file size (-s)
    # rather than mere presence; an empty file means score-everything (no restriction).
    FILTER_ARG=""
    ~{if defined(inference_filters_json) then "if [ -s " + select_first([inference_filters_json]) + " ]; then FILTER_ARG=\"--inference-filters-json " + select_first([inference_filters_json]) + "\"; fi" else "true"}

    cram_to_tensors \
      --cram ~{input_cram} \
      --featuremap-vcf ~{featuremap_vcf} \
      --label inference \
      --fold-idx ~{fold_idx} \
      --training-regions ~{training_interval_list} \
      --num-folds ~{effective_num_folds} \
      --random-seed ~{deep_srsnv_params.random_seed} \
      --holdout-chromosomes ~{holdout} \
      --reference ~{references.ref_fasta} \
      --output tensor_cache \
      --tensorizer ~{select_first([deep_srsnv_params.tensorizer, "rust"])} \
      --tensor-length ~{deep_srsnv_params.tensor_length} \
      --tensorize-output-rows ~{deep_srsnv_params.tensorize_output_rows} \
      --num-workers ~{cpus} \
      --compress \
      --channel-config ~{deep_srsnv_params.channel_registry} \
      --vocab-config ~{deep_srsnv_params.vocab_config} \
      ${FILTER_ARG} \
      ~{if defined(interval_bed) then "--interval-bed " + select_first([interval_bed]) else ""}

    ls -ltr tensor_cache/
  >>>

  runtime {
    preemptible: preemptible_effective
    # maxRetries covers NON-preemption transient failures (e.g. GCS localization/read IOExceptions) that
    # `preemptible` does not; a single flaky shard should not fail the whole scatter.
    maxRetries: 2
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    Array[File] tensor_shards = glob("tensor_cache/shard_*.pt.gz")
    File tensor_cache_index = "tensor_cache/index.json"
    File monitoring_log = "monitoring.log"
  }
}

task DNNFoldInference {
  input {
    Array[File] tensor_shards
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

  Float input_size = size(tensor_shards, "GiB")
  Int disk_size = ceil(input_size + 20)
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

    # Link tensor shards into a single directory for inference. Each source shard-tensorize task names
    # its outputs shard_00000.pt.gz, shard_00001.pt.gz, ... so files from DIFFERENT shard tasks share
    # basenames; linking by basename (ln -sf "$shard" tensor_cache/) would collide and silently keep only
    # the last one per name (dropping ~all tensors when many shards are flattened in). Renumber each linked
    # file to a globally-unique name. Inference scores every tensor independently and the merge joins by
    # (CHROM,POS,RN), so the on-disk ordering/names do not matter — only that ALL shards are present.
    mkdir -p tensor_cache
    idx=0
    while read -r shard; do
      ln -sf "$shard" "$(printf 'tensor_cache/shard_%08d.pt.gz' "$idx")"
      idx=$((idx + 1))
    done < ~{write_lines(tensor_shards)}
    echo "Linked $(ls tensor_cache/ | wc -l) shards into tensor_cache/ (from $idx input files)"

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
