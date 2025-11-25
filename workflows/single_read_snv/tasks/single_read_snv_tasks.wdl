version 1.0
# LICENSE
#   Copyright 2022 Ultima Genomics
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#


import "structs.wdl" as Structs
task PrepareFeatureMapForTraining {
  parameter_meta {
    memory_gb: {
      help: "Memory in GB to allocate for the PrepareFeatureMapForTraining task",
      type: "Int",
      category: "runtime"
    }
    cpus: {
      help: "Number of CPUs to allocate for the PrepareFeatureMapForTraining task",
      type: "Int",
      category: "runtime"
    }
  }
  input {
    File featuremap
    File featuremap_index
    File training_regions_interval_list
    File training_regions_interval_list_index
    FeatureMapParams featuremap_params
    SingleReadSNVParams single_read_snv_params
    Int train_set_size
    Float mean_coverage
    Array[String]? filters  # Additional filters to apply after pre_filters
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 2
    Int cpus = 1
  }

  Float featuremap_size = size(featuremap, "GB")
  Int disk_size = ceil(featuremap_size * 3 + 20)

  String base_file_name = basename(featuremap, ".vcf.gz")
  String intermediate_featuremap = base_file_name + ".training_regions.vcf.gz"
  String intermediate_parquet = base_file_name + ".training_regions.parquet"
  String filtered_parquet = base_file_name + ".filtered.parquet"
  String stats_out_json = base_file_name + ".stats.json"
  
  # Calculate coverage threshold (factor times mean coverage)
  Int coverage_threshold = ceil(mean_coverage * single_read_snv_params.max_coverage_factor)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Filter featuremap by training regions
    bcftools view ~{featuremap} -T ~{training_regions_interval_list} -Oz -o ~{intermediate_featuremap} 
    bcftools index -t ~{intermediate_featuremap}

    # Convert featuremap to parquet
    featuremap_to_dataframe \
    --input ~{intermediate_featuremap} \
    --output ~{intermediate_parquet} \
    --drop-format GT AD X_TCM

    # Run filtering with hardcoded and user-defined filters
    filter_featuremap \
    --in ~{intermediate_parquet} \
    --out ~{filtered_parquet} \
    --stats ~{stats_out_json} \
    --filter name=coverage_ge_min:field=DP:op=ge:value=~{single_read_snv_params.min_coverage_filter}:type=region \
    --filter name=coverage_le_max:field=DP:op=le:value=~{coverage_threshold}:type=region \
    ~{true="--filter " false="" defined(single_read_snv_params.pre_filters)}~{sep=" --filter " single_read_snv_params.pre_filters} \
    ~{true="--filter " false="" defined(filters)}~{sep=" --filter " filters} \
    --downsample random:~{train_set_size}:~{single_read_snv_params.random_seed}
    
    # Print stats JSON content
    echo "=== Filter Statistics ==="
    cat ~{stats_out_json}
    
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File filtered_featuremap_parquet = "~{filtered_parquet}"
    File stats_file = "~{stats_out_json}"
    File monitoring_log = "monitoring.log"
  }
}

task TrainModel {
  parameter_meta {
    memory_gb: {
      help: "Memory in GB to allocate for the TrainModel task",
      type: "Int",
      category: "runtime"
    }
    cpus: {
      help: "Number of CPUs to allocate for the TrainModel task",
      type: "Int",
      category: "runtime"
    }
  }
  input {
    File raw_filtered_featuremap_parquet
    File random_sample_filtered_featuremap_parquet
    File random_sample_positive_label_stats
    File random_sample_negative_label_stats
    File raw_featuremap_stats
    Float mean_coverage
    File training_regions_interval_list
    File xgboost_params_file
    SingleReadSNVParams single_read_snv_params
    Array[String] features
    String base_file_name
    String docker
    String pipeline_version
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 32
    Int cpus = ceil(memory_gb / 2)
  }

  Float featuremap_size = size(raw_filtered_featuremap_parquet, "GB") + size(random_sample_filtered_featuremap_parquet, "GB")
  Int disk_size = ceil(featuremap_size * 3 + size(training_regions_interval_list, "GB") + 20)

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    srsnv_training \
    --positive ~{random_sample_filtered_featuremap_parquet} \
    --negative ~{raw_filtered_featuremap_parquet} \
    --stats-positive ~{random_sample_positive_label_stats} \
    --stats-negative ~{random_sample_negative_label_stats} \
    --stats-featuremap ~{raw_featuremap_stats} \
    --mean-coverage ~{mean_coverage} \
    --training-regions ~{training_regions_interval_list} \
    --k-folds ~{single_read_snv_params.num_CV_folds} \
    --model-params ~{xgboost_params_file} \
    --output $PWD \
    --basename ~{base_file_name} \
    --features ~{sep=":" features} \
    --random-seed ~{single_read_snv_params.random_seed} \
    --metadata docker_image="~{docker}" \
    --metadata pipeline_version="~{pipeline_version}" \
    --verbose

    ls -ltr
    
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File featuremap_df = "~{base_file_name}.featuremap_df.parquet"
    File srsnv_metadata_json = "~{base_file_name}.srsnv_metadata.json"
    Array[File] model_files = glob("~{base_file_name}.model_fold_*.json")
    File monitoring_log = "monitoring.log"
  }
}

task CreateFeatureMap {
  parameter_meta {
    memory_gb: {
      help: "Memory in GB to allocate for the CreateFeatureMap task",
      type: "Int",
      category: "runtime"
    }
    cpus: {
      help: "Number of CPUs to allocate for the CreateFeatureMap task",
      type: "Int",
      category: "runtime"
    }
  }
  input {
    Array[File] input_cram_bam_list
    Array[File] input_cram_bam_index_list
    References references
    String base_file_name
    String total_aligned_bases
    Int random_sample_size
    File? random_sample_trinuc_freq_
    FeatureMapParams featuremap_params
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 2
    Int cpus = 1
  }

  Float input_size = size(input_cram_bam_list, "GB")
  Float ref_size = size(references.ref_fasta, "GB")
  Int disk_size = ceil(input_size * 3 + ref_size + 20)

  String out_vcf = "~{base_file_name}.raw.featuremap.vcf.gz"
  String out_vcf_random_sample = "~{base_file_name}.random_sample.featuremap.vcf.gz"
  String out_random_sample_trinuc_freq = "~{base_file_name}.random_sample.trinuc_freq.csv"

  command <<<
    lscpu
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Build downsampling rate
    TOTAL_ALIGNED_BASES="~{total_aligned_bases}"
    echo "Total aligned bases: $TOTAL_ALIGNED_BASES"
    DOWNSAMPLING_RATE=$(awk -v num=~{random_sample_size} -v den="$TOTAL_ALIGNED_BASES" 'BEGIN{printf "%.12f", num/den}')
    echo "Downsampling rate: $DOWNSAMPLING_RATE"

    # Build random sample args
    RANDOM_SAMPLE_TRINUC_ARGS=""
    if [ -n "~{default='' random_sample_trinuc_freq_}" ]; then
      RANDOM_SAMPLE_TRINUC_ARGS=",~{random_sample_trinuc_freq_},~{out_random_sample_trinuc_freq}"
    fi

    # Run snvfind with parameters
    snvfind \
      ~{sep="," input_cram_bam_list} \
      ~{references.ref_fasta} \
      -o ~{out_vcf} \
      -f ~{out_vcf_random_sample},${DOWNSAMPLING_RATE}${RANDOM_SAMPLE_TRINUC_ARGS} \
      ~{true="-p" false="" defined(featuremap_params.padding_size)}~{default="" featuremap_params.padding_size} \
      ~{true="-L" false="" defined(featuremap_params.score_limit)}~{default="" featuremap_params.score_limit} \
      ~{true="-X" false="" defined(featuremap_params.max_score_to_emit)}~{default="" featuremap_params.max_score_to_emit} \
      ~{true="-N" false="" defined(featuremap_params.min_score_to_emit)}~{default="" featuremap_params.min_score_to_emit} \
      ~{true="-n" false="" select_first([featuremap_params.exclude_nan_scores, true])} \
      ~{true="-d" false="" select_first([featuremap_params.include_dup_reads, false])} \
      ~{true="-k" false="" select_first([featuremap_params.keep_supplementary, false])} \
      ~{true="-Q" false="" defined(featuremap_params.surrounding_quality_size)}~{default="" featuremap_params.surrounding_quality_size} \
      ~{true="-r" false="" defined(featuremap_params.reference_context_size)}~{default="" featuremap_params.reference_context_size} \
      ~{true="-m" false="" defined(featuremap_params.min_mapq)}~{default="" featuremap_params.min_mapq} \
      ~{true="-c" false="" defined(featuremap_params.cram_tags_to_copy)} ~{default="" sep="," featuremap_params.cram_tags_to_copy} \
      ~{true="-C" false="" defined(featuremap_params.attributes_prefix)} ~{default="" featuremap_params.attributes_prefix} \
      ~{true="-b" false="" defined(featuremap_params.bed_file)} ~{default="" featuremap_params.bed_file}
   
    bcftools index -t ~{out_vcf_random_sample}
    bcftools index -t ~{out_vcf}

    printf '%s\n' "${DOWNSAMPLING_RATE}" > downsampling_rate.txt

    ls -ltr
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File featuremap = "~{out_vcf}"
    File featuremap_index = "~{out_vcf}.tbi"
    File featuremap_random_sample = "~{out_vcf_random_sample}"
    File featuremap_random_sample_index = "~{out_vcf_random_sample}.tbi"
    File? random_sample_trinuc_freq = "~{out_random_sample_trinuc_freq}"
    Float downsampling_rate = read_float("downsampling_rate.txt")
    File monitoring_log = "monitoring.log"
  }
}

task Inference {
  input {
    String base_file_name
    Array[File] model_files
    File srsnv_metadata_json
    File featuremap
    File monitoring_script
    String docker
    Int preemptible_tries
    Int cpus = 4
    Int memory_gb = 8
  }

  Float featuremap_size = size(featuremap, "GB")
  Int disk_size = ceil(featuremap_size * 2 + 20)

  String out_vcf = "~{base_file_name}.featuremap.vcf.gz"

  command <<<
      lscpu
      set -xeuo pipefail
      bash ~{monitoring_script} | tee monitoring.log >&2 &
 
      # prepare model files directory
      mkdir -p model_files
      cp ~{sep=" " model_files} model_files/
      cp ~{srsnv_metadata_json} model_files/srsnv_metadata.json

      # Run inference with snvqual
      snvqual \
        "~{featuremap}" \
        "~{out_vcf}" \
        model_files/srsnv_metadata.json \
        -v

      bcftools index -t "~{out_vcf}"

      ls -ltr
   >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File featuremap_out = "~{out_vcf}"
    File featuremap_out_index = "~{out_vcf}.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task CreateReport {
  input {
    File   featuremap_df
    File   srsnv_metadata_json
    Array[File] model_files
    String basename
    String docker
    Int    preemptible_tries
    File   monitoring_script
    Int    memory_gb = 32
    Int    cpus      = 16
  }

  Float fm_size = size(featuremap_df, "GB")
  Int   disk_size = ceil(fm_size + 10)

  command <<<
    set -euox pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Copy metadata and model files to PWD (Cromwell and Omics path support)
    cp ~{sep=" " model_files} .
    cp ~{srsnv_metadata_json} srsnv_metadata.json

    srsnv_report \
      --featuremap-df ~{featuremap_df} \
      --srsnv-metadata srsnv_metadata.json \
      --report-path . \
      --basename ~{basename} \
      --verbose

    ls -ltr
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_size} HDD"
  }

  output {
    File application_qc_h5 = "~{basename}.single_read_snv.applicationQC.h5"
    File report_html       = "~{basename}.report.html"
    File monitoring_log    = "monitoring.log"
  }
}