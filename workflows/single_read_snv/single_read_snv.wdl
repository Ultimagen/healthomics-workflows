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
#   See the License specific language governing permissions and
#   limitations under the License.
#

# DESCRIPTION
# Runs substitution error analysis for Ultima Genomics data. Includes three subworkflows:
# 1. FeatureMap generation - aggregates all the substitutions from a bam file using snvfind
# 2. Annotation with additional features
# 3. Single read ML model training and inference

# CHANGELOG in reverse chronological order
# 1.19.0 Replaced FeatureMap generation with snvfind from featuremap_docker
# 1.10.1 Added LA_7 adapter version
# 1.7.0 Initial implementation of single_read_snv wdl

import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/qc_tasks.wdl" as UGQCTasks
import "tasks/mrd.wdl" as UGMrdTasks
import "tasks/globals.wdl" as Globals


workflow SingleReadSNV {
input {
  File input_cram_bam
  File input_cram_bam_index
  File sorter_json_stats_file
  String base_file_name
  String pipeline_version = "1.23.2"
  References references

  FeatureMapParams featuremap_params
  SingleReadSNVParams? single_read_snv_params
  Array[String] features
  Int? featuremap_scatter_count_override

  # snv quality model include/exclude regions
  File training_regions_interval_list
  File training_regions_interval_list_index
  File xgboost_params_file

  Float min_coverage_to_train_model

  Array[File]? pre_trained_model_files
  File? pre_trained_srsnv_metadata_json

  Boolean raise_exceptions_in_report

  Int? process_featuremap_memory_gb_override
  Int? preemptible_tries
  Boolean? no_address_override
  String? cloud_provider_override
  # Used for running on other clouds (aws)
  File? monitoring_script_input

  Boolean create_md5_checksum_outputs = false
  # winval validations
  #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
  #@wv not("test" in base_file_name or "train" in base_file_name)
  #@wv suffix(input_cram_bam) in {".bam", ".cram"}
  #@wv suffix(input_cram_bam_index) in {".bai", ".crai"}
  #@wv prefix(basename(input_cram_bam_index)) == basename(input_cram_bam)
  #@wv prefix(basename(training_regions_interval_list_index)) == basename(training_regions_interval_list)
  #@wv suffix(training_regions_interval_list) in {".gz"}
  #@wv suffix(prefix(training_regions_interval_list)) in {".interval_list"}
  #@wv suffix(sorter_json_stats_file) == ".json"
  #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
  #@wv suffix(references['ref_dict']) == '.dict'
  #@wv suffix(references['ref_fasta_index']) == '.fai'
  #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
  #@wv not defined(pre_trained_model_files) -> defined(single_read_snv_params)
  #@wv defined(single_read_snv_params) -> single_read_snv_params['num_CV_folds'] > 1
  #@wv defined(pre_trained_model_files) <-> defined(pre_trained_srsnv_metadata_json)
  ##@wv defined(pre_trained_model_files) -> suffix(pre_trained_model_files) == ".json"
  #@wv defined(pre_trained_srsnv_metadata_json) -> suffix(pre_trained_srsnv_metadata_json) == ".json"
  
}

meta {
    description : "Single Read SNV Quality Recalibration workflow (single_read_snv wdl) assigns accurate quality scores to all SNV candidates. The output is a FeatureMap VCF file with the recalibrated SNV quality scores. The input cram file coverage must be over some minimal coverage for a new model to be trained and quality scores to be generated, otherwise a pre-trained model can be provided, or a FeatureMap with no scores is generated."
    author: "Ultima Genomics"
    WDL_AID: { exclude: [
        "pipeline_version",
        "cloud_provider_override",
        "no_address_override",
        "preemptible_tries",
        "monitoring_script_input",
        "process_featuremap_memory_gb_override",
        "featuremap_scatter_count_override",
        "CreateTrainingRegionsBed.disk_size",
        "CreateTrainingRegionsBed.cpus",
        "TrainModel.cpus",
        "TrainModel.memory_gb",
        "TrainModel.disk_size",
        "Inference.cpus",
        "Inference.memory_gb",
        "Inference.disk_size",
        "Inference.input_size",
        "Inference.out_vcf",
        "CreateFeatureMap.cpus",
        "CreateFeatureMap.memory_gb",
        "CreateReport.cpus",
        "CreateReport.memory_gb",
        "PrepareRawFeatureMap.cpus",
        "PrepareRawFeatureMap.memory_gb",
        "PrepareRandomSampleFeatureMap.cpus",
        "PrepareRandomSampleFeatureMap.memory_gb",
        "TrainModel.xgboost_params_file",
        "MergeMd5sToJson.output_json",
        "Globals.glob"
    ]}
}    

parameter_meta {
    base_file_name: {
        help: "Base file name for output files. The output files will be named [base_file_name].with_ml_qual.vcf.gz",
        type: "String", 
        category: "input_required"
    }
    input_cram_bam: {
        help: "Input CRAM or BAM file",
        type: "File",
        category: "input_required"
    }
    input_cram_bam_index: {
        help: "Input CRAM or BAM index file",
        type: "File",
        category: "input_required"
    }
    sorter_json_stats_file: {
        help: "Sorter json stats file provided by the Ultima Genomics pipeline (same base name as the input cram/bam file with a json extension)",
        type: "File",
        category: "input_required"
    }
    references: {
        help: "Reference files: fasta, dict and fai, recommended value set in the template",
        type: "References",
        category: "ref_required"
    }
    featuremap_params: {
        type: "FeatureMapParams",
        help: "FeatureMap parameters, recommended value set in the template.",
        category: "param_required"
    }
    single_read_snv_params: {
        type: "SingleReadSNVParams",
        help: "SingleReadSNV parameters, recommended value set in the template.",
        category: "param_optional"
    }
    features: {
        type: "Array[String]",
        help: "Features to be used for training the SNV quality model, should match pre trained model if one is given, recommended value set in the template.",
        category: "param_required"
    }
    xgboost_params_file: {
        type: "File",
        help: "XGBoost parameters file for training the SNV quality model, recommended value set in the template.",
        category: "param_required"
    }
    create_md5_checksum_outputs: {
        help: "Create md5 checksum for requested output files",
        type: "Boolean",
        category: "input_optional"
    }
    training_regions_interval_list: {
        type: "File",
        help: "Genomic regions to include in the training set, the recommended value is set in the template",
        category: "param_required"
    }
    training_regions_interval_list_index: {
        type: "File",
        help: "Index for genomic regions to exclude from the training set, the recommended value is set in the template",
        category: "optional"
    }
    min_coverage_to_train_model: {
        type: "Float",
        help: "Minimum coverage to train the ML model, needed as label assignment (true/false SNV) is unreliable at low coverages, the recommended value is set in the template",
        category: "param_required"
    }
    pre_trained_model_files: {
        type: "Array[File]",
        help: "Pre-trained ML model json files, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features",
        category: "optional"
    }
    pre_trained_srsnv_metadata_json: {
        type: "File",
        help: "Pre-trained SNV quality model metadata json file, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features",
        category: "optional"
    }
    raise_exceptions_in_report: {
        type: "Boolean",
        help: "Raise and exception and fail the pipeline if an error is raised in the QC report",
        category: "optional"
    }
    no_address_override: {
        type: "Boolean",
        help: "Should the instances used be without external IP address. Allows for more parallelization, but not supported with Dockerhub dockers Default: true",
        category: "optional"
    }
    preemptible_tries_override: {
        type: "Int",
        help: "Number of times to retry preemptible instances, default: 1",
        category: "optional"
    }
    featuremap: {
        type: "File",
        help: "FeatureMap VCF file with SNVQ values",
        category: "output"
    }
    featuremap_index: {
        type: "File",
        help: "FeatureMap VCF index file",
        category: "output"
    }
    snv_qualities_assigned: {
        type: "Boolean",
        help: "Indicates whether SNV qualities were assigned (coverage was sufficient to train the model, over min_coverage_to_train_model), otherwise a FeatureMap with no quality scores is generated",
        category: "output"
    }
    used_self_trained_model: {
        type: "Boolean",
        help: "Indicates whether a self-trained model was used for inference, otherwise a pre-trained model was used (if snv_qualities_assigned) or no model was used",
        category: "output"
    }
    md5_checksums_json: {
        help: "json file that will contain md5 checksums for requested output files",
        type: "File",
        category: "output"
    }
    report_html: {
        type: "File",
        help: "SRSNV QC report html file",
        category: "output_optional"
    }
    model_file: {
        type: "File",
        help: "ML model file, saved with joblib",
        category: "output_optional"
    }
    params_file: {
        type: "File",
        help: "ML model parameters json file",
        category: "output_optional"
    }
    featuremap_df_file: {
        type: "File",
        help: "FeatureMap DataFrame with X matrix (features), y (labels) and qual, in parquet format",
        category: "output_optional"
    }
    test_set_mrd_simulation_dataframe: {
        type: "File",
        help: "ML model test set MRD simulation DataFrame, parquet format",
        category: "output_optional"
    }
    application_qc_h5: {
        type: "File",
        help: "Application QC statistics h5 file",
        category: "output_optional"
    }
    aggregated_metrics_json: {
        type: "File",
        help: "ML model aggregated metrics json file",
        category: "output_optional"
    }
    raw_filtered_featuremap_parquet: {
        type: "File",
        help: "Filtered parquet dataframe of raw featuremap for training",
        category: "output_optional"
    }
    raw_featuremap_stats: {
        type: "File",
        help: "Statistics for raw featuremap filtering",
        category: "output_optional"
    }
    random_sample_filtered_featuremap_parquet: {
        type: "File",
        help: "Filtered parquet dataframe of random sample featuremap for training",
        category: "output_optional"
    }
    random_sample_featuremap_stats: {
        type: "File",
        help: "Statistics for random sample featuremap filtering",
        category: "output_optional"
    }
    test_report_file_notebook: {
        type: "File",
        help: "ML model test set report notebook file",
        category: "output_optional"
    }
    featuremap_random_sample: {
        type: "File",
        help: "Downsampled FeatureMap VCF file for training",
        category: "output"
    }
    featuremap_random_sample_index: {
        type: "File",
        help: "Downsampled FeatureMap VCF index file",
        category: "output"
    }
    downsampling_rate: {
        type: "Float",
        help: "The downsampling rate used to create the random sample featuremap",
        category: "output"
    }
    featuremap_df: {
        type: "File",
        help: "Parquet file with model data, including features and labels",
        category: "output_optional"
    }
    model_files: {
        type: "Array[File]",
        help: "Array of model files, each file should be used for a specific set of chromosome",
        category: "output_optional"
    }
    srsnv_metadata_json: {
        type: "File",
        help: "Metadata JSON file for the SNV quality model, containing information about the model, features, and training parameters",
        category: "output_optional"
    }
}


  Int preemptibles = select_first([preemptible_tries, 1])
  String base_file_name_sub = sub(base_file_name, "#", "")
  Boolean no_address = select_first([no_address_override, true])

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

  Boolean use_pre_trained_model = defined(pre_trained_model_files) && defined(pre_trained_srsnv_metadata_json)

  call UGGeneralTasks.GetMeanCoverageFromSorterStats as GetMeanCoverageFromSorterStats {
    input:
      sorter_json_stats_file = select_first([sorter_json_stats_file]),
      docker = global.ugbio_srsnv_docker,
      preemptible_tries = preemptibles,
      monitoring_script = monitoring_script,  #!FileCoercion
  }
  Float mean_coverage = GetMeanCoverageFromSorterStats.mean_coverage

  # Calculate random_sample_size based on tp_train_set_size and tp_train_set_size_sampling_overhead
  Int random_sample_size = if (defined(single_read_snv_params)) then
    ceil(select_first([single_read_snv_params]).tp_train_set_size * select_first([single_read_snv_params]).tp_train_set_size_sampling_overhead)
  else
    1000000

  # Run snvfind to create featuremap
  call CreateFeatureMap {
    input:
      input_cram_bam = input_cram_bam,
      input_cram_bam_index = input_cram_bam_index,
      references = references,
      base_file_name = base_file_name_sub,
      sorter_json_stats_file = sorter_json_stats_file,
      random_sample_size = random_sample_size,
      featuremap_params = featuremap_params,
      docker = global.featuremap_docker,
      preemptible_tries = preemptibles,
      monitoring_script = monitoring_script,
  }

  Boolean sufficient_coverage_to_train_model = (mean_coverage >= min_coverage_to_train_model)
  Boolean snv_qualities_can_be_assigned = (sufficient_coverage_to_train_model) || (use_pre_trained_model)

  if ((sufficient_coverage_to_train_model) && (!use_pre_trained_model)) {
    # Prepare the raw featuremap for training
    SingleReadSNVParams single_read_snv_params_ = select_first([single_read_snv_params])
    # Dynamically assign memory for PrepareRawFeatureMap by mean_coverage
    # Int memory_gb_PrepareRawFeatureMap_default = if (mean_coverage < 40.0) then 16 else if (mean_coverage < 80.0) then 128 else 256
    # Int memory_gb_PrepareRawFeatureMap        = select_first([process_featuremap_memory_gb_override,
    #                                                           memory_gb_PrepareRawFeatureMap_default])
    Int memory_gb_PrepareRawFeatureMap = 128
    Int cpu_PrepareRawFeatureMap              = ceil(memory_gb_PrepareRawFeatureMap / 2)

    call PrepareFeatureMapForTraining as PrepareRawFeatureMap {
        input:
            featuremap = CreateFeatureMap.featuremap,
            featuremap_index = CreateFeatureMap.featuremap_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.tp_train_set_size,
            mean_coverage = mean_coverage,
            filters = ["name=vaf_le_5perc:field=RAW_VAF:op=le:value=" + single_read_snv_params_.max_vaf_for_fp + ":type=label"],
            docker = global.ugbio_featuremap_docker,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_PrepareRawFeatureMap,
            cpus      = cpu_PrepareRawFeatureMap,
            monitoring_script = monitoring_script
    }

    # Prepare the random sample featuremap for training
    Int memory_gb_PrepareRandomSampleFeatureMap = select_first([process_featuremap_memory_gb_override, 16])
    Int cpu_PrepareRandomSampleFeatureMap = ceil(memory_gb_PrepareRandomSampleFeatureMap / 2)

    call PrepareFeatureMapForTraining as PrepareRandomSampleFeatureMap {
        input:
            featuremap = CreateFeatureMap.featuremap_random_sample,
            featuremap_index = CreateFeatureMap.featuremap_random_sample_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.fp_train_set_size,
            mean_coverage = mean_coverage,
            filters = ["name=ref_eq_alt:field=REF:op=eq:value_field=ALT:type=label"],
            docker = global.ugbio_featuremap_docker,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_PrepareRandomSampleFeatureMap,
            cpus = cpu_PrepareRandomSampleFeatureMap,
            monitoring_script = monitoring_script
    }

    call PrepareFeatureMapForTraining as RandomSampleFeatureMapApplyNegativeFilter {
        input:
            featuremap = CreateFeatureMap.featuremap_random_sample,
            featuremap_index = CreateFeatureMap.featuremap_random_sample_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.fp_train_set_size,
            mean_coverage = mean_coverage,
            filters = [
              "name=ref_ne_alt:field=REF:op=ne:value_field=ALT:type=label",
              "name=vaf_le_5perc:field=RAW_VAF:op=le:value=" + single_read_snv_params_.max_vaf_for_fp + ":type=label"
              ],
            docker = global.ugbio_featuremap_docker,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_PrepareRandomSampleFeatureMap,
            cpus = cpu_PrepareRandomSampleFeatureMap,
            monitoring_script = monitoring_script
    }

    call TrainModel {
        input:
            raw_filtered_featuremap_parquet = PrepareRawFeatureMap.filtered_featuremap_parquet,
            random_sample_filtered_featuremap_parquet = PrepareRandomSampleFeatureMap.filtered_featuremap_parquet,
            random_sample_negative_label_stats = RandomSampleFeatureMapApplyNegativeFilter.stats_file,
            random_sample_positive_label_stats = PrepareRandomSampleFeatureMap.stats_file,
            raw_featuremap_stats = PrepareRawFeatureMap.stats_file,
            mean_coverage = mean_coverage,
            training_regions_interval_list = training_regions_interval_list,
            xgboost_params_file = xgboost_params_file,
            single_read_snv_params = single_read_snv_params_,
            features=features,
            base_file_name = base_file_name_sub,
            docker = global.ugbio_srsnv_docker,
            pipeline_version = pipeline_version,
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script
    }

    call CreateReport {
      input:
        featuremap_df        = TrainModel.featuremap_df,
        srsnv_metadata_json  = TrainModel.srsnv_metadata_json,
        model_files          = TrainModel.model_files,
        basename             = base_file_name_sub,
        docker               = global.ugbio_srsnv_docker,
        preemptible_tries    = preemptibles,
        monitoring_script    = monitoring_script
    }

    # Wire report outputs to outer scope
    File featuremap_df_output      = TrainModel.featuremap_df
    File application_qc_h5_output  = CreateReport.application_qc_h5
    File report_html_output        = CreateReport.report_html
  }

  if (snv_qualities_can_be_assigned) {
    Array[File] model_files_ = select_all(select_first([pre_trained_model_files, TrainModel.model_files]))
    File srsnv_metadata_json_ = select_first([pre_trained_srsnv_metadata_json, TrainModel.srsnv_metadata_json])

    call Inference {
        input:
          base_file_name      = base_file_name_sub,
          model_files         = model_files_,
          srsnv_metadata_json = srsnv_metadata_json_,
          featuremap          = CreateFeatureMap.featuremap,
          monitoring_script   = monitoring_script,  #!FileCoercion
          docker              = global.featuremap_docker,
          preemptible_tries   = preemptibles,
    }
  }
  File featuremap_output = select_first([Inference.featuremap_out, CreateFeatureMap.featuremap])
  File featuremap_index_output = select_first([Inference.featuremap_out_index, CreateFeatureMap.featuremap_index])
  

  if (create_md5_checksum_outputs) {

        Array[File] output_files = select_all(
          flatten(
            [
              select_first([[featuremap_output], []]),
              select_first([[featuremap_index_output], []]),
              select_first([[featuremap_df_output], []]),
            ]
          )
        )

        scatter (file in output_files) {
            call UGGeneralTasks.ComputeMd5 as compute_md5 {
                input:
                    input_file = file,
                    docker = global.ubuntu_docker,
            }
        }

        call UGGeneralTasks.MergeMd5sToJson {
            input:
                md5_files = compute_md5.checksum,
                docker = global.ugbio_core_docker
        }
    }

  output {
    File featuremap = featuremap_output
    File featuremap_index = featuremap_index_output
    File? featuremap_random_sample = CreateFeatureMap.featuremap_random_sample
    File? featuremap_random_sample_index = CreateFeatureMap.featuremap_random_sample_index
    Float downsampling_rate = CreateFeatureMap.downsampling_rate
    Boolean snv_qualities_assigned = snv_qualities_can_be_assigned
    Boolean used_self_trained_model = snv_qualities_assigned && (!use_pre_trained_model)
    # File? report_html = TrainModel.report_html
    # File? model_file = TrainModel.model
    # File? combined_stats = TrainModel.combined_stats
    File? raw_filtered_featuremap_parquet = PrepareRawFeatureMap.filtered_featuremap_parquet
    File? raw_featuremap_stats = PrepareRawFeatureMap.stats_file
    File? random_sample_filtered_featuremap_parquet = PrepareRandomSampleFeatureMap.filtered_featuremap_parquet
    File? random_sample_featuremap_stats = PrepareRandomSampleFeatureMap.stats_file

    File? featuremap_df = featuremap_df_output
    File? application_qc_h5 = application_qc_h5_output
    File? report_html = report_html_output
    File? srsnv_metadata_json = srsnv_metadata_json_
    Array[File]? model_files = model_files_

    File? md5_checksums_json = MergeMd5sToJson.md5_json

   }
}


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
  String intermediate_featuremap = base_file_name + ".training_regions.parquet"
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
    --drop-format GT AD

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
    Int memory_gb = 16
    Int cpus = 8
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
    File input_cram_bam
    File input_cram_bam_index
    References references
    String base_file_name
    
    File sorter_json_stats_file
    
    # Required training set size for downsampling calculation
    Int random_sample_size
    
    # All snvfind parameters are in the featuremap_params struct
    FeatureMapParams featuremap_params
    
    String docker
    Int preemptible_tries
    File monitoring_script
    Int memory_gb = 2
    Int cpus = 1
  }

  Float input_size = size(input_cram_bam, "GB")
  Float ref_size = size(references.ref_fasta, "GB")
  Int disk_size = ceil(input_size * 3 + ref_size + 20)

  String out_vcf = "~{base_file_name}.raw.featuremap.vcf.gz"
  String out_vcf_random_sample = "~{base_file_name}.random_sample.featuremap.vcf.gz"

  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    # Extract total aligned bases from sorter stats
    TOTAL_ALIGNED_BASES=$(jq -re '.total_aligned_bases // .total_bases // error("missing total_aligned_bases and total_bases")' "~{sorter_json_stats_file}")
    echo "Total aligned bases: $TOTAL_ALIGNED_BASES"
    DOWNSAMPLING_RATE=$(awk -v num=~{random_sample_size} -v den="$TOTAL_ALIGNED_BASES" 'BEGIN{printf "%.12f", num/den}')
    echo "Downsampling rate: $DOWNSAMPLING_RATE"

    # Run snvfind with parameters
    snvfind \
      ~{input_cram_bam} \
      ~{references.ref_fasta} \
      -o ~{out_vcf} \
      -f ~{out_vcf_random_sample},$DOWNSAMPLING_RATE \
      -v \
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
    Int memory_gb = 4
  }

  Float featuremap_size = size(featuremap, "GB")
  Int disk_size = ceil(featuremap_size * 2 + 20)

  String out_vcf = "~{base_file_name}.featuremap.vcf.gz"

  command <<<
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
    Int    memory_gb = 16
    Int    cpus      = 2
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