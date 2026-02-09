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
import "tasks/single_read_snv_tasks.wdl" as SRSNVTasks
import "tasks/qc_tasks.wdl" as UGQCTasks
import "tasks/mrd.wdl" as UGMrdTasks
import "tasks/globals.wdl" as Globals


workflow SingleReadSNV {
input {
  Array[File] input_cram_bam_list
  Array[File] input_cram_bam_index_list
  Array[File]? sorter_json_stats_file_list
  String base_file_name
  String pipeline_version = "1.26.1"
  References references

  FeatureMapParams featuremap_params
  SingleReadSNVParams? single_read_snv_params
  Array[String] features
  Int? featuremap_scatter_count_override

  # snv quality model include/exclude regions
  File training_regions_interval_list
  File training_regions_interval_list_index
  File xgboost_params_file
  File? random_sample_trinuc_freq

  Float min_coverage_to_train_model

  Array[File]? pre_trained_model_files
  File? pre_trained_srsnv_metadata_json

  Boolean raise_exceptions_in_report

  Int? prepare_featuremap_for_training_gb
  Int? training_memory_gb
  Int? preemptible_tries
  Boolean? no_address_override
  String? cloud_provider_override
  # Used for running on other clouds (aws)
  File? monitoring_script_input

  Boolean create_md5_checksum_outputs = false

  Float? mean_coverage
  String? total_aligned_bases
  # winval validations
  #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
  #@wv not("test" in base_file_name or "train" in base_file_name)
  #@wv suffix(input_cram_bam_list) <= {".bam", ".cram"}
  #@wv suffix(input_cram_bam_index_list) <= {".bai", ".crai"}
  #@wv len(input_cram_bam_list) == len(input_cram_bam_index_list)
  ##@wv prefix(basename(input_cram_bam_index_list)) == basename(input_cram_bam_list)
  ##@wv prefix(basename(training_regions_interval_list_index)) == basename(training_regions_interval_list)
  #@wv suffix(training_regions_interval_list) in {".gz"}
  #@wv suffix(prefix(training_regions_interval_list)) in {".interval_list"}
  #@wv defined(sorter_json_stats_file_list) -> suffix(sorter_json_stats_file_list) <= {".json"}
  #@wv defined(sorter_json_stats_file_list) -> len(sorter_json_stats_file_list) == len(input_cram_bam_list)
  #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
  #@wv suffix(references['ref_dict']) == '.dict'
  #@wv suffix(references['ref_fasta_index']) == '.fai'
  #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
  #@wv not defined(pre_trained_model_files) -> defined(single_read_snv_params)
  #@wv defined(single_read_snv_params) -> single_read_snv_params['num_CV_folds'] > 1
  #@wv defined(pre_trained_model_files) <-> defined(pre_trained_srsnv_metadata_json)
  ##@wv defined(pre_trained_model_files) -> suffix(pre_trained_model_files) == ".json"
  #@wv defined(pre_trained_srsnv_metadata_json) -> suffix(pre_trained_srsnv_metadata_json) == ".json"
  #@wv (defined(mean_coverage) or defined(total_aligned_bases)) -> (defined(mean_coverage) and defined(total_aligned_bases) and not defined(sorter_json_stats_file_list))
  #@wv defined(sorter_json_stats_file_list) -> (not defined(mean_coverage) and not defined(total_aligned_bases))
  #@wv defined(random_sample_trinuc_freq) -> suffix(random_sample_trinuc_freq) in {".csv", ".tsv"}
}

meta {
    description : "The Single Read SNV (SRSNV) pipeline is a read-centric de-noising framework, developed to overcome the limitations of traditional locus-centric variant calling, particularly in scenarios where rare mutations may be supported by only a single read. These rare mutations need to be distinguished from artefactual SNVs, which can derive from sequencing, library or alignment errors. To achieve this, we employed a supervised machine learning model trained to classify actual SNVs (labelled True or TP) from noise (False or FP). First, a comprehensive dataset capturing every candidate SNV is generated, along with a rich suite of annotations that describe sequencing quality, local sequence motifs, fragment-specific features, and locus-specific information. Randomly selected bases in the data matching the reference genome are collected and annotated as True, while low VAF (<=5%) SNVs in high-coverage (>=20x) regions (SNVs with low support, indicating they are likely to be artifacts) are annotated as False SNVs. Using these curated sets, we train an XGBoost classifier to robustly distinguish between true and artifactual SNVs. Once trained, the classifier assigns a calibrated quality score to each SNV in the input CRAM, providing a precise estimate of the residual error rate. To avoid overfitting, an ensemble of models (3) are trained on different sets of chromosomes and applied using a cross-validation scheme. \n\nThe following input templates are available for different kinds of input data: \n\n1) `single_read_snv_template-ppmSeq.json` | Use this template for ppmSeq data. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. st, et). \n\n2) `single_read_snv_template-ppmSeq_legacy_v5.json` | Use this template for LEGACY v5 ppmSeq data. This is an older version of the ppmSeq adapters, generally not available since 2024. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. as, ts). \n\n3) `single_read_snv_template-Standard-WG.json` | Use this template for any non-ppmSeq data. "
    author: "Ultima Genomics"
    WDL_AID: { exclude: [
        "pipeline_version",
        "cloud_provider_override",
        "no_address_override",
        "preemptible_tries",
        "monitoring_script_input",
        "prepare_featuremap_for_training_gb",
        "featuremap_scatter_count_override",
        "CreateTrainingRegionsBed.disk_size",
        "CreateTrainingRegionsBed.cpus",
        "TrainModel.cpus",
        "TrainModel.memory_gb",
        "TrainModel.disk_size",
        "training_memory_gb",
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
    input_cram_bam_list: {
        help: "Input CRAM file/s",
        type: "Array[File]",
        category: "input_required"
    }
    input_cram_bam_index_list: {
        help: "Input CRAM index file/s",
        type: "Array[File]",
        category: "input_required"
    }
    sorter_json_stats_file_list: {
        help: "(Optional) Sorter json stats files. Provide EITHER these files OR both mean_coverage and total_aligned_bases.",
        type: "Array[File]",
        category: "input_optional"
    }
    mean_coverage: {
        help: "(Optional) Mean coverage value. Provide together with total_aligned_bases and without sorter_json_stats_file_list.",
        type: "Float",
        category: "input_optional"
    }
    total_aligned_bases: {
        help: "(Optional) Total aligned bases used for downsampling rate calculation. Provide together with mean_coverage and without sorter_json_stats_file_list.",
        type: "String",
        category: "input_optional"
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
    random_sample_trinuc_freq: {
        type: "File",
        help: "(Optional) CSV or TSV file with trinucleotide frequencies for the random sample. If provided, the random sample featuremap will be sampled according to the given trinucleotide frequency. If not provided, sampling is uniform.",
        category: "input_optional"
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
    monitoring_script_input: {
        type: "File",
        help: "Monitoring script override for AWS HealthOmics workflow templates multi-region support",
        category: "input_optional"
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
    random_sample_trinuc_freq_stats: {
        type: "File",
        help: "(Optional) CSV or TSV file with trinucleotide frequencies for the random sample.",
        category: "output_optional"
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

  if (defined(sorter_json_stats_file_list)) {
    Array[File] sorter_json_stats_file_list_ = select_first([sorter_json_stats_file_list])
    call UGGeneralTasks.ExtractSorterStatsMetrics {
      input:
        sorter_json_stats_files = sorter_json_stats_file_list_,
        docker = global.ugbio_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script,
    }
  }
  Float mean_coverage_used = select_first([mean_coverage, ExtractSorterStatsMetrics.mean_coverage])
  String total_aligned_bases_used = select_first([total_aligned_bases, ExtractSorterStatsMetrics.total_aligned_bases])

  # Calculate random_sample_size based on tp_train_set_size and tp_train_set_size_sampling_overhead
  Int random_sample_size = if (defined(single_read_snv_params)) then
    ceil(select_first([single_read_snv_params]).tp_train_set_size * select_first([single_read_snv_params]).tp_train_set_size_sampling_overhead)
  else
    1000000

  # Run snvfind to create featuremap
  call SRSNVTasks.CreateFeatureMap {
    input:
      input_cram_bam_list = input_cram_bam_list,
      input_cram_bam_index_list = input_cram_bam_index_list,
      references = references,
      base_file_name = base_file_name_sub,
      total_aligned_bases = total_aligned_bases_used,
      random_sample_size = random_sample_size,
      random_sample_trinuc_freq_ = random_sample_trinuc_freq,
      featuremap_params = featuremap_params,
      docker = global.featuremap_docker,
      preemptible_tries = preemptibles,
      monitoring_script = monitoring_script,
  }

  Boolean sufficient_coverage_to_train_model = (mean_coverage_used >= min_coverage_to_train_model)
  Boolean snv_qualities_can_be_assigned = (sufficient_coverage_to_train_model) || (use_pre_trained_model)

  if ((sufficient_coverage_to_train_model) && (!use_pre_trained_model)) {
    # Prepare the raw featuremap for training
    SingleReadSNVParams single_read_snv_params_ = select_first([single_read_snv_params])
    # Dynamically assign memory for PrepareRawFeatureMap by mean_coverage
    # Int memory_gb_PrepareRawFeatureMap_default = if (mean_coverage < 40.0) then 16 else if (mean_coverage < 80.0) then 128 else 256
    # Int memory_gb_PrepareRawFeatureMap        = select_first([prepare_featuremap_for_training_gb,
    #                                                           memory_gb_PrepareRawFeatureMap_default])
    Int memory_gb_PrepareRawFeatureMap = select_first([prepare_featuremap_for_training_gb, 128])
    Int cpu_PrepareRawFeatureMap_ = ceil(memory_gb_PrepareRawFeatureMap / 2)
    Int cpu_PrepareRawFeatureMap = if (cpu_PrepareRawFeatureMap_ > 64) then 64 else cpu_PrepareRawFeatureMap_

    call SRSNVTasks.PrepareFeatureMapForTraining as PrepareRawFeatureMap {
        input:
            featuremap = CreateFeatureMap.featuremap,
            featuremap_index = CreateFeatureMap.featuremap_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.tp_train_set_size,
            mean_coverage = mean_coverage_used,
            filters = ["name=low_vaf:field=RAW_VAF:op=le:value=" + single_read_snv_params_.max_vaf_for_fp + ":type=label"],
            docker = global.ugbio_featuremap_docker,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_PrepareRawFeatureMap,
            cpus      = cpu_PrepareRawFeatureMap,
            monitoring_script = monitoring_script
    }

    # Prepare the random sample featuremap for training
    Int memory_gb_PrepareRandomSampleFeatureMap = select_first([prepare_featuremap_for_training_gb, 16])
    Int cpu_PrepareRandomSampleFeatureMap_ = ceil(memory_gb_PrepareRandomSampleFeatureMap / 2)
    Int cpu_PrepareRandomSampleFeatureMap = if (cpu_PrepareRandomSampleFeatureMap_ > 8) then 8 else cpu_PrepareRandomSampleFeatureMap_

    call SRSNVTasks.PrepareFeatureMapForTraining as PrepareRandomSampleFeatureMap {
        input:
            featuremap = CreateFeatureMap.featuremap_random_sample,
            featuremap_index = CreateFeatureMap.featuremap_random_sample_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.fp_train_set_size,
            mean_coverage = mean_coverage_used,
            filters = ["name=ref_eq_alt:field=REF:op=eq:value_field=ALT:type=label"],
            docker = global.ugbio_featuremap_docker,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_PrepareRandomSampleFeatureMap,
            cpus = cpu_PrepareRandomSampleFeatureMap,
            monitoring_script = monitoring_script
    }

    call SRSNVTasks.PrepareFeatureMapForTraining as RandomSampleFeatureMapApplyNegativeFilter {
        input:
            featuremap = CreateFeatureMap.featuremap_random_sample,
            featuremap_index = CreateFeatureMap.featuremap_random_sample_index,
            training_regions_interval_list = training_regions_interval_list,
            training_regions_interval_list_index = training_regions_interval_list_index,
            featuremap_params = featuremap_params,
            single_read_snv_params = single_read_snv_params_,
            train_set_size = single_read_snv_params_.fp_train_set_size,
            mean_coverage = mean_coverage_used,
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

    Int memory_gb_TrainModel = select_first([training_memory_gb, 32])
    call SRSNVTasks.TrainModel {
        input:
            raw_filtered_featuremap_parquet = PrepareRawFeatureMap.filtered_featuremap_parquet,
            random_sample_filtered_featuremap_parquet = PrepareRandomSampleFeatureMap.filtered_featuremap_parquet,
            random_sample_negative_label_stats = RandomSampleFeatureMapApplyNegativeFilter.stats_file,
            random_sample_positive_label_stats = PrepareRandomSampleFeatureMap.stats_file,
            raw_featuremap_stats = PrepareRawFeatureMap.stats_file,
            mean_coverage = mean_coverage_used,
            training_regions_interval_list = training_regions_interval_list,
            xgboost_params_file = xgboost_params_file,
            single_read_snv_params = single_read_snv_params_,
            features=features,
            base_file_name = base_file_name_sub,
            docker = global.ugbio_srsnv_docker,
            pipeline_version = pipeline_version,
            preemptible_tries = preemptibles,
            memory_gb = memory_gb_TrainModel,
            monitoring_script = monitoring_script
    }

    call SRSNVTasks.CreateReport {
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

    call SRSNVTasks.Inference {
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
    File? random_sample_trinuc_freq_stats = CreateFeatureMap.random_sample_trinuc_freq

    File? featuremap_df = featuremap_df_output
    File? application_qc_h5 = application_qc_h5_output
    File? report_html = report_html_output
    File? srsnv_metadata_json = srsnv_metadata_json_
    Array[File]? model_files = model_files_

    File? md5_checksums_json = MergeMd5sToJson.md5_json

   }
}
