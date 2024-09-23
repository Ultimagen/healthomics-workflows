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

# DESCRIPTION
# Runs substitution error analysis for Ultima Genomics data. Includes three subworkflows:
# 1. FeatureMap generation - aggregates all the substitutions from a bam file
# 2. Annotation with additional features
# 3. Single read ML model training and inference

# CHANGELOG in reverse chronological order
# 1.10.1 Added LA_7 adapter version
# 1.7.0 Initial implementation of single_read_snv wdl

import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "featuremap.wdl" as FeaturemapSubWF
import "tasks/qc_tasks.wdl" as UGQCTasks
import "tasks/mrd.wdl" as UGMrdTasks
import "tasks/globals.wdl" as Globals


workflow SingleReadSNV {
input {
  File input_cram_bam
  File input_cram_bam_index
  File sorter_json_stats_file
  String base_file_name
  Array[File]? somatic_mutations_list
  String pipeline_version = "1.14.2" 
  References references

  File wgs_calling_interval_list  # TODO update this name to interval_list
  Int break_bands_at_multiples_of

  FeatureMapParams featuremap_params
  SingleReadSNVParams single_read_snv_params
  Map[String, Array[String]] categorical_features
  Int? featuremap_scatter_count_override

  # snv quality model include/exclude regions
  Array[File] training_include_regions
  Array[File]? tp_training_exclude_regions
  Array[File]? fp_training_exclude_regions

  Float min_coverage_to_train_model

  File? pre_trained_model_file

  Boolean raise_exceptions_in_report

  Int? process_featuremap_memory_gb_override
  Int? preemptible_tries
  Boolean? no_address_override
  String? cloud_provider_override
  # Used for running on other clouds (aws)
  File? monitoring_script_input

  # winval validations
  #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
  #@wv not("test" in base_file_name or "train" in base_file_name)
  #@wv suffix(input_cram_bam) in {".bam", ".cram"}
  #@wv suffix(input_cram_bam_index) in {".bai", ".crai"}
  #@wv prefix(input_cram_bam_index) == input_cram_bam
  #@wv suffix(sorter_json_stats_file) == ".json"
  #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
  #@wv suffix(references['ref_dict']) == '.dict'
  #@wv suffix(references['ref_fasta_index']) == '.fai'
  #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
  #@wv featuremap_params['motif_length_to_annotate'] <= 4 and featuremap_params['motif_length_to_annotate'] >= 1
  #@wv single_read_snv_params['ppmSeq_adapter_version'] in ["None", "v1", "legacy_v5", "legacy_v5_start", "legacy_v5_end", "dmbl"]
  #@wv single_read_snv_params['num_CV_folds'] >= 1
  #@wv single_read_snv_params['split_folds_by'] in ['random', 'chromosome', 'chrom']

  # pre-trained model parameters
  #@wv defined(pre_trained_model_file) -> suffix(pre_trained_model_file) == ".joblib"
  
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
        "FeatureMap.pipeline_version",
        "CreateHomSnvFeatureMap.disk_size",
        "CreateTpTrainingRegionsBed.disk_size",
        "CreateTpTrainingRegionsBed.cpus",
        "CreateFpTrainingRegionsBed.disk_size",
        "CreateFpTrainingRegionsBed.cpus",
        "TrainSnvQualityRecalibrationModel.cpus",
        "TrainSnvQualityRecalibrationModel.memory_gb",
        "TrainSnvQualityRecalibrationModel.disk_size",
        "InferenceSnvQualityRecalibrationModel.cpus",
        "InferenceSnvQualityRecalibrationModel.memory_gb",
        "InferenceSnvQualityRecalibrationModel.disk_size",
        "FeatureMap.Globals.glob",
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
    somatic_mutations_list: {
        help: "Somatic mutations to be excluded from FP training set, will be appended to the fp_training_exclude_regions optional",
        type: "Array[File]",
        category: "optional"
    }
    references: {
        help: "Reference files: fasta, dict and fai, recommended value set in the template",
        type: "References",
        category: "ref_required"
    }
    wgs_calling_interval_list: {
        type: "File",
        help: "interval list defining the region to perform variant calling on, recommended value set in the template",
        category: "param_required"
    }
    break_bands_at_multiples_of: {
        type: "Int",
        help: "Break wgs_calling_interval_list bands at multiples of this number, recommended value set in the template",
        category: "param_required"
    }
    featuremap_params: {
        type: "FeatureMapParams",
        help: "FeatureMap parameters, recommended value set in the template.",
        category: "param_required"
    }
    single_read_snv_params: {
        type: "SingleReadSNVParams",
        help: "SingleReadSNV parameters, recommended value set in the template.",
        category: "param_required"
    }
    categorical_features: {
        type: "Map[String, Array[String]]",
        help: "Categorical features in SingleReadSNV model, list of feature names with allowed values per feature. Separate from single_read_snv_params due to technical reasons. The recommended value is set in the template.",
        category: "param_required"
    }
    training_include_regions: {
        type: "Array[File]",
        help: "Genomic regions to include in the training set, the recommended value is set in the template",
        category: "param_required"
    }
    tp_training_exclude_regions: {
        type: "Array[File]",
        help: "Genomic regions to exclude from the training set TP examples, the recommended value is set in the template",
        category: "optional"
    }
    fp_training_exclude_regions: {
        type: "Array[File]",
        help: "Genomic regions to exclude from the training set FP examples, the recommended value is set in the template",
        category: "optional"
    }
    min_coverage_to_train_model: {
        type: "Float",
        help: "Minimum coverage to train the ML model, needed as label assignment (true/false SNV) is unreliable at low coverages, the recommended value is set in the template",
        category: "param_required"
    }
    pre_trained_model_file: {
        type: "File",
        help: "Pre-trained ML model file, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features",
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
    test_set_statistics_h5: {
        type: "File",
        help: "ML model test set statistics h5 file",
        category: "output_optional"
    }
    aggregated_metrics_json: {
        type: "File",
        help: "ML model aggregated metrics json file",
        category: "output_optional"
    }
    tp_training_regions_bed: {
        type: "File",
        help: "ML model training set TP regions bed file",
        category: "output_optional"
    }
    fp_training_regions_bed: {
        type: "File",
        help: "ML model training set FP regions bed file",
        category: "output_optional"
    }
    test_report_file_notebook: {
        type: "File",
        help: "ML model test set report notebook file",
        category: "output_optional"
    }
    flow_order: {
        type: "String",
        help: "Flow order for the sample",
        category: "output_optional"
    }
}


  Int preemptibles = select_first([preemptible_tries, 1])
  String base_file_name_sub = sub(base_file_name, "#", "")
  Boolean no_address = select_first([no_address_override, true ])

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

  Boolean use_pre_trained_model = defined(pre_trained_model_file)

  call UGGeneralTasks.ExtractSampleNameFlowOrder as ExtractSampleNameFlowOrder{
    input:
      input_bam = input_cram_bam,
      references = references,
      preemptible_tries = preemptibles,
      monitoring_script = monitoring_script,
      docker = global.broad_gatk_docker,
      no_address = no_address,
      cloud_provider_override = cloud_provider_override
  }

  call UGGeneralTasks.GetMeanCoverageFromSorterStats as GetMeanCoverageFromSorterStats {
    input:
      sorter_json_stats_file = sorter_json_stats_file,
      docker = global.ug_vc_docker,
      preemptible_tries = preemptibles,
      monitoring_script = monitoring_script,  #!FileCoercion
  }
  Float mean_coverage = GetMeanCoverageFromSorterStats.mean_coverage
  Int featuremap_shard_number_calc = select_first([featuremap_scatter_count_override, ceil(mean_coverage / 1.4)])
  # patch because there is no "max" function...
  if (featuremap_shard_number_calc<2) {
    Int featuremap_shard_number_min = 2
  }
  Int featuremap_shard_number = select_first([featuremap_shard_number_min, featuremap_shard_number_calc])
  if (single_read_snv_params.ppmSeq_adapter_version != "None") {
    String ppmSeq_adapter_version_ = single_read_snv_params.ppmSeq_adapter_version
  }
  call FeaturemapSubWF.FeatureMap {
    input:
      input_cram_bam = input_cram_bam,
      input_cram_bam_index = input_cram_bam_index,
      references = references,
      wgs_calling_interval_list = wgs_calling_interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      featuremap_params = featuremap_params,
      scatter_count = featuremap_shard_number,
      flow_order = ExtractSampleNameFlowOrder.flow_order,
      ppmSeq_adapter_version = ppmSeq_adapter_version_,
      base_file_name = base_file_name_sub,
      preemptible_tries = preemptibles,
      process_featuremap_memory_gb_override = process_featuremap_memory_gb_override,
  }

  Boolean sufficient_coverage_to_train_model = (mean_coverage >= min_coverage_to_train_model)
  Boolean snv_qualities_can_be_assigned = (sufficient_coverage_to_train_model) || (use_pre_trained_model)

  if ((sufficient_coverage_to_train_model) && (!use_pre_trained_model)) {
    call UGMrdTasks.CreateHomSnvFeatureMap as CreateHomSnvFeatureMap {
        input:
        featuremap =              FeatureMap.featuremap,
        featuremap_index =        FeatureMap.featuremap_index,
        sorter_json_stats_file =  sorter_json_stats_file,
        base_file_name =          base_file_name_sub,
        min_af =                  0.7,
        min_coverage =            20,
        memory_gb =               2,
        cpus =                    1,
        docker =                  global.ug_vc_docker,
        preemptibles =            preemptibles,
        monitoring_script =       monitoring_script,  #!FileCoercion
    }

    call UGMrdTasks.BedIntersectAndExclude as CreateTpTrainingRegionsBed {
        input:
            include_regions = training_include_regions,
            exclude_regions = tp_training_exclude_regions,
            output_basename = "TP",
            memory_gb = 4,
            docker = global.ug_vc_docker,
            monitoring_script = monitoring_script,  #!FileCoercion
            preemptibles = preemptibles
    }


    call UGMrdTasks.BedIntersectAndExclude as CreateFpTrainingRegionsBed {
        input:
            include_regions = training_include_regions,
            exclude_regions = flatten([select_first([fp_training_exclude_regions, []]), select_first([somatic_mutations_list, []])]),
            output_basename = "FP",
            memory_gb = 8,
            cpus = 4,
            docker = global.ug_vc_docker,
            monitoring_script = monitoring_script,  #!FileCoercion
            preemptibles = preemptibles
    }

  call UGMrdTasks.TrainSnvQualityRecalibrationModel as TrainSnvQualityRecalibrationModel {
    input:
      basename                        = base_file_name_sub,
      hom_snv_featuremap              = CreateHomSnvFeatureMap.hom_snv_featuremap,
      hom_snv_featuremap_index        = CreateHomSnvFeatureMap.hom_snv_featuremap_index,
      singleton_snv_featuremap        = FeatureMap.featuremap_single_substitutions,
      singleton_snv_featuremap_index  = FeatureMap.featuremap_single_substitutions_index,
      sorter_json_stats_file          = sorter_json_stats_file,
      single_read_snv_params          = single_read_snv_params,
      categorical_features            = categorical_features,
      hom_snv_regions_bed             = CreateTpTrainingRegionsBed.merged_bed,
      single_substitution_regions_bed = CreateFpTrainingRegionsBed.merged_bed,
      references                      = references,
      raise_exceptions_in_report      = raise_exceptions_in_report,
      flow_order                      = ExtractSampleNameFlowOrder.flow_order,
      monitoring_script               = monitoring_script,  #!FileCoercion
      docker                          = global.ug_vc_docker,
      pipeline_version                = pipeline_version,
      preemptible_tries               = preemptibles,
    }

  }
  if (snv_qualities_can_be_assigned) {
    File srsnv_model_file = select_first([pre_trained_model_file, TrainSnvQualityRecalibrationModel.model_file])

    call UGMrdTasks.InferenceSnvQualityRecalibrationModel as InferenceSnvQualityRecalibrationModel {
        input:
        output_file                             = base_file_name_sub + ".featuremap_with_qual.vcf.gz",
        model_file                              = srsnv_model_file,
        featuremap                              = FeatureMap.featuremap,
        featuremap_index		                = FeatureMap.featuremap_index,
        monitoring_script                       = monitoring_script,  #!FileCoercion
        docker                                  = global.ug_vc_docker,
        preemptible_tries                       = preemptibles,
    }
  }

  output {
    File featuremap = select_first([InferenceSnvQualityRecalibrationModel.featuremap_with_qual, FeatureMap.featuremap])
    File featuremap_index = select_first([InferenceSnvQualityRecalibrationModel.featuremap_with_qual_index, FeatureMap.featuremap_index])
    Boolean snv_qualities_assigned = snv_qualities_can_be_assigned
    Boolean used_self_trained_model = snv_qualities_can_be_assigned && (!use_pre_trained_model)
    File? report_html = TrainSnvQualityRecalibrationModel.test_report_file_html
    File? featuremap_df_file = TrainSnvQualityRecalibrationModel.featuremap_df_file
    File? model_file = TrainSnvQualityRecalibrationModel.model_file  
    File? params_file = TrainSnvQualityRecalibrationModel.params_file
    File? test_set_mrd_simulation_dataframe = TrainSnvQualityRecalibrationModel.test_set_mrd_simulation_dataframe
    File? test_set_statistics_h5 = TrainSnvQualityRecalibrationModel.test_set_statistics_h5
    File? aggregated_metrics_json = TrainSnvQualityRecalibrationModel.test_set_statistics_json

    File? tp_training_regions_bed = CreateTpTrainingRegionsBed.merged_bed
    File? fp_training_regions_bed =  CreateFpTrainingRegionsBed.merged_bed

    File? test_report_file_notebook =  TrainSnvQualityRecalibrationModel.test_report_file_notebook
   }
}

