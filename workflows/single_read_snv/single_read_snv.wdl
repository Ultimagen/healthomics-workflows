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
  String pipeline_version = "1.11.4" # !UnusedDeclaration
  References references

  File wgs_calling_interval_list  # TODO update this name to interval_list
  Int break_bands_at_multiples_of

  FeatureMapParams featuremap_params

  # featuremap annotation parameters
  String? balanced_strand_adapter_version
  Int motif_length_to_annotate
  Int max_hmer_length
  # snv quality recalibration params
  Int train_set_size
  Int test_set_size
  Array[File] tp_training_include_regions
  Array[File]? tp_training_exclude_regions
  Array[File] fp_training_include_regions
  Array[File]? fp_training_exclude_regions

  Array[String] numerical_features
  Array[String] categorical_features
  Array[String]? balanced_sampling_info_fields
  String? pre_filter
  Int random_seed

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
  #@wv defined(balanced_strand_adapter_version) -> balanced_strand_adapter_version in ['LA_v5', 'LA_v6', 'LA_v5and6', 'LA_v7']
  #@wv motif_length_to_annotate <= 4 and motif_length_to_annotate >= 1
}

meta {
    description : "Single Read SNV Quality Recalibration workflow (single_read_snv wdl) is a software tool that assigns accurate quality scores to all SNV candidates. The output is a FeatureMap VCF file with the recalibrated SNV quality scores."
    author: "Ultima Genomics"
    WDL_AID: { exclude: [
        "pipeline_version",
        "cloud_provider_override",
        "no_address_override",
        "preemptible_tries",
        "monitoring_script_input",
        "process_featuremap_memory_gb_override",
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
        help: "FeatureMap parameters, recommended value set in the template. Int scatter_count: number of scatter tasks to use. Int min_mapq, Int snv_identical_bases, Int snv_identical_bases_after, Int min_score, Int limit_score, String extra_args",
        category: "param_required"
    }
    balanced_strand_adapter_version: {
        type: "String",
        help: "ppmSeq adapter version, for ppmSeq data the recommended value is set in the template",
        category: "param_required"
    }
    motif_length_to_annotate: {
        type: "Int",
        help: "Length of the motif (-+N bp) to annotate in the FeatureMap, the recommended value is set in the template",
        category: "param_required"
    }
    max_hmer_length: {
        type: "Int",
        help: "Maximum length of the homopolymer to annotate in the FeatureMap, the recommended value is set in the template",
        category: "param_required"
    }
    train_set_size: {
        type: "Int",
        help: "Number of SNVs to use for the ML model training set, the recommended value is set in the template",
        category: "param_required"
    }
    test_set_size: {
        type: "Int",
        help: "Number of SNVs to use for the ML model test set, the recommended value is set in the template",
        category: "param_required"
    }
    tp_training_include_regions: {
        type: "Array[File]",
        help: "Genomic regions to include in the training set TP examples, the recommended value is set in the template",
        category: "param_required"
    }
    tp_training_exclude_regions: {
        type: "Array[File]",
        help: "Genomic regions to exclude from the training set TP examples, the recommended value is set in the template",
        category: "optional"
    }
    fp_training_include_regions: {
        type: "Array[File]",
        help: "Genomic regions to include in the training set FP examples, the recommended value is set in the template",
        category: "param_required"
    }
    fp_training_exclude_regions: {
        type: "Array[File]",
        help: "Genomic regions to exclude from the training set FP examples, the recommended value is set in the template",
        category: "optional"
    }
    numerical_features: {
        type: "Array[String]",
        help: "Numerical features to use in the ML model, the recommended value is set in the template",
        category: "param_required"
    }
    categorical_features: {
        type: "Array[String]",
        help: "Categorical features to use in the ML model, the recommended value is set in the template",
        category: "param_required"
    }
    balanced_sampling_info_fields: {
        type: "Array[String]",
        help: "Fields to use for balanced sampling of TP examples to remove the prior distribution of the homozygous SNVs, the pipeline will attempt for the distribution over these arguments to be uniform. The recommended value is set in the template",
        category: "optional"
    }
    pre_filter: {
        type: "String",
        help: "SNV filter to apply to the FeatureMap before model, any SNV not passing the pre_filter is assigned QUAL=0, the recommended value is set in the template",
        category: "optional"
    }
    random_seed: {
        type: "Int",
        help: "Random seed to use for the ML model, the recommended value is set in the template",
        category: "advanced"
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
    report_html: {
        type: "File",
        help: "SRSNV QC report html file",
        category: "output"
    }
    model_file: {
        type: "File",
        help: "ML model file, saved with joblib",
        category: "output"
    }
    X_test_file: {
        type: "File",
        help: "ML model test set features DataFrame, parquet format",
        category: "output"
    }
    y_test_file: {
        type: "File",
        help: "ML model test set labels DataFrame, parquet format",
        category: "output"
    }
    qual_test_file: {
        type: "File",
        help: "ML model test set qual (SNVQ) DataFrame, parquet format",
        category: "output"
    }
    X_train_file: {
        type: "File",
        help: "ML model training set features DataFrame, parquet format",
        category: "output"
    }
    y_train_file: {
        type: "File",
        help: "ML model training set labels DataFrame, parquet format",
        category: "output"
    }
    params_file: {
        type: "File",
        help: "ML model parameters json file",
        category: "output"
    }
    test_set_mrd_simulation_dataframe: {
        type: "File",
        help: "ML model test set MRD simulation DataFrame, parquet format",
        category: "output"
    }
    train_set_mrd_simulation_dataframe: {
        type: "File",
        help: "ML model training set MRD simulation DataFrame, parquet format",
        category: "output"
    }
    test_set_statistics_h5: {
        type: "File",
        help: "ML model test set statistics h5 file",
        category: "output"
    }
    train_set_statistics_h5: {
        type: "File",
        help: "ML model training set statistics h5 file",
        category: "output"
    }
    train_set_statistics_json: {
        type: "File",
        help: "ML model training set statistics json file",
        category: "output"
    }
    aggregated_metrics_json: {
        type: "File",
        help: "ML model aggregated metrics json file",
        category: "output"
    }
    tp_training_regions_bed: {
        type: "File",
        help: "ML model training set TP regions bed file",
        category: "output"
    }
    fp_training_regions_bed: {
        type: "File",
        help: "ML model training set FP regions bed file",
        category: "output"
    }
    train_report_file_notebook: {
        type: "File",
        help: "ML model training set report notebook file",
        category: "output"
    }
    train_report_file_html: {
        type: "File",
        help: "ML model training set report html file",
        category: "output"
    }
    test_report_file_notebook: {
        type: "File",
        help: "ML model test set report notebook file",
        category: "output"
    }
    flow_order: {
        type: "String",
        help: "Flow order for the sample",
        category: "output"
    }

        
    }


  Int preemptibles = select_first([preemptible_tries, 1])
  String base_file_name_sub = sub(base_file_name, "#", "")
  Boolean no_address = select_first([no_address_override, true ])

  call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

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
  Int mean_coverage = GetMeanCoverageFromSorterStats.mean_coverage
  Int featuremap_shard_number = ceil(mean_coverage / 1.4)   # overrides featuremap_params.scatter_count

  call FeaturemapSubWF.FeatureMap {
    input:
      input_cram_bam = input_cram_bam,
      input_cram_bam_index = input_cram_bam_index,
      references = references,
      wgs_calling_interval_list = wgs_calling_interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      featuremap_params = featuremap_params,
      featuremap_shard_number = featuremap_shard_number,
      flow_order = ExtractSampleNameFlowOrder.flow_order,
      balanced_strand_adapter_version = balanced_strand_adapter_version,
      motif_length_to_annotate = motif_length_to_annotate,
      max_hmer_length = max_hmer_length,
      base_file_name = base_file_name_sub,
      preemptible_tries = preemptibles,
      process_featuremap_memory_gb_override = process_featuremap_memory_gb_override,
  }

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
          include_regions = tp_training_include_regions,
          exclude_regions = tp_training_exclude_regions,
          output_basename = "TP",
          memory_gb = 8,
          docker = global.ug_vc_docker,
          monitoring_script = monitoring_script,  #!FileCoercion
          preemptibles = preemptibles
  }


  call UGMrdTasks.BedIntersectAndExclude as CreateFpTrainingRegionsBed {
      input:
          include_regions = fp_training_include_regions,
          exclude_regions = flatten([select_first([fp_training_exclude_regions, []]), select_first([somatic_mutations_list, []])]),
          output_basename = "FP",
          memory_gb = 16,
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
      train_set_size                  = train_set_size,
      test_set_size                   = test_set_size,
      hom_snv_regions_bed             = CreateTpTrainingRegionsBed.merged_bed,
      single_substitution_regions_bed = CreateFpTrainingRegionsBed.merged_bed,
      numerical_features              = numerical_features,
      categorical_features            = categorical_features,
      balanced_sampling_info_fields   = balanced_sampling_info_fields,
      pre_filter                      = pre_filter,
      random_seed                     = random_seed,
      balanced_strand_adapter_version = balanced_strand_adapter_version,
      references                      = references,
      flow_order                      = ExtractSampleNameFlowOrder.flow_order,
      monitoring_script               = monitoring_script,  #!FileCoercion
      docker                          = global.ug_vc_docker,
      preemptible_tries               = preemptibles,
    }


  call UGMrdTasks.InferenceSnvQualityRecalibrationModel as InferenceSnvQualityRecalibrationModel {
    input:
      output_file                             = base_file_name_sub + ".with_ml_qual.vcf.gz",
      model_file                              = TrainSnvQualityRecalibrationModel.model_file,
      params_file                             = TrainSnvQualityRecalibrationModel.params_file,
      featuremap                              = FeatureMap.featuremap,
      featuremap_index		                    = FeatureMap.featuremap_index,
      X_train_file                            = TrainSnvQualityRecalibrationModel.X_train_file,
      test_set_mrd_simulation_dataframe_file  = TrainSnvQualityRecalibrationModel.test_set_mrd_simulation_dataframe,
      monitoring_script                       = monitoring_script,  #!FileCoercion
      docker                                  = global.ug_vc_docker,
      preemptible_tries                       = preemptibles,
  }

  output {
    File featuremap = InferenceSnvQualityRecalibrationModel.featuremap_with_ml_qual
    File featuremap_index = InferenceSnvQualityRecalibrationModel.featuremap_with_ml_qual_index
    File report_html = TrainSnvQualityRecalibrationModel.test_report_file_html

    File model_file = TrainSnvQualityRecalibrationModel.model_file  
    #TODO save all the dataframes and statistics in one h5 file  
    File X_test_file = TrainSnvQualityRecalibrationModel.X_test_file
    File y_test_file = TrainSnvQualityRecalibrationModel.y_test_file
    File qual_test_file = TrainSnvQualityRecalibrationModel.qual_test_file
    File X_train_file = TrainSnvQualityRecalibrationModel.X_train_file
    File y_train_file = TrainSnvQualityRecalibrationModel.y_train_file
    File params_file = TrainSnvQualityRecalibrationModel.params_file
    File test_set_mrd_simulation_dataframe = TrainSnvQualityRecalibrationModel.test_set_mrd_simulation_dataframe
    File train_set_mrd_simulation_dataframe = TrainSnvQualityRecalibrationModel.train_set_mrd_simulation_dataframe
    File test_set_statistics_h5 = TrainSnvQualityRecalibrationModel.test_set_statistics_h5
    File train_set_statistics_h5 = TrainSnvQualityRecalibrationModel.train_set_statistics_h5
    File train_set_statistics_json = TrainSnvQualityRecalibrationModel.train_set_statistics_json
    File aggregated_metrics_json = TrainSnvQualityRecalibrationModel.test_set_statistics_json

    File tp_training_regions_bed = CreateTpTrainingRegionsBed.merged_bed
    File fp_training_regions_bed =  CreateFpTrainingRegionsBed.merged_bed

    File train_report_file_notebook= TrainSnvQualityRecalibrationModel.train_report_file_notebook
    File train_report_file_html = TrainSnvQualityRecalibrationModel.train_report_file_html
    File test_report_file_notebook =  TrainSnvQualityRecalibrationModel.test_report_file_notebook
    String flow_order = ExtractSampleNameFlowOrder.flow_order
   }
}

