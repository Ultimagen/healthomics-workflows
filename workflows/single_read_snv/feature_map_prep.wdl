version 1.0

# LICENSE
#   Copyright 2024 Ultima Genomics
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

# DESCRIPTION
# Sub-workflow for FeatureMap preparation: creates the annotated featuremap VCF and
# prepares positive/negative parquet training sets for both XGBoost and DNN pipelines.

import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/single_read_snv_tasks.wdl" as SRSNVTasks

workflow FeatureMapPrep {
input {
  Array[File] input_cram_bam_list
  Array[File] input_cram_bam_index_list
  String base_file_name

  # References
  References references
  File training_interval_list

  # Params
  FeatureMapParams featuremap_params
  SingleReadSNVParams single_read_snv_params

  # annotation files
  FeaturemapAnnotationFiles annotation_files

  # Coverage info
  Float mean_coverage
  String total_aligned_bases

  # Training params
  File? random_sample_trinuc_freq

  # Multi-VCF filtering field names
  String exclude_from_training_field_name
  String include_in_inference_field_name
  String pcawg_field_name
  String include_vcf_bcftools_filter_args

  # Memory overrides
  Int? override_memory_gb_CreateFeatureMap
  Int? override_memory_gb_PrepareRawFeatureMap
  Int? override_memory_gb_PrepareRandomSampleFeatureMap

  # Execution params
  Int preemptible_tries = 1
  File monitoring_script

  # Docker images
  String featuremap_docker
  String ugbio_featuremap_docker
}

  String base_file_name_sub = sub(base_file_name, "#", "")

  # Calculate random_sample_size based on tp_train_set_size and tp_train_set_size_sampling_overhead
  Int random_sample_size = ceil(single_read_snv_params.tp_train_set_size * single_read_snv_params.tp_train_set_size_sampling_overhead)

  # Coverage threshold for read filters
  Float max_coverage_factor = single_read_snv_params.max_coverage_factor
  if (defined(featuremap_params.read_filters)) {
    Int coverage_threshold = ceil(mean_coverage * max_coverage_factor)
  }

  # Multi-VCF annotation (all optional)
  Array[File] exclude_vcfs_raw = select_first([annotation_files.exclude_from_training_vcf_list, []])
  Array[File] exclude_vcf_idxs_raw = select_first([annotation_files.exclude_from_training_vcf_index_list, []])
  Array[File] include_vcfs_raw = select_first([annotation_files.include_in_inference_vcf_list, []])
  Array[File] include_vcf_idxs_raw = select_first([annotation_files.include_in_inference_vcf_index_list, []])

  # Filter include VCFs to PASS-only biallelic SNPs
  if (length(include_vcfs_raw) > 0) {
    scatter (idx in range(length(include_vcfs_raw))) {
      call UGGeneralTasks.FilterVcfWithBcftools as FilterIncludeVcf {
        input:
          input_vcf = include_vcfs_raw[idx],
          bcftools_extra_args = include_vcf_bcftools_filter_args,
          docker = featuremap_docker,
          monitoring_script = monitoring_script,
          preemptible_tries = preemptible_tries
      }
    }
  }
  Array[File] include_vcfs_filtered = select_first([FilterIncludeVcf.output_vcf, []])
  Array[File] include_vcf_idxs_filtered = select_first([FilterIncludeVcf.output_vcf_index, []])

  Boolean has_annotation_vcfs = (length(exclude_vcfs_raw) > 0 || length(include_vcfs_filtered) > 0 || defined(annotation_files.pcawg_vcf))

  # Prepare read_filters JSON with annotation exclusion filters + inference_filters.json
  if (has_annotation_vcfs) {
    call SRSNVTasks.PrepareAnnotationVcfs {
      input:
        exclude_field_name = exclude_from_training_field_name,
        include_field_name = include_in_inference_field_name,
        pcawg_field_name = pcawg_field_name,
        has_exclude = length(exclude_vcfs_raw) > 0,
        has_include = length(include_vcfs_filtered) > 0,
        has_pcawg = defined(annotation_files.pcawg_vcf),
        read_filters_json = featuremap_params.read_filters,
        coverage_threshold = coverage_threshold,
        docker = ugbio_featuremap_docker,
        preemptible_tries = preemptible_tries,
        monitoring_script = monitoring_script
    }
  }

  # Run snvfind to create featuremap
  call SRSNVTasks.CreateFeatureMap {
    input:
      input_cram_bam_list = input_cram_bam_list,
      input_cram_bam_index_list = input_cram_bam_index_list,
      references = references,
      base_file_name = base_file_name_sub,
      total_aligned_bases = total_aligned_bases,
      mean_coverage = mean_coverage,
      max_coverage_factor = max_coverage_factor,
      random_sample_size = random_sample_size,
      random_sample_trinuc_freq_ = random_sample_trinuc_freq,
      featuremap_params = featuremap_params,
      annotation_files = annotation_files,
      exclude_annotation_vcfs = exclude_vcfs_raw,
      exclude_annotation_vcf_indices = exclude_vcf_idxs_raw,
      exclude_annotation_field_name = exclude_from_training_field_name,
      include_annotation_vcfs = include_vcfs_filtered,
      include_annotation_vcf_indices = include_vcf_idxs_filtered,
      include_annotation_field_name = include_in_inference_field_name,
      pcawg_annotation_vcf = annotation_files.pcawg_vcf,
      pcawg_annotation_vcf_index = annotation_files.pcawg_vcf_index,
      pcawg_annotation_field_name = pcawg_field_name,
      augmented_read_filters = PrepareAnnotationVcfs.augmented_read_filters,
      docker = featuremap_docker,
      preemptible_tries = preemptible_tries,
      monitoring_script = monitoring_script,
      memory_gb = select_first([override_memory_gb_CreateFeatureMap, 4]),
      cpus = 2
  }

  # Prepare the raw featuremap for training
  # Dynamically assign memory for PrepareRawFeatureMap by mean_coverage
  Int memory_gb_PrepareRawFeatureMap_default = if (mean_coverage < 50.0) then 64 else 128
  Int memory_gb_PrepareRawFeatureMap = select_first([override_memory_gb_PrepareRawFeatureMap, memory_gb_PrepareRawFeatureMap_default])
  Int cpu_PrepareRawFeatureMap_ = ceil(memory_gb_PrepareRawFeatureMap / 2)
  Int cpu_PrepareRawFeatureMap = if (cpu_PrepareRawFeatureMap_ > 32) then 32 else cpu_PrepareRawFeatureMap_
  # Prefer augmented JSON (with annotation is_null filters + coverage threshold) from PrepareAnnotationVcfs
  File read_filters_with_max_coverage_ = select_first([PrepareAnnotationVcfs.augmented_read_filters, CreateFeatureMap.read_filters_with_max_coverage])

  # Generate FP training set
  call SRSNVTasks.PrepareFeatureMapForTraining as PrepareRawFeatureMap {
      input:
          featuremap = CreateFeatureMap.featuremap,
          featuremap_index = CreateFeatureMap.featuremap_index,
          single_read_snv_params = single_read_snv_params,
          train_set_size = single_read_snv_params.fp_train_set_size,
          filter_json_key = "filters_full_output",
          read_filters_with_max_coverage = read_filters_with_max_coverage_,
          docker = ugbio_featuremap_docker,
          preemptible_tries = preemptible_tries,
          memory_gb = memory_gb_PrepareRawFeatureMap,
          cpus      = cpu_PrepareRawFeatureMap,
          monitoring_script = monitoring_script
  }

  # Prepare the random sample featuremap for training
  Int memory_gb_PrepareRandomSampleFeatureMap = select_first([override_memory_gb_PrepareRandomSampleFeatureMap, 8])
  Int cpu_PrepareRandomSampleFeatureMap_ = ceil(memory_gb_PrepareRandomSampleFeatureMap / 2)
  Int cpu_PrepareRandomSampleFeatureMap = if (cpu_PrepareRandomSampleFeatureMap_ > 4) then 4 else cpu_PrepareRandomSampleFeatureMap_

  # Generate TP training set
  call SRSNVTasks.PrepareFeatureMapForTraining as PrepareRandomSampleFeatureMap {
      input:
          featuremap = select_first([CreateFeatureMap.featuremap_random_sample]),
          featuremap_index = select_first([CreateFeatureMap.featuremap_random_sample_index]),
          single_read_snv_params = single_read_snv_params,
          train_set_size = single_read_snv_params.tp_train_set_size,
          filter_json_key = "filters_random_sample",
          read_filters_with_max_coverage = read_filters_with_max_coverage_,
          docker = ugbio_featuremap_docker,
          preemptible_tries = preemptible_tries,
          memory_gb = memory_gb_PrepareRandomSampleFeatureMap,
          cpus = cpu_PrepareRandomSampleFeatureMap,
          monitoring_script = monitoring_script
  }

  output {
    File featuremap = CreateFeatureMap.featuremap
    File featuremap_index = CreateFeatureMap.featuremap_index
    File? featuremap_random_sample = CreateFeatureMap.featuremap_random_sample
    File? featuremap_random_sample_index = CreateFeatureMap.featuremap_random_sample_index
    Float downsampling_rate = CreateFeatureMap.downsampling_rate
    File? random_sample_trinuc_freq_stats = CreateFeatureMap.random_sample_trinuc_freq
    File model_filters_status_funnel = CreateFeatureMap.model_filters_status_funnel
    File positive_parquet = PrepareRandomSampleFeatureMap.filtered_featuremap_parquet
    File negative_parquet = PrepareRawFeatureMap.filtered_featuremap_parquet
    File? inference_filters = PrepareAnnotationVcfs.inference_filters
    File? augmented_read_filters = PrepareAnnotationVcfs.augmented_read_filters
  }
}
