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
# Deep Single Read SNV (DNN) pipeline with a `mode` selector.
# Trains a convolutional neural network on read tensors (pileup images) and/or runs inference
# to produce DNN quality-annotated featuremap VCF outputs. Modes: full (train + inference),
# train_only, inference_only (from a provided pre-trained model), and data_prep_only (tensor export).

# CHANGELOG in reverse chronological order
# 1.33.0 Merge DeepSingleReadSNVTrain via `mode` selector; add inference_only (provided model,
#        optional featuremap skip) and data_prep_only (tensor cache export) modes
# 1.31.0 Split from single_read_snv.wdl into standalone workflow

import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/single_read_snv_tasks.wdl" as SRSNVTasks
import "tasks/deep_single_read_snv_tasks.wdl" as DeepSRSNVTasks
import "tasks/globals.wdl" as Globals
import "tasks/genome_resources.wdl" as GenomeResourcesLib
import "feature_map_prep.wdl" as PrepareFeaturemap

workflow DeepSingleReadSNV {
input {
  File input_cram_bam
  File input_cram_bam_index
  Array[File]? sorter_json_stats_file_list
  String base_file_name
  String pipeline_version = "1.35.0"

  # Genome resources
  String reference_genome = "hg38"

  FeatureMapParams featuremap_params
  SingleReadSNVParams single_read_snv_params
  Array[String] features

  # annotation files
  FeaturemapAnnotationFiles annotation_files

  # model training parameters
  File? random_sample_trinuc_freq

  # DNN parameters (required)
  DeepSRSNVParams deep_srsnv_params

  # Mode selector: "full" (train + inference), "train_only", "inference_only", "data_prep_only".
  # Value comes from the input template (set to "full" there), not a WDL default.
  String mode

  # inference_only inputs: a pre-trained N-fold model, and (optionally) a pre-computed featuremap
  # VCF to skip CreateFeatureMap and go straight to vcf-to-parquet.
  DeepSRSNVModel? inference_models
  File? input_featuremap_vcf
  File? input_featuremap_vcf_index

  # Multi-VCF filtering field names
  String exclude_from_training_field_name
  String include_in_inference_field_name
  String pcawg_field_name
  String include_vcf_bcftools_filter_args

  Boolean raise_exceptions_in_report = false

  # Scatter configuration for snvfind parallelization
  Int num_shards_featuremap
  File scatter_interval_list

  Int? override_memory_gb_CreateFeatureMap
  Int? override_memory_gb_PrepareRawFeatureMap
  Int? override_memory_gb_PrepareRandomSampleFeatureMap
  Int? preemptible_tries
  Boolean? no_address_override
  String? cloud_provider_override
  File? monitoring_script_input

  Float? mean_coverage
  String? total_aligned_bases

  # winval validations
  #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
  #@wv not("test" in base_file_name or "train" in base_file_name)
  #@wv suffix(input_cram_bam) in {".bam", ".cram"}
  #@wv suffix(input_cram_bam_index) in {".bai", ".crai"}
  #@wv reference_genome in {"hg38"}
  #@wv defined(sorter_json_stats_file_list) -> suffix(sorter_json_stats_file_list) <= {".json"}
  #@wv single_read_snv_params['num_CV_folds'] > 1
  #@wv (defined(mean_coverage) or defined(total_aligned_bases)) -> (defined(mean_coverage) and defined(total_aligned_bases) and not defined(sorter_json_stats_file_list))
  #@wv defined(sorter_json_stats_file_list) -> (not defined(mean_coverage) and not defined(total_aligned_bases))
  #@wv defined(random_sample_trinuc_freq) -> suffix(random_sample_trinuc_freq) in {".csv", ".tsv"}
  #@wv mode in {"full", "train_only", "inference_only", "data_prep_only"}
  #@wv mode == "inference_only" -> defined(inference_models)
  #@wv mode != "inference_only" -> not defined(inference_models)
  #@wv defined(input_featuremap_vcf) -> (defined(input_featuremap_vcf_index) and mode == "inference_only")
}

meta {
    description : "Deep Single Read SNV (DeepSRSNV) pipeline. A convolutional neural network on read tensors (pileup images) with k-fold cross-validation, driven by a `mode` selector: 'full' (train + inference, default), 'train_only' (produce trained fold models + QC report, no inference), 'inference_only' (apply a provided pre-trained N-fold model to produce a DNN quality-annotated featuremap VCF, optionally skipping CreateFeatureMap when a featuremap VCF is supplied), and 'data_prep_only' (produce per-CRAM training tensor caches and stop, for later pooled multi-CRAM training). Requires GPU-enabled execution environment for training/inference modes."
    author: "Ultima Genomics"
    WDL_AID: { exclude: [
        "pipeline_version",
        "cloud_provider_override",
        "no_address_override",
        "preemptible_tries",
        "monitoring_script_input",
        "Globals.glob",
        "FeatureMapPrep.preemptible_tries",
        "FeatureMapPrep.monitoring_script",
        "FeatureMapPrep.featuremap_docker",
        "FeatureMapPrep.ugbio_featuremap_docker",
        "FeatureMapPrep.gatk_docker",
        "FeatureMapPrep.scatter_intervals_break",
        "FeatureMapPrep.model_files",
        "FeatureMapPrep.run_training_prep",
        "FeatureMapPrep.ConcatFeaturemapVcfs.disk_size",
        "FeatureMapPrep.ConcatRandomSampleVcfs.disk_size",
        "DNNCramToTensorsPos.memory_gb",
        "DNNCramToTensorsPos.cpus",
        "DNNCramToTensorsNeg.memory_gb",
        "DNNCramToTensorsNeg.cpus",
        "DNNCombineSplits.memory_gb",
        "DNNCombineSplits.cpus",
        "DNNTrainFold.cpus",
        "DNNTrainFold.override_memory_gb",
        "DNNRecalibrateFolds.memory_gb",
        "DNNRecalibrateFolds.cpus",
        "DNNRecalibrateFolds.featuremap_parquets",
        "DNNVcfToParquet.memory_gb",
        "DNNVcfToParquet.cpus",
        "DNNCramToTensorsInference.cpus",
        "DNNFoldInference.memory_gb",
        "DNNFoldInference.cpus",
        "DNNMergeAndAnnotate.memory_gb",
        "DNNMergeAndAnnotate.cpus",
        "DNNPrepareReport.memory_gb",
        "DNNPrepareReport.cpus",
        "DNNReport.memory_gb",
        "DNNReport.cpus"
    ]}
}

parameter_meta {
    base_file_name: {
        help: "Base file name for output files",
        type: "String",
        category: "input_required"
    }
    input_cram_bam: {
        help: "Input CRAM file",
        type: "File",
        category: "input_required"
    }
    input_cram_bam_index: {
        help: "Input CRAM index file",
        type: "File",
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
        help: "(Optional) Total aligned bases used for downsampling rate calculation.",
        type: "String",
        category: "input_optional"
    }
    reference_genome: {
        type: "String",
        help: "Genome type selector. The workflow currently supports only hg38.",
        category: "input_optional"
    }
    featuremap_params: {
        type: "FeatureMapParams",
        help: "FeatureMap parameters, recommended value set in the template.",
        category: "param_required"
    }
    single_read_snv_params: {
        type: "SingleReadSNVParams",
        help: "SingleReadSNV parameters for training set preparation.",
        category: "param_required"
    }
    features: {
        type: "Array[String]",
        help: "Features to be used for the report quality plots, should match the XGBoost feature set.",
        category: "param_required"
    }
    annotation_files: {
        type: "FeaturemapAnnotationFiles",
        help: "Annotation files for featuremap generation: dbSNP, gnomAD, and UG High Confidence Regions with their indices",
        category: "ref_required"
    }
    deep_srsnv_params: {
        type: "DeepSRSNVParams",
        help: "Deep SRSNV (DNN) parameters. Controls training hyperparameters, fold count, GPU resources, and inference backend.",
        category: "param_required"
    }
    mode: {
        type: "String",
        help: "Pipeline mode (set in the input template): 'full' (train + inference), 'train_only' (training + QC report), 'inference_only' (apply provided model, requires inference_models), or 'data_prep_only' (produce training tensor caches and stop).",
        category: "input_required"
    }
    inference_models: {
        type: "DeepSRSNVModel",
        help: "Pre-trained N-fold model for mode='inference_only'. fold_metadata (recalibrated), fold_checkpoints, and fold_onnx_models are REQUIRED (engine rebuilt from ONNX in-runtime); fold_engines optional. Fold count inferred from the arrays.",
        category: "input_optional"
    }
    input_featuremap_vcf: {
        type: "File",
        help: "(inference_only) Pre-computed FeatureMap VCF. If provided, CreateFeatureMap is skipped and inference runs directly on it.",
        category: "input_optional"
    }
    input_featuremap_vcf_index: {
        type: "File",
        help: "(inference_only) Index for input_featuremap_vcf. Required if input_featuremap_vcf is provided.",
        category: "input_optional"
    }
    random_sample_trinuc_freq: {
        type: "File",
        help: "(Optional) CSV or TSV file with trinucleotide frequencies for the random sample.",
        category: "input_optional"
    }
    exclude_from_training_field_name: {
        type: "String",
        help: "INFO field name for the exclude-from-training annotation in the featuremap VCF.",
        category: "param_required"
    }
    include_in_inference_field_name: {
        type: "String",
        help: "INFO field name for the include-in-inference annotation in the featuremap VCF.",
        category: "param_required"
    }
    pcawg_field_name: {
        type: "String",
        help: "INFO field name for the PCAWG annotation in the featuremap VCF.",
        category: "param_required"
    }
    include_vcf_bcftools_filter_args: {
        type: "String",
        help: "Bcftools filter arguments applied to include-in-inference VCFs before annotation.",
        category: "param_required"
    }
    num_shards_featuremap: {
        type: "Int",
        help: "Number of genomic shards to scatter the snvfind (CreateFeatureMap) step across. Higher values reduce wall-clock time but add scatter overhead.",
        category: "input_required"
    }
    scatter_interval_list: {
        type: "File",
        help: "Interval list defining the genomic regions to scatter snvfind across. Should match the regions in featuremap_params.bed_file.",
        category: "input_required"
    }
    override_memory_gb_CreateFeatureMap: {
        type: "Int",
        help: "Override memory in GB for the CreateFeatureMap task. If an out of memory error occurs, try increasing this value.",
        category: "optional"
    }
    override_memory_gb_PrepareRawFeatureMap: {
        type: "Int",
        help: "Override memory in GB for the PrepareRawFeatureMap task. If an out of memory error occurs, try increasing this value.",
        category: "optional"
    }
    override_memory_gb_PrepareRandomSampleFeatureMap: {
        type: "Int",
        help: "Override memory in GB for the PrepareRandomSampleFeatureMap task. If an out of memory error occurs, try increasing this value.",
        category: "optional"
    }
    raise_exceptions_in_report: {
        type: "Boolean",
        help: "Raise an exception and fail the pipeline if an error is raised in the QC report",
        category: "optional"
    }
    featuremap: {
        type: "File?",
        help: "FeatureMap VCF with all SNV candidates",
        category: "output"
    }
    featuremap_index: {
        type: "File?",
        help: "Index for the FeatureMap VCF",
        category: "output"
    }
    featuremap_random_sample: {
        type: "File?",
        help: "Downsampled FeatureMap VCF file for training",
        category: "output"
    }
    featuremap_random_sample_index: {
        type: "File?",
        help: "Downsampled FeatureMap VCF index file",
        category: "output"
    }
    downsampling_rate: {
        type: "Float?",
        help: "The downsampling rate used to create the random sample featuremap",
        category: "output"
    }
    random_sample_trinuc_freq_stats: {
        type: "File?",
        help: "Trinucleotide frequency statistics from the random sample featuremap",
        category: "output"
    }
    fold_metadata: {
        type: "Array[File]?",
        help: "Per-fold metadata JSON files from DNN training",
        category: "output"
    }
    fold_checkpoints: {
        type: "Array[File]?",
        help: "Per-fold model checkpoint files from DNN training",
        category: "output"
    }
    fold_onnx_models: {
        type: "Array[File]?",
        help: "Per-fold ONNX model files from DNN training",
        category: "output"
    }
    fold_engines: {
        type: "Array[File]?",
        help: "Per-fold TensorRT engine files from DNN training",
        category: "output"
    }
    fold_trt_timing_caches: {
        type: "Array[File]?",
        help: "Per-fold TensorRT timing caches from DNN training; feed back as DeepSRSNVModel.fold_timing_caches in a later inference_only run to rebuild bit-identical engines (same GPU + TRT version)",
        category: "output"
    }
    report_html: {
        type: "File?",
        help: "QC report HTML file",
        category: "output"
    }
    application_qc_h5: {
        type: "File?",
        help: "Application QC statistics h5 file",
        category: "output"
    }
    combined_metadata: {
        type: "File?",
        help: "Shared quality recalibration LUT metadata",
        category: "output"
    }
    featuremap_df: {
        type: "File?",
        help: "Combined featuremap DataFrame (parquet) with per-fold predictions",
        category: "output"
    }
    updated_fold_metadata: {
        type: "Array[File]?",
        help: "Per-fold metadata JSONs updated with the shared recalibration LUT (training modes only)",
        category: "output"
    }
    featuremap_vcf: {
        type: "File?",
        help: "FeatureMap VCF annotated with DNN quality scores (inference modes)",
        category: "output"
    }
    featuremap_vcf_index: {
        type: "File?",
        help: "Index for the annotated FeatureMap VCF (inference modes)",
        category: "output"
    }
    positive_parquet: {
        type: "File?",
        help: "Positive-label (random-sample) training featuremap parquet (training/data_prep modes)",
        category: "output"
    }
    negative_parquet: {
        type: "File?",
        help: "Negative-label (raw) training featuremap parquet (training/data_prep modes)",
        category: "output"
    }
    positive_tensor_cache_tar: {
        type: "File?",
        help: "Positive-label read tensor cache tar for pooled multi-CRAM training (data_prep_only mode only)",
        category: "output"
    }
    negative_tensor_cache_tar: {
        type: "File?",
        help: "Negative-label read tensor cache tar for pooled multi-CRAM training (data_prep_only mode only)",
        category: "output"
    }
    stats_funnel: {
        type: "File?",
        help: "snvfind model-filters status funnel JSON, needed as --stats-file for pooled multi-CRAM training (data_prep_only mode only)",
        category: "output"
    }
    mean_coverage_out: {
        type: "Float?",
        help: "Resolved mean coverage, needed as --mean-coverage for pooled multi-CRAM training (data_prep_only mode only)",
        category: "output"
    }
}

  Int preemptibles = select_first([preemptible_tries, 1])
  String base_file_name_sub = sub(base_file_name, "#", "")
  Boolean no_address = select_first([no_address_override, true])

  # Mode-derived control flags
  Boolean do_training  = mode == "full" || mode == "train_only"        # combine+split -> train -> recalibrate -> report
  Boolean do_data_prep = mode == "data_prep_only"                      # tensorize and stop (export raw tensor tars)
  Boolean do_tensorize = do_training || do_data_prep                   # DNNCramToTensorsPos/Neg (shared stage)
  Boolean do_inference = mode == "full" || mode == "inference_only"    # vcf->parquet -> fold inference -> annotate
  Boolean need_training_prep = do_training || do_data_prep             # positive/negative parquet generation in FeatureMapPrep
  # Skip CreateFeatureMap only when a featuremap VCF is provided in inference_only (bypasses FeatureMapPrep entirely).
  Boolean skip_featuremap_prep = defined(input_featuremap_vcf) && mode == "inference_only"
  Boolean run_featuremap_prep  = !skip_featuremap_prep

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

  call GenomeResourcesLib.GenomeResourcesWorkflow as GenomeResources

  References references = object {
    ref_fasta: GenomeResources.resources[reference_genome].ref_fasta,
    ref_fasta_index: GenomeResources.resources[reference_genome].ref_fasta_index,
    ref_dict: GenomeResources.resources[reference_genome].ref_dict
  }

  File training_interval_list = GenomeResources.resources[reference_genome].srsnv_training_interval_list

  # Coverage/sorter-stats resolution + FeatureMap preparation. Nested under run_featuremap_prep:
  # in inference_only-with-featuremap there is no coverage input, and select_first on an all-None
  # list is a runtime error, so these declarations must not be evaluated in that path.
  if (run_featuremap_prep) {
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

    # FeatureMap preparation via sub-workflow. run_training_prep gates positive/negative parquet
    # generation: needed for training and data_prep, not for inference_only.
    call PrepareFeaturemap.FeatureMapPrep {
      input:
        input_cram_bam_list = [input_cram_bam],
        input_cram_bam_index_list = [input_cram_bam_index],
        base_file_name = base_file_name_sub,
        references = references,
        training_interval_list = training_interval_list,
        scatter_interval_list = scatter_interval_list,
        featuremap_params = featuremap_params,
        single_read_snv_params = single_read_snv_params,
        annotation_files = annotation_files,
        mean_coverage = mean_coverage_used,
        total_aligned_bases = total_aligned_bases_used,
        random_sample_trinuc_freq = random_sample_trinuc_freq,
        exclude_from_training_field_name = exclude_from_training_field_name,
        include_in_inference_field_name = include_in_inference_field_name,
        pcawg_field_name = pcawg_field_name,
        include_vcf_bcftools_filter_args = include_vcf_bcftools_filter_args,
        override_memory_gb_CreateFeatureMap = override_memory_gb_CreateFeatureMap,
        override_memory_gb_PrepareRawFeatureMap = override_memory_gb_PrepareRawFeatureMap,
        override_memory_gb_PrepareRandomSampleFeatureMap = override_memory_gb_PrepareRandomSampleFeatureMap,
        run_training_prep = need_training_prep,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script,
        num_shards = num_shards_featuremap,
        featuremap_docker = global.featuremap_docker,
        ugbio_featuremap_docker = global.ugbio_featuremap_docker,
        gatk_docker = global.broad_gatk_docker
    }
  }

  # Featuremap used for inference: provided VCF (inference_only) or the freshly-created one.
  File? featuremap_for_inference       = if defined(input_featuremap_vcf) then input_featuremap_vcf else FeatureMapPrep.featuremap
  File? featuremap_for_inference_index = if defined(input_featuremap_vcf_index) then input_featuremap_vcf_index else FeatureMapPrep.featuremap_index

  # ============================================================
  # Tensorize stage (shared by full, train_only, data_prep_only). data_prep_only stops here.
  # do_tensorize => need_training_prep => FeatureMapPrep ran with parquets, so select_first is safe.
  # ============================================================
  if (do_tensorize) {
    call DeepSRSNVTasks.DNNCramToTensors as DNNCramToTensorsPos {
      input:
        input_cram = input_cram_bam,
        input_cram_index = input_cram_bam_index,
        featuremap_parquet = select_first([FeatureMapPrep.positive_parquet]),
        label = "positive",
        references = references,
        deep_srsnv_params = deep_srsnv_params,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }
    call DeepSRSNVTasks.DNNCramToTensors as DNNCramToTensorsNeg {
      input:
        input_cram = input_cram_bam,
        input_cram_index = input_cram_bam_index,
        featuremap_parquet = select_first([FeatureMapPrep.negative_parquet]),
        label = "negative",
        references = references,
        deep_srsnv_params = deep_srsnv_params,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }
  }

  # Expose the raw tensor tars ONLY in data_prep_only. In full/train_only they are consumed
  # internally by DNNCombineSplits and must not become workflow outputs.
  if (do_data_prep) {
    File? data_prep_positive_tensor_cache_tar = DNNCramToTensorsPos.tensor_cache_tar
    File? data_prep_negative_tensor_cache_tar = DNNCramToTensorsNeg.tensor_cache_tar
    # Also surface the snvfind stats funnel and resolved mean coverage — a later pooled
    # multi-CRAM training run needs these (DNNTrainFold/DNNRecalibrateFolds --stats-file /
    # --mean-coverage) in addition to the tensor caches and parquets.
    File? data_prep_stats_funnel = select_first([FeatureMapPrep.model_filters_status_funnel])
    Float? data_prep_mean_coverage = mean_coverage_used
  }

  # ============================================================
  # Training stage (full, train_only) — combine+split, train folds, recalibrate, QC report.
  # ============================================================
  if (do_training) {
    # 2. Combine + k-fold split
    call DeepSRSNVTasks.DNNCombineSplits {
      input:
        positive_tensor_cache_tar = select_first([DNNCramToTensorsPos.tensor_cache_tar]),
        negative_tensor_cache_tar = select_first([DNNCramToTensorsNeg.tensor_cache_tar]),
        training_interval_list = training_interval_list,
        deep_srsnv_params = deep_srsnv_params,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }

    # 3. Train each fold in parallel (scatter)
    scatter (fold_idx in range(deep_srsnv_params.num_folds)) {
      call DeepSRSNVTasks.DNNTrainFold {
        input:
          fold_tar = DNNCombineSplits.fold_tars[fold_idx],
          fold_idx = fold_idx,
          deep_srsnv_params = deep_srsnv_params,
          pretrained_checkpoint = deep_srsnv_params.pretrained_checkpoint,
          stats_file = select_first([FeatureMapPrep.model_filters_status_funnel]),
          mean_coverage = select_first([mean_coverage_used]),
          training_interval_list = training_interval_list,
          base_file_name = base_file_name_sub,
          docker = global.ugbio_deep_srsnv_docker,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script
      }
    }

    # 4. Build shared quality recalibration LUT
    call DeepSRSNVTasks.DNNRecalibrateFolds {
      input:
        fold_parquets = DNNTrainFold.dnn_featuremap_df,
        fold_metadata = DNNTrainFold.dnn_metadata,
        featuremap_parquets = [select_first([FeatureMapPrep.positive_parquet]), select_first([FeatureMapPrep.negative_parquet])],
        stats_file = select_first([FeatureMapPrep.model_filters_status_funnel]),
        training_interval_list = training_interval_list,
        mean_coverage = select_first([mean_coverage_used]),
        base_file_name = base_file_name_sub,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script,
    }

    # Report — training modes only (consumes training parquets/metadata; skipped in inference_only/data_prep_only)
    call DeepSRSNVTasks.DNNPrepareReport {
      input:
        positive_featuremap_df = select_first([FeatureMapPrep.positive_parquet]),
        negative_featuremap_df = select_first([FeatureMapPrep.negative_parquet]),
        dnn_combined_featuremap_df = DNNRecalibrateFolds.combined_featuremap_df,
        training_metadata = DNNRecalibrateFolds.updated_fold_metadata[0],
        dnn_fold_0_metadata = DNNTrainFold.dnn_metadata[0],
        dnn_fold_metadata = DNNRecalibrateFolds.updated_fold_metadata,
        features = features,
        base_file_name = base_file_name_sub + "_dnn",
        pipeline_version = pipeline_version,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }

    call DeepSRSNVTasks.DNNReport {
      input:
        report_featuremap_df = DNNPrepareReport.report_featuremap_df,
        report_metadata = DNNPrepareReport.report_metadata,
        base_file_name = base_file_name_sub + "_dnn",
        docker = global.ugbio_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }
  }

  # Fold count for inference. Computed from params (full) or the provided model (inference_only) —
  # NOT from any training output — so it never gates the inference data-prep on training. Kept inside
  # a do_inference guard so select_first([inference_models]) is not evaluated in data_prep_only.
  # ============================================================
  # Inference data-prep (full, inference_only): featuremap VCF -> parquet -> per-fold tensors.
  # Depends only on the featuremap VCF + CRAM (available before/independently of training), so in
  # full mode this CPU work runs in PARALLEL with training and data-prep. It is deliberately in a
  # separate conditional scope from the model-resolution block below, whose declarations reference
  # training outputs — keeping them together would make the engine gate this prep on training.
  # ============================================================
  if (do_inference) {
    Int num_folds_prep = if do_training then deep_srsnv_params.num_folds else length(select_first([inference_models]).fold_metadata)

    # Convert featuremap VCF to parquet
    call DeepSRSNVTasks.DNNVcfToParquet {
      input:
        featuremap_vcf = select_first([featuremap_for_inference]),
        featuremap_vcf_index = select_first([featuremap_for_inference_index]),
        inference_filters = FeatureMapPrep.inference_filters,
        base_file_name = base_file_name_sub,
        docker = global.ugbio_featuremap_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }

    # Pre-compute inference tensors per fold (CPU-only)
    scatter (tensor_fold_idx in range(num_folds_prep)) {
      call DeepSRSNVTasks.DNNCramToTensorsInference {
        input:
          input_cram = input_cram_bam,
          input_cram_index = input_cram_bam_index,
          featuremap_parquet = DNNVcfToParquet.featuremap_parquet,
          training_interval_list = training_interval_list,
          references = references,
          deep_srsnv_params = deep_srsnv_params,
          fold_idx = tensor_fold_idx,
          num_folds = num_folds_prep,
          base_file_name = base_file_name_sub,
          docker = global.ugbio_deep_srsnv_docker,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script
      }
    }
  }

  # ============================================================
  # Model resolution + per-fold GPU inference + merge/annotate. The model arrays come from the
  # training outputs (full) or the provided inference_models (inference_only). DNNFoldInference
  # genuinely needs the trained model, so this scope legitimately waits for training in full mode;
  # it consumes the tensors pre-computed above (optional across the conditional-scope boundary).
  # ============================================================
  if (do_inference) {
    # if/then/else branches are lazily evaluated, so select_first([inference_models]) is never
    # touched in full mode and the training outputs are never touched in inference_only mode.
    Array[File]  fold_checkpoints_used = if do_training then select_first([DNNTrainFold.dnn_checkpoint]) else select_first([inference_models]).fold_checkpoints
    Array[File]  fold_metadata_used    = if do_training then select_first([DNNRecalibrateFolds.updated_fold_metadata]) else select_first([inference_models]).fold_metadata
    Array[File]? fold_onnx_used   = if do_training then DNNTrainFold.dnn_onnx   else select_first([inference_models]).fold_onnx_models
    Array[File]? fold_engine_used = if do_training then DNNTrainFold.dnn_engine else select_first([inference_models]).fold_engines
    Array[File]? fold_timing_cache_used = if do_training then DNNTrainFold.dnn_trt_timing_cache else select_first([inference_models]).fold_timing_caches
    Int num_folds_for_inference = if do_training then deep_srsnv_params.num_folds else length(fold_metadata_used)
    # Tensors from the data-prep block above become optional outside their conditional scope.
    Array[Array[File]] inference_tensor_shards = select_first([DNNCramToTensorsInference.tensor_shards])

    # Per-fold GPU inference. onnx/engine fall back to the metadata file when the (optional)
    # model arrays are absent (non-trt backend), mirroring the inference-only reference workflow.
    scatter (infer_fold_idx in range(num_folds_for_inference)) {
      # Optional per-fold timing cache (File?): present only when the model carries one. Declared as a
      # conditional so it stays truly optional (no metadata-JSON fallback, which is not a timing cache).
      if (defined(fold_timing_cache_used)) {
        File fold_timing_cache_for_fold = select_first([fold_timing_cache_used])[infer_fold_idx]
      }
      call DeepSRSNVTasks.DNNFoldInference {
        input:
          tensor_shards = inference_tensor_shards[infer_fold_idx],
          fold_metadata = fold_metadata_used[infer_fold_idx],
          fold_checkpoint = fold_checkpoints_used[infer_fold_idx],
          # inference_only REQUIRES ONNX (the engine is rebuilt from it); select_first fails loudly
          # if fold_onnx_models was not supplied. full mode always has onnx from DNNTrainFold.
          fold_onnx_model = if !do_training then select_first([fold_onnx_used])[infer_fold_idx]
                            else (if defined(fold_onnx_used) then select_first([fold_onnx_used])[infer_fold_idx] else fold_metadata_used[infer_fold_idx]),
          # engine: rebuilt in-task for inference_only, so the provided/placeholder value is
          # overwritten; keep the metadata fallback so a missing engine array is harmless there.
          fold_engine = if defined(fold_engine_used) then select_first([fold_engine_used])[infer_fold_idx] else fold_metadata_used[infer_fold_idx],
          # optional timing cache: passed to the engine rebuild (inference_only) for a bit-identical
          # engine. Absent -> the rebuild times tactics fresh. full mode doesn't rebuild, so it's unused.
          fold_timing_cache = fold_timing_cache_for_fold,
          fold_idx = infer_fold_idx,
          deep_srsnv_params = deep_srsnv_params,
          # Rebuild the engine from ONNX only in inference_only AND when the trt backend is used
          # (the pytorch backend loads the .ckpt, not the engine, so a rebuild is pointless/wasteful
          # and could fail). Within do_inference, !do_training is exactly mode=="inference_only";
          # full keeps its in-env engine. Backend defaults to "trt".
          rebuild_engine_from_onnx = !do_training && select_first([deep_srsnv_params.inference_backend, "trt"]) == "trt",
          base_file_name = base_file_name_sub,
          docker = global.ugbio_deep_srsnv_docker,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script
      }
    }

    # Merge fold predictions + annotate VCF
    call DeepSRSNVTasks.DNNMergeAndAnnotate {
      input:
        featuremap_vcf = select_first([featuremap_for_inference]),
        featuremap_vcf_index = select_first([featuremap_for_inference_index]),
        fold_predictions = DNNFoldInference.fold_predictions_parquet,
        fold_metadata_0 = fold_metadata_used[0],
        deep_srsnv_params = deep_srsnv_params,
        base_file_name = base_file_name_sub,
        docker = global.ugbio_deep_srsnv_docker,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script
    }
  }

  output {
    # FeatureMap outputs (present whenever FeatureMapPrep ran)
    File? featuremap = FeatureMapPrep.featuremap
    File? featuremap_index = FeatureMapPrep.featuremap_index
    File? featuremap_random_sample = FeatureMapPrep.featuremap_random_sample
    File? featuremap_random_sample_index = FeatureMapPrep.featuremap_random_sample_index
    File? random_sample_trinuc_freq_stats = FeatureMapPrep.random_sample_trinuc_freq_stats
    Float? downsampling_rate = FeatureMapPrep.downsampling_rate

    # Training-set parquets (training/data_prep modes)
    File? positive_parquet = FeatureMapPrep.positive_parquet
    File? negative_parquet = FeatureMapPrep.negative_parquet

    # Raw tensor caches + training inputs for pooled multi-CRAM training (data_prep_only mode only)
    File? positive_tensor_cache_tar = data_prep_positive_tensor_cache_tar
    File? negative_tensor_cache_tar = data_prep_negative_tensor_cache_tar
    File? stats_funnel = data_prep_stats_funnel
    Float? mean_coverage_out = data_prep_mean_coverage

    # Annotated VCF (inference modes)
    File? featuremap_vcf = DNNMergeAndAnnotate.dnn_featuremap_vcf
    File? featuremap_vcf_index = DNNMergeAndAnnotate.dnn_featuremap_vcf_index

    # Per-fold model files + recalibration (training modes)
    Array[File]? fold_metadata = DNNTrainFold.dnn_metadata
    Array[File]? fold_checkpoints = DNNTrainFold.dnn_checkpoint
    Array[File]? fold_onnx_models = DNNTrainFold.dnn_onnx
    Array[File]? fold_engines = DNNTrainFold.dnn_engine
    # Per-fold TensorRT timing caches: feed these back as DeepSRSNVModel.fold_timing_caches in a later
    # inference_only run to rebuild bit-identical engines (same GPU + TRT version).
    Array[File]? fold_trt_timing_caches = DNNTrainFold.dnn_trt_timing_cache
    File? combined_metadata = DNNRecalibrateFolds.shared_lut_metadata
    Array[File]? updated_fold_metadata = DNNRecalibrateFolds.updated_fold_metadata
    File? featuremap_df = DNNRecalibrateFolds.combined_featuremap_df

    # QC report (training modes)
    File? report_html = DNNReport.dnn_report_html
    File? application_qc_h5 = DNNReport.dnn_application_qc_h5
  }
}
