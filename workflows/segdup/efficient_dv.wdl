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

# DESCRIPTION
# Runs efficientDV variant calling pipeline for Ultima Genomics data.
#
# CHANGELOG
# 1.11.0 Bug fixes in somatic efficient DV, added post-processing tasks
# 1.10.3 Externalized hard cutoff on the quality
# 1.9.0 Initial release

import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/efficient_dv_tasks.wdl" as UGDVTasks
import "tasks/globals.wdl" as Globals
import "tasks/single_sample_vc_tasks.wdl" as VCTasks
import "tasks/vcf_postprocessing_tasks.wdl" as PostProcesTasks

workflow EfficientDV {
  input {
    # Workflow args
    String pipeline_version = "1.24.2" # !UnusedDeclaration
    String base_file_name

    # Mandatory inputs
    Array[File] cram_files
    Array[File] cram_index_files
    References references

    Boolean make_gvcf
    Boolean recalibrate_vaf
    Boolean is_somatic  # Enable somatic calling mode
    Boolean show_bg_fields = is_somatic # Show background fields in the output vcf

    # Scatter interval list args
    Int num_shards = 40
    Int scatter_intervals_break = 10000000 # Maximal resolution for scattering intervals
    File? target_intervals
    String? intervals_string

    # Make examples args
    Float min_fraction_hmer_indels = 0.12
    Float min_fraction_non_hmer_indels = 0.06
    Float min_fraction_snps = 0.12
    Float min_fraction_single_strand_non_snps = 0.15
    Int min_read_count_snps = 2
    Int min_read_count_hmer_indels = 2
    Int min_read_count_non_hmer_indels = 2
    Int min_base_quality = 5
    Int pileup_min_mapping_quality = 5
    Int candidate_min_mapping_quality = 5
    Int max_reads_per_partition = 1500
    Int dbg_min_base_quality = 0 # Minimal base quality during the assembly process
    Boolean prioritize_alt_supporting_reads = false
    Float p_error = 0.005
    Array[Int] optimal_coverages = [ 50 ]
    Boolean cap_at_optimal_coverage = false
    Boolean output_realignment = false
    Boolean single_strand_filter = false
    Boolean keep_duplicates = true
    Boolean add_ins_size_channel = true
    String? ug_make_examples_extra_args
    Boolean log_make_examples_progress = false
    Boolean normalize_strand_bias
    Array[Float]? strand_bias_normalization_thresholds

    File? germline_vcf

    # Background files (for somatic calling)
    Array[File] background_cram_files = []
    Array[File] background_cram_index_files = []

    # Call variants args
    File model_onnx
    File? model_serialized
    Int? optimization_level
    Boolean output_call_variants_tfrecords = false

    # PostProcessing args
    Int min_variant_quality_hmer_indels = 5
    Int min_variant_quality_non_hmer_indels = 0
    Int min_variant_quality_snps = 0
    Int min_variant_quality_exome_hmer_indels = 20
    Int hard_qual_filter = 1
    Float? allele_frequency_ratio
    Float? h_indel_vaf_to_pass
    Float? h_indel_allele_frequency_ratio
    String ug_post_processing_extra_args = ""

    String dummy_input_for_call_caching = ""

    # Annotation args
    String? input_flow_order
    File exome_intervals
    Array[File]? annotation_intervals
    File ref_dbsnp
    File ref_dbsnp_index

    # Runtime args
    Float? ug_make_examples_memory_override
    Int? ug_make_examples_cpus_override
    Int preemptible_tries = 1
    Int? ug_call_variants_extra_mem
    String call_variants_gpu_type = "nvidia-l4" # For AWS
    Int call_variants_gpus = 1
    Int call_variants_cpus = 8
    Int call_variants_threads = 8
    Int call_variants_uncompr_buf_size_gb = 1
    Boolean? no_address_override

    # Used for running on other clouds (aws)
    String? cloud_provider_override
    File? monitoring_script_input


   # winval validations
   #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
   #@wv min_fraction_hmer_indels <= 1.1 and min_fraction_hmer_indels >= 0
   #@wv min_fraction_non_hmer_indels <= 1.1 and min_fraction_non_hmer_indels >= 0
   #@wv min_fraction_snps <= 1.1 and min_fraction_snps >= 0
   #@wv min_fraction_single_strand_non_snps <= 1.1 and min_fraction_single_strand_non_snps >= 0
   #@mv min_read_count_snps >= 1
   #@mv min_read_count_hmer_indels >= 1
   #@mv min_read_count_non_hmer_indels >= 1
   #@wv cloud_provider_override == "aws" -> suffix(cram_files) <= {".cram"}
   #@wv cloud_provider_override == "aws" -> suffix(cram_index_files) <= {".crai", ".csi"}
   #@wv cloud_provider_override == "gcp" -> suffix(cram_files) <= {".cram", ".bam"}
   #@wv cloud_provider_override == "gcp" -> suffix(cram_index_files) <= {".crai", ".bai", ".csi"}
   #@wv prefix(cram_index_files) == cram_files
   #@wv len(cram_files) >= 0
   #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
   #@wv suffix(references['ref_dict']) == '.dict'
   #@wv suffix(references['ref_fasta_index']) == '.fai'
   #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
   #@wv len(background_cram_files) == len(background_cram_index_files)
   #@wv cloud_provider_override == "aws" and len(background_cram_files) > 0 ->  suffix(background_cram_files) <= {".cram"}
   #@wv cloud_provider_override == "aws" and len(background_cram_files) > 0 ->  suffix(background_cram_index_files) <= {".crai", ".csi"}
   #@wv cloud_provider_override == "gcp" and len(background_cram_files) > 0 ->  suffix(background_cram_files) <= {".cram", ".bam"}
   #@wv cloud_provider_override == "gcp" and len(background_cram_files) > 0 ->  suffix(background_cram_index_files) <= {".crai", ".bai", ".csi"}
   #@wv len(optimal_coverages) == 1 + (len(background_cram_files) > 0)
   #@wv is_somatic -> len(background_cram_files) > 0
   #@wv is_somatic -> defined(allele_frequency_ratio)
   
  }
  meta {
      description:"Performs variant calling on an input cram, using a re-write of (DeepVariant)[https://www.nature.com/articles/nbt.4235] which is adapted for Ultima Genomics data. There are three stages to the variant calling: (1) make_examples - Looks for “active regions” with potential candidates. Within these regions, it performs local assembly (haplotypes), re-aligns the reads, and defines candidate variant. Images of the reads in the vicinity of the candidates are saved as protos in a tfrecord format. (2) call_variants - Collects the images from make_examples and uses a deep learning model to infer the statistics of each variant (i.e. quality, genotype likelihoods etc.). (3) post_process - Uses the output of call_variants to generate a vcf and annotates it."
      author: "Ultima Genomics"
      WDL_AID: {
          exclude: ["pipeline_version",
              "cloud_provider_override",
              "monitoring_script_input",
              "dummy_input_for_call_caching",
              "UGCallVariants.disk_size",
              "UGPostProcessing.disk_size",
              "MergeRealignedCrams.cache_tarball",
              "Globals.glob",
              "Sentieon.Globals.glob",
              "AnnotateVCF.Globals.glob",
              "SingleSampleQC.Globals.glob",
              "VariantCallingEvaluation.Globals.glob",
              "QCReport.disk_size",
              "UGMakeExamples.count_candidates_with_dvtools"
          ]}
  }

  parameter_meta {
    base_file_name: {
      type: "String",
      help: "Prefix for name of all output files",
      category: "input_required"
    }
    cram_files: {
      type: "Array[File]",
      help: "Input cram files. Multiple files are merged.",
      category: "input_required"
    }
    cram_index_files: {
      type: "Array[File]",
      help: "Input cram index files.",
      category: "input_required"
    }
    background_cram_files: {
      type: "Array[File]",
      help: "Background (normal sample) cram files for somatic calling",
      category: "input_optional"
    }
    background_cram_index_files: {
      type: "Array[File]",
      help: "Background (normal sample) cram index files for somatic calling",
      category: "input_optional"
    }

    references: {
      type: "References",
      help: "Reference files: fasta, dict and fai, recommended value set in the template",
      category: "ref_required"
    }
    make_gvcf: {
      type: "Boolean",
      help: "Whether to generate a gvcf. Default: False",
      category: "param_required"
    }
    recalibrate_vaf: {
      type: "Boolean",
      help: "Whether to recalculate the variant allele frequency on the PASS variants, improves over the naive VAF calculate of DeepVariant (somatic calling only)",
      category: "param_advanced"
    }
    is_somatic: {
      help: "Enable somatic calling mode, which enables somatic-specific post-processing options",
      category: "param_required"
    }
    num_shards: {
      type: "Int",
      help: "Maximal number of intervals the genome is broken into when parallelizing the make_examples step",
      category: "param_advanced"
    }
    scatter_intervals_break: {
      type: "Int",
      help: "The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals.",
      category: "param_optional"
    }
    target_intervals: {
      type: "File",
      help: "Limit calling to these regions. If target_intervals and intervals_string are not provided then entire genome is used.",
      category: "param_optional"
    }
    show_bg_fields: {
      help: "Show background fields in the output vcf. Default: false. Mostly relevant for somatic calling.",
      category: "param_optional"
    }
    intervals_string: {
      type: "String",
      help: "Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. hese regions. Takes precedence over target_intervals. If both are not provided then entire genome is used.",
      category: "param_optional"
    }
    min_fraction_snps: {
      type: "Float",
      help: "Minimal fraction of reads, that support a snp, required to  generate a candidate variant",
      category: "param_optional"
    }
    min_fraction_hmer_indels: {
      type: "Float",
      help: "Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant",
      category: "param_optional"
    }
    min_fraction_non_hmer_indels: {
      type: "Float",
      help: "Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant",
      category: "param_optional"
    }
    min_fraction_single_strand_non_snps: {
      help: "In active region detection, in case single_strand_filter is set to true, keep only non snps candidates that have at least min_fraction_single_strand_non_snps fraction",
      category: "param_advanced"
    }
    min_read_count_snps: {
      type: "Int",
      help: "Minimal number of reads, that support a snp, required to  generate a candidate variant",
      category: "param_optional"
    }
    min_read_count_hmer_indels: {
      type: "Int",
      help: "Minimal number of reads, that support an h-mer indel, required to generate a candidate variant",
      category: "param_optional"
    }
    min_read_count_non_hmer_indels: {
      type: "Int",
      help: "Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant",
      category: "param_optional"
    }
    min_base_quality: {
      type: "Int",
      help: "Minimal base quality for candidate generation",
      category: "param_optional"
    }
    pileup_min_mapping_quality: {
      type: "Int",
      help: "Minimal mapping quality to be included in image (the input to the CNN)",
      category: "param_optional"
    }
    candidate_min_mapping_quality: {
      type: "Int",
      help: "Minimal mapping quality for candidate generation",
      category: "param_optional"
    }
    max_reads_per_partition: {
      type: "Int",
      help: "Maximal number of reads that are stored in memory when analyzing an active region",
      category: "param_optional"
    }
    dbg_min_base_quality: {
      type: "Int",
      help: "Minimal base quality for local assembly of haplotypes",
      category: "param_optional"
    }
    prioritize_alt_supporting_reads: {
      type: "Boolean",
      help: "Generate an image with all available alt-supporting reads, and only then add non-supporting reads",
      category: "param_optional"
    }
    p_error: {
      type: "Float",
      help: "Basecalling error for reference confidence model in gvcf",
      category: "param_optional"
    }
    optimal_coverages: {
      type: "Array[Int]",
      help: "Each sample is downsampled to the \"optimal coverage\" (dictated by the coverage of the training set). Downsampling method is determined by cap_at_optimal_coverage.",
      category: "param_advanced"
    }
    cap_at_optimal_coverage: {
      type: "Boolean",
      help: "Defines downsampling behavior. When false, then the reads are downsampled such that the average coverage equals \"optimal coverage\". When true, each position is downsampled to \"optimal coverage\".",
      category: "param_advanced"
    }
    output_realignment: {
      type: "Boolean",
      help: "Output haplotypes and re-aligned reads to a bam file. Default: false.",
      category: "param_optional"
    }
    single_strand_filter : {
      type: "Boolean",
      help: "Whether to filter out non snp candidates that are on a single strand. Reduces the number of candidates and hence the cost. Most useful for somatic calling.",
      category: "param_advanced"
    }
    normalize_strand_bias: {
      type: "Boolean",
      help: "Whether to normalize the strand bias in the images. This is useful for WES calling, where target enrichment may introduce strand bias between alleles.",
      category: "param_advanced"
    }
    strand_bias_normalization_thresholds: {
      type: "Array[Float]",
      help: "Thresholds for strand bias normalization. The first value is the lowest strand bias (further from 1:1) to normalize, the second value is the highest strand bias (closest to 1:1) to normalize",
      category: "param_advanced"
    } 
    keep_duplicates : {
      type: "Boolean",
      help: "Keep duplicated reads in the images. Do not use in high depth samples (e.g. WES).",
      category: "param_advanced"
    }
    add_ins_size_channel : {
      type: "Boolean",
      help: "Use a channel with the insertion size (depends on the model).",
      category: "param_advanced"
    }
    ug_make_examples_extra_args: {
      help: "Additional arguments for make-examples tool",
      category: "param_optional"
    }
    log_make_examples_progress: {
      help: "Cause make_examples to output detailed progress information (for debugging)",
      category: "param_optional"
    }
    germline_vcf:{
      help: "Germline vcf file in order to generate haplotypes that incorporate germline variants",
      category: "param_optional"
    }
    model_onnx: {
      help: "TensorRT model for calling variants (onnx format)",
      category: "ref_required"
    }
    model_serialized: {
      help: "TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided)",
      category: "ref_optional"
    }
    optimization_level: {
      type: "Int",
      help: "Optimization level for TensorRT engine in call_variants",
      category: "param_optional"
    }
    output_call_variants_tfrecords: {
      help: "Output tfrecords from call_variants",
      category: "param_optional"
    }
    min_variant_quality_snps: {
      help: "Minimal snp variant quality in order to be labeled as PASS",
      category: "param_optional"
    }
    min_variant_quality_hmer_indels: {
      help: "Minimal h-mer indel quality in order to be labeled as PASS",
      category: "param_optional"
    }
    min_variant_quality_non_hmer_indels: {
      help: "Minimal non-h-mer indel quality in order to be labeled as PASS",
      category: "param_optional"
    }
    min_variant_quality_exome_hmer_indels: {
      help: "Minimal non-h-mer indel quality in order to be labeled as PASS",
      category: "param_optional"
    }
    allele_frequency_ratio: {
        type: "Float",
        help: "Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering",
        category: "param_optional"
    }
    h_indel_vaf_to_pass: {
        type: "Float",
        help: "Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio",
        category: "param_optional"
    }
    h_indel_allele_frequency_ratio: {
        type: "Float",
        help: "Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering",
        category: "param_optional"
    }
    hard_qual_filter: {
        help: "Any variant with QUAL < hard_qual_filter will be discarded from the VCF file",
        category: "param_optional"
    }
    input_flow_order: {
      help: "Flow order. If not provided, it will be extracted from the CRAM header",
      category: "param_optional"
    }
    exome_intervals: {
      help: "A bed file with exome intervals. Used at the post-processing step to annotate the vcf and modify the FILTER of variants in the exome.",
      category: "ref_required"
    }
    annotation_intervals: {
      help: "List of bed files for VCF annotation",
      category: "ref_optional"
    }
    ref_dbsnp: {
      help: "DbSNP vcf for the annotation of known variants",
      category: "ref_required"
    }
    ref_dbsnp_index: {
      help: "DbSNP vcf index",
      category: "ref_required"
    }
    ug_post_processing_extra_args: {
      help: "Additional arguments for post-processing",
      category: "param_optional"
    }
    ug_make_examples_memory_override: {
      help: "Memory override for make_examples step",
      category: "param_advanced"
    }
    ug_make_examples_cpus_override: {
      help: "CPU number override for make_examples step",
      category: "param_advanced"
    }
    preemptible_tries: {
      help: "Number of preemptible tries",
      category: "param_advanced"
    }
    no_address_override: {
      help: "Whether to disable assigning external IP addresses to VMs (relevant for Google)",
      category: "param_advanced"
    }
    call_variants_uncompr_buf_size_gb: {
      help: "Memory buffer allocated for each uncompression thread in calll_variants",
      category: "param_optional"
    }
    ug_call_variants_extra_mem: {
      help: "Extra memory for call_variants",
      category: "param_advanced"
    }
    call_variants_gpu_type: {
      help: "GPU type for call variants",
      category: "param_optional"
    }
    call_variants_gpus: {
      help: "Number of GPUs for call_variants",
      category: "param_optional"
    }
    call_variants_cpus: {
      help: "Number of CPUs for call_variants",
      category: "param_optional"
    }
    call_variants_threads: {
      help: "Number of decompression threads for call_variants",
      category: "param_optional"
    }
    nvidia_smi_log: {
      help: "Nvidia System Management (nvidia-smi) log",
      category: "output"
    }
    output_vcf: {
      help: "Called variants in vcf format",
      category: "output"
    }
    output_vcf_index: {
      help: "vcf index",
      category: "output"
    }
    vcf_no_ref_calls: {
      help: "Called variants without reference calls",
      type: "File",
      category: "output"
    }
    vcf_no_ref_calls_index: {
      help: "vcf without references calls index",
      type: "File",
      category: "output"
    }
    call_variants_output_tfrecords: {
      help: "The tfrecords that call_variants outputs",
      category: "output"
    }
    # output_model_serialized: {  # uncomment to save serialized model
    #   help: "Output model serialized",
    #   category: "output"
    # }
    output_gvcf: {
      help: "Variant in each position (gvcf file)",
      category: "output"
    }
    output_gvcf_index: {
      help: "gvcf index",
      category: "output"
    }
    output_gvcf_hcr: {
      help: "HCR file - callability regions BED file defined from the gVCF",
      category: "output"
    }
    realigned_cram: {
      help: "Realigned reads cram from make_examples",
      category: "output"
    }
    realigned_cram_index: {
      help: "Realigned CRAM index",

      category: "output"
    }
    flow_order: {
      help: "Flow order",
      type: "String",
      category: "output"
    }
    report_html: {
      help: "QC report html",
      type: "File",
      category: "output"
    }
    qc_h5: {
      help: "QC stats in h5 file format",
      type: "File",
      category: "output"
    }
    qc_metrics_h5: {
      help: "QC stats in specific format for UGDV workflow",
      type: "File",
      category: "output"
    }
    custom_annotation_names: {
      help: "List of custom annotations added to the VCF by postprocessing (currently only through the annotation_intervals)",
      type: "Array[String]",
      category: "output"
    }
    num_candidates: {
      help: "Number of candidates that call_variants processed",
      type: "File",
      category: "output"
    }
    num_candidates_as_int: {
      help: "Number of candidates that call_variants processed (as an integer)",
      type: "Int",
      category: "output"
    }
  }
  String cloud_provider = select_first([cloud_provider_override, 'gcp'])
  Float ug_make_examples_memory = select_first([ug_make_examples_memory_override, 4])
  Int ug_make_examples_cpus = select_first([ug_make_examples_cpus_override, 2])
  Boolean no_address = select_first([no_address_override, true])

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers

  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
  if (!defined(target_intervals)){
    call UGGeneralTasks.IntervalListOfGenome as IntervalListOfGenome{
      input:
        ref_fai = references.ref_fasta_index,
        ref_dict = references.ref_dict,
        disk_size = 1,
        preemptible_tries = preemptible_tries,
        docker = global.ubuntu_docker,
        monitoring_script = monitoring_script,
        no_address = no_address
    }
  }

  if (defined(intervals_string)){
    call UGGeneralTasks.IntervalListFromString{
      input:
        intervals_string = select_first([intervals_string]),
        ref_fai = references.ref_fasta_index,
        ref_dict = references.ref_dict,
        disk_size = 1,
        preemptible_tries = preemptible_tries,
        docker = global.ubuntu_docker,
        monitoring_script = monitoring_script,
        no_address = no_address
    }
  }

  File interval_list = select_first([IntervalListFromString.interval_list, IntervalListOfGenome.interval_list, target_intervals])

  String output_prefix = base_file_name #basename(cram_files[0], ".cram")

   call UGGeneralTasks.IntervalListTotalLength as CallingIntervalsLength {
      input:
        interval_list = interval_list,
        docker = global.ubuntu_docker,
        no_address = no_address,
        monitoring_script = monitoring_script
    }

  call UGGeneralTasks.FastaLengthFromIndex as GenomeSize {
      input:
        fasta_index = references.ref_fasta_index,
        docker = global.ubuntu_docker,
        no_address = no_address,
        monitoring_script = monitoring_script
  }

  call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList{
    input:
      interval_list = interval_list,
      scatter_count = num_shards,
      break_bands_at_multiples_of = scatter_intervals_break,
      dummy_input_for_call_caching = dummy_input_for_call_caching,
      docker = global.broad_gatk_docker,
      no_address = no_address,
      monitoring_script = monitoring_script
  }

  scatter (interval in ScatterIntervalList.out){
    call UGDVTasks.UGMakeExamples {
      input:
        interval = interval,
        total_number_of_shards = ScatterIntervalList.interval_count,
        overall_calling_regions_length = CallingIntervalsLength.interval_list_length,
        is_somatic = is_somatic,
        genome_length = GenomeSize.fasta_length,
        references = references,
        cram_files = cram_files,
        cram_index_files = cram_index_files,
        background_cram_files = background_cram_files,
        background_cram_index_files = background_cram_index_files,
        germline_vcf = germline_vcf,
        min_base_quality = min_base_quality,
        pileup_min_mapping_quality = pileup_min_mapping_quality,
        min_read_count_snps = min_read_count_snps,
        min_read_count_hmer_indels = min_read_count_hmer_indels,
        min_read_count_non_hmer_indels = min_read_count_non_hmer_indels,
        min_fraction_snps = min_fraction_snps,
        min_fraction_hmer_indels = min_fraction_hmer_indels,
        min_fraction_non_hmer_indels = min_fraction_non_hmer_indels,
        candidate_min_mapping_quality = candidate_min_mapping_quality,
        max_reads_per_partition = max_reads_per_partition,
        assembly_min_base_quality  = dbg_min_base_quality,
        make_gvcf = make_gvcf,
        p_error = p_error,
        single_strand_filter = single_strand_filter,
        min_fraction_single_strand_non_snps = min_fraction_single_strand_non_snps,
        keep_duplicates = keep_duplicates,
        add_ins_size_channel = add_ins_size_channel,
        extra_args = ug_make_examples_extra_args,
        prioritize_alt_supporting_reads = prioritize_alt_supporting_reads,
        normalize_strand_bias = normalize_strand_bias,
        strand_bias_normalization_thresholds = strand_bias_normalization_thresholds,
        log_progress = log_make_examples_progress,
        docker = global.ug_make_examples_docker,
        optimal_coverages = optimal_coverages,
        cap_at_optimal_coverage = cap_at_optimal_coverage,
        output_realignment = output_realignment,
        monitoring_script = monitoring_script,
        preemptible_tries = preemptible_tries,
        memory = ug_make_examples_memory,
        cpu = ug_make_examples_cpus,
        cloud_provider = cloud_provider,
        no_address = no_address
    }
  }

  Array[File] examples_array = flatten(UGMakeExamples.output_examples)

  call UGDVTasks.UGCallVariants {
    input:
      examples = examples_array,
      model_onnx = model_onnx,
      model_serialized = model_serialized,
      is_somatic = is_somatic,
      docker = global.ug_call_variants_docker,
      call_variants_uncompr_buf_size_gb = call_variants_uncompr_buf_size_gb,
      gpu_type = call_variants_gpu_type,
      num_gpus = call_variants_gpus,
      num_cpus = call_variants_cpus,
      num_threads = call_variants_threads,
      monitoring_script = monitoring_script,
      call_variants_extra_mem = ug_call_variants_extra_mem,
      optimization_level = optimization_level,
      no_address = no_address
  }

  if (make_gvcf){
    Array[File] gvcf_records = select_all(flatten(UGMakeExamples.gvcf_records))
  }

  if ( defined(input_flow_order) == false ) {
    call UGGeneralTasks.ExtractSampleNameFlowOrder as ExtractSampleNameFlowOrder{
        input:
        input_bam         = cram_files[0],
        monitoring_script = monitoring_script,
        preemptible_tries = preemptible_tries,
        docker            = global.broad_gatk_docker,
        references         = references,
        no_address        =  true,
        cloud_provider_override = cloud_provider_override
    }
  }

  String flow_order_ = select_first([input_flow_order, ExtractSampleNameFlowOrder.flow_order])
  Array[File] cram_files_for_post_processing = if recalibrate_vaf then cram_files else []
  Array[File] cram_index_files_for_post_processing = if recalibrate_vaf then cram_index_files else []
  Array[File] background_cram_files_for_post_processing = if  recalibrate_vaf then background_cram_files else []
  Array[File] background_cram_index_files_for_post_processing = if recalibrate_vaf then background_cram_index_files else []
  
  call UGDVTasks.UGPostProcessing {
    input:
      called_records = UGCallVariants.output_records,
      cram_files = cram_files_for_post_processing,
      cram_index_files = cram_index_files_for_post_processing,
      background_cram_files = background_cram_files_for_post_processing,
      background_cram_index_files = background_cram_index_files_for_post_processing,
      gvcf_records = gvcf_records,
      make_gvcf    = make_gvcf,
      recalibrate_vaf = recalibrate_vaf,
      is_somatic = is_somatic,
      ref = references.ref_fasta,
      ref_index = references.ref_fasta_index,
      flow_order = flow_order_,
      exome_intervals = exome_intervals,
      annotation_intervals = annotation_intervals,
      min_variant_quality_exome_hmer_indels = min_variant_quality_exome_hmer_indels,
      min_variant_quality_hmer_indels = min_variant_quality_hmer_indels,
      min_variant_quality_non_hmer_indels = min_variant_quality_non_hmer_indels,
      min_variant_quality_snps = min_variant_quality_snps,
      qual_filter = hard_qual_filter,
      show_bg_fields = show_bg_fields,
      docker = global.ug_make_examples_docker,
      output_prefix = output_prefix,
      dbsnp = ref_dbsnp,
      dbsnp_index = ref_dbsnp_index,
      extra_args = ug_post_processing_extra_args,
      monitoring_script = monitoring_script,
      no_address = no_address
  }

  if (output_realignment){
    call UGGeneralTasks.MergeCramFiles as MergeRealignedCrams{
      input:
        monitoring_script = monitoring_script,
        crams = UGMakeExamples.realigned_cram,
        output_base_name = output_prefix + "_realign",
        sample_name = "realigned",
        references = references,
        docker = global.broad_gatk_docker,
        no_address = no_address,
        cpus = 4,
        preemptible_tries = preemptible_tries
    }
  }

  if (output_call_variants_tfrecords){
    Array[File] call_variants_output_tfrecords_maybe = UGCallVariants.output_records
  }

  File raw_output_vcf = UGPostProcessing.vcf_file
  File raw_output_vcf_index = UGPostProcessing.vcf_index

  if (length(background_cram_files) > 0) {
    call PostProcesTasks.ApplyAlleleFrequencyRatioFilter {
      input:
          input_vcf = raw_output_vcf,
          af_ratio = select_first([allele_frequency_ratio]),
          h_indel_vaf_to_pass = h_indel_vaf_to_pass,
          h_indel_vaf_ratio_to_pass = h_indel_allele_frequency_ratio,
          final_vcf_base_name = base_file_name,
          monitoring_script = monitoring_script,
          ugbio_filtering_docker =  global.ugbio_filtering_docker,
          no_address = no_address
    } 
  }

  call PostProcesTasks.RemoveRefCalls as RemoveRefCalls {
       input:
          input_vcf = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf, raw_output_vcf]),
          final_vcf_base_name = base_file_name,
          monitoring_script = monitoring_script,
          bcftools_docker =  global.bcftools_docker,
          no_address = no_address
  }

  if (make_gvcf) {
    File gvcf_maybe = UGPostProcessing.gvcf_file
    File gvcf_index_maybe = UGPostProcessing.gvcf_file_index
    File gvcf_hcr_maybe = UGPostProcessing.gvcf_hcr
  }

  call UGDVTasks.QCReport {
    input:
      input_vcf = RemoveRefCalls.output_vcf,
      input_vcf_index = RemoveRefCalls.output_vcf_index,
      callable_bed  = gvcf_hcr_maybe,
      output_prefix = output_prefix,
      ref = references.ref_fasta,
      docker = global.ug_make_examples_docker,
      monitoring_script = monitoring_script,
      no_address = no_address
  }

  if (output_realignment) {
    File? realigned_cram_maybe = MergeRealignedCrams.output_cram
    File? realigned_cram_index_maybe = MergeRealignedCrams.output_cram_index
  }



  output 
  {
    File nvidia_smi_log     = UGCallVariants.nvidia_smi_log
    # File output_model_serialized   = UGCallVariants.output_model_serialized # uncomment to output the serilized model
    File output_vcf         = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf, raw_output_vcf])
    File output_vcf_index   = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf_index, raw_output_vcf_index])
    File vcf_no_ref_calls   = RemoveRefCalls.output_vcf
    File vcf_no_ref_calls_index = RemoveRefCalls.output_vcf_index
    Array[File]? call_variants_output_tfrecords = call_variants_output_tfrecords_maybe
    File? output_gvcf       = gvcf_maybe
    File? output_gvcf_index = gvcf_index_maybe
    File? output_gvcf_hcr   = gvcf_hcr_maybe
    File? realigned_cram    = realigned_cram_maybe
    File? realigned_cram_index = realigned_cram_index_maybe
    String flow_order       = flow_order_
    File report_html        = QCReport.qc_report
    File qc_h5              = QCReport.qc_h5
    File qc_metrics_h5      = QCReport.qc_metrics_h5
    Array[File] num_candidates   = UGCallVariants.num_candidates
    Int num_candidates_as_int    = UGCallVariants.num_candidates_as_int
  }
}
