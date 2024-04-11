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
    String pipeline_version = "1.11" # !UnusedDeclaration
    String base_file_name

    # Mandatory inputs
    Array[File] cram_files
    Array[File] cram_index_files
    References references

    Boolean make_gvcf

    # Scatter interval list args
    Int num_shards = 40
    Int scatter_intervals_break = 10000000 # Maximal resolution for scattering intervals
    File? target_intervals

    # Make examples args
    Float min_fraction_hmer_indels = 0.12
    Float min_fraction_non_hmer_indels = 0.06
    Float min_fraction_snps = 0.12
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
    String ug_make_examples_extra_args = "--add-ins-size-channel --add-proxy-support-to-non-hmer-insertion --pragmatic"

    # Background files (for somatic calling)
    Array[File] background_cram_files = []
    Array[File] background_cram_index_files = []

    # Call variants args
    File model_onnx
    File? model_serialized

    # PostProcessing args
    Int min_variant_quality_hmer_indels = 5
    Int min_variant_quality_non_hmer_indels = 0
    Int min_variant_quality_snps = 0
    Int min_variant_quality_exome_hmer_indels = 20
    Int hard_qual_filter = 1
    Float? allele_frequency_ratio
    Boolean show_bg_fields = false  # Show fields from background sample (for somatic calling)

    # Systematic error correction args
    Boolean annotate_systematic_errors = false
    File? hmer_runs_bed
    File? filtering_blacklist_file
    Array[File]? sec_models

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
    String call_variants_gpu_type = "nvidia-tesla-a10g" # For AWS
    Int call_variants_gpus = 1
    Int call_variants_cpus = 8
    Int call_variants_threads = 8
    Int call_variants_uncompr_buf_size_gb = 1

    # Used for running on other clouds (aws)
    String? cloud_provider_override
    File? monitoring_script_input


   # winval validations
   #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
   #@wv min_fraction_hmer_indels <= 1.1 and min_fraction_hmer_indels >= 0
   #@wv min_fraction_non_hmer_indels <= 1.1 and min_fraction_non_hmer_indels >= 0
   #@wv min_fraction_snps <= 1.1 and min_fraction_snps >= 0
   #@mv min_read_count_snps >= 1
   #@mv min_read_count_hmer_indels >= 1
   #@mv min_read_count_non_hmer_indels >= 1
   #@wv suffix(cram_files) <= {".cram"}
   #@wv suffix(cram_index_files) <= {".crai"}
   #@wv prefix(cram_index_files) == cram_files
   #@wv len(cram_files) >= 0
   #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
   #@wv suffix(references['ref_dict']) == '.dict'
   #@wv suffix(references['ref_fasta_index']) == '.fai'
   #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
   #@wv annotate_systematic_errors -> (defined(filtering_blacklist_file) and defined(sec_models) and defined(hmer_runs_bed))
   #@wv len(background_cram_files) == len(background_cram_index_files)
   #@wv len(background_cram_files) > 0 ->  suffix(background_cram_files) <= {".cram"}
   #@wv len(background_cram_files) > 0 ->  suffix(background_cram_index_files) <= {".crai"}
   #@wv len(optimal_coverages) == 1 + (len(background_cram_files) > 0)
   #@wv len(background_cram_files) > 0 -> defined(allele_frequency_ratio)
  }
  meta {
      description:"Performs variant calling on an input cram, using a re-write of (DeepVariant)[https://www.nature.com/articles/nbt.4235] which is adapted for Ultima Genomics data. There are three stages to the variant calling: (1) make_examples - Looks for “active regions” with potential candidates. Within these regions, it performs local assembly (haplotypes), re-aligns the reads, and defines candidate variant. Images of the reads in the vicinity of the candidates are saved as protos in a tfrecord format. (2) call_variants - Collects the images from make_examples and uses a deep learning model to infer the statistics of each variant (i.e. quality, genotype likelihoods etc.). (3) post_process - Uses the output of call_variants to generate a vcf and annotates it."
      author: "Ultima Genomics"
      WDL_AID: {
          exclude: ["pipeline_version",
              "cloud_provider_override",
              "monitoring_script_input",
              "no_address_override",
              "dummy_input_for_call_caching",
              "UGMakeExamples.interval_nreads",
              "UGCallVariants.disk_size",
              "UGPostProcessing.disk_size",
              "CreateSECBlacklist.disk_size",
              "FilterVCF.input_model",
              "FilterVCF.model_name",
              "FilterVCF.flow_order",
              "FilterVCF.annotation_intervals",
              "FilterVCF.is_mutect",
              "FilterVCF.disk_size",
              "MergeRealignedCrams.cache_tarball",
              "Globals.glob",
              "Sentieon.Globals.glob",
              "AnnotateVCF.Globals.glob",
              "SingleSampleQC.Globals.glob",
              "VariantCallingEvaluation.Globals.glob"
          ]}
  }

  parameter_meta {
    base_file_name: {
      help: "Prefix for name of all output files",
      type: "String",
      category: "input_required"
    }
    cram_files: {
      help: "Input cram file. At the moment only a single file is supported.",
      type: "Array[File]",
      category: "input_required"
    }
    cram_index_files: {
      help: "Input cram index file. At the moment only a single file is supported.",
      type: "Array[File]",
      category: "input_required"
    }
    background_cram_files: {
      help: "Background (normal sample) cram files for somatic calling",
      type: "Array[File]?",
      category: "input_optional"
    }
    background_cram_index_files: {
      help: "Background (normal sample) cram index files for somatic calling",
      type: "Array[File]?",
      category: "input_optional"
    }

    references: {
      help: "Reference files: fasta, dict and fai, recommended value set in the template",
      type: "References",
      category: "ref_required"
    }
    make_gvcf: {
      help: "Whether to generate a gvcf. Default: False",
      type: "Boolean",
      category: "param_required"
    }
    num_shards: {
      help: "Maximal number of intervals the genome is broken into when parallelizing the make_examples step",
      type: "Int",
      category: "param_advanced"
    }
    scatter_intervals_break: {
      help: "The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals.",
      type: "Int",
      category: "param_advanced"
    }
    target_intervals: {
      help: "Limit calling to these regions. If not provided then entire genome is used.",
      type: "File?",
      category: "param_optional"
    }
    min_fraction_snps: {
      help: "Minimal fraction of reads, that support a snp, required to  generate a candidate variant",
      type: "Float",
      category: "param_advanced"
    }
    min_fraction_hmer_indels: {
      help: "Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant",
      type: "Float",
      category: "param_advanced"
    }
    min_fraction_non_hmer_indels: {
      help: "Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant",
      type: "Float",
      category: "param_advanced"
    }
    min_read_count_snps: {
      help: "Minimal number of reads, that support a snp, required to  generate a candidate variant",
      type: "Int",
      category: "param_advanced"
    }
    min_read_count_hmer_indels: {
      help: "Minimal number of reads, that support an h-mer indel, required to generate a candidate variant",
      type: "Int",
      category: "param_advanced"
    }
    min_read_count_non_hmer_indels: {
      help: "Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant",
      type: "Int",
      category: "param_advanced"
    }
    min_base_quality: {
      help: "Minimal base quality for candidate generation",
      type: "Int",
      category: "param_advanced"
    }
    pileup_min_mapping_quality: {
      help: "Minimal mapping quality to be included in image (the input to the CNN)",
      type: "Int",
      category: "param_advanced"
    }
    candidate_min_mapping_quality: {
      help: "Minimal mapping quality for candidate generation",
      type: "Int",
      category: "param_advanced"
    }
    max_reads_per_partition: {
      help: "Maximal number of reads that are stored in memory when analyzing an active region",
      type: "Int",
      category: "param_advanced"
    }
    dbg_min_base_quality: {
      help: "Minimal base quality for local assembly of haplotypes",
      type: "Int",
      category: "param_advanced"
    }
    prioritize_alt_supporting_reads: {
      help: "Generate an image with all available alt-supporting reads, and only then add non-supporting reads",
      type: "Boolean",
      category: "param_advanced"
    }
    p_error: {
      help: "Basecalling error for reference confidence model in gvcf",
      type: "Float",
      category: "param_advanced"
    }
    optimal_coverages: {
      help: "Each sample is downsampled to the \"optimal coverage\" (dictated by the coverage of the training set). Downsampling method is determined by cap_at_optimal_coverage.",
      type: "Array[Int]",
      category: "param_advanced"
    }
    cap_at_optimal_coverage: {
      help: "Defines downsampling behavior. When false, then the reads are downsampled such that the average coverage equals \"optimal coverage\". When true, each position is downsampled to \"optimal coverage\".",
      type: "Boolean",
      category: "param_advanced"
    }
    output_realignment: {
      help: "Output haplotypes and re-aligned reads to a bam file. Default: false.",
      type: "Boolean",
      category: "param_optional"
    }
    ug_make_examples_extra_args: {
      help: "Additional arguments for make-examples tool",
      type: "String?",
      category: "param_advanced"
    }
    model_onnx: {
      help: "Machine-learning model for calling variants (onnx format)",
      type: "File",
      category: "ref_required"
    }
    model_serialized: {
      help: "Machine-learning model for calling variants, serialized for a specific platform (it is regenerated if not provided)",
      type: "File?",
      category: "ref_optional"
    }
    min_variant_quality_snps: {
      help: "Minimal snp variant quality in order to be labeled as PASS",
      type: "Int",
      category: "param_optional"
    }
    min_variant_quality_hmer_indels: {
      help: "Minimal h-mer indel quality in order to be labeled as PASS",
      type: "Int",
      category: "param_optional"
    }
    min_variant_quality_non_hmer_indels: {
      help: "Minimal non-h-mer indel quality in order to be labeled as PASS",
      type: "Int",
      category: "param_optional"
    }
    min_variant_quality_exome_hmer_indels: {
      help: "Minimal non-h-mer indel quality in order to be labeled as PASS",
      type: "Int",
      category: "param_optional"
    }
    allele_frequency_ratio: {
        type: "Float",
        help: "Minimal ratio between the allele frequency in tumor and normal, for vcf filtering",
        category: "param_optional"
    }
    hard_qual_filter: {
        type: "Int",
        help: "Any variant with QUAL < hard_qual_filter will be discarded from the VCF file",
        category: "param_optional"
    }
    show_bg_fields: {
      help: "Show background statistics BG_AD, BG_SB in the output VCF (relevant for somatic calling)",
      type: "Boolean",
      category: "param_optional"
    }
    annotate_systematic_errors: {
      help: "Should systematic errors be annotated from a database of common systematic errors",
      type: "Boolean",
      category: "param_advanced"
    }
    filtering_blacklist_file: {
      help: "Database of known positions with systematic errors",
      type: "File?",
      category: "param_advanced"
    }
    sec_models: {
      help: "Models to annotate systematic errors",
      type: "Array[File]?",
      category: "param_advanced"
    }
    hmer_runs_bed: {
      help: "Bed file annotating all homopolymer runs longer than 7 in the reference genome. Used to annotate potentially difficult to sequence regions",
      type: "File?",
      category: "param_optional"
    }
    input_flow_order: {
      help: "Flow order. If not provided, it will be extracted from the CRAM header",
      type: "String?",
      category: "param_optional"
    }
    exome_intervals: {
      help: "A bed file with exome intervals",
      type: "File",
      category: "ref_required"
    }
    annotation_intervals: {
      help: "List of bed files for VCF annotation",
      type: "Array[File]?",
      category: "ref_optional"
    }
    ref_dbsnp: {
      help: "DbSNP vcf for the annotation of known variants",
      type: "File",
      category: "ref_required"
    }
    ref_dbsnp_index: {
      help: "DbSNP vcf index",
      type: "File",
      category: "ref_required"
    }
    ug_make_examples_memory_override: {
      help: "Memory override for make_examples step",
      type: "Float?",
      category: "param_advanced"
    }
    ug_make_examples_cpus_override: {
      help: "CPU number override for make_examples step",
      type: "Int?",
      category: "param_advanced"
    }
    preemptible_tries: {
      help: "Number of preemptible tries",
      type: "Int",
      category: "param_advanced"
    }
    call_variants_uncompr_buf_size_gb: {
      help: "Memory buffer allocated for each uncompression thread in calll_variants",
      type: "Int",
      category: "param_advanced"
    }
    ug_call_variants_extra_mem: {
      help: "Extra memory for call_variants",
      type: "Int?",
      category: "param_advanced"
    }
    call_variants_gpu_type: {
      help: "GPU type for call variants",
      type: "String",
      category: "param_advanced"
    }
    call_variants_gpus: {
      help: "Number of GPUs for call_variants",
      type: "Int",
      category: "param_advanced"
    }
    call_variants_cpus: {
      help: "Number of CPUs for call_variants",
      type: "Int",
      category: "param_advanced"
    }
    call_variants_threads: {
      help: "Number of decompression threads for call_variants",
      type: "Int",
      category: "param_advanced"
    }
    nvidia_smi_log: {
      help: "Nvidia System Management (nvidia-smi) log",
      type: "File",
      category: "output"
    }
    vcf_file: {
      help: "Called variants in vcf format",
      type: "File",
      category: "output"
    }
    vcf_index: {
      help: "vcf index",
      type: "File",
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
    output_gvcf: {
      help: "Variant in each position (gvcf file)",
      type: "File?",
      category: "output"
    }
    output_gvcf_index: {
      help: "gvcf index",
      type: "File?",
      category: "output"
    }
    realigned_cram: {
      help: "Realigned reads cram from make_examples",
      type: "File?",
      category: "output"
    }
    realigned_cram_index: {
      help: "Realigned CRAM index",
      type: "File?",
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
  }
  String cloud_provider = select_first([cloud_provider_override, 'gcp'])
  Float ug_make_examples_memory = select_first([ug_make_examples_memory_override, 4])
  Int ug_make_examples_cpus = select_first([ug_make_examples_cpus_override, 2])

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
        no_address = true
    }
  }

  File interval_list = select_first([IntervalListOfGenome.interval_list, target_intervals])

  String output_prefix = base_file_name #basename(cram_files[0], ".cram")

   call UGGeneralTasks.IntervalListTotalLength as CallingIntervalsLength {
      input:
        interval_list = interval_list,
        docker = global.ubuntu_docker,
        no_address = true,
        monitoring_script = monitoring_script
    }

  call UGGeneralTasks.FastaLengthFromIndex as GenomeSize {
      input:
        fasta_index = references.ref_fasta_index,
        docker = global.ubuntu_docker,
        no_address = true,
        monitoring_script = monitoring_script
  }

  call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList{
    input:
      interval_list = interval_list,
      scatter_count = num_shards,
      break_bands_at_multiples_of = scatter_intervals_break,
      dummy_input_for_call_caching = dummy_input_for_call_caching,
      docker = global.gitc_docker,
      gitc_path = global.gitc_jar_path,
      no_address = true,
      monitoring_script = monitoring_script
  }

  scatter (interval in ScatterIntervalList.out){
    call UGDVTasks.UGMakeExamples {
      input:
        interval = interval,
        total_number_of_shards = ScatterIntervalList.interval_count,
        overall_calling_regions_length = CallingIntervalsLength.interval_list_length,
        genome_length = GenomeSize.fasta_length,
        references = references,
        cram_files = cram_files,
        cram_index_files = cram_index_files,
        background_cram_files = background_cram_files,
        background_cram_index_files = background_cram_index_files,
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
        ug_make_examples_extra_args = ug_make_examples_extra_args,
        prioritize_alt_supporting_reads = prioritize_alt_supporting_reads,
        docker = global.ug_make_examples_docker,
        optimal_coverages = optimal_coverages,
        cap_at_optimal_coverage = cap_at_optimal_coverage,
        output_realignment = output_realignment,
        monitoring_script = monitoring_script,
        preemptible_tries = preemptible_tries,
        memory = ug_make_examples_memory,
        cpu = ug_make_examples_cpus,
        cloud_provider = cloud_provider
    }
  }

  Array[File] examples_array = flatten(UGMakeExamples.output_examples)

  call UGDVTasks.UGCallVariants {
    input:
      examples = examples_array,
      model_onnx = model_onnx,
      model_serialized = model_serialized,
      docker = global.ug_call_variants_docker,
      call_variants_uncompr_buf_size_gb = call_variants_uncompr_buf_size_gb,
      gpu_type = call_variants_gpu_type,
      num_gpus = call_variants_gpus,
      num_cpus = call_variants_cpus,
      num_threads = call_variants_threads,
      monitoring_script = monitoring_script,
      call_variants_extra_mem = ug_call_variants_extra_mem,
  }

  if (make_gvcf){
    Array[File] gvcf_records = flatten(UGMakeExamples.gvcf_records)
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

  call UGDVTasks.UGPostProcessing {
    input:
      called_records = UGCallVariants.output_records,
      gvcf_records = gvcf_records,
      make_gvcf    = make_gvcf,
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
      monitoring_script = monitoring_script
  }

  if (output_realignment){
    call UGGeneralTasks.MergeCramFiles as MergeRealignedCrams{
      input:
        monitoring_script = monitoring_script,
        crams = UGMakeExamples.realigned_cram,
        output_base_name = output_prefix + "_realign",
        sample_name = "realigned",
        references = references,
        docker = global.ug_vc_docker,
        no_address = true,
        cpus = 4,
        preemptible_tries = preemptible_tries
    }
  }

  if (annotate_systematic_errors) {
    call VCTasks.CreateSECBlacklist {
      input:
       monitoring_script = monitoring_script,
       input_gvcf = UGPostProcessing.vcf_file,
       input_gvcf_index = UGPostProcessing.vcf_index,
       output_blacklist_path = base_file_name + '_sec.pickle',
       blacklist_file = filtering_blacklist_file,
       sec_models = sec_models,
       preemptible_tries = preemptible_tries,
       no_address = true,
       docker = global.ug_vc_docker,
    }

    call VCTasks.FilterVCF as FilterVCF {
      input:
        input_vcf = UGPostProcessing.vcf_file,
        runs_file = select_first([hmer_runs_bed]),
        references = references,
        filter_cg_insertions = false,
        blacklist_file = CreateSECBlacklist.output_blacklist,
        final_vcf_base_name = base_file_name,
        preemptible_tries = preemptible_tries,
        docker = global.ug_vc_docker,
        monitoring_script = monitoring_script,
        no_address = true
    }
  }

  call PostProcesTasks.CalibrateBridgingSnvs as CalibrateBridgingSnvs {
    input:
        input_vcf         = select_first([FilterVCF.output_vcf_filtered, UGPostProcessing.vcf_file]),
        input_vcf_index   = select_first([FilterVCF.output_vcf_filtered_index, UGPostProcessing.vcf_index]),
        references = references,
        final_vcf_base_name = base_file_name,
        monitoring_script = monitoring_script,
        ugvc_docker =  global.ug_vc_docker
  }
  if (length(background_cram_files) > 0) {
    call PostProcesTasks.ApplyAlleleFrequencyRatioFilter as ApplyAlleleFrequencyRatioFilter {
      input:
          input_vcf = CalibrateBridgingSnvs.output_vcf,
          af_ratio = select_first([allele_frequency_ratio]),
          final_vcf_base_name = base_file_name,
          monitoring_script = monitoring_script,
          bcftools_docker =  global.bcftools_docker
    } 
  }

   call PostProcesTasks.RemoveRefCalls as RemoveRefCalls {
        input:
            input_vcf = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf, CalibrateBridgingSnvs.output_vcf]),
            final_vcf_base_name = base_file_name,
            monitoring_script = monitoring_script,
            bcftools_docker =  global.bcftools_docker
    }
  
  if (make_gvcf) {
    File gvcf_maybe = UGPostProcessing.gvcf_file
    File gvcf_index_maybe = UGPostProcessing.gvcf_file_index
  }

  if (output_realignment) {
    File? realigned_cram_maybe = MergeRealignedCrams.output_cram
    File? realigned_cram_index_maybe = MergeRealignedCrams.output_cram_index
  }



  output 
  {
    File nvidia_smi_log     = UGCallVariants.nvidia_smi_log
    # File output_model_serialized   = UGCallVariants.output_model_serialized # uncomment to output the serilized model
    File vcf_file           = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf, CalibrateBridgingSnvs.output_vcf])
    File vcf_index          = select_first([ApplyAlleleFrequencyRatioFilter.output_vcf_index, CalibrateBridgingSnvs.output_vcf_index])
    File vcf_no_ref_calls   = RemoveRefCalls.output_vcf
    File vcf_no_ref_calls_index = RemoveRefCalls.output_vcf_index
    File? output_gvcf       = gvcf_maybe
    File? output_gvcf_index = gvcf_index_maybe
    File? realigned_cram    = realigned_cram_maybe
    File? realigned_cram_index = realigned_cram_index_maybe
    String flow_order       = flow_order_
    File report_html        = UGPostProcessing.qc_report
    File qc_h5              = UGPostProcessing.qc_h5
    File qc_metrics_h5      = UGPostProcessing.qc_metrics_h5
  }
}
