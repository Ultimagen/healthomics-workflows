version 1.0
# LICENSE
#   Copyright 2025 Ultima Genomics
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

# ------------------------------------------------------------------------------
# AUTO-GENERATED FILE - DO NOT EDIT.
#
# Generated from pipelines.yaml (pipeline key: germline_pipeline) by
#   python -m python_libs.wdl_pipeline_gen generate --pipeline germline_pipeline
#
# To change this workflow, edit pipelines.yaml and regenerate. CI runs
#   python -m python_libs.wdl_pipeline_gen check
# which fails if this file and the spec have drifted apart.
#
# generated-from: pipelines.yaml@b5c6e8fb0aa8
# ------------------------------------------------------------------------------

# DESCRIPTION
# Runs the Ultima Genomics germline analysis on a single aligned sample: SNV/indel calling with EfficientDV, CNV, structural variant, short tandem repeat, HLA, pharmacogenomics and segmental duplication analysis. Every stage is switched by its own run_* flag, all of which default to true.

import "efficient_dv.wdl" as EfficientDVSubWF
import "germline_CNV_pipeline.wdl" as CNVSubWF
import "structural_variant_pipeline.wdl" as SVSubWF
import "str_genotyper.wdl" as STRSubWF
import "hla_genotyping.wdl" as HLASubWF
import "pypgx.wdl" as PyPGxSubWF
import "segdup.wdl" as SegDupSubWF

workflow GermlinePipeline {
  input {
    String pipeline_version = "1.36.0" # !UnusedDeclaration

    # ---- GermlinePipeline inputs ----
    Array[File] input_cram_bam_list
    Array[File] input_cram_bam_index_list
    Boolean run_variant_calling = true
    Boolean run_cnv = true
    Boolean run_sv = true
    Boolean run_str = true
    Boolean run_hla = true
    Boolean run_pgx = true
    Boolean run_segdup = true

    # ---- shared across stages ----
    String base_file_name
    String reference_genome = "hg38"
    File? monitoring_script_input
    String? cloud_provider_override
    Boolean? no_address_override

    # ---- VariantCalling (EfficientDV, wdls/efficient_dv.wdl) ----
    Int active_areas_min_base_quality = 5
    Boolean add_ins_size_channel = true
    Float? allele_frequency_ratio
    Array[File]? annotation_intervals
    Array[File] background_cram_files = []
    Array[File] background_cram_index_files = []
    Int call_variants_cpus = 8
    String? call_variants_gpu_type_override
    Int call_variants_gpus = 1
    Int call_variants_threads = 8
    Int call_variants_uncompr_buf_size_gb = 1
    Boolean cap_at_optimal_coverage = false
    Int dbg_min_base_quality = 0
    Boolean? dv_diploid_sampling_in_haplotypes
    String dv_dummy_input_for_call_caching = ""
    Int dv_ensemble_size = 0
    Boolean? dv_include_reference_in_haplotypes
    Int dv_max_reads_per_partition = 1500
    Int dv_min_base_quality = 5
    Float dv_min_fraction_single_strand_non_snps = 0.15
    Int dv_min_mapping_quality = 5
    File? dv_model_onnx
    File? dv_model_serialized
    Int? dv_num_haplotypes
    Int dv_num_shards = 40
    Int dv_preemptible_tries = 1
    File? dv_ref_gbz_for_haplotypes
    File? dv_ref_hapl
    Boolean dv_run_haplotype_sampling = false
    Int dv_scatter_intervals_break = 10000000
    Boolean dv_single_strand_filter = false
    String? dv_ug_make_examples_extra_args
    Int ensemble_reference_rows = 5
    File? germline_vcf
    Array[Int]? gq_bins
    Int? gq_resolution_override
    Float? h_indel_allele_frequency_ratio
    Float? h_indel_vaf_to_pass
    Int hard_qual_filter = 1
    String? input_flow_order
    String? intervals_string
    Boolean keep_duplicates = true
    Boolean log_make_examples_progress = false
    Boolean? make_gvcf
    Float min_fraction_hmer_indels = 0.12
    Float min_fraction_non_hmer_indels = 0.06
    Float min_fraction_snps = 0.12
    Int min_hmer_plus_one_candidate = 7
    Int min_read_count_hmer_indels = 2
    Int min_read_count_non_hmer_indels = 2
    Int min_read_count_snps = 2
    Int min_variant_quality_exome_hmer_indels = 20
    Int min_variant_quality_hmer_indels = 5
    Int min_variant_quality_non_hmer_indels = 0
    Int min_variant_quality_snps = 0
    Boolean? normalize_strand_bias
    Array[Int] optimal_coverages = [50]
    Int? optimization_level
    Boolean output_call_variants_tfrecords = false
    Boolean output_realignment = false
    File? override_target_intervals
    Float p_error = 0.005
    File? pangenome_haplotypes
    File? pangenome_haplotypes_index
    Boolean prioritize_alt_supporting_reads = false
    Boolean prioritize_high_quality_reads = true
    Int random_seed = 42
    Boolean? recalibrate_vaf
    File? ref_dbsnp
    File? ref_dbsnp_index
    Float roh_af_default = 0.6
    File? roh_blacklist_override
    Boolean run_ploidy_estimation = false
    Boolean run_roh = false
    Array[String] sex_chromosomes = ["chrX", "chrY", "X", "Y"]
    Boolean shuffle_all_samples = false
    Array[Float]? strand_bias_normalization_thresholds
    Float strong_call_threshold = 0.995
    Boolean trim_soft_clips = true
    Int? ug_call_variants_extra_mem
    Int? ug_make_examples_cpus_override
    Float? ug_make_examples_memory_override
    String ug_post_processing_extra_args = ""
    Int v_gpu_tile_size = 4

    # ---- CNV (GermlineCNVPipeline, wdls/germline_CNV_pipeline.wdl) ----
    Array[File]? bed_graph
    String? chrX_name_override
    String? chrY_name_override
    Int? cnmops_mapq_override
    Int? cnmops_min_cnv_length_override
    Int? cnmops_min_width_value_override
    Int? cnmops_parallel_override
    Int? cnmops_window_length_override
    Boolean cnv_create_md5_checksum_outputs = false
    Int? cnv_preemptible_tries_override
    Array[Int]? cnvpytor_window_length_override
    File? cohort_reads_count_matrix_override
    Int? cushion_size
    Boolean? disable_mod_cnv
    File? filtering_model
    Int? filtering_model_decision_threshold
    File? ploidy_file
    Boolean? skip_figure_generation
    Boolean? skip_filtering
    File? sv_calls_vcf
    File? sv_calls_vcf_index

    # ---- SV (SVPipeline, wdls/structural_variant_pipeline.wdl) ----
    Int? annotate_variants_cpu_override
    Int? annotate_variants_memory_override
    File? blacklist_bed
    String? config_file_string
    Int? convert_vcf_format_memory_override
    Int? create_assembly_memory_override
    String? exclude_filters
    Int? germline_link_variants_memory_override
    GiraffeParameters? giraffe_parameters
    String? gridss_metrics_interval
    Int? homopolymer_length
    Array[File] input_tumor_crams = []
    Array[File] input_tumor_crams_indexes = []
    File? known_hotspot_file
    Int? max_num_haps
    Int? max_reads_per_working_area
    Int? min_base
    String? min_indel_sc_size_to_include
    Int? min_mapq
    String? min_mismatch_count_to_include
    Int? min_normal_coverage
    File? pon_sgl_file
    File? pon_sv_file
    String? prefilter_query
    Int? realign_mapq
    String? reference_name
    Int? rematching_memory_override
    File? repeat_mask_file
    Boolean? run_giraffe
    Boolean? run_ua
    Boolean sv_create_md5_checksum_outputs = false
    String sv_dummy_input_for_call_caching = ""
    Int? sv_max_reads_per_partition
    Boolean sv_no_address = true
    Int? sv_num_shards
    Int? sv_preemptible_tries_override
    Int? sv_scatter_intervals_break
    Boolean? symbolic_vcf_format
    UaParameters? ua_parameters
    File? wgs_calling_interval_list_override

    # ---- STR (STRGenotyper, wdls/str_genotyper.wdl) ----
    Boolean haploid = false
    Int max_repeat = 100
    Int? memory_gb_override
    Float micro_allele_consensus_ratio = 0.8
    Int micro_allele_min_reads = 10
    Int min_repeat = 1
    Float min_score_ratio = 0.85
    Boolean output_detailed_csv = true
    Boolean output_summary_csv = true
    Int ref_padding = 500
    Boolean report_micro_alleles = false
    Int spanning_flank_bases = 10
    Int str_min_mapping_quality = 1
    Boolean str_no_address = true
    Int str_preemptible_tries = 1
    Int threads = 2
    File? variant_catalog

    # ---- HLA (HLAGenotyping, wdls/hla_genotyping.wdl) ----
    File? graphs_files_tar
    String hla_genotyping_tool = "T1K"
    Int? hla_preemptible_tries_override
    File? t1k_index_tar

    # ---- PGx (PyPGx, wdls/pypgx.wdl) ----
    Array[String]? gene_symbols
    File? input_vcf_file
    File? input_vcf_index_file
    Boolean? pgx_diploid_sampling_in_haplotypes
    Int? pgx_ensemble_size
    Boolean? pgx_include_reference_in_haplotypes
    Int? pgx_min_base_quality
    Float? pgx_min_fraction_single_strand_non_snps
    File? pgx_model_onnx
    Int? pgx_num_haplotypes
    Int pgx_preemptible_tries = 1
    File? pgx_ref_gbz_for_haplotypes
    File? pgx_ref_hapl
    Boolean? pgx_run_haplotype_sampling
    Boolean? pgx_single_strand_filter
    String? pgx_ug_make_examples_extra_args
    Array[File]? ref_files_for_tarball

    # ---- SegDup (SegDupAnalysis, wdls/segdup.wdl) ----
    File? background_bed
    File? cn_model
    File? homology_table
    File? homology_table_index
    Int? n_threads
    File? segdup_model_onnx
    File? segdup_model_serialized
    Boolean segdup_no_address = true
    Int segdup_preemptible_tries = 3
    File? segdup_regions

    # winval validations
    #@wv run_variant_calling -> defined(make_gvcf)
    #@wv run_variant_calling -> defined(recalibrate_vaf)
    #@wv run_variant_calling -> defined(normalize_strand_bias)
    #@wv run_variant_calling -> defined(dv_model_onnx)
    #@wv run_cnv -> defined(skip_filtering)
    #@wv run_cnv -> defined(bed_graph)
    #@wv run_cnv -> defined(ploidy_file)
    #@wv run_cnv -> defined(cushion_size)
    #@wv run_sv -> defined(ua_parameters)
    #@wv run_sv -> defined(min_base)
    #@wv run_sv -> defined(min_mapq)
    #@wv run_sv -> defined(sv_max_reads_per_partition)
    #@wv run_sv -> defined(max_reads_per_working_area)
    #@wv run_sv -> defined(realign_mapq)
    #@wv run_sv -> defined(homopolymer_length)
    #@wv run_sv -> defined(config_file_string)
    #@wv run_sv -> defined(reference_name)
    #@wv run_sv -> defined(run_ua)
    #@wv run_sv -> defined(run_giraffe)
    #@wv run_sv -> defined(symbolic_vcf_format)
    #@wv run_sv -> defined(sv_num_shards)
    #@wv run_sv -> defined(sv_scatter_intervals_break)
    #@wv run_str -> defined(variant_catalog)
    #@wv run_pgx -> defined(gene_symbols)
    #@wv run_pgx -> defined(ref_files_for_tarball)
    #@wv run_segdup -> defined(homology_table)
    #@wv run_segdup -> defined(homology_table_index)
    #@wv run_segdup -> defined(segdup_regions)
    #@wv run_segdup -> defined(background_bed)
    #@wv run_segdup -> defined(cn_model)
    #@wv run_segdup -> defined(n_threads)
    #@wv run_segdup -> defined(segdup_model_onnx)
    #@wv not(" " in base_file_name or "#" in base_file_name or "," in base_file_name)
    #@wv len(input_cram_bam_list) == len(input_cram_bam_index_list)
    #@wv len(input_cram_bam_list) > 0
    #@wv suffix(input_cram_bam_list) <= {".bam", ".cram"}
    #@wv suffix(input_cram_bam_index_list) <= {".bai", ".crai"}
    #@wv reference_genome in {"hg38", "b37", "hg38_taps", "hg38_no_alt", "hg38_nist_v3_with_decoy"}
  }

  meta {
    description: "Runs the Ultima Genomics germline analysis on a single aligned sample: SNV/indel calling with EfficientDV, CNV, structural variant, short tandem repeat, HLA, pharmacogenomics and segmental duplication analysis. Every stage is switched by its own run_* flag, all of which default to true."
    author: "Ultima Genomics"
    WDL_AID: { exclude: [
      "cnv_preemptible_tries_override",
      "dv_dummy_input_for_call_caching",
      "hla_preemptible_tries_override",
      "pgx_preemptible_tries",
      "pipeline_version",
      "str_no_address",
      "str_preemptible_tries",
      "sv_dummy_input_for_call_caching",
      "sv_no_address",
      "sv_preemptible_tries_override"
    ]}
  }

  parameter_meta {
    input_cram_bam_list: {
      help: "Aligned, sorted, duplicate-marked CRAM/BAM file(s) of a single sample",
      type: "Array[File]",
      category: "input_required"
    }
    input_cram_bam_index_list: {
      help: "CRAI/BAI index file(s) matching input_cram_bam_list",
      type: "Array[File]",
      category: "input_required"
    }
    run_variant_calling: {
      help: "Run SNV/indel variant calling (EfficientDV). Needs a CRAM with real coverage: EfficientDV aborts when the measured median coverage is 0",
      type: "Boolean",
      category: "param_optional"
    }
    run_cnv: {
      help: "Run germline CNV calling (cn.MOPS + CNVpytor). Requires bed_graph and ploidy_file",
      type: "Boolean",
      category: "param_optional"
    }
    run_sv: {
      help: "Run structural variant calling",
      type: "Boolean",
      category: "param_optional"
    }
    run_str: {
      help: "Run short tandem repeat genotyping",
      type: "Boolean",
      category: "param_optional"
    }
    run_hla: {
      help: "Run HLA/KIR genotyping",
      type: "Boolean",
      category: "param_optional"
    }
    run_pgx: {
      help: "Run pharmacogenomics (PyPGx) analysis",
      type: "Boolean",
      category: "param_optional"
    }
    run_segdup: {
      help: "Run segmental duplication analysis (parascopy/LPA)",
      type: "Boolean",
      category: "param_optional"
    }
    base_file_name: {
      help: "(shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Prefix for name of all output files",
      type: "String",
      category: "input_required"
    }
    reference_genome: {
      help: "(shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Genome selector: hg38, b37, hg38_taps, hg38_nist_v3, hg38_nist_v3_with_decoy, hg38_no_alt, mm10, mm39. Default to hg38",
      type: "String",
      category: "input_optional"
    }
    monitoring_script_input: {
      help: "(shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Monitoring script override for AWS HealthOmics workflow templates multi-region support",
      type: "File",
      category: "input_optional"
    }
    cloud_provider_override: {
      help: "(shared: VariantCalling, CNV, SV, HLA, PGx, SegDup) cloud_provider_override",
      type: "String",
      category: "input_optional"
    }
    no_address_override: {
      help: "(shared: VariantCalling, CNV, HLA, PGx) Whether to disable assigning external IP addresses to VMs (relevant for Google)",
      type: "Boolean",
      category: "param_advanced"
    }
    active_areas_min_base_quality: {
      help: "Minimum base quality for active areas detection",
      type: "Int",
      category: "param_optional"
    }
    add_ins_size_channel: {
      help: "Use a channel with the insertion size (depends on the model).",
      type: "Boolean",
      category: "param_advanced"
    }
    allele_frequency_ratio: {
      help: "Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering",
      type: "Float",
      category: "param_optional"
    }
    annotation_intervals: {
      help: "List of bed files for VCF annotation",
      type: "Array[File]",
      category: "ref_optional"
    }
    background_cram_files: {
      help: "Background (normal sample) cram files for somatic calling",
      type: "Array[File]",
      category: "input_optional"
    }
    background_cram_index_files: {
      help: "Background (normal sample) cram index files for somatic calling",
      type: "Array[File]",
      category: "input_optional"
    }
    call_variants_cpus: {
      help: "Number of CPUs for call_variants",
      type: "Int",
      category: "param_optional"
    }
    call_variants_gpu_type_override: {
      help: "GPU type for call variants",
      type: "String",
      category: "param_optional"
    }
    call_variants_gpus: {
      help: "Number of GPUs for call_variants",
      type: "Int",
      category: "param_optional"
    }
    call_variants_threads: {
      help: "Number of decompression threads for call_variants",
      type: "Int",
      category: "param_optional"
    }
    call_variants_uncompr_buf_size_gb: {
      help: "Memory buffer allocated for each uncompression thread in calll_variants",
      type: "Int",
      category: "param_optional"
    }
    cap_at_optimal_coverage: {
      help: "Defines downsampling behavior. When false, then the reads are downsampled such that the average coverage equals \"optimal coverage\". When true, each position is downsampled to \"optimal coverage\".",
      type: "Boolean",
      category: "param_advanced"
    }
    dbg_min_base_quality: {
      help: "Minimal base quality for local assembly of haplotypes",
      type: "Int",
      category: "param_optional"
    }
    dv_diploid_sampling_in_haplotypes: {
      help: "(VariantCalling) Use diploid sampling strategy for haplotype selection",
      type: "Boolean",
      category: "param_optional"
    }
    dv_dummy_input_for_call_caching: {
      help: "(VariantCalling) dummy_input_for_call_caching (forwarded to VariantCalling.dummy_input_for_call_caching)",
      type: "String",
      category: "input_optional"
    }
    dv_ensemble_size: {
      help: "(VariantCalling) Number of augmented passes for ensemble inference. Values <= 1 disable ensemble entirely (no augmentation is applied); values >= 2 enable selective ensemble.",
      type: "Int",
      category: "param_optional"
    }
    dv_include_reference_in_haplotypes: {
      help: "(VariantCalling) Include the reference sequence in the sampled haplotypes",
      type: "Boolean",
      category: "param_optional"
    }
    dv_max_reads_per_partition: {
      help: "(VariantCalling) Maximal number of reads that are stored in memory when analyzing an active region",
      type: "Int",
      category: "param_optional"
    }
    dv_min_base_quality: {
      help: "(VariantCalling) Minimal base quality for candidate generation",
      type: "Int",
      category: "param_optional"
    }
    dv_min_fraction_single_strand_non_snps: {
      help: "(VariantCalling) In active region detection, in case single_strand_filter is set to true, keep only non snps candidates that have at least min_fraction_single_strand_non_snps fraction",
      type: "Float",
      category: "param_advanced"
    }
    dv_min_mapping_quality: {
      help: "(VariantCalling) Minimum mapping quality for reads to appear in pileup images (input to CNN) and to be considered as supporting an alt-allele in candidate generation",
      type: "Int",
      category: "param_optional"
    }
    dv_model_onnx: {
      help: "(VariantCalling) TensorRT model for calling variants (onnx format)",
      type: "File",
      category: "ref_optional"
    }
    dv_model_serialized: {
      help: "(VariantCalling) TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided)",
      type: "File",
      category: "ref_optional"
    }
    dv_num_haplotypes: {
      help: "(VariantCalling) Number of haplotypes in the pangenome haplotype CRAM. Also determines the haplotype band height in the pileup image.",
      type: "Int",
      category: "param_optional"
    }
    dv_num_shards: {
      help: "(VariantCalling) Maximal number of intervals the genome is broken into when parallelizing the make_examples step",
      type: "Int",
      category: "param_advanced"
    }
    dv_preemptible_tries: {
      help: "(VariantCalling) Number of preemptible tries",
      type: "Int",
      category: "param_advanced"
    }
    dv_ref_gbz_for_haplotypes: {
      help: "(VariantCalling) Pangenome GBZ index file for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided)",
      type: "File",
      category: "ref_optional"
    }
    dv_ref_hapl: {
      help: "(VariantCalling) Pre-computed haplotype index file (.hapl) for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided)",
      type: "File",
      category: "ref_optional"
    }
    dv_run_haplotype_sampling: {
      help: "(VariantCalling) Whether to run haplotype sampling to create pangenome haplotypes. Default: false",
      type: "Boolean",
      category: "param_optional"
    }
    dv_scatter_intervals_break: {
      help: "(VariantCalling) The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals.",
      type: "Int",
      category: "param_optional"
    }
    dv_single_strand_filter: {
      help: "(VariantCalling) Whether to filter out non snp candidates that are on a single strand. Reduces the number of candidates and hence the cost. Most useful for somatic calling.",
      type: "Boolean",
      category: "param_advanced"
    }
    dv_ug_make_examples_extra_args: {
      help: "(VariantCalling) Additional arguments for make-examples tool",
      type: "String",
      category: "param_optional"
    }
    ensemble_reference_rows: {
      help: "Number of reference rows for ensemble inference",
      type: "Int",
      category: "param_optional"
    }
    germline_vcf: {
      help: "Germline vcf file in order to generate haplotypes that incorporate germline variants",
      type: "File",
      category: "param_optional"
    }
    gq_bins: {
      help: "GQ bins to use instead of a fixed resolution (overrides gq_resolution)",
      type: "Array[Int]",
      category: "param_optional"
    }
    gq_resolution_override: {
      help: "Override for gq resolution (default: 5)",
      type: "Int",
      category: "param_optional"
    }
    h_indel_allele_frequency_ratio: {
      help: "Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering",
      type: "Float",
      category: "param_optional"
    }
    h_indel_vaf_to_pass: {
      help: "Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio",
      type: "Float",
      category: "param_optional"
    }
    hard_qual_filter: {
      help: "Any variant with QUAL < hard_qual_filter will be discarded from the VCF file",
      type: "Int",
      category: "param_optional"
    }
    input_flow_order: {
      help: "Flow order. If not provided, it will be extracted from the CRAM header",
      type: "String",
      category: "param_optional"
    }
    intervals_string: {
      help: "Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. Takes precedence over override_target_intervals.",
      type: "String",
      category: "param_optional"
    }
    keep_duplicates: {
      help: "Keep duplicated reads in the images. Do not use in high depth samples (e.g. WES).",
      type: "Boolean",
      category: "param_advanced"
    }
    log_make_examples_progress: {
      help: "Cause make_examples to output detailed progress information (for debugging)",
      type: "Boolean",
      category: "param_optional"
    }
    make_gvcf: {
      help: "Whether to generate a gvcf. Default: False",
      type: "Boolean",
      category: "param_optional"
    }
    min_fraction_hmer_indels: {
      help: "Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant",
      type: "Float",
      category: "param_optional"
    }
    min_fraction_non_hmer_indels: {
      help: "Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant",
      type: "Float",
      category: "param_optional"
    }
    min_fraction_snps: {
      help: "Minimal fraction of reads, that support a snp, required to generate a candidate variant",
      type: "Float",
      category: "param_optional"
    }
    min_hmer_plus_one_candidate: {
      help: "Minimal hmer length, above which more 1-bp insertion candidates are generated, provided they also meet allele frequency conditions",
      type: "Int",
      category: "param_optional"
    }
    min_read_count_hmer_indels: {
      help: "Minimal number of reads, that support an h-mer indel, required to generate a candidate variant",
      type: "Int",
      category: "param_optional"
    }
    min_read_count_non_hmer_indels: {
      help: "Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant",
      type: "Int",
      category: "param_optional"
    }
    min_read_count_snps: {
      help: "Minimal number of reads, that support a snp, required to generate a candidate variant",
      type: "Int",
      category: "param_optional"
    }
    min_variant_quality_exome_hmer_indels: {
      help: "Minimal non-h-mer indel quality in order to be labeled as PASS",
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
    min_variant_quality_snps: {
      help: "Minimal snp variant quality in order to be labeled as PASS",
      type: "Int",
      category: "param_optional"
    }
    normalize_strand_bias: {
      help: "Whether to normalize the strand bias in the images. This is useful for WES calling, where target enrichment may introduce strand bias between alleles.",
      type: "Boolean",
      category: "param_advanced"
    }
    optimal_coverages: {
      help: "Each sample is downsampled to the \"optimal coverage\" (dictated by the coverage of the training set). Downsampling method is determined by cap_at_optimal_coverage.",
      type: "Array[Int]",
      category: "param_advanced"
    }
    optimization_level: {
      help: "Optimization level for TensorRT engine in call_variants",
      type: "Int",
      category: "param_optional"
    }
    output_call_variants_tfrecords: {
      help: "Output tfrecords from call_variants",
      type: "Boolean",
      category: "param_optional"
    }
    output_realignment: {
      help: "Output haplotypes and re-aligned reads to a bam file. Default: false.",
      type: "Boolean",
      category: "param_optional"
    }
    override_target_intervals: {
      help: "Override default genome-specific target intervals. If not provided, uses genome-specific default intervals.",
      type: "File",
      category: "param_optional"
    }
    p_error: {
      help: "Basecalling error for reference confidence model in gvcf",
      type: "Float",
      category: "param_optional"
    }
    pangenome_haplotypes: {
      help: "Optional pangenome haplotypes cram file",
      type: "File",
      category: "param_optional"
    }
    pangenome_haplotypes_index: {
      help: "Optional pangenome haplotypes cram index file",
      type: "File",
      category: "param_optional"
    }
    prioritize_alt_supporting_reads: {
      help: "Generate an image with all available alt-supporting reads, and only then add non-supporting reads",
      type: "Boolean",
      category: "param_optional"
    }
    prioritize_high_quality_reads: {
      help: "When min-mapq=0, add mapq=0 reads last, only filling remaining image capacity after high-mapq reads",
      type: "Boolean",
      category: "param_optional"
    }
    random_seed: {
      help: "Random seed for ensemble inference",
      type: "Int",
      category: "param_optional"
    }
    recalibrate_vaf: {
      help: "Whether to recalculate the variant allele frequency on the PASS variants, improves over the naive VAF calculate of DeepVariant (somatic calling only)",
      type: "Boolean",
      category: "param_advanced"
    }
    ref_dbsnp: {
      help: "DbSNP vcf for the annotation of known variants",
      type: "File",
      category: "ref_optional"
    }
    ref_dbsnp_index: {
      help: "DbSNP vcf index",
      type: "File",
      category: "ref_optional"
    }
    roh_af_default: {
      help: "Alternate allele frequency assumed for every marker by the ROH caller, in place of a population frequency table",
      type: "Float",
      category: "param_optional"
    }
    roh_blacklist_override: {
      help: "BED of alignment-artefact regions to exclude from the reported runs of homozygosity, overriding the genome default (ENCODE blacklist v2)",
      type: "File",
      category: "ref_optional"
    }
    run_ploidy_estimation: {
      help: "Run VCF-based ploidy estimation and chrX/Y haploid conversion for germline samples. Default: false; enabled by germline use cases.",
      type: "Boolean",
      category: "param_optional"
    }
    run_roh: {
      help: "Whether to call runs of homozygosity (ROH). Enabled by default in the germline WGS use-cases, off otherwise. Requires a reference genome that has a roh_blacklist resource (the hg38 builds and b37) unless roh_blacklist_override is given",
      type: "Boolean",
      category: "param_optional"
    }
    sex_chromosomes: {
      help: "Sex chromosome names to exclude from autosomal ploidy baseline. Defaults support chr-prefixed and non-prefixed human references.",
      type: "Array[String]",
      category: "param_optional"
    }
    shuffle_all_samples: {
      help: "Whether to shuffle all samples during inference",
      type: "Boolean",
      category: "param_optional"
    }
    strand_bias_normalization_thresholds: {
      help: "Thresholds for strand bias normalization. The first value is the lowest strand bias (further from 1:1) to normalize, the second value is the highest strand bias (closest to 1:1) to normalize",
      type: "Array[Float]",
      category: "param_advanced"
    }
    strong_call_threshold: {
      help: "Probability threshold for selective ensemble inference. When ensemble_size >= 2, examples with max probability below this threshold are re-evaluated using ensemble inference; examples above it are accepted as-is.",
      type: "Float",
      category: "param_optional"
    }
    trim_soft_clips: {
      help: "Trim soft-clipped bases from pileup images",
      type: "Boolean",
      category: "param_optional"
    }
    ug_call_variants_extra_mem: {
      help: "Extra memory for call_variants",
      type: "Int",
      category: "param_advanced"
    }
    ug_make_examples_cpus_override: {
      help: "CPU number override for make_examples step",
      type: "Int",
      category: "param_advanced"
    }
    ug_make_examples_memory_override: {
      help: "Memory override for make_examples step",
      type: "Float",
      category: "param_advanced"
    }
    ug_post_processing_extra_args: {
      help: "Additional arguments for post-processing",
      type: "String",
      category: "param_optional"
    }
    v_gpu_tile_size: {
      help: "Virtual GPU tile size for call_variants",
      type: "Int",
      category: "param_optional"
    }
    bed_graph: {
      help: "Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data).",
      type: "Array[File]",
      category: "input_optional"
    }
    chrX_name_override: {
      help: "Name of the X chromosome in the cohort reads count matrix, needed for the ploidy correction. By default derived from reference_genome ('X' for b37, 'chrX' otherwise)",
      type: "String",
      category: "param_advanced"
    }
    chrY_name_override: {
      help: "Name of the Y chromosome in the cohort reads count matrix, needed for the ploidy correction. By default derived from reference_genome ('Y' for b37, 'chrY' otherwise)",
      type: "String",
      category: "param_advanced"
    }
    cnmops_mapq_override: {
      help: "Reads mapping-quality cutoff for coverage aggregation used in cn.mops, default value is 1",
      type: "Int",
      category: "param_advanced"
    }
    cnmops_min_cnv_length_override: {
      help: "Minimum length for reporting CNV. Default is: 0",
      type: "Int",
      category: "param_advanced"
    }
    cnmops_min_width_value_override: {
      help: "Minimum of consecutive windows with a significant signal to consider for CNV reporting. Default is: 2",
      type: "Int",
      category: "param_advanced"
    }
    cnmops_parallel_override: {
      help: "Number of cpus for cn.mops run. Default value is 4",
      type: "Int",
      category: "param_advanced"
    }
    cnmops_window_length_override: {
      help: "Window length on which the read counts will be aggregated, default value is 1000",
      type: "Int",
      category: "param_advanced"
    }
    cnv_create_md5_checksum_outputs: {
      help: "(CNV) Create md5 checksum for requested output files",
      type: "Boolean",
      category: "input_optional"
    }
    cnv_preemptible_tries_override: {
      help: "(CNV) preemptible_tries_override (forwarded to CNV.preemptible_tries_override)",
      type: "Int",
      category: "input_optional"
    }
    cnvpytor_window_length_override: {
      help: "Window length on which the read counts will be aggregated, default value is [500,2500]",
      type: "Array[Int]",
      category: "param_advanced"
    }
    cohort_reads_count_matrix_override: {
      help: "GenomicRanges object of the cohort reads count matrix in rds file format. By default the cohort matching reference_genome is taken from the genome resources.",
      type: "File",
      category: "input_optional"
    }
    cushion_size: {
      help: "Cushion size around CNV breakpoints for split-read analysis and jump alignment analysis",
      type: "Int",
      category: "param_optional"
    }
    disable_mod_cnv: {
      help: "whether to call moderate cnvs (Fold-Change~1.5 will be tagged as CN2.5 and Fold-Change~0.7 will be tagged as CN1.5). Default is: True",
      type: "Boolean",
      category: "param_advanced"
    }
    filtering_model: {
      help: "CNV filtering model, default in template, calls are not filtered if not provided",
      type: "File",
      category: "input_optional"
    }
    filtering_model_decision_threshold: {
      help: "Decision threshold for the filtering model, default is set in template. Lower- less stringent, Higher- more stringent",
      type: "Int",
      category: "param_optional"
    }
    ploidy_file: {
      help: "X chromosome ploidy of the cohort and the additional sample. Each sample is represented on a number on a separate row. Ploidy of the default cohort can be found in the template. The last row corresponds to the sample being called. Genome independent.",
      type: "File",
      category: "input_optional"
    }
    skip_figure_generation: {
      help: "Skip CNV calls figure generation. Default is: False",
      type: "Boolean",
      category: "param_optional"
    }
    skip_filtering: {
      help: "Whether to skip CNV filtering step, default is False",
      type: "Boolean",
      category: "param_optional"
    }
    sv_calls_vcf: {
      help: "SV calls in VCF format (MANTA-like, single record per SV call) to be used for annotation of combined CNV calls, default is empty and annotation is not performed.<br> The input tested is the output of structrual_variant_pipeline.wdl",
      type: "File",
      category: "input_optional"
    }
    sv_calls_vcf_index: {
      help: "Index file for the SV calls VCF",
      type: "File",
      category: "input_optional"
    }
    annotate_variants_cpu_override: {
      help: "cpu override for annotate_variants task",
      type: "Int",
      category: "advanced"
    }
    annotate_variants_memory_override: {
      help: "memory override for annotate_variants task",
      type: "Int",
      category: "advanced"
    }
    blacklist_bed: {
      help: "Gridss blacklist file. When not provided, resolved per genome from genome_resources (sv_blacklist). Provide empty file to disable",
      type: "File",
      category: "optional"
    }
    config_file_string: {
      help: "Gridss config file content",
      type: "String",
      category: "advanced"
    }
    convert_vcf_format_memory_override: {
      help: "memory override for convert_vcf_format task",
      type: "Int",
      category: "advanced"
    }
    create_assembly_memory_override: {
      help: "memory override for create_assembly task",
      type: "Int",
      category: "advanced"
    }
    exclude_filters: {
      help: "gripss paramter: Exclude filters from the output vcf, separated by ;",
      type: "String",
      category: "optional"
    }
    germline_link_variants_memory_override: {
      help: "memory override for germline_link_variants task",
      type: "Int",
      category: "advanced"
    }
    giraffe_parameters: {
      help: "vg giraffe index files to improve haplotype interpretation using population graphs",
      type: "GiraffeParameters",
      category: "optional"
    }
    gridss_metrics_interval: {
      help: "Interval for collecting gridss metrics",
      type: "String",
      category: "optional"
    }
    homopolymer_length: {
      help: "Realignment parameter: do realignment on homopolymeres longer than this value",
      type: "Int",
      category: "optional"
    }
    input_tumor_crams: {
      help: "Input CRAM file for the tumor (in case of matched T/N calling)",
      type: "Array[File]",
      category: "optional"
    }
    input_tumor_crams_indexes: {
      help: "Input CRAM index for the tumor (in case of matched T/N calling)",
      type: "Array[File]",
      category: "optional"
    }
    known_hotspot_file: {
      help: "gripss paramter: Known locations that are hot spot for SVs (see https://github.com/hartwigmedical/hmftools/tree/master/linx), filtered less stringently",
      type: "File",
      category: "optional"
    }
    max_num_haps: {
      help: "Assembly parameter: Maximum number of haplotypes showing an evidence of SV to report",
      type: "Int",
      category: "required"
    }
    max_reads_per_working_area: {
      help: "Rematching parameter: Maximal number of reads that are stored in memory when rematching reads to haplotypes (similar to max_reads_per_partition in assembly)",
      type: "Int",
      category: "advanced"
    }
    min_base: {
      help: "Assembly parameter: Minimum base quality for using in DeBruijn graph construction. Default value in template",
      type: "Int",
      category: "optional"
    }
    min_indel_sc_size_to_include: {
      help: "Assembly parameter: Minimum size of an indel and soft-clipping in the read to include the read in the assembly. ;-separated between samples",
      type: "String",
      category: "optional"
    }
    min_mapq: {
      help: "Assembly parameter: Minimum mapping quality. Default value in template",
      type: "Int",
      category: "optional"
    }
    min_mismatch_count_to_include: {
      help: "Assembly parameter: Minimal number of counts to require to include the read in the assembly. ;-separated between samples",
      type: "String",
      category: "optional"
    }
    min_normal_coverage: {
      help: "gripss paramter: Minimum coverage in the normal sample to determine somatic status. Default value:8",
      type: "Int",
      category: "optional"
    }
    pon_sgl_file: {
      help: "gripss paramter: Panel of normals for single end breakend (partially resolved) calls. Note that the default value is in template",
      type: "File",
      category: "optional"
    }
    pon_sv_file: {
      help: "gripss paramter: panel of normals for breakpoint (fully resolved) calls. Note that the default value is in template",
      type: "File",
      category: "optional"
    }
    prefilter_query: {
      help: "Expression (in bcftools view format) to filter the variants before annotation",
      type: "String",
      category: "optional"
    }
    realign_mapq: {
      help: "Realignment parameter: Below this value we skip realignment on the supplementary alignment",
      type: "Int",
      category: "optional"
    }
    reference_name: {
      help: "Can be 38 or 19",
      type: "String",
      category: "optional"
    }
    rematching_memory_override: {
      help: "memory override for rematching task",
      type: "Int",
      category: "advanced"
    }
    repeat_mask_file: {
      help: "gripss paramter: Repeat mask file. Note that the default value is in template",
      type: "File",
      category: "optional"
    }
    run_giraffe: {
      help: "Whether to run Giraffe haplotype aware alignment or not",
      type: "Boolean",
      category: "optional"
    }
    run_ua: {
      help: "Whether to run UA realignment on the output of the assembly (helps resolving some deletions) or not",
      type: "Boolean",
      category: "optional"
    }
    sv_create_md5_checksum_outputs: {
      help: "(SV) Create md5 checksum for requested output files",
      type: "Boolean",
      category: "input_optional"
    }
    sv_dummy_input_for_call_caching: {
      help: "(SV) dummy_input_for_call_caching (forwarded to SV.dummy_input_for_call_caching)",
      type: "String",
      category: "input_optional"
    }
    sv_max_reads_per_partition: {
      help: "(SV) Assembly parameter: Maximal number of reads that are stored in memory when analyzing an active region",
      type: "Int",
      category: "advanced"
    }
    sv_no_address: {
      help: "(SV) no_address (forwarded to SV.no_address)",
      type: "Boolean",
      category: "input_optional"
    }
    sv_num_shards: {
      help: "(SV) Relevant for scatter tasks, which are CreateAssembly and gridss.AnnotateVariants",
      type: "Int",
      category: "optional"
    }
    sv_preemptible_tries_override: {
      help: "(SV) preemptible_tries_override (forwarded to SV.preemptible_tries_override)",
      type: "Int",
      category: "input_optional"
    }
    sv_scatter_intervals_break: {
      help: "(SV) Maximal resolution for scattering intervals",
      type: "Int",
      category: "advanced"
    }
    symbolic_vcf_format: {
      help: "Whether to convert the output vcf to the region format or not, default True",
      type: "Boolean",
      category: "optional"
    }
    ua_parameters: {
      help: "UA parameters: v_aware_alignment_flag and ua_extra_args, recommended value set in the template",
      type: "UaParameters",
      category: "optional"
    }
    wgs_calling_interval_list_override: {
      help: "Optional override for the interval list defining the region to perform variant calling on. When not provided, resolved per genome from genome_resources (calling_interval_list_without_artefacts)",
      type: "File",
      category: "optional"
    }
    haploid: {
      help: "Enable haploid mode: report single allele instead of diploid pairs. Use for X/Y chromosomes in males or haploid organisms.",
      type: "Boolean",
      category: "input_optional"
    }
    max_repeat: {
      help: "Maximum number of repeat units to include in auxiliary reference sequences",
      type: "Int",
      category: "input_optional"
    }
    memory_gb_override: {
      help: "Optional memory allocation override in GB (default: 4)",
      type: "Int",
      category: "input_advanced"
    }
    micro_allele_consensus_ratio: {
      help: "Minimum fraction of supporting spanning reads required to report a micro-allele decimal. Only used when report_micro_alleles is true. Range 0.0-1.0.",
      type: "Float",
      category: "input_optional"
    }
    micro_allele_min_reads: {
      help: "Minimum number of spanning reads supporting the consensus tract length required to report a micro-allele. Guards against low-coverage indel artifacts. Only used when report_micro_alleles is true. Default: 10.",
      type: "Int",
      category: "input_optional"
    }
    min_repeat: {
      help: "Minimum number of repeat units to include in auxiliary reference sequences",
      type: "Int",
      category: "input_optional"
    }
    min_score_ratio: {
      help: "Minimum ratio of alignment score to the theoretical maximum score (read_length * match_score). Alignments below this threshold are filtered out. Range: 0.0-1.0, where 1.0 requires perfect alignment.",
      type: "Float",
      category: "input_optional"
    }
    output_detailed_csv: {
      help: "Whether to output detailed per-read CSV file. Set to false for large catalogs to reduce I/O.",
      type: "Boolean",
      category: "input_advanced"
    }
    output_summary_csv: {
      help: "Whether to output summary per-locus CSV file. Set to false for large catalogs to reduce I/O.",
      type: "Boolean",
      category: "input_optional"
    }
    ref_padding: {
      help: "Number of bases to extend around the STR repeat region when building auxiliary references for alignment. Larger values provide more flanking sequence context for accurate alignment.",
      type: "Int",
      category: "input_optional"
    }
    report_micro_alleles: {
      help: "Report micro-alleles (e.g. 15.3) for haploid loci when a partial-repeat insertion is present. REPCN stays the integer floor; the micro-allele decimal appears in the genotype (GT) and a new RCMA VCF/BED field. No-op for diploid loci. Default: off.",
      type: "Boolean",
      category: "input_optional"
    }
    spanning_flank_bases: {
      help: "Minimum number of bases that must align on each side of the STR repeat region for a read to be considered 'spanning' the locus",
      type: "Int",
      category: "input_optional"
    }
    str_min_mapping_quality: {
      help: "(STR) Minimum mapping quality for reads to be included in analysis",
      type: "Int",
      category: "input_optional"
    }
    str_no_address: {
      help: "(STR) no_address (forwarded to STR.no_address)",
      type: "Boolean",
      category: "input_optional"
    }
    str_preemptible_tries: {
      help: "(STR) Number of preemptible tries before running on non-preemptible",
      type: "Int",
      category: "input_advanced"
    }
    threads: {
      help: "Number of threads for parallel processing",
      type: "Int",
      category: "input_advanced"
    }
    variant_catalog: {
      help: "Variant catalog (json). Example: https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_with_offtargets.GRCh38.json",
      type: "File",
      category: "input_optional"
    }
    graphs_files_tar: {
      help: "HLA-LA graphs files tar (required if using HLA-LA)",
      type: "File",
      category: "optional"
    }
    hla_genotyping_tool: {
      help: "HLA genotyping tool to use. Options: 'HLA-LA' or 'T1K'",
      type: "String",
      category: "required"
    }
    hla_preemptible_tries_override: {
      help: "(HLA) preemptible_tries_override (forwarded to HLA.preemptible_tries_override)",
      type: "Int",
      category: "input_optional"
    }
    t1k_index_tar: {
      help: "T1K index tar.gz containing hlaidx/ and kiridx/ directories with all index files (_seq.fa and _coord.fa). Required if using T1K.",
      type: "File",
      category: "optional"
    }
    gene_symbols: {
      help: "List of gene symbols to analyze",
      type: "Array[String]",
      category: "input_optional"
    }
    input_vcf_file: {
      help: "Input VCF file with variants. Use of high quality variants (i.e. PASS). If not provided, Efficient DV will be run",
      type: "File",
      category: "input_optional"
    }
    input_vcf_index_file: {
      help: "Input VCF index file",
      type: "File",
      category: "input_optional"
    }
    pgx_diploid_sampling_in_haplotypes: {
      help: "(PGx) EfficientDV: use diploid sampling for haplotypes",
      type: "Boolean",
      category: "input_advanced"
    }
    pgx_ensemble_size: {
      help: "(PGx) EfficientDV: number of augmented passes for ensemble inference (0 disables)",
      type: "Int",
      category: "input_advanced"
    }
    pgx_include_reference_in_haplotypes: {
      help: "(PGx) EfficientDV: include the reference among sampled haplotypes",
      type: "Boolean",
      category: "input_advanced"
    }
    pgx_min_base_quality: {
      help: "(PGx) EfficientDV: minimum base quality",
      type: "Int",
      category: "input_advanced"
    }
    pgx_min_fraction_single_strand_non_snps: {
      help: "(PGx) EfficientDV: minimum fraction of single-strand support for non-SNP variants",
      type: "Float",
      category: "input_advanced"
    }
    pgx_model_onnx: {
      help: "(PGx) TensorRT model for calling variants (onnx format)",
      type: "File",
      category: "input_optional"
    }
    pgx_num_haplotypes: {
      help: "(PGx) EfficientDV: number of haplotypes to sample",
      type: "Int",
      category: "input_advanced"
    }
    pgx_preemptible_tries: {
      help: "(PGx) Number of preemptible tries",
      type: "Int",
      category: "param_advanced"
    }
    pgx_ref_gbz_for_haplotypes: {
      help: "(PGx) EfficientDV: pangenome graph (.gbz) file, required when run_haplotype_sampling is true",
      type: "File",
      category: "input_advanced"
    }
    pgx_ref_hapl: {
      help: "(PGx) EfficientDV: pangenome haplotypes (.hapl) file, required when run_haplotype_sampling is true",
      type: "File",
      category: "input_advanced"
    }
    pgx_run_haplotype_sampling: {
      help: "(PGx) EfficientDV: run haplotype sampling to create pangenome haplotypes (enabled by the PE use-case)",
      type: "Boolean",
      category: "input_advanced"
    }
    pgx_single_strand_filter: {
      help: "(PGx) EfficientDV: enable single-strand filtering (default set by template)",
      type: "Boolean",
      category: "input_advanced"
    }
    pgx_ug_make_examples_extra_args: {
      help: "(PGx) EfficientDV: extra args passed to make_examples (channels / max-ins-size for PE)",
      type: "String",
      category: "input_advanced"
    }
    ref_files_for_tarball: {
      help: "List of references for CreateReferenceCache task.",
      type: "Array[File]",
      category: "input_optional"
    }
    background_bed: {
      help: "Background regions (non-segmental duplicated) for CNV calling, see template and `parascopy`",
      type: "File",
      category: "input_advanced"
    }
    cn_model: {
      help: "CNV model file from parascopy, see template",
      type: "File",
      category: "input_advanced"
    }
    homology_table: {
      help: "Segmental duplication table (see parascopy), see template",
      type: "File",
      category: "input_advanced"
    }
    homology_table_index: {
      help: "Segmental duplication table index (see parascopy), see template",
      type: "File",
      category: "input_advanced"
    }
    n_threads: {
      help: "Number of threads to use",
      type: "Int",
      category: "input_optional"
    }
    segdup_model_onnx: {
      help: "(SegDup) DeepVariant model for variant calling on segmental duplications, see template",
      type: "File",
      category: "input_advanced"
    }
    segdup_model_serialized: {
      help: "(SegDup) Serialized model for variant calling",
      type: "File",
      category: "input_advanced"
    }
    segdup_no_address: {
      help: "(SegDup) Start instances with no public IP address",
      type: "Boolean",
      category: "input_advanced"
    }
    segdup_preemptible_tries: {
      help: "(SegDup) Number of preemptible tries",
      type: "Int",
      category: "input_optional"
    }
    segdup_regions: {
      help: "Segmental duplication regions. All reads will be remapped to `segdup_regions` and the calling will happen only on these regions, see template (BED file)",
      type: "File",
      category: "input_advanced"
    }
    dv_nvidia_smi_log: {
      help: "Nvidia System Management (nvidia-smi) log (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    snv_indel_vcf: {
      help: "Called variants in vcf format (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    snv_indel_vcf_index: {
      help: "vcf index (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    snv_indel_vcf_no_ref_calls: {
      help: "Called variants without reference calls (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    snv_indel_vcf_no_ref_calls_index: {
      help: "vcf without references calls index (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_roh_tsv: {
      help: "Runs of homozygosity, as the regions tsv of bcftools roh (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_call_variants_output_tfrecords: {
      help: "The tfrecords that call_variants outputs (only when run_variant_calling)",
      type: "Array[File]",
      category: "output"
    }
    gvcf: {
      help: "Variant in each position (gvcf file) (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    gvcf_index: {
      help: "gvcf index (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_output_gvcf_hcr: {
      help: "HCR file - callability regions BED file defined from the gVCF (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_realigned_cram: {
      help: "Realigned reads cram from make_examples (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_realigned_cram_index: {
      help: "Realigned CRAM index (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_flow_order: {
      help: "Flow order (only when run_variant_calling)",
      type: "String",
      category: "output"
    }
    qc_report_html: {
      help: "QC report html (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    qc_report_h5: {
      help: "QC stats in h5 file format (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_qc_metrics_h5: {
      help: "QC stats in specific format for UGDV workflow (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    dv_num_candidates: {
      help: "Number of candidates that call_variants processed (only when run_variant_calling)",
      type: "Array[File]",
      category: "output"
    }
    dv_num_candidates_as_int: {
      help: "Number of candidates that call_variants processed (as an integer) (only when run_variant_calling)",
      type: "Int",
      category: "output"
    }
    dv_ploidy_report: {
      help: "Human-readable genome ploidy report with sex karyotype, per-chromosome ploidy, and BAF summary when available (only when run_variant_calling)",
      type: "File",
      category: "output"
    }
    cnv_cnmops_cnv_calls_bed: {
      help: "CNMOPS CNV calls in bed format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_cnmops_vcf: {
      help: "CNMOPS CNV calls in VCF format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_cnmops_cnv_calls_vcf_index: {
      help: "Index file for the CNMOPS CNV calls VCF (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_cnvpytor_cnv_calls_bed: {
      help: "CNVpytor CNV calls in bed format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_cnvpytor_vcf: {
      help: "CNVpytor CNV calls in VCF format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_cnvpytor_cnv_calls_vcf_index: {
      help: "Index file for the CNVpytor CNV calls VCF (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_bed: {
      help: "Final (combined) CNV calls in bed format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_vcf: {
      help: "Combined CNV calls in vcf format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_vcf_index: {
      help: "Index of the combined CNV calls in vcf format (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_split_read_evidence: {
      help: "BAM file with split read evidence supporting combined CNV calls (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_split_read_evidence_index: {
      help: "Index file for the BAM with split read evidence supporting combined CNV calls (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_realign_read_evidence: {
      help: "BAM file with read evidence supporting combined CNV calls (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_realign_read_evidence_index: {
      help: "Index file for the BAM with read evidence supporting combined CNV calls (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_combine_read_scores_csv: {
      help: "CSV file with jalign scores for each read (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_combined_coverage_plot: {
      help: "CNV coverage plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_combined_dup_del_plot: {
      help: "Duplication and deletion calls plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_combined_copy_number_plot: {
      help: "Copy number calls plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)",
      type: "File",
      category: "output"
    }
    cnv_md5_checksums_json: {
      help: "json file that will contain md5 checksums for requested output files (only when run_cnv)",
      type: "File",
      category: "output"
    }
    sv_annotated_unlinked_vcf: {
      help: "Annotated VCF file, before GRIPSS or GermlineLinkVariants (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_annotated_unlinked_vcf_index: {
      help: "Annotated VCF index file (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_vcf: {
      help: "Final VCF (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_vcf_index: {
      help: "Final VCF index (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_assembly: {
      help: "Raw assembly - before the realignment (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_assembly_index: {
      help: "Raw assembly - before the realignment - index (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_realigned_assembly: {
      help: "Assembly output after UA realingment (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_realigned_assembly_index: {
      help: "Assembly output index after UA realingment (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_converted_vcf: {
      help: "Final VCF file in the region (non-breakend) format (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_converted_vcf_index: {
      help: "Final VCF index file in the region (non-breakend) format (only when run_sv)",
      type: "File",
      category: "output"
    }
    sv_md5_checksums_json: {
      help: "json file that will contain md5 checksums for requested output files (only when run_sv)",
      type: "File",
      category: "output"
    }
    str_detailed_csv_files: {
      help: "Detailed per-read alignment results in CSV format, containing alignment scores, repeat counts, and read metadata for each alignment. Empty array if output_detailed_csv=false. (only when run_str)",
      type: "Array[File]",
      category: "output"
    }
    str_summary_csv_files: {
      help: "Per-locus summary statistics in CSV format, aggregating alignment results across all reads for each STR locus. Empty array if output_summary_csv=false. (only when run_str)",
      type: "Array[File]",
      category: "output"
    }
    str_bed: {
      help: "Final genotype calls in BED format for visualization in genome browsers (IGV, UCSC). Contains chromosome, start, end, and genotype information (only when run_str)",
      type: "File",
      category: "output"
    }
    str_vcf: {
      help: "Final genotype calls in compressed VCF format with per-allele support counts (ADSP, ADFL). Compatible with standard VCF tools. (only when run_str)",
      type: "File",
      category: "output"
    }
    str_vcf_index: {
      help: "Tabix index for the genotypes VCF file (only when run_str)",
      type: "File",
      category: "output"
    }
    hla_genotypes: {
      help: "HLA genotyping output file (only when run_hla)",
      type: "File",
      category: "output"
    }
    kir_genotypes: {
      help: "KIR genotyping output file (only when run_hla)",
      type: "File",
      category: "output"
    }
    pgx_allele_fraction_profiles: {
      help: "Allele fraction profiles for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_alleles: {
      help: "Alleles for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_cnv_calls: {
      help: "CNV calls for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_consolidated_variants: {
      help: "Consolidated variants for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_copy_number_profiles: {
      help: "Copy number profiles for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_copy_numbers: {
      help: "Copy numbers for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_genotypes: {
      help: "Genotypes for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_imported_variants: {
      help: "Imported variants for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_phased_variants: {
      help: "Phased variants for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_phenotypes: {
      help: "Phenotypes for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_read_depths: {
      help: "Read depths for each gene (only when run_pgx)",
      type: "Array[File]",
      category: "output"
    }
    pgx_vcf: {
      help: "Output VCF file (either the input VCF file or the one produced by Efficient DV if no input VCF file was provided) (only when run_pgx)",
      type: "File",
      category: "output"
    }
    pgx_vcf_index: {
      help: "Output VCF index file (either the input VCF index file or the one produced by Efficient DV if no input VCF index file was provided) (only when run_pgx)",
      type: "File",
      category: "output"
    }
    pgx_results: {
      help: "Results for each gene (only when run_pgx)",
      type: "File",
      category: "output"
    }
    segdup_remap_bam: {
      help: "Remapped BAM file (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_remap_bam_index: {
      help: "Remapped BAM index file (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_acnv_calls: {
      help: "CNV calls (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_pcnv_calls: {
      help: "Paralog CNV calls (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_small_variants: {
      help: "Small variants (VCF) combining ParascopyCall output and LPA KIV-2 targeted small variants (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_small_variants_index: {
      help: "Small variants index (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_lpa_vcf: {
      help: "LPA KIV-2 targeted caller VCF (full: KIV-2 CNV symbolic record + LPA small variants) (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_lpa_vcf_index: {
      help: "LPA KIV-2 targeted caller VCF index (only when run_segdup)",
      type: "File",
      category: "output"
    }
    segdup_lpa_json: {
      help: "LPA KIV-2 targeted caller JSON report (only when run_segdup)",
      type: "File",
      category: "output"
    }
  }

  File aligned_cram = input_cram_bam_list[0]
  File aligned_cram_index = input_cram_bam_index_list[0]
  if (run_variant_calling) {
    call EfficientDVSubWF.EfficientDV as VariantCalling {
      input:
        cram_files = input_cram_bam_list,
        cram_index_files = input_cram_bam_index_list,
        is_somatic = false,
        show_bg_fields = false,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        no_address_override = no_address_override,
        make_gvcf = select_first([make_gvcf]),
        recalibrate_vaf = select_first([recalibrate_vaf]),
        run_haplotype_sampling = dv_run_haplotype_sampling,
        num_shards = dv_num_shards,
        scatter_intervals_break = dv_scatter_intervals_break,
        override_target_intervals = override_target_intervals,
        intervals_string = intervals_string,
        min_fraction_hmer_indels = min_fraction_hmer_indels,
        min_fraction_non_hmer_indels = min_fraction_non_hmer_indels,
        min_fraction_snps = min_fraction_snps,
        min_fraction_single_strand_non_snps = dv_min_fraction_single_strand_non_snps,
        min_read_count_snps = min_read_count_snps,
        min_read_count_hmer_indels = min_read_count_hmer_indels,
        min_read_count_non_hmer_indels = min_read_count_non_hmer_indels,
        min_base_quality = dv_min_base_quality,
        min_mapping_quality = dv_min_mapping_quality,
        min_hmer_plus_one_candidate = min_hmer_plus_one_candidate,
        max_reads_per_partition = dv_max_reads_per_partition,
        dbg_min_base_quality = dbg_min_base_quality,
        prioritize_alt_supporting_reads = prioritize_alt_supporting_reads,
        active_areas_min_base_quality = active_areas_min_base_quality,
        prioritize_high_quality_reads = prioritize_high_quality_reads,
        trim_soft_clips = trim_soft_clips,
        p_error = p_error,
        gq_resolution_override = gq_resolution_override,
        gq_bins = gq_bins,
        optimal_coverages = optimal_coverages,
        cap_at_optimal_coverage = cap_at_optimal_coverage,
        output_realignment = output_realignment,
        single_strand_filter = dv_single_strand_filter,
        keep_duplicates = keep_duplicates,
        add_ins_size_channel = add_ins_size_channel,
        ug_make_examples_extra_args = dv_ug_make_examples_extra_args,
        log_make_examples_progress = log_make_examples_progress,
        normalize_strand_bias = select_first([normalize_strand_bias]),
        strand_bias_normalization_thresholds = strand_bias_normalization_thresholds,
        germline_vcf = germline_vcf,
        pangenome_haplotypes = pangenome_haplotypes,
        pangenome_haplotypes_index = pangenome_haplotypes_index,
        ref_gbz_for_haplotypes = dv_ref_gbz_for_haplotypes,
        ref_hapl = dv_ref_hapl,
        num_haplotypes = dv_num_haplotypes,
        include_reference_in_haplotypes = dv_include_reference_in_haplotypes,
        diploid_sampling_in_haplotypes = dv_diploid_sampling_in_haplotypes,
        background_cram_files = background_cram_files,
        background_cram_index_files = background_cram_index_files,
        model_onnx = select_first([dv_model_onnx]),
        model_serialized = dv_model_serialized,
        optimization_level = optimization_level,
        output_call_variants_tfrecords = output_call_variants_tfrecords,
        run_ploidy_estimation = run_ploidy_estimation,
        sex_chromosomes = sex_chromosomes,
        strong_call_threshold = strong_call_threshold,
        ensemble_size = dv_ensemble_size,
        ensemble_reference_rows = ensemble_reference_rows,
        random_seed = random_seed,
        shuffle_all_samples = shuffle_all_samples,
        min_variant_quality_hmer_indels = min_variant_quality_hmer_indels,
        min_variant_quality_non_hmer_indels = min_variant_quality_non_hmer_indels,
        min_variant_quality_snps = min_variant_quality_snps,
        min_variant_quality_exome_hmer_indels = min_variant_quality_exome_hmer_indels,
        hard_qual_filter = hard_qual_filter,
        allele_frequency_ratio = allele_frequency_ratio,
        h_indel_vaf_to_pass = h_indel_vaf_to_pass,
        h_indel_allele_frequency_ratio = h_indel_allele_frequency_ratio,
        ug_post_processing_extra_args = ug_post_processing_extra_args,
        run_roh = run_roh,
        roh_blacklist_override = roh_blacklist_override,
        roh_af_default = roh_af_default,
        dummy_input_for_call_caching = dv_dummy_input_for_call_caching,
        input_flow_order = input_flow_order,
        annotation_intervals = annotation_intervals,
        ref_dbsnp = ref_dbsnp,
        ref_dbsnp_index = ref_dbsnp_index,
        ug_make_examples_memory_override = ug_make_examples_memory_override,
        ug_make_examples_cpus_override = ug_make_examples_cpus_override,
        preemptible_tries = dv_preemptible_tries,
        ug_call_variants_extra_mem = ug_call_variants_extra_mem,
        call_variants_gpu_type_override = call_variants_gpu_type_override,
        call_variants_gpus = call_variants_gpus,
        call_variants_cpus = call_variants_cpus,
        call_variants_threads = call_variants_threads,
        call_variants_uncompr_buf_size_gb = call_variants_uncompr_buf_size_gb,
        v_gpu_tile_size = v_gpu_tile_size
    }
  }

  if (run_cnv) {
    call CNVSubWF.GermlineCNVPipeline as CNV {
      input:
        input_bam_file = aligned_cram,
        input_bam_file_index = aligned_cram_index,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        no_address_override = no_address_override,
        skip_filtering = select_first([skip_filtering]),
        filtering_model = filtering_model,
        filtering_model_decision_threshold = filtering_model_decision_threshold,
        cnmops_mapq_override = cnmops_mapq_override,
        cnmops_window_length_override = cnmops_window_length_override,
        cnmops_parallel_override = cnmops_parallel_override,
        bed_graph = select_first([bed_graph]),
        cohort_reads_count_matrix_override = cohort_reads_count_matrix_override,
        ploidy_file = select_first([ploidy_file]),
        cnmops_min_width_value_override = cnmops_min_width_value_override,
        cnmops_min_cnv_length_override = cnmops_min_cnv_length_override,
        disable_mod_cnv = disable_mod_cnv,
        chrX_name_override = chrX_name_override,
        chrY_name_override = chrY_name_override,
        cnvpytor_window_length_override = cnvpytor_window_length_override,
        cushion_size = select_first([cushion_size]),
        sv_calls_vcf = sv_calls_vcf,
        sv_calls_vcf_index = sv_calls_vcf_index,
        skip_figure_generation = skip_figure_generation,
        preemptible_tries_override = cnv_preemptible_tries_override,
        create_md5_checksum_outputs = cnv_create_md5_checksum_outputs
    }
  }

  if (run_sv) {
    call SVSubWF.SVPipeline as SV {
      input:
        input_germline_crams = input_cram_bam_list,
        input_germline_crams_indexes = input_cram_bam_index_list,
        is_somatic = false,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        input_tumor_crams = input_tumor_crams,
        input_tumor_crams_indexes = input_tumor_crams_indexes,
        ua_parameters = select_first([ua_parameters]),
        giraffe_parameters = giraffe_parameters,
        wgs_calling_interval_list_override = wgs_calling_interval_list_override,
        min_base = select_first([min_base]),
        min_mapq = select_first([min_mapq]),
        max_reads_per_partition = select_first([sv_max_reads_per_partition]),
        max_reads_per_working_area = select_first([max_reads_per_working_area]),
        max_num_haps = max_num_haps,
        realign_mapq = select_first([realign_mapq]),
        min_indel_sc_size_to_include = min_indel_sc_size_to_include,
        min_mismatch_count_to_include = min_mismatch_count_to_include,
        homopolymer_length = select_first([homopolymer_length]),
        config_file_string = select_first([config_file_string]),
        blacklist_bed = blacklist_bed,
        reference_name = select_first([reference_name]),
        run_ua = select_first([run_ua]),
        run_giraffe = select_first([run_giraffe]),
        prefilter_query = prefilter_query,
        gridss_metrics_interval = gridss_metrics_interval,
        pon_sgl_file = pon_sgl_file,
        pon_sv_file = pon_sv_file,
        repeat_mask_file = repeat_mask_file,
        known_hotspot_file = known_hotspot_file,
        min_normal_coverage = min_normal_coverage,
        exclude_filters = exclude_filters,
        symbolic_vcf_format = select_first([symbolic_vcf_format]),
        num_shards = select_first([sv_num_shards]),
        no_address = sv_no_address,
        preemptible_tries_override = sv_preemptible_tries_override,
        create_assembly_memory_override = create_assembly_memory_override,
        rematching_memory_override = rematching_memory_override,
        annotate_variants_cpu_override = annotate_variants_cpu_override,
        annotate_variants_memory_override = annotate_variants_memory_override,
        convert_vcf_format_memory_override = convert_vcf_format_memory_override,
        germline_link_variants_memory_override = germline_link_variants_memory_override,
        scatter_intervals_break = select_first([sv_scatter_intervals_break]),
        dummy_input_for_call_caching = sv_dummy_input_for_call_caching,
        create_md5_checksum_outputs = sv_create_md5_checksum_outputs
    }
  }

  if (run_str) {
    call STRSubWF.STRGenotyper as STR {
      input:
        cram_file = aligned_cram,
        cram_index = aligned_cram_index,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        variant_catalog = select_first([variant_catalog]),
        ref_padding = ref_padding,
        min_repeat = min_repeat,
        max_repeat = max_repeat,
        min_score_ratio = min_score_ratio,
        spanning_flank_bases = spanning_flank_bases,
        min_mapping_quality = str_min_mapping_quality,
        threads = threads,
        memory_gb_override = memory_gb_override,
        output_detailed_csv = output_detailed_csv,
        output_summary_csv = output_summary_csv,
        haploid = haploid,
        report_micro_alleles = report_micro_alleles,
        micro_allele_consensus_ratio = micro_allele_consensus_ratio,
        micro_allele_min_reads = micro_allele_min_reads,
        preemptible_tries = str_preemptible_tries,
        no_address = str_no_address
    }
  }

  if (run_hla) {
    call HLASubWF.HLAGenotyping as HLA {
      input:
        input_cram_bam = aligned_cram,
        input_cram_bam_index = aligned_cram_index,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        no_address_override = no_address_override,
        graphs_files_tar = graphs_files_tar,
        hla_genotyping_tool = hla_genotyping_tool,
        t1k_index_tar = t1k_index_tar,
        preemptible_tries_override = hla_preemptible_tries_override
    }
  }

  if (run_pgx) {
    call PyPGxSubWF.PyPGx as PGx {
      input:
        cram_file = aligned_cram,
        cram_index_file = aligned_cram_index,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        no_address_override = no_address_override,
        gene_symbols = select_first([gene_symbols]),
        input_vcf_file = input_vcf_file,
        input_vcf_index_file = input_vcf_index_file,
        model_onnx = pgx_model_onnx,
        run_haplotype_sampling = pgx_run_haplotype_sampling,
        ref_hapl = pgx_ref_hapl,
        ref_gbz_for_haplotypes = pgx_ref_gbz_for_haplotypes,
        num_haplotypes = pgx_num_haplotypes,
        include_reference_in_haplotypes = pgx_include_reference_in_haplotypes,
        diploid_sampling_in_haplotypes = pgx_diploid_sampling_in_haplotypes,
        ensemble_size = pgx_ensemble_size,
        single_strand_filter = pgx_single_strand_filter,
        min_base_quality = pgx_min_base_quality,
        min_fraction_single_strand_non_snps = pgx_min_fraction_single_strand_non_snps,
        ug_make_examples_extra_args = pgx_ug_make_examples_extra_args,
        preemptible_tries = pgx_preemptible_tries,
        ref_files_for_tarball = select_first([ref_files_for_tarball])
    }
  }

  if (run_segdup) {
    call SegDupSubWF.SegDupAnalysis as SegDup {
      input:
        input_cram_bam = aligned_cram,
        input_crai_bai = aligned_cram_index,
        base_file_name = base_file_name,
        reference_genome = reference_genome,
        monitoring_script_input = monitoring_script_input,
        cloud_provider_override = cloud_provider_override,
        homology_table = select_first([homology_table]),
        homology_table_index = select_first([homology_table_index]),
        segdup_regions = select_first([segdup_regions]),
        background_bed = select_first([background_bed]),
        cn_model = select_first([cn_model]),
        n_threads = select_first([n_threads]),
        model_onnx = select_first([segdup_model_onnx]),
        model_serialized = segdup_model_serialized,
        preemptible_tries = segdup_preemptible_tries,
        no_address = segdup_no_address
    }
  }

  output {
    File? dv_nvidia_smi_log = VariantCalling.nvidia_smi_log
    File? snv_indel_vcf = VariantCalling.output_vcf
    File? snv_indel_vcf_index = VariantCalling.output_vcf_index
    File? snv_indel_vcf_no_ref_calls = VariantCalling.vcf_no_ref_calls
    File? snv_indel_vcf_no_ref_calls_index = VariantCalling.vcf_no_ref_calls_index
    File? dv_roh_tsv = VariantCalling.roh_tsv
    Array[File]? dv_call_variants_output_tfrecords = VariantCalling.call_variants_output_tfrecords
    File? gvcf = VariantCalling.output_gvcf
    File? gvcf_index = VariantCalling.output_gvcf_index
    File? dv_output_gvcf_hcr = VariantCalling.output_gvcf_hcr
    File? dv_realigned_cram = VariantCalling.realigned_cram
    File? dv_realigned_cram_index = VariantCalling.realigned_cram_index
    String? dv_flow_order = VariantCalling.flow_order
    File? qc_report_html = VariantCalling.report_html
    File? qc_report_h5 = VariantCalling.qc_h5
    File? dv_qc_metrics_h5 = VariantCalling.qc_metrics_h5
    Array[File]? dv_num_candidates = VariantCalling.num_candidates
    Int? dv_num_candidates_as_int = VariantCalling.num_candidates_as_int
    File? dv_ploidy_report = VariantCalling.ploidy_report
    File? cnv_cnmops_cnv_calls_bed = CNV.cnmops_cnv_calls_bed
    File? cnv_cnmops_vcf = CNV.cnmops_cnv_calls_vcf
    File? cnv_cnmops_cnv_calls_vcf_index = CNV.cnmops_cnv_calls_vcf_index
    File? cnv_cnvpytor_cnv_calls_bed = CNV.cnvpytor_cnv_calls_bed
    File? cnv_cnvpytor_vcf = CNV.cnvpytor_cnv_calls_vcf
    File? cnv_cnvpytor_cnv_calls_vcf_index = CNV.cnvpytor_cnv_calls_vcf_index
    File? cnv_bed = CNV.combined_cnv_calls_bed
    File? cnv_vcf = CNV.combined_cnv_calls_bed_vcf
    File? cnv_vcf_index = CNV.combined_cnv_calls_bed_vcf_index
    File? cnv_split_read_evidence = CNV.split_read_evidence
    File? cnv_split_read_evidence_index = CNV.split_read_evidence_index
    File? cnv_realign_read_evidence = CNV.realign_read_evidence
    File? cnv_realign_read_evidence_index = CNV.realign_read_evidence_index
    File? cnv_combine_read_scores_csv = CNV.combine_read_scores_csv
    File? cnv_combined_coverage_plot = CNV.combined_coverage_plot
    File? cnv_combined_dup_del_plot = CNV.combined_dup_del_plot
    File? cnv_combined_copy_number_plot = CNV.combined_copy_number_plot
    File? cnv_md5_checksums_json = CNV.md5_checksums_json
    File? sv_annotated_unlinked_vcf = SV.annotated_vcf_out
    File? sv_annotated_unlinked_vcf_index = SV.annotated_vcf_index_out
    File? sv_vcf = SV.output_vcf
    File? sv_vcf_index = SV.output_vcf_index
    File? sv_assembly = SV.assembly
    File? sv_assembly_index = SV.assembly_index
    File? sv_realigned_assembly = SV.realigned_assembly
    File? sv_realigned_assembly_index = SV.realigned_assembly_index
    File? sv_converted_vcf = SV.converted_vcf
    File? sv_converted_vcf_index = SV.converted_vcf_index
    File? sv_md5_checksums_json = SV.md5_checksums_json
    Array[File]? str_detailed_csv_files = STR.detailed_csv_files
    Array[File]? str_summary_csv_files = STR.summary_csv_files
    File? str_bed = STR.genotypes_bed
    File? str_vcf = STR.genotypes_vcf
    File? str_vcf_index = STR.genotypes_vcf_index
    File? hla_genotypes = HLA.output_hla
    File? kir_genotypes = HLA.output_kir
    Array[File]? pgx_allele_fraction_profiles = PGx.allele_fraction_profiles
    Array[File]? pgx_alleles = PGx.alleles
    Array[File]? pgx_cnv_calls = PGx.cnv_calls
    Array[File]? pgx_consolidated_variants = PGx.consolidated_variants
    Array[File]? pgx_copy_number_profiles = PGx.copy_number_profiles
    Array[File]? pgx_copy_numbers = PGx.copy_numbers
    Array[File]? pgx_genotypes = PGx.genotypes
    Array[File]? pgx_imported_variants = PGx.imported_variants
    Array[File]? pgx_phased_variants = PGx.phased_variants
    Array[File]? pgx_phenotypes = PGx.phenotypes
    Array[File]? pgx_read_depths = PGx.read_depths
    File? pgx_vcf = PGx.output_vcf
    File? pgx_vcf_index = PGx.output_vcf_index
    File? pgx_results = PGx.results
    File? segdup_remap_bam = SegDup.remap_bam
    File? segdup_remap_bam_index = SegDup.remap_bam_index
    File? segdup_acnv_calls = SegDup.acnv_calls
    File? segdup_pcnv_calls = SegDup.pcnv_calls
    File? segdup_small_variants = SegDup.small_variants
    File? segdup_small_variants_index = SegDup.small_variants_idx
    File? segdup_lpa_vcf = SegDup.lpa_vcf
    File? segdup_lpa_vcf_index = SegDup.lpa_vcf_index
    File? segdup_lpa_json = SegDup.lpa_json
  }
}
