# GermlinePipeline
Runs the Ultima Genomics germline analysis on a single aligned sample: SNV/indel calling with EfficientDV, CNV, structural variant, short tandem repeat, HLA, pharmacogenomics and segmental duplication analysis. Every stage is switched by its own run_* flag, all of which default to true.

## Inputs

### Required inputs
<p name="GermlinePipeline.max_num_haps">
        <b>GermlinePipeline.max_num_haps</b><br />
        <i>Int? </i> &mdash;
         Assembly parameter: Maximum number of haplotypes showing an evidence of SV to report <br />
</p>
<p name="GermlinePipeline.hla_genotyping_tool">
        <b>GermlinePipeline.hla_genotyping_tool</b><br />
        <i>String </i> &mdash;
         HLA genotyping tool to use. Options: 'HLA-LA' or 'T1K' <br />
</p>

### Required inputs
<p name="GermlinePipeline.input_cram_bam_list">
        <b>GermlinePipeline.input_cram_bam_list</b><br />
        <i>Array[File] </i> &mdash;
         Aligned, sorted, duplicate-marked CRAM/BAM file(s) of a single sample <br />
</p>
<p name="GermlinePipeline.input_cram_bam_index_list">
        <b>GermlinePipeline.input_cram_bam_index_list</b><br />
        <i>Array[File] </i> &mdash;
         CRAI/BAI index file(s) matching input_cram_bam_list <br />
</p>
<p name="GermlinePipeline.base_file_name">
        <b>GermlinePipeline.base_file_name</b><br />
        <i>String </i> &mdash;
         (shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Prefix for name of all output files <br />
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="GermlinePipeline.blacklist_bed">
        <b>GermlinePipeline.blacklist_bed</b><br />
        <i>File? &mdash; Default: None</i><br />
        Gridss blacklist file. When not provided, resolved per genome from genome_resources (sv_blacklist). Provide empty file to disable
</p>
<p name="GermlinePipeline.exclude_filters">
        <b>GermlinePipeline.exclude_filters</b><br />
        <i>String? &mdash; Default: None</i><br />
        gripss paramter: Exclude filters from the output vcf, separated by ;
</p>
<p name="GermlinePipeline.giraffe_parameters">
        <b>GermlinePipeline.giraffe_parameters</b><br />
        <i>GiraffeParameters? &mdash; Default: None</i><br />
        vg giraffe index files to improve haplotype interpretation using population graphs
</p>
<p name="GermlinePipeline.gridss_metrics_interval">
        <b>GermlinePipeline.gridss_metrics_interval</b><br />
        <i>String? &mdash; Default: None</i><br />
        Interval for collecting gridss metrics
</p>
<p name="GermlinePipeline.homopolymer_length">
        <b>GermlinePipeline.homopolymer_length</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Realignment parameter: do realignment on homopolymeres longer than this value
</p>
<p name="GermlinePipeline.input_tumor_crams">
        <b>GermlinePipeline.input_tumor_crams</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM file for the tumor (in case of matched T/N calling)
</p>
<p name="GermlinePipeline.input_tumor_crams_indexes">
        <b>GermlinePipeline.input_tumor_crams_indexes</b><br />
        <i>Array[File] &mdash; Default: []</i><br />
        Input CRAM index for the tumor (in case of matched T/N calling)
</p>
<p name="GermlinePipeline.known_hotspot_file">
        <b>GermlinePipeline.known_hotspot_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Known locations that are hot spot for SVs (see https://github.com/hartwigmedical/hmftools/tree/master/linx), filtered less stringently
</p>
<p name="GermlinePipeline.min_base">
        <b>GermlinePipeline.min_base</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Assembly parameter: Minimum base quality for using in DeBruijn graph construction. Default value in template
</p>
<p name="GermlinePipeline.min_indel_sc_size_to_include">
        <b>GermlinePipeline.min_indel_sc_size_to_include</b><br />
        <i>String? &mdash; Default: None</i><br />
        Assembly parameter: Minimum size of an indel and soft-clipping in the read to include the read in the assembly. ;-separated between samples
</p>
<p name="GermlinePipeline.min_mapq">
        <b>GermlinePipeline.min_mapq</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Assembly parameter: Minimum mapping quality. Default value in template
</p>
<p name="GermlinePipeline.min_mismatch_count_to_include">
        <b>GermlinePipeline.min_mismatch_count_to_include</b><br />
        <i>String? &mdash; Default: None</i><br />
        Assembly parameter: Minimal number of counts to require to include the read in the assembly. ;-separated between samples
</p>
<p name="GermlinePipeline.min_normal_coverage">
        <b>GermlinePipeline.min_normal_coverage</b><br />
        <i>Int? &mdash; Default: None</i><br />
        gripss paramter: Minimum coverage in the normal sample to determine somatic status. Default value:8
</p>
<p name="GermlinePipeline.pon_sgl_file">
        <b>GermlinePipeline.pon_sgl_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Panel of normals for single end breakend (partially resolved) calls. Note that the default value is in template
</p>
<p name="GermlinePipeline.pon_sv_file">
        <b>GermlinePipeline.pon_sv_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: panel of normals for breakpoint (fully resolved) calls. Note that the default value is in template
</p>
<p name="GermlinePipeline.prefilter_query">
        <b>GermlinePipeline.prefilter_query</b><br />
        <i>String? &mdash; Default: None</i><br />
        Expression (in bcftools view format) to filter the variants before annotation
</p>
<p name="GermlinePipeline.realign_mapq">
        <b>GermlinePipeline.realign_mapq</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Realignment parameter: Below this value we skip realignment on the supplementary alignment
</p>
<p name="GermlinePipeline.reference_name">
        <b>GermlinePipeline.reference_name</b><br />
        <i>String? &mdash; Default: None</i><br />
        Can be 38 or 19
</p>
<p name="GermlinePipeline.repeat_mask_file">
        <b>GermlinePipeline.repeat_mask_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        gripss paramter: Repeat mask file. Note that the default value is in template
</p>
<p name="GermlinePipeline.run_giraffe">
        <b>GermlinePipeline.run_giraffe</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Whether to run Giraffe haplotype aware alignment or not
</p>
<p name="GermlinePipeline.run_ua">
        <b>GermlinePipeline.run_ua</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Whether to run UA realignment on the output of the assembly (helps resolving some deletions) or not
</p>
<p name="GermlinePipeline.sv_num_shards">
        <b>GermlinePipeline.sv_num_shards</b><br />
        <i>Int? &mdash; Default: None</i><br />
        (SV) Relevant for scatter tasks, which are CreateAssembly and gridss.AnnotateVariants
</p>
<p name="GermlinePipeline.symbolic_vcf_format">
        <b>GermlinePipeline.symbolic_vcf_format</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Whether to convert the output vcf to the region format or not, default True
</p>
<p name="GermlinePipeline.ua_parameters">
        <b>GermlinePipeline.ua_parameters</b><br />
        <i>UaParameters? &mdash; Default: None</i><br />
        UA parameters: v_aware_alignment_flag and ua_extra_args, recommended value set in the template
</p>
<p name="GermlinePipeline.wgs_calling_interval_list_override">
        <b>GermlinePipeline.wgs_calling_interval_list_override</b><br />
        <i>File? &mdash; Default: None</i><br />
        Optional override for the interval list defining the region to perform variant calling on. When not provided, resolved per genome from genome_resources (calling_interval_list_without_artefacts)
</p>
<p name="GermlinePipeline.graphs_files_tar">
        <b>GermlinePipeline.graphs_files_tar</b><br />
        <i>File? &mdash; Default: None</i><br />
        HLA-LA graphs files tar (required if using HLA-LA)
</p>
<p name="GermlinePipeline.t1k_index_tar">
        <b>GermlinePipeline.t1k_index_tar</b><br />
        <i>File? &mdash; Default: None</i><br />
        T1K index tar.gz containing hlaidx/ and kiridx/ directories with all index files (_seq.fa and _coord.fa). Required if using T1K.
</p>

### Optional inputs
<p name="GermlinePipeline.reference_genome">
        <b>GermlinePipeline.reference_genome</b><br />
        <i>String </i> &mdash;
         (shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Genome selector: hg38, b37, hg38_taps, hg38_nist_v3, hg38_nist_v3_with_decoy, hg38_no_alt, mm10, mm39. Default to hg38 <br />
</p>
<p name="GermlinePipeline.monitoring_script_input">
        <b>GermlinePipeline.monitoring_script_input</b><br />
        <i>File? </i> &mdash;
         (shared: VariantCalling, CNV, SV, STR, HLA, PGx, SegDup) Monitoring script override for AWS HealthOmics workflow templates multi-region support <br />
</p>
<p name="GermlinePipeline.cloud_provider_override">
        <b>GermlinePipeline.cloud_provider_override</b><br />
        <i>String? </i> &mdash;
         (shared: VariantCalling, CNV, SV, HLA, PGx, SegDup) cloud_provider_override <br />
</p>
<p name="GermlinePipeline.background_cram_files">
        <b>GermlinePipeline.background_cram_files</b><br />
        <i>Array[File] </i> &mdash;
         Background (normal sample) cram files for somatic calling <br />
</p>
<p name="GermlinePipeline.background_cram_index_files">
        <b>GermlinePipeline.background_cram_index_files</b><br />
        <i>Array[File] </i> &mdash;
         Background (normal sample) cram index files for somatic calling <br />
</p>
<p name="GermlinePipeline.bed_graph">
        <b>GermlinePipeline.bed_graph</b><br />
        <i>Array[File]? </i> &mdash;
         Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data). <br />
</p>
<p name="GermlinePipeline.cnv_create_md5_checksum_outputs">
        <b>GermlinePipeline.cnv_create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         (CNV) Create md5 checksum for requested output files <br />
</p>
<p name="GermlinePipeline.cohort_reads_count_matrix_override">
        <b>GermlinePipeline.cohort_reads_count_matrix_override</b><br />
        <i>File? </i> &mdash;
         GenomicRanges object of the cohort reads count matrix in rds file format. By default the cohort matching reference_genome is taken from the genome resources. <br />
</p>
<p name="GermlinePipeline.filtering_model">
        <b>GermlinePipeline.filtering_model</b><br />
        <i>File? </i> &mdash;
         CNV filtering model, default in template, calls are not filtered if not provided <br />
</p>
<p name="GermlinePipeline.ploidy_file">
        <b>GermlinePipeline.ploidy_file</b><br />
        <i>File? </i> &mdash;
         X chromosome ploidy of the cohort and the additional sample. Each sample is represented on a number on a separate row. Ploidy of the default cohort can be found in the template. The last row corresponds to the sample being called. Genome independent. <br />
</p>
<p name="GermlinePipeline.sv_calls_vcf">
        <b>GermlinePipeline.sv_calls_vcf</b><br />
        <i>File? </i> &mdash;
         SV calls in VCF format (MANTA-like, single record per SV call) to be used for annotation of combined CNV calls, default is empty and annotation is not performed.<br> The input tested is the output of structrual_variant_pipeline.wdl <br />
</p>
<p name="GermlinePipeline.sv_calls_vcf_index">
        <b>GermlinePipeline.sv_calls_vcf_index</b><br />
        <i>File? </i> &mdash;
         Index file for the SV calls VCF <br />
</p>
<p name="GermlinePipeline.sv_create_md5_checksum_outputs">
        <b>GermlinePipeline.sv_create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         (SV) Create md5 checksum for requested output files <br />
</p>
<p name="GermlinePipeline.haploid">
        <b>GermlinePipeline.haploid</b><br />
        <i>Boolean </i> &mdash;
         Enable haploid mode: report single allele instead of diploid pairs. Use for X/Y chromosomes in males or haploid organisms. <br />
</p>
<p name="GermlinePipeline.max_repeat">
        <b>GermlinePipeline.max_repeat</b><br />
        <i>Int </i> &mdash;
         Maximum number of repeat units to include in auxiliary reference sequences <br />
</p>
<p name="GermlinePipeline.micro_allele_consensus_ratio">
        <b>GermlinePipeline.micro_allele_consensus_ratio</b><br />
        <i>Float </i> &mdash;
         Minimum fraction of supporting spanning reads required to report a micro-allele decimal. Only used when report_micro_alleles is true. Range 0.0-1.0. <br />
</p>
<p name="GermlinePipeline.micro_allele_min_reads">
        <b>GermlinePipeline.micro_allele_min_reads</b><br />
        <i>Int </i> &mdash;
         Minimum number of spanning reads supporting the consensus tract length required to report a micro-allele. Guards against low-coverage indel artifacts. Only used when report_micro_alleles is true. Default: 10. <br />
</p>
<p name="GermlinePipeline.min_repeat">
        <b>GermlinePipeline.min_repeat</b><br />
        <i>Int </i> &mdash;
         Minimum number of repeat units to include in auxiliary reference sequences <br />
</p>
<p name="GermlinePipeline.min_score_ratio">
        <b>GermlinePipeline.min_score_ratio</b><br />
        <i>Float </i> &mdash;
         Minimum ratio of alignment score to the theoretical maximum score (read_length * match_score). Alignments below this threshold are filtered out. Range: 0.0-1.0, where 1.0 requires perfect alignment. <br />
</p>
<p name="GermlinePipeline.output_summary_csv">
        <b>GermlinePipeline.output_summary_csv</b><br />
        <i>Boolean </i> &mdash;
         Whether to output summary per-locus CSV file. Set to false for large catalogs to reduce I/O. <br />
</p>
<p name="GermlinePipeline.ref_padding">
        <b>GermlinePipeline.ref_padding</b><br />
        <i>Int </i> &mdash;
         Number of bases to extend around the STR repeat region when building auxiliary references for alignment. Larger values provide more flanking sequence context for accurate alignment. <br />
</p>
<p name="GermlinePipeline.report_micro_alleles">
        <b>GermlinePipeline.report_micro_alleles</b><br />
        <i>Boolean </i> &mdash;
         Report micro-alleles (e.g. 15.3) for haploid loci when a partial-repeat insertion is present. REPCN stays the integer floor; the micro-allele decimal appears in the genotype (GT) and a new RCMA VCF/BED field. No-op for diploid loci. Default: off. <br />
</p>
<p name="GermlinePipeline.spanning_flank_bases">
        <b>GermlinePipeline.spanning_flank_bases</b><br />
        <i>Int </i> &mdash;
         Minimum number of bases that must align on each side of the STR repeat region for a read to be considered 'spanning' the locus <br />
</p>
<p name="GermlinePipeline.str_min_mapping_quality">
        <b>GermlinePipeline.str_min_mapping_quality</b><br />
        <i>Int </i> &mdash;
         (STR) Minimum mapping quality for reads to be included in analysis <br />
</p>
<p name="GermlinePipeline.variant_catalog">
        <b>GermlinePipeline.variant_catalog</b><br />
        <i>File? </i> &mdash;
         Variant catalog (json). Example: https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_with_offtargets.GRCh38.json <br />
</p>
<p name="GermlinePipeline.gene_symbols">
        <b>GermlinePipeline.gene_symbols</b><br />
        <i>Array[String]? </i> &mdash;
         List of gene symbols to analyze <br />
</p>
<p name="GermlinePipeline.input_vcf_file">
        <b>GermlinePipeline.input_vcf_file</b><br />
        <i>File? </i> &mdash;
         Input VCF file with variants. Use of high quality variants (i.e. PASS). If not provided, Efficient DV will be run <br />
</p>
<p name="GermlinePipeline.input_vcf_index_file">
        <b>GermlinePipeline.input_vcf_index_file</b><br />
        <i>File? </i> &mdash;
         Input VCF index file <br />
</p>
<p name="GermlinePipeline.pgx_model_onnx">
        <b>GermlinePipeline.pgx_model_onnx</b><br />
        <i>File? </i> &mdash;
         (PGx) TensorRT model for calling variants (onnx format) <br />
</p>
<p name="GermlinePipeline.ref_files_for_tarball">
        <b>GermlinePipeline.ref_files_for_tarball</b><br />
        <i>Array[File]? </i> &mdash;
         List of references for CreateReferenceCache task. <br />
</p>
<p name="GermlinePipeline.n_threads">
        <b>GermlinePipeline.n_threads</b><br />
        <i>Int? </i> &mdash;
         Number of threads to use <br />
</p>
<p name="GermlinePipeline.segdup_preemptible_tries">
        <b>GermlinePipeline.segdup_preemptible_tries</b><br />
        <i>Int </i> &mdash;
         (SegDup) Number of preemptible tries <br />
</p>
<p name="GermlinePipeline.VariantCalling.ScatterIntervalList.convert_to_bed">
        <b>GermlinePipeline.VariantCalling.ScatterIntervalList.convert_to_bed</b><br />
        <i>Boolean? </i> &mdash;
         If true, convert interval_list files to BED format in addition to interval_list format <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.input_bam_file">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.input_bam_file</b><br />
        <i>File? </i> &mdash;
         Input sample BAM/CRAM file. one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.input_bam_file_index">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.input_bam_file_index</b><br />
        <i>File? </i> &mdash;
         Input sample BAI/CRAI index file <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.input_sample_reads_count">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.input_sample_reads_count</b><br />
        <i>File? </i> &mdash;
         Inputs sample windowed coverage stored as GenomicRanges object in rds file. can be calculated using cn.mops::getReadCountsFromBAM R function.  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.SingleSampleReadsCount.ref_seq_names">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.SingleSampleReadsCount.ref_seq_names</b><br />
        <i>Array[String]? </i> &mdash;
         Chromosome names for which reads counts will be calculated. Mutually exclusive with genome_windows <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.SingleSampleReadsCount.window_length">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.SingleSampleReadsCount.window_length</b><br />
        <i>Int? </i> &mdash;
         Window lenght for which reads counts will be calculated for. Mutually exclusive with genome_windows <br />
</p>
<p name="GermlinePipeline.SV.ScatterIntervalList.convert_to_bed">
        <b>GermlinePipeline.SV.ScatterIntervalList.convert_to_bed</b><br />
        <i>Boolean? </i> &mdash;
         If true, convert interval_list files to BED format in addition to interval_list format <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.ScatterIntervalList.convert_to_bed">
        <b>GermlinePipeline.PGx.EfficientDV.ScatterIntervalList.convert_to_bed</b><br />
        <i>Boolean? </i> &mdash;
         If true, convert interval_list files to BED format in addition to interval_list format <br />
</p>
<p name="GermlinePipeline.SegDup.DV.background_cram_files">
        <b>GermlinePipeline.SegDup.DV.background_cram_files</b><br />
        <i>Array[File] </i> &mdash;
         Background (normal sample) cram files for somatic calling <br />
</p>
<p name="GermlinePipeline.SegDup.DV.background_cram_index_files">
        <b>GermlinePipeline.SegDup.DV.background_cram_index_files</b><br />
        <i>Array[File] </i> &mdash;
         Background (normal sample) cram index files for somatic calling <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ScatterIntervalList.convert_to_bed">
        <b>GermlinePipeline.SegDup.DV.ScatterIntervalList.convert_to_bed</b><br />
        <i>Boolean? </i> &mdash;
         If true, convert interval_list files to BED format in addition to interval_list format <br />
</p>

### Optional parameters
<p name="GermlinePipeline.run_variant_calling">
        <b>GermlinePipeline.run_variant_calling</b><br />
        <i>Boolean </i> &mdash;
         Run SNV/indel variant calling (EfficientDV). Needs a CRAM with real coverage: EfficientDV aborts when the measured median coverage is 0 <br />
</p>
<p name="GermlinePipeline.run_cnv">
        <b>GermlinePipeline.run_cnv</b><br />
        <i>Boolean </i> &mdash;
         Run germline CNV calling (cn.MOPS + CNVpytor). Requires bed_graph and ploidy_file <br />
</p>
<p name="GermlinePipeline.run_sv">
        <b>GermlinePipeline.run_sv</b><br />
        <i>Boolean </i> &mdash;
         Run structural variant calling <br />
</p>
<p name="GermlinePipeline.run_str">
        <b>GermlinePipeline.run_str</b><br />
        <i>Boolean </i> &mdash;
         Run short tandem repeat genotyping <br />
</p>
<p name="GermlinePipeline.run_hla">
        <b>GermlinePipeline.run_hla</b><br />
        <i>Boolean </i> &mdash;
         Run HLA/KIR genotyping <br />
</p>
<p name="GermlinePipeline.run_pgx">
        <b>GermlinePipeline.run_pgx</b><br />
        <i>Boolean </i> &mdash;
         Run pharmacogenomics (PyPGx) analysis <br />
</p>
<p name="GermlinePipeline.run_segdup">
        <b>GermlinePipeline.run_segdup</b><br />
        <i>Boolean </i> &mdash;
         Run segmental duplication analysis (parascopy/LPA) <br />
</p>
<p name="GermlinePipeline.active_areas_min_base_quality">
        <b>GermlinePipeline.active_areas_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimum base quality for active areas detection <br />
</p>
<p name="GermlinePipeline.allele_frequency_ratio">
        <b>GermlinePipeline.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering <br />
</p>
<p name="GermlinePipeline.call_variants_cpus">
        <b>GermlinePipeline.call_variants_cpus</b><br />
        <i>Int </i> &mdash;
         Number of CPUs for call_variants <br />
</p>
<p name="GermlinePipeline.call_variants_gpu_type_override">
        <b>GermlinePipeline.call_variants_gpu_type_override</b><br />
        <i>String? </i> &mdash;
         GPU type for call variants <br />
</p>
<p name="GermlinePipeline.call_variants_gpus">
        <b>GermlinePipeline.call_variants_gpus</b><br />
        <i>Int </i> &mdash;
         Number of GPUs for call_variants <br />
</p>
<p name="GermlinePipeline.call_variants_threads">
        <b>GermlinePipeline.call_variants_threads</b><br />
        <i>Int </i> &mdash;
         Number of decompression threads for call_variants <br />
</p>
<p name="GermlinePipeline.call_variants_uncompr_buf_size_gb">
        <b>GermlinePipeline.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash;
         Memory buffer allocated for each uncompression thread in calll_variants <br />
</p>
<p name="GermlinePipeline.dbg_min_base_quality">
        <b>GermlinePipeline.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for local assembly of haplotypes <br />
</p>
<p name="GermlinePipeline.dv_diploid_sampling_in_haplotypes">
        <b>GermlinePipeline.dv_diploid_sampling_in_haplotypes</b><br />
        <i>Boolean? </i> &mdash;
         (VariantCalling) Use diploid sampling strategy for haplotype selection <br />
</p>
<p name="GermlinePipeline.dv_ensemble_size">
        <b>GermlinePipeline.dv_ensemble_size</b><br />
        <i>Int </i> &mdash;
         (VariantCalling) Number of augmented passes for ensemble inference. Values <= 1 disable ensemble entirely (no augmentation is applied); values >= 2 enable selective ensemble. <br />
</p>
<p name="GermlinePipeline.dv_include_reference_in_haplotypes">
        <b>GermlinePipeline.dv_include_reference_in_haplotypes</b><br />
        <i>Boolean? </i> &mdash;
         (VariantCalling) Include the reference sequence in the sampled haplotypes <br />
</p>
<p name="GermlinePipeline.dv_max_reads_per_partition">
        <b>GermlinePipeline.dv_max_reads_per_partition</b><br />
        <i>Int </i> &mdash;
         (VariantCalling) Maximal number of reads that are stored in memory when analyzing an active region <br />
</p>
<p name="GermlinePipeline.dv_min_base_quality">
        <b>GermlinePipeline.dv_min_base_quality</b><br />
        <i>Int </i> &mdash;
         (VariantCalling) Minimal base quality for candidate generation <br />
</p>
<p name="GermlinePipeline.dv_min_mapping_quality">
        <b>GermlinePipeline.dv_min_mapping_quality</b><br />
        <i>Int </i> &mdash;
         (VariantCalling) Minimum mapping quality for reads to appear in pileup images (input to CNN) and to be considered as supporting an alt-allele in candidate generation <br />
</p>
<p name="GermlinePipeline.dv_num_haplotypes">
        <b>GermlinePipeline.dv_num_haplotypes</b><br />
        <i>Int? </i> &mdash;
         (VariantCalling) Number of haplotypes in the pangenome haplotype CRAM. Also determines the haplotype band height in the pileup image. <br />
</p>
<p name="GermlinePipeline.dv_run_haplotype_sampling">
        <b>GermlinePipeline.dv_run_haplotype_sampling</b><br />
        <i>Boolean </i> &mdash;
         (VariantCalling) Whether to run haplotype sampling to create pangenome haplotypes. Default: false <br />
</p>
<p name="GermlinePipeline.dv_scatter_intervals_break">
        <b>GermlinePipeline.dv_scatter_intervals_break</b><br />
        <i>Int </i> &mdash;
         (VariantCalling) The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br />
</p>
<p name="GermlinePipeline.dv_ug_make_examples_extra_args">
        <b>GermlinePipeline.dv_ug_make_examples_extra_args</b><br />
        <i>String? </i> &mdash;
         (VariantCalling) Additional arguments for make-examples tool <br />
</p>
<p name="GermlinePipeline.ensemble_reference_rows">
        <b>GermlinePipeline.ensemble_reference_rows</b><br />
        <i>Int </i> &mdash;
         Number of reference rows for ensemble inference <br />
</p>
<p name="GermlinePipeline.germline_vcf">
        <b>GermlinePipeline.germline_vcf</b><br />
        <i>File? </i> &mdash;
         Germline vcf file in order to generate haplotypes that incorporate germline variants <br />
</p>
<p name="GermlinePipeline.gq_bins">
        <b>GermlinePipeline.gq_bins</b><br />
        <i>Array[Int]? </i> &mdash;
         GQ bins to use instead of a fixed resolution (overrides gq_resolution) <br />
</p>
<p name="GermlinePipeline.gq_resolution_override">
        <b>GermlinePipeline.gq_resolution_override</b><br />
        <i>Int? </i> &mdash;
         Override for gq resolution (default: 5) <br />
</p>
<p name="GermlinePipeline.h_indel_allele_frequency_ratio">
        <b>GermlinePipeline.h_indel_allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering <br />
</p>
<p name="GermlinePipeline.h_indel_vaf_to_pass">
        <b>GermlinePipeline.h_indel_vaf_to_pass</b><br />
        <i>Float? </i> &mdash;
         Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio <br />
</p>
<p name="GermlinePipeline.hard_qual_filter">
        <b>GermlinePipeline.hard_qual_filter</b><br />
        <i>Int </i> &mdash;
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br />
</p>
<p name="GermlinePipeline.input_flow_order">
        <b>GermlinePipeline.input_flow_order</b><br />
        <i>String? </i> &mdash;
         Flow order. If not provided, it will be extracted from the CRAM header <br />
</p>
<p name="GermlinePipeline.intervals_string">
        <b>GermlinePipeline.intervals_string</b><br />
        <i>String? </i> &mdash;
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. Takes precedence over override_target_intervals. <br />
</p>
<p name="GermlinePipeline.log_make_examples_progress">
        <b>GermlinePipeline.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash;
         Cause make_examples to output detailed progress information (for debugging) <br />
</p>
<p name="GermlinePipeline.make_gvcf">
        <b>GermlinePipeline.make_gvcf</b><br />
        <i>Boolean? </i> &mdash;
         Whether to generate a gvcf. Default: False <br />
</p>
<p name="GermlinePipeline.min_fraction_hmer_indels">
        <b>GermlinePipeline.min_fraction_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_fraction_non_hmer_indels">
        <b>GermlinePipeline.min_fraction_non_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_fraction_snps">
        <b>GermlinePipeline.min_fraction_snps</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a snp, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_hmer_plus_one_candidate">
        <b>GermlinePipeline.min_hmer_plus_one_candidate</b><br />
        <i>Int </i> &mdash;
         Minimal hmer length, above which more 1-bp insertion candidates are generated, provided they also meet allele frequency conditions <br />
</p>
<p name="GermlinePipeline.min_read_count_hmer_indels">
        <b>GermlinePipeline.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_read_count_non_hmer_indels">
        <b>GermlinePipeline.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_read_count_snps">
        <b>GermlinePipeline.min_read_count_snps</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a snp, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.min_variant_quality_exome_hmer_indels">
        <b>GermlinePipeline.min_variant_quality_exome_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal non-h-mer indel quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.min_variant_quality_hmer_indels">
        <b>GermlinePipeline.min_variant_quality_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal h-mer indel quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.min_variant_quality_non_hmer_indels">
        <b>GermlinePipeline.min_variant_quality_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal non-h-mer indel quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.min_variant_quality_snps">
        <b>GermlinePipeline.min_variant_quality_snps</b><br />
        <i>Int </i> &mdash;
         Minimal snp variant quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.optimization_level">
        <b>GermlinePipeline.optimization_level</b><br />
        <i>Int? </i> &mdash;
         Optimization level for TensorRT engine in call_variants <br />
</p>
<p name="GermlinePipeline.output_call_variants_tfrecords">
        <b>GermlinePipeline.output_call_variants_tfrecords</b><br />
        <i>Boolean </i> &mdash;
         Output tfrecords from call_variants <br />
</p>
<p name="GermlinePipeline.output_realignment">
        <b>GermlinePipeline.output_realignment</b><br />
        <i>Boolean </i> &mdash;
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br />
</p>
<p name="GermlinePipeline.override_target_intervals">
        <b>GermlinePipeline.override_target_intervals</b><br />
        <i>File? </i> &mdash;
         Override default genome-specific target intervals. If not provided, uses genome-specific default intervals. <br />
</p>
<p name="GermlinePipeline.p_error">
        <b>GermlinePipeline.p_error</b><br />
        <i>Float </i> &mdash;
         Basecalling error for reference confidence model in gvcf <br />
</p>
<p name="GermlinePipeline.pangenome_haplotypes">
        <b>GermlinePipeline.pangenome_haplotypes</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram file <br />
</p>
<p name="GermlinePipeline.pangenome_haplotypes_index">
        <b>GermlinePipeline.pangenome_haplotypes_index</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram index file <br />
</p>
<p name="GermlinePipeline.prioritize_alt_supporting_reads">
        <b>GermlinePipeline.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash;
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br />
</p>
<p name="GermlinePipeline.prioritize_high_quality_reads">
        <b>GermlinePipeline.prioritize_high_quality_reads</b><br />
        <i>Boolean </i> &mdash;
         When min-mapq=0, add mapq=0 reads last, only filling remaining image capacity after high-mapq reads <br />
</p>
<p name="GermlinePipeline.random_seed">
        <b>GermlinePipeline.random_seed</b><br />
        <i>Int </i> &mdash;
         Random seed for ensemble inference <br />
</p>
<p name="GermlinePipeline.roh_af_default">
        <b>GermlinePipeline.roh_af_default</b><br />
        <i>Float </i> &mdash;
         Alternate allele frequency assumed for every marker by the ROH caller, in place of a population frequency table <br />
</p>
<p name="GermlinePipeline.run_ploidy_estimation">
        <b>GermlinePipeline.run_ploidy_estimation</b><br />
        <i>Boolean </i> &mdash;
         Run VCF-based ploidy estimation and chrX/Y haploid conversion for germline samples. Default: false; enabled by germline use cases. <br />
</p>
<p name="GermlinePipeline.run_roh">
        <b>GermlinePipeline.run_roh</b><br />
        <i>Boolean </i> &mdash;
         Whether to call runs of homozygosity (ROH). Enabled by default in the germline WGS use-cases, off otherwise. Requires a reference genome that has a roh_blacklist resource (the hg38 builds and b37) unless roh_blacklist_override is given <br />
</p>
<p name="GermlinePipeline.sex_chromosomes">
        <b>GermlinePipeline.sex_chromosomes</b><br />
        <i>Array[String] </i> &mdash;
         Sex chromosome names to exclude from autosomal ploidy baseline. Defaults support chr-prefixed and non-prefixed human references. <br />
</p>
<p name="GermlinePipeline.shuffle_all_samples">
        <b>GermlinePipeline.shuffle_all_samples</b><br />
        <i>Boolean </i> &mdash;
         Whether to shuffle all samples during inference <br />
</p>
<p name="GermlinePipeline.strong_call_threshold">
        <b>GermlinePipeline.strong_call_threshold</b><br />
        <i>Float </i> &mdash;
         Probability threshold for selective ensemble inference. When ensemble_size >= 2, examples with max probability below this threshold are re-evaluated using ensemble inference; examples above it are accepted as-is. <br />
</p>
<p name="GermlinePipeline.trim_soft_clips">
        <b>GermlinePipeline.trim_soft_clips</b><br />
        <i>Boolean </i> &mdash;
         Trim soft-clipped bases from pileup images <br />
</p>
<p name="GermlinePipeline.ug_post_processing_extra_args">
        <b>GermlinePipeline.ug_post_processing_extra_args</b><br />
        <i>String </i> &mdash;
         Additional arguments for post-processing <br />
</p>
<p name="GermlinePipeline.v_gpu_tile_size">
        <b>GermlinePipeline.v_gpu_tile_size</b><br />
        <i>Int </i> &mdash;
         Virtual GPU tile size for call_variants <br />
</p>
<p name="GermlinePipeline.cushion_size">
        <b>GermlinePipeline.cushion_size</b><br />
        <i>Int? </i> &mdash;
         Cushion size around CNV breakpoints for split-read analysis and jump alignment analysis <br />
</p>
<p name="GermlinePipeline.filtering_model_decision_threshold">
        <b>GermlinePipeline.filtering_model_decision_threshold</b><br />
        <i>Int? </i> &mdash;
         Decision threshold for the filtering model, default is set in template. Lower- less stringent, Higher- more stringent <br />
</p>
<p name="GermlinePipeline.skip_figure_generation">
        <b>GermlinePipeline.skip_figure_generation</b><br />
        <i>Boolean? </i> &mdash;
         Skip CNV calls figure generation. Default is: False <br />
</p>
<p name="GermlinePipeline.skip_filtering">
        <b>GermlinePipeline.skip_filtering</b><br />
        <i>Boolean? </i> &mdash;
         Whether to skip CNV filtering step, default is False <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.kmer_length">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.kmer_length</b><br />
        <i>Int </i> &mdash;
         K-mer length for KMC counting (default: 29) <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.min_kmer_count">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.min_kmer_count</b><br />
        <i>Int </i> &mdash;
         Minimum k-mer count threshold for sampling (default: 2) <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.window_size">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.window_size</b><br />
        <i>Int </i> &mdash;
         Sliding window size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.step_size">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.step_size</b><br />
        <i>Int </i> &mdash;
         Sliding window step size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.minimap2_preset">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.minimap2_preset</b><br />
        <i>String </i> &mdash;
         Minimap2 preset for alignment (default: asm5) <br />
</p>
<p name="GermlinePipeline.VariantCalling.HaplotypeSampling.minimap_extra_args">
        <b>GermlinePipeline.VariantCalling.HaplotypeSampling.minimap_extra_args</b><br />
        <i>String? </i> &mdash;
         Additional extra arguments to pass to minimap2 (default: empty) <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.cap_coverage_override">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.cap_coverage_override</b><br />
        <i>Boolean? </i> &mdash;
         whether to cap extremely high average coverage windows to 2*cohort's average coverage quantile 99.9% value <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.save_hdf_override">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.save_hdf_override</b><br />
        <i>Boolean? </i> &mdash;
         Whether to save sample reads counts/cohort including sample/cnmops output data in hdf5 format (additionally to RDS format). Default is: False. <br />
</p>
<p name="GermlinePipeline.CNV.CnmopsCNVCalling.save_csv_override">
        <b>GermlinePipeline.CNV.CnmopsCNVCalling.save_csv_override</b><br />
        <i>Boolean? </i> &mdash;
         Whether to save sample reads counts/cohort including sample/cnmops output data in csv format (additionally to RDS format). Default is: False. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.show_bg_fields">
        <b>GermlinePipeline.PGx.EfficientDV.show_bg_fields</b><br />
        <i>Boolean </i> &mdash;
         Show background fields in the output vcf. Default: false. Mostly relevant for somatic calling. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.scatter_intervals_break">
        <b>GermlinePipeline.PGx.EfficientDV.scatter_intervals_break</b><br />
        <i>Int </i> &mdash;
         The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.intervals_string">
        <b>GermlinePipeline.PGx.EfficientDV.intervals_string</b><br />
        <i>String? </i> &mdash;
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. Takes precedence over override_target_intervals. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_fraction_hmer_indels">
        <b>GermlinePipeline.PGx.EfficientDV.min_fraction_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_fraction_non_hmer_indels">
        <b>GermlinePipeline.PGx.EfficientDV.min_fraction_non_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_fraction_snps">
        <b>GermlinePipeline.PGx.EfficientDV.min_fraction_snps</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a snp, required to  generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_read_count_snps">
        <b>GermlinePipeline.PGx.EfficientDV.min_read_count_snps</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a snp, required to  generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_read_count_hmer_indels">
        <b>GermlinePipeline.PGx.EfficientDV.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_read_count_non_hmer_indels">
        <b>GermlinePipeline.PGx.EfficientDV.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.min_hmer_plus_one_candidate">
        <b>GermlinePipeline.PGx.EfficientDV.min_hmer_plus_one_candidate</b><br />
        <i>Int </i> &mdash;
         Minimal hmer length, above which more 1-bp insertion candidates are generated, provided they also meet allele frequency conditions <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.max_reads_per_partition">
        <b>GermlinePipeline.PGx.EfficientDV.max_reads_per_partition</b><br />
        <i>Int </i> &mdash;
         Maximal number of reads that are stored in memory when analyzing an active region <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.dbg_min_base_quality">
        <b>GermlinePipeline.PGx.EfficientDV.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for local assembly of haplotypes <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.prioritize_alt_supporting_reads">
        <b>GermlinePipeline.PGx.EfficientDV.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash;
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.active_areas_min_base_quality">
        <b>GermlinePipeline.PGx.EfficientDV.active_areas_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimum base quality for active areas detection <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.trim_soft_clips">
        <b>GermlinePipeline.PGx.EfficientDV.trim_soft_clips</b><br />
        <i>Boolean </i> &mdash;
         Trim soft-clipped bases from pileup images <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.p_error">
        <b>GermlinePipeline.PGx.EfficientDV.p_error</b><br />
        <i>Float </i> &mdash;
         Basecalling error for reference confidence model in gvcf <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.gq_resolution_override">
        <b>GermlinePipeline.PGx.EfficientDV.gq_resolution_override</b><br />
        <i>Int? </i> &mdash;
         Override for gq resolution (default: 5) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.gq_bins">
        <b>GermlinePipeline.PGx.EfficientDV.gq_bins</b><br />
        <i>Array[Int]? </i> &mdash;
         GQ bins to use instead of a fixed resolution (overrides gq_resolution) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.output_realignment">
        <b>GermlinePipeline.PGx.EfficientDV.output_realignment</b><br />
        <i>Boolean </i> &mdash;
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.log_make_examples_progress">
        <b>GermlinePipeline.PGx.EfficientDV.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash;
         Cause make_examples to output detailed progress information (for debugging) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.germline_vcf">
        <b>GermlinePipeline.PGx.EfficientDV.germline_vcf</b><br />
        <i>File? </i> &mdash;
         Germline vcf file in order to generate haplotypes that incorporate germline variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.pangenome_haplotypes">
        <b>GermlinePipeline.PGx.EfficientDV.pangenome_haplotypes</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram file <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.pangenome_haplotypes_index">
        <b>GermlinePipeline.PGx.EfficientDV.pangenome_haplotypes_index</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram index file <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.optimization_level">
        <b>GermlinePipeline.PGx.EfficientDV.optimization_level</b><br />
        <i>Int? </i> &mdash;
         Optimization level for TensorRT engine in call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.output_call_variants_tfrecords">
        <b>GermlinePipeline.PGx.EfficientDV.output_call_variants_tfrecords</b><br />
        <i>Boolean </i> &mdash;
         Output tfrecords from call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.run_ploidy_estimation">
        <b>GermlinePipeline.PGx.EfficientDV.run_ploidy_estimation</b><br />
        <i>Boolean </i> &mdash;
         Run VCF-based ploidy estimation and chrX/Y haploid conversion for germline samples. Default: false; enabled by germline use cases. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.sex_chromosomes">
        <b>GermlinePipeline.PGx.EfficientDV.sex_chromosomes</b><br />
        <i>Array[String] </i> &mdash;
         Sex chromosome names to exclude from autosomal ploidy baseline. Defaults support chr-prefixed and non-prefixed human references. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.strong_call_threshold">
        <b>GermlinePipeline.PGx.EfficientDV.strong_call_threshold</b><br />
        <i>Float </i> &mdash;
         Probability threshold for selective ensemble inference. When ensemble_size >= 2, examples with max probability below this threshold are re-evaluated using ensemble inference; examples above it are accepted as-is. <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.ensemble_reference_rows">
        <b>GermlinePipeline.PGx.EfficientDV.ensemble_reference_rows</b><br />
        <i>Int </i> &mdash;
         Number of reference rows for ensemble inference <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.random_seed">
        <b>GermlinePipeline.PGx.EfficientDV.random_seed</b><br />
        <i>Int </i> &mdash;
         Random seed for ensemble inference <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.shuffle_all_samples">
        <b>GermlinePipeline.PGx.EfficientDV.shuffle_all_samples</b><br />
        <i>Boolean </i> &mdash;
         Whether to shuffle all samples during inference <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.hard_qual_filter">
        <b>GermlinePipeline.PGx.EfficientDV.hard_qual_filter</b><br />
        <i>Int </i> &mdash;
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.allele_frequency_ratio">
        <b>GermlinePipeline.PGx.EfficientDV.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.h_indel_vaf_to_pass">
        <b>GermlinePipeline.PGx.EfficientDV.h_indel_vaf_to_pass</b><br />
        <i>Float? </i> &mdash;
         Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.h_indel_allele_frequency_ratio">
        <b>GermlinePipeline.PGx.EfficientDV.h_indel_allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.ug_post_processing_extra_args">
        <b>GermlinePipeline.PGx.EfficientDV.ug_post_processing_extra_args</b><br />
        <i>String </i> &mdash;
         Additional arguments for post-processing <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.run_roh">
        <b>GermlinePipeline.PGx.EfficientDV.run_roh</b><br />
        <i>Boolean </i> &mdash;
         Whether to call runs of homozygosity (ROH). Enabled by default in the germline WGS use-cases, off otherwise. Requires a reference genome that has a roh_blacklist resource (the hg38 builds and b37) unless roh_blacklist_override is given <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.roh_af_default">
        <b>GermlinePipeline.PGx.EfficientDV.roh_af_default</b><br />
        <i>Float </i> &mdash;
         Alternate allele frequency assumed for every marker by the ROH caller, in place of a population frequency table <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.input_flow_order">
        <b>GermlinePipeline.PGx.EfficientDV.input_flow_order</b><br />
        <i>String? </i> &mdash;
         Flow order. If not provided, it will be extracted from the CRAM header <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.call_variants_gpu_type_override">
        <b>GermlinePipeline.PGx.EfficientDV.call_variants_gpu_type_override</b><br />
        <i>String? </i> &mdash;
         GPU type for call variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.call_variants_gpus">
        <b>GermlinePipeline.PGx.EfficientDV.call_variants_gpus</b><br />
        <i>Int </i> &mdash;
         Number of GPUs for call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.call_variants_cpus">
        <b>GermlinePipeline.PGx.EfficientDV.call_variants_cpus</b><br />
        <i>Int </i> &mdash;
         Number of CPUs for call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.call_variants_threads">
        <b>GermlinePipeline.PGx.EfficientDV.call_variants_threads</b><br />
        <i>Int </i> &mdash;
         Number of decompression threads for call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.call_variants_uncompr_buf_size_gb">
        <b>GermlinePipeline.PGx.EfficientDV.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash;
         Memory buffer allocated for each uncompression thread in calll_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.v_gpu_tile_size">
        <b>GermlinePipeline.PGx.EfficientDV.v_gpu_tile_size</b><br />
        <i>Int </i> &mdash;
         Virtual GPU tile size for call_variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.kmer_length">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.kmer_length</b><br />
        <i>Int </i> &mdash;
         K-mer length for KMC counting (default: 29) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.min_kmer_count">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.min_kmer_count</b><br />
        <i>Int </i> &mdash;
         Minimum k-mer count threshold for sampling (default: 2) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.window_size">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.window_size</b><br />
        <i>Int </i> &mdash;
         Sliding window size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.step_size">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.step_size</b><br />
        <i>Int </i> &mdash;
         Sliding window step size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.minimap2_preset">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.minimap2_preset</b><br />
        <i>String </i> &mdash;
         Minimap2 preset for alignment (default: asm5) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.minimap_extra_args">
        <b>GermlinePipeline.PGx.EfficientDV.HaplotypeSampling.minimap_extra_args</b><br />
        <i>String? </i> &mdash;
         Additional extra arguments to pass to minimap2 (default: empty) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.show_bg_fields">
        <b>GermlinePipeline.SegDup.DV.show_bg_fields</b><br />
        <i>Boolean </i> &mdash;
         Show background fields in the output vcf. Default: false. Mostly relevant for somatic calling. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.run_haplotype_sampling">
        <b>GermlinePipeline.SegDup.DV.run_haplotype_sampling</b><br />
        <i>Boolean </i> &mdash;
         Whether to run haplotype sampling to create pangenome haplotypes. Default: false <br />
</p>
<p name="GermlinePipeline.SegDup.DV.scatter_intervals_break">
        <b>GermlinePipeline.SegDup.DV.scatter_intervals_break</b><br />
        <i>Int </i> &mdash;
         The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.intervals_string">
        <b>GermlinePipeline.SegDup.DV.intervals_string</b><br />
        <i>String? </i> &mdash;
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. Takes precedence over override_target_intervals. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_read_count_snps">
        <b>GermlinePipeline.SegDup.DV.min_read_count_snps</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a snp, required to  generate a candidate variant <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_read_count_hmer_indels">
        <b>GermlinePipeline.SegDup.DV.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_read_count_non_hmer_indels">
        <b>GermlinePipeline.SegDup.DV.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_base_quality">
        <b>GermlinePipeline.SegDup.DV.min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for candidate generation <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_mapping_quality">
        <b>GermlinePipeline.SegDup.DV.min_mapping_quality</b><br />
        <i>Int </i> &mdash;
         Minimum mapping quality for reads to appear in pileup images (input to CNN) and to be considered as supporting an alt-allele in candidate generation <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_hmer_plus_one_candidate">
        <b>GermlinePipeline.SegDup.DV.min_hmer_plus_one_candidate</b><br />
        <i>Int </i> &mdash;
         Minimal hmer length, above which more 1-bp insertion candidates are generated, provided they also meet allele frequency conditions <br />
</p>
<p name="GermlinePipeline.SegDup.DV.max_reads_per_partition">
        <b>GermlinePipeline.SegDup.DV.max_reads_per_partition</b><br />
        <i>Int </i> &mdash;
         Maximal number of reads that are stored in memory when analyzing an active region <br />
</p>
<p name="GermlinePipeline.SegDup.DV.dbg_min_base_quality">
        <b>GermlinePipeline.SegDup.DV.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for local assembly of haplotypes <br />
</p>
<p name="GermlinePipeline.SegDup.DV.prioritize_alt_supporting_reads">
        <b>GermlinePipeline.SegDup.DV.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash;
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br />
</p>
<p name="GermlinePipeline.SegDup.DV.active_areas_min_base_quality">
        <b>GermlinePipeline.SegDup.DV.active_areas_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimum base quality for active areas detection <br />
</p>
<p name="GermlinePipeline.SegDup.DV.prioritize_high_quality_reads">
        <b>GermlinePipeline.SegDup.DV.prioritize_high_quality_reads</b><br />
        <i>Boolean </i> &mdash;
         When min-mapq=0, add mapq=0 reads last, only filling remaining image capacity after high-mapq reads <br />
</p>
<p name="GermlinePipeline.SegDup.DV.trim_soft_clips">
        <b>GermlinePipeline.SegDup.DV.trim_soft_clips</b><br />
        <i>Boolean </i> &mdash;
         Trim soft-clipped bases from pileup images <br />
</p>
<p name="GermlinePipeline.SegDup.DV.p_error">
        <b>GermlinePipeline.SegDup.DV.p_error</b><br />
        <i>Float </i> &mdash;
         Basecalling error for reference confidence model in gvcf <br />
</p>
<p name="GermlinePipeline.SegDup.DV.gq_resolution_override">
        <b>GermlinePipeline.SegDup.DV.gq_resolution_override</b><br />
        <i>Int? </i> &mdash;
         Override for gq resolution (default: 5) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.gq_bins">
        <b>GermlinePipeline.SegDup.DV.gq_bins</b><br />
        <i>Array[Int]? </i> &mdash;
         GQ bins to use instead of a fixed resolution (overrides gq_resolution) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.output_realignment">
        <b>GermlinePipeline.SegDup.DV.output_realignment</b><br />
        <i>Boolean </i> &mdash;
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ug_make_examples_extra_args">
        <b>GermlinePipeline.SegDup.DV.ug_make_examples_extra_args</b><br />
        <i>String? </i> &mdash;
         Additional arguments for make-examples tool <br />
</p>
<p name="GermlinePipeline.SegDup.DV.log_make_examples_progress">
        <b>GermlinePipeline.SegDup.DV.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash;
         Cause make_examples to output detailed progress information (for debugging) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.germline_vcf">
        <b>GermlinePipeline.SegDup.DV.germline_vcf</b><br />
        <i>File? </i> &mdash;
         Germline vcf file in order to generate haplotypes that incorporate germline variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.pangenome_haplotypes">
        <b>GermlinePipeline.SegDup.DV.pangenome_haplotypes</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram file <br />
</p>
<p name="GermlinePipeline.SegDup.DV.pangenome_haplotypes_index">
        <b>GermlinePipeline.SegDup.DV.pangenome_haplotypes_index</b><br />
        <i>File? </i> &mdash;
         Optional pangenome haplotypes cram index file <br />
</p>
<p name="GermlinePipeline.SegDup.DV.num_haplotypes">
        <b>GermlinePipeline.SegDup.DV.num_haplotypes</b><br />
        <i>Int? </i> &mdash;
         Number of haplotypes in the pangenome haplotype CRAM. Also determines the haplotype band height in the pileup image. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.include_reference_in_haplotypes">
        <b>GermlinePipeline.SegDup.DV.include_reference_in_haplotypes</b><br />
        <i>Boolean? </i> &mdash;
         Include the reference sequence in the sampled haplotypes <br />
</p>
<p name="GermlinePipeline.SegDup.DV.diploid_sampling_in_haplotypes">
        <b>GermlinePipeline.SegDup.DV.diploid_sampling_in_haplotypes</b><br />
        <i>Boolean? </i> &mdash;
         Use diploid sampling strategy for haplotype selection <br />
</p>
<p name="GermlinePipeline.SegDup.DV.optimization_level">
        <b>GermlinePipeline.SegDup.DV.optimization_level</b><br />
        <i>Int? </i> &mdash;
         Optimization level for TensorRT engine in call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.output_call_variants_tfrecords">
        <b>GermlinePipeline.SegDup.DV.output_call_variants_tfrecords</b><br />
        <i>Boolean </i> &mdash;
         Output tfrecords from call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.run_ploidy_estimation">
        <b>GermlinePipeline.SegDup.DV.run_ploidy_estimation</b><br />
        <i>Boolean </i> &mdash;
         Run VCF-based ploidy estimation and chrX/Y haploid conversion for germline samples. Default: false; enabled by germline use cases. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.sex_chromosomes">
        <b>GermlinePipeline.SegDup.DV.sex_chromosomes</b><br />
        <i>Array[String] </i> &mdash;
         Sex chromosome names to exclude from autosomal ploidy baseline. Defaults support chr-prefixed and non-prefixed human references. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.strong_call_threshold">
        <b>GermlinePipeline.SegDup.DV.strong_call_threshold</b><br />
        <i>Float </i> &mdash;
         Probability threshold for selective ensemble inference. When ensemble_size >= 2, examples with max probability below this threshold are re-evaluated using ensemble inference; examples above it are accepted as-is. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ensemble_size">
        <b>GermlinePipeline.SegDup.DV.ensemble_size</b><br />
        <i>Int </i> &mdash;
         Number of augmented passes for ensemble inference. Values <= 1 disable ensemble entirely (no augmentation is applied); values >= 2 enable selective ensemble. <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ensemble_reference_rows">
        <b>GermlinePipeline.SegDup.DV.ensemble_reference_rows</b><br />
        <i>Int </i> &mdash;
         Number of reference rows for ensemble inference <br />
</p>
<p name="GermlinePipeline.SegDup.DV.random_seed">
        <b>GermlinePipeline.SegDup.DV.random_seed</b><br />
        <i>Int </i> &mdash;
         Random seed for ensemble inference <br />
</p>
<p name="GermlinePipeline.SegDup.DV.shuffle_all_samples">
        <b>GermlinePipeline.SegDup.DV.shuffle_all_samples</b><br />
        <i>Boolean </i> &mdash;
         Whether to shuffle all samples during inference <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_variant_quality_hmer_indels">
        <b>GermlinePipeline.SegDup.DV.min_variant_quality_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal h-mer indel quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_variant_quality_non_hmer_indels">
        <b>GermlinePipeline.SegDup.DV.min_variant_quality_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal non-h-mer indel quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.SegDup.DV.min_variant_quality_snps">
        <b>GermlinePipeline.SegDup.DV.min_variant_quality_snps</b><br />
        <i>Int </i> &mdash;
         Minimal snp variant quality in order to be labeled as PASS <br />
</p>
<p name="GermlinePipeline.SegDup.DV.hard_qual_filter">
        <b>GermlinePipeline.SegDup.DV.hard_qual_filter</b><br />
        <i>Int </i> &mdash;
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br />
</p>
<p name="GermlinePipeline.SegDup.DV.allele_frequency_ratio">
        <b>GermlinePipeline.SegDup.DV.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering <br />
</p>
<p name="GermlinePipeline.SegDup.DV.h_indel_vaf_to_pass">
        <b>GermlinePipeline.SegDup.DV.h_indel_vaf_to_pass</b><br />
        <i>Float? </i> &mdash;
         Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio <br />
</p>
<p name="GermlinePipeline.SegDup.DV.h_indel_allele_frequency_ratio">
        <b>GermlinePipeline.SegDup.DV.h_indel_allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ug_post_processing_extra_args">
        <b>GermlinePipeline.SegDup.DV.ug_post_processing_extra_args</b><br />
        <i>String </i> &mdash;
         Additional arguments for post-processing <br />
</p>
<p name="GermlinePipeline.SegDup.DV.run_roh">
        <b>GermlinePipeline.SegDup.DV.run_roh</b><br />
        <i>Boolean </i> &mdash;
         Whether to call runs of homozygosity (ROH). Enabled by default in the germline WGS use-cases, off otherwise. Requires a reference genome that has a roh_blacklist resource (the hg38 builds and b37) unless roh_blacklist_override is given <br />
</p>
<p name="GermlinePipeline.SegDup.DV.roh_af_default">
        <b>GermlinePipeline.SegDup.DV.roh_af_default</b><br />
        <i>Float </i> &mdash;
         Alternate allele frequency assumed for every marker by the ROH caller, in place of a population frequency table <br />
</p>
<p name="GermlinePipeline.SegDup.DV.input_flow_order">
        <b>GermlinePipeline.SegDup.DV.input_flow_order</b><br />
        <i>String? </i> &mdash;
         Flow order. If not provided, it will be extracted from the CRAM header <br />
</p>
<p name="GermlinePipeline.SegDup.DV.call_variants_gpu_type_override">
        <b>GermlinePipeline.SegDup.DV.call_variants_gpu_type_override</b><br />
        <i>String? </i> &mdash;
         GPU type for call variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.call_variants_gpus">
        <b>GermlinePipeline.SegDup.DV.call_variants_gpus</b><br />
        <i>Int </i> &mdash;
         Number of GPUs for call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.call_variants_cpus">
        <b>GermlinePipeline.SegDup.DV.call_variants_cpus</b><br />
        <i>Int </i> &mdash;
         Number of CPUs for call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.call_variants_threads">
        <b>GermlinePipeline.SegDup.DV.call_variants_threads</b><br />
        <i>Int </i> &mdash;
         Number of decompression threads for call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.call_variants_uncompr_buf_size_gb">
        <b>GermlinePipeline.SegDup.DV.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash;
         Memory buffer allocated for each uncompression thread in calll_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.v_gpu_tile_size">
        <b>GermlinePipeline.SegDup.DV.v_gpu_tile_size</b><br />
        <i>Int </i> &mdash;
         Virtual GPU tile size for call_variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.kmer_length">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.kmer_length</b><br />
        <i>Int </i> &mdash;
         K-mer length for KMC counting (default: 29) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.min_kmer_count">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.min_kmer_count</b><br />
        <i>Int </i> &mdash;
         Minimum k-mer count threshold for sampling (default: 2) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.window_size">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.window_size</b><br />
        <i>Int </i> &mdash;
         Sliding window size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.step_size">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.step_size</b><br />
        <i>Int </i> &mdash;
         Sliding window step size for seqkit (default: 50000) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.minimap2_preset">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.minimap2_preset</b><br />
        <i>String </i> &mdash;
         Minimap2 preset for alignment (default: asm5) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.HaplotypeSampling.minimap_extra_args">
        <b>GermlinePipeline.SegDup.DV.HaplotypeSampling.minimap_extra_args</b><br />
        <i>String? </i> &mdash;
         Additional extra arguments to pass to minimap2 (default: empty) <br />
</p>

### Optional reference files
<p name="GermlinePipeline.annotation_intervals">
        <b>GermlinePipeline.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash;
         List of bed files for VCF annotation <br />
</p>
<p name="GermlinePipeline.dv_model_onnx">
        <b>GermlinePipeline.dv_model_onnx</b><br />
        <i>File? </i> &mdash;
         (VariantCalling) TensorRT model for calling variants (onnx format) <br />
</p>
<p name="GermlinePipeline.dv_model_serialized">
        <b>GermlinePipeline.dv_model_serialized</b><br />
        <i>File? </i> &mdash;
         (VariantCalling) TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided) <br />
</p>
<p name="GermlinePipeline.dv_ref_gbz_for_haplotypes">
        <b>GermlinePipeline.dv_ref_gbz_for_haplotypes</b><br />
        <i>File? </i> &mdash;
         (VariantCalling) Pangenome GBZ index file for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided) <br />
</p>
<p name="GermlinePipeline.dv_ref_hapl">
        <b>GermlinePipeline.dv_ref_hapl</b><br />
        <i>File? </i> &mdash;
         (VariantCalling) Pre-computed haplotype index file (.hapl) for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided) <br />
</p>
<p name="GermlinePipeline.ref_dbsnp">
        <b>GermlinePipeline.ref_dbsnp</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf for the annotation of known variants <br />
</p>
<p name="GermlinePipeline.ref_dbsnp_index">
        <b>GermlinePipeline.ref_dbsnp_index</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf index <br />
</p>
<p name="GermlinePipeline.roh_blacklist_override">
        <b>GermlinePipeline.roh_blacklist_override</b><br />
        <i>File? </i> &mdash;
         BED of alignment-artefact regions to exclude from the reported runs of homozygosity, overriding the genome default (ENCODE blacklist v2) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.model_serialized">
        <b>GermlinePipeline.PGx.EfficientDV.model_serialized</b><br />
        <i>File? </i> &mdash;
         TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.roh_blacklist_override">
        <b>GermlinePipeline.PGx.EfficientDV.roh_blacklist_override</b><br />
        <i>File? </i> &mdash;
         BED of alignment-artefact regions to exclude from the reported runs of homozygosity, overriding the genome default (ENCODE blacklist v2) <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.annotation_intervals">
        <b>GermlinePipeline.PGx.EfficientDV.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash;
         List of bed files for VCF annotation <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.ref_dbsnp">
        <b>GermlinePipeline.PGx.EfficientDV.ref_dbsnp</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf for the annotation of known variants <br />
</p>
<p name="GermlinePipeline.PGx.EfficientDV.ref_dbsnp_index">
        <b>GermlinePipeline.PGx.EfficientDV.ref_dbsnp_index</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf index <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ref_gbz_for_haplotypes">
        <b>GermlinePipeline.SegDup.DV.ref_gbz_for_haplotypes</b><br />
        <i>File? </i> &mdash;
         Pangenome GBZ index file for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ref_hapl">
        <b>GermlinePipeline.SegDup.DV.ref_hapl</b><br />
        <i>File? </i> &mdash;
         Pre-computed haplotype index file (.hapl) for haplotype sampling (required if run_haplotype_sampling is true and pangenome_haplotypes is not provided) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.roh_blacklist_override">
        <b>GermlinePipeline.SegDup.DV.roh_blacklist_override</b><br />
        <i>File? </i> &mdash;
         BED of alignment-artefact regions to exclude from the reported runs of homozygosity, overriding the genome default (ENCODE blacklist v2) <br />
</p>
<p name="GermlinePipeline.SegDup.DV.annotation_intervals">
        <b>GermlinePipeline.SegDup.DV.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash;
         List of bed files for VCF annotation <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ref_dbsnp">
        <b>GermlinePipeline.SegDup.DV.ref_dbsnp</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf for the annotation of known variants <br />
</p>
<p name="GermlinePipeline.SegDup.DV.ref_dbsnp_index">
        <b>GermlinePipeline.SegDup.DV.ref_dbsnp_index</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf index <br />
</p>
</details>


### Advanced inputs
<details>
<summary> Show/Hide </summary>
<p name="GermlinePipeline.annotate_variants_cpu_override">
        <b>GermlinePipeline.annotate_variants_cpu_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         cpu override for annotate_variants task
</p>
<p name="GermlinePipeline.annotate_variants_memory_override">
        <b>GermlinePipeline.annotate_variants_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for annotate_variants task
</p>
<p name="GermlinePipeline.config_file_string">
        <b>GermlinePipeline.config_file_string</b><br />
        <i>String? &mdash; Default: None</i><br />
         Gridss config file content
</p>
<p name="GermlinePipeline.convert_vcf_format_memory_override">
        <b>GermlinePipeline.convert_vcf_format_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for convert_vcf_format task
</p>
<p name="GermlinePipeline.create_assembly_memory_override">
        <b>GermlinePipeline.create_assembly_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for create_assembly task
</p>
<p name="GermlinePipeline.germline_link_variants_memory_override">
        <b>GermlinePipeline.germline_link_variants_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for germline_link_variants task
</p>
<p name="GermlinePipeline.max_reads_per_working_area">
        <b>GermlinePipeline.max_reads_per_working_area</b><br />
        <i>Int? &mdash; Default: None</i><br />
         Rematching parameter: Maximal number of reads that are stored in memory when rematching reads to haplotypes (similar to max_reads_per_partition in assembly)
</p>
<p name="GermlinePipeline.rematching_memory_override">
        <b>GermlinePipeline.rematching_memory_override</b><br />
        <i>Int? &mdash; Default: None</i><br />
         memory override for rematching task
</p>
<p name="GermlinePipeline.sv_max_reads_per_partition">
        <b>GermlinePipeline.sv_max_reads_per_partition</b><br />
        <i>Int? &mdash; Default: None</i><br />
         (SV) Assembly parameter: Maximal number of reads that are stored in memory when analyzing an active region
</p>
<p name="GermlinePipeline.sv_scatter_intervals_break">
        <b>GermlinePipeline.sv_scatter_intervals_break</b><br />
        <i>Int? &mdash; Default: None</i><br />
         (SV) Maximal resolution for scattering intervals
</p>
</details>

## Outputs
<p name="GermlinePipeline.dv_nvidia_smi_log">
        <b>GermlinePipeline.dv_nvidia_smi_log</b><br />
        <i>File?</i><br />
        Nvidia System Management (nvidia-smi) log (only when run_variant_calling)
</p>
<p name="GermlinePipeline.snv_indel_vcf">
        <b>GermlinePipeline.snv_indel_vcf</b><br />
        <i>File?</i><br />
        Called variants in vcf format (only when run_variant_calling)
</p>
<p name="GermlinePipeline.snv_indel_vcf_index">
        <b>GermlinePipeline.snv_indel_vcf_index</b><br />
        <i>File?</i><br />
        vcf index (only when run_variant_calling)
</p>
<p name="GermlinePipeline.snv_indel_vcf_no_ref_calls">
        <b>GermlinePipeline.snv_indel_vcf_no_ref_calls</b><br />
        <i>File?</i><br />
        Called variants without reference calls (only when run_variant_calling)
</p>
<p name="GermlinePipeline.snv_indel_vcf_no_ref_calls_index">
        <b>GermlinePipeline.snv_indel_vcf_no_ref_calls_index</b><br />
        <i>File?</i><br />
        vcf without references calls index (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_roh_tsv">
        <b>GermlinePipeline.dv_roh_tsv</b><br />
        <i>File?</i><br />
        Runs of homozygosity, as the regions tsv of bcftools roh (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_call_variants_output_tfrecords">
        <b>GermlinePipeline.dv_call_variants_output_tfrecords</b><br />
        <i>Array[File]?</i><br />
        The tfrecords that call_variants outputs (only when run_variant_calling)
</p>
<p name="GermlinePipeline.gvcf">
        <b>GermlinePipeline.gvcf</b><br />
        <i>File?</i><br />
        Variant in each position (gvcf file) (only when run_variant_calling)
</p>
<p name="GermlinePipeline.gvcf_index">
        <b>GermlinePipeline.gvcf_index</b><br />
        <i>File?</i><br />
        gvcf index (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_output_gvcf_hcr">
        <b>GermlinePipeline.dv_output_gvcf_hcr</b><br />
        <i>File?</i><br />
        HCR file - callability regions BED file defined from the gVCF (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_realigned_cram">
        <b>GermlinePipeline.dv_realigned_cram</b><br />
        <i>File?</i><br />
        Realigned reads cram from make_examples (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_realigned_cram_index">
        <b>GermlinePipeline.dv_realigned_cram_index</b><br />
        <i>File?</i><br />
        Realigned CRAM index (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_flow_order">
        <b>GermlinePipeline.dv_flow_order</b><br />
        <i>String?</i><br />
        Flow order (only when run_variant_calling)
</p>
<p name="GermlinePipeline.qc_report_html">
        <b>GermlinePipeline.qc_report_html</b><br />
        <i>File?</i><br />
        QC report html (only when run_variant_calling)
</p>
<p name="GermlinePipeline.qc_report_h5">
        <b>GermlinePipeline.qc_report_h5</b><br />
        <i>File?</i><br />
        QC stats in h5 file format (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_qc_metrics_h5">
        <b>GermlinePipeline.dv_qc_metrics_h5</b><br />
        <i>File?</i><br />
        QC stats in specific format for UGDV workflow (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_num_candidates">
        <b>GermlinePipeline.dv_num_candidates</b><br />
        <i>Array[File]?</i><br />
        Number of candidates that call_variants processed (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_num_candidates_as_int">
        <b>GermlinePipeline.dv_num_candidates_as_int</b><br />
        <i>Int?</i><br />
        Number of candidates that call_variants processed (as an integer) (only when run_variant_calling)
</p>
<p name="GermlinePipeline.dv_ploidy_report">
        <b>GermlinePipeline.dv_ploidy_report</b><br />
        <i>File?</i><br />
        Human-readable genome ploidy report with sex karyotype, per-chromosome ploidy, and BAF summary when available (only when run_variant_calling)
</p>
<p name="GermlinePipeline.cnv_cnmops_cnv_calls_bed">
        <b>GermlinePipeline.cnv_cnmops_cnv_calls_bed</b><br />
        <i>File?</i><br />
        CNMOPS CNV calls in bed format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_cnmops_vcf">
        <b>GermlinePipeline.cnv_cnmops_vcf</b><br />
        <i>File?</i><br />
        CNMOPS CNV calls in VCF format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_cnmops_cnv_calls_vcf_index">
        <b>GermlinePipeline.cnv_cnmops_cnv_calls_vcf_index</b><br />
        <i>File?</i><br />
        Index file for the CNMOPS CNV calls VCF (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_cnvpytor_cnv_calls_bed">
        <b>GermlinePipeline.cnv_cnvpytor_cnv_calls_bed</b><br />
        <i>File?</i><br />
        CNVpytor CNV calls in bed format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_cnvpytor_vcf">
        <b>GermlinePipeline.cnv_cnvpytor_vcf</b><br />
        <i>File?</i><br />
        CNVpytor CNV calls in VCF format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_cnvpytor_cnv_calls_vcf_index">
        <b>GermlinePipeline.cnv_cnvpytor_cnv_calls_vcf_index</b><br />
        <i>File?</i><br />
        Index file for the CNVpytor CNV calls VCF (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_bed">
        <b>GermlinePipeline.cnv_bed</b><br />
        <i>File?</i><br />
        Final (combined) CNV calls in bed format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_vcf">
        <b>GermlinePipeline.cnv_vcf</b><br />
        <i>File?</i><br />
        Combined CNV calls in vcf format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_vcf_index">
        <b>GermlinePipeline.cnv_vcf_index</b><br />
        <i>File?</i><br />
        Index of the combined CNV calls in vcf format (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_split_read_evidence">
        <b>GermlinePipeline.cnv_split_read_evidence</b><br />
        <i>File?</i><br />
        BAM file with split read evidence supporting combined CNV calls (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_split_read_evidence_index">
        <b>GermlinePipeline.cnv_split_read_evidence_index</b><br />
        <i>File?</i><br />
        Index file for the BAM with split read evidence supporting combined CNV calls (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_realign_read_evidence">
        <b>GermlinePipeline.cnv_realign_read_evidence</b><br />
        <i>File?</i><br />
        BAM file with read evidence supporting combined CNV calls (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_realign_read_evidence_index">
        <b>GermlinePipeline.cnv_realign_read_evidence_index</b><br />
        <i>File?</i><br />
        Index file for the BAM with read evidence supporting combined CNV calls (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_combine_read_scores_csv">
        <b>GermlinePipeline.cnv_combine_read_scores_csv</b><br />
        <i>File?</i><br />
        CSV file with jalign scores for each read (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_combined_coverage_plot">
        <b>GermlinePipeline.cnv_combined_coverage_plot</b><br />
        <i>File?</i><br />
        CNV coverage plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_combined_dup_del_plot">
        <b>GermlinePipeline.cnv_combined_dup_del_plot</b><br />
        <i>File?</i><br />
        Duplication and deletion calls plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_combined_copy_number_plot">
        <b>GermlinePipeline.cnv_combined_copy_number_plot</b><br />
        <i>File?</i><br />
        Copy number calls plot for combined calls in JPEG format (only generated if skip_figure_generation is false) (only when run_cnv)
</p>
<p name="GermlinePipeline.cnv_md5_checksums_json">
        <b>GermlinePipeline.cnv_md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files (only when run_cnv)
</p>
<p name="GermlinePipeline.sv_annotated_unlinked_vcf">
        <b>GermlinePipeline.sv_annotated_unlinked_vcf</b><br />
        <i>File?</i><br />
        Annotated VCF file, before GRIPSS or GermlineLinkVariants (only when run_sv)
</p>
<p name="GermlinePipeline.sv_annotated_unlinked_vcf_index">
        <b>GermlinePipeline.sv_annotated_unlinked_vcf_index</b><br />
        <i>File?</i><br />
        Annotated VCF index file (only when run_sv)
</p>
<p name="GermlinePipeline.sv_vcf">
        <b>GermlinePipeline.sv_vcf</b><br />
        <i>File?</i><br />
        Final VCF (only when run_sv)
</p>
<p name="GermlinePipeline.sv_vcf_index">
        <b>GermlinePipeline.sv_vcf_index</b><br />
        <i>File?</i><br />
        Final VCF index (only when run_sv)
</p>
<p name="GermlinePipeline.sv_assembly">
        <b>GermlinePipeline.sv_assembly</b><br />
        <i>File?</i><br />
        Raw assembly - before the realignment (only when run_sv)
</p>
<p name="GermlinePipeline.sv_assembly_index">
        <b>GermlinePipeline.sv_assembly_index</b><br />
        <i>File?</i><br />
        Raw assembly - before the realignment - index (only when run_sv)
</p>
<p name="GermlinePipeline.sv_realigned_assembly">
        <b>GermlinePipeline.sv_realigned_assembly</b><br />
        <i>File?</i><br />
        Assembly output after UA realingment (only when run_sv)
</p>
<p name="GermlinePipeline.sv_realigned_assembly_index">
        <b>GermlinePipeline.sv_realigned_assembly_index</b><br />
        <i>File?</i><br />
        Assembly output index after UA realingment (only when run_sv)
</p>
<p name="GermlinePipeline.sv_converted_vcf">
        <b>GermlinePipeline.sv_converted_vcf</b><br />
        <i>File?</i><br />
        Final VCF file in the region (non-breakend) format (only when run_sv)
</p>
<p name="GermlinePipeline.sv_converted_vcf_index">
        <b>GermlinePipeline.sv_converted_vcf_index</b><br />
        <i>File?</i><br />
        Final VCF index file in the region (non-breakend) format (only when run_sv)
</p>
<p name="GermlinePipeline.sv_md5_checksums_json">
        <b>GermlinePipeline.sv_md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files (only when run_sv)
</p>
<p name="GermlinePipeline.str_detailed_csv_files">
        <b>GermlinePipeline.str_detailed_csv_files</b><br />
        <i>Array[File]?</i><br />
        Detailed per-read alignment results in CSV format, containing alignment scores, repeat counts, and read metadata for each alignment. Empty array if output_detailed_csv=false. (only when run_str)
</p>
<p name="GermlinePipeline.str_summary_csv_files">
        <b>GermlinePipeline.str_summary_csv_files</b><br />
        <i>Array[File]?</i><br />
        Per-locus summary statistics in CSV format, aggregating alignment results across all reads for each STR locus. Empty array if output_summary_csv=false. (only when run_str)
</p>
<p name="GermlinePipeline.str_bed">
        <b>GermlinePipeline.str_bed</b><br />
        <i>File?</i><br />
        Final genotype calls in BED format for visualization in genome browsers (IGV, UCSC). Contains chromosome, start, end, and genotype information (only when run_str)
</p>
<p name="GermlinePipeline.str_vcf">
        <b>GermlinePipeline.str_vcf</b><br />
        <i>File?</i><br />
        Final genotype calls in compressed VCF format with per-allele support counts (ADSP, ADFL). Compatible with standard VCF tools. (only when run_str)
</p>
<p name="GermlinePipeline.str_vcf_index">
        <b>GermlinePipeline.str_vcf_index</b><br />
        <i>File?</i><br />
        Tabix index for the genotypes VCF file (only when run_str)
</p>
<p name="GermlinePipeline.hla_genotypes">
        <b>GermlinePipeline.hla_genotypes</b><br />
        <i>File?</i><br />
        HLA genotyping output file (only when run_hla)
</p>
<p name="GermlinePipeline.kir_genotypes">
        <b>GermlinePipeline.kir_genotypes</b><br />
        <i>File?</i><br />
        KIR genotyping output file (only when run_hla)
</p>
<p name="GermlinePipeline.pgx_allele_fraction_profiles">
        <b>GermlinePipeline.pgx_allele_fraction_profiles</b><br />
        <i>Array[File]?</i><br />
        Allele fraction profiles for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_alleles">
        <b>GermlinePipeline.pgx_alleles</b><br />
        <i>Array[File]?</i><br />
        Alleles for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_cnv_calls">
        <b>GermlinePipeline.pgx_cnv_calls</b><br />
        <i>Array[File]?</i><br />
        CNV calls for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_consolidated_variants">
        <b>GermlinePipeline.pgx_consolidated_variants</b><br />
        <i>Array[File]?</i><br />
        Consolidated variants for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_copy_number_profiles">
        <b>GermlinePipeline.pgx_copy_number_profiles</b><br />
        <i>Array[File]?</i><br />
        Copy number profiles for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_copy_numbers">
        <b>GermlinePipeline.pgx_copy_numbers</b><br />
        <i>Array[File]?</i><br />
        Copy numbers for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_genotypes">
        <b>GermlinePipeline.pgx_genotypes</b><br />
        <i>Array[File]?</i><br />
        Genotypes for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_imported_variants">
        <b>GermlinePipeline.pgx_imported_variants</b><br />
        <i>Array[File]?</i><br />
        Imported variants for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_phased_variants">
        <b>GermlinePipeline.pgx_phased_variants</b><br />
        <i>Array[File]?</i><br />
        Phased variants for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_phenotypes">
        <b>GermlinePipeline.pgx_phenotypes</b><br />
        <i>Array[File]?</i><br />
        Phenotypes for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_read_depths">
        <b>GermlinePipeline.pgx_read_depths</b><br />
        <i>Array[File]?</i><br />
        Read depths for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_vcf">
        <b>GermlinePipeline.pgx_vcf</b><br />
        <i>File?</i><br />
        Output VCF file (either the input VCF file or the one produced by Efficient DV if no input VCF file was provided) (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_vcf_index">
        <b>GermlinePipeline.pgx_vcf_index</b><br />
        <i>File?</i><br />
        Output VCF index file (either the input VCF index file or the one produced by Efficient DV if no input VCF index file was provided) (only when run_pgx)
</p>
<p name="GermlinePipeline.pgx_results">
        <b>GermlinePipeline.pgx_results</b><br />
        <i>File?</i><br />
        Results for each gene (only when run_pgx)
</p>
<p name="GermlinePipeline.segdup_remap_bam">
        <b>GermlinePipeline.segdup_remap_bam</b><br />
        <i>File?</i><br />
        Remapped BAM file (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_remap_bam_index">
        <b>GermlinePipeline.segdup_remap_bam_index</b><br />
        <i>File?</i><br />
        Remapped BAM index file (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_acnv_calls">
        <b>GermlinePipeline.segdup_acnv_calls</b><br />
        <i>File?</i><br />
        CNV calls (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_pcnv_calls">
        <b>GermlinePipeline.segdup_pcnv_calls</b><br />
        <i>File?</i><br />
        Paralog CNV calls (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_small_variants">
        <b>GermlinePipeline.segdup_small_variants</b><br />
        <i>File?</i><br />
        Small variants (VCF) combining ParascopyCall output and LPA KIV-2 targeted small variants (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_small_variants_index">
        <b>GermlinePipeline.segdup_small_variants_index</b><br />
        <i>File?</i><br />
        Small variants index (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_lpa_vcf">
        <b>GermlinePipeline.segdup_lpa_vcf</b><br />
        <i>File?</i><br />
        LPA KIV-2 targeted caller VCF (full: KIV-2 CNV symbolic record + LPA small variants) (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_lpa_vcf_index">
        <b>GermlinePipeline.segdup_lpa_vcf_index</b><br />
        <i>File?</i><br />
        LPA KIV-2 targeted caller VCF index (only when run_segdup)
</p>
<p name="GermlinePipeline.segdup_lpa_json">
        <b>GermlinePipeline.segdup_lpa_json</b><br />
        <i>File?</i><br />
        LPA KIV-2 targeted caller JSON report (only when run_segdup)
</p>

<hr />

> Generated using WDL AID (1.0.1)