# SegDupAnalysis
Processes segmental duplications in the genome by collapsing all copies on a single copy of the segmental duplication

## Inputs

### Required inputs
<p name="SegDupAnalysis.base_file_name">
        <b>SegDupAnalysis.base_file_name</b><br />
        <i>String </i> &mdash; 
         Prefix of the output files <br /> 
</p>
<p name="SegDupAnalysis.input_cram_bam">
        <b>SegDupAnalysis.input_cram_bam</b><br />
        <i>File </i> &mdash; 
         Input CRAM/BAM file <br /> 
</p>
<p name="SegDupAnalysis.input_crai_bai">
        <b>SegDupAnalysis.input_crai_bai</b><br />
        <i>File </i> &mdash; 
         Input CRAM/BAM index file <br /> 
</p>
<p name="SegDupAnalysis.references">
        <b>SegDupAnalysis.references</b><br />
        <i>References </i> &mdash; 
         Reference genome files <br /> 
</p>
<p name="SegDupAnalysis.n_threads">
        <b>SegDupAnalysis.n_threads</b><br />
        <i>Int </i> &mdash; 
         Number of threads to use <br /> 
</p>
<p name="SegDupAnalysis.exome_intervals">
        <b>SegDupAnalysis.exome_intervals</b><br />
        <i>File </i> &mdash; 
         Exome intervals for variant calling (required for deepVariant, otherwise not important) <br /> 
</p>
<p name="SegDupAnalysis.dbsnp">
        <b>SegDupAnalysis.dbsnp</b><br />
        <i>File </i> &mdash; 
         dbSNP reference file (for annotation) <br /> 
</p>
<p name="SegDupAnalysis.dbsnp_index">
        <b>SegDupAnalysis.dbsnp_index</b><br />
        <i>File </i> &mdash; 
         dbSNP reference index file (for annotation) <br /> 
</p>

### Optional inputs
<p name="SegDupAnalysis.cloud_provider_override">
        <b>SegDupAnalysis.cloud_provider_override</b><br />
        <i>String? </i> &mdash; 
         Cloud provider override (for running on other clouds): gcp or aws <br /> 
</p>
<p name="SegDupAnalysis.preemptible_tries">
        <b>SegDupAnalysis.preemptible_tries</b><br />
        <i>Int </i> &mdash; 
         Number of preemptible tries <br /> 
</p>
<p name="SegDupAnalysis.DV.background_cram_files">
        <b>SegDupAnalysis.DV.background_cram_files</b><br />
        <i>Array[File] </i> &mdash; 
         Background (normal sample) cram files for somatic calling <br /> 
</p>
<p name="SegDupAnalysis.DV.background_cram_index_files">
        <b>SegDupAnalysis.DV.background_cram_index_files</b><br />
        <i>Array[File] </i> &mdash; 
         Background (normal sample) cram index files for somatic calling <br /> 
</p>

### Optional parameters
<p name="SegDupAnalysis.DV.show_bg_fields">
        <b>SegDupAnalysis.DV.show_bg_fields</b><br />
        <i>Boolean </i> &mdash; 
         Show background fields in the output vcf. Default: false. Mostly relevant for somatic calling. <br /> 
</p>
<p name="SegDupAnalysis.DV.scatter_intervals_break">
        <b>SegDupAnalysis.DV.scatter_intervals_break</b><br />
        <i>Int </i> &mdash; 
         The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br /> 
</p>
<p name="SegDupAnalysis.DV.intervals_string">
        <b>SegDupAnalysis.DV.intervals_string</b><br />
        <i>String? </i> &mdash; 
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. hese regions. Takes precedence over target_intervals. If both are not provided then entire genome is used. <br /> 
</p>
<p name="SegDupAnalysis.DV.min_read_count_snps">
        <b>SegDupAnalysis.DV.min_read_count_snps</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support a snp, required to  generate a candidate variant <br /> 
</p>
<p name="SegDupAnalysis.DV.min_read_count_hmer_indels">
        <b>SegDupAnalysis.DV.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="SegDupAnalysis.DV.min_read_count_non_hmer_indels">
        <b>SegDupAnalysis.DV.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="SegDupAnalysis.DV.min_base_quality">
        <b>SegDupAnalysis.DV.min_base_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal base quality for candidate generation <br /> 
</p>
<p name="SegDupAnalysis.DV.pileup_min_mapping_quality">
        <b>SegDupAnalysis.DV.pileup_min_mapping_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal mapping quality to be included in image (the input to the CNN) <br /> 
</p>
<p name="SegDupAnalysis.DV.candidate_min_mapping_quality">
        <b>SegDupAnalysis.DV.candidate_min_mapping_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal mapping quality for candidate generation <br /> 
</p>
<p name="SegDupAnalysis.DV.max_reads_per_partition">
        <b>SegDupAnalysis.DV.max_reads_per_partition</b><br />
        <i>Int </i> &mdash; 
         Maximal number of reads that are stored in memory when analyzing an active region <br /> 
</p>
<p name="SegDupAnalysis.DV.dbg_min_base_quality">
        <b>SegDupAnalysis.DV.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal base quality for local assembly of haplotypes <br /> 
</p>
<p name="SegDupAnalysis.DV.prioritize_alt_supporting_reads">
        <b>SegDupAnalysis.DV.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash; 
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br /> 
</p>
<p name="SegDupAnalysis.DV.p_error">
        <b>SegDupAnalysis.DV.p_error</b><br />
        <i>Float </i> &mdash; 
         Basecalling error for reference confidence model in gvcf <br /> 
</p>
<p name="SegDupAnalysis.DV.output_realignment">
        <b>SegDupAnalysis.DV.output_realignment</b><br />
        <i>Boolean </i> &mdash; 
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br /> 
</p>
<p name="SegDupAnalysis.DV.ug_make_examples_extra_args">
        <b>SegDupAnalysis.DV.ug_make_examples_extra_args</b><br />
        <i>String? </i> &mdash; 
         Additional arguments for make-examples tool <br /> 
</p>
<p name="SegDupAnalysis.DV.log_make_examples_progress">
        <b>SegDupAnalysis.DV.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash; 
         Cause make_examples to output detailed progress information (for debugging) <br /> 
</p>
<p name="SegDupAnalysis.DV.germline_vcf">
        <b>SegDupAnalysis.DV.germline_vcf</b><br />
        <i>File? </i> &mdash; 
         Germline vcf file in order to generate haplotypes that incorporate germline variants <br /> 
</p>
<p name="SegDupAnalysis.DV.optimization_level">
        <b>SegDupAnalysis.DV.optimization_level</b><br />
        <i>Int? </i> &mdash; 
         Optimization level for TensorRT engine in call_variants <br /> 
</p>
<p name="SegDupAnalysis.DV.output_call_variants_tfrecords">
        <b>SegDupAnalysis.DV.output_call_variants_tfrecords</b><br />
        <i>Boolean </i> &mdash; 
         Output tfrecords from call_variants <br /> 
</p>
<p name="SegDupAnalysis.DV.min_variant_quality_hmer_indels">
        <b>SegDupAnalysis.DV.min_variant_quality_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal h-mer indel quality in order to be labeled as PASS <br /> 
</p>
<p name="SegDupAnalysis.DV.min_variant_quality_non_hmer_indels">
        <b>SegDupAnalysis.DV.min_variant_quality_non_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal non-h-mer indel quality in order to be labeled as PASS <br /> 
</p>
<p name="SegDupAnalysis.DV.min_variant_quality_snps">
        <b>SegDupAnalysis.DV.min_variant_quality_snps</b><br />
        <i>Int </i> &mdash; 
         Minimal snp variant quality in order to be labeled as PASS <br /> 
</p>
<p name="SegDupAnalysis.DV.hard_qual_filter">
        <b>SegDupAnalysis.DV.hard_qual_filter</b><br />
        <i>Int </i> &mdash; 
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br /> 
</p>
<p name="SegDupAnalysis.DV.allele_frequency_ratio">
        <b>SegDupAnalysis.DV.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash; 
         Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering <br /> 
</p>
<p name="SegDupAnalysis.DV.h_indel_vaf_to_pass">
        <b>SegDupAnalysis.DV.h_indel_vaf_to_pass</b><br />
        <i>Float? </i> &mdash; 
         Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio <br /> 
</p>
<p name="SegDupAnalysis.DV.h_indel_allele_frequency_ratio">
        <b>SegDupAnalysis.DV.h_indel_allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash; 
         Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering <br /> 
</p>
<p name="SegDupAnalysis.DV.ug_post_processing_extra_args">
        <b>SegDupAnalysis.DV.ug_post_processing_extra_args</b><br />
        <i>String </i> &mdash; 
         Additional arguments for post-processing <br /> 
</p>
<p name="SegDupAnalysis.DV.input_flow_order">
        <b>SegDupAnalysis.DV.input_flow_order</b><br />
        <i>String? </i> &mdash; 
         Flow order. If not provided, it will be extracted from the CRAM header <br /> 
</p>
<p name="SegDupAnalysis.DV.call_variants_gpu_type">
        <b>SegDupAnalysis.DV.call_variants_gpu_type</b><br />
        <i>String </i> &mdash; 
         GPU type for call variants <br /> 
</p>
<p name="SegDupAnalysis.DV.call_variants_gpus">
        <b>SegDupAnalysis.DV.call_variants_gpus</b><br />
        <i>Int </i> &mdash; 
         Number of GPUs for call_variants <br /> 
</p>
<p name="SegDupAnalysis.DV.call_variants_cpus">
        <b>SegDupAnalysis.DV.call_variants_cpus</b><br />
        <i>Int </i> &mdash; 
         Number of CPUs for call_variants <br /> 
</p>
<p name="SegDupAnalysis.DV.call_variants_threads">
        <b>SegDupAnalysis.DV.call_variants_threads</b><br />
        <i>Int </i> &mdash; 
         Number of decompression threads for call_variants <br /> 
</p>
<p name="SegDupAnalysis.DV.call_variants_uncompr_buf_size_gb">
        <b>SegDupAnalysis.DV.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash; 
         Memory buffer allocated for each uncompression thread in calll_variants <br /> 
</p>

### Optional reference files
<p name="SegDupAnalysis.DV.annotation_intervals">
        <b>SegDupAnalysis.DV.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash; 
         List of bed files for VCF annotation <br /> 
</p>
</details>


## Outputs
<p name="SegDupAnalysis.remap_bam">
        <b>SegDupAnalysis.remap_bam</b><br />
        <i>File</i><br />
        Remapped BAM file
</p>
<p name="SegDupAnalysis.remap_bam_index">
        <b>SegDupAnalysis.remap_bam_index</b><br />
        <i>File</i><br />
        Remapped BAM index file
</p>
<p name="SegDupAnalysis.acnv_calls">
        <b>SegDupAnalysis.acnv_calls</b><br />
        <i>File</i><br />
        CNV calls
</p>
<p name="SegDupAnalysis.acnv_calls_index">
        <b>SegDupAnalysis.acnv_calls_index</b><br />
        <i>File</i><br />
        CNV calls index
</p>
<p name="SegDupAnalysis.small_variants">
        <b>SegDupAnalysis.small_variants</b><br />
        <i>File</i><br />
        Small variants (VCF)
</p>
<p name="SegDupAnalysis.small_variants_idx">
        <b>SegDupAnalysis.small_variants_idx</b><br />
        <i>File</i><br />
        Small variants index
</p>

<hr />

> Generated using WDL AID (1.0.1)
