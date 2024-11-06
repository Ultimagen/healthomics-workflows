# EfficientDV
Performs variant calling on an input cram, using a re-write of (DeepVariant)[https://www.nature.com/articles/nbt.4235] which is adapted for Ultima Genomics data. There are three stages to the variant calling: (1) make_examples - Looks for “active regions” with potential candidates. Within these regions, it performs local assembly (haplotypes), re-aligns the reads, and defines candidate variant. Images of the reads in the vicinity of the candidates are saved as protos in a tfrecord format. (2) call_variants - Collects the images from make_examples and uses a deep learning model to infer the statistics of each variant (i.e. quality, genotype likelihoods etc.). (3) post_process - Uses the output of call_variants to generate a vcf and annotates it.

## Inputs

### Required inputs
<p name="EfficientDV.base_file_name">
        <b>EfficientDV.base_file_name</b><br />
        <i>String </i> &mdash; 
         Prefix for name of all output files <br /> 
</p>
<p name="EfficientDV.cram_files">
        <b>EfficientDV.cram_files</b><br />
        <i>Array[File] </i> &mdash; 
         Input cram files. Multiple files are merged. <br /> 
</p>
<p name="EfficientDV.cram_index_files">
        <b>EfficientDV.cram_index_files</b><br />
        <i>Array[File] </i> &mdash; 
         Input cram index files. <br /> 
</p>

### Required parameters
<p name="EfficientDV.make_gvcf">
        <b>EfficientDV.make_gvcf</b><br />
        <i>Boolean </i> &mdash; 
         Whether to generate a gvcf. Default: False <br /> 
</p>

### Required references
<p name="EfficientDV.references">
        <b>EfficientDV.references</b><br />
        <i>References </i> &mdash; 
         Reference files: fasta, dict and fai, recommended value set in the template <br /> 
</p>
<p name="EfficientDV.model_onnx">
        <b>EfficientDV.model_onnx</b><br />
        <i>File </i> &mdash; 
         TensorRT model for calling variants (onnx format) <br /> 
</p>
<p name="EfficientDV.exome_intervals">
        <b>EfficientDV.exome_intervals</b><br />
        <i>File </i> &mdash; 
         A bed file with exome intervals <br /> 
</p>
<p name="EfficientDV.ref_dbsnp">
        <b>EfficientDV.ref_dbsnp</b><br />
        <i>File </i> &mdash; 
         DbSNP vcf for the annotation of known variants <br /> 
</p>
<p name="EfficientDV.ref_dbsnp_index">
        <b>EfficientDV.ref_dbsnp_index</b><br />
        <i>File </i> &mdash; 
         DbSNP vcf index <br /> 
</p>

### Optional inputs
<p name="EfficientDV.background_cram_files">
        <b>EfficientDV.background_cram_files</b><br />
        <i>Array[File] </i> &mdash; 
         Background (normal sample) cram files for somatic calling <br /> 
</p>
<p name="EfficientDV.background_cram_index_files">
        <b>EfficientDV.background_cram_index_files</b><br />
        <i>Array[File] </i> &mdash; 
         Background (normal sample) cram index files for somatic calling <br /> 
</p>

### Optional parameters
<p name="EfficientDV.scatter_intervals_break">
        <b>EfficientDV.scatter_intervals_break</b><br />
        <i>Int </i> &mdash; 
         The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br /> 
</p>
<p name="EfficientDV.target_intervals">
        <b>EfficientDV.target_intervals</b><br />
        <i>File? </i> &mdash; 
         Limit calling to these regions. If target_intervals and intervals_string are not provided then entire genome is used. <br /> 
</p>
<p name="EfficientDV.intervals_string">
        <b>EfficientDV.intervals_string</b><br />
        <i>String? </i> &mdash; 
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. hese regions. Takes precedence over target_intervals. If both are not provided then entire genome is used. <br /> 
</p>
<p name="EfficientDV.min_fraction_hmer_indels">
        <b>EfficientDV.min_fraction_hmer_indels</b><br />
        <i>Float </i> &mdash; 
         Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_fraction_non_hmer_indels">
        <b>EfficientDV.min_fraction_non_hmer_indels</b><br />
        <i>Float </i> &mdash; 
         Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_fraction_snps">
        <b>EfficientDV.min_fraction_snps</b><br />
        <i>Float </i> &mdash; 
         Minimal fraction of reads, that support a snp, required to  generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_read_count_snps">
        <b>EfficientDV.min_read_count_snps</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support a snp, required to  generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_read_count_hmer_indels">
        <b>EfficientDV.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_read_count_non_hmer_indels">
        <b>EfficientDV.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br /> 
</p>
<p name="EfficientDV.min_base_quality">
        <b>EfficientDV.min_base_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal base quality for candidate generation <br /> 
</p>
<p name="EfficientDV.pileup_min_mapping_quality">
        <b>EfficientDV.pileup_min_mapping_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal mapping quality to be included in image (the input to the CNN) <br /> 
</p>
<p name="EfficientDV.candidate_min_mapping_quality">
        <b>EfficientDV.candidate_min_mapping_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal mapping quality for candidate generation <br /> 
</p>
<p name="EfficientDV.max_reads_per_partition">
        <b>EfficientDV.max_reads_per_partition</b><br />
        <i>Int </i> &mdash; 
         Maximal number of reads that are stored in memory when analyzing an active region <br /> 
</p>
<p name="EfficientDV.dbg_min_base_quality">
        <b>EfficientDV.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash; 
         Minimal base quality for local assembly of haplotypes <br /> 
</p>
<p name="EfficientDV.prioritize_alt_supporting_reads">
        <b>EfficientDV.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash; 
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br /> 
</p>
<p name="EfficientDV.p_error">
        <b>EfficientDV.p_error</b><br />
        <i>Float </i> &mdash; 
         Basecalling error for reference confidence model in gvcf <br /> 
</p>
<p name="EfficientDV.output_realignment">
        <b>EfficientDV.output_realignment</b><br />
        <i>Boolean </i> &mdash; 
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br /> 
</p>
<p name="EfficientDV.ug_make_examples_extra_args">
        <b>EfficientDV.ug_make_examples_extra_args</b><br />
        <i>String </i> &mdash; 
         Additional arguments for make-examples tool <br /> 
</p>
<p name="EfficientDV.log_make_examples_progress">
        <b>EfficientDV.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash; 
         Cause make_examples to output detailed progress information (for debugging) <br /> 
</p>
<p name="EfficientDV.min_variant_quality_hmer_indels">
        <b>EfficientDV.min_variant_quality_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal h-mer indel quality in order to be labeled as PASS <br /> 
</p>
<p name="EfficientDV.min_variant_quality_non_hmer_indels">
        <b>EfficientDV.min_variant_quality_non_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal non-h-mer indel quality in order to be labeled as PASS <br /> 
</p>
<p name="EfficientDV.min_variant_quality_snps">
        <b>EfficientDV.min_variant_quality_snps</b><br />
        <i>Int </i> &mdash; 
         Minimal snp variant quality in order to be labeled as PASS <br /> 
</p>
<p name="EfficientDV.min_variant_quality_exome_hmer_indels">
        <b>EfficientDV.min_variant_quality_exome_hmer_indels</b><br />
        <i>Int </i> &mdash; 
         Minimal non-h-mer indel quality in order to be labeled as PASS <br /> 
</p>
<p name="EfficientDV.hard_qual_filter">
        <b>EfficientDV.hard_qual_filter</b><br />
        <i>Int </i> &mdash; 
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br /> 
</p>
<p name="EfficientDV.allele_frequency_ratio">
        <b>EfficientDV.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash; 
         Minimal ratio between the allele frequency in tumor and normal, for vcf filtering <br /> 
</p>
<p name="EfficientDV.show_bg_fields">
        <b>EfficientDV.show_bg_fields</b><br />
        <i>Boolean </i> &mdash; 
         Show background statistics BG_AD, BG_SB in the output VCF (relevant for somatic calling) <br /> 
</p>
<p name="EfficientDV.annotate_systematic_errors">
        <b>EfficientDV.annotate_systematic_errors</b><br />
        <i>Boolean </i> &mdash; 
         Should systematic errors be annotated from a database of common systematic errors <br /> 
</p>
<p name="EfficientDV.input_flow_order">
        <b>EfficientDV.input_flow_order</b><br />
        <i>String? </i> &mdash; 
         Flow order. If not provided, it will be extracted from the CRAM header <br /> 
</p>
<p name="EfficientDV.call_variants_gpu_type">
        <b>EfficientDV.call_variants_gpu_type</b><br />
        <i>String </i> &mdash; 
         GPU type for call variants <br /> 
</p>
<p name="EfficientDV.call_variants_gpus">
        <b>EfficientDV.call_variants_gpus</b><br />
        <i>Int </i> &mdash; 
         Number of GPUs for call_variants <br /> 
</p>
<p name="EfficientDV.call_variants_cpus">
        <b>EfficientDV.call_variants_cpus</b><br />
        <i>Int </i> &mdash; 
         Number of CPUs for call_variants <br /> 
</p>
<p name="EfficientDV.call_variants_threads">
        <b>EfficientDV.call_variants_threads</b><br />
        <i>Int </i> &mdash; 
         Number of decompression threads for call_variants <br /> 
</p>
<p name="EfficientDV.call_variants_uncompr_buf_size_gb">
        <b>EfficientDV.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash; 
         Memory buffer allocated for each uncompression thread in calll_variants <br /> 
</p>

### Optional reference files
<p name="EfficientDV.model_serialized">
        <b>EfficientDV.model_serialized</b><br />
        <i>File? </i> &mdash; 
         TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided) <br /> 
</p>
<p name="EfficientDV.annotation_intervals">
        <b>EfficientDV.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash; 
         List of bed files for VCF annotation <br /> 
</p>
</details>


## Outputs
<p name="EfficientDV.nvidia_smi_log">
        <b>EfficientDV.nvidia_smi_log</b><br />
        <i>File</i><br />
        Nvidia System Management (nvidia-smi) log
</p>
<p name="EfficientDV.vcf_file">
        <b>EfficientDV.vcf_file</b><br />
        <i>File</i><br />
        Called variants in vcf format
</p>
<p name="EfficientDV.vcf_index">
        <b>EfficientDV.vcf_index</b><br />
        <i>File</i><br />
        vcf index
</p>
<p name="EfficientDV.vcf_no_ref_calls">
        <b>EfficientDV.vcf_no_ref_calls</b><br />
        <i>File</i><br />
        Called variants without reference calls
</p>
<p name="EfficientDV.vcf_no_ref_calls_index">
        <b>EfficientDV.vcf_no_ref_calls_index</b><br />
        <i>File</i><br />
        vcf without references calls index
</p>
<p name="EfficientDV.output_gvcf">
        <b>EfficientDV.output_gvcf</b><br />
        <i>File?</i><br />
        Variant in each position (gvcf file)
</p>
<p name="EfficientDV.output_gvcf_index">
        <b>EfficientDV.output_gvcf_index</b><br />
        <i>File?</i><br />
        gvcf index
</p>
<p name="EfficientDV.output_gvcf_hcr">
        <b>EfficientDV.output_gvcf_hcr</b><br />
        <i>File?</i><br />
        HCR file - callability regions BED file defined from the gVCF
</p>
<p name="EfficientDV.realigned_cram">
        <b>EfficientDV.realigned_cram</b><br />
        <i>File?</i><br />
        Realigned reads cram from make_examples
</p>
<p name="EfficientDV.realigned_cram_index">
        <b>EfficientDV.realigned_cram_index</b><br />
        <i>File?</i><br />
        Realigned CRAM index
</p>
<p name="EfficientDV.flow_order">
        <b>EfficientDV.flow_order</b><br />
        <i>String</i><br />
        Flow order
</p>
<p name="EfficientDV.report_html">
        <b>EfficientDV.report_html</b><br />
        <i>File</i><br />
        QC report html
</p>
<p name="EfficientDV.qc_h5">
        <b>EfficientDV.qc_h5</b><br />
        <i>File</i><br />
        QC stats in h5 file format
</p>
<p name="EfficientDV.qc_metrics_h5">
        <b>EfficientDV.qc_metrics_h5</b><br />
        <i>File</i><br />
        QC stats in specific format for UGDV workflow
</p>

<hr />

> Generated using WDL AID (1.0.1)
