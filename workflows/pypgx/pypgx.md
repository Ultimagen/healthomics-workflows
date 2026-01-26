# PyPGx
Runs pharmacogenomics analysis on several genes.

## Inputs

### Required inputs
<p name="PyPGx.base_file_name">
        <b>PyPGx.base_file_name</b><br />
        <i>String </i> &mdash;
         Prefix for name of all output files <br />
</p>
<p name="PyPGx.cram_file">
        <b>PyPGx.cram_file</b><br />
        <i>File </i> &mdash;
         Input (bam/cram) file <br />
</p>
<p name="PyPGx.cram_index_file">
        <b>PyPGx.cram_index_file</b><br />
        <i>File </i> &mdash;
         Input (bai/crai) index file <br />
</p>
<p name="PyPGx.gene_symbols">
        <b>PyPGx.gene_symbols</b><br />
        <i>Array[String] </i> &mdash;
         List of gene symbols to analyze <br />
</p>
<p name="PyPGx.ref_files_for_tarball">
        <b>PyPGx.ref_files_for_tarball</b><br />
        <i>Array[File] </i> &mdash;
         List of references for CreateReferenceCache task. <br />
</p>

### Required references
<p name="PyPGx.references">
        <b>PyPGx.references</b><br />
        <i>References </i> &mdash;
         Reference files: fasta, dict and fai <br />
</p>

### Optional inputs
<p name="PyPGx.input_vcf_file">
        <b>PyPGx.input_vcf_file</b><br />
        <i>File? </i> &mdash;
         Input VCF file with variants. If not provided, Efficient DV will be run <br />
</p>
<p name="PyPGx.input_vcf_index_file">
        <b>PyPGx.input_vcf_index_file</b><br />
        <i>File? </i> &mdash;
         Input VCF index file <br />
</p>
<p name="PyPGx.model_onnx">
        <b>PyPGx.model_onnx</b><br />
        <i>File? </i> &mdash;
         TensorRT model for calling variants (onnx format) <br />
</p>
<p name="PyPGx.exome_intervals">
        <b>PyPGx.exome_intervals</b><br />
        <i>File? </i> &mdash;
         A bed file with exome intervals. Used at the post-processing step to annotate the vcf and modify the FILTER of variants in the exome. <br />
</p>
<p name="PyPGx.ref_dbsnp">
        <b>PyPGx.ref_dbsnp</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf for the annotation of known variants <br />
</p>
<p name="PyPGx.ref_dbsnp_index">
        <b>PyPGx.ref_dbsnp_index</b><br />
        <i>File? </i> &mdash;
         DbSNP vcf index <br />
</p>

### Optional parameters
<p name="PyPGx.EfficientDV.show_bg_fields">
        <b>PyPGx.EfficientDV.show_bg_fields</b><br />
        <i>Boolean </i> &mdash;
         Show background fields in the output vcf. Default: false. Mostly relevant for somatic calling. <br />
</p>
<p name="PyPGx.EfficientDV.scatter_intervals_break">
        <b>PyPGx.EfficientDV.scatter_intervals_break</b><br />
        <i>Int </i> &mdash;
         The length of the intervals for parallelization are multiples of scatter_intervals_break. This is also the maximal length of the intervals. <br />
</p>
<p name="PyPGx.EfficientDV.intervals_string">
        <b>PyPGx.EfficientDV.intervals_string</b><br />
        <i>String? </i> &mdash;
         Regions for variant calling, in the format chrom:start-end. Multiple regions are separated by semi-colon. hese regions. Takes precedence over target_intervals. If both are not provided then entire genome is used. <br />
</p>
<p name="PyPGx.EfficientDV.min_fraction_hmer_indels">
        <b>PyPGx.EfficientDV.min_fraction_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_fraction_non_hmer_indels">
        <b>PyPGx.EfficientDV.min_fraction_non_hmer_indels</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_fraction_snps">
        <b>PyPGx.EfficientDV.min_fraction_snps</b><br />
        <i>Float </i> &mdash;
         Minimal fraction of reads, that support a snp, required to  generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_read_count_snps">
        <b>PyPGx.EfficientDV.min_read_count_snps</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a snp, required to  generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_read_count_hmer_indels">
        <b>PyPGx.EfficientDV.min_read_count_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support an h-mer indel, required to generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_read_count_non_hmer_indels">
        <b>PyPGx.EfficientDV.min_read_count_non_hmer_indels</b><br />
        <i>Int </i> &mdash;
         Minimal number of reads, that support a non-h-mer indel, required to generate a candidate variant <br />
</p>
<p name="PyPGx.EfficientDV.min_base_quality">
        <b>PyPGx.EfficientDV.min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for candidate generation <br />
</p>
<p name="PyPGx.EfficientDV.min_hmer_plus_one_candidate">
        <b>PyPGx.EfficientDV.min_hmer_plus_one_candidate</b><br />
        <i>Int </i> &mdash;
         Minimal hmer length, above which more 1-bp insertion candidates are generated, provided they also meet allele frequency conditions <br />
</p>
<p name="PyPGx.EfficientDV.max_reads_per_partition">
        <b>PyPGx.EfficientDV.max_reads_per_partition</b><br />
        <i>Int </i> &mdash;
         Maximal number of reads that are stored in memory when analyzing an active region <br />
</p>
<p name="PyPGx.EfficientDV.dbg_min_base_quality">
        <b>PyPGx.EfficientDV.dbg_min_base_quality</b><br />
        <i>Int </i> &mdash;
         Minimal base quality for local assembly of haplotypes <br />
</p>
<p name="PyPGx.EfficientDV.prioritize_alt_supporting_reads">
        <b>PyPGx.EfficientDV.prioritize_alt_supporting_reads</b><br />
        <i>Boolean </i> &mdash;
         Generate an image with all available alt-supporting reads, and only then add non-supporting reads <br />
</p>
<p name="PyPGx.EfficientDV.p_error">
        <b>PyPGx.EfficientDV.p_error</b><br />
        <i>Float </i> &mdash;
         Basecalling error for reference confidence model in gvcf <br />
</p>
<p name="PyPGx.EfficientDV.output_realignment">
        <b>PyPGx.EfficientDV.output_realignment</b><br />
        <i>Boolean </i> &mdash;
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br />
</p>
<p name="PyPGx.EfficientDV.ug_make_examples_extra_args">
        <b>PyPGx.EfficientDV.ug_make_examples_extra_args</b><br />
        <i>String? </i> &mdash;
         Additional arguments for make-examples tool <br />
</p>
<p name="PyPGx.EfficientDV.log_make_examples_progress">
        <b>PyPGx.EfficientDV.log_make_examples_progress</b><br />
        <i>Boolean </i> &mdash;
         Cause make_examples to output detailed progress information (for debugging) <br />
</p>
<p name="PyPGx.EfficientDV.germline_vcf">
        <b>PyPGx.EfficientDV.germline_vcf</b><br />
        <i>File? </i> &mdash;
         Germline vcf file in order to generate haplotypes that incorporate germline variants <br />
</p>
<p name="PyPGx.EfficientDV.optimization_level">
        <b>PyPGx.EfficientDV.optimization_level</b><br />
        <i>Int? </i> &mdash;
         Optimization level for TensorRT engine in call_variants <br />
</p>
<p name="PyPGx.EfficientDV.output_call_variants_tfrecords">
        <b>PyPGx.EfficientDV.output_call_variants_tfrecords</b><br />
        <i>Boolean </i> &mdash;
         Output tfrecords from call_variants <br />
</p>
<p name="PyPGx.EfficientDV.hard_qual_filter">
        <b>PyPGx.EfficientDV.hard_qual_filter</b><br />
        <i>Int </i> &mdash;
         Any variant with QUAL < hard_qual_filter will be discarded from the VCF file <br />
</p>
<p name="PyPGx.EfficientDV.allele_frequency_ratio">
        <b>PyPGx.EfficientDV.allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for non h indels and snvs, for vcf filtering <br />
</p>
<p name="PyPGx.EfficientDV.h_indel_vaf_to_pass">
        <b>PyPGx.EfficientDV.h_indel_vaf_to_pass</b><br />
        <i>Float? </i> &mdash;
         Minimal variant allele frequency for h-indels to not filter out by allele frequency ratio <br />
</p>
<p name="PyPGx.EfficientDV.h_indel_allele_frequency_ratio">
        <b>PyPGx.EfficientDV.h_indel_allele_frequency_ratio</b><br />
        <i>Float? </i> &mdash;
         Minimal ratio between the allele frequency in tumor and normal for h-indels for vcf filtering <br />
</p>
<p name="PyPGx.EfficientDV.ug_post_processing_extra_args">
        <b>PyPGx.EfficientDV.ug_post_processing_extra_args</b><br />
        <i>String </i> &mdash;
         Additional arguments for post-processing <br />
</p>
<p name="PyPGx.EfficientDV.input_flow_order">
        <b>PyPGx.EfficientDV.input_flow_order</b><br />
        <i>String? </i> &mdash;
         Flow order. If not provided, it will be extracted from the CRAM header <br />
</p>
<p name="PyPGx.EfficientDV.call_variants_gpu_type">
        <b>PyPGx.EfficientDV.call_variants_gpu_type</b><br />
        <i>String </i> &mdash;
         GPU type for call variants <br />
</p>
<p name="PyPGx.EfficientDV.call_variants_gpus">
        <b>PyPGx.EfficientDV.call_variants_gpus</b><br />
        <i>Int </i> &mdash;
         Number of GPUs for call_variants <br />
</p>
<p name="PyPGx.EfficientDV.call_variants_cpus">
        <b>PyPGx.EfficientDV.call_variants_cpus</b><br />
        <i>Int </i> &mdash;
         Number of CPUs for call_variants <br />
</p>
<p name="PyPGx.EfficientDV.call_variants_threads">
        <b>PyPGx.EfficientDV.call_variants_threads</b><br />
        <i>Int </i> &mdash;
         Number of decompression threads for call_variants <br />
</p>
<p name="PyPGx.EfficientDV.call_variants_uncompr_buf_size_gb">
        <b>PyPGx.EfficientDV.call_variants_uncompr_buf_size_gb</b><br />
        <i>Int </i> &mdash;
         Memory buffer allocated for each uncompression thread in calll_variants <br />
</p>

### Optional reference files
<p name="PyPGx.EfficientDV.model_serialized">
        <b>PyPGx.EfficientDV.model_serialized</b><br />
        <i>File? </i> &mdash;
         TensorRT model for calling variants, serialized for a specific platform (it is regenerated if not provided) <br />
</p>
<p name="PyPGx.EfficientDV.annotation_intervals">
        <b>PyPGx.EfficientDV.annotation_intervals</b><br />
        <i>Array[File]? </i> &mdash;
         List of bed files for VCF annotation <br />
</p>
</details>


## Outputs
<p name="PyPGx.allele_fraction_profiles">
        <b>PyPGx.allele_fraction_profiles</b><br />
        <i>Array[File]</i><br />
        Allele fraction profiles for each gene
</p>
<p name="PyPGx.alleles">
        <b>PyPGx.alleles</b><br />
        <i>Array[File]</i><br />
        Alleles for each gene
</p>
<p name="PyPGx.cnv_calls">
        <b>PyPGx.cnv_calls</b><br />
        <i>Array[File]</i><br />
        CNV calls for each gene
</p>
<p name="PyPGx.consolidated_variants">
        <b>PyPGx.consolidated_variants</b><br />
        <i>Array[File]</i><br />
        Consolidated variants for each gene
</p>
<p name="PyPGx.copy_number_profiles">
        <b>PyPGx.copy_number_profiles</b><br />
        <i>Array[File]</i><br />
        Copy number profiles for each gene
</p>
<p name="PyPGx.copy_numbers">
        <b>PyPGx.copy_numbers</b><br />
        <i>Array[File]</i><br />
        Copy numbers for each gene
</p>
<p name="PyPGx.genotypes">
        <b>PyPGx.genotypes</b><br />
        <i>Array[File]</i><br />
        Genotypes for each gene
</p>
<p name="PyPGx.imported_variants">
        <b>PyPGx.imported_variants</b><br />
        <i>Array[File]</i><br />
        Imported variants for each gene
</p>
<p name="PyPGx.phased_variants">
        <b>PyPGx.phased_variants</b><br />
        <i>Array[File]</i><br />
        Phased variants for each gene
</p>
<p name="PyPGx.phenotypes">
        <b>PyPGx.phenotypes</b><br />
        <i>Array[File]</i><br />
        Phenotypes for each gene
</p>
<p name="PyPGx.read_depths">
        <b>PyPGx.read_depths</b><br />
        <i>Array[File]</i><br />
        Read depths for each gene
</p>
<p name="PyPGx.results">
        <b>PyPGx.results</b><br />
        <i>File</i><br />
        Results for each gene
</p>

<hr />

> Generated using WDL AID (1.0.1)