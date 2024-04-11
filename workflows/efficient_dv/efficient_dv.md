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
         Input cram file. At the moment only a single file is supported. <br /> 
</p>
<p name="EfficientDV.cram_index_files">
        <b>EfficientDV.cram_index_files</b><br />
        <i>Array[File] </i> &mdash; 
         Input cram index file. At the moment only a single file is supported. <br /> 
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
         Machine-learning model for calling variants (onnx format) <br /> 
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
<p name="EfficientDV.target_intervals">
        <b>EfficientDV.target_intervals</b><br />
        <i>File? </i> &mdash; 
         Limit calling to these regions. If not provided then entire genome is used. <br /> 
</p>
<p name="EfficientDV.output_realignment">
        <b>EfficientDV.output_realignment</b><br />
        <i>Boolean </i> &mdash; 
         Output haplotypes and re-aligned reads to a bam file. Default: false. <br /> 
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
<p name="EfficientDV.hmer_runs_bed">
        <b>EfficientDV.hmer_runs_bed</b><br />
        <i>File? </i> &mdash; 
         Bed file annotating all homopolymer runs longer than 7 in the reference genome. Used to annotate potentially difficult to sequence regions <br /> 
</p>
<p name="EfficientDV.input_flow_order">
        <b>EfficientDV.input_flow_order</b><br />
        <i>String? </i> &mdash; 
         Flow order. If not provided, it will be extracted from the CRAM header <br /> 
</p>

### Optional reference files
<p name="EfficientDV.model_serialized">
        <b>EfficientDV.model_serialized</b><br />
        <i>File? </i> &mdash; 
         Machine-learning model for calling variants, serialized for a specific platform (it is regenerated if not provided) <br /> 
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

> Generated using WDL AID (1.0.0)
