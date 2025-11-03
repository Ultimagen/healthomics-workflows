# SingleSampleCnmopsCNVCalling
Runs single sample germline CNV calling workflow based on [cn.mops](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)

The pipeline uses a given cohort's coverage profile for normalization.

The pipeline can recieve one of the following options as input:

&nbsp;&nbsp;1. Input CRAM/BAM file. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_template.json

&nbsp;&nbsp;2. A rds file which stores a GenomicRanges object with coverage collected in the same windows as the given cohort. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_skip_reads_count_template.json

&nbsp;&nbsp;3. A BedGraph holding the coverage per location. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_input_bedGraph_template.json

The pipeline calls CNVs for the given sample and filters them by length (>10,000b) and overlap with UG-CNV-LCR.

<b>When Running in AWS HealthOmics this pipeline should run with [dynamic storage](https://docs.omics.ai/products/workbench/engines/parameters/aws-healthomics#storage_type-dynamic-or-static)</b>

## Inputs

### Required inputs
<p name="SingleSampleCnmopsCNVCalling.base_file_name">
        <b>SingleSampleCnmopsCNVCalling.base_file_name</b><br />
        <i>String </i> &mdash; 
         Sample name <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.cohort_reads_count_matrix">
        <b>SingleSampleCnmopsCNVCalling.cohort_reads_count_matrix</b><br />
        <i>File </i> &mdash; 
         GenomicRanges object of the cohort reads count matrix in rds file format. default cohort can be found in the template. can be created by cn.mops::getReadCountsFromBAM R function  <br /> 
</p>

### Required parameters
<p name="SingleSampleCnmopsCNVCalling.mapq">
        <b>SingleSampleCnmopsCNVCalling.mapq</b><br />
        <i>Int </i> &mdash; 
         Reads mapping-quality cutoff for coverage aggregation, recommended value set in the template <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.ref_seq_names">
        <b>SingleSampleCnmopsCNVCalling.ref_seq_names</b><br />
        <i>Array[String] </i> &mdash; 
         Chromosome names for which coverage will be calculated <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.window_length">
        <b>SingleSampleCnmopsCNVCalling.window_length</b><br />
        <i>Int </i> &mdash; 
         Window length on which the read counts will be aggregated <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.min_cnv_length">
        <b>SingleSampleCnmopsCNVCalling.min_cnv_length</b><br />
        <i>Int </i> &mdash; 
         Minimum length for reporting CNV. Default is: 10,000 <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.intersection_cutoff">
        <b>SingleSampleCnmopsCNVCalling.intersection_cutoff</b><br />
        <i>Float </i> &mdash; 
         Intersection cutoff with UG-CNV-LCR regions to filter out CNV calls. Default is:  0.5 <br /> 
</p>

### Required references
<p name="SingleSampleCnmopsCNVCalling.reference_genome">
        <b>SingleSampleCnmopsCNVCalling.reference_genome</b><br />
        <i>File </i> &mdash; 
         Genome fasta file associated with the CRAM file <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.reference_genome_index">
        <b>SingleSampleCnmopsCNVCalling.reference_genome_index</b><br />
        <i>File </i> &mdash; 
         Fai index of the fasta file <br /> 
</p>

### Optional inputs
<p name="SingleSampleCnmopsCNVCalling.input_bam_file">
        <b>SingleSampleCnmopsCNVCalling.input_bam_file</b><br />
        <i>File? </i> &mdash; 
         Input sample BAM/CRAM file. one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.input_bam_file_index">
        <b>SingleSampleCnmopsCNVCalling.input_bam_file_index</b><br />
        <i>File? </i> &mdash; 
         Input sample BAI/CRAI index file <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.input_sample_reads_count">
        <b>SingleSampleCnmopsCNVCalling.input_sample_reads_count</b><br />
        <i>File? </i> &mdash; 
         Inputs sample windowed coverage stored as GenomicRanges object in rds file. can be calculated using cn.mops::getReadCountsFromBAM R function.  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.bed_graph">
        <b>SingleSampleCnmopsCNVCalling.bed_graph</b><br />
        <i>Array[File]? </i> &mdash; 
         Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data).  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.genome_windows">
        <b>SingleSampleCnmopsCNVCalling.genome_windows</b><br />
        <i>File? </i> &mdash; 
         Bed file of the genome binned to equal sized windows similar to the cohort_reads_count_matrix. if bed_graph input is set, this file must be given.  <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.merged_cohort_ploidy_file">
        <b>SingleSampleCnmopsCNVCalling.merged_cohort_ploidy_file</b><br />
        <i>File? </i> &mdash; 
         Cohort ploidy file indicating 1 for male and 2 for female, per sample. The number of lines should be the same as the number of samples in cohort + current_sample. if not given, defaults to 2 for all samples. <br /> 
</p>

### Optional parameters
<p name="SingleSampleCnmopsCNVCalling.chrX_name">
        <b>SingleSampleCnmopsCNVCalling.chrX_name</b><br />
        <i>String? </i> &mdash; 
         The name of the female sex chromosome in the genome. default is: chrX <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.chrY_name">
        <b>SingleSampleCnmopsCNVCalling.chrY_name</b><br />
        <i>String? </i> &mdash; 
         The name of the male sex chromosome in the genome. default is: chrY <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.cap_coverage_override">
        <b>SingleSampleCnmopsCNVCalling.cap_coverage_override</b><br />
        <i>Boolean? </i> &mdash; 
         whether to cap extremely high average coverage windows to 2*cohort's average coverage quantile 99.9% value <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.cnv_lcr_file">
        <b>SingleSampleCnmopsCNVCalling.cnv_lcr_file</b><br />
        <i>File? </i> &mdash; 
         UG-CNV-LCR bed file <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.enable_mod_cnv_override">
        <b>SingleSampleCnmopsCNVCalling.enable_mod_cnv_override</b><br />
        <i>Boolean? </i> &mdash; 
         whether to call moderate cnvs (Fold-Change~1.5 will be tagged as CN2.5 and Fold-Change~0.7 will be tagged as CN1.5). Default is: False <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.skip_figure_generation">
        <b>SingleSampleCnmopsCNVCalling.skip_figure_generation</b><br />
        <i>Boolean? </i> &mdash; 
         Whether to skip figure generation. set true when using reference genome different than hg38.  Default is: False <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.save_hdf_override">
        <b>SingleSampleCnmopsCNVCalling.save_hdf_override</b><br />
        <i>Boolean? </i> &mdash; 
         Whether to save sample reads counts/cohort including sample/cnmops output data in hdf5 format (additionally to RDS format). Default is: False. <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.save_csv_override">
        <b>SingleSampleCnmopsCNVCalling.save_csv_override</b><br />
        <i>Boolean? </i> &mdash; 
         Whether to save sample reads counts/cohort including sample/cnmops output data in csv format (additionally to RDS format). Default is: False. <br /> 
</p>
<p name="SingleSampleCnmopsCNVCalling.preemptible_tries_override">
        <b>SingleSampleCnmopsCNVCalling.preemptible_tries_override</b><br />
        <i>Int? </i> &mdash; 
         Number of preemptible tries,default is: 1 <br /> 
</p>
</details>


## Outputs
<p name="SingleSampleCnmopsCNVCalling.out_sample_reads_count">
        <b>SingleSampleCnmopsCNVCalling.out_sample_reads_count</b><br />
        <i>File</i><br />
        GenomicRanges object of the sample's reads count in rds file format
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_cnvs_bed">
        <b>SingleSampleCnmopsCNVCalling.out_sample_cnvs_bed</b><br />
        <i>Array[File]</i><br />
        Bed file with sample's called CNVs
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_cnvs_vcf">
        <b>SingleSampleCnmopsCNVCalling.out_sample_cnvs_vcf</b><br />
        <i>Array[File]</i><br />
        VCF file with sample's called CNVs
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_cnvs_vcf_index">
        <b>SingleSampleCnmopsCNVCalling.out_sample_cnvs_vcf_index</b><br />
        <i>Array[File]</i><br />
        Index file for the VCF file with sample's called CNVs
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_cnvs_filtered_bed">
        <b>SingleSampleCnmopsCNVCalling.out_sample_cnvs_filtered_bed</b><br />
        <i>Array[File]</i><br />
        Bed file with CNVs filtered by length and overlap with low confidence regions
</p>
<p name="SingleSampleCnmopsCNVCalling.out_coverage_plot_files">
        <b>SingleSampleCnmopsCNVCalling.out_coverage_plot_files</b><br />
        <i>Array[File]</i><br />
        List of coverage figures (file per sample)
</p>
<p name="SingleSampleCnmopsCNVCalling.out_dup_del_plot_files">
        <b>SingleSampleCnmopsCNVCalling.out_dup_del_plot_files</b><br />
        <i>Array[File]</i><br />
        List of duplication/deletion figures (file per sample)
</p>
<p name="SingleSampleCnmopsCNVCalling.out_copy_number_plot_files">
        <b>SingleSampleCnmopsCNVCalling.out_copy_number_plot_files</b><br />
        <i>Array[File]</i><br />
        List of copy number figures (file per sample)
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_merged_bedGraph">
        <b>SingleSampleCnmopsCNVCalling.out_sample_merged_bedGraph</b><br />
        <i>File?</i><br />
        Merged bedGraph file of the sample's coverage
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_reads_count_hdf5">
        <b>SingleSampleCnmopsCNVCalling.out_sample_reads_count_hdf5</b><br />
        <i>File?</i><br />
        GenomicRanges object of the sample's reads count in hdf5 file format
</p>

<hr />

> Generated using WDL AID (1.0.1)
