# SingleSampleCnmopsCNVCalling
Runs single sample germline CNV calling workflow based on \<a href="https://bioconductor.org/packages/release/bioc/html/cn.mops.html"\>cn.mops</a>

The pipeline uses a given cohort's coverage profile for normalization.

The pipeline can recieve one of the following options as input:

&nbsp;&nbsp;1. input CRAM/BAM file.

&nbsp;&nbsp;2. rds file which stores a GenomicRanges object with coverage collected in the same windows as the given cohort.

&nbsp;&nbsp;3. BedGraph holding the coverage per location.

The pipeline calls CNVs for the given sample and filters them by length (>10,000b) and overlap with UG-CNV-LCR.

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
<p name="SingleSampleCnmopsCNVCalling.cnv_lcr_file">
        <b>SingleSampleCnmopsCNVCalling.cnv_lcr_file</b><br />
        <i>File </i> &mdash; 
         UG-CNV-LCR bed file <br /> 
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
        <i>File? </i> &mdash; 
         Previously calculated input bedGraph holding the coverage per base (outputs with the sequencing data).  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br /> 
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
        <i>File</i><br />
        Bed file with sample's called CNVs
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_cnvs_filtered_bed">
        <b>SingleSampleCnmopsCNVCalling.out_sample_cnvs_filtered_bed</b><br />
        <i>File</i><br />
        Bed file with CNVs filtered by length and overlap with low confidence regions
</p>
<p name="SingleSampleCnmopsCNVCalling.out_sample_reads_count_hdf5">
        <b>SingleSampleCnmopsCNVCalling.out_sample_reads_count_hdf5</b><br />
        <i>File?</i><br />
        GenomicRanges object of the sample's reads count in hdf5 file format
</p>

<hr />

> Generated using WDL AID (1.0.0)
