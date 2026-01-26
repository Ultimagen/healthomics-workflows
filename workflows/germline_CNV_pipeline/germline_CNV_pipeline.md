# GermlineCNVPipeline
Runs: <br>1. single sample germline CNV calling workflow based on [cn.mops](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)<br>2. cnvpytor workflow<br> 3. combines results, verifies them using split reads and jump alignments<br>4. Applies ML model to estimate quality of the CNV</b>

## Inputs

### Required inputs
<p name="GermlineCNVPipeline.base_file_name">
        <b>GermlineCNVPipeline.base_file_name</b><br />
        <i>String </i> &mdash;
         Sample name <br />
</p>
<p name="GermlineCNVPipeline.input_bam_file">
        <b>GermlineCNVPipeline.input_bam_file</b><br />
        <i>File </i> &mdash;
         Input sample BAM/CRAM file <br />
</p>
<p name="GermlineCNVPipeline.input_bam_file_index">
        <b>GermlineCNVPipeline.input_bam_file_index</b><br />
        <i>File </i> &mdash;
         Input sample BAI/CRAI index file <br />
</p>
<p name="GermlineCNVPipeline.bed_graph">
        <b>GermlineCNVPipeline.bed_graph</b><br />
        <i>Array[File] </i> &mdash;
         Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data). <br />
</p>
<p name="GermlineCNVPipeline.genome_windows">
        <b>GermlineCNVPipeline.genome_windows</b><br />
        <i>File </i> &mdash;
         Bed file of the genome binned to equal sized windows similar to the cohort_reads_count_matrix. <br />
</p>
<p name="GermlineCNVPipeline.cohort_reads_count_matrix">
        <b>GermlineCNVPipeline.cohort_reads_count_matrix</b><br />
        <i>File </i> &mdash;
         GenomicRanges object of the cohort reads count matrix in rds file format. default cohort can be found in the template. <br />
</p>
<p name="GermlineCNVPipeline.ploidy_file">
        <b>GermlineCNVPipeline.ploidy_file</b><br />
        <i>File </i> &mdash;
         X chromosome ploidy of the cohort and the additional sample. Each sample is represented on a number on a separate row. Ploidy of the default cohort can be found in the template. The last row corresponds to the sample being called <br />
</p>

### Required parameters
<p name="GermlineCNVPipeline.ref_seq_names">
        <b>GermlineCNVPipeline.ref_seq_names</b><br />
        <i>Array[String] </i> &mdash;
         Chromosome names for which coverage will be calculated <br />
</p>
<p name="GermlineCNVPipeline.skip_filtering">
        <b>GermlineCNVPipeline.skip_filtering</b><br />
        <i>Boolean </i> &mdash;
         Whether to skip CNV filtering step, default is False <br />
</p>
<p name="GermlineCNVPipeline.cushion_size">
        <b>GermlineCNVPipeline.cushion_size</b><br />
        <i>Int </i> &mdash;
         Cushion size around CNV breakpoints for split-read analysis and jump alignment analysis <br />
</p>

### Required references
<p name="GermlineCNVPipeline.reference_genome">
        <b>GermlineCNVPipeline.reference_genome</b><br />
        <i>File </i> &mdash;
         Genome fasta file associated with the CRAM file <br />
</p>
<p name="GermlineCNVPipeline.reference_genome_index">
        <b>GermlineCNVPipeline.reference_genome_index</b><br />
        <i>File </i> &mdash;
         Fai index of the fasta file <br />
</p>

### Optional inputs
<p name="GermlineCNVPipeline.ug_cnv_lcr_file">
        <b>GermlineCNVPipeline.ug_cnv_lcr_file</b><br />
        <i>File? </i> &mdash;
         UG-CNV-LCR bed file <br />
</p>
<p name="GermlineCNVPipeline.filtering_model">
        <b>GermlineCNVPipeline.filtering_model</b><br />
        <i>File? </i> &mdash;
         CNV filtering model, default in template, calls are not filtered if not provided <br />
</p>
<p name="GermlineCNVPipeline.create_md5_checksum_outputs">
        <b>GermlineCNVPipeline.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.input_bam_file">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.input_bam_file</b><br />
        <i>File? </i> &mdash;
         Input sample BAM/CRAM file. one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.input_bam_file_index">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.input_bam_file_index</b><br />
        <i>File? </i> &mdash;
         Input sample BAI/CRAI index file <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.input_sample_reads_count">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.input_sample_reads_count</b><br />
        <i>File? </i> &mdash;
         Inputs sample windowed coverage stored as GenomicRanges object in rds file. can be calculated using cn.mops::getReadCountsFromBAM R function.  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set <br />
</p>

### Optional parameters
<p name="GermlineCNVPipeline.filtering_model_decision_threshold">
        <b>GermlineCNVPipeline.filtering_model_decision_threshold</b><br />
        <i>Int? </i> &mdash;
         Decision threshold for the filtering model, default is set in template. Lower- less stringent, Higher- more stringent <br />
</p>
<p name="GermlineCNVPipeline.skip_figure_generation">
        <b>GermlineCNVPipeline.skip_figure_generation</b><br />
        <i>Boolean? </i> &mdash;
         Skip CNV calls figure generation. please set to True if reference genome is not hg38. Default is: False <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.chrX_name">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.chrX_name</b><br />
        <i>String? </i> &mdash;
         The name of the female sex chromosome in the genome. default is: chrX <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.chrY_name">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.chrY_name</b><br />
        <i>String? </i> &mdash;
         The name of the male sex chromosome in the genome. default is: chrY <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.cap_coverage_override">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.cap_coverage_override</b><br />
        <i>Boolean? </i> &mdash;
         whether to cap extremely high average coverage windows to 2*cohort's average coverage quantile 99.9% value <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.save_hdf_override">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.save_hdf_override</b><br />
        <i>Boolean? </i> &mdash;
         Whether to save sample reads counts/cohort including sample/cnmops output data in hdf5 format (additionally to RDS format). Default is: False. <br />
</p>
<p name="GermlineCNVPipeline.CnmopsCNVCalling.save_csv_override">
        <b>GermlineCNVPipeline.CnmopsCNVCalling.save_csv_override</b><br />
        <i>Boolean? </i> &mdash;
         Whether to save sample reads counts/cohort including sample/cnmops output data in csv format (additionally to RDS format). Default is: False. <br />
</p>
</details>


## Outputs
<p name="GermlineCNVPipeline.cnmops_cnv_calls_bed">
        <b>GermlineCNVPipeline.cnmops_cnv_calls_bed</b><br />
        <i>File</i><br />
        CNMOPS CNV calls in bed format
</p>
<p name="GermlineCNVPipeline.cnmops_cnv_calls_vcf">
        <b>GermlineCNVPipeline.cnmops_cnv_calls_vcf</b><br />
        <i>File</i><br />
        CNMOPS CNV calls in VCF format
</p>
<p name="GermlineCNVPipeline.cnmops_cnv_calls_vcf_index">
        <b>GermlineCNVPipeline.cnmops_cnv_calls_vcf_index</b><br />
        <i>File</i><br />
        Index file for the CNMOPS CNV calls VCF
</p>
<p name="GermlineCNVPipeline.cnvpytor_cnv_calls_bed">
        <b>GermlineCNVPipeline.cnvpytor_cnv_calls_bed</b><br />
        <i>File</i><br />
        CNVpytor CNV calls in bed format
</p>
<p name="GermlineCNVPipeline.cnvpytor_cnv_calls_vcf">
        <b>GermlineCNVPipeline.cnvpytor_cnv_calls_vcf</b><br />
        <i>File</i><br />
        CNVpytor CNV calls in VCF format
</p>
<p name="GermlineCNVPipeline.cnvpytor_cnv_calls_vcf_index">
        <b>GermlineCNVPipeline.cnvpytor_cnv_calls_vcf_index</b><br />
        <i>File</i><br />
        Index file for the CNVpytor CNV calls VCF
</p>
<p name="GermlineCNVPipeline.combined_cnv_calls_bed">
        <b>GermlineCNVPipeline.combined_cnv_calls_bed</b><br />
        <i>File</i><br />
        Final (combined) CNV calls in bed format
</p>
<p name="GermlineCNVPipeline.combined_cnv_calls_bed_vcf">
        <b>GermlineCNVPipeline.combined_cnv_calls_bed_vcf</b><br />
        <i>File</i><br />
        Combined CNV calls in vcf format
</p>
<p name="GermlineCNVPipeline.combined_cnv_calls_bed_vcf_index">
        <b>GermlineCNVPipeline.combined_cnv_calls_bed_vcf_index</b><br />
        <i>File</i><br />
        Index of the combined CNV calls in vcf format
</p>
<p name="GermlineCNVPipeline.combine_read_evidence">
        <b>GermlineCNVPipeline.combine_read_evidence</b><br />
        <i>File</i><br />
        BAM file with read evidence supporting combined CNV calls
</p>
<p name="GermlineCNVPipeline.combine_read_evidence_index">
        <b>GermlineCNVPipeline.combine_read_evidence_index</b><br />
        <i>File</i><br />
        Index file for the BAM with read evidence supporting combined CNV calls
</p>
<p name="GermlineCNVPipeline.combine_read_scores_csv">
        <b>GermlineCNVPipeline.combine_read_scores_csv</b><br />
        <i>File</i><br />
        CSV file with jalign scores for each read
</p>
<p name="GermlineCNVPipeline.md5_checksums_json">
        <b>GermlineCNVPipeline.md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files
</p>

<hr />

> Generated using WDL AID (1.0.1)