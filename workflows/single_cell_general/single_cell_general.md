# SingleCell
Create simulated paired end fastq reads from Ultima single-ended CRAM or BAM.

 Intended for reads with a cell barcode, UMI and insert (e.g., 10x, Parse Biosciences, Fluent).

 In addition, runs fastqc on the insert and, optionally, runs FastQC and STAR or STARsolo.

## Inputs

### Required inputs
<p name="SingleCell.input_file">
        <b>SingleCell.input_file</b><br />
        <i>File </i> &mdash;
         Input CRAM or BAM file <br />
</p>
<p name="SingleCell.base_file_name">
        <b>SingleCell.base_file_name</b><br />
        <i>String </i> &mdash;
         Base file name for output files. The output files will be named <tt>[base_file_name]*.fastq.gz</tt> <br />
</p>
<p name="SingleCell.steps">
        <b>SingleCell.steps</b><br />
        <i>TrimAlignSortSteps </i> &mdash;
         The steps to run in the workflow (trim+sort) <br />
</p>
<p name="SingleCell.ref_fastas_cram">
        <b>SingleCell.ref_fastas_cram</b><br />
        <i>Array[File] </i> &mdash;
         Reference fasta files for the CRAM file <br />
</p>
<p name="SingleCell.references">
        <b>SingleCell.references</b><br />
        <i>References </i> &mdash;
         References for the workflow <br />
</p>
<p name="SingleCell.insert_rg">
        <b>SingleCell.insert_rg</b><br />
        <i>String </i> &mdash;
         Read group name for the insert reads, e.g. S1_L001_R2_001 <br />
</p>
<p name="SingleCell.barcode_rg">
        <b>SingleCell.barcode_rg</b><br />
        <i>String </i> &mdash;
         Read group name for the barcode reads, e.g. S1_L001_R1_001 <br />
</p>
<p name="SingleCell.qc_thresholds">
        <b>SingleCell.qc_thresholds</b><br />
        <i>SingleCellQcThresholds </i> &mdash;
         Thresholds for the single cell qc <br />
</p>

### Required parameters
<p name="SingleCell.no_address">
        <b>SingleCell.no_address</b><br />
        <i>Boolean </i> &mdash;
         Should the instances used be without external IP address. Allows for more parallelization, but not supported with Dockerhub dockers Default: true <br />
</p>
<p name="SingleCell.preemptible_tries">
        <b>SingleCell.preemptible_tries</b><br />
        <i>Int </i> &mdash;
         Number of preemptible tries <br />
</p>
<p name="SingleCell.cpu">
        <b>SingleCell.cpu</b><br />
        <i>Int </i> &mdash;
         Number of CPUs to use (for Trimmer, conversion to fastq, and alignment tasks) <br />
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="SingleCell.trimmer_parameters">
        <b>SingleCell.trimmer_parameters</b><br />
        <i>TrimmerParameters &mdash; Default: None</i><br />
        Parameters for Trimmer task.  See input template
</p>
<p name="SingleCell.sorter_params">
        <b>SingleCell.sorter_params</b><br />
        <i>SorterParams &mdash; Default: None</i><br />
        Parameters for Sorter task.  See input template
</p>
<p name="SingleCell.additional_rg">
        <b>SingleCell.additional_rg</b><br />
        <i>String? &mdash; Default: None</i><br />
        Additional read group name for a third read, e.g. S1_L001_I1_001
</p>
<p name="SingleCell.downstream_analysis">
        <b>SingleCell.downstream_analysis</b><br />
        <i>String? &mdash; Default: None</i><br />
        Can be either star_solo, star, or undefined (default)
</p>
<p name="SingleCell.star_solo_params">
        <b>SingleCell.star_solo_params</b><br />
        <i>StarSoloParams? &mdash; Default: None</i><br />
        Parameters for running the StarSolo task.  See input template
</p>
<p name="SingleCell.genome_generate_params">
        <b>SingleCell.genome_generate_params</b><br />
        <i>StarGenomeGenerateParams? &mdash; Default: None</i><br />
        Parameters that can be passed to the StarSolo or StarAlignment task, for formatting a STAR genome.  See input template
</p>
<p name="SingleCell.star_align_gtf_override">
        <b>SingleCell.star_align_gtf_override</b><br />
        <i>File? &mdash; Default: None</i><br />
        The gtf to use in the StarAlignment task.  See input template
</p>

### Optional inputs
<p name="SingleCell.create_md5_checksum_outputs">
        <b>SingleCell.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
<p name="SingleCell.TrimAlignSort.aligner">
        <b>SingleCell.TrimAlignSort.aligner</b><br />
        <i>String? </i> &mdash;
         Aligner to be used. Options are: ua, ua-meth, star. Mandatory if align step is selected. <br />
</p>
<p name="SingleCell.TrimAlignSort.ua_parameters">
        <b>SingleCell.TrimAlignSort.ua_parameters</b><br />
        <i>UaParameters? </i> &mdash;
         Parameters for the UA aligner. Mandatory if aligner is ua. <br />
</p>
<p name="SingleCell.TrimAlignSort.ua_meth_parameters">
        <b>SingleCell.TrimAlignSort.ua_meth_parameters</b><br />
        <i>UaMethParameters? </i> &mdash;
         Parameters for the UA meth aligner. Mandatory if aligner is ua-meth. <br />
</p>
<p name="SingleCell.TrimAlignSort.star_genome">
        <b>SingleCell.TrimAlignSort.star_genome</b><br />
        <i>File? </i> &mdash;
         Star genome file. If aligner is star, supllay either star genome file or generate new genome index. <br />
</p>
<p name="SingleCell.TrimAlignSort.star_genome_generate_params">
        <b>SingleCell.TrimAlignSort.star_genome_generate_params</b><br />
        <i>StarGenomeGenerateParams? </i> &mdash;
         Parameters for generating the star genome. Mandatory if aligner is star and not given star genome file. <br />
</p>
<p name="SingleCell.TrimAlignSort.star_align_extra_args">
        <b>SingleCell.TrimAlignSort.star_align_extra_args</b><br />
        <i>String? </i> &mdash;
         Extra arguments for the STAR aligner. <br />
</p>
<p name="SingleCell.TrimAlignSort.star_align_gtf_override">
        <b>SingleCell.TrimAlignSort.star_align_gtf_override</b><br />
        <i>File? </i> &mdash;
         GTF file to be used for STAR aligner. <br />
</p>
<p name="SingleCell.TrimAlignSort.create_md5_checksum_outputs">
        <b>SingleCell.TrimAlignSort.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
<p name="SingleCell.StarAlignSubSample.genome_generate_params">
        <b>SingleCell.StarAlignSubSample.genome_generate_params</b><br />
        <i>StarGenomeGenerateParams? </i> &mdash;
         Parameters for generating the reference genome. <br />
</p>
<p name="SingleCell.StarAlignSubSample.star_align_gtf_override">
        <b>SingleCell.StarAlignSubSample.star_align_gtf_override</b><br />
        <i>File? </i> &mdash;
         Override the GTF file used for STAR alignment. <br />
</p>

### Optional parameters
<p name="SingleCell.star_align_extra_args">
        <b>SingleCell.star_align_extra_args</b><br />
        <i>String? </i> &mdash;
         Extra parameters to pass to the StarAlignment task.  See input template <br />
</p>
</details>


## Outputs
<p name="SingleCell.trimmer_stats_output">
        <b>SingleCell.trimmer_stats_output</b><br />
        <i>File?</i><br />
        Trimmer output statistics
</p>
<p name="SingleCell.trimmer_failure_codes_csv">
        <b>SingleCell.trimmer_failure_codes_csv</b><br />
        <i>File?</i><br />
        Trimmer failure codes csv
</p>
<p name="SingleCell.trimmer_histogram">
        <b>SingleCell.trimmer_histogram</b><br />
        <i>Array[File?]?</i><br />
        Trimmer histograms
</p>
<p name="SingleCell.trimmer_histogram_extra">
        <b>SingleCell.trimmer_histogram_extra</b><br />
        <i>Array[File?]?</i><br />
        Trimmer extra histograms
</p>
<p name="SingleCell.sorter_stats_csv">
        <b>SingleCell.sorter_stats_csv</b><br />
        <i>Array[File?]?</i><br />
        Sorter statistics csv
</p>
<p name="SingleCell.sorter_stats_json">
        <b>SingleCell.sorter_stats_json</b><br />
        <i>Array[File?]?</i><br />
        Sorter statistics json
</p>
<p name="SingleCell.unmatched_cram">
        <b>SingleCell.unmatched_cram</b><br />
        <i>File?</i><br />
        Unmatched cram file output from sorter
</p>
<p name="SingleCell.output_barcodes_fastq">
        <b>SingleCell.output_barcodes_fastq</b><br />
        <i>File</i><br />
        The fastq with the barcodes portion of the read
</p>
<p name="SingleCell.output_insert_fastq">
        <b>SingleCell.output_insert_fastq</b><br />
        <i>File</i><br />
        The fastq with the insert portion of the read
</p>
<p name="SingleCell.output_additional_fastq">
        <b>SingleCell.output_additional_fastq</b><br />
        <i>File?</i><br />
        The additional fastq, if an additional_rg was specified
</p>
<p name="SingleCell.report_html">
        <b>SingleCell.report_html</b><br />
        <i>File</i><br />
        The report from the single cell qc
</p>
<p name="SingleCell.application_qc_h5">
        <b>SingleCell.application_qc_h5</b><br />
        <i>File</i><br />
        Single Cell application QC h5 file
</p>
<p name="SingleCell.aggregated_metrics_json">
        <b>SingleCell.aggregated_metrics_json</b><br />
        <i>File</i><br />
        Single Cell application QC json file
</p>
<p name="SingleCell.star_bam">
        <b>SingleCell.star_bam</b><br />
        <i>File?</i><br />
        The Star-aligned bam.  Created only if STAR is run.
</p>
<p name="SingleCell.star_reads_per_gene_file">
        <b>SingleCell.star_reads_per_gene_file</b><br />
        <i>File?</i><br />
        The Star reads per gene file
</p>
<p name="SingleCell.star_stats">
        <b>SingleCell.star_stats</b><br />
        <i>File?</i><br />
        The Star output statistics
</p>
<p name="SingleCell.md5_checksums_json">
        <b>SingleCell.md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files
</p>
<p name="SingleCell.star_solo_outputs">
        <b>SingleCell.star_solo_outputs</b><br />
        <i>StarSoloOutputs?</i><br />
        The outputs from the StarSolo workflow.  Created only if STARsolo is run
</p>

<hr />

> Generated using WDL AID (1.0.1)