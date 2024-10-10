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
        <i>TrimmerParameters? &mdash; Default: None</i><br />
        Parameters for Trimmer task.  See input template
</p>
<p name="SingleCell.trimmer_stats">
        <b>SingleCell.trimmer_stats</b><br />
        <i>File? &mdash; Default: None</i><br />
        If trimmer was run outside of this workflow, the stats can still be combined in the final report
</p>
<p name="SingleCell.barcode_fastq_file_suffix">
        <b>SingleCell.barcode_fastq_file_suffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Suffix to add to the name of the file with the barcode reads (in case downstream software has name requirements)
</p>
<p name="SingleCell.insert_fastq_file_suffix">
        <b>SingleCell.insert_fastq_file_suffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Suffix to add to the name of the file with the insert reads (in case downstream software has name requirements)
</p>
<p name="SingleCell.barcode_fastq_header_suffix">
        <b>SingleCell.barcode_fastq_header_suffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Suffix to add to the end of the fastq headers for the barcode reads (in case downstream software has fastq header requirements)
</p>
<p name="SingleCell.insert_fastq_header_suffix">
        <b>SingleCell.insert_fastq_header_suffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Suffix to add to the end of the fastq headers for the insert reads (in case downstream software has fastq header requirements)
</p>
<p name="SingleCell.fastqc_adapter">
        <b>SingleCell.fastqc_adapter</b><br />
        <i>String? &mdash; Default: None</i><br />
        Adapter that can be passed to fastqc with the --adapters options
</p>
<p name="SingleCell.fastqc_limits">
        <b>SingleCell.fastqc_limits</b><br />
        <i>File? &mdash; Default: None</i><br />
        Adapter that can be passed to fastqc with the --limits option
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

### Optional parameters
<p name="SingleCell.demux_extra_args">
        <b>SingleCell.demux_extra_args</b><br />
        <i>String? </i> &mdash; 
         Extra parameters to pass to SortTasks.Demux , when converting to fastq <br /> 
</p>
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
<p name="SingleCell.fastqc_reports">
        <b>SingleCell.fastqc_reports</b><br />
        <i>Array[File]</i><br />
        Fastqc output report for insert
</p>
<p name="SingleCell.combined_statistics">
        <b>SingleCell.combined_statistics</b><br />
        <i>File</i><br />
        A csv with the trimming and alignment statistics
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
<p name="SingleCell.star_solo_outputs">
        <b>SingleCell.star_solo_outputs</b><br />
        <i>StarSoloOutputs?</i><br />
        The outputs from the StarSolo workflow.  Created only if STARsolo is run
</p>

<hr />

> Generated using WDL AID (1.0.0)
