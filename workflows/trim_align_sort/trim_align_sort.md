# TrimAlignSort
Pipeline for trimming, aligning and sorting Ultima data in a fast, cost-effective, and easier-to-maintain way, similiar to the way it is done on-tool. It can be used to run the entire pipeline or just a subset of the steps.

## Inputs

### Required inputs
<p name="TrimAlignSort.input_cram_bam_list">
        <b>TrimAlignSort.input_cram_bam_list</b><br />
        <i>Array[File] </i> &mdash;
         List of input cram or bam files to be processed. <br />
</p>
<p name="TrimAlignSort.ref_fastas_cram">
        <b>TrimAlignSort.ref_fastas_cram</b><br />
        <i>Array[File] </i> &mdash;
         List of references for CreateReferenceCache task. <br />
</p>
<p name="TrimAlignSort.base_file_name">
        <b>TrimAlignSort.base_file_name</b><br />
        <i>String </i> &mdash;
         Base name for the output files. <br />
</p>
<p name="TrimAlignSort.steps">
        <b>TrimAlignSort.steps</b><br />
        <i>TrimAlignSortSteps </i> &mdash;
         Steps to be executed in the pipeline.Options are: trim, align, sort <br />
</p>
<p name="TrimAlignSort.references">
        <b>TrimAlignSort.references</b><br />
        <i>References </i> &mdash;
         References for merging inputs into one file, alignment, and sorting. <br />
</p>
<p name="TrimAlignSort.cpu">
        <b>TrimAlignSort.cpu</b><br />
        <i>Int </i> &mdash;
         Number of cpus to be used for the tasks. <br />
</p>

### Optional inputs
<p name="TrimAlignSort.trimmer_parameters">
        <b>TrimAlignSort.trimmer_parameters</b><br />
        <i>TrimmerParameters? </i> &mdash;
         Parameters for the trimmer task. Mandatory if trim step is selected. <br />
</p>
<p name="TrimAlignSort.aligner">
        <b>TrimAlignSort.aligner</b><br />
        <i>String? </i> &mdash;
         Aligner to be used. Options are: ua, ua-meth, star. Mandatory if align step is selected. <br />
</p>
<p name="TrimAlignSort.ua_parameters">
        <b>TrimAlignSort.ua_parameters</b><br />
        <i>UaParameters? </i> &mdash;
         Parameters for the UA aligner. Mandatory if aligner is ua. <br />
</p>
<p name="TrimAlignSort.ua_meth_parameters">
        <b>TrimAlignSort.ua_meth_parameters</b><br />
        <i>UaMethParameters? </i> &mdash;
         Parameters for the UA meth aligner. Mandatory if aligner is ua-meth. <br />
</p>
<p name="TrimAlignSort.star_genome">
        <b>TrimAlignSort.star_genome</b><br />
        <i>File? </i> &mdash;
         Star genome file. If aligner is star, supllay either star genome file or generate new genome index. <br />
</p>
<p name="TrimAlignSort.star_genome_generate_params">
        <b>TrimAlignSort.star_genome_generate_params</b><br />
        <i>StarGenomeGenerateParams? </i> &mdash;
         Parameters for generating the star genome. Mandatory if aligner is star and not given star genome file. <br />
</p>
<p name="TrimAlignSort.star_align_extra_args">
        <b>TrimAlignSort.star_align_extra_args</b><br />
        <i>String? </i> &mdash;
         Extra arguments for the STAR aligner. <br />
</p>
<p name="TrimAlignSort.star_align_gtf_override">
        <b>TrimAlignSort.star_align_gtf_override</b><br />
        <i>File? </i> &mdash;
         GTF file to be used for STAR aligner. <br />
</p>
<p name="TrimAlignSort.sorter_params">
        <b>TrimAlignSort.sorter_params</b><br />
        <i>SorterParams? </i> &mdash;
         Parameters for the sorter task. Mandatory if sort step is selected. <br />
</p>
<p name="TrimAlignSort.create_md5_checksum_outputs">
        <b>TrimAlignSort.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
</details>


## Outputs
<p name="TrimAlignSort.output_cram_bam">
        <b>TrimAlignSort.output_cram_bam</b><br />
        <i>File?</i><br />
        Output file after the pipeline is executed.
</p>
<p name="TrimAlignSort.output_cram_bam_index">
        <b>TrimAlignSort.output_cram_bam_index</b><br />
        <i>File?</i><br />
        Index file for the output cram file.
</p>
<p name="TrimAlignSort.fastq_file">
        <b>TrimAlignSort.fastq_file</b><br />
        <i>File?</i><br />
        Sorter output in Fastq format.
</p>
<p name="TrimAlignSort.sorter_stats_csv">
        <b>TrimAlignSort.sorter_stats_csv</b><br />
        <i>File?</i><br />
        Sorter stats in csv format.
</p>
<p name="TrimAlignSort.sorter_stats_json">
        <b>TrimAlignSort.sorter_stats_json</b><br />
        <i>File?</i><br />
        Sorter stats in json format.
</p>
<p name="TrimAlignSort.unmatched_cram">
        <b>TrimAlignSort.unmatched_cram</b><br />
        <i>File?</i><br />
        Unmatched cram file output from sorter (if defined).
</p>
<p name="TrimAlignSort.unmatched_sorter_stats_csv">
        <b>TrimAlignSort.unmatched_sorter_stats_csv</b><br />
        <i>File?</i><br />
        Unmatched output cram files sorter stats in csv format (if defined).
</p>
<p name="TrimAlignSort.unmatched_sorter_stats_json">
        <b>TrimAlignSort.unmatched_sorter_stats_json</b><br />
        <i>File?</i><br />
        Unmatched output cram files sorter stats in json format (if defined).
</p>
<p name="TrimAlignSort.bedgraph_mapq0">
        <b>TrimAlignSort.bedgraph_mapq0</b><br />
        <i>File?</i><br />
        Bedgraph mapq0 file.
</p>
<p name="TrimAlignSort.bedgraph_mapq1">
        <b>TrimAlignSort.bedgraph_mapq1</b><br />
        <i>File?</i><br />
        Bedgraph mapq1 file.
</p>
<p name="TrimAlignSort.downsampling_seed">
        <b>TrimAlignSort.downsampling_seed</b><br />
        <i>File?</i><br />
        Seed used for downsampling (if downsampling_frac is set).
</p>
<p name="TrimAlignSort.aggregated_metrics_h5">
        <b>TrimAlignSort.aggregated_metrics_h5</b><br />
        <i>File?</i><br />
        Aggregated metrics in h5 format.
</p>
<p name="TrimAlignSort.aggregated_metrics_json">
        <b>TrimAlignSort.aggregated_metrics_json</b><br />
        <i>File?</i><br />
        Aggregated metrics in json format.
</p>
<p name="TrimAlignSort.report_html">
        <b>TrimAlignSort.report_html</b><br />
        <i>File?</i><br />
        Report in html format.
</p>
<p name="TrimAlignSort.fastq_file_list">
        <b>TrimAlignSort.fastq_file_list</b><br />
        <i>Array[File?]?</i><br />
        Sorter output in Fastq format (allow multiple files for applications like single-cell).
</p>
<p name="TrimAlignSort.sub_sampled_output">
        <b>TrimAlignSort.sub_sampled_output</b><br />
        <i>Array[File?]?</i><br />
        Sub-sampling files as part of sorter output.
</p>
<p name="TrimAlignSort.output_cram_bam_list">
        <b>TrimAlignSort.output_cram_bam_list</b><br />
        <i>Array[File?]?</i><br />
        Output file list after the pipeline is executed in multiple outputs mode
</p>
<p name="TrimAlignSort.output_cram_bam_index_list">
        <b>TrimAlignSort.output_cram_bam_index_list</b><br />
        <i>Array[File?]?</i><br />
        Index files for the output cram files (when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.sorter_stats_csv_list">
        <b>TrimAlignSort.sorter_stats_csv_list</b><br />
        <i>Array[File?]?</i><br />
        Sorter stats in csv format (when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.sorter_stats_json_list">
        <b>TrimAlignSort.sorter_stats_json_list</b><br />
        <i>Array[File?]?</i><br />
        Sorter stats in json format (when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.unmatched_cram_list">
        <b>TrimAlignSort.unmatched_cram_list</b><br />
        <i>Array[File?]?</i><br />
        Unmatched cram file output from sorter (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.unmatched_sorter_stats_csv_list">
        <b>TrimAlignSort.unmatched_sorter_stats_csv_list</b><br />
        <i>Array[File?]?</i><br />
        Unmatched output cram files sorter stats in csv format (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.unmatched_sorter_stats_json_list">
        <b>TrimAlignSort.unmatched_sorter_stats_json_list</b><br />
        <i>Array[File?]?</i><br />
        Unmatched output cram files sorter stats in json format (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.bedgraph_mapq0_list">
        <b>TrimAlignSort.bedgraph_mapq0_list</b><br />
        <i>Array[File?]?</i><br />
        Bedgraph mapq0 files (when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.bedgraph_mapq1_list">
        <b>TrimAlignSort.bedgraph_mapq1_list</b><br />
        <i>Array[File?]?</i><br />
        Bedgraph mapq1 files (when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.aggregated_metrics_h5_list">
        <b>TrimAlignSort.aggregated_metrics_h5_list</b><br />
        <i>Array[File?]?</i><br />
        Aggregated metrics in h5 format (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.aggregated_metrics_json_list">
        <b>TrimAlignSort.aggregated_metrics_json_list</b><br />
        <i>Array[File?]?</i><br />
        Aggregated metrics in json format (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.report_html_list">
        <b>TrimAlignSort.report_html_list</b><br />
        <i>Array[File?]?</i><br />
        Reports in html format (if defined, when running in multiple outputs mode).
</p>
<p name="TrimAlignSort.trimmer_stats">
        <b>TrimAlignSort.trimmer_stats</b><br />
        <i>File?</i><br />
        Trimmer stats file.
</p>
<p name="TrimAlignSort.trimmer_failure_codes_csv">
        <b>TrimAlignSort.trimmer_failure_codes_csv</b><br />
        <i>File?</i><br />
        Trimmer failure codes in csv format.
</p>
<p name="TrimAlignSort.trimmer_histogram">
        <b>TrimAlignSort.trimmer_histogram</b><br />
        <i>Array[File?]?</i><br />
        Trimmer histogram files.
</p>
<p name="TrimAlignSort.trimmer_histogram_extra">
        <b>TrimAlignSort.trimmer_histogram_extra</b><br />
        <i>Array[File?]?</i><br />
        Trimmer histogram extra files.
</p>
<p name="TrimAlignSort.ua_stats_jsons">
        <b>TrimAlignSort.ua_stats_jsons</b><br />
        <i>Array[File]?</i><br />
        Statistics file in json format from UA.
</p>
<p name="TrimAlignSort.ua_index_build">
        <b>TrimAlignSort.ua_index_build</b><br />
        <i>File?</i><br />
        UA index files if UA alignment is used and index is built as part of the workflow.
</p>
<p name="TrimAlignSort.ua_index_c2t_build">
        <b>TrimAlignSort.ua_index_c2t_build</b><br />
        <i>File?</i><br />
        UA-meth index C2T files if UA-meth alignment is used and index is built as part of the workflow.
</p>
<p name="TrimAlignSort.ua_index_g2a_build">
        <b>TrimAlignSort.ua_index_g2a_build</b><br />
        <i>File?</i><br />
        UA-meth index G2A files if UA-meth alignment is used and index is built as part of the workflow.
</p>
<p name="TrimAlignSort.align_star_reads_per_gene_file">
        <b>TrimAlignSort.align_star_reads_per_gene_file</b><br />
        <i>File?</i><br />
        STAR reads per gene file.
</p>
<p name="TrimAlignSort.align_star_stats">
        <b>TrimAlignSort.align_star_stats</b><br />
        <i>File?</i><br />
        STAR stats file.
</p>
<p name="TrimAlignSort.md5_checksums_json">
        <b>TrimAlignSort.md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files
</p>

<hr />

> Generated using WDL AID (1.0.1)