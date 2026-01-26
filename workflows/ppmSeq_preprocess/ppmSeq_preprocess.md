# ppmSeqPreprocess
The ppmSeq Preprocess pipeline is designed to process untrimmed ppmSeq sequencing data through a complete preprocessing workflow including adapter trimming, alignment to reference genome, sorting, duplicate marking, and QC report generation. This analysis is generally done on the UG sequencer, this pipeline is intended for cases where it did not happen or was improperly configured. For more details on ppmSeq, see https://www.ultimagenomics.com/products/ppmseq-tm/ and https://www.biorxiv.org/content/10.1101/2025.08.11.669689v1. 

The following input templates are available for different input data: 

1) `ppmSeq_preprocess_template-ppmSeq.json` | Use this template for ppmSeq data. The input CRAM file should NOT be trimmed. If it contains the ppmSeq tags (e.g. st, et), it was trimmed. 

2) `ppmSeq_preprocess_template-ppmSeq_legacy_v5.json` | Use this template for LEGACY v5 ppmSeq data. This is an older version of the ppmSeq adapters, generally not available since 2024. The input CRAM file should NOT be trimmed. If it contains the ppmSeq tags (e.g. as, ts), it was trimmed. 

3) `ppmSeq_preprocess_template-ppmSeq_post_native_adapter_trimming.json` | Use this template only for the case where the UG native adapters were trimmed, but not the ppmSeq adapters and loop. This generally happens if the application_type is configured to be 'native' instead of 'ppmSeq', can be verified by the presence of an 'a3' tag in some reads but the absence of ppmSeq tags (e.g. st, et). 

## Inputs

### Required inputs
<p name="ppmSeqPreprocess.input_cram_bam_list">
        <b>ppmSeqPreprocess.input_cram_bam_list</b><br />
        <i>Array[File] </i> &mdash;
         Input CRAM or BAM file list <br />
</p>
<p name="ppmSeqPreprocess.base_file_name">
        <b>ppmSeqPreprocess.base_file_name</b><br />
        <i>String </i> &mdash;
         Base file name for output files. <br />
</p>
<p name="ppmSeqPreprocess.adapter_version">
        <b>ppmSeqPreprocess.adapter_version</b><br />
        <i>String </i> &mdash;
         ppmSeq adapter version <br />
</p>

### Required parameters
<p name="ppmSeqPreprocess.trimmer_parameters">
        <b>ppmSeqPreprocess.trimmer_parameters</b><br />
        <i>TrimmerParameters </i> &mdash;
         Trimmer parameters <br />
</p>
<p name="ppmSeqPreprocess.ua_parameters">
        <b>ppmSeqPreprocess.ua_parameters</b><br />
        <i>UaParameters </i> &mdash;
         UA alignment parameters <br />
</p>
<p name="ppmSeqPreprocess.sorter_params">
        <b>ppmSeqPreprocess.sorter_params</b><br />
        <i>SorterParams </i> &mdash;
         Sorter parameters <br />
</p>

### Required references
<p name="ppmSeqPreprocess.ref_fastas_cram">
        <b>ppmSeqPreprocess.ref_fastas_cram</b><br />
        <i>Array[File] </i> &mdash;
         Reference fasta file for cache tarball <br />
</p>
<p name="ppmSeqPreprocess.references">
        <b>ppmSeqPreprocess.references</b><br />
        <i>References </i> &mdash;
         Reference files <br />
</p>

### Optional inputs
<p name="ppmSeqPreprocess.ppmSeq_analysis_extra_args">
        <b>ppmSeqPreprocess.ppmSeq_analysis_extra_args</b><br />
        <i>String? </i> &mdash;
         Extra arguments for ppmSeq analysis <br />
</p>
<p name="ppmSeqPreprocess.create_md5_checksum_outputs">
        <b>ppmSeqPreprocess.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.ua_meth_parameters">
        <b>ppmSeqPreprocess.TrimAlignSort.ua_meth_parameters</b><br />
        <i>UaMethParameters? </i> &mdash;
         Parameters for the UA meth aligner. Mandatory if aligner is ua-meth. <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.star_genome">
        <b>ppmSeqPreprocess.TrimAlignSort.star_genome</b><br />
        <i>File? </i> &mdash;
         Star genome file. If aligner is star, supllay either star genome file or generate new genome index. <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.star_genome_generate_params">
        <b>ppmSeqPreprocess.TrimAlignSort.star_genome_generate_params</b><br />
        <i>StarGenomeGenerateParams? </i> &mdash;
         Parameters for generating the star genome. Mandatory if aligner is star and not given star genome file. <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.star_align_extra_args">
        <b>ppmSeqPreprocess.TrimAlignSort.star_align_extra_args</b><br />
        <i>String? </i> &mdash;
         Extra arguments for the STAR aligner. <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.star_align_gtf_override">
        <b>ppmSeqPreprocess.TrimAlignSort.star_align_gtf_override</b><br />
        <i>File? </i> &mdash;
         GTF file to be used for STAR aligner. <br />
</p>
<p name="ppmSeqPreprocess.TrimAlignSort.create_md5_checksum_outputs">
        <b>ppmSeqPreprocess.TrimAlignSort.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash;
         Create md5 checksum for requested output files <br />
</p>
</details>


## Outputs
<p name="ppmSeqPreprocess.output_cram_bam">
        <b>ppmSeqPreprocess.output_cram_bam</b><br />
        <i>File</i><br />
        Output CRAM or BAM file, trimmed aligned and sorted
</p>
<p name="ppmSeqPreprocess.output_cram_bam_index">
        <b>ppmSeqPreprocess.output_cram_bam_index</b><br />
        <i>File</i><br />
        Output CRAM or BAM index file
</p>
<p name="ppmSeqPreprocess.sorter_stats_csv">
        <b>ppmSeqPreprocess.sorter_stats_csv</b><br />
        <i>File</i><br />
        Sorter stats csv output
</p>
<p name="ppmSeqPreprocess.sorter_stats_json">
        <b>ppmSeqPreprocess.sorter_stats_json</b><br />
        <i>File</i><br />
        Sorter stats json output
</p>
<p name="ppmSeqPreprocess.unmatched_cram">
        <b>ppmSeqPreprocess.unmatched_cram</b><br />
        <i>File?</i><br />
        Unmatched cram file output from sorter (if defined).
</p>
<p name="ppmSeqPreprocess.unmatched_sorter_stats_csv">
        <b>ppmSeqPreprocess.unmatched_sorter_stats_csv</b><br />
        <i>File?</i><br />
        Unmatched output cram files sorter stats in csv format (if defined).
</p>
<p name="ppmSeqPreprocess.unmatched_sorter_stats_json">
        <b>ppmSeqPreprocess.unmatched_sorter_stats_json</b><br />
        <i>File?</i><br />
        Unmatched output cram files sorter stats in json format (if defined).
</p>
<p name="ppmSeqPreprocess.bedgraph_mapq0">
        <b>ppmSeqPreprocess.bedgraph_mapq0</b><br />
        <i>File?</i><br />
        Bedgraph mapq0 file.
</p>
<p name="ppmSeqPreprocess.bedgraph_mapq1">
        <b>ppmSeqPreprocess.bedgraph_mapq1</b><br />
        <i>File?</i><br />
        Bedgraph mapq1 file.
</p>
<p name="ppmSeqPreprocess.report_html">
        <b>ppmSeqPreprocess.report_html</b><br />
        <i>File</i><br />
        ppmSeq QC report html
</p>
<p name="ppmSeqPreprocess.application_qc_h5">
        <b>ppmSeqPreprocess.application_qc_h5</b><br />
        <i>File</i><br />
        ppmSeq QC aggregated metrics h5
</p>
<p name="ppmSeqPreprocess.aggregated_metrics_json">
        <b>ppmSeqPreprocess.aggregated_metrics_json</b><br />
        <i>File</i><br />
        ppmSeq QC aggregated metrics json
</p>
<p name="ppmSeqPreprocess.md5_checksums_json">
        <b>ppmSeqPreprocess.md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files
</p>

<hr />

> Generated using WDL AID (1.0.1)