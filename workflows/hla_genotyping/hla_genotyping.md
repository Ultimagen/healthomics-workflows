# HLAGenotyping
HLA Genotyping

## Inputs

### Required inputs
<p name="HLAGenotyping.base_file_name">
        <b>HLAGenotyping.base_file_name</b><br />
        <i>String </i> &mdash;
         Base file name for the output files (to be used as the prefix) <br />
</p>
<p name="HLAGenotyping.input_cram_bam">
        <b>HLAGenotyping.input_cram_bam</b><br />
        <i>File </i> &mdash;
         Input CRAM or BAM file for annalysing HLA genotyping <br />
</p>
<p name="HLAGenotyping.input_cram_bam_index">
        <b>HLAGenotyping.input_cram_bam_index</b><br />
        <i>File </i> &mdash;
         Input CRAM or BAM index file for annalysing HLA genotyping <br />
</p>
<p name="HLAGenotyping.hla_genotyping_tool">
        <b>HLAGenotyping.hla_genotyping_tool</b><br />
        <i>String </i> &mdash;
         HLA genotyping tool to use. Options: 'HLA-LA' or 'T1K' <br />
</p>
<p name="HLAGenotyping.reference_genome">
        <b>HLAGenotyping.reference_genome</b><br />
        <i>String </i> &mdash;
         Genome type selector (hg38 or hg38_nist_v3_with_decoy). Determines which reference files to use. <br />
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="HLAGenotyping.graphs_files_tar">
        <b>HLAGenotyping.graphs_files_tar</b><br />
        <i>File? &mdash; Default: None</i><br />
        HLA-LA graphs files tar (required if using HLA-LA)
</p>
<p name="HLAGenotyping.t1k_index_tar">
        <b>HLAGenotyping.t1k_index_tar</b><br />
        <i>File? &mdash; Default: None</i><br />
        T1K index tar.gz containing hlaidx/ and kiridx/ directories with all index files (_seq.fa and _coord.fa). Required if using T1K.
</p>
</details>


## Outputs
<p name="HLAGenotyping.output_hla">
        <b>HLAGenotyping.output_hla</b><br />
        <i>File</i><br />
        HLA genotyping output file
</p>
<p name="HLAGenotyping.output_kir">
        <b>HLAGenotyping.output_kir</b><br />
        <i>File?</i><br />
        KIR genotyping output file
</p>

<hr />

> Generated using WDL AID (1.0.1)