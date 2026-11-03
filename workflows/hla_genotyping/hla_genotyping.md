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
<p name="HLAGenotyping.t1k_index_tar">
        <b>HLAGenotyping.t1k_index_tar</b><br />
        <i>File? </i> &mdash;
         T1K index tar.gz containing hlaidx/ and kiridx/ directories with all index files (_seq.fa and _coord.fa) <br />
</p>
<p name="HLAGenotyping.reference_genome">
        <b>HLAGenotyping.reference_genome</b><br />
        <i>String </i> &mdash;
         Genome type selector (hg38 or hg38_nist_v3_with_decoy). Determines which reference files to use. <br />
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