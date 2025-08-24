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
<p name="HLAGenotyping.references">
        <b>HLAGenotyping.references</b><br />
        <i>References </i> &mdash; 
         Reference files: fasta, dict and fai, recommended value set in the template <br /> 
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="HLAGenotyping.graphs_files_tar">
        <b>HLAGenotyping.graphs_files_tar</b><br />
        <i>File &mdash; Default: None</i><br />
        HLA-LA graphs files tar
</p>
</details>


## Outputs
<p name="HLAGenotyping.output_hla">
        <b>HLAGenotyping.output_hla</b><br />
        <i>File</i><br />
        HLA genotyping output file
</p>

<hr />

> Generated using WDL AID (1.0.1)
