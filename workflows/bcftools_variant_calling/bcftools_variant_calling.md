# BcftoolsVariantCalling
BcftoolsVariantCalling: bcftools-based variant calling for CRAM files

## Inputs

### Required inputs
<p name="BcftoolsVariantCalling.base_file_name">
        <b>BcftoolsVariantCalling.base_file_name</b><br />
        <i>String </i> &mdash;
         Base name for all output files <br />
</p>
<p name="BcftoolsVariantCalling.input_cram">
        <b>BcftoolsVariantCalling.input_cram</b><br />
        <i>File </i> &mdash;
         Input CRAM file with aligned reads <br />
</p>
<p name="BcftoolsVariantCalling.input_cram_index">
        <b>BcftoolsVariantCalling.input_cram_index</b><br />
        <i>File </i> &mdash;
         CRAM index file (.crai) <br />
</p>
<p name="BcftoolsVariantCalling.autosomes_bed">
        <b>BcftoolsVariantCalling.autosomes_bed</b><br />
        <i>File </i> &mdash;
         BED file with autosomal targets <br />
</p>
<p name="BcftoolsVariantCalling.chrx_par_bed">
        <b>BcftoolsVariantCalling.chrx_par_bed</b><br />
        <i>File </i> &mdash;
         BED file with X chromosome PAR (pseudoautosomal) regions <br />
</p>
<p name="BcftoolsVariantCalling.chrx_nonpar_bed">
        <b>BcftoolsVariantCalling.chrx_nonpar_bed</b><br />
        <i>File </i> &mdash;
         BED file with X chromosome non-PAR regions <br />
</p>
<p name="BcftoolsVariantCalling.chry_nonpar_bed">
        <b>BcftoolsVariantCalling.chry_nonpar_bed</b><br />
        <i>File </i> &mdash;
         BED file with Y chromosome non-PAR regions (used for XY samples) <br />
</p>
<p name="BcftoolsVariantCalling.reference_genome">
        <b>BcftoolsVariantCalling.reference_genome</b><br />
        <i>String </i> &mdash;
         Genome type for reference selection (hg38, b37). Reference files are automatically selected based on this value. <br />
</p>

### Optional inputs
<p name="BcftoolsVariantCalling.sex">
        <b>BcftoolsVariantCalling.sex</b><br />
        <i>String? </i> &mdash;
         Sample sex: XX (female) or XY (male). If not provided and auto_determine_sex=true, will be determined automatically from CRAM <br />
</p>
<p name="BcftoolsVariantCalling.auto_determine_sex">
        <b>BcftoolsVariantCalling.auto_determine_sex</b><br />
        <i>Boolean </i> &mdash;
         Automatically determine sex from CRAM file if sex parameter not provided (default: true) <br />
</p>
<p name="BcftoolsVariantCalling.run_filtering">
        <b>BcftoolsVariantCalling.run_filtering</b><br />
        <i>Boolean </i> &mdash;
         Whether to generate filtered gVCF using depth and GQ thresholds <br />
</p>
<p name="BcftoolsVariantCalling.num_autosome_shards">
        <b>BcftoolsVariantCalling.num_autosome_shards</b><br />
        <i>Int </i> &mdash;
         Number of shards to split autosomes for parallel processing <br />
</p>
<p name="BcftoolsVariantCalling.memory_gb">
        <b>BcftoolsVariantCalling.memory_gb</b><br />
        <i>Int </i> &mdash;
         Memory in GB for variant calling tasks <br />
</p>
<p name="BcftoolsVariantCalling.cpus">
        <b>BcftoolsVariantCalling.cpus</b><br />
        <i>Int </i> &mdash;
         Number of CPUs for variant calling tasks <br />
</p>
</details>


## Outputs
<p name="BcftoolsVariantCalling.final_gvcf">
        <b>BcftoolsVariantCalling.final_gvcf</b><br />
        <i>File</i><br />
        Final merged gVCF file
</p>
<p name="BcftoolsVariantCalling.final_gvcf_index">
        <b>BcftoolsVariantCalling.final_gvcf_index</b><br />
        <i>File</i><br />
        Index for final gVCF
</p>
<p name="BcftoolsVariantCalling.filtered_gvcf">
        <b>BcftoolsVariantCalling.filtered_gvcf</b><br />
        <i>File?</i><br />
        Filtered gVCF with depth and GQ filtering applied
</p>
<p name="BcftoolsVariantCalling.filtered_gvcf_index">
        <b>BcftoolsVariantCalling.filtered_gvcf_index</b><br />
        <i>File?</i><br />
        Index for filtered gVCF
</p>
<p name="BcftoolsVariantCalling.determined_sex">
        <b>BcftoolsVariantCalling.determined_sex</b><br />
        <i>String</i><br />
        Final sex used for variant calling (provided or auto-determined): XX or XY
</p>
<p name="BcftoolsVariantCalling.sex_metrics">
        <b>BcftoolsVariantCalling.sex_metrics</b><br />
        <i>File?</i><br />
        Sex determination metrics file (only present if sex was auto-determined)
</p>

<hr />

> Generated using WDL AID (1.0.1)