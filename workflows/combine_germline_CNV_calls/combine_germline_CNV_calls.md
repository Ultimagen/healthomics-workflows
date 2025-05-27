# CombineGermlineCNVCalls
Combine cn.mops and cnvpytor CNV calls for a single sample.

## Inputs

### Required inputs
<p name="CombineGermlineCNVCalls.base_file_name">
        <b>CombineGermlineCNVCalls.base_file_name</b><br />
        <i>String </i> &mdash; 
         Sample name <br /> 
</p>
<p name="CombineGermlineCNVCalls.cnmops_cnvs_bed">
        <b>CombineGermlineCNVCalls.cnmops_cnvs_bed</b><br />
        <i>File </i> &mdash; 
         cn.mops CNV calls in bed format <br /> 
</p>
<p name="CombineGermlineCNVCalls.cnvpytor_cnvs_bed">
        <b>CombineGermlineCNVCalls.cnvpytor_cnvs_bed</b><br />
        <i>File </i> &mdash; 
         cnvpytor CNV calls in bed format <br /> 
</p>
<p name="CombineGermlineCNVCalls.input_bam_file">
        <b>CombineGermlineCNVCalls.input_bam_file</b><br />
        <i>File </i> &mdash; 
         Input sample BAM/CRAM file. <br /> 
</p>
<p name="CombineGermlineCNVCalls.input_bam_file_index">
        <b>CombineGermlineCNVCalls.input_bam_file_index</b><br />
        <i>File </i> &mdash; 
         Input sample BAI/CRAI index file <br /> 
</p>

### Required references
<p name="CombineGermlineCNVCalls.reference_genome">
        <b>CombineGermlineCNVCalls.reference_genome</b><br />
        <i>File </i> &mdash; 
         Genome fasta file associated with the CRAM file <br /> 
</p>
<p name="CombineGermlineCNVCalls.reference_genome_index">
        <b>CombineGermlineCNVCalls.reference_genome_index</b><br />
        <i>File </i> &mdash; 
         Fai index of the fasta file <br /> 
</p>

### Optional inputs
<p name="CombineGermlineCNVCalls.cnv_lcr_file">
        <b>CombineGermlineCNVCalls.cnv_lcr_file</b><br />
        <i>File? </i> &mdash; 
         UG-CNV-LCR bed file <br /> 
</p>

### Optional parameters
<p name="CombineGermlineCNVCalls.preemptible_tries_override">
        <b>CombineGermlineCNVCalls.preemptible_tries_override</b><br />
        <i>Int? </i> &mdash; 
         Number of preemptible tries,default is: 1 <br /> 
</p>
</details>


## Outputs
<p name="CombineGermlineCNVCalls.out_jalign_del_bed">
        <b>CombineGermlineCNVCalls.out_jalign_del_bed</b><br />
        <i>File</i><br />
        Output file with DEL candidates after jalign
</p>
<p name="CombineGermlineCNVCalls.out_sample_cnvs_bed">
        <b>CombineGermlineCNVCalls.out_sample_cnvs_bed</b><br />
        <i>File</i><br />
        Bed file with sample's called CNVs
</p>
<p name="CombineGermlineCNVCalls.out_sample_cnvs_vcf">
        <b>CombineGermlineCNVCalls.out_sample_cnvs_vcf</b><br />
        <i>File</i><br />
        VCF file with sample's called CNVs
</p>
<p name="CombineGermlineCNVCalls.out_sample_cnvs_vcf_index">
        <b>CombineGermlineCNVCalls.out_sample_cnvs_vcf_index</b><br />
        <i>File</i><br />
        Index file for the VCF file with sample's called CNVs
</p>

<hr />

> Generated using WDL AID (1.0.1)
