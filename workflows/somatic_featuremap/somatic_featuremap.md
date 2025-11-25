# SomaticFeaturemap
The SomaticFeaturemap workflow generates a merged somatic featuremap VCF file from tumor and normal featuremap VCF files. This workflow is designed to combine information from tumor-normal pairs and enrich the data with additional features for somatic variant analysis.

## Inputs

### Required inputs
<p name="SomaticFeaturemap.base_file_name">
        <b>SomaticFeaturemap.base_file_name</b><br />
        <i>String </i> &mdash; 
         Sample name <br /> 
</p>
<p name="SomaticFeaturemap.tumor_featuremap_vcf">
        <b>SomaticFeaturemap.tumor_featuremap_vcf</b><br />
        <i>File? </i> &mdash; 
         Input tumor featuremap VCF file. <br /> 
</p>
<p name="SomaticFeaturemap.tumor_featuremap_vcf_index">
        <b>SomaticFeaturemap.tumor_featuremap_vcf_index</b><br />
        <i>File? </i> &mdash; 
         Input tumor featuremap VCF index file. <br /> 
</p>
<p name="SomaticFeaturemap.normal_featuremap_vcf">
        <b>SomaticFeaturemap.normal_featuremap_vcf</b><br />
        <i>File? </i> &mdash; 
         Input normal featuremap VCF file. <br /> 
</p>
<p name="SomaticFeaturemap.normal_featuremap_vcf_index">
        <b>SomaticFeaturemap.normal_featuremap_vcf_index</b><br />
        <i>File? </i> &mdash; 
         Input normal featuremap VCF index file. <br /> 
</p>
<p name="SomaticFeaturemap.ref_tandem_repeats">
        <b>SomaticFeaturemap.ref_tandem_repeats</b><br />
        <i>File </i> &mdash; 
         Reference tandem repeats file. <br /> 
</p>
<p name="SomaticFeaturemap.input_tumor_cram_bam">
        <b>SomaticFeaturemap.input_tumor_cram_bam</b><br />
        <i>Array[File]+? </i> &mdash; 
         Input tumor CRAM/BAM file. <br /> 
</p>
<p name="SomaticFeaturemap.input_tumor_cram_bam_index">
        <b>SomaticFeaturemap.input_tumor_cram_bam_index</b><br />
        <i>Array[File]+? </i> &mdash; 
         Input tumor CRAM/BAM index file. <br /> 
</p>
<p name="SomaticFeaturemap.input_normal_cram_bam">
        <b>SomaticFeaturemap.input_normal_cram_bam</b><br />
        <i>Array[File]+? </i> &mdash; 
         Input normal CRAM/BAM file. <br /> 
</p>
<p name="SomaticFeaturemap.input_normal_cram_bam_index">
        <b>SomaticFeaturemap.input_normal_cram_bam_index</b><br />
        <i>Array[File]+? </i> &mdash; 
         Input normal CRAM/BAM index file. <br /> 
</p>
<p name="SomaticFeaturemap.interval_list">
        <b>SomaticFeaturemap.interval_list</b><br />
        <i>File? </i> &mdash; 
         Input interval list file. <br /> 
</p>
<p name="SomaticFeaturemap.references">
        <b>SomaticFeaturemap.references</b><br />
        <i>References? </i> &mdash; 
         Input reference genome file, index file and dict. <br /> 
</p>

### Optional inputs
<p name="SomaticFeaturemap.somatic_featuremap_given">
        <b>SomaticFeaturemap.somatic_featuremap_given</b><br />
        <i>File? </i> &mdash; 
         If provided, use this somatic featuremap file instead of generating a new one. <br /> 
</p>
<p name="SomaticFeaturemap.somatic_featuremap_given_index">
        <b>SomaticFeaturemap.somatic_featuremap_given_index</b><br />
        <i>File? </i> &mdash; 
         If provided, use this somatic featuremap index file instead of generating a new one. <br /> 
</p>
<p name="SomaticFeaturemap.keep_non_pass_tumor_candidates_override">
        <b>SomaticFeaturemap.keep_non_pass_tumor_candidates_override</b><br />
        <i>Boolean? </i> &mdash; 
         If true, keep non-PASS tumor candidates in the somatic featuremap. <br /> 
</p>
<p name="SomaticFeaturemap.run_mpileup_override">
        <b>SomaticFeaturemap.run_mpileup_override</b><br />
        <i>Boolean? </i> &mdash; 
         If true, run the mpileup step. <br /> 
</p>
<p name="SomaticFeaturemap.sfm_xgb_model">
        <b>SomaticFeaturemap.sfm_xgb_model</b><br />
        <i>File? </i> &mdash; 
         XGBoost model file for somatic featuremap pileup scoring. <br /> 
</p>
<p name="SomaticFeaturemap.pad_size">
        <b>SomaticFeaturemap.pad_size</b><br />
        <i>Int </i> &mdash; 
         Size of padding to add to the positions of variants defined as tumor-PASS in somatic-pileup-featuremap. <br /> 
</p>
<p name="SomaticFeaturemap.min_mapq_for_samtools_mpileup">
        <b>SomaticFeaturemap.min_mapq_for_samtools_mpileup</b><br />
        <i>Int </i> &mdash; 
         Minimum mapping quality for samtools mpileup. <br /> 
</p>
<p name="SomaticFeaturemap.num_shards">
        <b>SomaticFeaturemap.num_shards</b><br />
        <i>Int? </i> &mdash; 
         Number of shards for splitting the interval list. <br /> 
</p>
<p name="SomaticFeaturemap.scatter_intervals_break">
        <b>SomaticFeaturemap.scatter_intervals_break</b><br />
        <i>Int? </i> &mdash; 
         Break scatter intervals at multiples of this value. <br /> 
</p>
<p name="SomaticFeaturemap.dummy_input_for_call_caching">
        <b>SomaticFeaturemap.dummy_input_for_call_caching</b><br />
        <i>String </i> &mdash; 
         Dummy input for call caching. <br /> 
</p>
</details>


## Outputs
<p name="SomaticFeaturemap.somatic_featuremap">
        <b>SomaticFeaturemap.somatic_featuremap</b><br />
        <i>File?</i><br />
        Output somatic featuremap pileup file.
</p>
<p name="SomaticFeaturemap.somatic_featuremap_index">
        <b>SomaticFeaturemap.somatic_featuremap_index</b><br />
        <i>File?</i><br />
        Output somatic featuremap pileup index file.
</p>
<p name="SomaticFeaturemap.tumor_mpileup">
        <b>SomaticFeaturemap.tumor_mpileup</b><br />
        <i>File?</i><br />
        Output tumor mpileup file.
</p>
<p name="SomaticFeaturemap.normal_mpileup">
        <b>SomaticFeaturemap.normal_mpileup</b><br />
        <i>File?</i><br />
        Output normal mpileup file.
</p>
<p name="SomaticFeaturemap.somatic_featuremap_mpileup_vcf">
        <b>SomaticFeaturemap.somatic_featuremap_mpileup_vcf</b><br />
        <i>File?</i><br />
        Output somatic featuremap pileup VCF with integrated mpileup info.
</p>
<p name="SomaticFeaturemap.somatic_featuremap_mpileup_vcf_index">
        <b>SomaticFeaturemap.somatic_featuremap_mpileup_vcf_index</b><br />
        <i>File?</i><br />
        somatic_featuremap_mpileup_vcf index file.
</p>
<p name="SomaticFeaturemap.final_out_sfm_vcf">
        <b>SomaticFeaturemap.final_out_sfm_vcf</b><br />
        <i>File</i><br />
        Final output somatic featuremap pileup VCF file.
</p>
<p name="SomaticFeaturemap.final_out_sfm_vcf_index">
        <b>SomaticFeaturemap.final_out_sfm_vcf_index</b><br />
        <i>File</i><br />
        Final output somatic featuremap pileup VCF index file.
</p>
<p name="SomaticFeaturemap.somatic_featuremap_padded_bed_file">
        <b>SomaticFeaturemap.somatic_featuremap_padded_bed_file</b><br />
        <i>File?</i><br />
        Output padded BED file for somatic featuremap variants.
</p>

<hr />

> Generated using WDL AID (1.0.1)
