# SomaticSNVfind
The SomaticSNVfind workflow performs somatic variant calling, producing an annotated somatic FeatureMap VCF file containing information on a tumor and a normal sample. It supports two execution modes: (1) Full pipeline mode, which takes matched tumor and normal CRAM files, divides the genome into regions, runs snvfind jointly on the tumor–normal pair per region using model metadata, applies somatic classification per region, and merges the results into a final VCF. (2) Classifier-only mode, which takes a pre-existing somatic FeatureMap VCF and applies the classification step only.

## Inputs

### Required inputs
<p name="SomaticSNVfind.base_file_name">
        <b>SomaticSNVfind.base_file_name</b><br />
        <i>String </i> &mdash;
         Sample name <br />
</p>
<p name="SomaticSNVfind.references">
        <b>SomaticSNVfind.references</b><br />
        <i>References </i> &mdash;
         Input reference genome file, index file and dict. <br />
</p>
<p name="SomaticSNVfind.interval_list">
        <b>SomaticSNVfind.interval_list</b><br />
        <i>File </i> &mdash;
         Input interval list file for scattering. <br />
</p>
<p name="SomaticSNVfind.num_shards">
        <b>SomaticSNVfind.num_shards</b><br />
        <i>Int </i> &mdash;
         Number of shards for splitting the interval list. <br />
</p>
<p name="SomaticSNVfind.tandem_repeats_bed">
        <b>SomaticSNVfind.tandem_repeats_bed</b><br />
        <i>File </i> &mdash;
         Reference tandem repeats file. <br />
</p>
<p name="SomaticSNVfind.xgb_model">
        <b>SomaticSNVfind.xgb_model</b><br />
        <i>File </i> &mdash;
         XGBoost model file for somatic featuremap classification. <br />
</p>

### Optional inputs
<p name="SomaticSNVfind.xgb_proba_threshold">
        <b>SomaticSNVfind.xgb_proba_threshold</b><br />
        <i>Float </i> &mdash;
         XGBoost probability threshold for somatic featuremap classification. <br />
</p>
<p name="SomaticSNVfind.scatter_intervals_break">
        <b>SomaticSNVfind.scatter_intervals_break</b><br />
        <i>Int? </i> &mdash;
         Break scatter intervals at multiples of this value. <br />
</p>
<p name="SomaticSNVfind.tumor_crams">
        <b>SomaticSNVfind.tumor_crams</b><br />
        <i>Array[File]? </i> &mdash;
         Array of input tumor CRAM files. All files must have the same sample name. Required for full pipeline mode, not needed for classifier-only mode. <br />
</p>
<p name="SomaticSNVfind.tumor_cram_index_list">
        <b>SomaticSNVfind.tumor_cram_index_list</b><br />
        <i>Array[File]? </i> &mdash;
         Array of input tumor CRAM index files, matching the order of tumor_crams. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.normal_crams">
        <b>SomaticSNVfind.normal_crams</b><br />
        <i>Array[File]? </i> &mdash;
         Array of input normal CRAM files. All files must have the same sample name. Required for full pipeline mode, not needed for classifier-only mode. <br />
</p>
<p name="SomaticSNVfind.normal_cram_index_list">
        <b>SomaticSNVfind.normal_cram_index_list</b><br />
        <i>Array[File]? </i> &mdash;
         Array of input normal CRAM index files, matching the order of normal_crams. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.tumor_single_read_snv_model">
        <b>SomaticSNVfind.tumor_single_read_snv_model</b><br />
        <i>SingleReadSNVModel? </i> &mdash;
         Tumor SingleReadSNVModel struct containing model_metadata and model_fold_files. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.normal_single_read_snv_model">
        <b>SomaticSNVfind.normal_single_read_snv_model</b><br />
        <i>SingleReadSNVModel? </i> &mdash;
         Normal SingleReadSNVModel struct containing model_metadata and model_fold_files. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.featuremap_params">
        <b>SomaticSNVfind.featuremap_params</b><br />
        <i>FeatureMapParams? </i> &mdash;
         FeatureMap parameters for snvfind. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.override_memory_gb_CreateFeatureMap">
        <b>SomaticSNVfind.override_memory_gb_CreateFeatureMap</b><br />
        <i>Int? </i> &mdash;
         Override memory in GB for the CreateFeatureMap task, default: 2 (GiB). If an out of memory error occurs in the CreateFeatureMap task, try increasing this value, e.g. double it. <br />
</p>
<p name="SomaticSNVfind.annotation_files">
        <b>SomaticSNVfind.annotation_files</b><br />
        <i>FeaturemapAnnotationFiles? </i> &mdash;
         Annotation files for featuremap generation: dbSNP, gnomAD, and UG High Confidence Regions with their indices. Required for full pipeline mode. <br />
</p>
<p name="SomaticSNVfind.input_somatic_snvfind_vcf">
        <b>SomaticSNVfind.input_somatic_snvfind_vcf</b><br />
        <i>File? </i> &mdash;
         Pre-existing somatic featuremap VCF file. When provided, runs classifier-only mode (skips featuremap generation). <br />
</p>
<p name="SomaticSNVfind.input_somatic_snvfind_vcf_index">
        <b>SomaticSNVfind.input_somatic_snvfind_vcf_index</b><br />
        <i>File? </i> &mdash;
         Index file for pre-existing somatic featuremap VCF. Required when input_somatic_snvfind_vcf is provided. <br />
</p>
</details>


## Outputs
<p name="SomaticSNVfind.somatic_snvfind_vcf">
        <b>SomaticSNVfind.somatic_snvfind_vcf</b><br />
        <i>File</i><br />
        Output final somatic snvfind VCF file.
</p>
<p name="SomaticSNVfind.somatic_snvfind_vcf_index">
        <b>SomaticSNVfind.somatic_snvfind_vcf_index</b><br />
        <i>File</i><br />
        Output final somatic snvfind VCF index file.
</p>
<p name="SomaticSNVfind.aggregated_parquet">
        <b>SomaticSNVfind.aggregated_parquet</b><br />
        <i>Array[File]</i><br />
        Array of aggregated parquet files, one per genomic region/shard, primarily intended for debugging and troubleshooting. Each parquet file contains structured data from the somatic featuremap classifier including features, classification scores, and results in a columnar format. These files facilitate debugging by providing detailed intermediate results that can be inspected, analyzed, and visualized more easily than VCF format. Useful for investigating classification behavior, feature distributions, and pipeline issues.
</p>

<hr />

> Generated using WDL AID (1.0.1)