# SingleReadSNV
Single Read SNV Quality Recalibration workflow (single_read_snv wdl) assigns accurate quality scores to all SNV candidates. The output is a FeatureMap VCF file with the recalibrated SNV quality scores. The input cram file coverage must be over some minimal coverage for a new model to be trained and quality scores to be generated, otherwise a pre-trained model can be provided, or a FeatureMap with no scores is generated.

## Inputs

### Required inputs
<p name="SingleReadSNV.input_cram_bam">
        <b>SingleReadSNV.input_cram_bam</b><br />
        <i>File </i> &mdash; 
         Input CRAM or BAM file <br /> 
</p>
<p name="SingleReadSNV.input_cram_bam_index">
        <b>SingleReadSNV.input_cram_bam_index</b><br />
        <i>File </i> &mdash; 
         Input CRAM or BAM index file <br /> 
</p>
<p name="SingleReadSNV.sorter_json_stats_file">
        <b>SingleReadSNV.sorter_json_stats_file</b><br />
        <i>File </i> &mdash; 
         Sorter json stats file provided by the Ultima Genomics pipeline (same base name as the input cram/bam file with a json extension) <br /> 
</p>
<p name="SingleReadSNV.base_file_name">
        <b>SingleReadSNV.base_file_name</b><br />
        <i>String </i> &mdash; 
         Base file name for output files. The output files will be named [base_file_name].with_ml_qual.vcf.gz <br /> 
</p>

### Required parameters
<p name="SingleReadSNV.featuremap_params">
        <b>SingleReadSNV.featuremap_params</b><br />
        <i>FeatureMapParams </i> &mdash; 
         FeatureMap parameters, recommended value set in the template. <br /> 
</p>
<p name="SingleReadSNV.features">
        <b>SingleReadSNV.features</b><br />
        <i>Array[String] </i> &mdash; 
         Features to be used for training the SNV quality model, should match pre trained model if one is given, recommended value set in the template. <br /> 
</p>
<p name="SingleReadSNV.training_regions_interval_list">
        <b>SingleReadSNV.training_regions_interval_list</b><br />
        <i>File </i> &mdash; 
         Genomic regions to include in the training set, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.xgboost_params_file">
        <b>SingleReadSNV.xgboost_params_file</b><br />
        <i>File </i> &mdash; 
         XGBoost parameters file for training the SNV quality model, recommended value set in the template. <br /> 
</p>
<p name="SingleReadSNV.min_coverage_to_train_model">
        <b>SingleReadSNV.min_coverage_to_train_model</b><br />
        <i>Float </i> &mdash; 
         Minimum coverage to train the ML model, needed as label assignment (true/false SNV) is unreliable at low coverages, the recommended value is set in the template <br /> 
</p>

### Required references
<p name="SingleReadSNV.references">
        <b>SingleReadSNV.references</b><br />
        <i>References </i> &mdash; 
         Reference files: fasta, dict and fai, recommended value set in the template <br /> 
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="SingleReadSNV.training_regions_interval_list_index">
        <b>SingleReadSNV.training_regions_interval_list_index</b><br />
        <i>File &mdash; Default: None</i><br />
        Index for genomic regions to exclude from the training set, the recommended value is set in the template
</p>
<p name="SingleReadSNV.pre_trained_model_files">
        <b>SingleReadSNV.pre_trained_model_files</b><br />
        <i>Array[File]? &mdash; Default: None</i><br />
        Pre-trained ML model json files, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features
</p>
<p name="SingleReadSNV.pre_trained_srsnv_metadata_json">
        <b>SingleReadSNV.pre_trained_srsnv_metadata_json</b><br />
        <i>File? &mdash; Default: None</i><br />
        Pre-trained SNV quality model metadata json file, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features
</p>
<p name="SingleReadSNV.raise_exceptions_in_report">
        <b>SingleReadSNV.raise_exceptions_in_report</b><br />
        <i>Boolean &mdash; Default: None</i><br />
        Raise and exception and fail the pipeline if an error is raised in the QC report
</p>

### Optional inputs
<p name="SingleReadSNV.create_md5_checksum_outputs">
        <b>SingleReadSNV.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash; 
         Create md5 checksum for requested output files <br /> 
</p>

### Optional parameters
<p name="SingleReadSNV.single_read_snv_params">
        <b>SingleReadSNV.single_read_snv_params</b><br />
        <i>SingleReadSNVParams? </i> &mdash; 
         SingleReadSNV parameters, recommended value set in the template. <br /> 
</p>
</details>


## Outputs
<p name="SingleReadSNV.featuremap">
        <b>SingleReadSNV.featuremap</b><br />
        <i>File</i><br />
        FeatureMap VCF file with SNVQ values
</p>
<p name="SingleReadSNV.featuremap_index">
        <b>SingleReadSNV.featuremap_index</b><br />
        <i>File</i><br />
        FeatureMap VCF index file
</p>
<p name="SingleReadSNV.featuremap_random_sample">
        <b>SingleReadSNV.featuremap_random_sample</b><br />
        <i>File?</i><br />
        Downsampled FeatureMap VCF file for training
</p>
<p name="SingleReadSNV.featuremap_random_sample_index">
        <b>SingleReadSNV.featuremap_random_sample_index</b><br />
        <i>File?</i><br />
        Downsampled FeatureMap VCF index file
</p>
<p name="SingleReadSNV.downsampling_rate">
        <b>SingleReadSNV.downsampling_rate</b><br />
        <i>Float</i><br />
        The downsampling rate used to create the random sample featuremap
</p>
<p name="SingleReadSNV.snv_qualities_assigned">
        <b>SingleReadSNV.snv_qualities_assigned</b><br />
        <i>Boolean</i><br />
        Indicates whether SNV qualities were assigned (coverage was sufficient to train the model, over min_coverage_to_train_model), otherwise a FeatureMap with no quality scores is generated
</p>
<p name="SingleReadSNV.used_self_trained_model">
        <b>SingleReadSNV.used_self_trained_model</b><br />
        <i>Boolean</i><br />
        Indicates whether a self-trained model was used for inference, otherwise a pre-trained model was used (if snv_qualities_assigned) or no model was used
</p>
<p name="SingleReadSNV.md5_checksums_json">
        <b>SingleReadSNV.md5_checksums_json</b><br />
        <i>File?</i><br />
        json file that will contain md5 checksums for requested output files
</p>
<p name="SingleReadSNV.raw_filtered_featuremap_parquet">
        <b>SingleReadSNV.raw_filtered_featuremap_parquet</b><br />
        <i>File?</i><br />
        Filtered parquet dataframe of raw featuremap for training
</p>
<p name="SingleReadSNV.raw_featuremap_stats">
        <b>SingleReadSNV.raw_featuremap_stats</b><br />
        <i>File?</i><br />
        Statistics for raw featuremap filtering
</p>
<p name="SingleReadSNV.random_sample_filtered_featuremap_parquet">
        <b>SingleReadSNV.random_sample_filtered_featuremap_parquet</b><br />
        <i>File?</i><br />
        Filtered parquet dataframe of random sample featuremap for training
</p>
<p name="SingleReadSNV.random_sample_featuremap_stats">
        <b>SingleReadSNV.random_sample_featuremap_stats</b><br />
        <i>File?</i><br />
        Statistics for random sample featuremap filtering
</p>
<p name="SingleReadSNV.featuremap_df">
        <b>SingleReadSNV.featuremap_df</b><br />
        <i>File?</i><br />
        Parquet file with model data, including features and labels
</p>
<p name="SingleReadSNV.application_qc_h5">
        <b>SingleReadSNV.application_qc_h5</b><br />
        <i>File?</i><br />
        Application QC statistics h5 file
</p>
<p name="SingleReadSNV.report_html">
        <b>SingleReadSNV.report_html</b><br />
        <i>File?</i><br />
        SRSNV QC report html file
</p>
<p name="SingleReadSNV.srsnv_metadata_json">
        <b>SingleReadSNV.srsnv_metadata_json</b><br />
        <i>File?</i><br />
        Metadata JSON file for the SNV quality model, containing information about the model, features, and training parameters
</p>
<p name="SingleReadSNV.model_files">
        <b>SingleReadSNV.model_files</b><br />
        <i>Array[File]?</i><br />
        Array of model files, each file should be used for a specific set of chromosome
</p>

<hr />

> Generated using WDL AID (1.0.1)
