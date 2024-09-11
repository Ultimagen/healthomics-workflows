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
<p name="SingleReadSNV.wgs_calling_interval_list">
        <b>SingleReadSNV.wgs_calling_interval_list</b><br />
        <i>File </i> &mdash; 
         interval list defining the region to perform variant calling on, recommended value set in the template <br /> 
</p>
<p name="SingleReadSNV.break_bands_at_multiples_of">
        <b>SingleReadSNV.break_bands_at_multiples_of</b><br />
        <i>Int </i> &mdash; 
         Break wgs_calling_interval_list bands at multiples of this number, recommended value set in the template <br /> 
</p>
<p name="SingleReadSNV.featuremap_params">
        <b>SingleReadSNV.featuremap_params</b><br />
        <i>FeatureMapParams </i> &mdash; 
         FeatureMap parameters, recommended value set in the template. <br /> 
</p>
<p name="SingleReadSNV.single_read_snv_params">
        <b>SingleReadSNV.single_read_snv_params</b><br />
        <i>SingleReadSNVParams </i> &mdash; 
         SingleReadSNV parameters, recommended value set in the template. <br /> 
</p>
<p name="SingleReadSNV.categorical_features">
        <b>SingleReadSNV.categorical_features</b><br />
        <i>Map[String,Array[String]] </i> &mdash; 
         Categorical features in SingleReadSNV model, list of feature names with allowed values per feature. Separate from single_read_snv_params due to technical reasons. The recommended value is set in the template. <br /> 
</p>
<p name="SingleReadSNV.training_include_regions">
        <b>SingleReadSNV.training_include_regions</b><br />
        <i>Array[File] </i> &mdash; 
         Genomic regions to include in the training set, the recommended value is set in the template <br /> 
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
<p name="SingleReadSNV.somatic_mutations_list">
        <b>SingleReadSNV.somatic_mutations_list</b><br />
        <i>Array[File]? &mdash; Default: None</i><br />
        Somatic mutations to be excluded from FP training set, will be appended to the fp_training_exclude_regions optional
</p>
<p name="SingleReadSNV.tp_training_exclude_regions">
        <b>SingleReadSNV.tp_training_exclude_regions</b><br />
        <i>Array[File]? &mdash; Default: None</i><br />
        Genomic regions to exclude from the training set TP examples, the recommended value is set in the template
</p>
<p name="SingleReadSNV.fp_training_exclude_regions">
        <b>SingleReadSNV.fp_training_exclude_regions</b><br />
        <i>Array[File]? &mdash; Default: None</i><br />
        Genomic regions to exclude from the training set FP examples, the recommended value is set in the template
</p>
<p name="SingleReadSNV.pre_trained_model_file">
        <b>SingleReadSNV.pre_trained_model_file</b><br />
        <i>File? &mdash; Default: None</i><br />
        Pre-trained ML model file, if provided the model will be used for inference and no self-trained model will be created. Use with care, the model must be trained on the same data type with the same features
</p>
<p name="SingleReadSNV.raise_exceptions_in_report">
        <b>SingleReadSNV.raise_exceptions_in_report</b><br />
        <i>Boolean &mdash; Default: None</i><br />
        Raise and exception and fail the pipeline if an error is raised in the QC report
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
<p name="SingleReadSNV.report_html">
        <b>SingleReadSNV.report_html</b><br />
        <i>File?</i><br />
        SRSNV QC report html file
</p>
<p name="SingleReadSNV.featuremap_df_file">
        <b>SingleReadSNV.featuremap_df_file</b><br />
        <i>File?</i><br />
        FeatureMap DataFrame with X matrix (features), y (labels) and qual, in parquet format
</p>
<p name="SingleReadSNV.model_file">
        <b>SingleReadSNV.model_file</b><br />
        <i>File?</i><br />
        ML model file, saved with joblib
</p>
<p name="SingleReadSNV.params_file">
        <b>SingleReadSNV.params_file</b><br />
        <i>File?</i><br />
        ML model parameters json file
</p>
<p name="SingleReadSNV.test_set_mrd_simulation_dataframe">
        <b>SingleReadSNV.test_set_mrd_simulation_dataframe</b><br />
        <i>File?</i><br />
        ML model test set MRD simulation DataFrame, parquet format
</p>
<p name="SingleReadSNV.test_set_statistics_h5">
        <b>SingleReadSNV.test_set_statistics_h5</b><br />
        <i>File?</i><br />
        ML model test set statistics h5 file
</p>
<p name="SingleReadSNV.aggregated_metrics_json">
        <b>SingleReadSNV.aggregated_metrics_json</b><br />
        <i>File?</i><br />
        ML model aggregated metrics json file
</p>
<p name="SingleReadSNV.tp_training_regions_bed">
        <b>SingleReadSNV.tp_training_regions_bed</b><br />
        <i>File?</i><br />
        ML model training set TP regions bed file
</p>
<p name="SingleReadSNV.fp_training_regions_bed">
        <b>SingleReadSNV.fp_training_regions_bed</b><br />
        <i>File?</i><br />
        ML model training set FP regions bed file
</p>
<p name="SingleReadSNV.test_report_file_notebook">
        <b>SingleReadSNV.test_report_file_notebook</b><br />
        <i>File?</i><br />
        ML model test set report notebook file
</p>

<hr />

> Generated using WDL AID (1.0.0)
