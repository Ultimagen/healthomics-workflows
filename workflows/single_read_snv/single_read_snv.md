# SingleReadSNV
Single Read SNV Quality Recalibration workflow (single_read_snv wdl) is a software tool that assigns accurate quality scores to all SNV candidates. The output is a FeatureMap VCF file with the recalibrated SNV quality scores.

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
         FeatureMap parameters, recommended value set in the template. Int scatter_count: number of scatter tasks to use. Int min_mapq, Int snv_identical_bases, Int snv_identical_bases_after, Int min_score, Int limit_score, String extra_args <br /> 
</p>
<p name="SingleReadSNV.balanced_strand_adapter_version">
        <b>SingleReadSNV.balanced_strand_adapter_version</b><br />
        <i>String? </i> &mdash; 
         ppmSeq adapter version, for ppmSeq data the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.motif_length_to_annotate">
        <b>SingleReadSNV.motif_length_to_annotate</b><br />
        <i>Int </i> &mdash; 
         Length of the motif (-+N bp) to annotate in the FeatureMap, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.max_hmer_length">
        <b>SingleReadSNV.max_hmer_length</b><br />
        <i>Int </i> &mdash; 
         Maximum length of the homopolymer to annotate in the FeatureMap, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.train_set_size">
        <b>SingleReadSNV.train_set_size</b><br />
        <i>Int </i> &mdash; 
         Number of SNVs to use for the ML model training set, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.test_set_size">
        <b>SingleReadSNV.test_set_size</b><br />
        <i>Int </i> &mdash; 
         Number of SNVs to use for the ML model test set, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.tp_training_include_regions">
        <b>SingleReadSNV.tp_training_include_regions</b><br />
        <i>Array[File] </i> &mdash; 
         Genomic regions to include in the training set TP examples, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.fp_training_include_regions">
        <b>SingleReadSNV.fp_training_include_regions</b><br />
        <i>Array[File] </i> &mdash; 
         Genomic regions to include in the training set FP examples, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.numerical_features">
        <b>SingleReadSNV.numerical_features</b><br />
        <i>Array[String] </i> &mdash; 
         Numerical features to use in the ML model, the recommended value is set in the template <br /> 
</p>
<p name="SingleReadSNV.categorical_features">
        <b>SingleReadSNV.categorical_features</b><br />
        <i>Array[String] </i> &mdash; 
         Categorical features to use in the ML model, the recommended value is set in the template <br /> 
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
<p name="SingleReadSNV.balanced_sampling_info_fields">
        <b>SingleReadSNV.balanced_sampling_info_fields</b><br />
        <i>Array[String]? &mdash; Default: None</i><br />
        Fields to use for balanced sampling of TP examples to remove the prior distribution of the homozygous SNVs, the pipeline will attempt for the distribution over these arguments to be uniform. The recommended value is set in the template
</p>
<p name="SingleReadSNV.pre_filter">
        <b>SingleReadSNV.pre_filter</b><br />
        <i>String? &mdash; Default: None</i><br />
        SNV filter to apply to the FeatureMap before model, any SNV not passing the pre_filter is assigned QUAL=0, the recommended value is set in the template
</p>
</details>


### Advanced inputs
<details>
<summary> Show/Hide </summary>
<p name="SingleReadSNV.random_seed">
        <b>SingleReadSNV.random_seed</b><br />
        <i>Int &mdash; Default: None</i><br />
         Random seed to use for the ML model, the recommended value is set in the template 
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
<p name="SingleReadSNV.report_html">
        <b>SingleReadSNV.report_html</b><br />
        <i>File</i><br />
        SRSNV QC report html file
</p>
<p name="SingleReadSNV.model_file">
        <b>SingleReadSNV.model_file</b><br />
        <i>File</i><br />
        ML model file, saved with joblib
</p>
<p name="SingleReadSNV.X_test_file">
        <b>SingleReadSNV.X_test_file</b><br />
        <i>File</i><br />
        ML model test set features DataFrame, parquet format
</p>
<p name="SingleReadSNV.y_test_file">
        <b>SingleReadSNV.y_test_file</b><br />
        <i>File</i><br />
        ML model test set labels DataFrame, parquet format
</p>
<p name="SingleReadSNV.qual_test_file">
        <b>SingleReadSNV.qual_test_file</b><br />
        <i>File</i><br />
        ML model test set qual (SNVQ) DataFrame, parquet format
</p>
<p name="SingleReadSNV.X_train_file">
        <b>SingleReadSNV.X_train_file</b><br />
        <i>File</i><br />
        ML model training set features DataFrame, parquet format
</p>
<p name="SingleReadSNV.y_train_file">
        <b>SingleReadSNV.y_train_file</b><br />
        <i>File</i><br />
        ML model training set labels DataFrame, parquet format
</p>
<p name="SingleReadSNV.params_file">
        <b>SingleReadSNV.params_file</b><br />
        <i>File</i><br />
        ML model parameters json file
</p>
<p name="SingleReadSNV.test_set_mrd_simulation_dataframe">
        <b>SingleReadSNV.test_set_mrd_simulation_dataframe</b><br />
        <i>File</i><br />
        ML model test set MRD simulation DataFrame, parquet format
</p>
<p name="SingleReadSNV.train_set_mrd_simulation_dataframe">
        <b>SingleReadSNV.train_set_mrd_simulation_dataframe</b><br />
        <i>File</i><br />
        ML model training set MRD simulation DataFrame, parquet format
</p>
<p name="SingleReadSNV.test_set_statistics_h5">
        <b>SingleReadSNV.test_set_statistics_h5</b><br />
        <i>File</i><br />
        ML model test set statistics h5 file
</p>
<p name="SingleReadSNV.train_set_statistics_h5">
        <b>SingleReadSNV.train_set_statistics_h5</b><br />
        <i>File</i><br />
        ML model training set statistics h5 file
</p>
<p name="SingleReadSNV.train_set_statistics_json">
        <b>SingleReadSNV.train_set_statistics_json</b><br />
        <i>File</i><br />
        ML model training set statistics json file
</p>
<p name="SingleReadSNV.aggregated_metrics_json">
        <b>SingleReadSNV.aggregated_metrics_json</b><br />
        <i>File</i><br />
        ML model aggregated metrics json file
</p>
<p name="SingleReadSNV.tp_training_regions_bed">
        <b>SingleReadSNV.tp_training_regions_bed</b><br />
        <i>File</i><br />
        ML model training set TP regions bed file
</p>
<p name="SingleReadSNV.fp_training_regions_bed">
        <b>SingleReadSNV.fp_training_regions_bed</b><br />
        <i>File</i><br />
        ML model training set FP regions bed file
</p>
<p name="SingleReadSNV.train_report_file_notebook">
        <b>SingleReadSNV.train_report_file_notebook</b><br />
        <i>File</i><br />
        ML model training set report notebook file
</p>
<p name="SingleReadSNV.train_report_file_html">
        <b>SingleReadSNV.train_report_file_html</b><br />
        <i>File</i><br />
        ML model training set report html file
</p>
<p name="SingleReadSNV.test_report_file_notebook">
        <b>SingleReadSNV.test_report_file_notebook</b><br />
        <i>File</i><br />
        ML model test set report notebook file
</p>
<p name="SingleReadSNV.flow_order">
        <b>SingleReadSNV.flow_order</b><br />
        <i>String</i><br />
        Flow order for the sample
</p>

<hr />

> Generated using WDL AID (1.0.0)
