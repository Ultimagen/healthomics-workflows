# SingleReadSNV
The Single Read SNV (SRSNV) pipeline is a read-centric de-noising framework, developed to overcome the limitations of traditional locus-centric variant calling, particularly in scenarios where rare mutations may be supported by only a single read. These rare mutations need to be distinguished from artefactual SNVs, which can derive from sequencing, library or alignment errors. To achieve this, we employed a supervised machine learning model trained to classify actual SNVs (labelled True or TP) from noise (False or FP). First, a comprehensive dataset capturing every candidate SNV is generated, along with a rich suite of annotations that describe sequencing quality, local sequence motifs, fragment-specific features, and locus-specific information. Randomly selected bases in the data matching the reference genome are collected and annotated as True, while low VAF (<=5%) SNVs in high-coverage (>=20x) regions (SNVs with low support, indicating they are likely to be artifacts) are annotated as False SNVs. Using these curated sets, we train an XGBoost classifier to robustly distinguish between true and artifactual SNVs. Once trained, the classifier assigns a calibrated quality score to each SNV in the input CRAM, providing a precise estimate of the residual error rate. To avoid overfitting, an ensemble of models (3) are trained on different sets of chromosomes and applied using a cross-validation scheme. 

The following input templates are available for different kinds of input data: 

1) `single_read_snv_template-ppmSeq.json` | Use this template for ppmSeq data. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. st, et). 

2) `single_read_snv_template-ppmSeq_legacy_v5.json` | Use this template for LEGACY v5 ppmSeq data. This is an older version of the ppmSeq adapters, generally not available since 2024. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. as, ts). 

3) `single_read_snv_template-Standard-WG.json` | Use this template for any non-ppmSeq data. 

## Inputs

### Required inputs
<p name="SingleReadSNV.input_cram_bam_list">
        <b>SingleReadSNV.input_cram_bam_list</b><br />
        <i>Array[File] </i> &mdash; 
         Input CRAM file/s <br /> 
</p>
<p name="SingleReadSNV.input_cram_bam_index_list">
        <b>SingleReadSNV.input_cram_bam_index_list</b><br />
        <i>Array[File] </i> &mdash; 
         Input CRAM index file/s <br /> 
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
<p name="SingleReadSNV.sorter_json_stats_file_list">
        <b>SingleReadSNV.sorter_json_stats_file_list</b><br />
        <i>Array[File]? </i> &mdash; 
         (Optional) Sorter json stats files. Provide EITHER these files OR both mean_coverage and total_aligned_bases. <br /> 
</p>
<p name="SingleReadSNV.random_sample_trinuc_freq">
        <b>SingleReadSNV.random_sample_trinuc_freq</b><br />
        <i>File? </i> &mdash; 
         (Optional) CSV or TSV file with trinucleotide frequencies for the random sample. If provided, the random sample featuremap will be sampled according to the given trinucleotide frequency. If not provided, sampling is uniform. <br /> 
</p>
<p name="SingleReadSNV.create_md5_checksum_outputs">
        <b>SingleReadSNV.create_md5_checksum_outputs</b><br />
        <i>Boolean </i> &mdash; 
         Create md5 checksum for requested output files <br /> 
</p>
<p name="SingleReadSNV.mean_coverage">
        <b>SingleReadSNV.mean_coverage</b><br />
        <i>Float? </i> &mdash; 
         (Optional) Mean coverage value. Provide together with total_aligned_bases and without sorter_json_stats_file_list. <br /> 
</p>
<p name="SingleReadSNV.total_aligned_bases">
        <b>SingleReadSNV.total_aligned_bases</b><br />
        <i>String? </i> &mdash; 
         (Optional) Total aligned bases used for downsampling rate calculation. Provide together with mean_coverage and without sorter_json_stats_file_list. <br /> 
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
<p name="SingleReadSNV.random_sample_trinuc_freq_stats">
        <b>SingleReadSNV.random_sample_trinuc_freq_stats</b><br />
        <i>File?</i><br />
        (Optional) CSV or TSV file with trinucleotide frequencies for the random sample.
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
