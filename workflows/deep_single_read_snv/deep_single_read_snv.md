# DeepSingleReadSNV
Deep Single Read SNV (DeepSRSNV) pipeline. A convolutional neural network on read tensors (pileup images) with k-fold cross-validation, driven by a `mode` selector: 'full' (train + inference, default), 'train_only' (produce trained fold models + QC report, no inference), 'inference_only' (apply a provided pre-trained N-fold model to produce a DNN quality-annotated featuremap VCF, optionally skipping CreateFeatureMap when a featuremap VCF is supplied), and 'data_prep_only' (produce per-CRAM training tensor caches and stop, for later pooled multi-CRAM training). Requires GPU-enabled execution environment for training/inference modes.

## Inputs

### Required inputs
<p name="DeepSingleReadSNV.input_cram_bam">
        <b>DeepSingleReadSNV.input_cram_bam</b><br />
        <i>File </i> &mdash;
         Input CRAM file <br />
</p>
<p name="DeepSingleReadSNV.input_cram_bam_index">
        <b>DeepSingleReadSNV.input_cram_bam_index</b><br />
        <i>File </i> &mdash;
         Input CRAM index file <br />
</p>
<p name="DeepSingleReadSNV.base_file_name">
        <b>DeepSingleReadSNV.base_file_name</b><br />
        <i>String </i> &mdash;
         Base file name for output files <br />
</p>
<p name="DeepSingleReadSNV.mode">
        <b>DeepSingleReadSNV.mode</b><br />
        <i>String </i> &mdash;
         Pipeline mode (set in the input template): 'full' (train + inference), 'train_only' (training + QC report), 'inference_only' (apply provided model, requires inference_models), or 'data_prep_only' (produce training tensor caches and stop). <br />
</p>
<p name="DeepSingleReadSNV.num_shards_featuremap">
        <b>DeepSingleReadSNV.num_shards_featuremap</b><br />
        <i>Int </i> &mdash;
         Number of genomic shards to scatter the snvfind (CreateFeatureMap) step across. Higher values reduce wall-clock time but add scatter overhead. <br />
</p>
<p name="DeepSingleReadSNV.scatter_interval_list">
        <b>DeepSingleReadSNV.scatter_interval_list</b><br />
        <i>File </i> &mdash;
         Interval list defining the genomic regions to scatter snvfind across. Should match the regions in featuremap_params.bed_file. <br />
</p>

### Required parameters
<p name="DeepSingleReadSNV.featuremap_params">
        <b>DeepSingleReadSNV.featuremap_params</b><br />
        <i>FeatureMapParams </i> &mdash;
         FeatureMap parameters, recommended value set in the template. <br />
</p>
<p name="DeepSingleReadSNV.single_read_snv_params">
        <b>DeepSingleReadSNV.single_read_snv_params</b><br />
        <i>SingleReadSNVParams </i> &mdash;
         SingleReadSNV parameters for training set preparation. <br />
</p>
<p name="DeepSingleReadSNV.features">
        <b>DeepSingleReadSNV.features</b><br />
        <i>Array[String] </i> &mdash;
         Features to be used for the report quality plots, should match the XGBoost feature set. <br />
</p>
<p name="DeepSingleReadSNV.deep_srsnv_params">
        <b>DeepSingleReadSNV.deep_srsnv_params</b><br />
        <i>DeepSRSNVParams </i> &mdash;
         Deep SRSNV (DNN) parameters. Controls training hyperparameters, fold count, GPU resources, and inference backend. <br />
</p>
<p name="DeepSingleReadSNV.exclude_from_training_field_name">
        <b>DeepSingleReadSNV.exclude_from_training_field_name</b><br />
        <i>String </i> &mdash;
         INFO field name for the exclude-from-training annotation in the featuremap VCF. <br />
</p>
<p name="DeepSingleReadSNV.include_in_inference_field_name">
        <b>DeepSingleReadSNV.include_in_inference_field_name</b><br />
        <i>String </i> &mdash;
         INFO field name for the include-in-inference annotation in the featuremap VCF. <br />
</p>
<p name="DeepSingleReadSNV.pcawg_field_name">
        <b>DeepSingleReadSNV.pcawg_field_name</b><br />
        <i>String </i> &mdash;
         INFO field name for the PCAWG annotation in the featuremap VCF. <br />
</p>
<p name="DeepSingleReadSNV.include_vcf_bcftools_filter_args">
        <b>DeepSingleReadSNV.include_vcf_bcftools_filter_args</b><br />
        <i>String </i> &mdash;
         Bcftools filter arguments applied to include-in-inference VCFs before annotation. <br />
</p>

### Required references
<p name="DeepSingleReadSNV.annotation_files">
        <b>DeepSingleReadSNV.annotation_files</b><br />
        <i>FeaturemapAnnotationFiles </i> &mdash;
         Annotation files for featuremap generation: dbSNP, gnomAD, and UG High Confidence Regions with their indices <br />
</p>

### Optional inputs
<details>
<summary> Show/Hide </summary>
<p name="DeepSingleReadSNV.raise_exceptions_in_report">
        <b>DeepSingleReadSNV.raise_exceptions_in_report</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Raise an exception and fail the pipeline if an error is raised in the QC report
</p>
<p name="DeepSingleReadSNV.override_memory_gb_CreateFeatureMap">
        <b>DeepSingleReadSNV.override_memory_gb_CreateFeatureMap</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Override memory in GB for the CreateFeatureMap task. If an out of memory error occurs, try increasing this value.
</p>
<p name="DeepSingleReadSNV.override_memory_gb_PrepareRawFeatureMap">
        <b>DeepSingleReadSNV.override_memory_gb_PrepareRawFeatureMap</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Override memory in GB for the PrepareRawFeatureMap task. If an out of memory error occurs, try increasing this value.
</p>
<p name="DeepSingleReadSNV.override_memory_gb_PrepareRandomSampleFeatureMap">
        <b>DeepSingleReadSNV.override_memory_gb_PrepareRandomSampleFeatureMap</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Override memory in GB for the PrepareRandomSampleFeatureMap task. If an out of memory error occurs, try increasing this value.
</p>

### Optional inputs
<p name="DeepSingleReadSNV.sorter_json_stats_file_list">
        <b>DeepSingleReadSNV.sorter_json_stats_file_list</b><br />
        <i>Array[File]? </i> &mdash;
         (Optional) Sorter json stats files. Provide EITHER these files OR both mean_coverage and total_aligned_bases. <br />
</p>
<p name="DeepSingleReadSNV.reference_genome">
        <b>DeepSingleReadSNV.reference_genome</b><br />
        <i>String </i> &mdash;
         Genome type selector. The workflow currently supports only hg38. <br />
</p>
<p name="DeepSingleReadSNV.random_sample_trinuc_freq">
        <b>DeepSingleReadSNV.random_sample_trinuc_freq</b><br />
        <i>File? </i> &mdash;
         (Optional) CSV or TSV file with trinucleotide frequencies for the random sample. <br />
</p>
<p name="DeepSingleReadSNV.inference_models">
        <b>DeepSingleReadSNV.inference_models</b><br />
        <i>DeepSRSNVModel? </i> &mdash;
         Pre-trained N-fold model for mode='inference_only'. fold_metadata (recalibrated), fold_checkpoints, and fold_onnx_models are REQUIRED (engine rebuilt from ONNX in-runtime); fold_engines optional. Fold count inferred from the arrays. <br />
</p>
<p name="DeepSingleReadSNV.input_featuremap_vcf">
        <b>DeepSingleReadSNV.input_featuremap_vcf</b><br />
        <i>File? </i> &mdash;
         (inference_only) Pre-computed FeatureMap VCF. If provided, CreateFeatureMap is skipped and inference runs directly on it. <br />
</p>
<p name="DeepSingleReadSNV.input_featuremap_vcf_index">
        <b>DeepSingleReadSNV.input_featuremap_vcf_index</b><br />
        <i>File? </i> &mdash;
         (inference_only) Index for input_featuremap_vcf. Required if input_featuremap_vcf is provided. <br />
</p>
<p name="DeepSingleReadSNV.mean_coverage">
        <b>DeepSingleReadSNV.mean_coverage</b><br />
        <i>Float? </i> &mdash;
         (Optional) Mean coverage value. Provide together with total_aligned_bases and without sorter_json_stats_file_list. <br />
</p>
<p name="DeepSingleReadSNV.total_aligned_bases">
        <b>DeepSingleReadSNV.total_aligned_bases</b><br />
        <i>String? </i> &mdash;
         (Optional) Total aligned bases used for downsampling rate calculation. <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.base_file_name">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.base_file_name</b><br />
        <i>String? </i> &mdash;
         Base file name for output files. <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.exclude_regions">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.exclude_regions</b><br />
        <i>Array[File]? </i> &mdash;
         Regions to exclude from the output vcf. Supported formats are bed, bed.gz, vcf, vcf.gz. VCF exclusion is done using bcftools view -T ^regions, by position and not by ref and alt. <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.include_regions">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.include_regions</b><br />
        <i>Array[File]? </i> &mdash;
         Regions to include in the output vcf. Supported formats are bed, bed.gz. Inclusion is done using 'bcftools view -T'. <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.disk_size">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.disk_size</b><br />
        <i>Int </i> &mdash;
         Size of the local disk to use for this task, in GB. By default it is calculated from the input file sizes. <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.memory_gb">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.memory_gb</b><br />
        <i>Int </i> &mdash;
         Amount of memory to use for this task, in GB. Default is 4 (GB). <br />
</p>
<p name="DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.cpus">
        <b>DeepSingleReadSNV.FeatureMapPrep.FilterIncludeVcf.cpus</b><br />
        <i>Int </i> &mdash;
         Number of cpus to use for this task. Default is 4. <br />
</p>
</details>


## Outputs
<p name="DeepSingleReadSNV.featuremap">
        <b>DeepSingleReadSNV.featuremap</b><br />
        <i>File?</i><br />
        FeatureMap VCF with all SNV candidates
</p>
<p name="DeepSingleReadSNV.featuremap_index">
        <b>DeepSingleReadSNV.featuremap_index</b><br />
        <i>File?</i><br />
        Index for the FeatureMap VCF
</p>
<p name="DeepSingleReadSNV.featuremap_random_sample">
        <b>DeepSingleReadSNV.featuremap_random_sample</b><br />
        <i>File?</i><br />
        Downsampled FeatureMap VCF file for training
</p>
<p name="DeepSingleReadSNV.featuremap_random_sample_index">
        <b>DeepSingleReadSNV.featuremap_random_sample_index</b><br />
        <i>File?</i><br />
        Downsampled FeatureMap VCF index file
</p>
<p name="DeepSingleReadSNV.random_sample_trinuc_freq_stats">
        <b>DeepSingleReadSNV.random_sample_trinuc_freq_stats</b><br />
        <i>File?</i><br />
        Trinucleotide frequency statistics from the random sample featuremap
</p>
<p name="DeepSingleReadSNV.downsampling_rate">
        <b>DeepSingleReadSNV.downsampling_rate</b><br />
        <i>Float?</i><br />
        The downsampling rate used to create the random sample featuremap
</p>
<p name="DeepSingleReadSNV.positive_parquet">
        <b>DeepSingleReadSNV.positive_parquet</b><br />
        <i>File?</i><br />
        Positive-label (random-sample) training featuremap parquet (training/data_prep modes)
</p>
<p name="DeepSingleReadSNV.negative_parquet">
        <b>DeepSingleReadSNV.negative_parquet</b><br />
        <i>File?</i><br />
        Negative-label (raw) training featuremap parquet (training/data_prep modes)
</p>
<p name="DeepSingleReadSNV.positive_tensor_cache_tar">
        <b>DeepSingleReadSNV.positive_tensor_cache_tar</b><br />
        <i>File?</i><br />
        Positive-label read tensor cache tar for pooled multi-CRAM training (data_prep_only mode only)
</p>
<p name="DeepSingleReadSNV.negative_tensor_cache_tar">
        <b>DeepSingleReadSNV.negative_tensor_cache_tar</b><br />
        <i>File?</i><br />
        Negative-label read tensor cache tar for pooled multi-CRAM training (data_prep_only mode only)
</p>
<p name="DeepSingleReadSNV.stats_funnel">
        <b>DeepSingleReadSNV.stats_funnel</b><br />
        <i>File?</i><br />
        snvfind model-filters status funnel JSON, needed as --stats-file for pooled multi-CRAM training (data_prep_only mode only)
</p>
<p name="DeepSingleReadSNV.mean_coverage_out">
        <b>DeepSingleReadSNV.mean_coverage_out</b><br />
        <i>Float?</i><br />
        Resolved mean coverage, needed as --mean-coverage for pooled multi-CRAM training (data_prep_only mode only)
</p>
<p name="DeepSingleReadSNV.featuremap_vcf">
        <b>DeepSingleReadSNV.featuremap_vcf</b><br />
        <i>File?</i><br />
        FeatureMap VCF annotated with DNN quality scores (inference modes)
</p>
<p name="DeepSingleReadSNV.featuremap_vcf_index">
        <b>DeepSingleReadSNV.featuremap_vcf_index</b><br />
        <i>File?</i><br />
        Index for the annotated FeatureMap VCF (inference modes)
</p>
<p name="DeepSingleReadSNV.fold_metadata">
        <b>DeepSingleReadSNV.fold_metadata</b><br />
        <i>Array[File]?</i><br />
        Per-fold metadata JSON files from DNN training
</p>
<p name="DeepSingleReadSNV.fold_checkpoints">
        <b>DeepSingleReadSNV.fold_checkpoints</b><br />
        <i>Array[File]?</i><br />
        Per-fold model checkpoint files from DNN training
</p>
<p name="DeepSingleReadSNV.fold_onnx_models">
        <b>DeepSingleReadSNV.fold_onnx_models</b><br />
        <i>Array[File]?</i><br />
        Per-fold ONNX model files from DNN training
</p>
<p name="DeepSingleReadSNV.fold_engines">
        <b>DeepSingleReadSNV.fold_engines</b><br />
        <i>Array[File]?</i><br />
        Per-fold TensorRT engine files from DNN training
</p>
<p name="DeepSingleReadSNV.fold_trt_timing_caches">
        <b>DeepSingleReadSNV.fold_trt_timing_caches</b><br />
        <i>Array[File]?</i><br />
        Per-fold TensorRT timing caches from DNN training; feed back as DeepSRSNVModel.fold_timing_caches in a later inference_only run to rebuild bit-identical engines (same GPU + TRT version)
</p>
<p name="DeepSingleReadSNV.combined_metadata">
        <b>DeepSingleReadSNV.combined_metadata</b><br />
        <i>File?</i><br />
        Shared quality recalibration LUT metadata
</p>
<p name="DeepSingleReadSNV.updated_fold_metadata">
        <b>DeepSingleReadSNV.updated_fold_metadata</b><br />
        <i>Array[File]?</i><br />
        Per-fold metadata JSONs updated with the shared recalibration LUT (training modes only)
</p>
<p name="DeepSingleReadSNV.featuremap_df">
        <b>DeepSingleReadSNV.featuremap_df</b><br />
        <i>File?</i><br />
        Combined featuremap DataFrame (parquet) with per-fold predictions
</p>
<p name="DeepSingleReadSNV.report_html">
        <b>DeepSingleReadSNV.report_html</b><br />
        <i>File?</i><br />
        QC report HTML file
</p>
<p name="DeepSingleReadSNV.application_qc_h5">
        <b>DeepSingleReadSNV.application_qc_h5</b><br />
        <i>File?</i><br />
        Application QC statistics h5 file
</p>

<hr />

> Generated using WDL AID (1.0.1)