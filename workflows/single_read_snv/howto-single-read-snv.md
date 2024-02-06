# Single Read SNV (SRSNV) pipeline

## Introduction
This document describes the UG Single Read SNV (SRSNV) calling pipeline, which is a tool for assessing the quality of individual base substitutions, denoted SNVs for convenience, compared to the reference genome. Each SNV reported from a CRAM file is reported in a custom VCF file denoted as FeatureMap, using the gatk FlowFeatureMapper tool. The FeatureMap is a VCF file that contains a record for each SNV in each read, so that multiple entries per locus are possible, with additional information about the SNV and about the read encoded as INFO fields. Additionally, a machine learning model is trained on these features to assign an SNV quality score, saved as the QUAL field of each SNV in the FeatureMap. 

SNVQ is a more precise quality metric than BQ, calculated per alt rather than aggregated across the three options. BQ, the standard substitution quality metric, measures the likelihood of any ALT at a given locus, which is inherently lossy. SNVQ, on the other hand, measures the likelihood of each ALT at a given locus, which is more precise.

The SRSNV pipeline is composed of the following stages:
1. Featuremap is created (potentially done in parallel over genomic intervals and merged) using gatk FlowFeatureMapper
2. Featuremap is annotated with additional features (e.g. softclip length, reference context, etc.) 
3. Training the ML model on the annotated FeatureMap:
  
    a. FeatureMap is split into two groups - SNV supporting homozygous SNVs used as TP, and SNV supported by 1 read only in a high coverage locus used as FP
  
    b. ML classifier model trained, by default XGBoost is used
  
    c. Predicted probabilities from model are calibrated to an SNV quality score (estimated residual SNV rate)
4. Applying the model to all the entries, creating the final output FeatureMap with SNVQ values

This documents describes the pipeline, its inputs and outputs, and how to run it. 


## Naming convention
Arguments annotated with curly brackets below indicate an internal variable name, e.g. {ref_fasta} indicates the 
reference fasta file and {input_cram_bam} indicates the input CRAM (or BAM) file.

For a list of input separated by spaces, {sep=" " include_regions} is used.


## Reference data files

### hg38 reference genome
{ref_fasta}: gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta

{ref_fasta_fai}: gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai

{ref_dict}: gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict

or 

{ref_fasta}: s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta

{ref_fasta_fai}: s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai

{ref_dict}: s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dict

### Interval list - region in which FeatureMap is generated
{interval_list}: gs://concordanz/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list

or

{interval_list}: s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list

### Model training bed and vcf files

SNVs are collected and annotated in as either ground truth True Positives (TP) or False Positives (FP), and the SRSNV model is then trained to distinguish between the two. Genomic regions used for training data collection can be specified for TP and FP separately, allowing for some flexibility in application. By default both models are limited to the UG HCR and hmers of length 7 or more are discarded. Additionally, a curated subset of the PCAWG mutation database is excluded to allow an evalutation of the results over these regions, often used to demonstrate the error rate in a human WG sample [Cheng 2022]. Additionally, common population variants from the dbsnp and gnomAD databases are excluded from the FP training set, to avoid misidentified germline variants or contamination reads from skewing the ground truth data.

Genomic region filtering can be modified by the user (e.g. by excluding a specific chromosome) or effectively disabled by removing the exclude regions and setting the include regions to span the entire genome, e.g. gs://concordanz/hg38/wgs_calling_regions.hg38.bed or s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38.bed.

#### Regions used for TP training data

{include_regions_tp}:
- gs://concordanz/hg38/UG-High-Confidence-Regions/v2.1/ug_hcr.bed

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v2.1/ug_hcr.bed

{exclude_regions_tp}:
- gs://concordanz/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- gs://concordanz/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- s3://ultimagen-workflow-resources-us-east-1/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz

#### Regions used for FP training data

{include_regions_fp}:
- gs://concordanz/hg38/UG-High-Confidence-Regions/v2.1/ug_hcr.bed

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed

{exclude_regions_fp}:
- gs://concordanz/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- gs://concordanz/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz
- gs://concordanz/hg38/somatic/Homo_sapiens_assembly38.dbsnp138.chr1-22XY.snps.vcf.gz
- gs://concordanz/hg38/somatic/af-only-gnomad.hg38.snps.AF_over_1e-3.vcf.gz

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- s3://ultimagen-workflow-resources-us-east-1/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz
- s3://ultimagen-workflow-resources-us-east-1/hg38/somatic/Homo_sapiens_assembly38.dbsnp138.chr1-22XY.snps.vcf.gz
- s3://ultimagen-workflow-resources-us-east-1/hg38/somatic/af-only-gnomad.hg38.snps.AF_over_1e-3.vcf.gz

***Note 1 - for cfDNA samples from cancer patients, it's recommended to add the somatic mutation vcf (signature) to the fp exclude regions to avoid true cancer mutations present in the cfDNA being used as FP training examples***

***Note 2 - you can add any bed or vcf.gz file with loci of interest to exclude from model training***

***Note 3 - some files appear in two lists***


## Test data files

The files below are provided for the user to run the pipeline in their local environment to ensure it is working properly. A ppmSeq sample is provided for this purpose.

### Small test file - chr20 only
{input_cram_bam}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram
{input_cram_bam_index}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram.crai
{sorter_json_stats_file}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.json
{interval_list}: gs://ug-cromwell-tests/balanced_strand/wgs_calling_regions.hg38.chr20.interval_list

or

{input_cram_bam}: s3://ultimagen-workflow-resources-us-east-1/test_data/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram
{input_cram_bam_index}: s3://ultimagen-workflow-resources-us-east-1/test_data/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram.crai
{sorter_json_stats_file}: s3://ultimagen-workflow-resources-us-east-1/test_data/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.json
{interval_list}: s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38.chr20.interval_list

### Full test file
{input_cram_bam}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.cram
{input_cram_bam_index}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.cram.crai
{sorter_json_stats_file}: gs://ug-cromwell-tests/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.json
{interval_list}: gs://concordanz/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list

or

{input_cram_bam}: s3://ultimagen-feb-2024-mrd/cfDNA_ppmSeq/Pa_46.333_LuNgs_08.Lb_744.cram
{input_cram_bam_index}: s3://ultimagen-feb-2024-mrd/cfDNA_ppmSeq/Pa_46.333_LuNgs_08.Lb_744.cram.crai
{sorter_json_stats_file}: s3://ultimagen-feb-2024-mrd/cfDNA_ppmSeq/Pa_46.333_LuNgs_08.Lb_744.json
{interval_list}: s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list

## Running using separate command lines
### Installation

#### GATK

GATK can be downloaded or built according to the instructions on https://github.com/broadinstitute/gatk
Using a built jar file, denoted {gatk_jar}, running java with a prescribed memory is recommended, e.g.
{gatk}="java -Xms4g -jar {gatk_jar}"

For 4GB of memory. In the remainder of the document, {gatk} refers to the command running gatk.


#### UGVC repository

1. Clone the https://github.com/Ultimagen/VariantCalling repository (e.g. to `software/VariantCalling`)
2. Create conda environment:
`conda env create -f software/VariantCalling/setup/environment.yml`
This will create an environment called `genomics.py3`
3. Install the ugvc package:
```
conda activate genomics.py3
cd software/VariantCalling
pip install .
```
4. To run, cd to the installation path and use "python ugvc"

In the remainder of the document, {ugvc} refers to the command running ugvc from the correct environment, e.g.
{ugvc} = "conda run -n genomics.py3 python /path/to/ugvc"

## Variables

### Input files and names
* {input_cram_bam}: input cram file
* {input_cram_bam_index}: index of input cram file
* {sorter_json_stats_file}: json file with UG cram statistics generated jointly, same basename with a json extension
* {base_file_name}: the base file name of the output files and plot titles

### Main outputs
* {output_featuremap}: FeatureMap vcf.gz file with all the annotations and SNVQ values, and a respective .tbi index file
* {model_file}: the ML model file saved with joblib
* {params_file}: the ML model parameters file saved as json
* {test_set_statistics_h5}: statistics of the test set used for model evaluation, saved as h5
* {test_report_file_html}: html report of the test set used for model evaluation

### Intermediate/advanced outputs
* {featuremap}: FeatureMap vcf.gz file with a respective .tbi index file
* {annotated_featuremap}: FeatureMap vcf.gz file with additional annotations and a respective .tbi index file
* {single_substitutions_featuremap}: FeatureMap vcf.gz file of single substitutions with a respective .tbi index file
* {hom_snv_featuremap}: FeatureMap vcf.gz file of homozygous SNVs with a respective .tbi index file
* {training_regions_tp}/{training_regions_fp} - TP/FP bed files with the intersected and excluded regions to train the model

## Model training features and filter

These are provided either as WDL inputs in the json template or given as command line arguments to the ugvc srsnv_training command described in TrainSnvQualityRecalibrationModel.
* {pre_filter} - "(X_SCORE>4) && (X_EDIST<10)"
* {numerical_features} - numerical features used for training the model:
  - "X_SCORE" (Sequencing error likelihood)
  - "X_EDIST" (Levenshtein distance between the read and reference)
  - "X_FC1" (Number of SNVs in the read)
  - "X_LENGTH" (Read length)
  - "X_INDEX" (Index of the SNV in the read)
  - "rq" (Read quality)
  - "max_softclip_length" (Maximum softclip length)
  - "hmer_context_ref" (homopolymer length in the reference allele)
  - "hmer_context_alt" (homopolymer length in the alternative allele)
* {categorical_features} - categorical features used for training the model:
  - "is_cycle_skip", (Is the SNV a cycle skip)
  - "is_forward", (Is the read aligned to the forward strand, required to interpret reference related features)
  - "ref", (Reference base)
  - "alt", (Alternative base)
  - "next_1", (Reference base after the SNV)
  - "next_2", (Reference base 2bp after the SNV)
  - "next_3", (Reference base 3bp after the SNV)
  - "prev_1", (Reference base before the SNV)
  - "prev_2", (Reference base 2bp before the SNV)
  - "prev_3" (Reference base 3bp before the SNV)
  - "strand_ratio_category_start" (ppmSeq read category measured in the read start)
  - "strand_ratio_category_end" (ppmSeq read category measured in the read end)

***Note - the last two, strand_ratio_category_start and strand_ratio_category_end, are only required for ppmSeq data, remove otherwise***

## Running the SRSNV pipline

The SRSNV pipeline includes 4 modules:
1. FeatureMap - create featuremap from cram file, along with single substitution featuremap (FP)
2. BedIntersectAndExclude - intersect include and exclude regions, run for TP and FP each
3. TrainSnvQualityRecalibrationModel - training of ML model 
4. InferenceSnvQualityRecalibrationModel - inference on the featuremap (stage 1), using the ML model (stage 3).

### Featuremap

The Featuremap stage includes the following stages:
1. (optional) Split IntervalList 
2. Create FeatureMap
3. Annotate FeatureMap
4. Merge FeatureMap parts

#### Split IntervalList
Using gatk IntervalListTools, the interval list is split into smaller intervals, to allow parallel processing on so called "shards".
```
{gatk} \
  IntervalListTools \
  SCATTER_COUNT=50 \
  SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
  UNIQUE=true \
  SORT=true \
  BREAK_BANDS_AT_MULTIPLES_OF=10000 \
  INPUT={interval_list} \
  OUTPUT=out
```

*Recommended hardware - 1 CPU, 2GB RAM*

***Note - a scatter count of 50 is used, but can be changed according to the desired number of parallel tasks.***

##### Output files:
  out/\*/\*.interval_list - the split interval lists


#### Create FeatureMap:
The gatk FlowFeatureMapper tool is used to create the FeatureMap, a file that contains a record for each SNV in each read, along with additional information about the SNV or the read saved in the INFO field. Since each entry represents a single substitution with respect to the reference genome in a specific read, multiple entries per locus are possible, and a specific read can appear multiple times for multiple SNVs. 

For each interval_list file used, whether the full interval list or one of the shards (split interval lists), the following command is run: 
```
{gatk} \
FlowFeatureMapper -I {input_cram_bam} -O {featuremap} -R {ref_fasta}  \
--intervals {interval_list} \
--snv-identical-bases 5 \
--snv-identical-bases-after 5 \
--min-score 0 \
--limit-score 10 \
--read-filter MappingQualityReadFilter --minimum-mapping-quality 60 \
--flow-use-t0-tag --flow-fill-empty-bins-value 0.0001 --surrounding-median-quality-size 20  \
--copy-attr tm --copy-attr a3 --copy-attr rq \
--copy-attr as --copy-attr ts --copy-attr ae --copy-attr te --copy-attr s3 --copy-attr s2
```

*Recommended hardware - 1 CPU, 2GB RAM*

***Note 1 - {interval_list} is either the full interval list or one of the split interval lists***

***Note 2 - the last line (--copy-attr as[...]) is only required for Native Duplex data***

***Note 3 - A minimum mapping quality of 60 is used, but can be changed according to the desired mapping quality threshold***

##### Output files:
  {featuremap}: FeatureMap vcf.gz file with a respective .tbi index file

#### Annotate FeatureMap:
After creating the FeatureMap file, additional annotations are added to the INFO field of each entry. This is done per shard. The following command is run:
```
{ugvc} annotate_featuremap \
  -i {featuremap} \
  -o {annotated_featuremap} \
  --ref_fasta {ref_fasta} \
  --flow_order TGCA \
  --motif_length_to_annotate 3 \
  --max_hmer_length 20 \
  --balanced_strand_adapter_version "LA_v5and6" 
```

*Recommended hardware - 1 CPU, 4GB RAM*

***Note - the "balanced_strand_adapter_version" is required for ppmSeq data, remove it when running on different data. ***

##### Output files:
  {annotated_featuremap}: FeatureMap vcf.gz file with additional annotations and a respective .tbi index file

#### Merge FeatureMap parts
Only relevant if the interval list was split into shards in the previous stages.
```
bcftools concat --threads {cpus} -a -Oz -o {annotated_featuremap} {sep=" " featuremap_parts}
bcftools index -t {annotated_featuremap}
```

*Recommended hardware - 4 CPU, 4GB RAM*

***Note 1 - {featuremap_parts} is the list of annotated FeatureMap parts***

***Note 2 - {cpus} is the number of cpus to use for merging***

##### Output files:
  {annotated_featuremap}: FeatureMap vcf.gz file with additional annotations and a respective .tbi index file

#### Generate a FeatureMap of single substitutions:
Ae explained in the opening section, SNV supported by 1 read only in a high coverage locus used as FP in the model training data. To find positions where only one substitution is observed in the FeatureMap, and create a FeatureMap of these positions, this command is used:
```
bcftools view {annotated_featuremap} -H | vcf2bed | bedtools groupby -c 3 -o count | awk '($4==1) {print }' > single_sub_loci.bed
bcftools view -R single_sub_loci.bed -Oz -o {single_substitutions_featuremap} {annotated_featuremap}
bcftools index -t {single_substitutions_featuremap}
```

Filtering positions by coverage is done in the TrainSnvQualityRecalibrationModel stage when sampling FeatureMap entries.

*Recommended hardware - 1 CPU, 4GB RAM*

***Note - vcf2bed is a bedops tool that is installed as part of the genomics.py3 conda environment***

##### Output files:
  {single_substitutions_featuremap}: FeatureMap vcf.gz file of single substitutions with a respective .tbi index file

### CreateHomSnvFeatureMap

As explained in the opening section, SNV supporting homozygous SNVs are used as TP in the model training data. To find positions where only homozygous substitutions are observed in the FeatureMap, and create a FeatureMap of these positions, this command is used:
```
{ugvc} create_hom_snv_featuremap \
  --featuremap {annotated_featuremap} \
  --sorter_stats_json {sorter_json_stats_file} \
  --hom_snv_featuremap {hom_snv_featuremap} \
  --requested_min_coverage 20 \
  --min_af 0.7
```
*Recommended hardware - 1 CPU, 2GB RAM*

{requested_min_coverage} Is the minimum coverage requested for locus to be propagated to the output. If the median coverage is lower than this value, the median coverage (deduced from the {sorter_json_stats_file}) will be used as the minimum coverage instead. 
{min_af} is the minimum allele frequency required to consider a locus for hom SNV filtering. The default is chosen as 0.7 and not higher because some SNVs are pre-filtered from the FeatureMap due to MAPQ<60 or due to adjacent hmers.

##### Output files:
  {hom_snv_featuremap}: FeatureMap vcf.gz file of homozygous SNVs with a respective .tbi index file


### BedIntersectAndExclude
Prepare regions for training and inference by intersecting and excluding regions of interest. One region for FP entries and one region for TP entries are created. Note that these two commands can run in parallel.
```
{ugvc} intersect_bed_regions \
  --include-regions {sep=" " include_regions_tp} \
  --exclude-regions {sep=" " exclude_regions_tp} \
  --output-bed {training_regions_tp}

{ugvc} intersect_bed_regions \
  --include-regions {sep=" " include_regions_fp} \
  --exclude-regions {sep=" " exclude_regions_fp} \
  --output-bed {training_regions_fp}
```
*Recommended hardware - 1 CPU, 8GB RAM*

##### Output files:
  {training_regions_tp}, {training_regions_fp} - bed files with the intersected and excluded regions


### TrainSnvQualityRecalibrationModel
This code trains the ML model on the annotated FeatureMap, and produces a report and a model file.
```
{ugvc} srsnv_training \
--hom_snv_featuremap {hom_snv_featuremap} \
--single_substitution_featuremap {single_substitutions_featuremap} \
--output output_dir \
--basename {base_file_name} \
--flow_order "TGCA" \
--reference_fasta {ref_fasta}  \
--cram_stats_file {sorter_json_stats_file} \
--train_set_size 3000000 \
--test_set_size 300000 \
--hom_snv_regions {training_regions_tp} \
--single_sub_regions {training_regions_fp} \
--numerical_features {sep=" " numerical_features} \
--categorical_features {sep=" " categorical_features} \
--pre_filter {pre_filter} \
--random_seed 0 \
--balanced_sampling_info_fields trinuc_context_with_alt is_forward \
--balanced_strand_adapter_version "LA_v5and6" \
```
*Recommended hardware - 1 CPU, 16GB RAM*

##### Output files:
  The training code produces a model and parameters file to be used downstream, the data used, and a report. {test_report_file_html} and {test_set_statistics_h5} can be used to QC the model results, and the model files are used for inference.

  - {model_file}: "{base_file_name}.model.joblib"
  - {params_file}: "{base_file_name}.params.json"
  - {X_test_file}: "{base_file_name}.X_test.parquet"
  - {y_test_file}: "{base_file_name}.y_test.parquet"
  - {qual_test_file}: "{base_file_name}.qual_test.parquet"   
  - {X_train_file}: "{base_file_name}.X_train.parquet"
  - {y_train_file}: "{base_file_name}.y_train.parquet"    
  - {test_set_mrd_simulation_dataframe}: "{base_file_name}.test.df_mrd_simulation.parquet"
  - {train_set_mrd_simulation_dataframe}: "{base_file_name}.train.df_mrd_simulation.parquet"
  - {test_set_statistics_h5}: "{base_file_name}.test.statistics.h5"
  - {train_set_statistics_h5}: "{base_file_name}.train.statistics.h5"
  - {test_set_statistics_json}: "{base_file_name}.test.statistics.json"
  - {train_set_statistics_json}: "{base_file_name}.train.statistics.json"
  - {train_report_file_notebook}: "{base_file_name}.train_report.ipynb"
  - {train_report_file_html}: "{base_file_name}.train_report.html"
  - {test_report_file_notebook}: "{base_file_name}.test_report.ipynb"
  - {test_report_file_html}: "{base_file_name}.test_report.html"


### InferenceSnvQualityRecalibrationModel
This code applies the ML model to the annotated FeatureMap, and produces a FeatureMap with SNVQ values.
```
{ugvc} srsnv_inference \
--featuremap_path "{annotated_featuremap}" \
--params_path "{params_file}" \
--model_path "{model_file}" \
--X_train_path "{X_train_file}" \
--output_path "{output_featuremap}" \
--test_set_mrd_simulation_dataframe_file "{test_set_mrd_simulation_dataframe_file}" \
--process_number {cpus}
```

*Recommended hardware - 10 CPU, 8GB RAM*

***Note - {cpus} is the number of cpus to use***

##### Output files:
  {output_featuremap}: FeatureMap vcf.gz file with all the annotations and SNVQ values, and a respective .tbi index file


# References
Cheng, Alexandre Pellan, et al. "Whole genome error-corrected sequencing for sensitive circulating tumor DNA cancer monitoring." bioRxiv (2022): 2022-11.