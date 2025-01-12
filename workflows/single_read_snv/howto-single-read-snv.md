# Single Read SNV (SRSNV) pipeline

## Table of Contents
- [Single Read SNV (SRSNV) pipeline](#single-read-snv-srsnv-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Naming convention](#naming-convention)
  - [Reference data files](#reference-data-files)
    - [hg38 reference genome](#hg38-reference-genome)
    - [Interval list - region in which FeatureMap is generated](#interval-list---region-in-which-featuremap-is-generated)
    - [Model training bed and vcf files](#model-training-bed-and-vcf-files)
      - [Regions used for TP training data](#regions-used-for-tp-training-data)
      - [Regions used for FP training data](#regions-used-for-fp-training-data)
  - [Variables](#variables)
    - [Input files and names](#input-files-and-names)
    - [Main outputs](#main-outputs)
    - [Intermediate/advanced outputs](#intermediateadvanced-outputs)
  - [Running the SRSNV pipline without the WDL](#running-the-srsnv-pipline-without-the-wdl)
    - [Fetching dockers](#fetching-dockers)
    - [Manual installation](#manual-installation)
      - [GATK and Picard](#gatk-and-picard)
      - [ugbio-utils repository](#ugbio-utils-repository)
    - [Pipeline structure](#pipeline-structure)
    - [Featuremap](#featuremap)
      - [Split IntervalList](#split-intervallist)
        - [Output files:](#output-files)
      - [Create FeatureMap:](#create-featuremap)
        - [Output files:](#output-files-1)
      - [Annotate FeatureMap:](#annotate-featuremap)
        - [Output files:](#output-files-2)
      - [Merge FeatureMap parts](#merge-featuremap-parts)
        - [Output files:](#output-files-3)
      - [Generate a FeatureMap of single substitutions:](#generate-a-featuremap-of-single-substitutions)
        - [Output files:](#output-files-4)
    - [CreateHomSnvFeatureMap](#createhomsnvfeaturemap)
        - [Output files:](#output-files-5)
    - [BedIntersectAndExclude](#bedintersectandexclude)
        - [Output files:](#output-files-6)
    - [TrainSnvQualityRecalibrationModel](#trainsnvqualityrecalibrationmodel)
        - [Output files:](#output-files-7)
    - [InferenceSnvQualityRecalibrationModel](#inferencesnvqualityrecalibrationmodel)
        - [Output files:](#output-files-8)
  - [Detailed explanation of keys output files](#detailed-explanation-of-keys-output-files)
    - [output\_featuremap](#output_featuremap)
    - [featuremap\_df\_file](#featuremap_df_file)
    - [Model joblib file](#model-joblib-file)
  - [Model training features and filter](#model-training-features-and-filter)
  - [Details of ML model Cross-Validation scheme](#details-of-ml-model-cross-validation-scheme)
    - [Training with train/test split](#training-with-traintest-split)
    - [Training with k-fold CV](#training-with-k-fold-cv)
  - [References](#references)


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
{interval_list}: gs://concordanz/hg38/wgs_calling_regions.without_encode_blacklist.interval_list

or

{interval_list}: s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.without_encode_blacklist.interval_list

### Model training bed and vcf files

SNVs are collected and annotated in as either ground truth True Positives (TP) or False Positives (FP), and the SRSNV model is then trained to distinguish between the two. Genomic regions used for training data collection can be specified for TP and FP separately, allowing for some flexibility in application. By default both models are limited to the UG HCR and hmers of length 7 or more are discarded. Additionally, a curated subset of the PCAWG mutation database is excluded to allow an evalutation of the results over these regions, often used to demonstrate the error rate in a human WG sample [Cheng 2022]. Additionally, common population variants from the dbsnp and gnomAD databases are excluded from the FP training set, to avoid misidentified germline variants or contamination reads from skewing the ground truth data.

Genomic region filtering can be modified by the user (e.g. by excluding a specific chromosome) or effectively disabled by removing the exclude regions and setting the include regions to span the entire genome, e.g. gs://concordanz/hg38/wgs_calling_regions.hg38.bed or s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38.bed.

#### Regions used for TP training data

{include_regions_tp}:
- gs://concordanz/hg38/UG-High-Confidence-Regions/v2.1.2/ug_hcr_autosomal.bed

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v2.1.2/ug_hcr_autosomal.bed

{exclude_regions_tp}:
- gs://concordanz/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- gs://concordanz/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz

or

- s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/hmers_7_and_higher.chr1-22XY.bed
- s3://ultimagen-workflow-resources-us-east-1/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz

#### Regions used for FP training data

{include_regions_fp}:
- gs://concordanz/hg38/UG-High-Confidence-Regions/v2.1.2/ug_hcr.bed

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

***Note 1 - for cfDNA samples from cancer patients, it is recommended to add the somatic mutation vcf (signature) to the fp exclude regions to avoid true cancer mutations present in the cfDNA being used as FP training examples***

***Note 2 - you can add any bed or vcf.gz file with loci of interest to exclude from model training***

***Note 3 - some files appear in two lists***


## Variables

### Input files and names
* {input_cram_bam}: input cram file
* {input_cram_bam_index}: index of input cram file
* {sorter_json_stats_file}: json file with UG cram statistics generated jointly, same basename with a json extension
* {base_file_name}: the base file name of the output files and plot titles

### Main outputs
* {output_featuremap}: FeatureMap vcf.gz file with all the annotations and SNVQ values, and a respective .tbi index file
* {model_file}: the ML model(s) and metadata required for inference, saved with joblib
* {params_file}: the ML model parameters file saved as json
* {test_set_statistics_h5}: statistics of the test set used for model evaluation, saved as h5
* {test_report_file_html}: html report of the test set used for model evaluation

### Intermediate/advanced outputs
* {featuremap}: FeatureMap vcf.gz file with a respective .tbi index file
* {annotated_featuremap}: FeatureMap vcf.gz file with additional annotations and a respective .tbi index file
* {single_substitutions_featuremap}: FeatureMap vcf.gz file of single substitutions with a respective .tbi index file
* {hom_snv_featuremap}: FeatureMap vcf.gz file of homozygous SNVs with a respective .tbi index file
* {training_regions_tp}/{training_regions_fp} - TP/FP bed files with the intersected and excluded regions to train the model


## Running the SRSNV pipline without the WDL
It is possible to run the same code as in the SingleReadSNV workflow (https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_read_snv/single_read_snv.wdl) without using the WDL. This section details the required steps and installation methods, lifted from the code in the WDL tasks. Readers fluent in the WDL language might prefer to use the code in the WDL itself and the various imported tasks in https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/single_read_snv/tasks, to ensure the exact same code runs in either case. 
Running the pipeline is possible either using the dockers (recommended), or by installing the required environments manually, as detailed below. 

### Fetching dockers
The docker required in this workflow are:
- ug_gatk_picard_docker
- broad_gatk_docker
- ugbio_srsnv_docker
- ugbio_featuremap_docker
- ugbio_vcflite_docker
- gitc_docker
Pull each docker according to the latest version referred to in https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_read_snv/tasks/globals.wdl

### Manual installation

#### GATK and Picard
For code running in "broad_gatk_docker", GATK can be downloaded or built according to the instructions on https://github.com/broadinstitute/gatk
Using a built jar file, denoted {gatk_jar}, running java with a prescribed memory is recommended, e.g.
{gatk}="java -Xms4g -jar {gatk_jar}"
For 4GB of memory.

For code running in "ug_gatk_picard_docker", refer to the UG gatk fork https://github.com/Ultimagen/gatk instead, and proceed using the same instructions.

For code running in "gitc_docker", follow the same instructions in https://github.com/broadinstitute/picard

#### ugbio-utils repository
ugbio-utils
1. Clone the https://github.com/Ultimagen/ugbio-utils repository
2. Install the environment using uv according to the repository instructions
3. To run code relying on the "ugbio_featuremap_docker", run "uv sync --package ugbio_featuremap" before code execution
4. To run code relying on the "ugbio_srsnv_docker", run "uv sync --package ugbio_srsnv" before code execution
5. To run code relying on the "vcflite_docker", run "uv sync --package ugbio_vcflite" before code execution


### Pipeline structure
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
Docker = broad_gatk_docker
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

For each interval_list file used, whether the full interval list or one of the shards (split interval lists), the following command is run. Note that the second command is needed when running in parallel across many shards to avoid duplication of entries around the interval edges: 
Docker = ug_gatk_picard_docker
```
{gatk} \
  FlowFeatureMapper -I {input_cram_bam} -O tmp.vcf.gz -R {ref_fasta}  \
  --intervals {interval_list} \
  --snv-identical-bases 5 \
  --snv-identical-bases-after 5 \
  --min-score 0 \
  --limit-score 10 \
  --read-filter MappingQualityReadFilter --minimum-mapping-quality 60 \
  --flow-use-t0-tag --flow-fill-empty-bins-value 0.0001 --surrounding-median-quality-size 20  \
  --copy-attr tm --copy-attr a3 --copy-attr rq --copy-attr st --copy-attr et \
  --copy-attr as --copy-attr ts --copy-attr ae --copy-attr te --copy-attr s3 --copy-attr s2

{gatk} \
  VariantFiltration \
  -V tmp.vcf.gz \
  -O {featuremap} \
  -R {ref_fasta}  \
  --intervals {interval_list}


java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{memory_gb-2}g -jar ~{gitc_path}GATK_ultima.jar  \
    FlowFeatureMapper \
    -I "~{input_cram_bam}" \
    -O tmp.vcf.gz \
    -R "~{references.ref_fasta}"  \
    --intervals "~{interval_list}" \
    --snv-identical-bases ~{featuremap_params.snv_identical_bases} \
    --snv-identical-bases-after ~{featuremap_params.snv_identical_bases_after} \
    --min-score ~{featuremap_params.min_score} \
    --limit-score ~{featuremap_params.limit_score} \
    --read-filter MappingQualityReadFilter --minimum-mapping-quality ~{featuremap_params.min_mapq} \
    ~{featuremap_params.extra_args} 

  echo "***************************** Filtering on region *****************************"
  java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{memory_gb-2}g -jar ~{gitc_path}GATK_ultima.jar  \
    VariantFiltration \
    -V tmp.vcf.gz \
    -O "~{output_basename}.vcf.gz" \
    -R "~{references.ref_fasta}"  \
    --intervals "~{interval_list}"
```

*Recommended hardware - 1 CPU, 2GB RAM*

***Note 1 - {interval_list} is either the full interval list or one of the split interval lists***

***Note 2 - the last line (--copy-attr as[...]) is only required for Native Duplex data***

***Note 3 - A minimum mapping quality of 60 is used, but can be changed according to the desired mapping quality threshold***

##### Output files:
  {featuremap}: FeatureMap vcf.gz file with a respective .tbi index file

#### Annotate FeatureMap:
After creating the FeatureMap file, additional annotations are added to the INFO field of each entry. This is done per shard. The following command is run:
Docker = ugbio_featuremap_docker
```
annotate_featuremap \
  -i {featuremap} \
  -o {annotated_featuremap} \
  --ref_fasta {ref_fasta} \
  --flow_order TGCA \
  --motif_length_to_annotate 3 \
  --max_hmer_length 20 \
  --ppmSeq_adapter_version "ppmSeq" 
```

*Recommended hardware - 1 CPU, 4GB RAM*

***Note - the "ppmSeq_adapter_version" is required for ppmSeq data, remove it when running on non-ppmSeq data.***

##### Output files:
  {annotated_featuremap}: FeatureMap vcf.gz file with additional annotations and a respective .tbi index file

#### Merge FeatureMap parts
Only relevant if the interval list was split into shards in the previous stages.
Docker = ugbio_featuremap_docker
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
Docker = ugbio_vcflite_docker
```
vcflite import --vcf-in {annotated_featuremap}
vcflite query --group-by "chrom, pos" --having "COUNT(*) = 1" --vcf-out {single_substitutions_featuremap}
tabix -p vcf {single_substitutions_featuremap}
```

Filtering positions by coverage is done in the TrainSnvQualityRecalibrationModel stage when sampling FeatureMap entries.

*Recommended hardware - 1 CPU, 4GB RAM*

##### Output files:
  {single_substitutions_featuremap}: FeatureMap vcf.gz file of single substitutions with a respective .tbi index file

### CreateHomSnvFeatureMap

As explained in the opening section, SNV supporting homozygous SNVs are used as TP in the model training data. To find positions where only homozygous substitutions are observed in the FeatureMap, and create a FeatureMap of these positions, this command is used:
Docker = ugbio_srsnv_docker
```
create_hom_snv_featuremap \
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
Docker = ugbio_srsnv_docker
```
intersect_bed_regions \
  --include-regions {sep=" " include_regions_fp} \
  --exclude-regions {sep=" " exclude_regions_fp} \
  --output-bed {training_regions_fp}
```
*Recommended hardware - 1 CPU, 8GB RAM*

##### Output files:
  {training_regions_tp}, {training_regions_fp} - bed files with the intersected and excluded regions


### TrainSnvQualityRecalibrationModel
This code trains the ML model on the annotated FeatureMap, and produces a report and a model file. 
Docker = ugbio_srsnv_docker
```
# Create a json file with parameters for the model
srsnv_training \
  --hom_snv_featuremap {hom_snv_featuremap} \
  --single_substitution_featuremap {single_substitutions_featuremap} \
  --dataset_params_json_path {single_read_snv_params} \
  --flow_order "TGCA" \
  --reference_fasta "~{references.ref_fasta}" \
  --reference_dict "~{references.ref_dict}" \
  --cram_stats_file "~{sorter_json_stats_file}" \
  --hom_snv_regions "~{hom_snv_regions_bed}" \
  --single_sub_regions "~{single_substitution_regions_bed}" \
  --output "$PWD" \
  --basename "~{basename}"
```
*Recommended hardware - 1 CPU, 16GB RAM*

##### Output files:
  The training code produces a model and parameters file to be used downstream, the data used, and a report. {test_report_file_html} and {test_set_statistics_h5} can be used to QC the model results, and the model files are used for inference.

  - {model_file}: "{base_file_name}.model.joblib"
  - {params_file}: "{base_file_name}.params.json"
  - {featuremap_df_file}: "{base_file_name}.featuremap_df.parquet"    
  - {test_set_statistics_h5}: "{base_file_name}.test.statistics.h5"
  - {test_set_statistics_json}: "{base_file_name}.test.statistics.json"
  - {test_report_file_notebook}: "{base_file_name}.test_report.ipynb"
  - {test_report_file_html}: "{base_file_name}.test_report.html"


### InferenceSnvQualityRecalibrationModel
This code applies the ML model to the annotated FeatureMap, and produces a FeatureMap with SNVQ values.
Docker = ugbio_srsnv_docker
```
srsnv_inference \
--featuremap_path "{annotated_featuremap}" \
--model_joblib_path "{model_file}" \
--output_path "{output_featuremap}" \
--process_number {cpus}
```

*Recommended hardware - 10 CPU, 8GB RAM*

***Note - {cpus} is the number of cpus to use***

##### Output files:
  {output_featuremap}: FeatureMap vcf.gz file with all the annotations and SNVQ values, and a respective .tbi index file


## Detailed explanation of keys output files

### output_featuremap
The FeatureMap is a VCF file that contains a record for each SNV in each read, so that multiple entries per locus are possible, with additional information about the SNV and about the read encoded as INFO fields. Additionally, a machine learning model is trained on these features to assign an SNVQ quality score, saved as the QUAL field of each SNV in the FeatureMap. This value indicates the predicted aggregate error rate, the rate of calling false SNVs in individual reads, when filtering for a given threshold. For example, setting a QUAL60 threshold is expected to yield a 1ppm SNV error rate.

Note that unlike standard vcf files, in a FeatureMap file multiple entries per locus with the same ref and alt are possible and indeed, common. 

### featuremap_df_file
- chrom: object, Chromosome (SNV coordinates)
- pos: int64, Position (SNV coordinates)
- ref: category, Reference base
- alt: category, Alt base
- qual: float64, SNVQ 
- filter: object, PASS for high quality, PreFiltered for SNVs not meeting the assigned pre_filter criteria, LowQual for SNVQ<40
- X_CIGAR: object, CIGAR string of the read, propagated from the input cram file
- X_EDIST: int64, Reference edit distance of the read from the reference, propagated from the input CRAM file
- X_FC1: int64, Edit distance of the read counting only SNVs
- X_FC2: int64, Edit distance of the read counting only SNVs that pass the adjacent base filter
- X_READ_COUNT: int64, Number of reads containing this locus (coverage)
- X_FILTERED_COUNT: int64, Number of reads containing this locus that agree with reference and pass the adjacent base filter
- X_FLAGS: int64, FLAGS propagated from the CRAM file
- X_INDEX: int64, Ordinal index, from start of the read, where the feature was found
- X_LENGTH: int64, Read length after adapter trimming
- X_MAPQ: int64, Mapping quality of the read, propagated from the input cram file
- X_RN: object, The name of the read, propagated from the input cram file
- X_SCORE: float64, Base calling quality, the likelihood of the SNV being a base calling error, evaluated from the read quality, Phred scaled
- X_SMQ_LEFT: int64, Median quality of N bases to the left of the feature
- X_SMQ_LEFT_MEAN: int64, Mean quality of N bases to the left of the feature
- X_SMQ_RIGHT: int64, Median quality of N bases to the right of the feature
- X_SMQ_RIGHT_MEAN: int64, Mean quality of N bases to the right of the feature
- rq: float64, read quality (lower is better), propagated from the input CRAM file
- st: int64, "MIXED" if ppmSeq tag in the start of the read called the read as mixed, named "strand_ratio_category_start" in ppmSeq_legacy_v5
- et: int64, "MIXED" if ppmSeq tag in the end of the read called the read as mixed, named "strand_ratio_category_end" in ppmSeq_legacy_v5
- tm: object, UG trimming tag - "A" indicates 3' adapter was trimmed, Q and/or Z indicate quality trimming, propagated from the CRAM file
- is_forward: bool, interpreted from X_FLAGS
- is_duplicate: bool, interpreted from X_FLAGS
- max_softclip_length: int64, maximal softclip length in either end of the read
- prev_1: category, reference base 1bp before the SNV
- next_1: category, reference base 1bp after the SNV
- prev_2: category, reference base 2bp before the SNV
- next_2: category, reference base 2bp after the SNV
- prev_3: category, reference base 3bp before the SNV
- next_3: category, reference base 3bp after the SNV
- hmer_context_ref: int64, Length of homopolymer the base is contained, in the reference allele
- hmer_context_alt: int64, Length of homopolymer the base is contained, in the alt allele (with the called base)
- is_cycle_skip: bool, is the SNV a cycle skip
- fold_id: int64, The Cross-Validation fold the SNV was contained in, used to match the correct model to use for inference
- label: bool, 0 for False SNVs (unique in locus), 1 for True SNVs (supporting germline variants)
- ML_prob_1_test: float64, Raw probability of the classifier for the SNV to be true, before calibration
- ML_prob_0_test: float64, Raw probability of the classifier for the SNV to be false, before calibration
- ML_qual_1_test: float64, Phred score of the classifier for the SNV to be true, before calibration
- ML_qual_0_test: float64, Phred score of the classifier for the SNV to be false, before calibration
- ML_prob_1_train: float64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training
- ML_prob_0_train: float64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training
- ML_qual_1_train: float64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training
- ML_qual_0_train: float64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training
- ML_prediction_1_train: int64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training
- ML_prediction_0_train: int64, same as above, score for the SNV when it was used in the training set, some values might be missing if SNV was not used for training

### Model joblib file
A dictionary that contains all components needed to run inference on a FeatureMap. Contains 3 keys: 
* "models" - a list of {num_CV_folds} models, one per fold. 
* "params" - a dictionary of all training parameters. Noteworthy parameters include
  - "chroms_to_folds" - a mapping of chromosome values to the corresponding `fold_id` values. Needed to decide which of the k-fold models is used for each SNV. Equalsl None when {split_folds_by}=="random".
  - "categorical_features_dict" - dictionary of all categorical variables and the corresponding category values. Note that the order of the category values is important: the same order should be used during inference as was used during training. 
* "quality_interpolation_function" - a function that maps `ML_qual_1` values (model outputs) to `qual` values. 


## Model training features and filter

Parameters for the srsnv_training command described in TrainSnvQualityRecalibrationModel are given as a json file referred to as {single_read_snv_params} with the values below. Notice the comments inline.
```json
{
  "ppmSeq_adapter_version": "v1",
  "SingleReadSNV.categorical_features": {
    "st": ["MIXED", "MINUS", "PLUS", "END_UNREACHED", "UNDETERMINED"],  # (for ppmSeq_legacy_v5 adapters change key to "strand_ratio_category_start", for non-ppmSeq remove)
    "et": ["MIXED", "MINUS", "PLUS", "END_UNREACHED", "UNDETERMINED"],  # (for ppmSeq_legacy_v5 adapters change key to "strand_ratio_category_end", for non-ppmSeq remove)
    "ref": ["A", "C", "G", "T"],
    "alt": ["A", "C", "G", "T"],
    "next_1": ["A", "C", "G", "T"],
    "next_2": ["A", "C", "G", "T"],
    "next_3": ["A", "C", "G", "T"],
    "prev_1": ["A", "C", "G", "T"],
    "prev_2": ["A", "C", "G", "T"],
    "prev_3": ["A", "C", "G", "T"]
  },
  "numerical_features": [
    "X_SCORE",
    "X_EDIST",
    "X_LENGTH",
    "X_INDEX",
    "X_FC1",
    "rq",
    "max_softclip_length",
    "hmer_context_ref",
    "hmer_context_alt"
  ],
  "boolean_features": [
    "is_cycle_skip",
    "is_forward"
  ],
  "balanced_sampling_info_fields": [
    "trinuc_context_with_alt",
    "is_forward"
  ],
  "pre_filter": "(X_SCORE>4) && (X_EDIST<10)",
  "random_seed": 0,
  "num_CV_folds": 5,
  "split_folds_by": "chrom",
  "train_set_size": 3000000,
  "test_set_size": 0
}
```

Explanation of key features:
* {pre_filter} - cutoff criteria for SNVs to be included in the model
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
  - "st" or "strand_ratio_category_start" (ppmSeq read category measured in the read start)
  - "et" or "strand_ratio_category_end" (ppmSeq read category measured in the read end)

  ***Note - the last two are only required for ppmSeq data***
* {num_CV_folds} - Number of Cross-Validation folds to use. 
* {split_folds_by} - when using k-fold Cross-Validation, method by which SNVs are assigned `fold_id` values. Two methods are possible: "chrom" and "random". 
  
  ***Note - See section [Details of ML model Cross-Validation scheme](#details-of-ml-model-cross-validation-scheme) below for more details on Cross-Validation scheme and parameter values.***


## Details of ML model Cross-Validation scheme
The ML model can be trained either by employing a test/train split, or by employing a k-fold Cross-Validation (CV) scheme. 
* Train/test split is appropriate when inference is done on SNVs in regions that are excluded from the ML model training set. 
* CV should be used for use-cases that require `qual` values for SNVs in regions that are included in the training set. A common example is when `qual` values are required for all SNVs in the FeatureMap. 

This section describes how `fold_id` values are assigned in each case, and how `fold_id` values are treated during inference. 

### Training with train/test split
When {num_CV_folds}==1, train/test split is employed. To avoid data leakage from training to test set, SNVs are split into a train and test set by position: {train_set_size}+{test_set_size} SNVs are chosen at random and sorted by position. Then the first {train_set_size} SNVs are assigned to the train set, and the remaining {test_set_size} SNVs are assigned to the test set. 

In the featuremap_df_file, SNVs belonging to the training set are assigned a `fold_id` value of -1, and SNVs belonging to the test set are assigned a `fold_id` value of 0. A single model is trained on the training set. During inference, this model is used for all SNVs. 

### Training with k-fold CV
When {num_CV_folds} >= 2, CV is employed. Data is split into {num_CV_folds} folds in one of two methods: by chromosome (preferred method), or randomly. 
* When {split_folds_by}=="chrom", SNVs are split into folds by chromosome. This is the preferred method, as it makes sure that there is no data leakage from training to validation/inference. Chromosomes are divided into {num_CV_folds} groups of approximately equal size (measured in bases). All chromosomes of the same group are assinged the same `fold_id` value in the featuremap_df_file. 

  Chromosomes X, Y, M are excluded from training. SNVs in these chromosomes are assigned a `fold_id` value of `nan`. 
* When {split_folds_by}=="random", all SNVs are assigned random `fold_id` values in the featuremap_df_file (between 0 and {num_CV_folds}-1). 

{num_CV_folds} models are trained, one per fold: model `i` is trained on all SNVs whose `fold_id` values are numerical (not `nan`) and not equal to `i`. Then, train and test model predictions are evaluated:
* Test values in the featuremap_df_file (e.g., ML_prob_0_test, ML_prob_1_test, ML_qual_1_test, etc.) are obtained by evaluating model `i` on all SNVs whose `fold_id` value is `i`, and evaluating all models on SNVs whose `fold_id` value is `nan` (model predictions are then averaged). 
* Train values in the featuremap_df_file (e.g., ML_prob_0_train, ML_prob_1_train, ML_qual_1_train, etc.) are obtained by evaluating all models except model `i` on all SNVs whose `fold_id` value is `i` (model predictions are then averaged). SNVs whose `fold_id` value is `nan` are not assigned train values. 

During inference, each SNV is evaluated using the model(s) of the appropriate fold_id. 
* When {split_folds_by}=="chrom", an SNV is assigned a `fold_id` value by its chromosome, and the corresponding model is used. Chromosomes excluded from training (`fold_id` is `nan`) are evaluated by averaging the predictions of all {num_CV_folds} models. 
* When {split_folds_by}=="random", all SNVs are evaluated by averaging the predictions of all models.

---

## References
Cheng, Alexandre Pellan, et al. "Whole genome error-corrected sequencing for sensitive circulating tumor DNA cancer monitoring." bioRxiv (2022): 2022-11.