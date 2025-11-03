<div style="font-size: 250%; font-weight: bold; text-align: center;">
Single Read SNV (SRSNV) pipeline (v1.23.0)
</div>

<div style="font-size: 150%; font-weight: bold; text-align: left;">
Table of Contents
</div>

- [Introduction](#introduction)
  - [SNV denoising and quality recalibration](#snv-denoising-and-quality-recalibration)
  - [Training set preparation](#training-set-preparation)
  - [Data filtering](#data-filtering)
    - [Genomic region filters](#genomic-region-filters)
    - [Alignment filters](#alignment-filters)
    - [Quality filters](#quality-filters)
  - [List of FeatureMap fields](#list-of-featuremap-fields)
- [Variables (WDL inputs/outputs)](#variables-wdl-inputsoutputs)
  - [Primary inputs](#primary-inputs)
  - [Model / training control inputs](#model--training-control-inputs)
  - [Key outputs](#key-outputs)
- [Pipeline overview](#pipeline-overview)
  - [FeatureMap generation with snvfind](#featuremap-generation-with-snvfind)
  - [Coverage gate and model decision](#coverage-gate-and-model-decision)
  - [Preparing training datasets](#preparing-training-datasets)
  - [Model training (srsnv\_training)](#model-training-srsnv_training)
  - [Inference (snvqual)](#inference-snvqual)
  - [Reporting (srsnv\_report)](#reporting-srsnv_report)
- [Template descriptions](#template-descriptions)
- [Manual execution (outside WDL)](#manual-execution-outside-wdl)
  - [Environment / dockers](#environment--dockers)
  - [Step-by-step commands](#step-by-step-commands)
- [FeatureMap content \& quality scores](#featuremap-content--quality-scores)
- [Training feature set \& filters](#training-feature-set--filters)
  - [single\_read\_snv\_params JSON example (updated)](#single_read_snv_params-json-example-updated)
  - [Explanation of important parameters](#explanation-of-important-parameters)
- [Cross-validation scheme](#cross-validation-scheme)
  - [Train/test split (num\_CV\_folds=1)](#traintest-split-num_cv_folds1)
  - [K-fold CV (num\_CV\_folds\>=2)](#k-fold-cv-num_cv_folds2)
- [Model files](#model-files)
- [When qualities are NOT assigned](#when-qualities-are-not-assigned)

## Introduction

### SNV denoising and quality recalibration
The Single Read SNV (SRSNV) pipeline is a read-centric de-noising framework, developed to overcome the limitations of traditional locus-centric variant calling, particularly in scenarios where rare mutations may be supported by only a single read. These rare mutations need to be distinguished from artefactual SNVs, which can derive from sequencing, library or alignment errors. To achieve this, we employed a supervised machine learning model trained to classify actual SNVs (labelled True or TP) from noise (False or FP). First, a comprehensive dataset capturing every candidate SNV is generated, along with a rich suite of annotations that describe sequencing quality, local sequence motifs, fragment-specific features, and locus-specific information. Randomly selected bases in the data matching the reference genome are collected and annotated as True, while low VAF (≤5%) SNVs in high-coverage (≥20×) regions (SNVs with low support, indicating they are likely to be artifacts) are annotated as False SNVs. Using these curated sets, we train an XGBoost classifier to robustly distinguish between true and artifactual SNVs. Once trained, the classifier assigns a calibrated quality score to each SNV in the input CRAM, providing a precise estimate of the residual error rate. To avoid overfitting, an ensemble of models (3) are trained on different sets of chromosomes and applied using a cross-validation scheme. 

We now provide detailed information on the denoising procedure outlined below.

### Training set preparation
First, all the SNVs in the input data (CRAM file) are exported to a vcf file denoted a FeatureMap (f1), where multiple annotation for every read supporting each SNV are collected. Additionally, a random sample of bases is collected and exported in the same format to a separate file (f2). Then, a series of filters is applied to each file, as detailed below, leaving in only high quality data. Finally, the randomly sample bases in f2 are filtered for bases matching the reference genome, which and annotated as True, and the SNVs in f1 are filtered for low VAF (≤5%) and high-coverage (≥20×) and annotated as False. 

The randomly sampled data can be sampled so that it follows a provided distribution of trinucleotide+alt motifs (see https://cancer.sanger.ac.uk/signatures/sbs/ for example of this kind of distribution). By default, the trinucleotide distribution of the reference genome inside the UG HCR is provided as input, with a uniform distribution across the 3 alts (e.g. ACG->AAG,ACG->AGG,ACG->ATG all have the same probability). 

### Data filtering
The filtering cascade described above incorporates both **hard** and **soft** filters. **Hard filters** remove SNVs entirely from the output VCF file. **Soft filters**, by contrast, retain SNVs in the VCF and assign them a quality score, but annotate them as filtered (FORMAT/FILT=0). Downstream analyses should therefore treat soft-filtered SNVs with caution. The effects of all filters on the randomly sampled bases are summarized in an output metadata file to support quality control (QC) and normalization. The default filters are described below.

#### Genomic region filters
Training data was filtered for SNVs included in the Ultima Genomics High Confidence Region in chromosomes 1-22. The UG high-confidence region (HCR) covers >99% of the GIAB v4.2.1 HCR and excludes genomic areas where UG performance is consistently of lower confidence, such as regions of low complexity (see [ug_hcr.md](https://github.com/Ultimagen/healthomics-workflows/blob/main/docs/ug_hcr.md)). Homopolymer regions of length >7 bp are also excluded. Additionally, SNVs represented in dbsnp or in gnomAD with AF>0.001 were excluded to avoid risk of germline contamination introducing labeling noise.

#### Alignment filters
Only reads with a mapping quality of 60 (maximal value in the aligner) are used in the model training set. From these reads, all SNVs with respect to the reference genome, where the adjacent bases (5) from the respective SNV match the reference genome, are extracted to a separate VCF file. The adjacent base filter is intended to avoid calling compound homopolymer errors as SNVs. For example, consider the following read: \
`REF : ACCGT (1A 2C 1G 1T)` \
`READ: ACGTT (1A 1C 1G 2T)` \
The difference from the reference can be interpreted as either two homopolymer (hmer) errors (2->1C & 1->2T) or as two SNVs (C->G & G->T). While the former is significantly more prevalent in Ultima Genomics data, the latter is generally preferred by most aligners, including the Ultima Aligner. To avoid these alignment errors being interpreted as true SNVs, the adjacent base filter is applied. Additionally, this filter is applied as a hard filter by default, leaving those mismatches out of the FeatureMap vcf, as they otherwise account for the majority of entries and bloat the file size.

#### Quality filters
Reads with a low base calling scores (BCSQ<40) are soft-filtered and not included in the training set. Additionally, reads with high edit distances (Levenshtein) with respect to the reference genome (EDIST≥10) are excluded. Lastly, entries where the alt allele create a long homopolymer (>7) are excluded. For example, this SNV would be soft-filtered: \
`REF : CTTTTGTTTTA` \
`READ: CTTTTTTTTTA` \
Because the ALT allele is contained in a 9T homopolymer. Note that the equivalent filter for reference homopolymers is applied as a genomic region filter.


### List of FeatureMap fields
Multiple values are reported in the FeatureMap output, comprised of the following types:
1. Variant features in INFO. These describe the variant itself regardless of sample and specific reads. Example - the reference base preceding the variant.
2. Sample features in FORMAT, Number=1/A (see below). These describe features that are specific to a given sample, but do not describe specific reads. Example - read depth.
3. Read features in FORMAT, Number=. (see below). These describe features that are specific to every read supporting the alt, given as lists that are all the same size within a reported variant. Example - read length. Can contain values copied from CRAM tags.

Below is a detailed table describing the various fields, as also described in the header of every output FeatureMap vcf.

| Name        | Legacy (<V1.23.0) equivalent  | Where   | Number*| Type    | Description                                                                                            |
|-------------|-------------------------------|---------|--------|---------|--------------------------------------------------------------------------------------------------------|
| X_PREV1     | prev_1                        | INFO    | 1      | String  | Reference base at position POS-1 {A,C,G,T}**                                                           |
| X_NEXT1     | next_1                        | INFO    | 1      | String  | Reference base at position POS+1 {A,C,G,T}**                                                           |
| X_PREV2     | prev_2                        | INFO    | 1      | String  | Reference base at position POS-2 {A,C,G,T}**                                                           |
| X_NEXT2     | next_2                        | INFO    | 1      | String  | Reference base at position POS+2 {A,C,G,T}**                                                           |
| X_PREV3     | prev_3                        | INFO    | 1      | String  | Reference base at position POS-3 {A,C,G,T}**                                                           |
| X_NEXT3     | next_3                        | INFO    | 1      | String  | Reference base at position POS+3 {A,C,G,T}**                                                           |
| X_HMER_REF  | hmer_context_ref              | INFO    | 1      | Integer | Homopolymer context in the ref allele (up to length 20)                                                |
| X_HMER_ALT  | hmer_context_alt              | INFO    | 1      | Integer | Homopolymer context in the alt allele (up to length 20)                                                |
| DP          | X_READ_COUNT                  | FORMAT  | 1      | Integer | Number of reads containing this location                                                               |
| DP_FILT     | X_FILTERED_COUNT              | FORMAT  | 1      | Integer | Number of reads containing this location that pass the adjacent base filter                            |
| DP_MAPQ60   |                               | FORMAT  | 1      | Integer | Number of reads with mapping quality ≥ 60                                                              |
| RAW_VAF     |                               | FORMAT  | 1      | Float   | Raw VAF := N_alt_reads/N_total_reads                                                                   |
| VAF         |                               | FORMAT  | 1      | Float   | VAF := N_alt_reads/(N_ref_reads+N_alt_reads)                                                           |
| AD          |                               | FORMAT  | A      | Integer | Number of reads supporting the reference allele in locus                                               |
| AD_A        |                               | FORMAT  | 1      | Integer | Number of reads supporting the base A in locus                                                         |
| AD_C        |                               | FORMAT  | 1      | Integer | Number of reads supporting the base C in locus                                                         |
| AD_G        |                               | FORMAT  | 1      | Integer | Number of reads supporting the base G in locus                                                         |
| AD_T        |                               | FORMAT  | 1      | Integer | Number of reads supporting the base T in locus                                                         |
| AD_DEL      |                               | FORMAT  | 1      | Integer | Number of reads supporting a deletion in locus                                                         |
| AD_INS      |                               | FORMAT  | 1      | Integer | Number of reads supporting an adjacent insertion in locus                                              |
| BCSQ        | X_SCORE                       | FORMAT  | .      | Integer | Base calling error likelihood in calling the SNV, in Phred scale                                       |
| BCSQCSS     |                               | FORMAT  | .      | Integer | Cycle Skip Size when computing base calling error likelihood                                           |
| RL          | X_LENGTH                      | FORMAT  | .      | Integer | Read length (post adapter trimming)                                                                    |
| INDEX       | X_INDEX                       | FORMAT  | .      | Integer | Position in the read of the SNV relative to read start                                                 |
| RN          | X_RN                          | FORMAT  | .      | String  | Query (read) name                                                                                      |
| DUP         | is_duplicate                  | FORMAT  | .      | Integer | Is the read a duplicate, interpreted from CRAM flag                                                    |
| REV         | is_forward                    | FORMAT  | .      | Integer | Is the read mapped to the reverse strand, interpreted from CRAM flag                                   |
| SCST        | ~max_softclip_length          | FORMAT  | .      | Integer | Softclip length in the start of the read (synthesis direction)                                         |
| SCED        | ~max_softclip_length          | FORMAT  | .      | Integer | Softclip length in the end of the read (synthesis direction)                                           |
| MAPQ        | X_MAPQ                        | FORMAT  | .      | Integer | Read mapping quality                                                                                   |
| EDIST       | X_EDIST                       | FORMAT  | .      | Integer | Read Levenshtein edit distance from the reference                                                      |
| HAMDIST     | X_FC1                         | FORMAT  | .      | Integer | Hamming distance (SNVs only, disregarding indels) from the reference                                   |
| HAMDIST_FILT| X_FC2                         | FORMAT  | .      | Integer | Filtered Hamming distance (SNVs passing adjacent base filter)                                          |
| SMQ_BEFORE  | X_SMQ_LEFT_MEAN               | FORMAT  | .      | Integer | Mean quality of 20 bases before the locus                                                              |
| SMQ_AFTER   | X_SMQ_RIGHT_MEAN              | FORMAT  | .      | Integer | Mean quality of 20 bases after the locus                                                               |
| ADJ_REF_DIFF|                               | FORMAT  | .      | Integer | The 5 adjacent bases to the locus do not fully match reference genome (if 0 the SNV passes)            |
| tm          | tm                            | FORMAT  | .      | String  | Trimming reasons {A,AQ,AQZ,AZ,Q,QZ,Z}                                                                  |
| a3          | a3                            | FORMAT  | .      | Integer | Start position in input of segment "A-tailing and native adapter"                                      |
| rq          | rq                            | FORMAT  | .      | Float   | Read quality (copied from CRAM tag)                                                                    |
| st          | st                            | FORMAT  | .      | String  | Name of pattern matched in segment "Start_loop" {MIXED,MINUS,PLUS,UNDETERMINED} (copied from CRAM tag) |
| et          | et                            | FORMAT  | .      | String  | Name of pattern matched in segment "End_loop" {MIXED,MINUS,PLUS,UNDETERMINED} (copied from CRAM tag)   |
| MI          |                               | FORMAT  | .      | String  | MI (copied from CRAM tag)                                                                              |
| DS          |                               | FORMAT  | .      | Integer | DS (copied from CRAM tag)                                                                              |
| FILT        |                               | FORMAT  | .      | Integer | Pre-filter status for SNV reads (1=pass, 0=fail)                                                       |
| FILT_BITMAP |                               | FORMAT  | .      | String  | Filter bitmaps                                                                                         |
| MQUAL       | ML_QUAL                       | FORMAT  | .      | Float   | SingleReadSNV model inferred raw Phred scaled quality                                                  |
| SNVQ        | QUAL                          | FORMAT  | .      | Float   | SingleReadSNV model inferred Phred scaled quality, recalibrated to SNVQ                                |
| sd          |                               | FORMAT  | .      | Integer  | Matching distance of "Start_loop" ***                                |
| ed          |                               | FORMAT  | .      | Integer  | Matching distance of "End_loop" ***                                |
| l1          |                               | FORMAT  | .      | Integer  | Length of insertion before "Start_loop" ***                                |
| l2          |                               | FORMAT  | .      | Integer  | Length "Start_loop" ***                                |
| l3          |                               | FORMAT  | .      | Integer  | Length "Stem_start" ***                                |
| l4          |                               | FORMAT  | .      | Integer  | Length of insert ***                                |
| l5          |                               | FORMAT  | .      | Integer  | Length of "Stem_end" ***                                |
| l6          |                               | FORMAT  | .      | Integer  | Length of "End_loop" ***                                |
| l7          |                               | FORMAT  | .      | Integer  | Length of adapter following "End_loop" ***                                |
| q2          |                               | FORMAT  | .      | Integer  | Average quality of "Start_loop" ***                                |
| q3          |                               | FORMAT  | .      | Integer  | Average quality of "Stem_start" ***                                |
| q4          |                               | FORMAT  | .      | Integer  | Average quality of insert ***                                |
| q5          |                               | FORMAT  | .      | Integer  | Average quality of "Stem_end" ***                                |
| q6          |                               | FORMAT  | .      | Integer  | Average quality of "End_loop" ***                                |

  `* Number according to the VCF format, 1="single value", .="multiple values", A="One value per ALT allele". In the FeatureMap output, by convention all the "." values are lists of the same length corresponding to the read supporting the reported SNV.  
  
  ** For categorical values, the list of allowed values is given in curly brackets.

  *** This tag is generated by the Trimmer software running on the UG tool (length=length of segment, distance=distance of matched sequence from expected sequence, quality=average quality of segment), only available in CRAM files produced by UG SWPKG1.9+. In older data this tag does not exist and will result in a missing value in all the entries, and will not affect model performance.

## Variables (WDL inputs/outputs)
Naming convention - Curly braces denote variable names as passed to / produced by the WDL (e.g. `{featuremap}`, `{model_files}`, `{srsnv_metadata_json}`).

### Primary inputs
- `{input_cram_bam}`, `{input_cram_bam_index}`
- `{sorter_json_stats_file}` (provides coverage & total aligned bases)
- `{base_file_name}`
- `{references.ref_fasta}`, `{references.ref_fasta_index}`, `{references.ref_dict}`
- `{featuremap_params}` (controls `snvfind` emission thresholds / context)
- `{training_regions_interval_list}` (+ its index)
- `{ref_trinuc_freq}` trinucleotide distribution to sample the random bases by
- `{min_coverage_to_train_model}`
- Optional: `{random_sample_trinuc_freq}` (trinucleotide frequency distribution for random sampling)

### Model / training control inputs
- `{single_read_snv_params}` (sampling sizes, filters, CV config)
- `{features}` (ordered list of feature names used in training)
- `{xgboost_params_file}` (model hyperparameters)
- Optional: `{pre_trained_model_files}`, `{pre_trained_srsnv_metadata_json}`

### Key outputs
- Output FeatureMap: `{featuremap}` + `{featuremap_index}`
- Random sample FeatureMap: `{featuremap_random_sample}`, `{featuremap_random_sample_index}`, `{downsampling_rate}`
- Optional: `{random_sample_trinuc_freq_stats}` (actual trinucleotide frequencies in random sample, if trinucleotide-aware sampling was used)
- Training datasets (if self-trained):
  - `{raw_filtered_featuremap_parquet}`, `{raw_featuremap_stats}`
  - `{random_sample_filtered_featuremap_parquet}`, `{random_sample_featuremap_stats}`
  - `{featuremap_df}` (combined labeled dataset with folds & predictions)
- Model files: `{model_files}` (array), `{srsnv_metadata_json}`
- Report & QC: `{report_html}`, `{application_qc_h5}`
- Flags: `{snv_qualities_assigned}`, `{used_self_trained_model}`

## Pipeline overview

### FeatureMap generation with snvfind
`snvfind` emits:
- Full raw FeatureMap: `{base}.raw.featuremap.vcf.gz`
- Randomly sampled bases FeatureMap: `{base}.random_sample.featuremap.vcf.gz`
A precise `downsampling_rate` is computed as:
```
random_sample_size / total_aligned_bases
```
(derived from sorter stats). This enables consistent positive sampling sizes independent of coverage.

**Trinucleotide-aware sampling**: If `{random_sample_trinuc_freq}` is provided (CSV file with trinucleotide frequencies), the random sample will be drawn according to the specified trinucleotide distribution rather than uniformly. By default, the trinucleotide distribution of the reference genome is used, to create a training set that is representative of the WGS data. The actual trinucleotide frequencies achieved in the random sample are saved to `{random_sample_trinuc_freq_stats}`.

### Coverage gate and model decision
If mean coverage `< {min_coverage_to_train_model}`:
- If a pre-trained model (metadata + models) is provided → use it.
- Else → emit raw FeatureMap without SNVQ (flag `{snv_qualities_assigned}=false`).

### Preparing training datasets
`PrepareFeatureMapForTraining`:
- Restricts to `{training_regions_interval_list}`
- Converts VCF → parquet (`featuremap_to_dataframe`)
- Applies sequential filters:
  - Coverage floor: `DP >= min_coverage_filter`
  - Coverage ceiling: `DP <= mean_coverage * max_coverage_factor`
  - Optional `pre_filters` (from `single_read_snv_params`)
  - Additional label filters (e.g. VAF, REF≠ALT) depending on positive/negative set construction
- Downsamples to target sizes (`tp_train_set_size`, `fp_train_set_size`)
- Produces stats JSON per filtered set

Labeling strategy:
- **Positive dataset**: Random sample filtered for REF == ALT (bases matching reference)
- **Negative datasets**: Raw featuremap filtered for low VAF (≤5%)
- Stats for both sets are passed to training for auditability

### Model training (srsnv_training)
Command consumes:
- `--positive` parquet file (random sample, REF==ALT)
- `--negative` parquet file (raw featuremap, low VAF)
- `--stats-positive`: statistics from positive filtering
- `--stats-negative`: statistics from random sample negative filtering (REF≠ALT, low VAF)
- `--stats-featuremap`: statistics from raw featuremap filtering
- Training regions (for metadata)
- Feature list (ordering preserved)
- K-fold CV (`num_CV_folds`) with chromosome-based fold assignment (preferred)
Outputs:
- `*.model_fold_{i}.json`
- `{base}.srsnv_metadata.json` (fold mapping, features, params, coverage, filters)
- `{base}.featuremap_df.parquet` (includes fold_id, train/test predictions)

### Inference (snvqual)
`snvqual` loads `srsnv_metadata.json` plus model fold JSON files; assigns SNVQ into QUAL for each record:
- Chromosome → model fold (if fold-defined)
- Chromosomes excluded from training → average predictions across folds
Outputs final `{featuremap}`.

### Reporting (srsnv_report)
Generates:
- `{base}.report.html`
- `{base}.single_read_snv.applicationQC.h5`

## Template descriptions

The following input templates are available for different kinds of input data:

| Template File | Description |
|---------------|-------------|
| `single_read_snv_template-ppmSeq.json` | Use this template for ppmSeq data. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. st, et). |
| `single_read_snv_template-ppmSeq_legacy_v5.json` | Use this template for LEGACY v5 ppmSeq data. This is an older version of the ppmSeq adapters, generally not available since 2024. The input CRAM file should be trimmed, aligned and sorted, and contain the ppmSeq tags (e.g. as, ts). |
| `single_read_snv_template-Standard-WG.json` | Use this template for any non-ppmSeq data. |

## Manual execution (outside WDL)

### Environment / dockers
Pull (matching versions referenced in workflows/single_read_snv/tasks/globals.wdl in this repository):
- `featuremap_docker` (snvfind + snvqual)
- `ugbio_featuremap_docker` (featuremap_to_dataframe, filter_featuremap)
- `ugbio_srsnv_docker` (srsnv_training, srsnv_report)


### Step-by-step commands
Variable names mirror WDL inputs; example corresponds to ppmSeq template. See [wdls/input_templates/single_read_snv_template-ppmSeq.json](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_read_snv/input_templates/single_read_snv_template-ppmSeq.json) for full filenames, only base names are used below for clarity.
```bash
# Set base name & inputs
BASE=sample
CRAM=input.cram
CRAM_INDEX=input.cram.crai
SORTER_STATS=sorter_stats.json
REF=Homo_sapiens_assembly38.fasta
TRAINING_REGIONS=ug_rare_variant_hcr.Homo_sapiens_assembly38.interval_list.gz
TRAINING_REGIONS_INDEX=${TRAINING_REGIONS}.tbi
XGBOOST_PARAMS=xgboost_model_params.json
BED=wgs_calling_regions.without_encode_blacklist.hg38.bed
TRINUC_FREQ=ref_trinuc_freq.csv  # Optional: reference genome trinucleotide frequency file for random sampling

# single_read_snv_params (ppmSeq template)
TP_TRAIN_SET_SIZE=1500000
FP_TRAIN_SET_SIZE=1500000
TP_OVERHEAD=10.0              # tp_train_set_size_sampling_overhead
MAX_VAF_FOR_FP=0.05
MIN_COV_FILTER=20             # min_coverage_filter
MAX_COV_FACTOR=2.0            # max_coverage_factor
RANDOM_SEED=0
NUM_FOLDS=3                   # num_CV_folds

# Derived random sample size (ceil(tp_train_set_size * overhead))
RANDOM_SAMPLE_SIZE=$(( TP_TRAIN_SET_SIZE * 10 ))  # 15000000

# 1. Compute downsampling rate from sorter stats
TOTAL_ALIGNED_BASES=$(jq -re '.total_aligned_bases // .total_bases // error("missing total_aligned_bases")' "$SORTER_STATS")
DOWNSAMPLING_RATE=$(awk -v num=$RANDOM_SAMPLE_SIZE -v den=$TOTAL_ALIGNED_BASES 'BEGIN{printf "%.12f", num/den}')
echo "Downsampling rate: $DOWNSAMPLING_RATE"

# Mean coverage extraction
# run inside ugbio_srsnv docker
MEAN_COVERAGE_FILE=${BASE}.mean_coverage.txt
sorter_stats_to_mean_coverage \
  --sorter-stats-json "$SORTER_STATS" \
  --output-file "$MEAN_COVERAGE_FILE"

MEAN_COVERAGE=$(cat "$MEAN_COVERAGE_FILE")
echo "Mean coverage: $MEAN_COVERAGE"
COVERAGE_CEIL=$(printf "%.0f" "$(echo "$MEAN_COVERAGE * $MAX_COV_FACTOR" | bc -l)")
echo "Coverage ceiling: $COVERAGE_CEIL"

# 2. snvfind (raw + random sample)
# FeatureMap params (ppmSeq): min_mapq=60 padding=5 score_limit=100 exclude_nan_scores=true include_dup_reads=true
# surrounding_quality_size=20 reference_context_size=3 keep_supplementary=false
# cram_tags_to_copy joined by commas below
CRAM_TAGS="tm:Z:A:AQ:AQZ:AZ:Q:QZ:Z,a3:i,rq:f,st:Z:MIXED:MINUS:PLUS:UNDETERMINED,et:Z:MIXED:MINUS:PLUS:UNDETERMINED,MI:Z,DS:i,sd:i,ed:i,l1:i,l2:i,l3:i,l4:i,l5:i,l6:i,l7:i,q2:i,q3:i,q4:i,q5:i,q6:i"

snvfind "$CRAM" "$REF" \
  -o ${BASE}.raw.featuremap.vcf.gz \
  -f ${BASE}.random_sample.featuremap.vcf.gz,${DOWNSAMPLING_RATE}${TRINUC_FREQ} \
  -v \
  -p 5 -L 100 -n -d -Q 20 -r 3 -m 60 -c "$CRAM_TAGS" -b "$BED"

bcftools index -t ${BASE}.raw.featuremap.vcf.gz
bcftools index -t ${BASE}.random_sample.featuremap.vcf.gz



# 3. Prepare RAW (negative / FP labeling set from raw featuremap)
#   a) Restrict to training regions
bcftools view ${BASE}.raw.featuremap.vcf.gz -T "$TRAINING_REGIONS" -Oz -o ${BASE}.raw.training_regions.vcf.gz
bcftools index -t ${BASE}.raw.training_regions.vcf.gz

#   b) Convert to parquet
featuremap_to_dataframe \
  --input ${BASE}.raw.training_regions.vcf.gz \
  --output ${BASE}.raw.training_regions.parquet \
  --drop-format GT AD X_TCM

#   c) Filter + label (RAW_VAF <= MAX_VAF_FOR_FP) + downsample to FP_TRAIN_SET_SIZE
filter_featuremap \
  --in  ${BASE}.raw.training_regions.parquet \
  --out ${BASE}.raw.filtered.parquet \
  --stats ${BASE}.raw.stats.json \
  --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
  --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
  --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
  --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
  --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
  --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
  --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
  --filter name=low_vaf:field=RAW_VAF:op=le:value=${MAX_VAF_FOR_FP}:type=label \
  --downsample random:${FP_TRAIN_SET_SIZE}:${RANDOM_SEED}

# 4. Prepare RANDOM SAMPLE (positive / TP labeling set)
#   a) Restrict to training regions
bcftools view ${BASE}.random_sample.featuremap.vcf.gz -T "$TRAINING_REGIONS" -Oz -o ${BASE}.rs.training_regions.vcf.gz
bcftools index -t ${BASE}.rs.training_regions.vcf.gz

#   b) Convert to parquet
featuremap_to_dataframe \
  --input ${BASE}.rs.training_regions.vcf.gz \
  --output ${BASE}.rs.training_regions.parquet \
  --drop-format GT AD X_TCM

#   c) Filter + label (REF == ALT) + downsample to TP_TRAIN_SET_SIZE
filter_featuremap \
  --in  ${BASE}.rs.training_regions.parquet \
  --out ${BASE}.rs.filtered.parquet \
  --stats ${BASE}.rs.stats.json \
  --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
  --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
  --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
  --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
  --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
  --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
  --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
  --filter name=ref_eq_alt:field=REF:op=eq:value_field=ALT:type=label \
  --downsample random:${TP_TRAIN_SET_SIZE}:${RANDOM_SEED}

#   d) Filter + negative label (RAW_VAF <= MAX_VAF_FOR_FP) + downsample to TP_TRAIN_SET_SIZE
filter_featuremap \
  --in  ${BASE}.rs.training_regions.parquet \
  --out ${BASE}.rs_neg.filtered.parquet \
  --stats ${BASE}.rs_neg.stats.json \
  --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
  --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
  --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
  --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
  --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
  --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
  --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
  --filter name=ref_ne_alt:field=REF:op=ne:value_field=ALT:type=label \
  --filter name=low_vaf:field=RAW_VAF:op=le:value=${MAX_VAF_FOR_FP}:type=label \
  --downsample random:${TP_TRAIN_SET_SIZE}:${RANDOM_SEED}

# 5. Train (NUM_FOLDS=3 per ppmSeq template)
FEATURES="REF:ALT:X_PREV1:X_NEXT1:X_PREV2:X_NEXT2:X_PREV3:X_NEXT3:X_HMER_REF:X_HMER_ALT:BCSQ:BCSQCSS:RL:INDEX:REV:SCST:SCED:SMQ_BEFORE:SMQ_AFTER:tm:rq:st:et:EDIST:HAMDIST:HAMDIST_FILT:l1,l2,l3,l4,l5,l6,l7,q2,q3,q4,q5,q6"

srsnv_training \
  --positive ${BASE}.rs.filtered.parquet \
  --negative ${BASE}.raw.filtered.parquet \
  --stats-positive ${BASE}.rs.stats.json \
  --stats-negative ${BASE}.rs_neg.stats.json \
  --stats-featuremap ${BASE}.raw.stats.json \
  --mean-coverage ${MEAN_COVERAGE} \
  --training-regions $TRAINING_REGIONS \
  --k-folds ${NUM_FOLDS} \
  --model-params $XGBOOST_PARAMS \
  --features $FEATURES \
  --basename $BASE \
  --output . \
  --random-seed ${RANDOM_SEED} \
  --verbose

# 6. Inference
mkdir -p model_files
cp ${BASE}.model_fold_*.json model_files/
cp ${BASE}.srsnv_metadata.json model_files/srsnv_metadata.json

snvqual ${BASE}.raw.featuremap.vcf.gz ${BASE}.featuremap.vcf.gz model_files/srsnv_metadata.json -v
bcftools index -t ${BASE}.featuremap.vcf.gz

# 7. Report
srsnv_report \
  --featuremap-df ${BASE}.featuremap_df.parquet \
  --srsnv-metadata model_files/srsnv_metadata.json \
  --report-path . \
  --basename ${BASE} \
  --verbose
```

## FeatureMap content & quality scores
Each record = one read-level substitution. QUAL = calibrated SNVQ (expected error rate phred-scaled). Multiple identical REF/ALT entries at same locus are normal (different reads).

## Training feature set & filters
Features are supplied as ordered list (`features` input). Categorical ordering must match training-time ordering.

### single_read_snv_params JSON example (updated)
```json
{
  "tp_train_set_size": 3000000,
  "tp_train_set_size_sampling_overhead": 1.05,
  "fp_train_set_size": 3000000,
  "max_vaf_for_fp": 0.05,
  "min_coverage_filter": 20,
  "max_coverage_factor": 10.0,
  "pre_filters": [
    "name=vaf_defined:field=RAW_VAF:op=ge:value=0:type=region"
  ],
  "random_seed": 0,
  "num_CV_folds": 5
}
```

### Explanation of important parameters
- `tp_train_set_size`, `fp_train_set_size`: target counts after filtering & downsampling
- `tp_train_set_size_sampling_overhead`: inflates random sample rate so enough positives survive filtering
- `max_vaf_for_fp`: VAF threshold used in labeling filters (legacy FP concept folded into new labeling logic)
- `max_coverage_factor`: coverage ceiling = mean_coverage * factor
- `min_coverage_filter`: floor for reliable labeling
- `pre_filters`: extra reusable filters (same syntax as used in WDL)
- `num_CV_folds`: ≥2 enables chromosome-based CV (preferred)

## Cross-validation scheme
(Logic unchanged from legacy; terminology updated.)

### Train/test split (num_CV_folds=1)
All SNVs use single model; `fold_id = -1` (train) or `0` (test) inside `featuremap_df`.

### K-fold CV (num_CV_folds>=2)
Chromosomes partitioned into ~balanced groups → `fold_id` per group.
Excluded chroms (e.g. X/Y/MT if configured) get `nan` and are inferred by averaging all folds.

## Model files
- `*.model_fold_X.json`: one JSON per fold (XGBoost serialized parameters)
- `*.srsnv_metadata.json`:
  - feature ordering
  - fold → chromosome map
  - filtering & sampling metadata
  - coverage statistics
- `*.featuremap_df.parquet`: includes:
  - raw & calibrated probabilities
  - fold assignments
  - per-fold train/test predictions

## When qualities are NOT assigned
`snv_qualities_assigned=false` if:
- Coverage < `{min_coverage_to_train_model}` AND no pre-trained model supplied.
In that case the output VCF = raw FeatureMap (no model is applied).
