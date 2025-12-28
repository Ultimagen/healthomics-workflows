# MRD WG Analysis

## Introduction

The UG pipeline for tumor informed MRD measures the tumor fraction in cfDNA from the presence of tumor-specific SNVs. The input data is generally 3 aligned cram files:
- cfDNA (plasma)
- Tumor tissue (FFPE / FF)
- Normal tissue (buffy coat / PBMCs)

It is possible to provide the cfDNA cram file only, with an existing somatic vcf file.

The analysis is composed of three parts:
1. Tumor signature mutation calling, where the tumor and normal tissues are used for finding the tumor somatic mutations signature with somatic variant calling (by default UG Somatic Efficient DeepVariant [efficient_dv.wdl], though these can be provided from other callers).
2. Single Read SNV pipeline, where all the SNV candidates compared to the reference genome are extracted from the cfDNA cram file to a FeatureMap vcf, annotated and assigned a quality score (SNVQ).
3. Intersection and MRD data analysis, where the FeatureMap and signature are intersected and filtered, then reads supporting the tumor mutations are counted and a circulating tumor variant allele fraction (ctDNA VAF) is measured. Control signatures can be added to estimate the background noise, e.g. from other cohort patients, and in addition control signatures are generated from a somatic mutation database. 

<img src="mrd_pipeline_scheme.png" width="800"/>

This pipeline describes step #3, the intersection and MRD data analysis, once #1 and #2 are completed.

## Template descriptions

The following input templates are available for different kinds of input data:

| Template File | Description |
|---------------|-------------|
| `mrd_featuremap_template-Matched-signature-with-cohort-with-quality-filtering.json` | Use this template to run MRD using a matched signature, including cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and control signatures (suitable for EfficientDV output). |
| `mrd_featuremap_template-Matched-signature-without-cohort-with-quality-filtering.json` | Use this template to run MRD using a matched signature, without cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and signature (suitable for EfficientDV output). |
| `mrd_featuremap_template-Matched-signature-with-cohort-without-quality-filtering.json` | Use this template to run MRD using a matched signature, with cohort controls (non-matched mutation signature vcf). Quality filtering is not applied to matched and signature (suitable for vcf files without a QUAL field). |
| `mrd_featuremap_template-Healthy-without-matched-with-cohort-with-quality-filtering.json` | Use this template to run MRD on a healthy control plasma without a matched signature, with cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and signature (suitable for EfficientDV output). |

## Running the pipeline
### Step 1: Somatic variant calling
Please refer to instruction in one of:
1. WDL - efficient_dv.wdl

2. Standalone - howto-somatic-calling-efficient-dv.md

The following outputs from this workflow are needed:
1. vcf_file - A vcf file containing the mutations found in the tumor tissue, with a quality score (QUAL) assigned to each mutation. Mutations suspected as germline appear in this vcf and are filtered out and marked RefCall in the FILTER column. The output vcf file is used as input to the next step.

### Step 2: Single Read SNV (SRSNV) pipeline
Please refer to instruction in one of:
1. WDL - single_read_snv.wdl

2. Standalone - howto-single-read-snv.md

The following outputs from this workflow are needed:
1. featuremap - Output FeatureMap, a VCF file that contains a record per variant (SNV), with aggregated information per read in the FORMAT fields. Additional information about variant is encoded in the INFO fields. Additionally, a machine learning model is trained on these features to assign an SNV quality score, saved in the SNVQ FORMAT field per read. The maximal SNVQ per variant is saved in the QUAL field. 
2. featuremap_index - index of the FeatureMap vcf file
3. application_qc_h5 - A file containing statistics about the test set used to train the machine learning model, saved as a h5 file.
4. featuremap_df - a parquet file containing the training dataset with labels and predictions
5. srsnv_metadata_json - A metadata JSON file for the SNV quality model, containing information about the model, features, and training parameters


### Step 3: Intersection and MRD data analysis 
In this stage the FeatureMap and signature are intersected, reads supporting the tumor mutations are counted and the circulating tumor variant allele fraction (ctDNA VAF) is measured. ctDNA VAF can be used to estimate the tumor fraction in plasma. 
In this stage control signatures can (and should) be added to estimate the background noise. In addition a mutation database is used for estimating background noise (see below).
The WDL used in this stage is: mrd_featuremap.wdl
Standalone - howto-mrd-featuremap.md

Either when using the WDL or running as standalone, the following inputs are needed:
1. General parameters

  a. base_file_name - Sample's base file name / sample name, will be used in the output files' name

  b. include_regions - region/s to which the analysis will be limited, multiple regions would be intersected. Default:
    
    [
      "gs://concordanz/hg38/UG-High-Confidence-Regions/v1.3/ug_hcr.bed"
    ]
    or 
    [
      "s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v1.3/ug_hcr.bed"
    ]

  c. exclude_regions - regions that will be excluded from the analysis, supporting bed and vcf formats. Default:

    1. [GNOMAD](https://gnomad.broadinstitute.org/): common population variants.

    2. [db_snp](https://www.ncbi.nlm.nih.gov/snp/): common population variants.

    3. MRD_blacklist: loci with high error rate in an internal HapMap project in UG

    Optionaly, one can add a list of germline variants of the corresponding sample, to exclude from MRD analysis.
    
    [
      "gs://concordanz/hg38/somatic/af-only-gnomad.hg38.snps.AF_over_1e-3.vcf.gz",
      "gs://concordanz/hg38/somatic/Homo_sapiens_assembly38.dbsnp138.chr1-22XY.snps.vcf.gz",
      "gs://concordanz/hg38/annotation_intervals/UG_MRD_blacklist_v0.bed.gz"
    ]
    or
    [
      "s3://ultimagen-workflow-resources-us-east-1/hg38/somatic/af-only-gnomad.hg38.snps.AF_over_1e-3.vcf.gz",
      "s3://ultimagen-workflow-resources-us-east-1/hg38/somatic/Homo_sapiens_assembly38.dbsnp138.chr1-22XY.snps.vcf.gz",
      "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/UG_MRD_blacklist_v0.bed.gz"
    ]

  d. references - Reference genome. Default:

    {
      "ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
      "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
      "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    }

2. Input signatures

  a. external_matched_signatures - a list of somatic vcf files matching the plasma sample. 
  An option to use multiple inputs is supported (e.g. signatures from different callers/tissues) and intersections are calculated for all options, but the first entry is used as the signature to sample for database controls. All loci in this signature will be excluded from the control signatures, implemented with bcftools_extra_args string (see below) before intersection with the FeatureMap.

  b. external_control_signatures - A list of control signatures that can be added to estimate the background noise, e.g. from other cohort patients. The control signatures are filtered with bcftools_extra_args.

  c. bcftools_extra_args - A string used by bcftools for filtering the matched and control signatures. Default:

    "-f PASS --type snps -m2 -M2 -i 'QUAL>10'"

  This expression filters out all non-SNVs, non-biallelic SNVs and SNVs with a quality score lower than 10. To disable filtering on quality, use:
    
    "-f PASS --type snps -m2 -M2"

  d. snv_database - a large database of whole-genome somatic cancer mutations from which variants for synthetic control signatures (also called database controls) will be drawn. Default is the PCAWG database (Nature 2020). The synthetic signatures (default: 5 synthetic signatures) are generated based on the matched signature: they have the same size and same trinucleotide motif distribution as the first matched signature. In case matched signatures are not part of the input, the sythetic signatures mimic the first control signature. The synthetic signatures appear as "db_control" signatures in the output ctdna_vaf.h5. snv_database default:
  
    "gs://concordanz/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz"
    or
    "s3://ultimagen-workflow-resources-us-east-1/hg38/pcawg/pancan_pcawg_2020.mutations_hg38_GNOMAD_dbsnp_beds.sorted.Annotated.HMER_LEN.edited.vcf.gz"

  e. n_synthetic_signatures - number of synthetic signatures to generate from the database. Default: 5

  f. diluent_germline_vcfs - optional argument. A list of vcf files which are output of germline calling of the diluent's DNA, in case of an experiment where patient's cfDNA was diluted into a cfDNA coming from a healthy donor, or a similar mixing experiment. Default: empty array []
  

3. Inputs from Single Read SNV pipeline

  a. cfdna_featuremap - SRSNV output: featuremap

  b. cfdna_featuremap_index - SRSNV output: featuremap_index

  c. featuremap_df_file - SRSNV output: featuremap_df

  d. srsnv_metadata_json - SRSNV output: srsnv_metadata_json
    
4. Analysis filters

  a. mapping_quality_threshold - minimum mapping quality of reads supporting a mutation, used only when estimating effective coverage. The FeatureMap is assumed to be filtered with the same value in the SRSNV pipeline. Default: 0
  
  b. mrd_analysis_params - filters used in the final analysis steps.
    1. "signature_filter_query": the default is: "(norm_coverage <= 2.5) and (norm_coverage >= 0.6)"
    It filters out variants found in regions of extreme coverage of the cfDNA sample
    2. "signature_filter_query": the default is: "filt>0 and snvq>60 and mapq>=60"
    Taking only featuremap entries with pass filter, have high SNVQ and maximal mapping quality.
    In order to take only mixed reads from a ppmSeq data, the following read_filter_query should be applied:
    "read_filter_query" : "(st == 'MIXED') and (et == 'MIXED')"

## Manual execution (outside WDL)

### Environment / dockers
Pull (matching versions referenced in workflows/single_read_snv/tasks/globals.wdl in this repository):

1. ugbio_core_docker
2. bcftools_docker
3. ugbio_mrd_docker
4. mosdepth_docker

### Step-by-step commands
Below is a bash script emulating the mrd_featuremap.wdl. Variable names mirror WDL inputs; example corresponds to test data provided (see file paths for test data in the last section of this deocument). See input templates for full filenames, only base names are used below for clarity.

```bash
# Set base name & inputs
BASE_FILE_NAME="Pa_46_333_LuNgs_08"
CFDNA_FEATUREMAP="Pa_46_333_LuNgs_08.featuremap.chr20.vcf.gz"
CFDNA_FEATUREMAP_INDEX="Pa_46_333_LuNgs_08.featuremap.chr20.vcf.gz.tbi"
CFDNA_CRAM_BAM="Pa_46.333_LuNgs_08.Lb_744.chr20.cram"
CFDNA_CRAM_BAM_INDEX="Pa_46.333_LuNgs_08.Lb_744.chr20.cram.crai"
FEATUREMAP_DF_FILE="Pa_46_333_LuNgs_08.featuremap_df.parquet"
SRSNV_METADATA_JSON="Pa_46_333_LuNgs_08.srsnv_metadata.json"

# External signatures (arrays)
EXTERNAL_MATCHED_SIGNATURES=("Pa_46_FreshFrozen.ann.chr20.vcf.gz")
EXTERNAL_CONTROL_SIGNATURES=("Pa_67_FFPE.ann.chr20.vcf.gz")
DILUENT_GERMLINE_VCFS=()  # Optional, empty in this example

# Reference files
REF_FASTA="Homo_sapiens_assembly38.fasta"
REF_FASTA_INDEX="Homo_sapiens_assembly38.fasta.fai"
REF_DICT="Homo_sapiens_assembly38.dict"

# Regions and database
INCLUDE_REGIONS=("ug_hcr.bed")
EXCLUDE_REGIONS=(
    "af-only-gnomad.hg38.snps.AF_over_1e-3.vcf.gz"
    "Homo_sapiens_assembly38.dbsnp138.chr1-22XY.snps.vcf.gz"
    "UG_MRD_blacklist_v0.bed.gz"
)
SNV_DATABASE="pancan_pcawg_2020.chr20.vcf.gz"

# Parameters
MAPPING_QUALITY_THRESHOLD=0
N_SYNTHETIC_SIGNATURES=5
BCFTOOLS_EXTRA_ARGS="-f PASS --type snps -m2 -M2 -i 'QUAL>10'"
SIGNATURE_FILTER_QUERY="(norm_coverage <= 2.5) and (norm_coverage >= 0.6)"
READ_FILTER_QUERY="filt>0 and snvq>60 and mapq>=60"

# Initialize arrays for tracking filtered signatures
filtered_matched_signatures=()
filtered_control_signatures=()
db_signatures=()
padded_diluent_files=()

# Define signature filtering function
filter_signatures() {
    local signature_type="$1"
    local array_name="$2"
    local bcftools_extra_args="${3:-$BCFTOOLS_EXTRA_ARGS}"  # Optional, defaults to global value
    shift 3
    local exclude_regions=("$@")
    
    # Extract signature files array
    eval "local signature_files=(\"\${${array_name}[@]}\")"
    
    echo "Filtering ${signature_type} signatures..."
    # run inside bcftools_docker
    
    for input_vcf in "${signature_files[@]}"; do
        basename_vcf=$(basename "$input_vcf" .vcf.gz)
        output_vcf="${basename_vcf}.filtered.vcf.gz"
        
        # Start with basic filtering and include regions
        cmd="bcftools view --threads 4 $bcftools_extra_args \"$input_vcf\""
        
        # Add include regions (pipe through each one)
        for include_region in "${INCLUDE_REGIONS[@]}"; do
            cmd="$cmd | bcftools view - -T \"$include_region\""
        done
        
        # Add exclude regions (pipe through each one)
        for exclude_region in "${exclude_regions[@]}"; do
            if [[ -n "$exclude_region" ]]; then
                cmd="$cmd | bcftools view - -T ^\"$exclude_region\""
            fi
        done
        
        # Finalize command with output
        cmd="$cmd -Oz -o \"$output_vcf\""
        
        # Execute the command
        echo "Executing command: $cmd"
        eval "$cmd"
        bcftools index -t "$output_vcf"
    done
}

# 1. OPTIONAL: Pad diluent VCF files (if provided)
# Skip in this example since DILUENT_GERMLINE_VCFS is empty
if [[ ${#DILUENT_GERMLINE_VCFS[@]} -gt 0 ]]; then
    echo "Processing diluent germline VCFs..."
    # run inside ugbio_core_docker
    for input_vcf in "${DILUENT_GERMLINE_VCFS[@]}"; do
        basename_vcf=$(basename "$input_vcf" .vcf.gz)
        
        # Extract variant positions and pad them
        bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' "$input_vcf" | \
        awk -F'\t' 'BEGIN{OFS="\t"}{
            s=$2;
            m=length($4);
            n=split($5,a,",");
            for(i=1;i<=n;i++) if(length(a[i])>m) m=length(a[i]);
            print $1, s, s+m
        }' | gzip > variants.bed.gz
        
        # Create genome file for bedtools
        cut -f1,2 "$REF_FASTA_INDEX" > genome.txt
        
        # Pad variants using bedtools
        zcat variants.bed.gz | bedtools slop -i stdin -g genome.txt -b 2 | \
        gzip > "${basename_vcf}.padded.bed.gz"
        
        padded_diluent_files+=("${basename_vcf}.padded.bed.gz")
    done
fi

# 2. Filter matched signatures
all_exclude_regions=("${EXCLUDE_REGIONS[@]}" "${padded_diluent_files[@]}")
if [[ ${#EXTERNAL_MATCHED_SIGNATURES[@]} -gt 0 ]]; then
    filter_signatures "matched" "EXTERNAL_MATCHED_SIGNATURES" "$BCFTOOLS_EXTRA_ARGS" "${all_exclude_regions[@]}"
    # Collect filtered outputs
    for input_vcf in "${EXTERNAL_MATCHED_SIGNATURES[@]}"; do
        basename_vcf=$(basename "$input_vcf" .vcf.gz)
        filtered_matched_signatures+=("${basename_vcf}.filtered.vcf.gz")
    done
fi

# 3. Filter control signatures (exclude matched signature loci too)
control_exclude_regions=("${EXCLUDE_REGIONS[@]}" "${EXTERNAL_MATCHED_SIGNATURES[@]}" "${padded_diluent_files[@]}")
if [[ ${#EXTERNAL_CONTROL_SIGNATURES[@]} -gt 0 ]]; then
    filter_signatures "control" "EXTERNAL_CONTROL_SIGNATURES"  "$BCFTOOLS_EXTRA_ARGS" "${control_exclude_regions[@]}"
    # Collect filtered outputs
    for input_vcf in "${EXTERNAL_CONTROL_SIGNATURES[@]}"; do
        basename_vcf=$(basename "$input_vcf" .vcf.gz)
        filtered_control_signatures+=("${basename_vcf}.filtered.vcf.gz")
    done
fi

# 4. Filter database and generate synthetic controls
echo "Filtering SNV database and generating synthetic controls..."
# run inside bcftools_docker
if [[ -n "$SNV_DATABASE" ]]; then
    # Create a temporary array with just the database file
    temp_db_array=("$SNV_DATABASE")
    filter_signatures "database" "temp_db_array" " " "${control_exclude_regions[@]}"
    
    # Get the filtered database file
    basename_db=$(basename "$SNV_DATABASE" .vcf.gz)
    filtered_db_vcf="${basename_db}.filtered.vcf.gz"
else
    filtered_db_vcf="filtered_db.vcf.gz"
fi

# run inside ugbio_mrd_docker
# Use first filtered signature as reference (prefer matched over control)
reference_signature="${filtered_matched_signatures[0]:-${filtered_control_signatures[0]}}"

generate_synthetic_signatures \
  --signature_vcf "$reference_signature" \
  --db_vcf "$filtered_db_vcf" \
  --n_synthetic_signatures $N_SYNTHETIC_SIGNATURES \
  --ref_fasta "$REF_FASTA" \
  --output_dir ./

# Collect generated database signatures
db_signatures=(syn*.vcf.gz)

# 5. Extract coverage over all signatures
echo "Extracting coverage over signature loci..."
# Combine all signature VCF files for coverage calculation
all_signature_files=("${filtered_matched_signatures[@]}" "${filtered_control_signatures[@]}" "${db_signatures[@]}")

# run inside ugbio_core_docker
# Convert VCF loci to BED format
echo "Combining all VCF loci into one BED file..."
for vcf in "${all_signature_files[@]}"; do
    zcat "$vcf" | grep -v "^#" | awk '{print $1"\t"($2-1)"\t"$2}' >> combined_loci.bed
done

# Sort and merge overlapping regions
echo "Sorting and merging the combined BED..."
sort -k1,1 -k2,2n combined_loci.bed | bedtools merge > merged_loci.bed

# run inside mosdepth_docker
# Extract coverage using mosdepth
echo "Extracting coverage from CRAM for the specified loci..."
mosdepth --by merged_loci.bed -f "$REF_FASTA" -Q $MAPPING_QUALITY_THRESHOLD --fast-mode \
"$BASE_FILE_NAME" "$CFDNA_CRAM_BAM"

coverage_bed="${BASE_FILE_NAME}.per-base.bed.gz"
coverage_bed_index="${BASE_FILE_NAME}.per-base.bed.gz.csi"

# 6. Intersect FeatureMap with signatures
echo "Intersecting FeatureMap with signatures..."
# run inside ugbio_mrd_docker
intersection_files=()
intersection_parquet_files=()

# Process matched signatures
for signature in "${filtered_matched_signatures[@]}"; do
    featuremap_base=$(basename "$CFDNA_FEATUREMAP")
    signature_base=$(basename "$signature")
    signature_type="matched"
    output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
    
    # Perform intersection
    bcftools isec -n=2 -w1 "$CFDNA_FEATUREMAP" "$signature" -Oz -o "$output_vcf"
    bcftools index -t "$output_vcf"
    intersection_files+=("$output_vcf")
    
    # Convert to parquet if not empty
    if [[ $(bcftools view "$output_vcf" -H | wc -l) -gt 0 ]]; then
        output_parquet="${output_vcf%.vcf.gz}.parquet"
        featuremap_to_dataframe --in "$output_vcf" --out "$output_parquet" --jobs 4 --drop-format AD GT
        intersection_parquet_files+=("$output_parquet")
    fi
done

# Process control signatures
for signature in "${filtered_control_signatures[@]}"; do
    featuremap_base=$(basename "$CFDNA_FEATUREMAP")
    signature_base=$(basename "$signature")
    signature_type="control"
    output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
    
    bcftools isec -n=2 -w1 "$CFDNA_FEATUREMAP" "$signature" -Oz -o "$output_vcf"
    bcftools index -t "$output_vcf"
    intersection_files+=("$output_vcf")
    
    if [[ $(bcftools view "$output_vcf" -H | wc -l) -gt 0 ]]; then
        output_parquet="${output_vcf%.vcf.gz}.parquet"
        featuremap_to_dataframe --in "$output_vcf" --out "$output_parquet" --jobs 4 --drop-format AD GT
        intersection_parquet_files+=("$output_parquet")
    fi
done

# Process database control signatures
for signature in "${db_signatures[@]}"; do
    featuremap_base=$(basename "$CFDNA_FEATUREMAP")
    signature_base=$(basename "$signature")
    signature_type="db_control"
    output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
    
    bcftools isec -n=2 -w1 "$CFDNA_FEATUREMAP" "$signature" -Oz -o "$output_vcf"
    bcftools index -t "$output_vcf"
    intersection_files+=("$output_vcf")
    
    if [[ $(bcftools view "$output_vcf" -H | wc -l) -gt 0 ]]; then
        output_parquet="${output_vcf%.vcf.gz}.parquet"
        featuremap_to_dataframe --in "$output_vcf" --out "$output_parquet" --jobs 4 --drop-format AD GT
        intersection_parquet_files+=("$output_parquet")
    fi
done

# 7. MRD data analysis - Generate final report and ctDNA VAF calculations
echo "Running MRD data analysis..."
# run inside ugbio_mrd_docker

# Build arguments for different signature types
matched_sigs_args=""
if [[ ${#filtered_matched_signatures[@]} -gt 0 ]]; then
    matched_sigs_args="--matched-signatures-vcf $(IFS=' '; echo "${filtered_matched_signatures[*]}")"
fi

control_sigs_args=""
if [[ ${#filtered_control_signatures[@]} -gt 0 ]]; then
    control_sigs_args="--control-signatures-vcf $(IFS=' '; echo "${filtered_control_signatures[*]}")"
fi

db_sigs_args=""
if [[ ${#db_signatures[@]} -gt 0 ]]; then
    db_sigs_args="--db-control-signatures-vcf $(IFS=' '; echo "${db_signatures[*]}")"
fi

# Run final MRD data analysis
generate_report \
    --intersected-featuremaps $(IFS=' '; echo "${intersection_parquet_files[*]}") \
    --coverage-bed "$coverage_bed" \
    $matched_sigs_args \
    $control_sigs_args \
    $db_sigs_args \
    --output-dir "$PWD" \
    --output-basename "$BASE_FILE_NAME" \
    --signature-filter-query "$SIGNATURE_FILTER_QUERY" \
    --read-filter-query "$READ_FILTER_QUERY" \
    --featuremap-file "$FEATUREMAP_DF_FILE" \
    --srsnv-metadata-json "$SRSNV_METADATA_JSON"

echo "MRD analysis complete!"
echo "Output files:"
echo "  - ${BASE_FILE_NAME}.features.parquet"
echo "  - ${BASE_FILE_NAME}.signatures.parquet"
echo "  - ${BASE_FILE_NAME}.mrd_data_analysis.html"
echo "  - ${BASE_FILE_NAME}.ctdna_vaf.h5"
```

## Summary of the important MRD output files
- report_html (automated analysis of the results in html format)
- features_dataframe (python pandas dataframe of all the substitutions in the cfDNA sample after intersection with all the signatures, parquet format)
- signatures_dataframe (python pandas dataframe of all the variants in all the signatures, parquet format)
- ctdna_vaf.h5 (dataframes of ctDNA VAF and supporting reads per locus)

## Summary of the relevant files
- **WDLs:**
  - mrd_featuremap.wdl
  - single_read_snv.wdl
  - efficient_dv.wdl

## test files
```
{cfdna_featuremap}:s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_46_333_LuNgs_08.featuremap.chr20.vcf.gz
{cfdna_featuremap_index}:s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_46_333_LuNgs_08.featuremap.chr20.vcf.gz.tbi
{cfdna_cram_bam}:s3://ultimagen-workflow-resources-us-east-1/test_data/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram
{cfdna_cram_bam_index}:s3://ultimagen-workflow-resources-us-east-1/test_data/single_read_snv/Pa_46.333_LuNgs_08.Lb_744.chr20.cram.crai
{external_matched_signatures}: ["s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_46_FreshFrozen.ann.chr20.vcf.gz"]
{external_control_signatures}: ["s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_67_FFPE.ann.chr20.vcf.gz"]
{featuremap_df_file}:s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_46_333_LuNgs_08.featuremap_df.parquet
{snv_database}:s3://ultimagen-workflow-resources-us-east-1/hg38/pcawg/pancan_pcawg_2020.chr20.vcf.gz
{srsnv_metadata_json}:s3://ultimagen-workflow-resources-us-east-1/test_data/mrd/Pa_46_333_LuNgs_08.srsnv_metadata.json
```