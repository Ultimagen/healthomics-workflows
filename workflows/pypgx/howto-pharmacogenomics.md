# Pharmacogenomics Analysis of WGS Data Using PyPGx

## Overview

The PyPGx workflow performs comprehensive pharmacogenomics analysis on whole-genome sequencing (WGS) data. It processes CRAM files, calls variants using Efficient DV, and generates star-allele representations and phenotypes for specified genes.

**Key Resources:**
- **PyPGx Documentation**: [GitHub](https://github.com/sbslee/pypgx) | [ReadTheDocs](https://pypgx.readthedocs.io/en/latest/readme.html)
- **Efficient DV Documentation**: [How-to Guide](howto-germline-calling-efficient-dv.md)

## Requirements

### Input Files

1. **CRAM File and Index**
   - WGS reads aligned to the GRCh38 reference genome
   - Must be sorted and duplicate-marked
   - Requires corresponding CRAI index file

2. **Gene Symbols**
   - List of gene symbols to analyze
   - The workflow supports 87 genes detailed [here](https://github.com/sbslee/pypgx)

3. **Reference Genome Files**
   - FASTA file with GRCh38 reference genome
   - FAI index file
   - Sequence dictionary file
   
   **Example (publicly available hg38 reference):**
   ```
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
   ```

4. **Optional VCF File**
   - It is preferable to use EfficientDV for variant calling, but any other pre-existing variant calls can be used instead

### Docker Images

The workflow uses multiple Docker images:

1. **PyPGx Docker:**
   ```
   ultimagenomics/ugbio_pypgx:0.26.0-r2
   ```

2. **Efficient DV Dockers:**
   
   a) **make_examples Docker:**
   ```
  ultimagenomics/make_examples
   ```
   
   b) **call_variants Docker:**
   ```
   ultimagenomics/make_examples
   ```

### Hardware Requirements

- **CPU**: Minimum 2 cores for most tasks
- **Memory**: 4-16 GB depending on task complexity
- **Disk**: Sufficient space for intermediate files and CRAM processing

## Workflow Steps

The main steps of the PyPGx workflow:

### 1. Variant Calling with Efficient DV

If no input VCF is provided, the workflow uses Efficient DV to call variants. For pharmacogenomics analysis, two parameters are modified to ensure variants are called in highly homologous regions (e.g., CYP2D6):

| Parameter name in WDL         | Argument in tool command line | Value |
| ----------------------------- | ----------------------------- | ----- |
| candidate_min_mapping_quality | cgp-min-mapping-quality       | 0     |
| pileup_min_mapping_quality    | min-map                       | 0     |

**Extracting Gene Regions (Recommended):**
To reduce run time, limit variant calling to relevant genomic regions. A bed file with these regions can be produced by PyPGx:

```bash
pypgx create-regions-bed --assembly GRCh38 --add-chr-prefix \
  --genes ABCG2 CYP2D6 VDR ... > gene_regions.bed
```

Use the bed in the --bed argument of the `tool` command line, or convert it to interval_list and use as target_intervals argument in the WDL.

### 2. Depth of Coverage Analysis

Coverage statistics are computed for genes that may exhibit structural variations:

**Primary Coverage Analysis:**
```bash
pypgx prepare-depth-of-coverage --assembly GRCh38 \
  sample.coverage.zip sample.cram
```

**Control Gene Analysis:**
For normalization purposes, compute coverage of a control gene (e.g., VDR):
```bash
pypgx compute-control-statistics --assembly GRCh38 VDR \
  sample_VDR_coverage.zip sample.cram
```

### 3. PyPGx Pipeline Execution

The main pharmacogenomics pipeline performs phasing, genotyping, and phenotype prediction:

```bash
pypgx run-ngs-pipeline CYP2D6 outputs \
  --assembly GRCh38 \
  --variants sample.vcf.gz \
  --depth-of-coverage sample.coverage.zip \
  --control-statistics sample_VDR_coverage.zip
```

### 4. Results Interpretation

The pipeline generates several [archive files](https://pypgx.readthedocs.io/en/latest/readme.html#archive-file-semantic-type-and-metadata). Of those the key output files are:

- **results.zip**: Main results file containing genotypes and phenotypes
- **copy-number-profile.png**: Visual representation of copy number variations
- **allele-fraction-profile.png**: Visual representation of allele fractions

**Viewing Results:**
```bash
pypgx print-data results.zip
```

For detailed interpretation of results, see the [PyPGx documentation](https://pypgx.readthedocs.io/en/latest/readme.html#results-interpretation).

