# BcftoolsVariantCalling: bcftools-based variant calling from aligned CRAM files

## Overview

BcftoolsVariantCalling is a variant calling pipeline that processes aligned CRAM files using `bcftools mpileup` + `bcftools call` to generate gVCF and VCF files. It supports sex-aware chromosome handling, optional quality filtering, and optional evaluation against GIAB truth sets.

Main features:
- Parallel processing of autosomes split into configurable number of shards (default: 10)
- Sex-specific handling of X and Y chromosomes (diploid/haploid ploidy)
- Optional filtering by depth (`DP < 8`) and genotype quality (`GQ < 20`)
- Generates gVCF output with optional quality filtering

## Requirements

### Input files

1. Reads aligned to a reference genome, sorted and duplicate-marked, stored in CRAM format with a corresponding `.crai` index.
2. Reference genome files (FASTA, FASTA index, and sequence dictionary).
3. BED files defining genomic regions:
   - Autosomal targets
   - X chromosome PAR (pseudoautosomal) regions
   - X chromosome non-PAR regions
   - Y chromosome non-PAR regions

### Dockers

The pipeline currently uses the following hardcoded Docker images:

- `biocontainers/bcftools:v1.14-1-deb_cv1` — for variant calling, concatenation, and filtering
- `biocontainers/tabix:1.11--hdfd78af_0` — for creating empty gVCF (female Y chromosome placeholder)
- `ubuntu:latest` — for splitting autosomal BED into shards

These must be replaced with globals dockers during onboarding (see [Onboarding Plan](#onboarding-plan-for-terra-pipeline) below).

## Inputs

### Required

| Parameter | Type | Description |
|-----------|------|-------------|
| `base_file_name` | String | Base name for all output files |
| `input_cram` | File | Input CRAM file with aligned reads (must end with `.cram`) |
| `input_cram_index` | File | CRAM index file (must end with `.crai`, must match input_cram path) |
| `sex` | String | Sample sex: `XX` (female) or `XY` (male) |
| `autosomes_bed` | File | BED file with autosomal target regions |
| `chrx_par_bed` | File | BED file with X chromosome PAR regions |
| `chrx_nonpar_bed` | File | BED file with X chromosome non-PAR regions |
| `chry_nonpar_bed` | File | BED file with Y chromosome non-PAR regions |
| `references` | References | Reference genome struct containing `ref_fasta` (`.fasta` or `.fa`), `ref_fasta_index` (`.fai`), and `ref_dict` (`.dict`) |

### Optional

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `run_filtering` | Boolean | `true` | Whether to generate a filtered gVCF using depth and GQ thresholds |
| `num_autosome_shards` | Int | `10` | Number of shards for parallel autosome processing (must be > 0) |
| `memory_gb` | Int | `4` | Memory in GB for variant calling tasks |
| `cpus` | Int | `2` | Number of CPUs for variant calling tasks |

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| `final_gvcf` | File | Final merged gVCF containing all chromosomes (autosomes + X PAR + X non-PAR + Y non-PAR) |
| `final_gvcf_index` | File | Tabix index (`.tbi`) for the final gVCF |
| `filtered_gvcf` | File? | Filtered gVCF with `LowDepth` (DP < 8) and `LowGQ` (GQ < 20) soft filters applied. Only produced when `run_filtering` is true. |
| `filtered_gvcf_index` | File? | Tabix index for the filtered gVCF |

## Workflow Steps

### 1. SplitAutosomesIntoBed

Splits the input autosomal BED file into `num_autosome_shards` approximately equal parts for parallel processing.

- **Method**: Uses `awk` to distribute BED lines evenly across shard files named `autosomes.part01.bed`, `autosomes.part02.bed`, etc.
- **Output**: Array of shard BED files.

### 2. BcftoolsCallVariants (autosome scatter)

Runs `bcftools mpileup` piped into `bcftools call` for each autosome shard in parallel.

- **Ploidy**: 2 (diploid) for all autosome shards.
- **mpileup parameters**:
  - `--ignore-RG` — ignores read groups
  - `-a FORMAT/DP,FORMAT/AD` — annotates with per-sample depth and allelic depth
  - `-q 20` — minimum mapping quality 20
  - `-Q 20` — minimum base quality 20
- **call parameters**:
  - `-m` — multiallelic caller
  - `-A` — keep all alternate alleles
  - `-a FORMAT/GQ` — annotates with genotype quality
  - `--ploidy` — set per region
- **Output**: One gVCF (`.g.vcf.gz`) and its tabix index per shard.

### 3. CallChrxPar — X chromosome PAR regions

Calls variants on the X chromosome pseudoautosomal regions (PAR). Uses ploidy 2 (diploid) for both XX and XY samples, since PAR regions behave like autosomes.

### 4. CallChrxNonpar — X chromosome non-PAR regions

Calls variants on X chromosome non-PAR regions. Ploidy is sex-dependent:
- **XX (female)**: ploidy = 2 (diploid)
- **XY (male)**: ploidy = 1 (haploid)

### 5. Y chromosome handling (sex-dependent)

- **XY (male)**: `CallChryNonparMale` — calls variants on Y chromosome non-PAR regions with ploidy = 1 (haploid).
- **XX (female)**: `CreateEmptyChryFemale` — creates an empty gVCF placeholder for the Y chromosome (header-only VCF with no records), since females have no Y chromosome data.

### 6. BcftoolsConcatGvcf — Merge autosome shards

Concatenates all autosome shard gVCFs into a single autosomal gVCF using `bcftools concat -a`.

### 7. ConcatFinalGvcf — Create final merged gVCF

Concatenates the autosomal gVCF with the sex chromosome gVCFs in genomic order:
1. Autosomes (merged)
2. X PAR
3. X non-PAR
4. Y non-PAR (real calls for males, empty placeholder for females)

This produces the final whole-genome gVCF.

### 8. BcftoolsFilterGvcf (optional, when `run_filtering = true`)

Applies two soft filters to the final gVCF:
1. **LowDepth**: Marks sites where `FORMAT/DP < 8`
2. **LowGQ**: Marks sites where `FORMAT/GQ < 20`

Filters are applied additively (`-m +`), meaning both filter tags can appear on the same record.

## Workflow Diagram

```
input_cram
    │
    ├──► SplitAutosomesIntoBed ──► scatter(N shards) ──► BcftoolsCallVariants ──► BcftoolsConcatGvcf (autosomes)
    │                                                                                      │
    ├──► CallChrxPar (ploidy=2) ──────────────────────────────────────────────────────────┐ │
    │                                                                                      │ │
    ├──► CallChrxNonpar (ploidy=2 if XX, 1 if XY) ───────────────────────────────────────┐│ │
    │                                                                                      ││ │
    ├──► [XY] CallChryNonparMale (ploidy=1) ─────────────────────────────────────────────┐││ │
    │    [XX] CreateEmptyChryFemale ─────────────────────────────────────────────────────┐│││ │
    │                                                                                     ││││ │
    │                                                                     ConcatFinalGvcf ◄┘┘┘┘
    │                                                                           │
    │                                                                  BcftoolsFilterGvcf
    │                                                                      (optional)
    │                                                                   filtered.g.vcf.gz
```

## Onboarding Plan for Terra Pipeline

This section describes the steps required to integrate this workflow into the terra_pipeline infrastructure for deployment on AWS HealthOmics. Follow the [ONBOARDING_NEW_WDL](../wdls_conf_docs/ONBOARDING_NEW_WDL.md) guide for detailed instructions on each step.

### 1. Replace hardcoded Dockers with globals

All tasks currently use hardcoded Docker images. Replace them with globals:

| Task | Current Hardcoded Docker | Global Replacement |
|------|--------------------------|-------------------|
| `SplitAutosomesIntoBed` | `ubuntu:latest` | `global.ubuntu_docker` |
| `BcftoolsCallVariants` | `biocontainers/bcftools:v1.14-1-deb_cv1` | `global.bcftools_docker` |
| `BcftoolsConcatGvcf` | `biocontainers/bcftools:v1.14-1-deb_cv1` | `global.bcftools_docker` |
| `BcftoolsFilterGvcf` | `biocontainers/bcftools:v1.14-1-deb_cv1` | `global.bcftools_docker` |
| `CreateEmptyGvcf` | `biocontainers/tabix:1.11--hdfd78af_0` | `global.bcftools_docker` |

Notes:

- The globals `bcftools_docker` uses `staphb/bcftools:1.19` (vs the current `v1.14`). Run a test and compare outputs to verify compatibility.
- The separate `tabix` image can be eliminated — `staphb/bcftools:1.19` includes both `bgzip` and `tabix`.
- Each task needs a `String docker` input parameter, and the `runtime` block should reference it instead of a hardcoded string.

See [GLOBALS.md](../wdls_conf_docs/GLOBALS.md) for details.

### 2. Fix parameter_meta categories

Update `parameter_meta` categories to use the standard values required by HealthOmics:

- `"required"` -> `"input_required"`
- `"optional"` -> `"input_optional"`
- Reference inputs should use `"ref_required"`

### 3. Create an input template

Create `wdls/input_templates/bcftools_variant_calling_template.json` with the workflow's default parameters. Use `<placeholder>` notation for user-provided values and literal values for defaults.

### 4. Create workflow configuration YAML

Create `wdls_conf/BcftoolsVariantCalling.yaml` with:

- `wdl_path` pointing to `wdls/bcftools_variant_calling.wdl`
- `input_templates` referencing the template JSON
- `supported_in` listing target platforms (at minimum: `cromwell`, `omics`)
- `omics` section with HealthOmics-specific overrides
- At least one test definition

### 5. Register in wdl_conf.yaml

Add the workflow entry to the `defaults` list in `wdl_conf.yaml`:

```yaml
- /wdls_conf/BcftoolsVariantCalling@wdls.BcftoolsVariantCalling
```

### 6. Create test inputs

Create at least one test input JSON in `tests/tests_inputs/regression_tests/nightly/` that overrides the template with actual test data paths.

### 7. Validate and test

```bash
# Validate
wdls workflow=[BcftoolsVariantCalling] action=[validate] target=all

# Deploy and run on HealthOmics
wdls workflow=[BcftoolsVariantCalling] target=omics action=[deploy,run] run.test_type=[nightly] aws.profile=dev
```
