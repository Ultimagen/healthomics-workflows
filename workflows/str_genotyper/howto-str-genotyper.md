# STR Genotyping with Alignment-Based Caller

## Overview

The STR Genotyper workflow performs Short Tandem Repeat (STR) genotyping on whole-genome sequencing (WGS) data using an alignment-based approach. It uses Smith-Waterman alignment to accurately determine repeat lengths at targeted STR loci and generates genotype calls in multiple output formats.

**Key Features:**
- Alignment-based STR length calling using Smith-Waterman alignment
- Support for clinically relevant STR loci (e.g., AFF2, AR, ATN1, ATXN1-10, C9ORF72, etc.)
- Multiple output formats: VCF, BED, and detailed CSV reports
- Configurable quality filtering via alignment score thresholds

## Requirements

### Input Files

1. **CRAM File and Index**
   - WGS reads aligned to the GRCh38 reference genome
   - Must be sorted and duplicate-marked
   - Requires corresponding CRAI index file

2. **Variant Catalog**
   - JSON file containing STR locus definitions
   - Each locus entry includes LocusId, LocusStructure (repeat motif), ReferenceRegion, and VariantType

3. **Reference Genome Files**
   - FASTA file with GRCh38 reference genome
   - FAI index file
   
   **Example (publicly available hg38 reference):**
   ```
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
   gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
   ```

### Docker Image

```
337532070941.dkr.ecr.us-east-1.amazonaws.com/str-genotyper:1.0.0
```

### Hardware Requirements

- **CPU**: 8 cores recommended (configurable via `threads` parameter)
- **Memory**: 16 GB default (can be overridden via `memory_gb_override`)
- **Disk**: Automatically calculated based on input file sizes plus 10 GB overhead

## Variant Catalog Format

The variant catalog is a JSON file containing an array of STR locus definitions. Each locus entry requires the following fields:

| Field | Description | Example |
| ----- | ----------- | ------- |
| LocusId | Unique identifier for the STR locus | "ATN1" |
| LocusStructure | Repeat motif pattern using regex syntax | "(CAG)*" |
| ReferenceRegion | Genomic coordinates (chr:start-end) | "chr12:6936716-6936773" |
| VariantType | Type of variant | "Repeat" |

**Example variant catalog entry:**
```json
[
    {
        "LocusId": "ATN1",
        "LocusStructure": "(CAG)*",
        "ReferenceRegion": "chr12:6936716-6936773",
        "VariantType": "Repeat"
    },
    {
        "LocusId": "ATXN7",
        "LocusStructure": "(GCA)*(GCC)+",
        "ReferenceRegion": [
            "chr3:63912684-63912714",
            "chr3:63912714-63912726"
        ],
        "VariantId": [
            "ATXN7",
            "ATXN7_GCC"
        ],
        "VariantType": [
            "Repeat",
            "Repeat"
        ]
    }
]
```

## Workflow Steps

### Running via WDL

The STR genotyper is available as a WDL workflow (`str_genotyper.wdl`). Create an inputs JSON file:

```json
{
    "STRGenotyper.base_file_name": "sample_name",
    "STRGenotyper.cram_file": "s3://path/to/sample.cram",
    "STRGenotyper.cram_index": "s3://path/to/sample.cram.crai",
    "STRGenotyper.variant_catalog": "s3://path/to/variant_catalog.json",
    "STRGenotyper.references": {
        "ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    },
    "STRGenotyper.ref_padding": 500,
    "STRGenotyper.min_repeat": 1,
    "STRGenotyper.max_repeat": 100,
    "STRGenotyper.min_score_ratio": 0.85,
    "STRGenotyper.spanning_flank_bases": 5,
    "STRGenotyper.threads": 8
}
```

### Running via Command Line

The underlying tool can be run directly using:

```bash
python -m alignment_str_len_caller.main \
    --cram-file sample.cram \
    --cram-index sample.cram.crai \
    --reference Homo_sapiens_assembly38.fasta \
    --reference-index Homo_sapiens_assembly38.fasta.fai \
    --variant-catalog variant_catalog.json \
    --ref-padding 500 \
    --min-repeat 1 \
    --max-repeat 100 \
    --min-score-ratio 0.85 \
    --spanning-flank-bases 5 \
    --threads 8 \
    --output-dir ./results \
    --output-prefix sample_name
```

## Parameters

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| ref_padding | 500 | Number of bases to extend around the STR repeat region when building auxiliary references for alignment |
| min_repeat | 1 | Minimum number of repeat units to include in auxiliary reference sequences |
| max_repeat | 100 | Maximum number of repeat units to include in auxiliary reference sequences |
| min_score_ratio | 0.85 | Minimum ratio of alignment score to theoretical maximum (0.0-1.0). Lower values allow more alignments; 1.0 requires perfect alignment |
| spanning_flank_bases | 5 | Minimum number of bases that must align on each side of the STR repeat region for a read to be considered "spanning" |
| threads | 8 | Number of threads for parallel processing |

## Output Files

The workflow generates the following output files:

| File | Format | Description |
| ---- | ------ | ----------- |
| `{sample}_detailed.csv` | CSV | Per-read alignment results containing alignment scores, repeat counts, and read metadata |
| `{sample}_summary.csv` | CSV | Per-locus summary statistics aggregated from detailed results |
| `{sample}_genotypes.bed` | BED | Final genotype calls for visualization in genome browsers (IGV, UCSC) |
| `{sample}_genotypes.vcf.gz` | VCF | Final genotype calls with per-allele support counts (ADSP, ADFL) |
| `{sample}_genotypes.vcf.gz.tbi` | TBI | Tabix index for the VCF file |

### Interpreting Results

The VCF output contains genotype calls with the following INFO/FORMAT fields:
- **ADSP**: Allele depth from spanning reads (reads that fully span the STR region)
- **ADFL**: Allele depth from flanking reads

The BED file can be loaded directly into genome browsers for visual inspection of STR calls alongside aligned reads.

## Clinically Relevant STR Loci

The variant catalog can include loci associated with various repeat expansion disorders:

| Gene | Associated Condition |
| ---- | -------------------- |
| AFF2 | Fragile XE syndrome |
| AR | Spinal and bulbar muscular atrophy |
| ATN1 | Dentatorubral-pallidoluysian atrophy |
| ATXN1-10 | Spinocerebellar ataxias |
| C9ORF72 | ALS/Frontotemporal dementia |
| FMR1 | Fragile X syndrome |
| HTT | Huntington disease |
| DMPK | Myotonic dystrophy type 1 |

