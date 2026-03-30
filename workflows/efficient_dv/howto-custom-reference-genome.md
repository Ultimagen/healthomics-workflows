# Using Custom Reference Genomes


## Overview

This version introduces reference bundles. Reference bundles are pre-defined sets of reference-specific files that are used by one or more workflows. Reference bundles make switching reference genome much easier.
When running the workflow, you select a reference genome by setting the `reference_genome` input parameter to the name of the desired genome type. Currently, only `trim_align_sort` and `efficient_dv` workflows
support reference bundles, but we are working on expanding this support.

## Built-in Reference Genomes

The workflow supports a subset of the following reference genomes (check `tasks/genome_resources.wdl` to see which are included):

| Genome Type      | Description                                               |
| ---------------- | --------------------------------------------------------- |
| `hg38`           | Human genome build 38 with alt contigs (b38) - default    |
| `b37`            | Human genome build 37 (b19)                               |
| `hg38_taps`      | hg38 with Lambda/pUC19 methylation control references     |
| `hg38_nist_v3`   | hg38 NIST reference                                       |
| `hg38_rna_seq`   | hg38 without ALT/HLA/Decoy contigs (for RNA-seq)         |
| `mm10`           | Mouse genome build 10 (GRCm38)                           |
| `mm10_methyl`    | Mouse genome build 10 (with methylation indices)          |
| `mm39`           | Mouse genome build 39 (GRCm39)                           |

To use a built-in genome, simply set the `reference_genome` parameter in your workflow input JSON:

```json
{
    "MyWorkflow.reference_genome": "hg38"
}
```

## How reference bundles work

Each workflow uses a file called `genome_resources.wdl` to manage reference genome paths. This file defines a `GenomeResources` struct containing paths to reference files (FASTA, index, dictionary, etc.) and a map that associates each genome type name (e.g., `"hg38"`, `"b37"`) with its resource paths.

The workflow's `genome_resources.wdl` contains a `GenomeResources` struct that defines the available resource fields. The struct looks like this:

```wdl
struct GenomeResources {
  File ref_fasta          # Reference FASTA file
  File ref_fasta_index    # Reference FASTA index (.fai)
  File ref_dict           # Sequence dictionary (.dict)
  File? ref_alt           # ALT contig file (optional)
  File? ua_index          # UA aligner index (optional)
  # ... additional workflow-specific fields
}
```

Fields marked as `File` (required) must be provided for every genome entry. Fields marked as `File?` (optional) can be omitted.

## Adding a Custom Reference Genome

To run the workflow with a reference genome that is not built in, edit the `genome_resources.wdl` file in the workflow directory and add your genome entry.

### Step 1: Prepare Your Reference Files

At minimum, you need the following files uploaded to an S3 bucket accessible from AWS HealthOmics:

- **Reference FASTA** (`.fasta` or `.fa`) - the genome sequence
- **FASTA index** (`.fasta.fai` or `.fa.fai`) - created with `samtools faidx`
- **Sequence dictionary** (`.dict`) - created with `samtools dict` or Picard `CreateSequenceDictionary`

Depending on the workflow, you may also need some of the following:

- **UA index** (`.uai`) - required for alignment workflows. See `howto-ua-align-sort.md` for instructions on building a UA index.
- **ALT contig file** (`.alt`) - for genomes with ALT contigs
- **Interval files** (`.bed`, `.interval_list`) - for targeted variant calling or coverage calculations

### Step 2: Locate genome_resources.wdl

In the workflow directory, find the `genome_resources.wdl` file under the `tasks/` subdirectory:

```
<workflow_name>/
    <workflow_name>.wdl
    tasks/
        genome_resources.wdl
        globals.wdl
        structs.wdl
        ...
```

### Step 3: Add Your Genome Entry

Open `genome_resources.wdl` and find the `GenomeResourcesWorkflow` workflow. It contains a `resources` output with a map of genome entries. Add a new entry for your custom genome inside the map.

For example, to add a genome called `my_custom_genome`:

```wdl
workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "hg38": {
        ...
      },
      "b37": {
        ...
      },
      "my_custom_genome": {
        "ref_fasta": "s3://my-bucket/my-ref/reference.fasta",
        "ref_fasta_index": "s3://my-bucket/my-ref/reference.fasta.fai",
        "ref_dict": "s3://my-bucket/my-ref/reference.dict",
        "ref_alt": "s3://my-bucket/my-ref/reference.fasta.alt",
        "ua_index": "s3://my-bucket/my-ref/reference.fasta.uai"
      }
    }
  }
}
```

Include all the fields that appear in the `GenomeResources` struct. Required fields (`File`) must be provided. Optional fields (`File?`) can be omitted from your entry.

### Step 4: Set the Workflow Input

In your workflow input JSON, set `reference_genome` to the name of your custom genome:

```json
{
    "MyWorkflow.reference_genome": "my_custom_genome"
}
```

## Important Notes

- All resource paths must be fully qualified S3 URIs (e.g., `s3://my-bucket/path/to/file`).
- Your S3 bucket must be accessible from the AWS HealthOmics service in the region where you run the workflow.
- Only include the resource fields that are present in the `GenomeResources` struct. Adding fields that are not in the struct will cause a validation error.
- If the workflow requires a UA index (`ua_index` field), you must build it for your custom reference before running the workflow. See `howto-ua-align-sort.md` for build instructions.
- The `reference_genome` input value must exactly match the key you added to the map (case-sensitive).

## Example: Adding a Custom Genome

This example walks through adding a custom genome and running the workflow with it, using TrimAlignSort as an illustration.

### 1. Upload reference files to S3

```bash
aws s3 cp my_reference.fasta s3://my-bucket/custom-ref/
aws s3 cp my_reference.fasta.fai s3://my-bucket/custom-ref/
aws s3 cp my_reference.dict s3://my-bucket/custom-ref/
aws s3 cp my_reference.fasta.uai s3://my-bucket/custom-ref/
```

### 2. Edit tasks/genome_resources.wdl

Add the following entry to the resources map:

```wdl
      "my_organism": {
        "ref_fasta": "s3://my-bucket/custom-ref/my_reference.fasta",
        "ref_fasta_index": "s3://my-bucket/custom-ref/my_reference.fasta.fai",
        "ref_dict": "s3://my-bucket/custom-ref/my_reference.dict",
        "ua_index": "s3://my-bucket/custom-ref/my_reference.fasta.uai"
      }
```

### 3. Set the input JSON

```json
{
    "TrimAlignSort.reference_genome": "my_organism",
    "TrimAlignSort.input_cram_bam_list": "s3://my-bucket/samples/sample1.cram"
}
```

### 4. Run the workflow

Deploy and run the workflow as usual through AWS HealthOmics. The workflow will use your custom reference files for alignment.
