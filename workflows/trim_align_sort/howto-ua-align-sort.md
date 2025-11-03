---
geometry: margin=1cm
output: pdf_document
---

# Realignment with UA and Sort

## Overview

The complete workflow for UA alignment and sorting consists of three main steps:

1. **Build UA Index** (if not already available): Create the UA index file from the reference genome. This part only needs to be done once per reference genome.
2. **Align with UA**: Align reads to the reference genome using the UA aligner.
3. **Process with Demux and Sorter**: Process the aligned BAM file using demux and sorter tools.

This workflow provides high-performance alignment and sorting capabilities, with options for variant-aware alignment, duplicate marking, UMI handling, and more. The commands provided in this document can be used to run these tools manually, outside of the WDL pipeline environment. The steps mirror some of the functionality of the workflow workflows/trim_align_sort/trim_align_sort.wdl in the UG HealthOmics repo (https://github.com/Ultimagen/healthomics-workflows).


## Docker Images

The most up-to-date docker images can be found in the `workflows/trim_align_sort/tasks/globals.wdl` file in the repository. As of the latest update, the following docker images are used:

### Prefix
    ultimagenomics/ (hosted on DockerHub)

### UA docker:
    ultimagenomics/alignment:3.0.4

### Sorter docker:
    ultimagenomics/sorter:1.4.15

## System Requirements
    1. CPUs: 32-40 
    2. Memory: 
       1. UA: UA_index_size + 10 GiB (60GiB for hg38)
       2. Demux: 16GiB
       3. Sorter: 64GiB (likely requires tuning)
    3. Platform : Intel Skylake

### Memory Requirements

Note that memory requirements can vary significantly based on the input data:

- UA alignment memory requirements depend on the UA index size.
- Demux memory requirements are generally fixed.
- Sorter memory requirements are data-dependent, especially when marking duplicates, even more so when UMIs are used. The UG WDL pipeline calculates the required memory based on the demux output, this part is left out for simplicity, but note that memory might need to be properly calibrated. When memory is insufficient, Sorter will raise an error immediately.

Always monitor memory usage and adjust resources accordingly for your specific dataset.


## Reference Files

The following reference files are needed for UA alignment and sorting. These files are specified in the `wdls/input_templates/trim_align_sort_template.json` file:

### Publicly available reference files - download locally

* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta` - Reference FASTA
* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai` - Reference FASTA index
* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict` - Reference dictionary
* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt` - Reference ALT file (optional)

### Pre-built UA index file 

* `gs://concordanz/hg38/UA/Homo_sapiens_assembly38.fasta.uai` - Pre-built UA index


### Step 1 : Build UA Index
This step creates the UA index file needed for alignment. Note that this file is typically very large, e.g. ~51GiB for the standard hg38 reference genome. This step only needs to be done once per reference genome, and a pre-built index file for hg38 is already provided at `gs://concordanz/hg38/UA/Homo_sapiens_assembly38.fasta.uai`.

This command is implemented in the `BuildUaIndex` task in `wdls/tasks/alignment_tasks.wdl`.

### Command:

```
ua \
    --build \
    --ref /path/to/reference.fasta \
    --seed 20,200,5 \
    --index /path/to/output.uai \
    --progress
```

### Step 2 : Align with UA
This step aligns reads to the reference genome using the UA aligner.

This command is implemented in the `AlignWithUA` task in `wdls/tasks/alignment_tasks.wdl`.

### Command:

```
# When working with a single input BAM/CRAM file
samtools view -h -@ 32 /path/to/input.cram -T /path/to/reference.fasta -F 2048 | \
ua \
    --index /path/to/reference.uai \
    --align true \
    --progress \
    --tp reference \
    --alt=/path/to/ref.alt \
    --stat=output_basename.%s.json \
    --nthread max \
    --sam-input - \
    --sam-output - \
    --seed-score-ratio 0.5 --vector --soft-clipping | \
samtools view -@ 32 -o output_basename.bam -

# When working with multiple input BAM files
samtools merge -@ 32 -c -O SAM /dev/stdout input1.bam input2.bam input3.bam | \
ua \
    --index /path/to/reference.uai \
    --align true \
    --progress \
    --tp reference \
    --alt=/path/to/ref.alt \
    --stat=output_basename.%s.json \
    --nthread max \
    --sam-input - \
    --sam-output - | \
samtools view -@ 32 -o output_basename.bam -
```

#### Key Parameters:
- `--index`: Path to the UA index file
- `--tp reference`: Sets the tp/t0 tags direction to reference
- `--alt`: Path to the reference alt file (optional, can be removed if not needed)
- `--stat`: Output statistics file pattern
- `--nthread max`: Use maximum available threads

### Step 3 : Sort UA Aligned BAM with Demux and Sorter
This step processes the aligned BAM file using demux and sorter tools.

### Demux Command:
Demux processes the aligned BAM file and prepares it for sorting.

This command is implemented in the `Demux` task in `wdls/tasks/sorting_tasks.wdl`.

```
# Basic usage
samtools view -h /path/to/aligned.bam -@ 32 | \
demux \
    --input=- \
    --output-dir=demux_output/ \
    --runid=output_basename \
    --nthreads 32 \
    --progress \
    --reference /path/to/reference.fasta \
    --mark-duplicates=true \
    --output-group=output_basename \
    --output-path={outputGroup}/{outputGroup} \
    --align=true


# With UMI tags
samtools view -h /path/to/aligned.bam -@ 32 | \
demux \
    --input=- \
    --output-dir=demux_output/ \
    --runid=output_basename \
    --nthreads 32 \
    --progress \
    --reference /path/to/reference.fasta \
    --mark-duplicates=true \
    --umi=RX \
    --output-group=output_basename \
    --output-path={outputGroup}/{outputGroup} \
    --align=true
```

#### Key Parameters:
- `--input`: Input BAM/SAM file or stream
- `--output-dir`: Directory for demux output
- `--runid`: Run identifier
- `--nthreads`: Number of threads to use
- `--reference`: Path to reference FASTA file
- `--mark-duplicates`: Whether to mark duplicates (true/false)
- `--output-group`: Output group name
- `--output-path`: Output path pattern
- `--align`: Whether the input is aligned (true/false)
- `--umi`: UMI tag name, comma-separated if there are multiple UMIs (e.g., RX or RX,RY)

### Sorter Command:
Sorter processes the output from demux to produce the final sorted output. It organizes reads and generates statistics.

This command is implemented in the `Sorter` task in `wdls/tasks/sorting_tasks.wdl`.

```
sorter \
    --runid=output_basename \
    --input-dir=demux_output/ \
    --output-dir=sorter_output/ \
    --nthreads 32 \
    --progress \
    --timestamp=000  # Used to override sorter output timestamp

```

#### Key Parameters:
- `--runid`: Run identifier
- `--input-dir`: Directory containing demux output
- `--output-dir`: Directory for sorter output
- `--nthreads`: Number of threads to use
- `--timestamp`: Timestamp for output organization



## Output Structure:
After running sorter, the output will be organized as follows:
```
sorter_output/
    output_basename-000/  # timestamp is appended to the run ID
        output_basename/output_basename.cram  # Main output CRAM file
        output_basename/output_basename.cram.crai  # Index file
        output_basename/output_basename.csv  # Statistics in CSV format
        output_basename/output_basename.json  # Statistics in JSON format
```
The other files that were output from UA and demux (demux_output/ + output_basename.bam) can be deleted once the sorter outputs are created.