---
geometry: margin=1cm
output: pdf_document
---

# Realignment with UA and Sort

### System Requirements
    1.cpus : 32 
    2.Memory : 200 GB
    3.Platform : Intel Skylake

Publicly available reference files - download locally

* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta`
* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai`
* `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict`

Pre-built UA index file 

* `gs://concordanz/hg38/UA/Homo_sapiens_assembly38.fasta.uai`

### UA docker : 
    us-central1-docker.pkg.dev/ganymede-331016/ultimagen/ua:master_c553451

### step 1 : Build UA Index

### Command:

```
/ua/ua \
    --build \
    --ref ~{references.ref_fasta} \
    --seed 20,200,5 \
    --index ~{ua_index}.uai \
    --progress
```
### Step 2 : Align with UA

### Command:

```
samtools view -h -@ 32 ~{input_cram} -T ~{references.ref_fasta} -F 2048 | \
/ua/ua \
    --index ~{ua_index}.uai \
    --align true \
    --progress \
    --tp reference \
    --alt=~{ref_alt} \
    --json=~{output_bam_basename}.%s.json \
    --nthread max \
    --sam-input - \
    --sam-output - \
    --seed-score-ratio 0.5 --vector --huge --soft-clipping | \
samtools view -@ 32 -o ~{output_bam_basename}.bam -
```

### Step 3 : Sort UA Aligned BAM with Sorter 

### Sorter docker :
    us-central1-docker.pkg.dev/ganymede-331016/ultimagen/sorter:master_51c8e93
### Command: 

```
samtools view -h ~{output_bam_basename}.bam -@ 32 | 
demux \
    --input=- \
    --output-dir=demux_output/ \
    --runid=output \
    --nthreads 32 \
    --progress \
    --reference ~{references.ref_fasta} \
    --mark-duplicates=true \
    --basename=~{output_bam_basename} \
    --align=true

sorter \
    --runid=output \
    --input-dir=demux_output/ \
    --output-dir=. \
    --nthreads 32 \
    --progress \
    --timestamp= 
```