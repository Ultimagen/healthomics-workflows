# Segmental Duplication Analysis

Segmental duplication analysis pipeline uses parascopy tool to analyze copy number variations and small variants in genomic regions with segmental duplications. The pipeline takes a CRAM/BAM file as input and produces remapped CRAM file, CNV calls, and small variant calls.

## Workflow Overview

The pipeline takes as input a CRAM/BAM file and outputs remapped reads, CNV calls, and small variants specific to segmental duplication regions.
The steps of the pipeline are as follows: pool reads from segmental duplication regions, call copy number variations, and call small variants using DeepVariant preprocessing.

## Requirements

1. Input CRAM/BAM files and respective indexes
2. Pull parascopy Docker image
3. Install samtools

### Files required for the analysis (download locally)

```bash
# Homology table and index
s3://ultimagen-workflow-resources-us-east-1/hg38/segmental_duplications/hg38.bed.gz
s3://ultimagen-workflow-resources-us-east-1/hg38/segmental_duplications/hg38.bed.gz.tbi

# Segmental duplication regions
s3://ultimagen-workflow-resources-us-east-1/hg38/segmental_duplications/hom.regions.bed

# Background regions for CNV calling
s3://ultimagen-workflow-resources-us-east-1/hg38/segmental_duplications/background_regions/hg38.bg.bed.gz

# Trained CNV model
s3://ultimagen-workflow-resources-us-east-1/hg38/segmental_duplications/parascopy_model_250205/model.tar.gz

# DeepVariant model for segmental duplications
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/segdup_1.0/model_dyn_1500.onnx

# Reference genome files
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
```

#### Step 1: Pool reads from segmental duplication regions

Run the following command to remap reads from each segmental duplication region:

```bash
# Pull the Docker image
docker pull ultimagenomics/parascopy:1.2.0_f42c9e4

# Process each region in the segdup_regions BED file
while IFS= read -r line; do 
    l1=$(echo "$line" | awk '{print $1":"$2+1"-"$3}')
    
    docker run --rm -v $(pwd):/data ultimagenomics/parascopy:1.2.0_f42c9e4 \
    parascopy pool -i /data/input.cram \
                   -t /data/hg38.bed.gz \
                   -f /data/Homo_sapiens_assembly38.fasta \
                   -o /data/sample."$l1".remap.bam \
                   -m 0 \
                   --tags_to_reverse t0 tp \
                   --tags_to_retain XA XB \
                   -r "$l1"
done < hom.regions.bed

# Merge all remapped BAM files
find . -name "sample.*:*.bam" > file.lst
samtools cat -o sample.remap.tmp.bam -b file.lst

# Sort and convert to CRAM
samtools sort --reference Homo_sapiens_assembly38.fasta \
              -o sample.remap.cram sample.remap.tmp.bam \
              --output-fmt-option embed_ref=1
samtools index sample.remap.cram
```

#### Step 2: Call copy number variations

Extract and use the trained CNV model to call copy numbers:

```bash
# Extract the CNV model
tar --no-same-owner --no-same-permissions -xvf model.tar.gz

# Calculate depth coverage
docker run --rm -v $(pwd):/data ultimagenomics/parascopy:1.2.0_f42c9e4 \
parascopy depth --input /data/input.cram \
                --bed-regions /data/hg38.bg.bed.gz \
                --fasta-ref /data/Homo_sapiens_assembly38.fasta \
                -o /data/sample.depth -@20 \
                --clipped-perc 20 --unpaired-perc 120

# Call copy numbers using the trained model
docker run --rm -v $(pwd):/data ultimagenomics/parascopy:1.2.0_f42c9e4 \
parascopy "cn-using" model \
          --input /data/input.cram \
          --fasta-ref /data/Homo_sapiens_assembly38.fasta \
          -d /data/sample.depth \
          -t /data/hg38.bed.gz \
          -o /data/sample.cn -@20
```

#### Step 3: Call small variants using DeepVariant

Before running parascopy call, you need to generate precalled variants using DeepVariant with the specialized segmental duplication model. See [howto-germline-calling-efficient-dv.md](howto-germline-calling-efficient-dv.md) for detailed DeepVariant instructions.

For segmental duplication analysis, use the specialized model:
`s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/segdup_1.0/model_dyn_1500.onnx`

After running DeepVariant to generate precalled variants, run parascopy call:

```bash
# Call small variants using copy number results and DeepVariant preprocessing
docker run --rm -v $(pwd):/data ultimagenomics/parascopy:1.2.0_f42c9e4 \
parascopy call -p /data/sample.cn \
               -i /data/input.cram \
               --fasta-ref /data/Homo_sapiens_assembly38.fasta \
               -t /data/hg38.bed.gz \
               --freebayes /usr/local/bin/freebayes \
               --precalled-variants /data/deepvariant_variants.vcf.gz \
               -o /data/sample.calls -@20
```

## Outputs

The pipeline produces the following outputs:

1. **sample.remap.cram** - Remapped CRAM file with reads from segmental duplication regions
2. **sample.remap.cram.crai** - Index for the remapped CRAM file  
3. **sample.cn/res.samples.bed.gz** - Copy number variation calls in BED format
4. **sample.cn/res.samples.bed.gz.tbi** - Index for CNV calls
5. **sample.calls/variants.vcf.gz** - Small variant calls in VCF format
6. **sample.calls/variants.vcf.gz.tbi** - Index for small variant calls

The output files contain:
- **Remapped CRAM**: Reads from segmental duplication regions collapsed and remapped to representative sequences
- **CNV calls**: Copy number estimates for segmental duplication regions
- **Small variants**: SNPs and indels called specifically in segmental duplication regions