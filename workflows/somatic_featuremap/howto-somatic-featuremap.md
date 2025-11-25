# SomaticFeaturemap Workflow

## Overview

The SomaticFeaturemap workflow generates a merged somatic featuremap VCF file from tumor and normal featuremap VCF files. This workflow is designed to combine information from tumor-normal pairs and enrich the data with additional features for somatic variant analysis.

## Description

Given tumor and normal featuremap VCF files, this pipeline:

1. **Merges** tumor and normal featuremap VCF files while preserving all tumor sample original information.
2. **Enriches** the merged file with additional information:
   - Tandem repeat information: proximity of variants to closest tandem repeats and their details.
   - Ref/non-ref counts: counts in positions around variant loci for both tumor and normal samples.
3. **Scores variants** using a pre-trained XGBoost model to infer probability of candidates being true variants.
4. **Outputs** the merged somatic featuremap VCF with added information and probability scores.

## Requirements

### Input Files Required

#### Required Files
- **Tumor featuremap VCF** and its index file 
- **Normal featuremap VCF** and its index file
- **Tumor CRAM/BAM file** and its index file
- **Normal CRAM/BAM file** and its index file
- **Reference genome** files (FASTA, FAI, and DICT)
- **Interval list** file defining genomic regions to process
- **Tandem repeats reference file** (tsv format)

#### Optional Files
- **XGBoost model file** for variant scoring (if not provided, default model is used)
default model: 
```
gs://concordanz/somatic_featuremap/sfm.xgb_model.V1.6.json
```
- **Pre-computed somatic featuremap** (if available, skips merging step)

### Reference Files

The following reference files are publicly available:

```
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
gs://concordanz/hg38/wgs_calling_regions.without_encode_blacklist.hg38.interval_list
gs://concordanz/hg38/tandem_repeats.hg38.bed
```




### Installation
* Use docker: <br>
	Pull ugbio_featuremap, ugbio_core_docker, ugbio_freec_docker docker images :
	```
	docker pull ultimagenomics/ugbio_featuremap:<release-tag>
	docker pull ultimagenomics/ugbio_freec:<release-tag>
    docker pull ultimagenomics/ugbio_core:<release-tag>
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data <docker-image> /bin/bash
	```
	for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/902c0def79e17c71ef810f7cdd887e06e736c5b4/workflows/single_read_snv/tasks/globals.wdl)<br>

## Workflow Steps

1. **Featuremap Merging** 
Combines tumor and normal featuremap VCFs

from inside ugbio_featuremap docker image: 
```
create_somatic_featuremap \
    --tumor_vcf  {tumor_featuremap_vcf} \
    --normal_vcf {normal_featuremap_vcf} \
    --sample_name {sample_name} \
    --out_directory . 
```
2. **Variant Padding** 
Creates padded regions around somatic featuremap (tumor-PASS) variants
from inside ugbio_core_docker image: 
```
    # Extract CHROM, POS0, END, REF, ALT and calculate the maximum length
    bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' {input_vcf} | awk -F'\t' 'BEGIN{OFS="\t"}{
        s=$2;
        m=length($4);
        n=split($5,a,",");
        for(i=1;i<=n;i++) if(length(a[i])>m) m=length(a[i]);
        print $1, s, s+m
    }' | gzip > variants.bed.gz
    
    # Create a genome file for bedtools slop (chromosome sizes)
    cut -f1,2 Homo_sapiens_assembly38.fasta.fai > genome.txt
    
    # Pad the bed file using bedtools slop"
    zcat variants.bed.gz | bedtools slop -i stdin -g genome.txt -b {pad_size} | gzip > sfm.padded.bed.gz
   ``` 
3. **Mpileup Generation** 
Generates pileup information for tumor and normal samples separately. 
from inside ugbio_freec docker image: 
```
gunzip sfm.padded.bed.gz
samtools mpileup -f Homo_sapiens_assembly38.fasta \
    -d 8000 \
    -Q 0 \
    -q 1 \
    -l sfm.padded.bed \
    {sample_name}.bam \
    >  {sample_name}_minipileup.pileup
```

4. **Feature Integration** 
Integrates mpileup data into the somatic featuremap
from inside ugbio_core_docker image: 
```
integrate_mpileup_to_sfm \
    --sfm_vcf {somatic_featuremap_vcf} \
    --tumor_mpileup  {tumor_mpileup} \
    --normal_mpileup  {normal_mpileup} \
    --distance_start_to_center 2 \
    --out_directory .
```

6. **Transformation and Variant Scoring** 
Transforms data to features and applies XGBoost model for variant probability scoring
from inside ugbio_core_docker image: 
```
#convert interval list to bed file: 
cat wgs_calling_regions.without_encode_blacklist.hg38.interval_list | grep -v "^@" | awk '{print $1"\t"$2"\t"$3}' > interval_list.bed

somatic_featuremap_fields_transformation \
    -sfm {somatic_featuremap_vcf} \
    -o ~{out_somatic_featuremap_vcf} \
    -i interval_list.bed \
    -ref_tr ~{ref_tandem_repeats} \
    -xgb_model sfm.xgb_model.V1.6.json
```

## Output Files

The workflow generates the following output files:

### Primary Outputs
- **`final_out_sfm_vcf`** - Final somatic featuremap VCF with all enhancements and scores
- **`final_out_sfm_vcf_index`** - Index file for the final output VCF

### Conditional Outputs (when generated)
- **`somatic_featuremap`** - Initial merged somatic featuremap VCF
- **`somatic_featuremap_index`** - Index for initial merged VCF
- **`tumor_mpileup`** - Tumor mpileup file
- **`normal_mpileup`** - Normal mpileup file
- **`somatic_featuremap_mpileup_vcf`** - Somatic featuremap with integrated mpileup info
- **`somatic_featuremap_padded_bed_file`** - Padded BED file for variant regions
