# Somatic CNV Calling

Somatic CNV calling pipeline detects CNVs for an paired tumor-normal sample 
based on coverage changes and B-allelic frequencies using the [controlFREEC](https://boevalab.inf.ethz.ch/FREEC/).

## Workflow Overview
The pipeline takes the aligned tumor and germline files (bam/cram) 
and calculates the coverage profile along the genome and 
allele frequencies of predefined locations. 
Finally, it runs controlFREEC to call CNVs and filters them by length and low confidence regions.

## Requirements

The workflow assumes that you have coverage representation of the tumor-normal sample in one of the following formats: 
1. BAM/CRAM - an aligned, sorted, duplicate marked UG BAM/CRAM file, seperatly for tumor and normal. 
2. BedGraph - previously calculated input bedGraph holding the coverage per base (outputs with the sequencing data), seperatly for tumor and normal.
3. CPN - pre-calculated binned coverage for the normal sample in the format needed for controlFREEC (cpn), seperatly for tumor and normal.

The workflow assumes that you have allele frequency representation of the tumor-normal sample in one of the following formats: 
1. BAM/CRAM - an aligned, sorted, duplicate marked UG BAM/CRAM file, seperatly for tumor and normal. 
2. MPILEUP - previously calculated mpileup, seperatly for tumor and normal.

### Files required for the analysis (download locally)
The following files are publicly available

    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
    gs://concordanz/hg38/af-only-gnomad.hg38.AF_gt0.35.CHR1-24.vcf.gz
    gs://concordanz/hg38/af-only-gnomad.hg38.AF_gt0.35.CHR1-24.vcf.gz.tbi
    gs://concordanz/hg38/Homo_sapiens_assembly38.w1000.chr1-23.bed
    gs://concordanz/hg38/out100m2_hg38.gem
    gs://concordanz/hg38/UG-High-Confidence-Regions/v1.4/ug_cnv_lcr.bed

## Generating Germline CNV calls for a single sample

### Installation
* Use docker: <br>
	Pull ugbio_freec and ugbio_cnv docker images :
	```
	docker pull ultimagenomics/ugbio_freec:1.5.5
	docker pull ultimagenomics/ugbio_cnv:1.5.5
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_freec:1.5.5 /bin/bash
	docker run -it -v /data:/data ultimagenomics/ugbio_cnv:1.5.5 /bin/bash
	```
	for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/902c0def79e17c71ef810f7cdd887e06e736c5b4/workflows/single_read_snv/tasks/globals.wdl#L68C31-L68C64)<br>
* manual installation: 
if you would like to manually install the enviorment for UG-germline-CNV-calling you can follow the following Dockerfiles:
    1. [ugbio_freec Dockerfile](https://github.com/Ultimagen/FREEC/blob/master/Dockerfile) used to build the ugbio_freec docker image.
    2. [ugbio_cnv Dockerfile](https://github.com/Ultimagen/ugbio-utils/blob/main/src/cnv/Dockerfile) used to build the ugbio_cnv docker image.

### Create mpileup file for tumor and normal seperatly. 
from inside ugbio_freec docker image: 
```
samtools mpileup -f Homo_sapiens_assembly38.fasta \
    -d 8000 \
    -Q 0 \
    -q 1 \
    -l out.vcf.gz \
    input.bam \
    >  {sample_name}_minipileup.pileup
```

### Collect coverage for tumor and normal seperatly
from inside ugbio_freec docker image: 
```
samtools depth \
	-J \
	-Q 1 \
	--reference Homo_sapiens_assembly38.fasta \
	{input_cram_bam} | \
	awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' > {sample_name}.bedgraph
```

### convert BedGraph to CPN format for tumor and normal seperatly. 
from inside ugbio_freec docker image: 
```
#unzip input bedgraph file if needed:
if [[ $bedgraph =~ \.gz$ ]];
then
	gzip -d -c $bedgraph > {sample_name}.bedgraph;

bedtools map -g Homo_sapiens_assembly38.fasta.fai \
	-a Homo_sapiens_assembly38.w1000.chr1-23.bed \
	-b {sample_name}.BedGraph \
	-c 4 -o mean | \
	awk '{if($4=="."){print $1"\t"$2"\t"0}else{print $1"\t"$2"\t"$4}}' | \
	grep -v "chrY" | \
	sed 's/^chr//' > \
	{sample_name}.cpn
```

### runControlFREEC
from inside ugbio_freec docker image: 
```
#split reference to file per chromosome
mkdir chrFiles_dir
cd chrFiles_dir
faidx -x ../Homo_sapiens_assembly38.fasta
cd ../

#create controlFREEC config file
generate_controlFREEC_config \
	--sample_name {sample_name} \
	--BedGraphOutput TRUE \
	--chrLenFile Homo_sapiens_assembly38.fasta.fai \
	--contaminationAdjustment TRUE \
	--maxThreads {maxThreads} \
	--window 1000 \
	--chrFiles chrFiles_dir \
	--degree 3 \
	--gemMappabilityFile out100m2_hg38.gem \
	--sample_mateFile {tumor}.bam \
	--sample_mateCopyNumberFile {tumor}.cpn \
	--sample_miniPileupFile {tumor}.mpileup \
	--sample_inputFormat BAM \
	--sample_mateOrientation 0 \
	--control_mateFile {normal}.bam \
	--control_mateCopyNumberFile {normal}.cpn \
	--control_miniPileupFile {normal}.mpileup \
	--control_inputFormat BAM \
	--control_mateOrientation 0 \
	--baf_makePileup af-only-gnomad.hg38.AF_gt0.35.CHR1-24.vcf.gz \
	--baf_fastaFile Homo_sapiens_assembly38.fasta \
	--baf_SNPfile af-only-gnomad.hg38.AF_gt0.35.CHR1-24.vcf.gz

#run controlFREEC
/freec -conf {sample_name}.config.txt
```

#### Filter ControlFREEC called CNVs by length and low confidence regions
from inside ugbio_cnv docker image:
```
#convert to bedfile
cat {tumor}_CNVs | sed 's/^/chr/' | cut -f1-4 > {tumor}.cnvs.bed

#annotate cnvs bed file
filter_sample_cnvs \
	--input_bed_file {tumor}.cnvs.bed \
	--intersection_cutoff 0.5 \
	--cnv_lcr_file ug_cnv_lcr.bed \
	--min_cnv_length 10000
```