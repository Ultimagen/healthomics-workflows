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
* Install ugvc package:
    1. Clone the `VariantCalling` repository (e.g. to `software/VariantCalling`)
    2. Create a clean conda environment defined by `software/VariantCalling/setup/environment.yml`
    3. Create genomics.py3 conda environment:
    `conda env create -f software/VariantCalling/setup/environment.yml`
    This will create an environment called `genomics.py3`
    4. Create genomics.py3 conda environment:
    `conda env create -f software/VariantCalling/setup/other_envs/ucsc.yml`
    This will create an environment called `ucsc`
    5. Install ugvc package:
        ```
        conda activate genomics.py3
        cd software/VariantCalling
        pip install .
        ```
* Download UG-controlFREEC docker image from: 
    `us-central1-docker.pkg.dev/ganymede-331016/ultimagen/ug_control_freec:1679a9`

### Create mpileup file for tumor and normal seperatly. 
```
conda activate genomics.py3
samtools mpileup -f Homo_sapiens_assembly38.fasta \
    -d 8000 \
    -Q 0 \
    -q 1 \
    -l out.vcf.gz \
    input.bam \
    >  {sample_name}_minipileup.pileup
```

### Collect coverage for tumor and normal seperatly
```
conda activate genomics.py3

COVERAGE_ANALYSIS="coverage_analysis.py collect_coverage"
OUTPUT=coverage
mkdir $OUTPUT

$COVERAGE_ANALYSIS \
	-i {input_cram_bam} \
	-o $OUTPUT \
	-Q 1 \
	--reference Homo_sapiens_assembly38.fasta
```

### Convert BigWig format to CPN format for tumor and normal seperatly
```
conda activate ucsc

for file in coverage/*.depth.bw; do \
	bigWigToBedGraph $file /dev/stdout | \
	bedtools map -g Homo_sapiens_assembly38.fasta.fai \
	-a Homo_sapiens_assembly38.w1000.chr1-23.bed \
	-b /dev/stdin \
	-c 4 -o mean | \
	awk '{if($4=="."){print $1"\t"$2"\t"0}else{print $1"\t"$2"\t"$4}}' | \
	grep -v "chrY" | \
	sed 's/^chr//' > \
	$file.cpn; \
	done
	
cat *.cpn > {sample_name}.cpn
```

### in case input format is BedGraph, convert BedGraph to CPN format for tumor and normal seperatly. 
```
conda activate genomics.py3

gzip -d {sample_name}.BedGraph.gz
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
run the ug_control_freec docker image in interactive mode and mount the volume with your data. 
for example: 
`docker run -it -v /data:/data us-central1-docker.pkg.dev/ganymede-331016/ultimagen/ug_control_freec:1679a9`

```
#split reference to file per chromosome
mkdir chrFiles_dir
cd chrFiles_dir
faidx -x ../Homo_sapiens_assembly38.fasta
cd ../

#create controlFREEC config file
python /generate_controlFREEC_config.py \
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
```
conda activate genomics.py3

#convert to bedfile
cat {tumor}_CNVs | sed 's/^/chr/' | cut -f1-4 > {tumor}.cnvs.bed

#annotate cnvs bed file
python /VariantCalling/ugvc filter_sample_cnvs \
	--input_bed_file {tumor}.cnvs.bed \
	--intersection_cutoff 0.5 \
	--cnv_lcr_file ug_cnv_lcr.bed \
	--min_cnv_length 10000
```