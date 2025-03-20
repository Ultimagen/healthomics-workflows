# Single sample germline CNV Calling

CNV calling pipeline to detect CNVs using multiple tools.
It runs: 
1. cn.mops pipeline - using [cn.mops](https://github.com/Ultimagen/cn.mops)
2. cnvpytor pipeline - using [cnvpytor](https://github.com/abyzovlab/CNVpytor)
3. deletion candidates validation - using [jalign](https://github.com/Ultimagen/jalign)
4. results combination

Outputs: 
- cnmops_cnv_calls_bed : cn.mops cnv calls in bed format
- cnvpytor_cnv_calls_bed : cnvpytor cnv calls in tsv format
- combined_cnv_calls_bed : combined cnv calls in bed format
- combined_cnv_calls_bed_vcf : combined cnv calls in vcf format

## Workflow Overview
The pipeline takes a single sample aligned file (bam/cram) alongside with pre-calculated coverage in bedgraph format. 
In the cn.mops pipeline, the sample’s reads counts are added to a predefined cohort reads count matrix. The combined matrix is used as input for cn.mops algorithm to detect CNVs for the entire cohort. The sample’s CNVs are extracted and reported.
In the cnvpytor pipeline, copy number variation detection is performed from whole-genome sequencing read depth.
The CNV calls are then merged and the deletion CNV candidates are validated using jalign, that searches for alignments of the reads that support the CNV breakpoints. 
The CNV calls are returned in bed and vcf format. 

## Requirements

The workflow assume that you have the following files: 
1. BAM/CRAM - an aligned, sorted, duplicate marked UG BAM/CRAM file.
2. BedGraph - Previously calculated input bedGraph holding the coverage per base (outputs with the sequencing data).

### Files required for the analysis (download locally)
The following required files are publicly available:

    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
    
	s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v2.0/HapMap2_65samples_cohort_v2.0.hg38.ReadsCount.rds
    s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v2.0/HapMap2_65samples_cohort_v2.0.plus_female.ploidy
	s3://ultimagen-workflow-resources-us-east-1/hg38/Homo_sapiens_assembly38.chr1-24.w1000.bed
    s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v1.4/ug_cnv_lcr.bed

	
## Generating Germline CNV calls for a single sample

### Installation 
* for running cnmops, cnvpytor, and post processing use ugbio_cnv docker: <br>
	Pull **ugbio_cnv** docker image :
	```
	docker pull ultimagenomics/ugbio_cnv:1.8.0
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_cnv:1.8.0 /bin/bash
	```	
* for running jalign use jalign docker: <br>
    Pull **jalign** docker image :
	```
	docker pull ultimagenomics/jalign:1.0.0
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/jalign:1.0.0 /bin/bash
	```
    for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_sample_cnmops_CNV_calling/tasks/globals.wdl)<br><br>

* manual installation: 
if you would like to manually install the enviorment for UG-germline-CNV-calling you can follow the corresponding Dockerfiles:
[ugbio_cnv Dockerfile](https://github.com/Ultimagen/ugbio-utils/blob/main/src/cnv/Dockerfile) used to build the ugbio_cnv docker image.
[jalign Dockerfile](https://github.com/Ultimagen/jalign/blob/main/Dockerfile) used to build the jalign docker image.

## Run cn.mops

use ugbio_cnv docker

### cn.mops: Single Sample coverage collection for the case input format is BedGraph
```
bedtools map \
	-g Homo_sapiens_assembly38.fasta.fai \
	-a Homo_sapiens_assembly38.chr1-24.w1000.bed \
	-b {input_bed_graph} \
	-c 4 -o mean | \
	awk '{if($4=="."){print $1"\t"$2"\t"$3"\t"0}else{print $1"\t"$2"\t"$3"\t"$4}}' \
	> {sample_name}.win.bedGraph


Rscript --vanilla /src/cnv/cnmops/convert_bedGraph_to_Granges.R \
	-i {sample_name}.win.bedGraph \
	-sample_name {sample_name}
```

### cn.mops: Add sample coverage profile to the cohort
```
Rscript --vanilla /src/cnv/cnmops/merge_reads_count_sample_to_cohort.R \
	-cohort_rc HapMap2_65samples_cohort_v2.0.hg38.ReadsCount.rds \
	-sample_rc {sample_name}.ReadCounts.rds
```

#### cn.mops: Run cn.mops to call CNVs
```
Rscript --vanilla /src/cnv/cnmops/normalize_reads_count.R \
	--cohort_reads_count_file merged_cohort_reads_count.rds \
    --ploidy HapMap2_65samples_cohort_v2.0.plus_female.ploidy
	
Rscript --vanilla /src/cnv/cnmops/cnv_calling_using_cnmops.R \
	-cohort_rc cohort_reads_count.norm.rds \
	-minWidth 1000 \
	--save_csv
```

#### cn.mops: Filter sample CNVs
```
grep "{sample_name}" cohort.cnmops.cnvs.csv > {sample_name}.cnvs.csv
awk -F "," '{print $1"\t"$2-1"\t"$3"\t"$NF}' {sample_name}.cnvs.csv > {sample_name}.cnvs.bed;

filter_sample_cnvs \
	--input_bed_file {sample_name}.cnvs.bed \
	--intersection_cutoff 0.5 \
	--cnv_lcr_file ug_cnv_lcr.bed \
	--min_cnv_length 10000;
```

## Run cnvpytor
use ugbio_cnv docker
```
run_cnvpytor --input_bam_cram_file {input_bam} \
    --sample_name {sample_name} \
    --ref_fasta Homo_sapiens_assembly38.fasta \
    --bin_size 500 \
    --chr_list chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --out_directory .
```

## Run jalign
use jalign docker
```
echo "fetch cnmops DEL candidates"
cat {cnmops_cnvs_bed} | sed 's/UG-CNV-LCR//g' | sed 's/LEN//g' | sed 's/|//g' | grep -E "CN0|CN1" | \
    bedtools merge -c 4 -o distinct -d 1500 -i - \
    > {sample_name}.cnmops.DEL.merged.bed

echo "fetch cnvpytor DEL candidates"
cat {cnvpytor_cnvs_bed} | grep "deletion" | awk '{print $2"\t"$3"\t"$4"\t"$1","$5}' | \
    bedtools merge -c 4 -o distinct -d 1500 -i - \
    > {sample_name}.cnvpytor.DEL.merged.bed

echo "concat DEL candidates"
cat {sample_name}.cnmops.DEL.merged.bed {sample_name}.cnvpytor.DEL.merged.bed | \
    bedtools sort -i - \
    > {sample_name}.cnmops_cnvpytor.DEL.bed

echo "split DEL candidates for parallel processing"
mkdir cnmops500mod_cnvpytor500_DEL_split
cat {sample_name}.cnmops_cnvpytor.DEL.bed | split -l 10 -d --additional-suffix .bed
mv x*.bed cnmops500mod_cnvpytor500_DEL_split/

echo "Running jalign for DEL candidates"
python /jalign/scripts/parallel_run_cnv_realign.py \
    --folder_with_cnv_del_bed_files cnmops500mod_cnvpytor500_DEL_split \
    --input_cram {input_bam_file} \
    --ref_fasta Homo_sapiens_assembly38.fasta \
    --out_folder out_jalign \
    --sample_name {sample_name} \
    --num_jobs 36
```

## Combine results
use ugbio_cnv docker
```
combine_cnmops_cnvpytor_cnv_calls \
    --cnmops_cnv_calls {cnmops_cnvs_bed} \
    --cnvpytor_cnv_calls {cnvpytor_cnvs_tsv} \
    --del_jalign_merged_results {jalign_del_candidates} \
    --ug_cnv_lcr ug_cnv_lcr.bed \
    --fasta_index Homo_sapiens_assembly38.fasta.fai \
    --out_directory . \
    --sample_name {sample_name}
```
 output files:
- **bed file**: {sample_name}.cnmops_cnvpytor.cnvs.combined.bed
  	shows combined CNV calls called by cn.mops and cnvpytor.
- **annotated bed file**: {sample_name}.cnmops_cnvpytor.cnvs.combined.UG-CNV-LCR_annotate.bed
	shows combined CNV calls with UG-CNV-LCR annotation.
- **vcf file**: {sample_name}.cnv.vcf.gz
	shows combined CNV calls with UG-CNV-LCR annotation.
- **vcf index file**: {sample_name}.cnv.vcf.gz.tbi
	corresponding index to output vcf.