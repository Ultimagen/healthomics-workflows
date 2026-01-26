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
The CNV calls are then merged and CNV candidates are validated using jalign, that searches for alignments of the reads that support the CNV breakpoints. 
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
    s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v2.1.2/ug_cnv_lcr.bed

**Note:** For male samples use 
    s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v2.0/HapMap2_65samples_cohort_v2.0.plus_male.ploidy
    The ploidy file contains single number per row that indicates the expected ploidy of the X chromosome. The sample being 
    called corresponds to the **last** number in the ploidy file. 
	
ML filtering model can be found here: 
    
    s3://ultimagen-workflow-resources-us-east-1/filtering_models/cnv_filtering_model_v1.5.0.pkl

## Generating Germline CNV calls for a single sample

In the pipeline we run two CNV callers: cn.mops and CNVPytor, combine their results and validate candidates using our tool `jalign` - that checks if there are reads that can be realigned accross the breakpoint. 

Finally, we apply simple ML model to filter candidates according to the length, copy number, the tool that proposed the candidate and the jump reads support.

### Installation 
* for running cnmops, cnvpytor, and post processing use ugbio_cnv docker: <br>
	Pull **ugbio_cnv** docker image :
	```
	docker pull ultimagenomics/ugbio_cnv:1.18.0
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_cnv:1.18.0 /bin/bash
	```	
for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/germline_CNV_pipeline/tasks/globals.wdl)

* for running jalign use jalign docker: <br>
    Pull **jalign** docker image :
	```
	docker pull ultimagenomics/jalign:1.1.0
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/jalign:1.1.0 /bin/bash
	```
for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/germline_CNV_pipeline/tasks/globals.wdl)

* For running ML-based filtering - use **ugbio_filtering** docker image. 
	```
	docker pull ultimagenomics/ugbio_filtering:1.18.0
	```
	Run docker in interactive mode: 
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_filtering:1.18.0 /bin/bash
	```	


### Run cn.mops

use ugbio_cnv docker

#### cn.mops: Single Sample coverage collection for the case input format is BedGraph
```
bedtools map \
	-g Homo_sapiens_assembly38.fasta.fai \
	-a Homo_sapiens_assembly38.chr1-24.w1000.bed \
	-b {input_bed_graph} \
	-c 4 -o mean |\
	awk '{if($4=="."){print $1"\t"$2"\t"$3"\t"0}else{print $1"\t"$2"\t"$3"\t"$4}}' \
	> {sample_name}.win.bedGraph


Rscript --vanilla /home/ugbio/src/cnv/cnmops/convert_bedGraph_to_Granges.R \
	-i {sample_name}.win.bedGraph \
	-sample_name {sample_name}
```

#### cn.mops: Add sample coverage profile to the cohort
```
Rscript --vanilla /home/ugbio/src/cnv/cnmops/merge_reads_count_sample_to_cohort.R \
	-cohort_rc HapMap2_65samples_cohort_v2.0.hg38.ReadsCount.rds \
	-sample_rc {sample_name}.ReadCounts.rds \
	--save_hdf
```

#### cn.mops: Run cn.mops to call CNVs
```
# normalize coverage (takes care of ploidies differences)
Rscript --vanilla /home/ugbio/src/cnv/cnmops/normalize_reads_count.R \
     --cohort_reads_count_file merged_cohort_reads_count.rds \
     --ploidy {ploidy_file} \
     --chrX_name chrX \
     --chrY_name chrY \
     --cap_coverage

Rscript --vanilla /home/ugbio/src/cnv/cnmops/cnv_calling_using_cnmops.R \
	-cohort_rc cohort_reads_count.norm.rds \
	-minWidth 2 \
	-p {parallel} \
	--save_csv \
	--mod_cnv
```

#### cn.mops: Extract normalized read counts

This step is done to annotate the output VCF with the cohort mean and the sample mean/std coverage. 

```
# Extract normalized read counts
Rscript --vanilla /home/ugbio/src/cnv/cnmops/export_cohort_matrix_to_bed.R \
    cohort_reads_count.norm.rds --sample {sample_name}
Rscript --vanilla /home/ugbio/src/cnv/cnmops/export_cohort_matrix_to_bed.R \
    cohort_reads_count.norm.rds --mean
```

These scripts generate files: `{sample_name}.cov.bed` and `coverage.cohort.bed`

#### cn.mops: Process CNV outputs of cn.mops - convert them to the CNV VCF file

```
grep "{sample_name}" cohort.cnmops.cnvs.csv > {sample_name}.cnvs.csv
awk -F "," '{print $1"\t"$2-1"\t"$3"\t"$NF}' {sample_name}.cnvs.csv > {sample_name}.cnvs.bed

process_cnvs \
	--sample_name {sample_name} \
	--input_bed_file {sample_name}.cnvs.bed \
	--intersection_cutoff 0.5 \
	--cnv_lcr_file ug_cnv_lcr.bed \
	--min_cnv_length 10000 \
	--out_directory . \
	--sample_norm_coverage_file {sample_name}.cov.bed \
	--cohort_avg_coverage_file coverage.cohort.bed \
	--fasta_index_file Homo_sapiens_assembly38.fasta.fai
```

This produces the result: `{sample_name}.cnvs.annotate.vcf.gz`

### Run cnvpytor
To improve robustness of the results it is recommended to run CNVPytor with two window sizes: 500 and 2500 bases. 

This is an example of one run, should be run twice. Use `ugbio_cnv` docker
```
        cnvpytor -root {sample_name}.pytor \
            -rd {input_bam} \
            -chrom chr1 chr2 .... chrY \
            -T Homo_sapiens_assembly38.fasta \ 
        
        cnvpytor -root {sample_name}.pytor \
            -his <window_size>

        cnvpytor -root {sample_name}.pytor \
            -partition <window_size>
        
		## output VCF
        cnvpytor -root {sample_name}.pytor  -view <window_size> <<EOF
                set print_filename {sample_name}.<window_size>.CNV.vcf
                print calls
                quit 
        EOF
```

The callsets from the two window sizes are combined using `ugbio_cnv` docker. 
```
        combine_cnmops_cnvpytor_cnv_calls concat \
            --cnvpytor_vcf vcf1 vcf2 \
            --output_vcf {sample_name}.cnvpytor.cnv_calls.vcf.gz \
            --fasta_index Homo_sapiens_assembly38.fasta.fai
```

### Combine, annotate and filter cn.mops and cnvpytor callsets

#### Annotate candidates: 
The following annotations are applied: 
1. cn.mops duplications called shorter than 10Kb are marked as SHORT_CNMOPS_DUPLICATION
2. Close cn.mops duplication calls (distance less than 1,500 bp) are combined
3. Percentage of N bases in the CNV is recorded
4. UG-CNV-LCR regions annotations are added 
   
   ```
    # Concatenate cn.mops and cnvpytor VCF files
	combine_cnmops_cnvpytor_cnv_calls concat \
		--cnmops_vcf {cnmops_vcf} \
		--cnvpytor_vcf {cnvpytor_vcf} \
		--output_vcf {sample_name}.step1.vcf.gz \
		--fasta_index Homo_sapiens_assembly38.fasta.fai

	# Filter short cn.mops duplications, merge adjacent ones
	combine_cnmops_cnvpytor_cnv_calls filter_cnmops_dups \
		--combined_calls {sample_name}.step1.vcf.gz \
		--combined_calls_annotated {sample_name}.step2.vcf.gz \
		--filtered_length 10000 \
		--distance_threshold 1500

	# Annotate with gaps information
	combine_cnmops_cnvpytor_cnv_calls annotate_gaps \
		--calls_vcf {sample_name}.step2.vcf.gz \
		--output_vcf {sample_name}.step3.vcf.gz \
		--ref_fasta Homo_sapiens_assembly38.fasta

	
#### Add annotations of the split read support 

```
analyze_cnv_breakpoint_reads \
    --bam-file {input_cram_bam_file} \
    --vcf-file {sample_name}.combined.vcf.gz \
    --reference-fasta Homo_sapiens_assembly38.fasta \
    --cushion 1500 \
    --output-file {sample_name}.split.annotated.vcf.gz
```

#### Run `jalign`

This task re-aligns reads that are close to the CNV breakpoints to check how many reads can support the breakpoint. For deletion, it looks for the reads
that span the whole CNV. For the duplication, this looks for the reads that can align with the pattern consistent with the duplication.

For this task - use `jalign` docker. Note that this tasks requires a 
large machine (~32 CPU recommended), which is expected to take about 1-2 hours. 

```
run_jalign --tool-path /opt/para_jalign \
    {input_bam_file} \
    {sample_name}.split.annotated.vcf.gz \
    Homo_sapiens_assembly38.fasta \
    {sample_name} \
    --max-reads-per-cnv 300 \
    --threads 32
bcftools index -tf {sample_name}.jalign.vcf.gz
samtools sort -o {sample_name}.jalign.sort.bam {sample_name}.jalign.bam
samtools index {sample_name}.jalign.sort.bam
```
Output files:
- **VCF file**: {sample_name}.jalign.vcf.gz - Final combined CNV calls with all annotations
- **VCF index file**: {sample_name}.jalign.vcf.gz.tbi - Corresponding index to output VCF
- **BAM file**: {sample_name}.jalign.sort.bam - Sorted BAM file with read evidence supporting CNV calls
- **BAM index file**: {sample_name}.jalign.sort.bam.bai - Index for the BAM file
- **CSV file**: {sample_name}.jalign.csv - CSV file with jalign scores for each read

#### Filter calls

The filtering process will annotate the callset with QUAL, and set LOW_SCORE/PASS filter string
Use `ugbio_filtering` docker to run it. 

The following annotations list should be provided for the current model: 

     [  "CNV_SOURCE",
        "RoundedCopyNumber",
        "CopyNumber",
        "CNMOPS_SAMPLE_STDEV",
        "CNMOPS_SAMPLE_MEAN",
        "CNMOPS_COHORT_STDEV",
        "CNMOPS_COHORT_MEAN",
        "pytorQ0",
        "pytorP2",
        "pytorRD",
        "pytorP1",
        "pytorP3",
        "CN",
        "GAP_PERCENTAGE",
        "CNV_DUP_READS",
        "CNV_DEL_READS",
        "CNV_DUP_FRAC",
        "CNV_DEL_FRAC",
        "JALIGN_DUP_SUPPORT",
        "JALIGN_DEL_SUPPORT",
        "JALIGN_DUP_SUPPORT_STRONG",
        "JALIGN_DEL_SUPPORT_STRONG",
        "SVTYPE",
        "SVLEN",
    ]

They should be provided in the following form: 

	--custom_annotation CNV_SOURCE --custom_annotation RoundedCopyNumber .... --custom_annotation SVLEN
	
	filter_variants_pipeline --input_file {sample_name}.jalign.vcf.gz \
			--model_file {input_model} \
            {custom_annotations}
			--decision_threshold 26
			--overwrite_qual_tag
			--output_file {sample_name}.filtered.vcf.gz

#### Combine redundant calls and produce final VCF

Use `ugbio_cnv` docker

    combine_cnmops_cnvpytor_cnv_calls merge_records \
            --input_vcf {sample_name}.filtered.vcf.gz \
            --output_vcf {sample_name}.mrg.vcf.gz \
            --distance 0 
