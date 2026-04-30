# Single sample germline CNV Calling

CNV calling pipeline to detect CNVs using multiple tools.
It runs:
1. cn.mops pipeline - using [cn.mops](https://github.com/Ultimagen/cn.mops)
2. cnvpytor pipeline - using [cnvpytor](https://github.com/abyzovlab/CNVpytor)
3. CNV candidates validation - using [jalign](https://github.com/Ultimagen/jalign)
4. results combination
5. (Optional) SV call integration - can merge with SV calls from [GRIDSS pipeline](https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/structural_variant_pipeline)

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
    s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v2.0/HapMap2_65samples_cohort_v2.0.plus_male.ploidy
	s3://ultimagen-workflow-resources-us-east-1/hg38/Homo_sapiens_assembly38.chr1-24.w1000.bed
    s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v2.1.2/ug_cnv_lcr.bed (used for cn.mops annotation only)

**Note:** For male samples use 
    s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v2.0/HapMap2_65samples_cohort_v2.0.plus_male.ploidy
    The ploidy file contains single number per row that indicates the expected ploidy of the X chromosome. The sample being 
    called corresponds to the **last** number in the ploidy file. 
	
ML filtering model can be found here: 
    
    s3://ultimagen-workflow-resources-us-east-1/filtering_models/cnv_filtering_model_v1.7.3.pkl

## Generating Germline CNV calls for a single sample

In the pipeline we run two CNV callers: cn.mops and CNVPytor, combine their results and validate candidates using our tool `jalign` - that checks if there are reads that can be realigned accross the breakpoint. 

Finally, we apply simple ML model to filter candidates according to the length, copy number, the tool that proposed the candidate and the jump reads support.

### Installation
* for running cnmops, cnvpytor, and post processing use ugbio_cnv docker: <br>
	Pull **ugbio_cnv** docker image :
	```
	docker pull ultimagenomics/ugbio_cnv:1.22.0
	```
	Run docker in interactive mode:
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_cnv:1.22.0 /bin/bash
	```
for latest docker version please see : (https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/germline_CNV_pipeline/tasks/globals.wdl)


* For running ML-based filtering - use **ugbio_filtering** docker image.
	```
	docker pull ultimagenomics/ugbio_filtering:1.22.0
	```
	Run docker in interactive mode:
	```
	docker run -it -v /data:/data ultimagenomics/ugbio_filtering:1.22.0 /bin/bash
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
    
    cnvpytor -root ~{sample_name}.pytor \
        -call <window_size> > {sample_name}.pytor.bin{window_size}.CNVs.1based.tsv

    ## output VCF
    cnvpytor -root {sample_name}.pytor  -view <window_size> <<EOF
            set print_filename {sample_name}.<window_size>.CNV.vcf
            print calls
            quit 
    EOF

    bcftools view -Oz -o {sample_name}.<window_size>.CNV.vcf.gz {sample_name}.<window_size>.CNV.vcf
    bcftools index -tf {sample_name}.<window_size>.CNV.vcf.gz
```

The callsets from the two window sizes are combined using `ugbio_cnv` docker. 
```
        combine_cnmops_cnvpytor_cnv_calls concat \
            --cnvpytor_vcf vcf1 vcf2 \
            --output_vcf {sample_name}.cnvpytor.cnv_calls.vcf.gz \
            --fasta_index Homo_sapiens_assembly38.fasta.fai \
            --make_ids_unique
```

### Combine, annotate and filter cn.mops and cnvpytor callsets

#### Annotate candidates: 
The following annotations are applied: 
1. Percentage of N bases in the CNV is recorded
2. Split read support is added
   
   ```
    # Concatenate cn.mops and cnvpytor VCF files
	combine_cnmops_cnvpytor_cnv_calls concat \
		--cnmops_vcf {cnmops_vcf} \
		--cnvpytor_vcf {cnvpytor_vcf} \
		--output_vcf {sample_name}.step1.vcf.gz \
		--fasta_index Homo_sapiens_assembly38.fasta.fai \
        --make_ids_unique

	# Annotate with gaps information
	combine_cnmops_cnvpytor_cnv_calls annotate_gaps \
		--calls_vcf {sample_name}.step1.vcf.gz \
		--output_vcf {sample_name}.combined.vcf.gz \
		--ref_fasta Homo_sapiens_assembly38.fasta

	
    # Add annotations of the split read support 
    analyze_cnv_breakpoint_reads \
        --bam-file {input_cram_bam_file} \
        --vcf-file {sample_name}.combined.vcf.gz \
        --reference-fasta Homo_sapiens_assembly38.fasta \
        --cushion 1500 \
        --output-file {sample_name}.split.annotated.vcf.gz \
        --output-bam {sample_name}.split.evidence.bam
    samtools sort {sample_name}.split.evidence.bam > {sample_name}.split.evidence.sort.bam
    samtools index {sample_name}.split.evidence.sort.bam
    ```

#### Run jump alignment

This task re-aligns reads that are close to the CNV breakpoints to check how many reads can be aligned in a way that supports the breakpoint. For deletion, it looks for the reads
that span the whole CNV. For the duplication, this looks for the reads that can align with the pattern consistent with the duplication.

For this task - use `ugbio_cnv` docker. Note that this tasks requires a 
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

#### Refine CNV breakpoints

This step refines CNV breakpoints using the evidence BAM files from both the split read analysis and jalign. It adds CIPOS (confidence interval around position) annotations to improve breakpoint accuracy.

Use `ugbio_cnv` docker to run it.

```
refine_cnv_breakpoints \
    --input-vcf {sample_name}.jalign.vcf.gz \
    --output-vcf {sample_name}.refined.tmp.vcf.gz \
    --bam-files {sample_name}.jalign.sort.bam {sample_name}.split.evidence.sort.bam
bcftools sort -Oz -o {sample_name}.refined.vcf.gz {sample_name}.refined.tmp.vcf.gz
bcftools index -tf {sample_name}.refined.vcf.gz
```

Output files:
- **VCF file**: {sample_name}.refined.vcf.gz - CNV calls with refined breakpoints and CIPOS annotations
- **VCF index file**: {sample_name}.refined.vcf.gz.tbi - Corresponding index to output VCF

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
        "DUP_READS_MEDIAN_INSERT_SIZE",
        "DEL_READS_MEDIAN_INSERT_SIZE",
        "SVTYPE",
        "SVLEN",
        "CIPOS"
    ]

They should be provided in the following form:

	--custom_annotation CNV_SOURCE --custom_annotation RoundedCopyNumber .... --custom_annotation SVLEN --custom_annotation CIPOS

	filter_variants_pipeline --input_file {sample_name}.refined.vcf.gz \
			--model_file {input_model} \
            {custom_annotations}
			--decision_threshold 31
			--overwrite_qual_tag
			--output_file {sample_name}.filtered.vcf.gz

#### Combine redundant calls and produce final VCF

In this step, overlapping or adjacent PASS-filtered CNV calls are consolidated into a single call, with aggregation of their INFO annotations.
Next, we merge nearby CNVs using a size-dependent criterion. Two CNVs are merged when the gap between them is smaller than `gap_scale_fraction × max(len(CNV1), len(CNV2))` and
less than `max_gap_absolute`.
Finally, filtered CNVs that overlap merged CNV calls are removed from the VCF.

Use `ugbio_cnv` docker.

    combine_cnmops_cnvpytor_cnv_calls merge_records \
            --input_vcf {sample_name}.filtered.vcf.gz \
            --output_vcf {sample_name}.mrg.vcf.gz \
            --distance 0 \
            --enable_smoothing \
            --max_gap_absolute 50000 \
            --gap_scale_fraction 0.05 \
            --cipos_threshold 50

#### Optional: Integrate SV calls

If you have structural variant (SV) calls available (e.g., from the GRIDSS pipeline: https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/structural_variant_pipeline), you can merge them with the CNV calls to improve breakpoint resolution.

This step should be performed after the collapse callset step above. Use `ugbio_cnv` docker.

```
combine_cnmops_cnvpytor_cnv_calls merge_cnv_sv \
    --cnv_vcf {sample_name}.mrg.vcf.gz \
    --sv_vcf {sv_calls_vcf} \
    --output_vcf {sample_name}.sv.cnv.merge.vcf.gz \
    --min_sv_length {min_cnv_length} \
    --fasta_index Homo_sapiens_assembly38.fasta.fai \
    --max_sv_length 5000000 \
    --min_sv_qual 600
```

**Parameters:**
- `--min_sv_length`: Minimum SV length to consider for merging with CNV calls
- `--max_sv_length`: Maximum SV length to consider (default: 5000000)
- `--min_sv_qual`: Minimum SV quality score to consider (default: 600)

The resulting VCF will contain both CNV and SV calls, with improved breakpoint positions where CNV and SV calls overlap.

## Output Files and Interpretation

The germline CNV calling pipeline generates the following final output files:

### 1. Final Merged VCF (`{sample_name}.mrg.vcf.gz` or `{sample_name}.sv.cnv.merge.vcf.gz`)

Final deduplicated CNV callset where overlapping/redundant calls have been merged. If SV calls were provided, the output will be `{sample_name}.sv.cnv.merge.vcf.gz` which includes integrated SV information. **This is the recommended file for downstream analysis and interpretation.**

**Processing:**
- Adjacent or overlapping CNV calls from the same caller are merged
- Keeps the most informative annotations from merged calls
- FILTER annotates high confidence calls after ML filtering (FILTER=PASS)
- If SV calls are provided, CNV breakpoints are refined using overlapping SV calls

**Key INFO Fields:**
- `SVTYPE`: DEL (deletion) or DUP (duplication)
- `SVLEN`: Length of the CNV in base pairs
- `CIPOS`: Confidence interval around the breakpoint (estimated from the CNV, split reads and the SV integration (optional))
- `CopyNumber`: Estimated copy number (float)
- `RoundedCopyNumber`: Rounded copy number (integer)
- `CNV_SOURCE`: Tool(s) that called the CNV (cn.mops, cnvpytor, or both)
- `CNMOPS_SAMPLE_MEAN`: Mean coverage in CNV region for the sample
- `CNMOPS_SAMPLE_STDEV`: Standard deviation of coverage in CNV region
- `CNMOPS_COHORT_MEAN`: Mean coverage in CNV region across cohort
- `CNMOPS_COHORT_STDEV`: Standard deviation of cohort coverage
- `pytorQ0`, `pytorP1-P3`, `pytorRD`: CNVpytor quality metrics
- `GAP_PERCENTAGE`: Fraction of N bases in the CNV region
- `CNV_DUP_READS`, `CNV_DEL_READS`: Originally aligned split reads supporting duplication/deletion
- `JALIGN_DUP_SUPPORT`, `JALIGN_DEL_SUPPORT`: Reads realigning across breakpoint consistent with duplication / deletion pattern
- `JALIGN_DUP_SUPPORT_STRONG`, `JALIGN_DEL_SUPPORT_STRONG`: Reads realigning across breakpoint consistent with duplication / deletion pattern
- `QUAL`: ML-based quality score, PHRED scale (0-100), above 30 is usually good

**FILTER Fields:**
- `PASS`: High-confidence call (QUAL ≥ 31 with default threshold)
- `LOW_SCORE`: Lower confidence call, may be false positive

**FORMAT Fields:**
- `GT`: Genotype
  - `0/1`: Heterozygous (single copy gain/loss)
  - `1/1`: Homozygous (two copy deletion or multi-copy gain)
  - `./1`: Unknown/variable state

### 2. Output Plots

The pipeline generates visualization plots to help interpret CNV calls. the plots are of high resolution and might not present short CNVs. :

#### Coverage Plot (`{sample_name}.CNV.coverage.jpeg`)

Shows normalized (log2 scale) coverage along the genome for the germline sample to visually identify regions with coverage changes.

**Interpretation:**
- X-axis: Location along the genome (separated by chromosome)
- Y-axis: Log2-transformed normalized coverage
- Blue dots: Germline sample coverage
- Deviations from baseline (y=0) indicate coverage changes
- Vertical black lines separate chromosomes

#### Duplication and Deletion Plot (`{sample_name}.dup_del.calls.jpeg`)

Shows called duplications and deletions as horizontal lines along the genome for a quick overview of the CNV distribution.

**Interpretation:**
- X-axis: Location along the genome (separated by chromosome)
- Green lines: Duplication calls
- Red lines: Deletion calls
- Line length represents CNV size
- Vertical black lines separate chromosomes

#### Copy Number Plot (`{sample_name}.CNV.calls.jpeg`)

Shows the estimated copy number along the genome for all CNV calls.

**Interpretation:**
- X-axis: Location along the genome (separated by chromosome)
- Y-axis: Copy number (0-8+)
- Horizontal gray line marks as the neutral ploidy (2)
- Black horizontal lines represent CNV segments with their copy number
- Vertical black lines separate chromosomes

### 3. Additional Output Files

**Combined CNV Calls BED** (`{sample_name}.combined_cnv_calls.bed`)
- BED format file with all CNVs from the final VCF
- Format: `chr`, `start`, `end`, `cnv_info`
- Useful for downstream tools requiring BED format

**Filtered BED Files** (optional outputs from individual callers)
- `{sample_name}.cnmops.bed`: cn.mops-only calls
- `{sample_name}.cnvpytor.bed`: CNVpytor-only calls
- Useful for comparing caller-specific results
