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
2. BedGraph - previously calculated input bedGraph holding the coverage per base (outputs with UG sequencing data), seperatly for tumor and normal.
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
    1. [ugbio_freec_docker Dockerfile](https://github.com/Ultimagen/ugbio-utils/blob/main/src/freec/Dockerfile) used to build the ugbio_freec_docker image.
    2. [ugbio_cnv_docker Dockerfile](https://github.com/Ultimagen/ugbio-utils/blob/main/src/cnv/Dockerfile) used to build the ugbio_cnv_docker image.

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
	--forceGCcontentNormalization 1 \
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

#### Process ControlFREEC called CNVs - convert them to VCF and filter by length and low confidence regions
from inside ugbio_cnv docker image:
```
#convert to bedfile
cat {tumor}_CNVs | sed 's/^/chr/' | cut -f1-4 > {tumor}.cnvs.bed

#annotate cnvs bed file and convert to VCF
process_cnvs \
	--sample_name {sample_name} \
	--input_bed_file {tumor}.cnvs.bed \
	--intersection_cutoff 0.5 \
	--cnv_lcr_file ug_cnv_lcr.bed \
	--min_cnv_length 10000 \
	--out_directory . \
	--fasta_index_file Homo_sapiens_assembly38.fasta.fai
```

### Generate CNV visualization plots
from inside ugbio_cnv docker image:
```
#prepare coverage files in BED format
awk '{print $1"\t"$2"\t"$2+999"\t"$NF}' {normal}.cpn | sed 's/^/chr/' > normal_coverage.cpn.bed
awk '{print $1"\t"$2"\t"$2+999"\t"$NF}' {tumor}.cpn | sed 's/^/chr/' > tumor_coverage.cpn.bed

#extract PASS CNVs and separate by type
bcftools view -f "PASS" {sample_name}.cnvs.annotate.vcf.gz | \
	bcftools query -f '%CHROM\t%POS0\t%END\t%CopyNumber\n' > helper.bed

#separate duplications and deletions based on ploidy
awk '$4>{ploidy}' helper.bed > {sample_name}.cnvs.filter.DUP.bed
awk '$4<{ploidy}' helper.bed > {sample_name}.cnvs.filter.DEL.bed

#create output directory for plots
mkdir CNV_figures

#generate CNV coverage, duplication/deletion, and copy number plots
plot_cnv_results \
	--germline_coverage normal_coverage.cpn.bed \
	--tumor_coverage tumor_coverage.cpn.bed \
	--duplication_cnv_calls {sample_name}.cnvs.filter.DUP.bed \
	--deletion_cnv_calls {sample_name}.cnvs.filter.DEL.bed \
	--sample_name {sample_name} \
	--out_directory CNV_figures \
	--neutral_ploidy {ploidy}

#generate neutral allele frequency plot
plot_FREEC_neutral_AF \
	--mpileup {tumor}_minipileup.pileup \
	--cnvs_file helper.bed \
	--sample_name {sample_name} \
	--out_directory CNV_figures
```

## Output Files and Interpretation

The pipeline generates the following final output files:

### 1. Final VCF Output (`{sample_name}.cnvs.annotate.vcf.gz`)

The final output is a VCF file containing filtered and annotated CNV calls. This is the primary output for downstream analysis.

#### VCF Header Fields

**ALT Types:**
- `<DEL>`: Deletion (copy number < 2)
- `<DUP>`: Duplication (copy number > 2)
- `<CNV>`: General copy number variant

**INFO Fields:**
- `SVTYPE`: Type of structural variant (DEL or DUP)
- `SVLEN`: Length of the CNV in base pairs
- `CopyNumber`: Estimated copy number (float)
- `RoundedCopyNumber`: Rounded copy number (integer)
- `CNV_SOURCE`: Tool that called this CNV (e.g., "controlFREEC")
**FILTER Fields:**
- `PASS`: CNV passes all quality filters
- `UG-CNV-LCR`: Overlaps with Ultima Genomics CNV low-complexity regions (>50% overlap by default)

**FORMAT Fields:**
- `GT`: Genotype
  - `0/1`: Heterozygous (copy number = 1 for deletion)
  - `1/1`: Homozygous (copy number = 0 for deletion)
  - `./1`: Unknown/variable (for duplications or uncertain copy states)

#### Example VCF Record

```
chr1	934	cnmops_dup_1	N	<DUP>	.	PASS	SVTYPE=DUP;SVLEN=3962;CopyNumber=3.2;RoundedCopyNumber=3	GT	./1
```

**Interpretation:**
- A duplication on chr1 starting at position 934
- Length: 3,962 bp
- Estimated copy number: 3.2 (approximately 3 copies)
- Passes all quality filters
- Genotype is uncertain/variable (typical for duplications)

### 2. Output Plots

The pipeline generates visualization plots to help interpret CNV calls:

#### Coverage Plot (`{sample_name}.CNV.coverage.jpeg`)

Shows normalized (log2 scale) coverage along the genome for germline (blue) and tumor (orange) samples.

**Interpretation:**
- X-axis: Location along the genome (separated by chromosome)
- Y-axis: Log2-transformed normalized coverage
- Deviations from baseline (y=0) indicate coverage changes
- Blue dots: Germline sample coverage
- Orange dots: Tumor sample coverage
- Vertical black lines separate chromosomes
- Coverage changes in tumor but not germline suggest somatic CNVs


#### Duplication and Deletion Plot (`{sample_name}.dup_del.calls.jpeg`)

Shows called duplications and deletions as horizontal lines along the genome.

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
- Horizontal gray line at y=3 marks threshold above diploid
- Black horizontal lines represent CNV segments with their copy number
- Vertical black lines separate chromosomes


#### Neutral Allele Frequency Plot (`{mpileup_basename}.freq.SNP.neutral.hist.jpeg`)

Histogram showing allele frequency distribution in neutral (non-CNV) regions.

**Interpretation:**
- X-axis: Allele frequency (0.0 to 1.0)
- Y-axis: Count of SNPs (log scale)
- Expected peak at ~0.5 for heterozygous SNPs in diploid regions
- Distribution shape indicates sample quality and tumor purity

**Usage:**
- Assess tumor purity (peak shift from 0.5 indicates contamination or purity issues)
- Quality control for allele frequency calculations
- Validate that neutral regions show expected diploid AF distribution
- Deviations from expected distribution may indicate:
  - Low tumor purity (peak < 0.5)
  - contamination
  - Technical artifacts in allele frequency calculation

### 3. Additional Output Files

**Filtered CNV BED File** (`{sample_name}.cnvs.filter.bed`)
- BED format file containing only high-confidence CNVs (FILTER=PASS)
- Format: `chr`, `start`, `end`, `copy_number`
- Recommended for downstream analysis requiring BED format

**Coverage Files** (`{tumor/normal}_coverage.cpn`)
- Binned coverage files in CPN format (input to controlFREEC)
- Can be reused for additional analyses

**Mpileup Files** (`{tumor/normal}_mpileup.pileup`)
- Allele frequency data at common variant positions
- Can be used for additional allele frequency analyses

**controlFREEC_info** (`*_info.txt`)
holds information for : 
- controlFREEC run parameters
- estimated tumor sample ploidy
- estimated tumor sample purity

### Quality Considerations

**High-confidence CNV calls should:**
1. Have `FILTER=PASS` (no quality flags)
2. Have length ≥ 10kb (unless using custom thresholds)
3. Not overlap with low-complexity regions (UG-CNV-LCR)
4. Show consistent copy number ratio across the region (low `sd_ratio`)

**Calls with filters should be interpreted cautiously:**
- `UG-CNV-LCR`: May be technical artifacts due to mapping issues in low-complexity regions
- `LEN`: Short CNVs (<10kb) have higher false positive rates
- Multiple filters indicate lower confidence

### Recommended Analysis Workflow

1. **Start with PASS calls**: Filter VCF for `FILTER=PASS` to get high-confidence CNVs
2. **Check CNV size**: Larger CNVs (>50kb) are generally more reliable than smaller ones
3. **Visual validation**: Use plots from `plot_FREEC_fold_change` to visualize fold-change across the genome
4. **Compare with germline**: For somatic analysis, ensure CNVs are present in tumor but not in matched normal