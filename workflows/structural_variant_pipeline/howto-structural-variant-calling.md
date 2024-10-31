# Germline/Somatic Structural Variant (SV) Calling 

SV calling pipeline detects SVs using an in-house assembly tool, GRIDSS tool for SV detection and GRIPSS tool for filtering the output callset. The pipeline supports both germline and somatic modes. 


## Workflow Overview

The pipeline takes as input aligned BAM/CRAM files and outputs a filtered VCF containing annotated SV calls. The steps of the pipeline are as follows: assembly, realignment with Ultima Aligner, reverting secondary low MAPQ alignments, GRIDSS identification of SVs, GRIDSS annotation and SV filtering using GRIPSS for somatic runs and GermlineLinkVariants which links GRIDSS breakends to create an SV for germline runs. 


## Requirements

1. Input CRAM files and respective indexes for germline or matched normal sample.
2. Input CRAM files and respective indexes for the tumor (in case of matched T/N calling).
3. For tumor only structural variant calling, Input CRAM and respective index from the tumor. 
4. Installation of picard and GATK tools.

### Files required for the analysis (download locally)
The following files are publicly available:

    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
    gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
    gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list

The following file is required as input for UA realignment:
    
    gs://concordanz/hg38/UA/Homo_sapiens_assembly38.fasta.uai

#### Generate scattered intervals to run SV calling on a single interval

hg38 interval file can be downloaded from: 
```
gs://concordanz/sv/temp/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list
```
Assembly is run on a small interval in the genome determined by an input bed file. This can be used to scatter the genome into small intervals, and then parallelize the workflow across these intervals. A convenient tool to generate scattered intervals is picard IntervalListTools. A typical command to scatter intervals provided below. Note that the input to this command is Picard's interval_list format, and not the simple bed format:

```
mkdir out && \
picard \
  IntervalListTools \
  SCATTER_COUNT=400 (for somatic, or 200 for germline) \
  SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
  UNIQUE=true \
  SORT=true \
  BREAK_BANDS_AT_MULTIPLES_OF=10000000 \
  INPUT=wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list \
  OUTPUT=out
```

#### Run Assembly on a single resulting interval

Assembly docker:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/haplotype:ebbe0bd
```

Generate interval bed : 
```
gatk IntervalListToBed -I out/001.scattered.interval_list -O interval.bed
```

T/N somatic calling: 
```
tool \
	--input input_tumor.cram\;input_germline.cram \
	--cram-index input_tumor.cram.crai\;input_germline.cram.crai \
	--somatic \
	--output output_basename_assembly \
	--ref Homo_sapiens_assembly38.fasta \
	--bed interval.bed \
	--min-base 0 \
	--min-mapq 5 \
	--min-feature-length '20;0' \
	--max-num-haps 10 \
	--prog \
	--interval-nreads 10000 \
	--sv \
	--single-strand-filter 	

samtools view -bS output_basename_assembly_hap_out.sam > output_basename_assembly_hap_out.bam
samtools sort output_basename_assembly_hap_out.bam -o output_basename_assembly_hap_out_sorted.bam
samtools index output_basename_assembly_hap_out_sorted.bam
```

Somatic calling (unmatched tumor)

```
tool \
	--input input_tumor.cram \
	--cram-index input_tumor.cram.crai \
	--output output_basename_assembly \
	--somatic \
	--ref Homo_sapiens_assembly38.fasta \
	--bed interval.bed \
	--min-base 0 \
	--min-mapq 5 \
	--min-feature-length '10' \
	--max-num-haps 10 \
	--prog \
	--interval-nreads 10000 \
	--sv \
	--single-strand-filter
samtools view -bS output_basename_assembly_hap_out.sam > output_basename_assembly_hap_out.bam
samtools sort output_basename_assembly_hap_out.bam -o output_basename_assembly_hap_out_sorted.bam
samtools index output_basename_assembly_hap_out_sorted.bam
```

Germline calling:

```
tool \
	--input input_germline.cram \
	--cram-index input_germline.cram.crai \
	--output output_basename_assembly \
	--ref Homo_sapiens_assembly38.fasta \
	--bed interval.bed \
	--min-base 0 \
	--min-mapq 5 \
	--min-feature-length '10' \
	--max-num-haps 10 \
	--prog \
	--interval-nreads 10000 \
	--sv
samtools view -bS output_basename_assembly_hap_out.sam > output_basename_assembly_hap_out.bam
samtools sort output_basename_assembly_hap_out.bam -o output_basename_assembly_hap_out_sorted.bam
samtools index output_basename_assembly_hap_out_sorted.bam
```

#### Realign with UA to imporves recall for deletions (optional)

UA Docker:

```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/ua:master_c553451
```
UA realignment command: (Realignment is done on a merged bam consisting of a merge of all the BAMs produced in the scattered assembly)
```
samtools view -h -@ 40 output_basename_assembly_hap_out_sorted.bam | \
/ua/ua \
    --index Homo_sapiens_assembly38.fasta.uai \
    --align true \
    --progress \
    --tp reference \
    --alt=Homo_sapiens_assembly38.fasta.64.alt \
    --json=output_basename_assembly_file_ua_aligned.%s.json \
    --nthread max \
    --sam-input - \
    --sam-output - \
    --seed-score-ratio 0.5 --vector --huge --soft-clipping --realignment-tag re --mismatch-cost -6 | \
samtools view -@ 40 -o output_basename_assembly_file_ua_aligned.bam -

samtools sort output_basename_assembly_file_ua_aligned.bam -o output_basename_assembly_file_ua_aligned_sorted.bam 
samtools index output_basename_assembly_file_ua_aligned_sorted.bam
```

GRIDSS Docker :

UG modifications to the GRIDSS pipeline are here:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/gridss:c0113ac
```
Reverting secondary low MAPQ alignments (on GRIDSS docker):

```
python3 /opt/gridss/revert_sup_low_mapq_ua_alignment.py
	--before output_basename_merged_assembly.bam \
	--after output_basename_assembly_file_ua_aligned.bam_sorted.bam \
	--output output_basename_assembly_ua_realigned_unsorted.bam \
	--min_mapping_quality 60

samtools sort -o output_basename_assembly_ua_realigned.bam output_basename_assembly_ua_realigned_unsorted.bam
samtools index output_basename_assembly_ua_realigned.bam 
```

##### GRIDSS Identify Variants

The gridss.config file we used was:

```
variantcalling.breakpointLowQuality = 120.0
variantcalling.breakendLowQuality = 120.0
variantcalling.writeFiltered = true
variantcalling.requireReadPair = false
variantcalling.requireSplitRead = false
variantcalling.callUnassembledBreakends = true
variantcalling.callUnassembledBreakpoints = true
variantcalling.breakendMaxAssemblySupportBias = 100
variantcalling.breakpointMaxAssemblySupportBias = 100
variantcalling.requiredReadAndAssemblyBreakpointOverlap = 5
``` 

Blacklist file (low confidence regions to avoid) can be downloaded from https://github.com/PapenfussLab/gridss/: 
```
example/ENCFF356LFX.bed
```

T/N IdentifyVariants:

```
java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
    SAMPLE_NAMES=~{tumor_sample_name} \
    SAMPLE_NAMES=~{germline_sample_name} \
    O={output_basename}.vcf \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=output_basename_assembly_ua_realigned.bam \ 
    C=gridss.config

bcftools view ~{output_basename}.vcf -Oz ~{output_basename}.vcf.gz
bcftools index ~{output_basename}.vcf.gz 
```

Tumor only IdentifyVariants:

```
java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
    SAMPLE_NAMES=~{tumor_sample_name} \
    O=~{output_basename}.vcf \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=output_basename_assembly_ua_realigned.bam \ 
    C=gridss.config

bcftools view ~{output_basename}.vcf -Oz ~{output_basename}.vcf.gz
bcftools index ~{output_basename}.vcf.gz 
```

Germline IdentifyVariants:

```
java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
    SAMPLE_NAMES=~{germline_sample_name} \
    O=output_basename.vcf \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=output_basename_assembly_ua_realigned.bam \ 
    C=gridss.config

bcftools view ~{output_basename}.vcf -Oz ~{output_basename}.vcf.gz
bcftools index ~{output_basename}.vcf.gz 

```

#### GRIDSS AnnotateVariants

Prefilter Candidates - Removal of variants with QUAL below 120 : 

```
# Prefilter VCF  

bcftools view -Oz -i "QUAL>120" ~{output_basename}.vcf.gz -o ~{output_basename.filtered}.vcf.gz

bcftools index -t ~{output_basename.filtered}.vcf.gz

bcftools view -R interval.bed ~{output_basename.filtered}.vcf.gz -o output_region.vcf 
input_vcf=output_region.vcf
```

Use the same gridss.config file for this step

T/N  AnnotateVariants 

Note that this step is run in a scatter mode using the intervals generated in the IntervalListTools task above

```
samtools view -b -L interval.bed input_tumor.cram -o input_partial_tumor.bam
samtools index input_partial_tumor.cram
samtools view -b -L interval.bed input_germline.cram -o input_partial_germline.bam
samtools index input_partial_germline.cram
samtools view -b -L interval.bed ~{output_basename_assembly}_ua_realigned.bam -o ~{output_basename_assembly}_partial_ua_realigned.bam
samtools index ~{output_basename_assembly}_partial_ua_realigned.bam
```


```
java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \ 
    INPUT_VCF=output_region.vcf \
    R=Homo_sapiens_assembly38.fasta \
    I=input_partial_tumor.bam \
    I=input_partial_germline.bam \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=~{output_basename_assembly}_partial_ua_realigned.bam \
    C=gridss.config \
    OUTPUT_VCF=~{output_basename}.ann.vcf
```
Tumor only  AnnotateVariants:

```

java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \ 
    INPUT_VCF=output_region.vcf \
    R=Homo_sapiens_assembly38.fasta \
    I=input_partial_tumor.bam \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=~{output_basename_assembly}_partial_ua_realigned.bam \
    C=gridss.config \
    OUTPUT_VCF=~{output_basename}.ann.vcf
```


Germline AnnotateVariants: 

```
java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \ 
    INPUT_VCF=output_basename.vcf \
    R=Homo_sapiens_assembly38.fasta \
    I=input_partial_germline.bam \
    BLACKLIST=ENCFF356LFX.bed \
    ASSEMBLY=~{output_basename_assembly}_partial_ua_realigned.bam \
    C=gridss.config \
    OUTPUT_VCF=output_basename.ann.vcf

bcftools view output_basename.ann.vcf -oZ output_basename.ann.vcf.gz
bcftools index -t output_basename.ann.vcf.gz
```
#### Germline filtering:

For germline mode, run GermlineLinkVariants (on the GRIDSS docker) using :

```
Rscript /opt/gridss/link_breakpoints \
    --input  output_basename.ann.vcf.gz \
    --fulloutput output_basename_linked.vcf \
    --ref BSgenome.Hsapiens.UCSC.hg38 \
    --scriptdir /opt/gridss/
```
#### GRIPSS filtering for somatic runs:

Download directory with necessary reference files with : 

```
gsutil -m cp -r gs://concordanz/sv/gripss/ .
```

GRIPSS docker: (run on a merge of the vcfs produced in scatter by AnnotateVariants above)
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/gripss:58cba04b7b
```

T/N GRIPSS:

```
mkdir gripss_output

java -jar gripss.jar \
	-vcf output_basename.ann.vcf \
	-sample vcf_tumor_sample_name \
	-reference vcf_normal_sample_name \
	-ref_genome_version 38 \
	-ref_genome Homo_sapiens_assembly38.fasta \
	-pon_sgl_file sv/gripss/sgl_pon.38.bed.gz \
	-pon_sv_file sv/gripss/sv_pon.38.bedpe.gz \
	-known_hotspot_file sv/gripss/known_fusions.38.sorted.bedpe \
	-repeat_mask_file sv/gripss/repeat_mask_data.38.fa.gz \
	-exclude_filters SHORT_SR_NORMAL \
	-output_dir gripss_output
```
Tumor only GRIPSS:
```
java -jar gripss.jar \
	-vcf output_basename.ann.vcf \
	-sample vcf_tumor_sample_name \
	-ref_genome_version 38 \
	-ref_genome Homo_sapiens_assembly38.fasta \
	-pon_sgl_file sv/gripss/sgl_pon.38.bed.gz \
	-pon_sv_file sv/gripss/sv_pon.38.bedpe.gz \
	-known_hotspot_file sv/gripss/known_fusions.38.sorted.bedpe \
	-repeat_mask_file sv/gripss/repeat_mask_data.38.fa.gz \
	-min_normal_coverage 8 \
	-exclude_filters SHORT_SR_NORMAL \
	-output_dir gripss_output
```

