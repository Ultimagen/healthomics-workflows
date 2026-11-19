# Calling variants from aligned cram files using efficient DV

## Overview

Efficient DV is an analysis pipeline designed to call variants from aligned cram files using the method of  [DeepVariant](https://www.nature.com/articles/nbt.4235) , adapted for Ultima Genomics data. There are three stages to the variant calling:
1. make_examples - Looks for “active regions” with potential candidates. Within these regions, it performs local assembly (haplotypes), re-aligns the reads, and defines candidate variant. Images of the reads in the vicinity of the candidates are saved as protos in a tfrecord format.
2. call_variants - Collects the images from make_examples and uses a deep learning model to infer the statistics of each variant (i.e. quality, genotype likelihoods etc.), using the TensorRT framework. It outputs sorted results in a CallVariantsOutput protos saved as tfrecords files.
3. post_process - Uses the call_variants tfrecord output. It resolves multi-allelic records and variants that overlap with indels. It annotates the variants with information such as variant type, cycle skip status and genomic intervals it appears in (e.g. exome). It filters the variants based on defined thresholds. Finally, it writes out a VCF file. Post_process can also output a gVCF file.

## Requirements


### Input files

The workflow takes three major inputs:

1. Reads aligned to a reference genome, sorted and duplicate-marked, stored in CRAM format, and indexed with CRAI.
2. Reference genome files (the same reference used to align the reads). Typically hg38 is used, and is publicly available in:
```
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list
or
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/wgs_calling_regions.hg38.interval_list
```
3. A model in ONNX format. For WGS data use:
```
gs://concordanz/deepvariant/model/germline/wgs/v2.0/ultimagen-germline-wgs-solaris2-hg38-regnet-v2.0.gpu-ens.onnx
or
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/wgs/v2.0/ultimagen-germline-wgs-solaris2-hg38-regnet-v2.0.gpu-ens.onnx
```
For paired-end (PE) WGS data use:
```
gs://concordanz/deepvariant/model/germline/wgs/pe/v1.2/ultimagen-germline-wgs-pe-v1.2.gpu-ens.onnx
or
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/wgs/pe/v1.2/ultimagen-germline-wgs-pe-v1.2.gpu-ens.onnx
```

**NOTE:** the exact model may differ between use-cases, we recommend consulting with the parameters templates provided to find the exact model

### Dockers and hardware requirements

The Efficient DV analysis pipeline is split into two docker images:

1. `make_examples` docker - contains binaries for the make_examples and post_process steps. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/make_examples:3.5.0
or
ultimagenomics/make_examples:3.5.0
```
2. `call_variants` docker - contains the binaries for the call_variants step. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/ncv:5.0.0
or
ultimagenomics/ncv:5.0.0
```

The make_examples and post_process steps are run on a single CPU. make_examples requires up to 2 GB of memory for each thread. post_process requires 8 GB of memory and runs on a single thread.

call_variants runs on a machine with a single NVIDIA GPU (e.g. T4, A10G or L4) and 8 CPUs, which are used for reading the input tfrecord files. It requires 24 GB of memory.

## Workflow details

The workflow is composed of three steps, as described above: make_examples, call_variants and post_process.

### Ploidy estimation

Optionally, sample sex/ploidy can be estimated from the filtered VCF after variant calling. The `estimate_ploidy` command uses PASS biallelic SNP depth and B-allele fractions to estimate per-chromosome ploidy and a sex karyotype:

```bash
estimate_ploidy \
  --vcf input.vcf.gz \
  --sample-id sample_name \
  --sex-chromosomes chrX chrY X Y \
  --output-dir .
```

The karyotype is reported in `sample_name.ploidy_report.txt`. For XY samples, run `convert_haploid_regions` on the VCF to convert non-PAR chrX/Y genotypes to haploid while preserving pseudoautosomal-region calls as diploid:

```bash
convert_haploid_regions \
  --input_vcf input.vcf.gz \
  --par_regions par_regions.bed \
  --output_vcf sample_name.haploid_chrXY.vcf.gz
```

Index the output VCF with `bcftools index -t` and use it as the final VCF for XY samples. For XX or undetermined karyotypes, keep the original VCF.

The PAR regions for hg38 are available in:
```
gs://concordanz/hg38/v0/par_regions.hg38.bed
or
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed
```

Ploidy can also be estimated from the CRAM alone, without a VCF, using a `mosdepth` coverage summary:
```bash
mosdepth --fasta Homo_sapiens_assembly38.fasta --no-per-base --fast-mode sample_name input_reads.cram
estimate_ploidy \
  --mosdepth-summary sample_name.mosdepth.summary.txt \
  --sample-id sample_name \
  --sex-chromosomes chrX chrY X Y \
  --output-dir .
```

### Running make_examples

make_examples is typically run on a small interval in the genome determined by an input bed file. This can be used to scatter the genome into small intervals, and then parallelize the workflow across these intervals. A convenient tool to generate scattered intervals is [picard IntervalListTools](https://gatk.broadinstitute.org/hc/en-us/articles/360036897212-IntervalListTools-Picard-). A typical command to scatter intervals provided below. Note that the input to this command is Picard's interval_list format, and not the simple bed format:

```
mkdir out && \
picard \
  IntervalListTools \
  SCATTER_COUNT=40 \
  SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
  UNIQUE=true \
  SORT=true \
  BREAK_BANDS_AT_MULTIPLES_OF=100000 \
  INPUT=wgs_calling_regions.hg38.interval_list \
  OUTPUT=out
```

The output interval_list can be converted to bed files using:

```
cat interval001.interval_list | grep -v @ | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3}' > interval001.bed
```

make_examples is invoked using `tool` command in the `make_examples` docker. A typical command line (from within the docker) will look like:

```
tool \
  --input input_reads.cram \
  --cram-index input_reads.cram.crai \
  --bed interval001.bed \
  --output 001 \
  --reference Homo_sapiens_assembly38.fasta \
  --min-base-quality 5 \
  --min-mapq 5 \
  --cgp-min-count-snps 2 \
  --cgp-min-count-hmer-indels 2 \
  --cgp-min-count-non-hmer-indels 2 \
  --cgp-min-fraction-snps 0.12 \
  --cgp-min-fraction-hmer-indels 0.12 \
  --cgp-min-fraction-non-hmer-indels 0.06 \
  --cgp-min-fraction-single-strand-non-snps 0.15 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --optimal-coverages 50 \
  --median-coverage <median_coverage> \
  --add-ins-size-channel \
  --single-strand-filter
```

The input cram files and the corresponding index files are provided to `--input` and `--cram-index`, respectively. Multiple cram files can be provided as a comma separated list.

The `--output` argument is the prefix for the output files (including tfrecords).

The `--optimal-coverages` and `--add-ins-size-channel` are parameters which affect the way images are produced, and should be aligned with the model. `optimal-coverages` is related to how reads are internally downsampled before the image is produced, and `add-ins-size-channel` adds a channel with the length of the insertion. `--median-coverage` provides `make_examples` with the actual sequencing depth of the sample, used alongside `--optimal-coverages` to control read downsampling.

The program will output a sam file with the re-aligned reads unless the argument `--no-realigned-sam` is provided. Note that these files are very large, so provide a large disk space if you want to save the re-aligned reads.

### Running call_variants

call_variants combines the tfrecords from all make_examples jobs. List the tfrecord files, one per line, in a text file (do not include the `gvcf.tfrecord.gz` files), and run from within the call_variants docker:

```bash
printf '%s\n' 001.tfrecord.gz 002.tfrecord.gz ... > examples.txt

call_variants \
  --model model.onnx \
  -T @examples.txt \
  -n 8 \
  --fp16 \
  --second-pass-threshold 0.995 \
  --second-pass-model-ensemble-k 5 \
  --second-pass-always-replace
```

`-n` sets the number of threads used for reading the input files. call_variants writes its output, `*_cvo.gz` files, to the current directory (use `-o <folder>` to change it). These files are the input to post_process.

All candidates are first called by the model. Candidates whose score is below `--second-pass-threshold` are called again by an ensemble of 5 (`--second-pass-model-ensemble-k 5`), and the ensemble result replaces the first-pass result (`--second-pass-always-replace`). To call with a single pass only, omit the `--second-pass-*` arguments.

### Running post_process:

post_process uses the `*_cvo.gz` output files of call_variants and generates a VCF file. A typical post_process command from within the docker will look like (same docker as make_examples):
```
ug_postproc \
  --infile shard001_cvo.gz,shard002_cvo.gz \
  --ref Homo_sapiens_assembly38.fasta \
  --outfile output_prefix.vcf.gz \
  --consider_strand_bias \
  --flow_order TGCA \
  --annotate \
  --bed_annotation_files exome.twist.bed,ug_hcr.bed,... \
  --qual_filter 1 \
  --filter \
  --filters_file filters.txt \
  --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf
```

Using `--annotate` together with `--bed_annotation_files` adds annotations to the VCF based on the provided BED files (given as a comma-separated list). A description of each BED file is recorded in the INFO field of the VCF. The description can be specified in one of two ways:

Using a `##INFO` in the header of the bed file. For example:
```
##INFO=<ID=EXOME,Number=1,Type=String,Description="Genomic Region Annotation: In the exome (s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed)">
chr1    69090   70008   TRUE
chr1    450739  451678  TRUE
```

If `##INFO` is not present in the bed file, then a json file with the same name as the bed file (e.g. exome.twist.bed.json) should be provided. The json will contain the ID, Type and Description:
```
{
 "ID": "EXOME",
 "Type": "String",
 "Description": "Genomic Region Annotation: In the exome (s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed)"
}
```

If `--filter` argument is used, the vcf will be filtered based on the criteria in `--filters_file`. In this file, each filter is composed of two lines, the first is the filter name (which will appear in FILTER column of the vcf), and the second is the expression for the filter. The syntax of the expression follows [JEXL filtering expressions](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions). Below is an example of a filters_file that is typically used. Note that these filters use the EXOME attribute, which is an annotation that was added using a bed file, if there is no annotation named EXOME this filter will fail.
```
LowQualInExome
QUAL < 16 and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and vc.hasAttribute('EXOME')
LowQual
QUAL < 16 and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and not vc.hasAttribute('EXOME')
LowQual
QUAL < 16 and VARIANT_TYPE=='non-h-indel' and not vc.isFiltered()
LowQual
QUAL < 14 and VARIANT_TYPE=='snp' and not vc.isFiltered()
LargeDeletion
REFLEN > 220 and vc.isFiltered()
```

dbSNP data can be downloaded from: gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
or s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

### Workflow That Includes Haplotype Data from Pangenomes

As demonstrated by [Asri et al.](https://www.biorxiv.org/content/10.1101/2025.06.05.657102v1), accuracy of variant calling can be significantly improved by using pangenome-aligned reads as input to DeepVariant, together with pangenome-derived personalized haplotype data.

Pangenome-aligned reads can be generated using [`vg giraffe`](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe). Haplotypes are sample-specific sequences selected from a pangenome graph based on the k-mer content of the input reads, and aligned to the reference genome. For detailed instructions on generating haplotypes from pangenome graphs, see [howto-haplotype-sampling.md](howto-haplotype-sampling.md).

*Important note:* The number of haplotypes in the image should fit the model that is being used in the call_variants step.

#### Running `make_examples` with Pangenome-aligned Haplotype Data

Below is the command to generate images of pangenome-aligned data that include haplotype data. The key differences from standard variant calling are the addition of `--exp-pangenome-haps haplotypes.cram` and `--min-mapq 1`:

```bash
tool \
  --input input_reads.cram \
  --cram-index input_reads.cram.crai \
  --bed interval001.bed \
  --output 001 \
  --reference Homo_sapiens_assembly38.fasta \
  --min-base-quality 5 \
  --min-mapq 1 \
  --cgp-min-count-snps 2 \
  --cgp-min-count-hmer-indels 2 \
  --cgp-min-count-non-hmer-indels 2 \
  --cgp-min-fraction-snps 0.12 \
  --cgp-min-fraction-hmer-indels 0.12 \
  --cgp-min-fraction-non-hmer-indels 0.06 \
  --cgp-min-fraction-single-strand-non-snps 0.15 \
  --cgp-min-hmer-plus-one-candidate 7 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --gvcf --p-error 0.005 \
  --optimal-coverages 50 \
  --median-coverage <median_coverage> \
  --cycle-examples-min 100000 \
  --keep-duplicates \
  --add-ins-size-channel \
  --exp-pangenome-haps haplotypes.cram
```

#### Running call_variants with Haplotype Data

Variant calling with haplotype data requires a specific model:

```bash
gs://concordanz/deepvariant/model/germline/wgs/v2.1/ultimagen-germline-wgs-solaris2-pan-hprc-v1.0-regnet-v2.1.onnx
```

or

```bash
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/wgs/v2.1/ultimagen-germline-wgs-solaris2-pan-hprc-v1.0-regnet-v2.1.onnx
```

The post-processing step remains unchanged.

### Producing a GVCF file
In case a GVCF is desired, then the commands should be modified in the following ways:

1. Add the arguments `--gvcf` and `--p-error 0.005` to the make_examples step. The p-error is an estimation of the probability of *any* base-calling error in the reads, and is used in the calculation of the reference confidence model. When these arguments are added, make_examples will output extra files with the `gvcf.tfrecord.gz` suffix.
2. When running post_process, add the argument `--gvcf_outfile output_prefix.g.vcf.gz` and provide the `gvcf.tfrecord.gz` files as input using the `--nonvariant_site_tfrecord_path` argument. The `gvcf.tfrecord.gz` files can be provided to `--nonvariant_site_tfrecord_path` either as a comma-separated list, or a text file that contains all the paths. In the latter case use the name of the ```--nonvariant_site_tfrecord_path @gvcf_records.txt```.
3. Optionally, use the `--gq-resolution` or `--gq-thresholds` arguments to reduce the output gvcf size, by binning intervals with similar GQ values together. `--gq-resolution` sets a constant difference between the bins, and `--gq-thresholds` accepts a list of specific bin thresholds, e.g. `--gq-thresholds 0,1,8,15,22` (the rounding is downwards). A value of 0 is special and results in a bin of 0.

### Calling runs of homozygosity (ROH)

Runs of homozygosity are long stretches in which the sample carries no heterozygous variant, reported as the regions tsv of `bcftools roh` from the VCF that post_process produced. This is meaningful for germline whole-genome data only.

```bash
marker_filter='FILTER="PASS" && TYPE="snp" && N_ALT=1'

# The calling regions used for variant calling, minus the blacklist.
gatk IntervalListToBed -I wgs_calling_regions.hg38.interval_list -O calling_regions.bed
bedtools subtract -a calling_regions.bed -b hg38-blacklist.v2.bed | cut -f1-3 > regions.bed

first=1
: > output_prefix.roh.tsv
while read -r line; do
  echo "$line" > region.bed
  bcftools roh --AF-dflt 0.6 --include "$marker_filter" \
      --regions-file region.bed --output-type r output_prefix.vcf.gz \
    | awk -v first=$first 'first && /^# RG/ {print} $1=="RG" {print}' >> output_prefix.roh.tsv
  first=0
done < regions.bed
```

The output is the output of `bcftools roh`: the literal `RG`, sample, chromosome, start, end, length in bp, number of markers, and the mean forward-backward phred quality. Coordinates are 1-based inclusive, so the length equals `end - start + 1`.

```
# RG	[2]Sample	[3]Chromosome	[4]Start	[5]End	[6]Length (bp)	[7]Number of markers	[8]Quality (average fwd-bwd phred score)
RG	sample_name	chr17	10274129	11420961	1146833	1119	82.4
```

### Calling variants in the mitochondrial genome (chrM)

Running on mitochondrial contig is done in the following way:

**Restrict the calling region to the mitochondrial contig.** Note that `wgs_calling_regions.hg38.interval_list` contains **no chrM interval lines** — so running it unchanged calls no mitochondrial variants at all. The same is true of the hg19/b37 calling regions.
Generate the interval list from the reference index instead of using a pre-built file, since the contig is named `chrM` in hg38 and `MT` in b37:

```bash
awk 'BEGIN{OFS="\t"} $1=="chrM" {print $1, 0, $2}' Homo_sapiens_assembly38.fasta.fai > chrM.bed
```

The contig is only 16,569 bp, so there is no need to scatter it (`num_shards` is 1): this single bed file is passed to one make_examples job.

**make_examples** is then run with lowered candidate fractions, so that low-heteroplasmy alleles become candidates:

```
tool \
  --input input_reads.cram \
  --cram-index input_reads.cram.crai \
  --bed chrM.bed \
  --output 001 \
  --reference Homo_sapiens_assembly38.fasta \
  --min-base-quality 5 \
  --min-mapq 5 \
  --cgp-min-count-snps 2 \
  --cgp-min-count-hmer-indels 2 \
  --cgp-min-count-non-hmer-indels 2 \
  --cgp-min-fraction-snps 0.05 \
  --cgp-min-fraction-hmer-indels 0.05 \
  --cgp-min-fraction-non-hmer-indels 0.05 \
  --cgp-min-fraction-single-strand-non-snps 0.06 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --optimal-coverages 70 \
  --median-coverage <median_coverage> \
  --keep-duplicates \
  --add-ins-size-channel \
  --single-strand-filter
```

**call_variants** uses a different model, with a single pass:

```bash
call_variants \
  --model model_dyn_1500.onnx \
  -T @examples.txt \
  -n 8 \
  --fp16
```

The model is available in:

```
gs://concordanz/deepvariant/model/germline/segdup_1.0/model_dyn_1500.onnx
or
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/segdup_1.0/model_dyn_1500.onnx
```

**post_process** is run with the same command as above.

**Interpreting the output.** EfficientDV genotypes diploid, so a homoplasmic mitochondrial allele is reported as `0/1` with an allele fraction close to 1.0. Assess heteroplasmy from the allele fraction rather than the genotype: it is the `FORMAT/VAF` field ("Variant allele fractions", one value per ALT allele).

### Legacy models

#### Solaris 1.0
The models described above apply to Solaris 2.0 data. For processing WGS Solaris 1.0 data we recommend the following changes to the workflow:

1. Use the following model: `s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/wgs/v1.9/ultima-usb4-amp_pcrfree-germline-model-v1.9.ckpt-420000.batch1500.onnx`
2. Run call_variants with a single pass (omit all `--second-pass-*` arguments).
