# Calling somatic variants from aligned cram files using efficient DV

## Overview

Efficient DV is an analysis pipeline designed to call variants from aligned cram files using the method of (DeepVariant)[https://www.nature.com/articles/nbt.4235], adapted for Ultima Genomics data. The pipeline was extended to call somatic variants.

There are three stages to the variant calling:
1. make_examples - Looks for “active regions” with potential candidates. Within these regions, it performs local assembly (haplotypes), re-aligns the reads, and defines candidate variant. Images of the reads in the vicinity of the candidates are saved as protos in a tfrecord format.
2. call_variants - Collects the images from make_examples and uses a deep learning model to infer the statistics of each variant (i.e. quality, genotype likelihoods etc.), using the TensorRT framework. It outputs sorted results in a CallVariantsOutput protos saved as tfrecords files.
3. post_process - Uses the call_variants tfrecord output. It resolves multi-allelic records and variants that overlap with indels. It annotates the variants with information such as variant type, cycle skip status and genomic intervals it appears in (e.g. exome). It filters the variants based on defined thresholds. Finally, it writes out a VCF file.

The workflow description below focuses on whole-genome somatic sequencing (tumor coverage 40-150x, normal coverage 40-100x). However, Efficient DV can support more applications (e.g. deep whole exome sequencing) with some modifications (most importantly - the model). For details on using Efficient DV for more sample preps, see the last sections of this document.

## Requirements

### Input files

The workflow takes three major inputs:

1. Reads aligned to a reference genome, for the tumor sample and a background "normal" sample (e.g. blood sample from the same patient). The reads need to be sorted and duplicate-marked, stored in CRAM format, and indexed with CRAI.
2. Reference genome files (the same reference used to align the reads). Typically hg38 is used, and is publicly available in:
```
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
or
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai
s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict
```
The calling regions interval list used in the WGS somatic pipeline excludes centromeres:
```
s3://ultimagen-workflow-resources-us-east-1/hg38/wgs_calling_regions.hg38_no_centromeres.interval_list
```
3. A model checkpoint in ONNX format
```
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/somatic/fresh_frozen/matched_normal/v1.3.1/deepvariant-ultima-somatic-wgs-model-v1.3.1-04-06.onnx
```
This model is intended for whole genome sequencing at a coverage of 40-150x for tumor and 40-100x for normal. For other applications, see the last sections of this document.

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

3. There are some [optional filtering steps](#additional-filtering-steps-recommended) that are implemented in the `ugbio_filtering` docker.

## Workflow details

The workflow is composed of three steps, as described above: make_examples, call_variants and post_process.

### Running make_examples

make_examples is typically run on a small interval in the genome determined by an input bed file. This can be used to scatter the genome into small intervals, and then parallelize the workflow across these intervals. A convenient tool to generate scattered intervals is [picard IntervalListTools](https://gatk.broadinstitute.org/hc/en-us/articles/360036897212-IntervalListTools-Picard-). A typical command to scatter intervals provided below. Note that the input to this command is Picard's interval_list format, and not the simple bed format:

```
mkdir out && \
picard \
  IntervalListTools \
  SCATTER_COUNT=100 \
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
  --input "input_reads.cram;background_reads.cram" \
  --cram-index "input_reads.cram.crai;background_reads.cram.crai" \
  --output output001 \
  --reference Homo_sapiens_assembly38.fasta \
  --bed interval001.bed \
  --somatic \
  --min-base-quality 5 \
  --min-mapq 5 \
  --cgp-min-count-snps 2 \
  --cgp-min-count-hmer-indels 2 \
  --cgp-min-count-non-hmer-indels 2 \
  --cgp-min-fraction-snps 0.03 \
  --cgp-min-fraction-hmer-indels 0.03 \
  --cgp-min-fraction-non-hmer-indels 0.03 \
  --cgp-min-fraction-single-strand-non-snps 0.15 \
  --cgp-min-hmer-plus-one-candidate 7 \
  --cgp-min-mapping-quality 5 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --optimal-coverages '200;200' \
  --median-coverage <tumor_median_coverage> \
  --background-median-coverage <background_median_coverage> \
  --single-strand-filter \
  --ignore-multi-allelic-cvos \
  --keep-duplicates \
  --add-ins-size-channel
```

The input cram files and the corresponding index files are provided to `--input` and `--cram-index`, respectively. The tumor and background crams are separated by a semicolon (note that the semicolon requires to quote the argument, in order for linux to interpret it correctly). Multiple cram files can be provided, separated by commas.

The `--output` argument is the prefix for the output files (including tfrecords).

The argument `optimal-coverages` (and the related argument `cap-at-optimal-coverage`) determine how reads are internally downsampled before image generation. The same values are used in the training of the model and the inference. Hence, their values are tightly linked to which model is used. The argument `add-ins-size-channel` adds a channel with the length of the insertion, and should also be aligned with the model. `--median-coverage` and `--background-median-coverage` provide `make_examples` with the actual sequencing depth of the tumor and normal samples, respectively, used alongside `--optimal-coverages` to control read downsampling.

The `single-strand-filter` reduces the number of candidates, therby reducing compute costs, at a negligible effect on recall.

#### Running somatic-on-germline
In somatic variant calling, complex somatic variants near germline variants can pose challenges during read alignment and variant calling. To address this, it can be beneficial to "correct" the reference genome by incorporating germline variants during the read-alignment step of `make_examples`. This simplifies the generated image and improves the accuracy of likelihood calculations during `call_variants`.

To enable this, you can provide a germline VCF to the workflow using the `--region-haplotypes-vcf` argument. This VCF can be generated through a separate germline variant-calling workflow on the normal sample.

It’s important to note that using `--region-haplotypes-vcf` influences only the likelihood computations for the variant calls. The REF and ALT fields in the output VCF remain unchanged.

### Running call_variants

The call_variants step combines the tfrecords from all make_examples jobs. List the tfrecord files, one per line, in a text file, and run from within the call_variants docker:
```bash
printf '%s\n' 001.tfrecord.gz 002.tfrecord.gz ... > examples.txt

call_variants \
  --model model.onnx \
  -T @examples.txt \
  -n 8 \
  --fp16
```

`model.onnx` is the somatic model listed above. `-n` sets the number of threads used for reading the input files. call_variants writes its output, `*_cvo.gz` files, to the current directory (use `-o <folder>` to change it). These files are the input to post_process.

### Running post_process:

post_process uses the output of call_variants (the `*_cvo.gz` files), and generates a vcf file. A typical post_process command from within the docker will look like (same docker as make_examples):
```
ug_postproc \
  --infile shard001_cvo.gz,shard002_cvo.gz,... \
  --ref Homo_sapiens_assembly38.fasta \
  --outfile output_prefix.vcf.gz \
  --consider_strand_bias \
  --flow_order TGCA \
  --annotate \
  --bed_annotation_files LCR-hs38.bed,mappability.0.bed,hmers_7_and_higher.bed,ug_hcr.bed,gnomad.common.indel.annotation.sort.bed \
  --ignore_multi_allelic_cvos \
  --qual_filter 1 \
  --filter \
  --filters_file filters.txt \
  --consider_bg_fields \
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

The `--ignore_multi_allelic_cvos` flag instructs post-processing to skip resolving multi-allelic records from call_variants output, which is the standard behavior for somatic calling.

If `--filter` argument is used, the vcf will be filtered based on the criteria in `--filters_file`. In this file, each filter is composed of two lines, the first is the filter name (which will appear in FILTER column of the vcf), and the second is the expression for the filter. The syntax of the expression follows [JEXL filtering expressions](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions). Below is an example of a filters_file that is typically used for WGS somatic calling:
```
LowQualInExome
QUAL < 20 and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and vc.hasAttribute('EXOME')
LowQual
QUAL < 20 and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and not vc.hasAttribute('EXOME')
LowQual
QUAL < 20 and VARIANT_TYPE=='non-h-indel' and not vc.isFiltered()
LowQual
QUAL < 14 and VARIANT_TYPE=='snp' and not vc.isFiltered()
LargeDeletion
REFLEN > 220 and vc.isFiltered()
```

dbSNP data can be downloaded from: gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
or s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

`ug_postproc` can also be used to improve the VAF estimates of INDELs. Note: this is not recommended when the callset contains more than 30K indels. To enable it, append the following parameters to `ug_postproc`:

```
--fix_allele_coverage 
--fix_allele_indels_only 
--fix_allele_crams <tumor.cram>;<normal.cram>
```


### Allele-frequency ratio filtering

The pipeline applies a filter on the allele-frequency ratio between tumor and normal. This is a safe measure to hard-filter variants with significant appearance in the normal sample..

Specifically, the filter removes any SNV or Non-homopolymer indel variants where the ratio of allele frequency between the tumor and the normal variant is below the threshold (default: 10).

To apply the filter, run:

```
  # Run in the ugbio_filtering docker
  filter_low_af_ratio_to_background --af_ratio_threshold 10 \
                                    --new_filter "LowAFRatioToBackground"  \
                                    <input_vcf> <output_vcf>
```

## More applications

Efficient DV can support somatic calling in more scenarios, with the following modifications:

### WGS somatic calling from FFPE 
Calling somatic variants from a tumor FFPE sample at a coverage of 100x and a normal sample at a coverage of 40-80x. The modifications required to use this application are:
1. Change the model of call_variants to:
```
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/somatic/wgs/ffpe/v1.8/deepvariant-ultima-somatic-wgs-ffpe-model-v1.8.onnx
```
2. Add the following argument to the make_examples step:
```
--channels hmer_deletion_quality,hmer_insertion_quality,non_hmer_insertion_quality
```
3. Change the value of the SNP quality for filtering:
```
LowQual
QUAL < 10 and VARIANT_TYPE=='snp' and not vc.isFiltered()
```
4. Update `--median-coverage` and `--background-median-coverage` in the make_examples step to match the actual tumor and normal sample coverages.

### WGS somatic calling for SNV signature detection 
Calling somatic variants for SNV signature detection from a coverage of 40x of tumor and normal samples. The same model as in WGS somatic calling is used, but with slight modifications to the coverages and candidate generation thresholds in the make_examples step:
```
  --cgp-min-fraction-snps 0.075 \
  --cgp-min-fraction-hmer-indels 1.1 \
  --cgp-min-fraction-non-hmer-indels 1.1 \
  --optimal-coverages "40;40" \
  --median-coverage <tumor_median_coverage> \
  --background-median-coverage <background_median_coverage> \
```

### Somatic calling from deep sequencing of whole exome
Calling from deep whole exome sequencing, at a coverage of 500x for tumor and at least >120x for normal. The required modifications are:
1. To the make_examples step:
```
  --cgp-min-fraction-snps 0.02 \
  --cgp-min-fraction-hmer-indels 0.02 \
  --cgp-min-fraction-non-hmer-indels 0.02 \
  --optimal-coverages "500;120" \
  --median-coverage <tumor_median_coverage> \
  --background-median-coverage <background_median_coverage> \
  --max-reads-per-region 6500 \
  --prioritize-alt-supporting-reads \
  --normalize-strand-bias \
  --strand-bias-threshold-to-normalize 0.1,0.35
```
Also, make sure you read variants only in exomic intervals.
2. To the call_variants model:
```
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/somatic/wes/deepvariant-ultima-somatic-wes-model-v0.1.ckpt-120000.onnx
```

3. We recommend the following additional parameters to the [optional filtering steps](#additional-filtering-steps-recommended): 
```
      --af_ratio_threshold_h_indels 2 \
      --tumor_vaf_threshold_h_indels 0.15
```


This addition filters out homopolymer indels where the ratio of allele frequency is below 2 and the tumor VAF of the variant is below 0.15.  
