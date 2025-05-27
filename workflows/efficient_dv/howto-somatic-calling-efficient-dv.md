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
gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list
```
3. A model checkpoint in ONNX format
```
onnxFileName = gs://concordanz/deepvariant/model/somatic/wgs/model_dyn_1000_200_221_9_2950000.onnx
```
This model is intended for whole genome sequencing at a coverage of 40-150x for tumor and 40-100x for normal. For other applications, see the last sections of this document.

### Dockers and hardware requirements

The Efficient DV analysis pipeline is split into two docker images:

1. `make_examples` docker - contains binaries for the make_examples and post_process steps. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/make_examples:3.1.2
or
337532070941.dkr.ecr.us-east-1.amazonaws.com/make_examples:3.1.2
```
2. `call_variants` docker - contains binaries for the call_variants step. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/call_variants:2.2.2
or
337532070941.dkr.ecr.us-east-1.amazonaws.com/call_variants:2.2.2
```

The make_examples and post_process steps are run on a single CPU. make_examples requires up to 2 GB of memory for each thread. post_process requires 8 GB of memory and runs on a single thread.

call_variants runs on a machine which contains a single GPU, such as nvidia-p100, or nvidia-v100. It also uses multiple CPUs for multi-threaded decompression of input tfrecord files. The required memory is 8 GB plus 1 GB for each decompression thread.


## Workflow details

The workflow is composed of three steps, as described above: make_examples, call_variants and post_process.

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
  --cgp-min-mapping-quality 5 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --optimal-coverages '200;200' \
  --single-strand-filter \
  --keep-duplicates \
  --add-ins-size-channel
```

The input cram files and the corresponding index files are provided to `--input` and `--cram-index`, respectively. The tumor and background crams are separated by a semicolon (note that the semicolon requires to quote the argument, in order for linux to interpret it correctly). Multiple cram files can be provided, separated by commas.

The `--output` argument is the prefix for the output files (including tfrecords).

The argument `optimal-coverages` (and the related argument `cap-at-optimal-coverage`) determine how reads are internally downsampled before image generation. The same values are used in the training of the model and the inference. Hence, their values are tightly linked to which model is used. The argument `add-ins-size-channel` adds a channel with the length of the insertion, and should also be aligned with the model.

The `single-strand-filter` reduces the number of candidates, therby reducing compute costs, at a negligible effect on recall.

#### Running somatic-on-germline
In somatic variant calling, complex somatic variants near germline variants can pose challenges during read alignment and variant calling. To address this, it can be beneficial to "correct" the reference genome by incorporating germline variants during the read-alignment step of `make_examples`. This simplifies the generated image and improves the accuracy of likelihood calculations during `call_variants`.

To enable this, you can provide a germline VCF to the workflow using the `--region-haplotypes-vcf` argument. This VCF can be generated through a separate germline variant-calling workflow on the normal sample.

It’s important to note that using `--region-haplotypes-vcf` influences only the likelihood computations for the variant calls. The REF and ALT fields in the output VCF remain unchanged.

### Running call_variants

The call_variants step combines the tfrecords from all make_examples jobs. The arguments to the call_variants step are provided as an `.ini`-formatted file. A typical file will look like:
```
[RT classification]
onnxFileName = model/somatic/fresh_frozen/matched_normal/v1.3/wgs_somatic_matched_normal_v1.3.onnx
useSerializedModel = 1
trtWorkspaceSizeMB = 2000
numInferTreadsPerGpu = 2
useGPUs = 1
gpuid = 0

[debug]
logFileFolder = .

[general]
tfrecord = 1
compressed = 1
outputInOneFile = 0
numUncomprThreads = 8
uncomprBufSizeGB = 1
outputFileName = call_variants
numConversionThreads = 2
numExampleFiles = 40

exampleFile 1 = input_dir/001.tfrecord.gz
exampleFile 2 = input_dir/002.tfrecord.gz
exampleFile 3 = input_dir/003.tfrecord.gz
...
```

The last part of the `ini` file is a list with the paths to all tfrecord files, and the first `numExampleFiles` files are used.
If `useSerializedModel` is set to 1, then the programs searches for a onnx-serialized file, with the same name as the onnx file and .serialized suffix. If the file is not found, it generates it. The serialized file can be re-used across different runs on the same platform (TRT version, GPU type etc.). The number of GPUs and the number of tfrecord.gz uncompression threads (runs on CPU) can be modified using the `useGPUs` and `numUncomprThreads` arguments, respectively.

Once the `ini` file is ready, call_variants can be invoked from within the docker using:
```
call_variants --param params.ini --fp16
```


### Running post_process:

post_process uses the output of call_variants, `call_variants.tfrecord.gz`, and generates a vcf file. A typical post_process command from within the docker will look like:
```
ug_postproc \
  --infile call_variants.1.gz,call_variants.2.gz,... \
  --ref Homo_sapiens_assembly38.fasta \
  --outfile output_prefix.vcf.gz \
  --consider_strand_bias \
  --flow_order TGCA \
  --annotate \
  --bed_annotation_files exome.twist.bed,ug_hcr.bed,... \
  --qual_filter 1 \
  --filter \
  --filters_file filters.txt \
  --consider_bg_fields \
  --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf
```

The `bed_annotation_files` are comma-separated list of bed files that are used to annotate the vcf. Description of the bed file is recorded in the INFO field of the vcf. The description can be provided in one of two ways:

Using a `##INFO` in the header of the bed file. For example:
```
##INFO=<ID=EXOME,Number=1,Type=String,Description="Genomic Region Annotation: In the exome (gs://concordanz/hg38/annotation_intervals/exome.twist.bed)">
chr1    69090   70008   TRUE
chr1    450739  451678  TRUE
```

If `##INFO` is not present in the bed file, then a json file with the same name as the bed file (e.g. exome.twist.bed.json) should be provided. The json will contain the ID, Type and Description:
```
{
 "ID": "EXOME",
 "Type": "String",
 "Description": "Genomic Region Annotation: In the exome (gs://concordanz/hg38/annotation_intervals/exome.twist.bed)"
}
```

If `--filter` argument is used, the vcf will be filtered based on the criteria in `--filters_file`. In this file, each filter is composed of two lines, the first is the filter name (which will appear in FILTER column of the vcf), and the second is the expression for the filter. The syntax of the expression follows [JEXL filtering expressions](https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions). Below is an example of a filters_file that is typically used. Note that these filters use the EXOME attribute, which is an annotation that was added using a bed file.
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

`ug_post_processing` can also be used to improve the VAF estimates of INDELs. This is recommended only run if the callset contains less than 30K indels. 
To do this, append the following parameters

```
--fix_allele_coverage 
--fix_allele_indels_only 
--fix_allele_crams <tumor.cram>;<normal.cram>
```


### Additional filtering steps recommended. 

We recommend applying an additional post-processes filter. The filter is applied to the allele-frequency ratio between tumor and normal. It filters out loci with germline variants that changed allele-frequency in the tumor, as the model was not trained to filter them, making this inconsistent with some other somatic calling pipelines.

The process to applying this filter is:

```
bcftools filter input.vcf.gz -e '(VARIANT_TYPE="snp" || VARIANT_TYPE="non-h-indel") && (AD[0:1]/DP)/(BG_AD[0:1]/BG_DP) < 10' -s "LowAFRatioToBackground" -m "+" -Oz -o  output.vcf.gz
bcftools index -t output.vcf.gz
```

## More applications

Efficient DV can support somatic calling in more scenarios, with the following modifications:

### WGS somatic calling from FFPE 
Calling somatic variants from a tumor FFPE sample at a coverage of 100x and a normal sample at a coverage of 40-80x. The modifications required to use this application are:
1. Add the argument `--channels` to the make_examples steps with the value of:
```
--channels hmer_deletion_quality,hmer_insertion_quality,non_hmer_insertion_quality,soft_clips
```
This adds a channel with information about soft-clipped reads.
2. Change the model of call_variants to:
```
gs://concordanz/deepvariant/model/somatic/wgs/ffpe/deepvariant-ultima-somatic-wgs-ffpe-model-v1.3.ckpt-890000.onnx
```
3. Change the value of the SNP quality for filtering:
```
LowQual
QUAL < 10 and VARIANT_TYPE=='snp' and not vc.isFiltered()
```

### WGS somatic calling for SNV signature detection 
Calling somatic variants for SNV signature detection from a coverage of 40x of tumor and normal samples. The same model as in WGS somatic calling is used, but with slight modifications to the coverages and candidate generation thresholds in the make_examples step:
```
  --cgp-min-fraction-snps 0.075 \
  --cgp-min-fraction-hmer-indels 1.1 \
  --cgp-min-fraction-non-hmer-indels 1.1 \
  --optimal-coverages "40;40" \
```

### Somatic calling from deep sequencing of whole exome
Calling from deep whole exome sequencing, at a coverage of 500x for tumor and at least >120x for normal. The required modifications are:
1. To the make_examples step:
```
  --cgp-min-fraction-snps 0.02 \
  --cgp-min-fraction-hmer-indels 0.02 \
  --cgp-min-fraction-non-hmer-indels 0.02 \
  --optimal-coverages "500;120" \
  --max-reads-per-region 6500 \
  --prioritize-alt-supporting-reads
```
Also, make sure you read variants only in exomic intervals.
2. To the call_variants model:
```
gs://concordanz/deepvariant/model/somatic/wes/deepvariant-ultima-somatic-wes-model-v0.1.ckpt-120000.onnx
```