# Calling variants from aligned cram files using efficient DV

## Overview

Efficient DV is an analysis pipeline designed to call variants from aligned cram files using the method of (DeepVariant)[https://www.nature.com/articles/nbt.4235], adapted for Ultima Genomics data. There are three stages to the variant calling:
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
```
3. A model checkpoint in ONNX format:
```
gs://concordanz/deepvariant/model/germline/v1.3/model.ckpt-890000.dyn_1500.onnx
```

or 
```
s3://ultimagen-workflow-resources-us-east-1/deepvariant/model/germline/v1.3/model.ckpt-890000.dyn_1500.onnx
```

### Dockers and hardware requirements

The Efficient DV analysis pipeline is split into two docker images:

1. `make_examples` docker - contains binaries for the make_examples and post_process steps. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/make_examples:edv_2.1.1_b0ca4ece
or
337532070941.dkr.ecr.us-east-1.amazonaws.com/make_examples:edv_2.1.1_b0ca4ece
```
2. `call_variants` docker - contains binaries for the call_variants step. Can be found in:
```
us-central1-docker.pkg.dev/ganymede-331016/ultimagen/call_variants:edv_2.1.1_b0ca4ece
or
337532070941.dkr.ecr.us-east-1.amazonaws.com/make_examples:edv_2.1.1_b0ca4ece
```

The make_examples and post_process steps are run on a single CPU. make_examples requires up to 4 GB of memory for each thread. post_process requires 8 GB of memory and runs on a single thread.

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
  --input input_reads.cram \
  --cram-index input_reads.cram.crai \
  --bed interval001.bed \
  --output 001 \
  --reference Homo_sapiens_assembly38.fasta \
  --min-base-quality 5 \
  --min-mapq 5 \
  --progress \
  --cgp-min-count-snps 2 \
  --cgp-min-count-hmer-indels 2 \
  --cgp-min-count-non-hmer-indels 2 \
  --cgp-min-fraction-snps 0.12 \
  --cgp-min-fraction-hmer-indels 0.12 \
  --cgp-min-fraction-non-hmer-indels 0.06 \
  --cgp-min-mapping-quality 5 \
  --max-reads-per-region 1500 \
  --assembly-min-base-quality 0 \
  --gzip-output \
  --no-realigned-sam \
  --interval-nreads 10000 \
  --optimal-coverages "50" --cap-at-optimal-coverage \
  --add-ins-size-channel --add-proxy-support-to-non-hmer-insertion --pragmatic
```

The input cram files and the corresponding index files are provided to `--input` and `--cram-index`, respectively. Multiple cram files can be provided as a comma separated list.

The `--output` argument is the prefix for the output files (including tfrecords).

The program will output a sam file with the re-aligned reads unless the argument `--no-realigned-sam` is provided. Note that these files are very large, so provide a large disk space if you want to save the re-aligned reads.

### Running call_variants

The call_variants step combines the tfrecords from all make_examples jobs. The arguments to the call_variants step are provided as an `.ini`-formatted file. A typical file will look like:
```
[RT classification]
onnxFileName = model/germline/v1.3/model.ckpt-890000.dyn_1500.onnx
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
call_variants --param params.ini
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


### Producing a GVCF file
In case a GVCF is desired, then the commands should be modified in the following ways:

1. Add the arguments `--gvcf` and `--p-error 0.005` to the make_examples step. The p-error is an estimation of the probability of *any* base-calling error in the reads, and is used in the calculation of the reference confidence model. When these arguments are added, make_examples will output extra files with the `gvcf.tfrecord.gz` suffix.
2. When running post_process, add the argument `--gvcf_outfile output_prefix.g.vcf.gz` and provide the `gvcf.tfrecord.gz` files as input using the `--nonvariant_site_tfrecord_path` argument. The `gvcf.tfrecord.gz` files can be provided to `--nonvariant_site_tfrecord_path` either as a comma-separated list, or a text file that contains all the paths. In the latter case use the name of the ```--nonvariant_site_tfrecord_path @gvcf_records.txt```.


## Debugging tfrecords using dvtools
The make_examples code also has a handy utility called `dvtools` to view the data in the tfrecord files (without the images). It can accept a tfrecord.gz file, and output a vcf with the records. You can use it in the following way:
```    
docker run -v <path mapping> <docker name> \
  dvtools --infile debug.tfrecord.gz \
  --filetype dv --op vcf \
  --outfile debug.dvtools.vcf
```
