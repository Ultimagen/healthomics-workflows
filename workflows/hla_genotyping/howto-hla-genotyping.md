# Human Leukocyte Antigen (HLA) Genotyping

HLA geontyping pipline uses HLA-LA tool that identifies HLA alleles from whole genome sequencing data.
The pipeline takes a CRAM file as input and produces a txt file of the HLA alleles.

## Workflow Overview

The pipeline takes as an input BAM/CRAM file and outputs a txt file of the HLA alleles.
The steps of the pipeline are as follows: extract sample name form CRAM/BAM file, running HLA-LA.

## Requirements

1. Input CRAM/BAM files and respective indexes
2. Install HLA-LA from https://github.com/Ultimagen/HLA-LA
3. Install samtools

### Files required for the analysis (download locally)

gs://concordanz/hla/PRG_MHC_GRCh38_withIMGT_indexed.tar.gz
This is the required indexed graphs that HLA-LA creates, just to save time and skip this step which requires time and resources

#### Get sample name out of the cram file
Run the following command to extract the sample name from the CRAM file:
```
samtools view -H input.cram | grep "^@RG"

The sample name is the value of the SM tag in the read group line.
```

#### Running HLA-LA
after installing and setting up HLA-LA, run the following command:
```
optinal:
Instead of running 
cd /usr/local/bin/HLA-LA/src
../bin/HLA-LA --action prepareGraph --PRG_graph_dir ../graphs/PRG_MHC_GRCh38_withIMGT

you can download the pre-built graphs from the following link:
gs://concordanz/hla/PRG_MHC_GRCh38_withIMGT_indexed.tar.gz
extract and paste it in the graphs directory /usr/local/bin/HLA-LA/graphs/
```

Run the following command to run HLA-LA:
```
/usr/local/bin/HLA-LA/src/HLA-LA.pl \
--BAM input.cram \
--workingDir /path_to_output_dir/ \
--graph PRG_MHC_GRCh38_withIMGT \
--sampleID <sample_name> \
--maxThreads 7 \
--longReads ultimagen
```


#### Output
The output of the pipeline is a txt file containing the HLA alleles.
The output file is located at /path_to_output_dir/hla/R1_bestguess_G.txt