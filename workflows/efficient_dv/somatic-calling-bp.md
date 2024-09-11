# Somatic mutation calling - best practices - v1.11.2

Last updated: Jun 04, 20242 

## Table of contents
1. [Running somatic variant calling pipeline](#running-somatic-variant-calling-pipeline)
2. [Available use-cases and requirements](#available-use-cases-and-requirements)
   
   a. [Sample characterization, fresh frozen, PCR-free WGS libraries](#sample-characterization-fresh-frozen-pcr-free-wgs-libraries)
   
   b. [Signature detection for MRD, fresh frozen, PCR-free WGS libraries](#signature-detection-for-mrd-fresh-frozen-pcr-free-wgs-libraries)
   
   c. [Sample characterization, FFPE, amplified WGS libraries](#sample-characterization-ffpe-amplified-wgs-libraries)
   
   d. [Sample characterization, amplified WES libraries](#sample-characterization-amplified-wes-libraries)
3. [Pipeline output](#pipeline-output)
4. [Evaluation process](#evaluation-process)

<a name="running-somatic-variant-calling-pipeline"></a>
## Running somatic variant calling pipeline
Ultima Genomics variant calling pipeline is based on our re-write of DeepVariant pipeline and is freely available [here](https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/efficient_dv) 

We provide the pipeline in the WDL format and sets of input parameters for different use cases (see below). See [efficient_dv.md](efficient_dv.md)  for the parameter documentation.

We recommend using AWS Healthomics For easy integration of the pipeline use these [instructions](https://github.com/Ultimagen/healthomics-workflows/blob/main/README.md).

The WDL files we provide should be also compatible with Cromwell engine running on GCP.

For the users that do not use WDL-based pipelines we provide a [description of the pipeline](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/efficient_dv/howto-somatic-calling-efficient-dv.md). We also list the necessary dockers and commands required to implement the pipeline on HPC of choice. Note that GPU instances are required for efficient inference. 

<a name="available-use-cases-and-requirements"></a>
## Available use-cases and requirements
<a name="sample-characterization-fresh-frozen-pcr-free-wgs-libraries"></a>
### Sample characterization, fresh frozen, PCR-free WGS libraries
Samples required: tumor sample at coverage 40x-150x, normal sample at coverage 40x-100x

Variants called at high confidence: SNVs with allele frequency > 5%, Indels  with allele frequency > 10%

Parameter set: [workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-T_N-for-sample-characterization-v1_0-BC-4_9--40-150x_40-100x-.json](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-T_N-for-sample-characterization-v1_0-BC-4_9--40-150x_40-100x-.json)

<a name="signature-detection-for-mrd-fresh-frozen-pcr-free-wgs-libraries"></a>
### Signature detection for MRD, fresh frozen, PCR-free WGS libraries
Samples required: tumor sample at coverage 40x, normal sample at coverage 40x

Variants called at high confidence: SNVs with allele frequency > 10%

Parameter set: [workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-T_N-for-SNV-signature-detection-v1_0--40x_40x-.json](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-T_N-for-SNV-signature-detection-v1_0--40x_40x-.json)

The pipeline can be used for tumor samples with coverage higher than 40x as a much faster and cheaper alternative to the full sample characterization

<a name="sample-characterization-ffpe-amplified-wgs-libraries"></a>
### Sample characterization, FFPE, amplified WGS libraries
Samples required: tumor sample at coverage >100x, normal sample at coverage 40x-80x. Only tumor sample is expected to be FFPE, normal sample is expected to be PCR-free WGS. 

Variants called: SNVs with allele frequency > 5%, indels with allele frequency > 10%

Parameter set: [workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-FFPE-T_N-for-sample-characterization-v1_3--100x_40-80x-.json](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/efficient_dv/input_templates/efficient_dv_template-WGS-somatic-FFPE-T_N-for-sample-characterization-v1_3--100x_40-80x-.json) 

<a name="sample-characterization-amplified-wes-libraries"></a>
### Sample characterization, amplified WES libraries
Samples required: tumor sample at coverage >500x, normal sample at coverage >120x.  

Variants called: SNVs with allele frequency > 5%, indels with allele frequency > 5%

Parameter set: [workflows/efficient_dv/input_templates/efficient_dv_template-WES-somatic-T_N-for-deep-sample-characterization-v0_1--500x_gt120x-.json](https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/efficient_dv/input_templates/efficient_dv_template-WES-somatic-T_N-for-deep-sample-characterization-v0_1--500x_gt120x-.json)

<a name="pipeline-output"></a>
## Pipeline output
The pipeline produces output in the VCF format. Note that by default, all variant candidates will appear in the VCF including germline candidates and sequencing artefacts. The called somatic variants will be labeled by PASS in the FILTER field of the VCF. 

**Notes:**

1. The pipeline does not reliably call long indels (longer than 30 bases). SV calling pipeline should be used to call those variants. 

2. The estimates of AD and VAF are naive and are biased towards lower values for indels. 

3. We do not try to call variants with AF < 3% for WGS and AF < 2% for WES 

<a name="evaluation-process"></a>
## Evaluation process
We suggest using standard tool [rtg vcfeval](https://github.com/RealTimeGenomics/rtg-tools) to perform comparison of the callset to the ground truth.  

We suggest to exclude from the evaluation variants that fall in the UG-LCR, the region where we do not call variants confidently [UG-LCR](https://github.com/Ultimagen/healthomics-workflows/blob/main/docs/ug_hcr.md). ug_lcr as well as the complement ug_hcr annotation bed files can be downloaded from `s3://ultimagen-workflow-resources-us-east-1/hg38/UG-High-Confidence-Regions/v2.1.2/`. 

### Example command: 


    rtg vcfeval
    -b <baseline>.vcf.gz
    -c <UG>.vcf.gz
    -o <output_dir>
    -t Homo_sapiens_assembly38.fasta.sdf
    --decompose
    --squash-ploidy
    --sample=<gt_sample>,<ug_sample>
    --bed-regions ug_hcr.bed
    -f QUAL
    