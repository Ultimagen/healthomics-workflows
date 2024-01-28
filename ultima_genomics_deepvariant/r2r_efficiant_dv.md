
## To convert the Ultima Genomics DeepVarinat Ready2Run workflow to a private workflow in AWS HealthOmics, follow the steps below:
i. Pre requisites: 
1. Ultima Genomics DeepVariant uses several ECR containers:
   - ultimagen/make_examples:latest
   - ultimagen/call_variants:latest
   - amazonlinux:2023
   - ous.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698
   - broadinstitute/gatk:4.2.6.1
2. Pull and push the Ultima Genomics DeepVariant public containers to your private ECR by following the steps here:

    a. Pull from docker hub or broad gcr into your local ecr
    b. Grant AWS HealthOmics permission to access your private ECR by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/permissions-resource.html#permissions-resource-ecr).
   
3. Import your input files into a S3 bucket.
4. Create an OmicsService role to access your resources by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/setting-up-workflows.html).

ii. Download [the workflow folder](UltimaGenomicsDV.zip) as a zipped file. You can save this zipped file locally or in a S3 bucket. 

iii. Download [the parameter template file](ultima_genomics_deepvarinat/r2r_efficient_dv_parameter_template.json). You can save this zipped file locally or in a S3 bucket.

iv. Modify and save the workflow scripts and parameter templates to meet your needs:
   - Update **tasks/globals.wdl** with the urls of you private ECR images.
   - Update the fixed parameter values in **r2r_efficient_dv.wdl** with their real s3 paths. Current links are to public buckets and transfer cost might incur by using them
