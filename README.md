# ultimagenomics-healthomics-workflows
UltimaGenomics repository for workflows compatible with AWS HealthOmics

## Introduction
i. UltimaGenomics offers our pipelines as Ready2Run workflows on AWS HealthOmics. Ready2Run workflows enable you to run these pipelines on AWS HealthOmics by simply bringing your data.

ii. For more flexibility or to modify the workflow, you can convert the Ready2Run workflow to a private workflow that can be run on AWS HealthOmics.

iii. The cost to run the workflow will now be based on the compute and run storage used during the private workflow run instead of a fixed price.

iv. For more questions about this workflow, please contact healthomics.support@ultimagen.com.

## To convert the Ultima Genomics DeepVarinat Ready2Run workflow to a private workflow in AWS HealthOmics, follow the steps below:
i. Pre requisites: 
1. Ultima Genomics DeepVariant uses several ECR containers:
   - ultimagenomics/make_examples:latest
   - ultimagenomics/call_variants:latest
   - omics-shared-amazonlinux:2023
   - omics-gitc:latest
   - omics-gatk:latest
2. Pull and push the Ultima Genomics DeepVariant public containers to your private ECR by following the steps here:

    a. Grant AWS HealthOmics permission to access your private ECR by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/permissions-resource.html#permissions-resource-ecr).
   
3. Import your input files into a S3 bucket.
5. Create an OmicsService role to access your resources by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/setting-up-workflows.html).

ii. Download [the workflow folder](ultima_genomics_deepvarinat/UltimaGenomicsDV.zip) as a zipped file. You can save this zipped file locally or in a S3 bucket. 

iii. Download [the parameter template file](ultima_genomics_deepvarinat/r2r_efficient_dv_parameter_template.json). You can save this zipped file locally or in a S3 bucket.

iv. Modify and save the workflow scripts and parameter templates to meet your needs:
   - Update **tasks/globals.wdl** with the urls of you private ECR images.
   - Update the fixed parameter values in **r2r_efficient_dv.wdl** with their real s3 paths.
   - Do other modifications as needed.

v. Create a private workflow in HealthOmics by following one of the two options below:
1. From the CLI:
 ~~~
$ aws omics create-workflow \
--name <workflow_name> \
--main r2r_efficient_dv.wdl \
--definition-zip <s3_uri_for_zipped_file> \
--parameter-template <s3_uri_for_parameter_template_file> \
--accelerators GPU
 ~~~
2. From the console:
    
    a. Click on **Private Workflows** from the left pane.
    
    b. Click on **Create Workflow** on the Workflows list.
    
    c. Follow the instructions on the console to create your workflow.

vi. Run your workflow by following one of the two options below:
1. From the CLI:
 ~~~
$ aws omics start-run \
--workflow-id <workflow_id> \
--role-arn <service_role_arn> \
 --output-uri <s3_uri_for_output_folder> \
 --parameters <s3_uri_or_local_parameters_file>
 --name <run_name> \
 --retention-mode REMOVE
 ~~~
2. From the console:
   
   a. Click **Private Workflows** from the left pane.

   b. Click the **Workflow ID** from the Workflows list.

   c. Click **Create Run** and enter the run information.
