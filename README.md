# healthomics-workflows
UltimaGenomics repository for workflows compatible with AWS HealthOmics

## Introduction

1.	Ultima Genomics offers pipelines as Ready2Run workflows on AWS HealthOmics. Ready2Run workflows enable you to run these pipelines on AWS HealthOmics by simply bringing your data. For more flexibility such as the use of larger file sizes or changing the reference genome, you can convert Ready2Run workflows to private workflows by following the steps in this repository. Once the Ready2Run workflow is converted to a private workflow, the cost to run the workflow will now be based on the compute and run storage used during the private workflow.

2.	Ultima Genomics also shares pipelines that has been modified to run as private workflows on AWS HealthOmics in this repository. You can follow the directions in this repository to create and run a private workflow on AWS HealthOmics. The instractions includes localizing resources, deploying workflow and creating a run.

3.	For more questions about these workflows, please contact healtomics.support@ultimagen.com.

## To localize workflow resources to create a private workflow in AWS HealthOmics, follow the steps below:
i. Pre requisites: 
1. Ultima Genomics workflows uses several ECR containers, they are listed under each workflow tasks\globals.wdl.

2. Pull and push the required public containers to your private ECR by following the steps here:

    a. Pull from docker hub or broad gcr into your local ecr
    b. Grant AWS HealthOmics permission to access your private ECR by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/permissions-resource.html#permissions-resource-ecr).
   
3. Import your input files into a S3 bucket.
4. Create an OmicsService role to access your resources by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/setting-up-workflows.html).

ii. Download the workflow folder as a zipped file, this should include main wdl file on the top level folder, tasks folders and <workflow_name>_params.json . You can save this zipped file locally or in a S3 bucket. 

iii. Download the parameter template for your desired use case from input_templates folder. You can save this file locally or in a S3 bucket.

iv. Modify and save the workflow scripts and parameter templates to meet your needs:
   - Update **tasks/globals.wdl** with the urls of you private ECR images.
   - Update the template's parameter values with their local s3 paths. Current links are to public buckets and transfer cost might incur by using them

## Deploying Private Workflow
Once the workflow resources have been deployed into locally (see instructions per workflow), user can create private workflow on AWS HealthOmics
1. Create a private workflow in HealthOmics by following one of the two options below:

i. From the CLI:
 ~~~
$ aws omics create-workflow \
    --name <workflow_name> \
    --main <main_wdl_file> \  # in case there is more than one wdl file, the main one is the one named after the directory
    --definition-zip fileb://<path_to_local_zip> \
    --parameter-template file://<path_to_parameters_definition_json> \
    --accelerators GPU
 ~~~
ii. From the console:
    
    a. Click on **Private Workflows** from the left pane.
    
    b. Click on **Create Workflow** on the Workflows list.
    
    c. Follow the instructions on the console to create your workflow.
        Note: if top level has more than one wdl, you should define "Main workflow definition file path" - r2r_efficient_dv.wdl


2. Run your workflow by following one of the two options below:
i. From the CLI:
 ~~~
$ aws omics start-run \
    --workflow-id <workflow_id> \
    --role-arn <service_role_arn> \
    --output-uri <s3_uri_for_output_folder> \
    --parameters file://<path_to_local_parameters_file> \
    --name <run_name> \
    --retention-mode REMOVE
 ~~~
ii. From the console:
   
   a. Click **Private Workflows** from the left pane.

   b. Click the **Workflow ID** from the Workflows list.

   c. Click **Create Run** and enter the run information.
