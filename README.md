# healthomics-workflows
UltimaGenomics repository for workflows compatible with AWS HealthOmics

# Table of Contents
1. [Introduction](#Introduction)
2. [Terraform Modules for Environment Setup](#Terraform-Modules-for-Environment-Setup)
3. [Deploying Private Workflow](#Deploying-Private-Workflow)
4. [Running Private Workflow](#Running-Private-Workflow)
5. [Support Tools](#Support-Tools)
## Introduction

1.	Ultima Genomics offers pipelines as Ready2Run workflows on AWS HealthOmics. Ready2Run workflows enable you to run these pipelines on AWS HealthOmics by simply bringing your data. For more flexibility such as the use of larger file sizes or changing the reference genome, you can convert Ready2Run workflows to private workflows by following the steps in this repository. Once the Ready2Run workflow is converted to a private workflow, the cost to run the workflow will now be based on the compute and run storage used during the private workflow.

2.	Ultima Genomics also shares pipelines that has been modified to run as private workflows on AWS HealthOmics in this repository. You can follow the directions in this repository to create and run a private workflow on AWS HealthOmics.
  
3.	Each workflow folder contains the following:
    - required wdl file\s
    - HowTo documentation that details the workflow flow and how to run it externally of wdl
    - documentation of the wdl inputs and outputs
    - json that list the parameters for creating workflow
    - folder with optional input templates with default parameters for the wdl
    - folder with the different tasks the wdl is running

4. The instructions below include localizing resources, deploying workflow and creating a run.
5. For more questions about these workflows, please contact healthomics.support@ultimagen.com.

## Terraform Modules for Environment Setup
This repository includes [Terraform modules](terraform/healthomics-wdl-env)
to simplify the deployment of AWS resources required to run Ultima Genomics WDL workflows on AWS HealthOmics.
These modules automate the creation of necessary infrastructure, such as S3 buckets, IAM roles, and permissions,
ensuring a consistent and secure environment for workflow execution.
We strongly recommend using it for building your AWS omics env.
Feel free to make customizations as per your request

Key features:
- Automated provisioning of AWS resources for HealthOmics workflows
- Easy integration with existing workflow deployment steps and version management

Refer to the [terraform/healthomics-wdl-env/README.md](terraform/healthomics-wdl-env/README.md) for detailed usage instructions and module configuration options.

## Deploying Private Workflow
### To localize workflow resources and create a private workflow in AWS HealthOmics you can:
- [use this script](scripts/healthomics_wf/create_healthomics_workflow.py) that localizes the necessary resources and create a private workflow on AWS HealthOmics. After running it your workflow [is ready to run](#running-private-workflow).
- Alternatively, follow the steps below to do it manually:
#### To localize workflow resources to create a private workflow in AWS HealthOmics, follow the steps below:
i. Pre requisites: 
1. Ultima Genomics workflows uses several ECR containers, they are listed under each workflow tasks\globals.wdl.

2. Pull and push the required public containers to your private ECR by following the steps:

    a. Pull from docker hub or broad gcr into your local ecr
     ~~~
     docker pull <hub_username>/<image_name>:<tag> #the docker as it appear on globals.wdl
     docker tag <hub_username>/<image_name>:<tag> <your_aws_account_id>.dkr.ecr.<region>.amazonaws.com/<repository_name>:<tag>
     aws ecr get-login-password --region <region> | docker login --username AWS --password-stdin <your_aws_account_id>.dkr.ecr.<region>.amazonaws.com
     docker push <your_aws_account_id>.dkr.ecr.<region>.amazonaws.com/<repository_name>:<tag> #if repository doesn't exist, you will need to create it first
     ~~~
   
    b. Grant AWS HealthOmics permission to access your private ECR by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/permissions-resource.html#permissions-resource-ecr).
   
3. Import your input files into a S3 bucket.
4. Create an OmicsService role to access your resources by following the instructions [here](https://docs.aws.amazon.com/omics/latest/dev/setting-up-workflows.html).

ii. Download the workflow folder as a zipped file, this should include main wdl file on the top level folder, tasks folders and <workflow_name>_params.json . You can save this zipped file locally or in a S3 bucket. 

iii. Download locally the parameter template for your desired use case from input_templates folder. 

iv. Modify and save the workflow scripts and parameter templates to meet your needs:
   - Update **tasks/globals.wdl** with the urls of you private ECR images.
   - Update the template's parameter values with their local s3 paths. Current links are to public buckets and transfer cost might incur by using them

Once the workflow resources have been deployed into locally (see instructions per workflow), user can create private workflow on AWS HealthOmics
#### Create a private workflow in HealthOmics by following one of the two options below:

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
       - Define "Main workflow definition file path" as <workflow_name>.wdl file

## Running Private Workflow
In case you chose to build your environment with the healthomics-wdl-env terraform module, and you've deployed the workflow by running the script, you can run your workflow easily using [this script](scripts/healthomics_wf/invoke_healthomics_run.py) that invokes the StartOmicsRun lambda.

Alternatively, follow the steps below to do it manually:
#### Run your workflow by following one of the two options below:
   
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
ii. From the console (current omics versoin doesn't work well with wdl scoped parameters, cli is preferred):
   
   a. Click **Private Workflows** from the left pane.

   b. Click the **Workflow ID** from the Workflows list.

   c. Click **Create Run** and enter the run information.

## Support Tools

In case your private workflow's run failed, you can [use this script](scripts/healthomics_support/README.md) to extract information and logs from AWS HealthOmics run to ease failures
debugging. Please attach the tar file generated by the script in any support call.
