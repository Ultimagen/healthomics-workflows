# healthomics-workflows
UltimaGenomics repository for workflows compatible with AWS HealthOmics

## Introduction

1.	Ultima Genomics offers pipelines as Ready2Run workflows on AWS HealthOmics. Ready2Run workflows enable you to run these pipelines on AWS HealthOmics by simply bringing your data. For more flexibility such as the use of larger file sizes or changing the reference genome, you can convert Ready2Run workflows to private workflows by following the steps in this repository. Once the Ready2Run workflow is converted to a private workflow, the cost to run the workflow will now be based on the compute and run storage used during the private workflow.

2.	Ultima Genomics also shares pipelines that has been modified to run as private workflows on AWS HealthOmics in this repository. You can follow the directions in this repository to create and run a private workflow on AWS HealthOmics.

3.	For more questions about these workflows, please contact healtomics.support@ultimagen.com.

## Deploying Private Workflow
Once the workflow resources have been deployed into locally (see instructions per workflow), user can create private workflow on AWS HealthOmics
1. Create a private workflow in HealthOmics by following one of the two options below:

i. From the CLI:
 ~~~
$ aws omics create-workflow \
--name <workflow_name> \
--main <main wdl file listed in the folder> \
--definition-zip <s3_uri_for_zipped_file> \
--parameter-template <s3_uri_for_parameter_template_file> \
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
 --parameters <s3_uri_or_local_parameters_file>
 --name <run_name> \
 --retention-mode REMOVE
 ~~~
ii. From the console:
   
   a. Click **Private Workflows** from the left pane.

   b. Click the **Workflow ID** from the Workflows list.

   c. Click **Create Run** and enter the run information.
