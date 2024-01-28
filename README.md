# healthomics-workflows
UltimaGenomics repository for workflows compatible with AWS HealthOmics

## Introduction
1. UltimaGenomics offers pipelines as Ready2Run workflows on AWS HealthOmics. Ready2Run workflows enable you to run these pipelines on AWS HealthOmics by simply bringing your data.

2. Utimagenomics also shares community workflows through this repository to allow more flexibility and possibility to modify the workflow. The user can convert the Ready2Run workflow and community workflows to a private workflow on AWS HealthOmics infrastructure.

3. The cost to run the workflow will now be based on the compute and run storage used during the private workflow run instead of a fixed price.

4. For more questions about this workflow, please contact healthomics.support@ultimagen.com.

## Deploying Private Workflow
Once the workflow resources have been deployed into locally (see instructions per workflow), user can create private workflow on AWS HealthOmics
1. Create a private workflow in HealthOmics by following one of the two options below:

i. From the CLI:
 ~~~
$ aws omics create-workflow \
--name <workflow_name> \
--main r2r_efficient_dv.wdl \
--definition-zip <s3_uri_for_zipped_file> \
--parameter-template <s3_uri_for_parameter_template_file> \
--accelerators GPU
 ~~~
ii. From the console:
    
    a. Click on **Private Workflows** from the left pane.
    
    b. Click on **Create Workflow** on the Workflows list.
    
    c. Follow the instructions on the console to create your workflow.

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
