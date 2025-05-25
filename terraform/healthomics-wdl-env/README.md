Terraform modules to simplify the deployment of AWS resources required to run Ultima Genomics WDL workflows on AWS HealthOmics. These modules automate the creation of necessary infrastructure, such as S3 buckets, IAM roles, and permissions, ensuring a consistent and secure environment for workflow execution.

## Apply Terraform on a region on an AWS account
- Make sure that the AWS healthomics is supported in the required region.
- Under terraform/healthomics-wdl-env run:
  1. Create bucket for terraform state
  2. Run:
     1. `export AWS_REGION=<aws region you want the modules to be deployed on>`
     2. `./init.sh <terraform bucket name created on step 1>` 
     3. `terraform plan/apply`
  3. [optional] Run `./scripts/request_omics_quotas.sh --help` and follow the instructions to create quota request as per your needs.
  