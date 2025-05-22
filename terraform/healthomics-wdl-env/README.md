Terraform modules to simplify the deployment of AWS resources required to run Ultima Genomics WDL workflows on AWS HealthOmics. These modules automate the creation of necessary infrastructure, such as S3 buckets, IAM roles, and permissions, ensuring a consistent and secure environment for workflow execution.

## Apply Terraform on a region on an AWS account
- Make sure that the AWS healthomics is supported in the required region.
- Under terraform/healthomics-wdl-env run:
  - create bucket for terraform state
  - Run ./scripts/request_omics_quotas.sh --help and follow the instructions to create quota request as per your needs.
  - run export AWS_REGION=<aws-region-you-want-the-modules-to-be-deployed>
  - ./init.sh <terraform-bucket-name> 
  - terraform plan/apply
  