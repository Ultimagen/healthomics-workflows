# Create HealthOmics Workflow

## Overview
The purpose of `create_healthomics_workflow.py` script is to generate AWS HealthOmics private workflows in a custom AWS account from a workflow under [this folder](../../workflows).
The script performs the following tasks:

1. Parses a given WDL to extract Docker image references and copies them from their source registries to a specified destination Amazon Elastic Container Registry (ECR).
2. Parses a JSON template file/s to extract S3 file references and copies them from their source locations to a specified destination S3 bucket.
3. Uses workflow files to create an AWS HealthOmics private workflow.

## Usage

### Prerequisites
- `Python 3.8+`
-  Use your favorite python env to install the needed packages listed in: `requirements.txt`

### Command Line Arguments

The script accepts the following command line arguments:

- `workflow`: Workflow to create (folder under workflows).
- `aws_region`: AWS region for ECR images.
- `s3_bucket`: Bucket name to copy resources files to. This bucket must be accessed by the service role that will be used to run the workflow.
#### Optional Arguments
- `--input_template`: Input template JSON file name to localize (optional). If empty, it will localize all input templates.
- `--omics_workflow_name`: Name for the generated omics workflow (optional). If empty, it will use the `workflow` argument.
- `--aws_profile`: AWS CLI profile (optional). If empty, it will use the current session's credentials.

### Example Command

```bash
python create_healthomics_workflow.py \
  single_rad_snv \
  us-east-1 \
  my-s3-bucket \
  --input_template single_read_snv_template-ppmSeq.json \
  --omics_workflow_name single_rad_snv_wf
```

## AWS Permissions

To successfully run this script, you need the following AWS permissions:

### ECR Permissions
`ecr:DescribeRepositories`
`ecr:CreateRepository`
`ecr:BatchCheckLayerAvailability`
`ecr:BatchGetImage`
`ecr:GetDownloadUrlForLayer`
`ecr:PutImage`
`ecr:CompleteLayerUpload`
`ecr:InitiateLayerUpload`
`ecr:UploadLayerPart`
### S3 Permissions
`s3:ListBucket`
`s3:GetObject`
`s3:PutObject`
### HealthOmics Permissions
`omics:CreateWorkflow`
`omics:GetWorkflow`
`omics:ListWorkflows`


## Example IAM Policy
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "ecr:DescribeRepositories",
        "ecr:CreateRepository",
        "ecr:BatchCheckLayerAvailability",
        "ecr:BatchGetImage",
        "ecr:GetDownloadUrlForLayer",
        "ecr:PutImage",
        "ecr:CompleteLayerUpload",
        "ecr:InitiateLayerUpload",
        "ecr:UploadLayerPart",
        "s3:ListBucket",
        "s3:GetObject",
        "s3:PutObject",
        "omics:CreateWorkflow",
        "omics:GetWorkflow",
        "omics:ListWorkflows"
      ],
      "Resource": "*"
    }
  ]
}
```
