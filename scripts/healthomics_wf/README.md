# Table of Contents

1. [Prerequisites](#prerequisites)
2. [Create HealthOmics Workflow](#create-healthomics-workflow)
3. [Invoke HealthOmics Run](#invoke-healthomics-run)


## Prerequisites
- `Python 3.8+`
- Docker
- Use your favorite python env to install the needed packages listed in: `requirements.txt`

## Create HealthOmics Workflow

### Overview
The purpose of `create_healthomics_workflow.py` script is to generate AWS HealthOmics private workflows in a custom AWS account from a workflow under [this folder](../../workflows).
The script performs the following tasks:

1. Parses a given WDL to extract Docker image references and copies them from their source registries to a specified destination Amazon Elastic Container Registry (ECR).
2. Parses a JSON template file/s to extract S3 file references and copies them from their source locations to a specified destination S3 bucket.
3. Uses workflow files to create an AWS HealthOmics private workflow.

### Usage


#### Command Line Arguments

The script accepts the following command line arguments:

- `workflow`: Workflow to create (folder under workflows).
- `--aws-region`: AWS region for ECR images.
- `--s3-bucket`: Bucket name to copy resources files to. This bucket must be accessible by the OmicsRole. If using the [Terraform module](../../terraform/healthomics-wdl-env), ensure this bucket is included in the `omics_accessible_buckets` variable.

##### Optional Arguments
- `--input-template`: Input template JSON file path to localize (optional). If empty, it will localize all input templates.
- `--use-dynamodb`: Add this flag for storing the version workflow in dynamodb table (Relevant when omics environment was created using the terraform code that shared in this repository)
- `--omics-workflow-name`: Name for the generated omics workflow (optional). If empty, it will use the `workflow` argument.
- `--aws-profile`: AWS CLI profile (optional). If empty, it will use the current session's credentials.

#### Example Command

```bash
python create_healthomics_workflow.py \
  single_read_snv \
  --aws-region us-east-1 \
  --s3-bucket my-s3-bucket \
  --input-template single_read_snv_template-ppmSeq.json \
  --use-dynamodb
```

### AWS Permissions
To successfully run this script, you need the following AWS permissions:

#### ECR Permissions
`ecr:DescribeRepositories`
`ecr:CreateRepository`
`ecr:BatchCheckLayerAvailability`
`ecr:BatchGetImage`
`ecr:GetDownloadUrlForLayer`
`ecr:PutImage`
`ecr:CompleteLayerUpload`
`ecr:InitiateLayerUpload`
`ecr:UploadLayerPart`
`ecr:GetAuthorizationToken`
`ecr:DescribeImages`
#### S3 Permissions
`s3:ListBucket`
`s3:GetObject`
`s3:PutObject`
#### HealthOmics Permissions
`omics:CreateWorkflow`
`omics:CreateWorkflowVersion`
`omics:GetWorkflow`
`omics:GetWorkflowVersion`
`omics:ListWorkflows`

#### DynamoDB Permissions
`dynamodb:PutItem`
`dynamodb:GetItem`

### Example IAM Policy
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
        "ecr:GetAuthorizationToken",
        "ecr:DescribeImages",
        "s3:ListBucket",
        "s3:GetObject",
        "s3:PutObject",
        "omics:CreateWorkflow",
        "omics:CreateWorkflowVersion",
        "omics:GetWorkflow",
        "omics:GetWorkflowVersion",
        "omics:ListWorkflows"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "dynamodb:PutItem",
        "dynamodb:GetItem"
      ],
      "Resource": "arn:aws:dynamodb:AWS_REGION:AWS_ACCOUNT:table/OmicsWorkflows"
    }
  ]
}
```

## Invoke HealthOmics Run

### Overview

The `invoke_healthomics_run.py` script starts an AWS HealthOmics private workflow
that was deployed using `create_healthomics_workflow.py` with the `--use-dynamodb` flag.

### Usage

#### Command Line Arguments

##### Required Arguments

- `--omics-workflow-name`: Workflow name to run.
- `--run-id`: Run ID for the run.
- `--input-params-file`: Path to the JSON file for InputParams.

##### Optional Arguments

- `--workflow-version`: Workflow version. If not provided, the latest version from DynamoDB will be used automatically.
- `--long-run`: Mark the run as long and run it under 'long' run group. By default, will run under 'standard' run group.
- `--aws-region`: AWS region (default: us-east-1).
- `--cache-behavior`: HealthOmics cache behavior (`CACHE_ALWAYS` or `CACHE_ON_FAILURE`). Default: `CACHE_ON_FAILURE`.
- `--aws-profile`: AWS CLI profile. If empty, it will use the current session's credentials.

#### Example Commands

Using latest version (recommended):

```bash
python invoke_healthomics_run.py \
  --omics-workflow-name germline_CNV_pipeline \
  --run-id DATA-8079 \
  --input-params-file germline_CNV_pipeline_input.json \
  --aws-region eu-west-1
```

Using specific version:

```bash
python invoke_healthomics_run.py \
  --omics-workflow-name germline_CNV_pipeline \
  --workflow-version 1.18.3 \
  --run-id DATA-8079 \
  --input-params-file germline_CNV_pipeline_input.json \
  --aws-region eu-west-1
```

### AWS Permissions

To successfully run this script, you need the following AWS permissions:

`lambda:InvokeFunction`
`dynamodb:Scan` (required when `--workflow-version` is not provided, to find the latest version)
`dynamodb:GetItem`
