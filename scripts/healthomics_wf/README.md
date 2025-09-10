# Table of Contents
1. [Prerequisites](#Prerequisites)
2. [Create HealthOmics Workflow](#Create-HealthOmics-Workflow)
3. [Create HealthOmics Workflow](#Invoke-HealthOmics-Run)


## Prerequisites
- `Python 3.8+`
-  Use your favorite python env to install the needed packages listed in: `requirements.txt`

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
- `--s3-bucket`: Bucket name to copy resources files to. This bucket must be accessed by the service role that will be used to run the workflow.
##### Optional Arguments
- `--s3-prefix`: S3 key prefix to place files under in the destination bucket (optional). If empty, the resources files keys will be similar to source keys.
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
#### S3 Permissions
`s3:ListBucket`
`s3:GetObject`
`s3:PutObject`
#### HealthOmics Permissions
`omics:CreateWorkflow`
`omics:GetWorkflow`
`omics:ListWorkflows`
#### Dynamodb Permissions
`dynamodb:PutItem`

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
        "s3:ListBucket",
        "s3:GetObject",
        "s3:PutObject",
        "omics:CreateWorkflow",
        "omics:GetWorkflow",
        "omics:ListWorkflows"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "dynamodb:PutItem"
      ],
      "Resource": "arn:aws:dynamodb:AWS_REGION:AWS_ACCOUNT:table/OmicsWorkflows"
    }
  ]
}
```

## Invoke HealthOmics Run

### Overview
The purpose of `invoke_healthomics_run.py` script is
to start a specific version of AWS HealthOmics private workflow
that deployed using the `create_healthomics_workflow.py` script,
and it's versioning is managed using dynamodb.

### Usage


#### Command Line Arguments

The script accepts the following command line arguments:
   parser.add_argument("--omics-workflow-name", required=True, help="Workflow name")
    parser.add_argument("--workflow-version", required=True, help="Workflow version")
    parser.add_argument("--run-id", required=True, help="Run ID")
    parser.add_argument("--input-params-file", required=True, help="Path to the JSON file for InputParams")
    parser.add_argument("--long-run", action='store_true', help="mark the run as long and run it under 'long' run group. By default, will run under 'standard' run group")
    parser.add_argument("--cache-behavior", help="HealthOmics cache behavior")
    parser.add_argument("--aws-region", help="AWS region", default="us-east-1")
    parser.add_argument("--aws-profile", help="AWS CLI profile", required=False)
- `--omics-workflow-name`: Workflow name to run.
- `--workflow-version`: Workflow version (one of healthomics-workflows release [tags](https://github.com/Ultimagen/healthomics-workflows/tags)).
- `--run-id`: Run ID for the run.
- `--input-params-file`: Path to the JSON file for InputParams.
##### Optional Arguments
- `--long-run`: mark the run as long and run it under 'long' run group. By default, will run under 'standard' run group.
- `--aws-region`: AWS region (default=us-east-1).
- `--cache-behavior`: HealthOmics cache behavior (CACHE_ALWAYS/CACHE_ON_FAILURE) (By default, will use CACHE_ON_FAILURE)
- `--aws-profile`: AWS CLI profile (optional). If empty, it will use the current session's credentials.

#### Example Command

```bash
python invoke_healthomics_run.py \
--workflow-name germline_CNV_pipeline \
--workflow-version 1.18.3 \
--run-id DATA-8079 \
--input-params-file germline_CNV_pipeline_input.json \
--aws-region eu-west-1
```

### AWS Permissions

To successfully run this script, you need the following AWS permissions:

`lambda:InvokeFunction`
`dynamodb:GetItem`
