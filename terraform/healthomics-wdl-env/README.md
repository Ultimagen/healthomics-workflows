# AWS HealthOmics WDL Environment

Terraform module to deploy AWS infrastructure for running Ultima Genomics WDL workflows on AWS HealthOmics. This module automates the creation of S3 buckets, IAM roles, HealthOmics run groups, DynamoDB tables, and Lambda functions required for workflow execution.

## Architecture

This module creates the following resources:

### S3 Buckets

Bucket names follow the pattern: `{project}-{account_id}-{region}-{suffix}`

| Bucket | Suffix | Purpose |
|--------|--------|---------|
| Inputs | `runs-data` | Workflow input files |
| Outputs | `outputs` | Workflow output files |
| Cache | `wdl-cache` | WDL execution cache (CACHE_ALWAYS, CACHE_ON_FAILURE) |
| Bioinfo Resources | `bioinfo-resources` | Bioinformatics reference data |

All buckets include lifecycle policies:
- Transition to INTELLIGENT_TIERING after 7 days
- Archive to GLACIER after 60 days (inputs/outputs, configurable)
- Cache expiration policies for CACHE_ALWAYS (60 days) and CACHE_ON_FAILURE (7 days)

### HealthOmics Run Groups
| Run Group | Max Duration | Purpose |
|-----------|--------------|---------|
| Standard | 5 days (7200 min) | Regular workflow runs |
| Long | No limit | Extended workflow runs |

### IAM Resources
- **OmicsRole**: IAM role assumed by HealthOmics and EC2 (for miniwdl local execution)
- **Instance Profile**: For local miniwdl execution on EC2

### DynamoDB Table
- **OmicsWorkflows**: Stores workflow metadata and version-to-workflow-ID mappings
  - Partition key: `version` (String)
  - Sort key: `workflow` (String)

### Lambda Function
- **startOmicsRun**: Triggers HealthOmics workflow runs
  - Queries DynamoDB for workflow IDs
  - Creates/retrieves run caches
  - Starts HealthOmics runs with appropriate configurations

## Prerequisites

- AWS CLI configured and authenticated
- Terraform >= 1.0
- AWS region with HealthOmics support
- **S3 bucket for Terraform state must be created manually before running this module** (used by `init.sh` to store remote state)

## Required IAM Permissions (Deployer)

The user or role running Terraform must have permissions to create the following resources:

| Service | Required Permissions |
|---------|---------------------|
| S3 | `s3:CreateBucket`, `s3:PutBucketPolicy`, `s3:PutBucketPublicAccessBlock`, `s3:PutLifecycleConfiguration`, `s3:GetBucketPolicy`, `s3:GetLifecycleConfiguration` |
| IAM | `iam:CreateRole`, `iam:CreatePolicy`, `iam:AttachRolePolicy`, `iam:CreateInstanceProfile`, `iam:AddRoleToInstanceProfile`, `iam:GetRole`, `iam:GetPolicy`, `iam:PassRole` |
| Lambda | `lambda:CreateFunction`, `lambda:GetFunction`, `lambda:UpdateFunctionCode`, `lambda:UpdateFunctionConfiguration` |
| DynamoDB | `dynamodb:CreateTable`, `dynamodb:DescribeTable`, `dynamodb:UpdateTable` |
| HealthOmics | `omics:CreateRunGroup`, `omics:GetRunGroup`, `omics:TagResource` |
| CloudWatch Logs | `logs:CreateLogGroup`, `logs:PutRetentionPolicy`, `logs:DescribeLogGroups` |
| ECR | `ecr:GetAuthorizationToken` (for validation) |

## Created Resources and Permissions

### OmicsRole Permissions

The OmicsRole created by this module grants HealthOmics workflows the following permissions:

| Permission | Resources |
|------------|-----------|
| S3 Read | Input, output, cache, bioinfo-resources buckets, configured accessible buckets, broad-references |
| S3 Write | Output and cache buckets |
| CloudWatch Logs | Create log groups/streams in `/aws/omics/WorkflowLog` |
| ECR | Pull images from current and shared account repositories |
| Sentieon | Access to Sentieon license bucket |

### Lambda Execution Role Permissions

| Permission | Resources |
|------------|-----------|
| HealthOmics | Full access (`omics:*`) |
| IAM | PassRole to OmicsRole |
| S3 | Full access to cache bucket |
| DynamoDB | GetItem from OmicsWorkflows table |
| CloudWatch Logs | Create and write logs |

### Cross-Account Access (Optional)

If `cross_aws_account` is provided, the input bucket is configured to allow read/write access from the specified external AWS account.

## Usage

### 1. Initialize Terraform Backend

The `init.sh` script configures the Terraform backend with S3 state storage.

```bash
./init.sh <terraform_bucket_name> [project_name]
```

| Parameter               | Required | Default                 | Description                                  |
|-------------------------|----------|-------------------------|----------------------------------------------|
| `terraform_bucket_name` | Yes      | -                       | S3 bucket name for storing Terraform state   |
| `project_name`          | No       | `healthomics-workflows` | Project name used in state key path          |

The script:
1. Retrieves the current AWS account ID
2. Uses `AWS_REGION` environment variable (or detects from AWS config)
3. Initializes Terraform with backend: `s3://{bucket}/{project}/{region}/terraform.tfstate`

### 2. Configure Variables

Create a `terraform.tfvars` file:

```hcl
aws_region = "us-east-1"
project    = "healthomics-workflows"

custom_tags = {
  Environment = "production"
  Team        = "bioinformatics"
}

omics_accessible_buckets = [
  "my-reference-data-bucket",
  "another-data-bucket"
]

# Optional: Account with shared ECR images and reference data
# aws_shared_account_id = "123456789012"

# Optional: Account with read/write access to input bucket
# cross_aws_account = "987654321098"
```

### 3. Deploy

```bash
export AWS_REGION=us-east-1
./init.sh my-terraform-state-bucket
terraform plan
terraform apply
```

### 4. (Optional) Request HealthOmics Quotas

```bash
./scripts/request_omics_quotas.sh --help
```

## Variables

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `aws_region` | string | Yes | - | AWS region for deployment |
| `project` | string | Yes | - | Project name used as prefix for S3 bucket names (bucket naming: `{project}-{account_id}-{region}-{suffix}`) |
| `custom_tags` | map(string) | Yes | - | Key-value tags to apply to all created resources |
| `omics_accessible_buckets` | list(string) | Yes | - | List of additional S3 bucket names for OmicsRole read access |
| `aws_shared_account_id` | string | No | `""` | AWS account ID storing shared ECR images (enables cross-account ECR access if provided) |
| `cross_aws_account` | string | No | `""` | AWS account ID granted read/write access to the input bucket (enables cross-account S3 access if provided) |

## Outputs

| Name | Description |
|------|-------------|
| `omics_inputs_bucket` | S3 bucket name for workflow input files |
| `omics_outputs_bucket` | S3 bucket name for workflow output files |
| `omics_cache_bucket` | S3 bucket name for WDL execution cache |
| `omics_role_arn` | ARN of the IAM role used by HealthOmics workflow runs |
| `omics_standard_run_group` | ARN of the run group with 5-day maximum duration |
| `omics_long_run_group` | ARN of the run group with no duration limit |
| `dynamodb_table` | Name of the DynamoDB table storing workflow metadata |

## Important: S3 Bucket Access for Workflow Resources

When deploying workflows using `create_healthomics_workflow.py`, you specify an S3 bucket via the `--s3-bucket` parameter where workflow resources (reference files, etc.) are copied.

For HealthOmics runs to access these resources, **this bucket must be included in the `omics_accessible_buckets` variable** when deploying this Terraform module.

**Terraform deployment** - Include your resources bucket:

```hcl
omics_accessible_buckets = [
  "my-workflow-resources-bucket"
]
```

**Workflow deployment** - Use the same bucket:

```bash
python create_healthomics_workflow.py my-workflow \
  --s3-bucket my-workflow-resources-bucket \
  --aws-region us-east-1 \
  --use-dynamodb
```

## Next Steps

After deploying the infrastructure, refer to the [main README](../../README.md) for:

- **Deploying workflows** - See [Deploying Private Workflow](../../README.md#deploying-private-workflow)
- **Running workflows** - See [Running Private Workflow](../../README.md#running-private-workflow)
