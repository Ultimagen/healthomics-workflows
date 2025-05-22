#!/bin/bash

set -eo pipefail

ACCOUNT=$(aws sts get-caller-identity | jq .Account | cut -d '"' -f 2)

REGION=${AWS_REGION:-$(aws ec2 describe-availability-zones --output text --query 'AvailabilityZones[0].[RegionName]')}
DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
TERRAFORM_BUCKET=$1
PROJECT=${3:-healthomics-workflows}
terraform init -reconfigure -backend-config="bucket=$TERRAFORM_BUCKET" -backend-config="key=$PROJECT/$REGION/terraform.tfstate"
