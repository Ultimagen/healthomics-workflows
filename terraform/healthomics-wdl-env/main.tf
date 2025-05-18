provider "aws" {
  region = var.aws_region
}

data "aws_caller_identity" "current" {}

module "omics_s3_buckets" {
  source = "./s3"

  aws_account_id = data.aws_caller_identity.current.account_id
  aws_region     = var.aws_region
  project        = var.project
  custom_tags    = var.custom_tags
}


module "omics_iam" {
  source = "./iam/omics"

  aws_account_id           = data.aws_caller_identity.current.account_id
  project                  = var.project
  aws_shared_account_id    = var.aws_shared_account_id
  aws_region               = var.aws_region
  omics_inputs_bucket      = module.omics_s3_buckets.omics_inputs_bucket
  omics_outputs_bucket     = module.omics_s3_buckets.omics_outputs_bucket
  omics_cache_bucket       = module.omics_s3_buckets.omics_cache_bucket
  bioinfo_resources_bucket = module.omics_s3_buckets.bioinfo_resources_bucket
}

module "omics" {
  source      = "./health-omics"
  custom_tags = var.custom_tags
}

# module "omics_deployment_iam" {
#   source = "./iam/omics-deployment"
#
#   aws_account_id      = data.aws_caller_identity.current.account_id
#   aws_region          = var.aws_region
#   dynamodb_table_name = module.dynamodb.table.name
# }
