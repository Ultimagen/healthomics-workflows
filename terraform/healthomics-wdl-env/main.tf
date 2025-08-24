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
  # cross_aws_account_id = var.cross_aws_account
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
  accessible_buckets_list = var.omics_accessible_buckets
}

# module "omics" {
#   source      = "./health-omics"
#   custom_tags = var.custom_tags
# }

module "dynamodb" {
  source = "./dynamodb"
}

module "start_omics_run_lambda" {
  source = "./lambda/start_omics_run"
  aws_region = var.aws_region
  custom_tags = var.custom_tags
  project = var.project
  aws_account_id = data.aws_caller_identity.current.account_id
  omics_cache_bucket = module.omics_s3_buckets.omics_cache_bucket
  # omics_long_run_group_id = module.omics.long_run_group_id
  omics_outputs_bucket = module.omics_s3_buckets.omics_outputs_bucket
  omics_role_arn = module.omics_iam.omics_role_arn
  # omics_standard_run_group_id = module.omics.standard_run_group_id
  dynamodb_table = module.dynamodb.table.name
}

