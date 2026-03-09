output "omics_inputs_bucket" {
  description = "S3 bucket name for workflow input files"
  value       = module.omics_s3_buckets.omics_inputs_bucket
}

output "omics_outputs_bucket" {
  description = "S3 bucket name for workflow output files"
  value       = module.omics_s3_buckets.omics_outputs_bucket
}

output "omics_cache_bucket" {
  description = "S3 bucket name for WDL execution cache"
  value       = module.omics_s3_buckets.omics_cache_bucket
}

output "omics_role_arn" {
  description = "ARN of the IAM role used by HealthOmics workflow runs"
  value       = module.omics_iam.omics_role_arn
}

output "omics_standard_run_group" {
  description = "ARN of the run group with 5-day maximum duration"
  value       = module.omics.standard_run_group_arn
}

output "omics_long_run_group" {
  description = "ARN of the run group with no duration limit"
  value       = module.omics.long_run_group_arn
}

output "dynamodb_table" {
  description = "Name of the DynamoDB table storing workflow metadata"
  value       = module.dynamodb.table.name
}