output "omics_inputs_bucket" {
  value = module.omics_s3_buckets.omics_inputs_bucket
}

output "omics_outputs_bucket" {
  value = module.omics_s3_buckets.omics_outputs_bucket
}

output "omics_cache_bucket" {
  value = module.omics_s3_buckets.omics_cache_bucket
}

output "omics_role_arn" {
  value = module.omics_iam.omics_role_arn
}

# output "omics_standard_run_group" {
#   value = module.omics.standard_run_group_arn
# }
#
# output "omics_long_run_group" {
#   value = module.omics.long_run_group_arn
# }

output "dynamodb_table" {
  value = module.dynamodb.table.name
}