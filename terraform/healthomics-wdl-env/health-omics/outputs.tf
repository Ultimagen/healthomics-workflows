output "standard_run_group_id" {
  description = "ID of the standard run group (5-day max duration)"
  value       = awscc_omics_run_group.omics_standard_run_group.id
}

output "long_run_group_id" {
  description = "ID of the long run group (no duration limit)"
  value       = awscc_omics_run_group.omics_long_run_group.id
}

output "standard_run_group_arn" {
  description = "ARN of the standard run group (5-day max duration)"
  value       = awscc_omics_run_group.omics_standard_run_group.arn
}

output "long_run_group_arn" {
  description = "ARN of the long run group (no duration limit)"
  value       = awscc_omics_run_group.omics_long_run_group.arn
}
