output "omics_role_arn" {
  description = "ARN of the IAM role assumed by HealthOmics for workflow execution"
  value       = aws_iam_role.omics_role.arn
}

