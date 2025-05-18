output "omics_deploy_role_arn" {
  value = aws_iam_role.omics_deployment_role.arn
}

output "omics_test_ec2_runner_role_arn" {
  value = aws_iam_role.omics_test_ec2_runner_role.arn
}
