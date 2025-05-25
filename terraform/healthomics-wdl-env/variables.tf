variable "aws_shared_account_id" {
  description = "account ID of AWS account that will store the S3 reference data and ECR images"
  type        = string
}
variable "aws_region" {
  description = "The AWS region where the modules will be deployed."
  type        = string
}
variable "project" {
  description = "Project Name"
  type        = string
}
variable "custom_tags" {
  description = "Map of key->value tags to be added to AWS resources"
  type = map(string)
}
variable "omics_accessible_buckets" {
  description = "List of s3 bucket to have read access from omics"
  type = list(string)
}
variable "cross_aws_account" {
  description = "AWS account id that should have read access to omics input bucket"
  type        = string
}