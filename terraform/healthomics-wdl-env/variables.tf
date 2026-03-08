variable "aws_shared_account_id" {
  description = "Account ID of AWS account that stores shared ECR images. Optional - if empty, cross-account ECR access is not configured."
  type        = string
  default     = ""
}
variable "aws_region" {
  description = "The AWS region where the modules will be deployed"
  type        = string
}
variable "project" {
  description = "Project Name. Will be used as prefix for s3 buckets"
  type        = string
}
variable "custom_tags" {
  description = "Map of key->value tags to be added to AWS resources"
  type        = map(string)
}
variable "omics_accessible_buckets" {
  description = "List of s3 buckets to have read access from omics"
  type        = list(string)
}
variable "cross_aws_account" {
  description = "AWS account ID that should have read/write access to the input bucket. Optional - if empty, cross-account S3 access is not configured."
  type        = string
  default     = ""
}
