variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}
variable "aws_region" {
  description = "The AWS region where the resource will be deployed."
  type        = string
}
variable "omics_role_arn" {
  description = "ARN of the IAM role to be assumed by HealthOmics runs"
  type        = string
}
variable "omics_outputs_bucket" {
  description = "Name of the S3 bucket for workflow outputs"
  type        = string
}

variable "omics_cache_bucket" {
  description = "Name of the S3 bucket for WDL execution cache"
  type        = string
}
variable "project" {
  description = "Project name used in Lambda function naming"
  type        = string
}

variable "custom_tags" {
  description = "Map of key-value tags to apply to Lambda resources"
  type        = map(string)
  default     = {}
}

variable "omics_standard_run_group_id" {
  description = "ID of the standard HealthOmics run group (5-day max duration)"
  type        = string
}

variable "omics_long_run_group_id" {
  description = "ID of the long HealthOmics run group (no duration limit)"
  type        = string
}

variable "dynamodb_table" {
  description = "Name of the DynamoDB table containing workflow metadata"
  type        = string
}