variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}
variable "aws_region" {
  description = "The AWS region where the resource will be deployed."
  type        = string
}
variable "omics_role_arn" {
  type = string
}
variable "omics_outputs_bucket" {
  type = string
}

variable "omics_cache_bucket" {
  type = string
}
variable "project" {}

variable "custom_tags" {
  type = map(string)
}

variable "omics_standard_run_group_id" {
  type = string
}

variable "omics_long_run_group_id" {
  type = string
}

variable "dynamodb_table" {
  type = string
}