variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}

variable "project" {
  type = string
}

variable "aws_region" {
  description = "The AWS region where the module will be deployed."
  type        = string
}

variable "omics_inputs_bucket" {
  type = string
}

variable "omics_outputs_bucket" {
  type = string
}

variable "omics_cache_bucket" {
  type = string
}

variable "bioinfo_resources_bucket" {
  type = string
}

variable "broad_references_public_bucket" {
  type    = string
  default = "broad-references"
}

variable "accessible_buckets_list" {
  type = list(string)
}

