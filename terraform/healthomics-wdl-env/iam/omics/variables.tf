variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}

variable "project" {
  description = "Project name used in IAM role naming"
  type        = string
}

variable "aws_shared_account_id" {
  description = "AWS account ID with shared ECR images. Optional - if empty, cross-account ECR access is not configured."
  type        = string
  default     = ""
}

variable "aws_region" {
  description = "The AWS region where the module will be deployed."
  type        = string
}

variable "omics_inputs_bucket" {
  description = "Name of the S3 bucket for workflow inputs (granted read access)"
  type        = string
}

variable "omics_outputs_bucket" {
  description = "Name of the S3 bucket for workflow outputs (granted read/write access)"
  type        = string
}

variable "omics_cache_bucket" {
  description = "Name of the S3 bucket for WDL cache (granted read/write access)"
  type        = string
}

variable "bioinfo_resources_bucket" {
  description = "Name of the S3 bucket for bioinformatics resources (granted read access)"
  type        = string
}

variable "broad_references_public_bucket" {
  description = "Name of the Broad Institute public references bucket (granted read access)"
  type        = string
  default     = "broad-references"
}

variable "accessible_buckets_list" {
  description = "List of additional S3 bucket names to grant read access to the OmicsRole"
  type        = list(string)
}

