variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}

variable "aws_region" {
  description = "The AWS region where the module will be deployed."
  type        = string
}

variable "project" {
  description = "Project name used as prefix for S3 bucket names"
  type        = string
}
variable "inputs_bucket_suffix" {
  description = "Suffix for the workflow inputs S3 bucket name"
  type        = string
  default     = "runs-data"
}

variable "outputs_bucket_suffix" {
  description = "Suffix for the workflow outputs S3 bucket name"
  type        = string
  default     = "outputs"
}

variable "bioinfo_resources_bucket_suffix" {
  description = "Suffix for the bioinformatics resources S3 bucket name"
  type        = string
  default     = "bioinfo-resources"
}

variable "cache_bucket_suffix" {
  description = "Suffix for the WDL cache S3 bucket name"
  type        = string
  default     = "wdl-cache"
}

variable "input_bucket_archive_days" {
  description = "Archive objects in input bucket after X days"
  type        = number
  default     = 60
}

variable "output_bucket_archive_days" {
  description = "Archive objects in output bucket after X days"
  type        = number
  default     = 60
}

variable "cache_always_expiration_days" {
  description = "Expire CACHE_ALWAYS files in cache bucket after X days"
  type        = number
  default     = 60
}

variable "cache_always_prefix" {
  description = "CACHE_ALWAYS prefix"
  type        = string
  default     = "cache_always/"
}

variable "cache_on_failure_expiration_days" {
  description = "Expire CACHE_ON_FAILURE files in cache bucket after X days"
  type        = number
  default     = 7
}

variable "cache_on_failure_prefix" {
  description = "CACHE_ON_FAILURE prefix"
  type        = string
  default     = "cache_on_failure/"
}

variable "custom_tags" {
  description = "Map of key-value tags to apply to all S3 buckets"
  type        = map(string)
}

variable "cross_aws_account_id" {
  description = "AWS account ID granted read/write access to the input bucket. Optional - if empty, cross-account access is not configured."
  type        = string
  default     = ""
}