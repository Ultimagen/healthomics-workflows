variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}

variable "aws_region" {
  description = "The AWS region where the module will be deployed."
  type        = string
}

variable "project" {
  type = string
}
variable "inputs_bucket_suffix" {
  type    = string
  default = "runs-data"
}

variable "outputs_bucket_suffix" {
  type    = string
  default = "outputs"
}

variable "bioinfo_resources_bucket_suffix" {
  type    = string
  default = "bioinfo-resources"
}

variable "cache_bucket_suffix" {
  type    = string
  default = "wdl-cache"
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
  type = map(string)
}

variable "cross_account_id" {
  type    = string
  default = "525048827230" #todo now - pass from main
}