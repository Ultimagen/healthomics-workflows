variable "aws_account_id" {
  description = "aws account ID"
  type        = string
}

variable "aws_region" {
  description = "aws region"
  type        = string
}

variable "dynamodb_table_name" {
  description = "dynamo DB table name"
}

variable "github_org" {
  description = "Github organization"
  type        = string
  default     = "Ultimagen"
}

variable "github_actions_instance_profile_name" {
  default = "github-actions-omics-role"
}
