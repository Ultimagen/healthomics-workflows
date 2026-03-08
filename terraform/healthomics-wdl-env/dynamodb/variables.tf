variable "table_name" {
  description = "Name of the DynamoDB table for storing workflow metadata (partition key: version, sort key: workflow)"
  type        = string
  default     = "OmicsWorkflows"
}
