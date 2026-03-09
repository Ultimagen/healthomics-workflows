variable "custom_tags" {
  description = "Map of key-value tags to apply to HealthOmics run groups"
  type        = map(string)
  default     = {}
}
variable "standard_run_group_name" {
  description = "Name of run group for running standard omics runs"
  type        = string
  default     = "standard"
}

variable "long_run_group_name" {
  description = "Name of run group for running long omics runs"
  type        = string
  default     = "long"
}

variable "standard_run_group_max_duration" {
  description = "Standard run group max duration in minutes"
  type        = number
  default     = 7200 # 5 days
}
