variable "custom_tags" {
  type = map(string)
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
