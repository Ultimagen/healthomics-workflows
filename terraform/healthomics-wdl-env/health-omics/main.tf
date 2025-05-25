resource "awscc_omics_run_group" "omics_standard_run_group" {
  name         = var.standard_run_group_name
  max_duration = var.standard_run_group_max_duration
  tags         = var.custom_tags
}

resource "awscc_omics_run_group" "omics_long_run_group" {
  name = var.long_run_group_name
  # This run group has no max_duration limitation - runs will be limited by the service Workflows - Maximum run duration quota

  tags = var.custom_tags
}