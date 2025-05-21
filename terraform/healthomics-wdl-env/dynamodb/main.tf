resource "aws_dynamodb_table" "omics_table" {
  name                        = "OmicsWorkflows"
  billing_mode                = "PAY_PER_REQUEST"
  deletion_protection_enabled = false

  hash_key  = "version"
  range_key = "workflow"

  attribute {
    name = "version"
    type = "S"
  }

  attribute {
    name = "workflow"
    type = "S"
  }
}
