output "table" {
  description = "The DynamoDB table object for storing workflow metadata"
  value       = aws_dynamodb_table.omics_table
}

