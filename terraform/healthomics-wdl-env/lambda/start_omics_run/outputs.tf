output "lambda_function_arn" {
  value       = aws_lambda_function.lambda_function.arn
  description = "The ARN of the deployed Lambda function."
}
output "lambda_function_name" {
  value       = aws_lambda_function.lambda_function.function_name
  description = "The function_name of the deployed Lambda function."
}
output "lambda_logs_group" {
  value       = aws_cloudwatch_log_group.lambda_log_group.name
  description = "The log_group name of the deployed Lambda function."
}