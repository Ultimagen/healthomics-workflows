provider "aws" {
  region = var.aws_region
}
locals {
  resource_name_prefix = "${var.aws_region}-start-omics-run-lambda"
  function_name        = "startOmicsRun"
}

data "archive_file" "lambda_zip" {
  type             = "zip"
  source_dir       = "${path.module}/lambda"
  output_path      = "${path.module}/dist/lambda.zip"
  output_file_mode = "0644"
}
resource "aws_lambda_function" "lambda_function" {
  filename         = data.archive_file.lambda_zip.output_path
  source_code_hash = data.archive_file.lambda_zip.output_base64sha256
  function_name    = local.function_name
  role             = aws_iam_role.lambda_role.arn
  handler          = "lambda_function.lambda_handler"
  runtime          = "python3.12"
  timeout          = 60
  depends_on = [aws_cloudwatch_log_group.lambda_log_group]
  environment {
    variables = {
      OMICS_ROLE_ARN           = var.omics_role_arn
      OMICS_OUTPUTS_BUCKET     = var.omics_outputs_bucket
      OMICS_CACHE_BUCKET       = var.omics_cache_bucket
      # STANDARD_RUN_GROUP_ID    = var.omics_standard_run_group_id
      # LONG_RUN_GROUP_ID        = var.omics_long_run_group_id
      OMICS_WORKFLOW_DDB_TABLE = var.dynamodb_table

    }
  }

  tags = merge(
    var.custom_tags,
    {
      Project = var.project
      Name    = "start-omics-run"
    }
  )
}

resource "aws_cloudwatch_log_group" "lambda_log_group" {
  name              = "/aws/lambda/${local.function_name}"
  retention_in_days = 90
}

resource "aws_iam_role_policy_attachment" "lambda_basic_execution_policy_attachment" {
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole"
  role       = aws_iam_role.lambda_role.name
}