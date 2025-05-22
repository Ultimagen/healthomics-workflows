resource "aws_iam_role" "lambda_role" {
  name = "${local.resource_name_prefix}_role"

  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "lambda.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF
}

resource "aws_iam_policy" "lambda_policy" {
  name   = "${local.resource_name_prefix}_policy"
  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": "logs:CreateLogGroup",
      "Resource": "arn:aws:logs:${var.aws_region}:${var.aws_account_id}:*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogStream",
        "logs:PutLogEvents"
      ],
      "Resource": "arn:aws:logs:${var.aws_region}:${var.aws_account_id}:log-group:/aws/lambda/${aws_lambda_function.lambda_function.function_name}:*"
    },
    {
      "Effect": "Allow",
      "Action": [
         "omics:*",
         "tag:TagResources"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
         "iam:PassRole"
      ],
      "Resource": "${var.omics_role_arn}"
    },
    {
      "Effect": "Allow",
      "Action": [
         "s3:*"
      ],
       "Resource": [
            "arn:aws:s3:::${var.omics_cache_bucket}/*",
            "arn:aws:s3:::${var.omics_cache_bucket}"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
         "dynamodb:GetItem"
      ],
       "Resource": [
            "arn:aws:dynamodb:${var.aws_region}:${var.aws_account_id}:table/${var.dynamodb_table}"
      ]
    }
  ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "lambda_policy_attachment" {
  policy_arn = aws_iam_policy.lambda_policy.arn
  role       = aws_iam_role.lambda_role.name
}