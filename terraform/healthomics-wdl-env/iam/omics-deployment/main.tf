resource "aws_iam_role" "omics_deployment_role" {
  name                 = "UltimagenOmicsDeploymentTest-role-${var.aws_region}"
  assume_role_policy = data.aws_iam_policy_document.assume_role_policy.json # (not shown)
  max_session_duration = 43200  # 12 h
}

resource "aws_iam_role" "omics_test_ec2_runner_role" {
  name                 = "GithubActionsEc2Runner-role-${var.aws_region}"
  assume_role_policy = data.aws_iam_policy_document.assume_role_policy.json # (not shown)
  max_session_duration = 43200 # 12 h
}

resource "aws_iam_policy" "iam_policy" {
  name = "UltimagenOmicsDeployment-iam-policy-${var.aws_region}"

  policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Effect" : "Allow",
        "Action" : [
          "s3:*",
          "omics:*",
          "lambda:InvokeFunction",
          "iam:PassRole",
          "tag:TagResources"
        ],
        "Resource" : "*"
      },
      {
        "Sid" : "DynamoDBIndexAndStreamAccess",
        "Effect" : "Allow",
        "Action" : [
          "dynamodb:GetShardIterator",
          "dynamodb:Scan",
          "dynamodb:Query",
          "dynamodb:DescribeStream",
          "dynamodb:GetRecords",
          "dynamodb:ListStreams",
          "dynamodb:DescribeLimits"
        ],
        "Resource" : [
          "arn:aws:dynamodb:*:${var.aws_account_id}:table/${var.dynamodb_table_name}/index/*",
          "arn:aws:dynamodb:*:${var.aws_account_id}:table/${var.dynamodb_table_name}/stream/*"
        ]
      },
      {
        "Sid" : "DynamoDBTableAccess"
        "Effect" : "Allow",
        "Action" : [
          "dynamodb:BatchGetItem",
          "dynamodb:BatchWriteItem",
          "dynamodb:ConditionCheckItem",
          "dynamodb:PutItem",
          "dynamodb:DescribeTable",
          "dynamodb:DeleteItem",
          "dynamodb:GetItem",
          "dynamodb:Scan",
          "dynamodb:Query",
          "dynamodb:UpdateItem"
        ],
        "Resource" : [
          "arn:aws:dynamodb:*:${var.aws_account_id}:table/${var.dynamodb_table_name}"
        ]
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "ecr_role_policy_attachment" {
  role       = aws_iam_role.omics_deployment_role.name
  policy_arn = aws_iam_policy.iam_policy.arn
}

resource "aws_iam_role_policy_attachment" "ssm_role_policy_attachment" {
  # when used as an instance profile ssm allows to connect to debug the runner
  role       = aws_iam_role.omics_deployment_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_policy" "ec2_runner_iam_policy" {
  name = "GithubActionsEc2Runner-iam-policy-${var.aws_region}"

  policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Effect" : "Allow",
        "Action" : [
          "ec2:RunInstances",
          "ec2:TerminateInstances",
          "ec2:Describe*",
          "ec2:AuthorizeSecurityGroupEgress",
          "ec2:AuthorizeSecurityGroupIngress",
          "ec2:RevokeSecurityGroupEgress",
          "ec2:RevokeSecurityGroupIngress",
          "ec2:CreateSecurityGroup",
          "ec2:DeleteSecurityGroup"
        ],
        "Resource" : "*"
      },
      {
        "Effect" : "Allow",
        "Action" : [
          "ec2:ReplaceIamInstanceProfileAssociation",
          "ec2:AssociateIamInstanceProfile"
        ],
        "Resource" : "*"
      },
      {
        "Effect" : "Allow",
        "Action" : "iam:PassRole",
        "Resource" : "*"
      },
      {
        "Effect" : "Allow",
        "Action" : [
          "ec2:CreateTags"
        ],
        "Resource" : "*",
        "Condition" : {
          "StringEquals" : {
            "ec2:CreateAction" : "RunInstances"
          }
        }
      }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "ec2_runner_role_policy_attachment" {
  role       = aws_iam_role.omics_test_ec2_runner_role.name
  policy_arn = aws_iam_policy.ec2_runner_iam_policy.arn
}

resource "aws_iam_role_policy_attachment" "ec2_runner_role_policy_attachment_for_omics" {
  role       = aws_iam_role.omics_test_ec2_runner_role.name
  policy_arn = aws_iam_policy.iam_policy.arn
}

data "aws_iam_openid_connect_provider" "github" {
  url = "https://token.actions.githubusercontent.com"
}

data "aws_iam_policy_document" "assume_role_policy" {
  statement {
    sid = "AssumeRole"
    actions = ["sts:AssumeRole"]

    principals {
      type = "Service"
      identifiers = ["ec2.amazonaws.com"]
    }
  }
  statement {
    sid = "GithubOID"

    actions = [
      "sts:AssumeRoleWithWebIdentity"
    ]

    effect = "Allow"

    principals {
      type = "Federated"
      identifiers = [data.aws_iam_openid_connect_provider.github.arn]
    }

    condition {
      test     = "StringLike"
      variable = "token.actions.githubusercontent.com:sub"
      values = ["repo:${var.github_org}/*"]
    }

    condition {
      test     = "ForAllValues:StringEquals"
      variable = "token.actions.githubusercontent.com:iss"
      values = ["https://token.actions.githubusercontent.com"]
    }

    condition {
      test     = "ForAllValues:StringEquals"
      variable = "token.actions.githubusercontent.com:aud"
      values = ["sts.amazonaws.com"]
    }
  }
}

resource "aws_iam_instance_profile" "github_actions_instance_profile" {
  count = var.aws_region == "us-east-1" ? 1 : 0 # create only once on primary region
  name = var.github_actions_instance_profile_name
  role = aws_iam_role.omics_deployment_role.name
}
