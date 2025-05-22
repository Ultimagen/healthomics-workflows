resource "aws_iam_role" "omics_role" {
  name               = "${var.project}-omics-role-${var.aws_region}"
  assume_role_policy = data.aws_iam_policy_document.assume_role_policy.json
}

data "aws_iam_policy_document" "assume_role_policy" {
  statement {
    sid = "OmicsAssumeRole"
    actions = ["sts:AssumeRole"]

    principals {
      type = "Service"
      identifiers = ["omics.amazonaws.com"]
    }

    condition {
      test     = "StringEquals"
      variable = "aws:SourceAccount"

      values = [
        var.aws_account_id
      ]
    }

    condition {
      test     = "ArnLike"
      variable = "aws:SourceArn"

      values = [
        "arn:aws:omics:${var.aws_region}:${var.aws_account_id}:run/*"
      ]
    }
  }

  statement {
    sid = "EC2AssumeRole"
    actions = ["sts:AssumeRole"]

    principals {
      type = "Service"
      identifiers = ["ec2.amazonaws.com"]
    }
  }
}


data "aws_iam_policy_document" "read_only_buckets" {
  statement {
    sid = "LisBuckets"

    actions = [
      "s3:ListBucket"
    ]

    resources = formatlist("arn:aws:s3:::%s",
      flatten([
        var.omics_inputs_bucket,
        var.omics_outputs_bucket,
        var.omics_cache_bucket,
        var.bioinfo_resources_bucket,
        var.accessible_buckets_list,
        var.broad_references_public_bucket
      ]))
  }
  statement {
    sid = "GetObjects"

    actions = [
      "s3:GetObject"
    ]

    resources = formatlist("arn:aws:s3:::%s/*",
      flatten([
        var.omics_inputs_bucket,
        var.omics_outputs_bucket,
        var.omics_cache_bucket,
        var.bioinfo_resources_bucket,
        var.accessible_buckets_list,
        var.broad_references_public_bucket
      ]))
  }
}

resource "aws_iam_policy" "read_only_buckets_policy" {
  name   = "UltimagenOmicsReadBuckets-policy-${var.aws_region}"
  path   = "/"
  policy = data.aws_iam_policy_document.read_only_buckets.json
}

resource "aws_iam_policy" "omics_policy" {
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:PutObject"
            ],
            "Resource": [
                "arn:aws:s3:::${var.omics_outputs_bucket}/*",
                "arn:aws:s3:::${var.omics_cache_bucket}/*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "logs:DescribeLogStreams",
                "logs:CreateLogStream",
                "logs:PutLogEvents"
            ],
            "Resource": [
                "arn:aws:logs:${var.aws_region}:${var.aws_account_id}:log-group:/aws/omics/WorkflowLog:log-stream:*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "logs:CreateLogGroup"
            ],
            "Resource": [
                "arn:aws:logs:${var.aws_region}:${var.aws_account_id}:log-group:/aws/omics/WorkflowLog:*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "ecr:BatchGetImage",
                "ecr:GetDownloadUrlForLayer",
                "ecr:BatchCheckLayerAvailability"
            ],
            "Resource": [
                "arn:aws:ecr:${var.aws_region}:${var.aws_account_id}:repository/*",
                "arn:aws:ecr:${var.aws_region}:${var.aws_shared_account_id}:repository/*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "ecr:GetAuthorizationToken"
            ],
            "Resource": [
                "*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject"
            ],
            "Resource": [
                "arn:aws:s3:::sentieon-omics-license-${var.aws_region}",
                "arn:aws:s3:::sentieon-omics-license-${var.aws_region}/*",
                "arn:aws:s3:::omics-${var.aws_region}",
                "arn:aws:s3:::omics-${var.aws_region}/*"
            ]
        }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "omics_policy_attachment" {
  policy_arn = aws_iam_policy.omics_policy.arn
  role       = aws_iam_role.omics_role.name
}

resource "aws_iam_role_policy_attachment" "read_only_buckets_policy_attachment" {
  policy_arn = aws_iam_policy.read_only_buckets_policy.arn
  role       = aws_iam_role.omics_role.name
}

resource "aws_iam_instance_profile" "omics_role_instance_profile" {
  # for local execution with miniwdl
  name = aws_iam_role.omics_role.name
  role = aws_iam_role.omics_role.name
}