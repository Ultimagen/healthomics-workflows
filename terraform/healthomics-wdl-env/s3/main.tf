resource "aws_s3_bucket" "pipelines-input-bucket" {
  bucket = "${var.project}-${var.aws_account_id}-${var.aws_region}-${var.inputs_bucket_suffix}"

  tags = merge(
    var.custom_tags, {
      Project = var.project
    }
  )
}

resource "aws_s3_bucket" "pipelines-output-bucket" {
  bucket = "${var.project}-${var.aws_account_id}-${var.aws_region}-${var.outputs_bucket_suffix}"
  tags = merge(
    var.custom_tags, {
      Project = var.project
    }
  )
}

resource "aws_s3_bucket" "bioinfo-resources-bucket" {
  bucket = "${var.project}-${var.aws_account_id}-${var.aws_region}-${var.bioinfo_resources_bucket_suffix}"

  tags = merge(
    var.custom_tags, {
      Project = var.project
    }
  )
}

resource "aws_s3_bucket" "pipelines-cache-bucket" {
  bucket = "${var.project}-${var.aws_account_id}-${var.aws_region}-${var.cache_bucket_suffix}"

  tags = merge(
    var.custom_tags, {
      Project = var.project
    }
  )
}

locals {
  omics_buckets = {
    inputs            = aws_s3_bucket.pipelines-input-bucket.id,
    outputs           = aws_s3_bucket.pipelines-output-bucket.id,
    bioinfo_resources = aws_s3_bucket.bioinfo-resources-bucket.id
    cache             = aws_s3_bucket.pipelines-cache-bucket.id
  }
}

resource "aws_s3_bucket_public_access_block" "buckets_acl" {
  for_each = local.omics_buckets

  bucket = each.value

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_lifecycle_configuration" "inputs_buckets_lifecycle" {

  bucket = aws_s3_bucket.pipelines-input-bucket.id



  rule {
    id = "intelligent_tiering_7d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = 7
      storage_class = "INTELLIGENT_TIERING"
    }

  }


  rule {
    id = "glacier_${var.input_bucket_archive_days}d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = var.input_bucket_archive_days
      storage_class = "GLACIER"
    }
  }


  rule {
    id = "cleanup_incomplete_uploads"

    status = "Enabled"
    filter {
      prefix = ""
    }
    abort_incomplete_multipart_upload {
      days_after_initiation = 2
    }
  }

}

resource "aws_s3_bucket_lifecycle_configuration" "outputs_buckets_lifecycle" {

  bucket = aws_s3_bucket.pipelines-output-bucket.id

  rule {
    id = "intelligent_tiering_7d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = 7
      storage_class = "INTELLIGENT_TIERING"
    }
  }

  rule {
    id = "glacier_${var.output_bucket_archive_days}d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = var.output_bucket_archive_days
      storage_class = "GLACIER"
    }
  }

}

resource "aws_s3_bucket_lifecycle_configuration" "bioinfo_resources_buckets_lifecycle" {

  bucket = aws_s3_bucket.bioinfo-resources-bucket.id

  rule {
    id = "intelligent_tiering_7d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = 7
      storage_class = "INTELLIGENT_TIERING"
    }
  }

  rule {
    id = "cleanup_incomplete_uploads"

    status = "Enabled"
    filter {
      prefix = ""
    }
    abort_incomplete_multipart_upload {
      days_after_initiation = 2
    }
  }
}

resource "aws_s3_bucket_lifecycle_configuration" "cache_buckets_lifecycle" {

  bucket = aws_s3_bucket.pipelines-cache-bucket.id

  rule {
    id = "intelligent_tiering_7d"

    status = "Enabled"
    filter {
      prefix = ""
    }
    transition {
      days          = 7
      storage_class = "INTELLIGENT_TIERING"
    }
  }

  rule {
    id = "cleanup_incomplete_uploads"

    status = "Enabled"
    filter {
      prefix = ""
    }
    abort_incomplete_multipart_upload {
      days_after_initiation = 2
    }
  }

  rule {
    id = "expire_cache_on_failure_objects_after_${var.cache_on_failure_expiration_days}_days"

    status = "Enabled"
    filter {
      prefix = var.cache_on_failure_prefix
    }
    expiration {
      days = var.cache_on_failure_expiration_days
    }
  }

  rule {
    id = "expire_cache_always_objects_after_${var.cache_always_expiration_days}_days"

    status = "Enabled"

    filter {
      prefix = var.cache_always_prefix
    }
    expiration {
      days = var.cache_always_expiration_days
    }
  }
}

# resource "aws_s3_bucket_policy" "allow_write_access_from_another_account" {
#   bucket = aws_s3_bucket.pipelines-input-bucket.id
#   policy = data.aws_iam_policy_document.allow_write_access_from_another_account.json
# }
#
# data "aws_iam_policy_document" "allow_write_access_from_another_account" {
#   statement {
#     principals {
#       type = "AWS"
#       identifiers = [var.cross_aws_account_id]
#     }
#
#     actions = [
#       "s3:GetObject",
#       "s3:GetObjectTagging",
#       "s3:PutObject",
#       "s3:ListBucket",
#       "s3:ListBucketMultipartUploads",
#       "s3:AbortMultipartUpload",
#       "s3:PutObjectVersionAcl",
#       "s3:DeleteObject",
#       "s3:PutObjectAcl",
#       "s3:ListMultipartUploadParts"
#
#     ]
#
#     resources = [
#       aws_s3_bucket.pipelines-input-bucket.arn,
#       "${aws_s3_bucket.pipelines-input-bucket.arn}/*",
#     ]
#   }
# }

