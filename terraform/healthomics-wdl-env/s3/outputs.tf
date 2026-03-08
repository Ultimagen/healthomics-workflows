output "omics_inputs_bucket" {
  description = "Name of the S3 bucket for workflow input files"
  value       = aws_s3_bucket.pipelines-input-bucket.id
}

output "omics_outputs_bucket" {
  description = "Name of the S3 bucket for workflow output files"
  value       = aws_s3_bucket.pipelines-output-bucket.id
}

output "omics_cache_bucket" {
  description = "Name of the S3 bucket for WDL execution cache"
  value       = aws_s3_bucket.pipelines-cache-bucket.id
}

output "bioinfo_resources_bucket" {
  description = "Name of the S3 bucket for bioinformatics reference data"
  value       = aws_s3_bucket.bioinfo-resources-bucket.id
}