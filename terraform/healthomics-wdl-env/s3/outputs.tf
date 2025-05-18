output "omics_inputs_bucket" {
  value = aws_s3_bucket.pipelines-input-bucket.id
}

output "omics_outputs_bucket" {
  value = aws_s3_bucket.pipelines-output-bucket.id
}

output "omics_cache_bucket" {
  value = aws_s3_bucket.pipelines-cache-bucket.id
}

output "bioinfo_resources_bucket" {
  value = aws_s3_bucket.bioinfo-resources-bucket.id
}