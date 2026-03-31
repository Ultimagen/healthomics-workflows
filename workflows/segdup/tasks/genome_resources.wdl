version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File efficient_dv_target_intervals
  File exome_intervals
  File ref_dict
  File ref_fasta
  File ref_fasta_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "hg38": {
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/wgs_calling_regions.hg38.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
}
    }
  }
}
