version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File srsnv_training_interval_list
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "hg38": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "srsnv_training_interval_list": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/wgs_calling_regions.without_encode_blacklist.hg38.chr1_22.interval_list"
}
    }
  }
}
