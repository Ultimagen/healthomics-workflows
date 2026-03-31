version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File trim_align_sort_coverage_intervals
  File ua_index
  File? ua_meth_index_c2t
  File? ua_meth_index_g2a
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "hg38": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
}
    }
  }
}
