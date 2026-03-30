version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ua_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "b37": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/UA/b19-v45-79372c0.uai"
},
      "hg38": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
}
    }
  }
}
