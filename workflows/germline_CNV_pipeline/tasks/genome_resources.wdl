version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File cnv_normalization_cohort
  File ref_dict
  File ref_fasta
  File ref_fasta_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "b37": {
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg19/germline_CNV_cohort/pe_v1.0/solaris_max_pe_66samples_cohort_v1.0.b37.ReadsCount.rds",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
},
      "hg38": {
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/pe_v1.0/solaris_max_pe_66samples_cohort_v1.0.hg38.ReadsCount.rds",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
},
      "hg38_nist_v3_with_decoy": {
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/germline_CNV_cohort/pe_v1.0/solaris_max_pe_66samples_cohort_v1.0.GIABv3_hs38d1.ReadsCount.rds",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta.fai"
},
      "hg38_no_alt": {
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/germline_CNV_cohort/pe_v1.0/solaris_max_pe_66samples_cohort_v1.0.hg38_no_alt.ReadsCount.rds",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
}
    }
  }
}
