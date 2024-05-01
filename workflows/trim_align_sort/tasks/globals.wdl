version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_docker
  String gitc_jar_path  
  String ug_vc_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String bcftools_docker
  String ug_control_freec_docker
  String ug_gatk_picard_docker
  String monitoring_script
  String perl_docker
  String ref_cache_script
  String star_docker
  String ua_docker
  String trimmer_docker
  String fastqc_docker
  String sorter_docker
}

workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_docker": "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698",
        "gitc_jar_path": "/usr/gitc/", 
        "ug_vc_docker": "ultimagenomics/ugvc:v0.21_e373679",
        "broad_gatk_docker": "broadinstitute/gatk:4.2.6.1",
        "ug_call_variants_docker": "ultimagenomics/call_variants:edv_2.1.1_b0ca4ece",
        "ug_make_examples_docker": "ultimagenomics/make_examples:edv_2.1.1_b0ca4ece",
        "bcftools_docker": "staphb/bcftools:1.19",
        "ug_control_freec_docker": "ultimagenomics/ug_control_freec:1679a9",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.12.2.1",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitoring.sh",
        "perl_docker": "perl:5.38",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "ua_docker": "ultimagenomics/alignment:1.0.1",
        "trimmer_docker": "ultimagenomics/trimmer:1.0.1",
        "fastqc_docker": "quay.io/biocontainers/fastqc:0.11.9--0",
        "sorter_docker": "ultimagenomics/sorter:1.0.1",
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}
