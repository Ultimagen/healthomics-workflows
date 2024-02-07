version 1.0

task global {
  command {}
  runtime { docker: "ubuntu:focal" }
  output {
    String ubuntu_docker = "ubuntu:focal"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698"
    String gitc_jar_path = "/usr/gitc/"
    String ug_vc_docker = "ultimagenomics/ugvc:v0.19.1_ef2ec20"
    String ug_gatk_picard_docker = "ultimagenomics/ug_gatk_picard:0.12.2.1"
    String broad_gatk_docker = "broadinstitute/gatk:4.2.6.1"
    String ug_call_variants_docker = "ultimagenomics/call_variants:edv_2.0.0.1_08cd993"
    String ug_make_examples_docker = "ultimagenomics/make_examples:edv_2.0.0.1_08cd993"
    String bcftools_docker = "staphb/bcftools:1.18"    
    String monitoring_script = "s3://ultimagen-workflow-resources-us-east-1/monitoring.sh"
  }
}
