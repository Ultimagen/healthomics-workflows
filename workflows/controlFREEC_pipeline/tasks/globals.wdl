version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_docker
  String gitc_jar_path
  String ug_vc_docker
  String ug_gatk_picard_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String perl_docker
  String bcftools_docker
  String monitoring_script
  String ref_cache_script
  String ua_docker
  String trimmer_docker
  String fastqc_docker
  String star_docker
  String sorter_docker
  String ug_control_freec_docker
  String pigz_docker
  String assembly_docker
  String gridss_docker
  String gripss_docker
  String reads_transformer_docker
  String segdup_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_docker": "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698",
        "gitc_jar_path": "/usr/gitc/",
        "ug_vc_docker": "ultimagenomics/ugvc:0.22_21bbfe8",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.13",
        "broad_gatk_docker": "broadinstitute/gatk:4.2.6.1",
        "ug_call_variants_docker": "ultimagenomics/call_variants:edv_2.2.0_173435b8",
        "ug_make_examples_docker": "ultimagenomics/make_examples:2.2.1",
        "perl_docker": "perl:5.38",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ua_docker": "ultimagenomics/alignment:1.1.2",
        "trimmer_docker": "ultimagenomics/trimmer:1.0.1",
        "fastqc_docker": "quay.io/biocontainers/fastqc:0.11.9--0",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.1.5",
        "ug_control_freec_docker": "ultimagenomics/ug_control_freec:1679a9",
        "pigz_docker": "nsheff/pigz:latest",
        "assembly_docker": "ultimagenomics/haplotype:1.0_db79395",
        "gridss_docker": "ultimagenomics/gridss:5ce69a8",
        "gripss_docker": "ultimagenomics/gripss:ug_2.3.5_8da5ab3",
        "reads_transformer_docker": "ultimagenomics/ug-reads-transformer:1.0",
        "segdup_docker": "ultimagenomics/parascopy:cf76c4b"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}