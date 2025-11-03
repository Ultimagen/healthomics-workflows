version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String bcftools_docker
  String monitoring_script
  String segdup_docker
  String ugbio_filtering_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.4",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.1.9",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "segdup_docker": "ultimagenomics/parascopy:1.2.0_f42c9e4",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}