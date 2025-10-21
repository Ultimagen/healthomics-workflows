version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_jar_path
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String bcftools_docker
  String monitoring_script
  String ref_cache_script
  String ugbio_filtering_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_jar_path": "/usr/gitc/",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.3",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.1.6",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.12.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}