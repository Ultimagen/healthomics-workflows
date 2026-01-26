version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String perl_docker
  String bcftools_docker
  String monitoring_script
  String ref_cache_script
  String ugbio_filtering_docker
  String pypgx_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.4",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.1.10",
        "perl_docker": "perl:5.38",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.20.0",
        "pypgx_docker": "ultimagenomics/ugbio_pypgx:0.26.0-r3"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}