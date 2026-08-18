version 1.0

struct GlobalVariables {
  String broad_gatk_docker
  String monitoring_script
  String ugbio_featuremap_docker
  String ugbio_srsnv_docker
  String ugbio_deep_srsnv_docker
  String featuremap_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.28.0",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.29.0",
        "ugbio_deep_srsnv_docker": "ultimagenomics/ugbio_deep_srsnv:1.5.0",
        "featuremap_docker": "ultimagenomics/featuremap:1.2.0_d2dbc55"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}