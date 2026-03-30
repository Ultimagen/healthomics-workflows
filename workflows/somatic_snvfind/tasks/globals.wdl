version 1.0

struct GlobalVariables {
  String broad_gatk_docker
  String monitoring_script
  String ugbio_featuremap_docker
  String featuremap_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.22.2",
        "featuremap_docker": "ultimagenomics/featuremap:master_b51ac53"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}