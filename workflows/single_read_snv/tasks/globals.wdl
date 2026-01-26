version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String monitoring_script
  String ugbio_core_docker
  String ugbio_featuremap_docker
  String ugbio_srsnv_docker
  String featuremap_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.18.0",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.20.0",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.18.0",
        "featuremap_docker": "ultimagenomics/featuremap:master_75f0535"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}