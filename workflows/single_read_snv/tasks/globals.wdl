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
        "ugbio_core_docker": "ultimagenomics/ugbio_core:main_9b97e48",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.17.2",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.16.1",
        "featuremap_docker": "ultimagenomics/featuremap:master_54894f4"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}