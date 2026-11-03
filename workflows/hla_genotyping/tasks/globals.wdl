version 1.0

struct GlobalVariables {
  String monitoring_script
  String ugbio_core_docker
  String t1k_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.33.0",
        "t1k_docker": "ultimagenomics/ugbio_t1k:1.30.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}