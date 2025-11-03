version 1.0

struct GlobalVariables {
  String monitoring_script
  String ugbio_cnv_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}