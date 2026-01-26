version 1.0

struct GlobalVariables {
  String monitoring_script
  String str_genotyper_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "str_genotyper_docker": "ultimagenomics/str_genotyper:1.0.1_f7ae986"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}