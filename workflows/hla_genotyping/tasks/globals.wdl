version 1.0

struct GlobalVariables {
  String broad_gatk_docker
  String monitoring_script
  String hla_la_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "hla_la_docker": "ultimagenomics/hla_la:f02c77c"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}