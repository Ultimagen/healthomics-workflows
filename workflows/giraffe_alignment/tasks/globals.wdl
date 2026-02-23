version 1.0

struct GlobalVariables {
  String ug_gatk_picard_docker
  String monitoring_script
  String giraffe_docker
  String crammer_docker
  String sorter_docker
  String ugbio_core_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.16",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "giraffe_docker": "ultimagenomics/giraffe:1.71.1",
        "crammer_docker": "ultimagenomics/crammer:1.3.1_1f194312_1f1943",
        "sorter_docker": "ultimagenomics/sorter:1.5.13",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.18.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}