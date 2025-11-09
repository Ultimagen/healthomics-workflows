version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String monitoring_script
  String ugbio_core_docker
  String ugbio_cnv_docker
  String ug_jalign_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.16.1",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:main_3d8a76b",
        "ug_jalign_docker": "ultimagenomics/jalign:1.2.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}