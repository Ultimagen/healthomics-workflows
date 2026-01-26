version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String bcftools_docker
  String monitoring_script
  String ugbio_core_docker
  String ugbio_cnv_docker
  String ugbio_filtering_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.18.0",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.20.0",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.20.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}