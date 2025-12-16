version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String bcftools_docker
  String monitoring_script
  String ugbio_core_docker
  String ugbio_mrd_docker
  String mosdepth_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:main_7d72af6",
        "ugbio_mrd_docker": "ultimagenomics/ugbio_mrd:1.17.0",
        "mosdepth_docker": "quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}