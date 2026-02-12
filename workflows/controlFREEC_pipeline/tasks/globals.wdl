version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String broad_gatk_docker
  String monitoring_script
  String ugbio_cnv_docker
  String ugbio_freec_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.20.1",
        "ugbio_freec_docker": "ultimagenomics/ugbio_freec:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}