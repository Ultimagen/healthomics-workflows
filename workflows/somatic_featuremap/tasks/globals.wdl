version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String broad_gatk_docker
  String monitoring_script
  String ugbio_core_docker
  String ugbio_featuremap_docker
  String ugbio_freec_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.18.0",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.20.0",
        "ugbio_freec_docker": "ultimagenomics/ugbio_freec:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}