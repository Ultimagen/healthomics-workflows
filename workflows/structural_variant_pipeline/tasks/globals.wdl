version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String ug_gatk_picard_docker
  String broad_gatk_docker
  String ug_make_examples_docker
  String monitoring_script
  String ua_docker
  String giraffe_docker
  String rematching_docker
  String gridss_docker
  String gripss_docker
  String ugbio_core_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.16",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.1.10",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ua_docker": "ultimagenomics/alignment:3.0.6",
        "giraffe_docker": "ultimagenomics/giraffe:1.47.1",
        "rematching_docker": "ultimagenomics/rematcher:main_28d1f2e",
        "gridss_docker": "ultimagenomics/gridss:0c97dd1",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.1_165b492",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}