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
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.16_fixup",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.3.4",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ua_docker": "ultimagenomics/alignment:3.0.8",
        "giraffe_docker": "ultimagenomics/giraffe:1.74.0-r1",
        "rematching_docker": "ultimagenomics/rematcher:1.1.2_08f0df1",
        "gridss_docker": "ultimagenomics/gridss:1.0.2",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.1_165b492",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.28.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}