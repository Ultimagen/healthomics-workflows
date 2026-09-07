version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String monitoring_script
  String ua_docker
  String trimmer_docker
  String star_docker
  String sorter_docker
  String single_cell_qc_docker
  String ugbio_core_docker
  String ug_vg_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ua_docker": "ultimagenomics/alignment:3.0.8",
        "trimmer_docker": "ultimagenomics/trimmer:2.3.7",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.6.11",
        "single_cell_qc_docker": "ultimagenomics/ugbio_single_cell:1.28.0",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.30.2",
        "ug_vg_docker": "ultimagenomics/vg:1.74.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}