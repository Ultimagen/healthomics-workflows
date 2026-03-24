version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String monitoring_script
  String ua_docker
  String trimmer_docker
  String star_docker
  String sorter_docker
  String ugbio_core_docker
  String ugbio_ppmseq_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ua_docker": "ultimagenomics/alignment:3.0.6",
        "trimmer_docker": "ultimagenomics/trimmer:2.3.4",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.5.18",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.18.0",
        "ugbio_ppmseq_docker": "ultimagenomics/ugbio_ppmseq:1.16.2"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}