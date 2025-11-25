version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String perl_docker
  String monitoring_script
  String ref_cache_script
  String ua_docker
  String trimmer_docker
  String star_docker
  String sorter_docker
  String single_cell_qc_docker
  String ugbio_core_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "perl_docker": "perl:5.38",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ua_docker": "ultimagenomics/alignment:3.0.6",
        "trimmer_docker": "ultimagenomics/trimmer:2.3.4",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.4.15",
        "single_cell_qc_docker": "ultimagenomics/ugbio_single_cell:1.14.0",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.16.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}