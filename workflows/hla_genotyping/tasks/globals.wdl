version 1.0

struct GlobalVariables {
  String gitc_jar_path
  String broad_gatk_docker
  String monitoring_script
  String ref_cache_script
  String hla_la_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "gitc_jar_path": "/usr/gitc/",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "hla_la_docker": "ultimagenomics/ugbio_hla_la:1.12.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}