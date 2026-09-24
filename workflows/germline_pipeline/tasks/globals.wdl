version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String ug_gatk_picard_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String bcftools_docker
  String monitoring_script
  String ua_docker
  String giraffe_docker
  String rematching_docker
  String gridss_docker
  String gripss_docker
  String segdup_docker
  String ugbio_core_docker
  String ugbio_cnv_docker
  String hla_la_docker
  String t1k_docker
  String ugbio_filtering_docker
  String pypgx_docker
  String str_genotyper_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.16_fixup",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:4.1.2",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.4.1",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ua_docker": "ultimagenomics/alignment:4.1.4",
        "giraffe_docker": "ultimagenomics/giraffe:1.74.0-r1",
        "rematching_docker": "ultimagenomics/rematcher:1.1.2_08f0df1",
        "gridss_docker": "ultimagenomics/gridss:1.0.2",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.1_165b492",
        "segdup_docker": "ultimagenomics/parascopy:1.3.2_f55b07e",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.31.0",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.30.0",
        "hla_la_docker": "ultimagenomics/hla_la:f02c77c",
        "t1k_docker": "ultimagenomics/ugbio_t1k:1.30.0",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.31.0",
        "pypgx_docker": "ultimagenomics/ugbio_pypgx:0.27.0-r2",
        "str_genotyper_docker": "ultimagenomics/str_genotyper:1.1.0_da93d4a"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}