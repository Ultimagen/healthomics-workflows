version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_docker
  String gitc_jar_path
  String ug_vc_docker
  String ug_gatk_picard_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String perl_docker
  String bcftools_docker
  String monitoring_script
  String ref_cache_script
  String ua_docker
  String trimmer_docker
  String fastqc_docker
  String star_docker
  String sorter_docker
  String ug_control_freec_docker
  String pigz_docker
  String gridss_docker
  String gripss_docker
  String single_cell_qc_docker
  String segdup_docker
  String arriba_docker
  String starfusion_docker
  String ugbio_core_docker
  String ugbio_cnv_docker
  String ugbio_vcflite_docker
  String ugbio_mrd_docker
  String ugbio_featuremap_docker
  String ugbio_srsnv_docker
  String ugbio_ppmseq_docker
  String ugbio_freec_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_docker": "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698",
        "gitc_jar_path": "/usr/gitc/",
        "ug_vc_docker": "ultimagenomics/ugvc:0.26.1_3d5a467",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.14",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.2",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.0.2",
        "perl_docker": "perl:5.38",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ua_docker": "ultimagenomics/alignment:2.1.2",
        "trimmer_docker": "ultimagenomics/trimmer:2.1.3",
        "fastqc_docker": "quay.io/biocontainers/fastqc:0.11.9--0",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.2.14",
        "ug_control_freec_docker": "ultimagenomics/ug_control_freec:26fe532",
        "pigz_docker": "nsheff/pigz:latest",
        "gridss_docker": "ultimagenomics/gridss:fb1cbae",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.0_cb91bca",
        "single_cell_qc_docker": "ultimagenomics/ugbio_single_cell:1.5.1",
        "segdup_docker": "ultimagenomics/parascopy:1.0_30e2e98",
        "arriba_docker": "uhrigs/arriba:2.4.0",
        "starfusion_docker": "trinityctat/starfusion:1.13.0",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.5.2",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.4.1",
        "ugbio_vcflite_docker": "ultimagenomics/ugbio_vcflite:1.4.1",
        "ugbio_mrd_docker": "ultimagenomics/ugbio_mrd:1.5.1",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.4.1",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.5.1",
        "ugbio_ppmseq_docker": "ultimagenomics/ugbio_ppmseq:1.5.1",
        "ugbio_freec_docker": "ultimagenomics/ugbio_freec:1.4.1"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}