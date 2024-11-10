version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_docker
  String gitc_jar_path
  String conda_docker
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
  String subread_docker
  String starfusion_docker
  String ugbio_cnv_docker
  String vcflite_docker
  String ugbio_mrd_docker
  String ugbio_featuremap_docker
  String ugbio_srsnv_docker
  String ugbio_ppmseq_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_docker": "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.6-1599252698",
        "gitc_jar_path": "/usr/gitc/",
        "conda_docker": "continuumio/anaconda3:2019.10",
        "ug_vc_docker": "ultimagenomics/ugvc:0.25_760275b",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.14",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.2",
        "ug_make_examples_docker": "ultimagenomics/make_examples:2.2.5",
        "perl_docker": "perl:5.38",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ua_docker": "ultimagenomics/alignment:1.1.2",
        "trimmer_docker": "ultimagenomics/trimmer:2.1.3",
        "fastqc_docker": "quay.io/biocontainers/fastqc:0.11.9--0",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.2.6",
        "ug_control_freec_docker": "ultimagenomics/ug_control_freec:26fe532",
        "pigz_docker": "nsheff/pigz:latest",
        "gridss_docker": "ultimagenomics/gridss:test_f822ea2",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.0_cb91bca",
        "single_cell_qc_docker": "ultimagenomics/ugbio_single_cell:1.3.2",
        "segdup_docker": "ultimagenomics/parascopy:1.0_30e2e98",
        "arriba_docker": "uhrigs/arriba:2.4.0",
        "subread_docker": "us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1689097353",
        "starfusion_docker": "trinityctat/starfusion:1.13.0",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.3.0",
        "vcflite_docker": "ultimagenomics/ugbio_vcflite:1.3.0",
        "ugbio_mrd_docker": "ultimagenomics/ugbio_mrd:1.3.3",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.3.0",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.3.2",
        "ugbio_ppmseq_docker": "ultimagenomics/ugbio_ppmseq:1.3.0"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}