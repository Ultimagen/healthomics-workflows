version 1.0

struct GlobalVariables {
  String ubuntu_docker
  String gitc_jar_path
  String ug_gatk_picard_docker
  String broad_gatk_docker
  String ug_call_variants_docker
  String ug_make_examples_docker
  String perl_docker
  String bcftools_docker
  String monitoring_script
  String ref_cache_script
  String ua_docker
  String giraffe_docker
  String rematching_docker
  String trimmer_docker
  String star_docker
  String sorter_docker
  String gridss_docker
  String gripss_docker
  String single_cell_qc_docker
  String segdup_docker
  String ugbio_core_docker
  String ugbio_cnv_docker
  String hla_la_docker
  String ugbio_mrd_docker
  String ugbio_featuremap_docker
  String ugbio_srsnv_docker
  String ugbio_ppmseq_docker
  String ugbio_freec_docker
  String ugbio_filtering_docker
  String ugbio_comparison_docker
  String featuremap_docker
  String ug_jalign_docker
  String mosdepth_docker
}
workflow Globals {
  input {
  GlobalVariables glob ={
        "ubuntu_docker": "ubuntu:focal",
        "gitc_jar_path": "/usr/gitc/",
        "ug_gatk_picard_docker": "ultimagenomics/ug_gatk_picard:0.15",
        "broad_gatk_docker": "broadinstitute/gatk:4.6.0.0",
        "ug_call_variants_docker": "ultimagenomics/call_variants:2.2.3",
        "ug_make_examples_docker": "ultimagenomics/make_examples:3.1.8",
        "perl_docker": "perl:5.38",
        "bcftools_docker": "staphb/bcftools:1.19",
        "monitoring_script": "s3://ultimagen-workflow-resources-us-east-1/monitor_1.0.sh",
        "ref_cache_script": "s3://ultimagen-workflow-resources-us-east-1/scripts/seq_cache_populate.pl",
        "ua_docker": "ultimagenomics/alignment:3.0.3",
        "giraffe_docker": "ultimagenomics/giraffe:1.47.1",
        "rematching_docker": "ultimagenomics/rematcher:v1.0.0",
        "trimmer_docker": "ultimagenomics/trimmer:2.3.3",
        "star_docker": "ultimagenomics/star:2.7.10a",
        "sorter_docker": "ultimagenomics/sorter:1.4.14",
        "gridss_docker": "ultimagenomics/gridss:0c97dd1",
        "gripss_docker": "ultimagenomics/gripss:ug_2.4.1_165b492",
        "single_cell_qc_docker": "ultimagenomics/ugbio_single_cell:1.14.0",
        "segdup_docker": "ultimagenomics/parascopy:1.1.2_3f14bd4",
        "ugbio_core_docker": "ultimagenomics/ugbio_core:1.12.0",
        "ugbio_cnv_docker": "ultimagenomics/ugbio_cnv:1.14.0",
        "hla_la_docker": "ultimagenomics/ugbio_hla_la:1.14.0",
        "ugbio_mrd_docker": "ultimagenomics/ugbio_mrd:1.15.0",
        "ugbio_featuremap_docker": "ultimagenomics/ugbio_featuremap:1.15.0",
        "ugbio_srsnv_docker": "ultimagenomics/ugbio_srsnv:1.15.0",
        "ugbio_ppmseq_docker": "ultimagenomics/ugbio_ppmseq:1.14.0",
        "ugbio_freec_docker": "ultimagenomics/ugbio_freec:1.14.0",
        "ugbio_filtering_docker": "ultimagenomics/ugbio_filtering:1.14.0",
        "ugbio_comparison_docker": "ultimagenomics/ugbio_comparison:1.12.0",
        "featuremap_docker": "ultimagenomics/featuremap:1.0.0_8797fc3",
        "ug_jalign_docker": "ultimagenomics/jalign:1.2.1",
        "mosdepth_docker": "quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2"
}
}

  output {
    GlobalVariables global_dockers = glob
  }
}