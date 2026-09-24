version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File calling_interval_list_without_artefacts
  File cnv_normalization_cohort
  File coverage_collection_interval_list
  File efficient_dv_target_intervals
  File exome_intervals
  File par_regions
  File ploidy_exclude_regions
  File ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File roh_blacklist
  File? srsnv_training_interval_list
  File sv_blacklist
  File trim_align_sort_coverage_intervals
  File ua_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "hg38": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/sv/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v3.0/solaris2_60samples_cohort_v3.0.hg38.ReadsCount.rds",
        "coverage_collection_interval_list": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/wgs_coverage_regions.hg38.interval_list",
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/wgs_calling_regions.hg38.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed",
        "ploidy_exclude_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/ploidy/chrY_low_confidence_regions.hg38.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "roh_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38-blacklist.v2.bed",
        "srsnv_training_interval_list": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/wgs_calling_regions.without_encode_blacklist.hg38.chr1_22.interval_list",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
},
      "hg38_nist_v3_with_decoy": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/germline_CNV_cohort/v3.0/solaris2_60samples_cohort_v3.0.GIABv3_hs38d1.ReadsCount.rds",
        "coverage_collection_interval_list": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/wgs_coverage_regions.hg38.interval_list",
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/wgs_calling_regions.hg38.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/annotation_intervals/exome.twist.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed",
        "ploidy_exclude_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/ploidy/chrY_low_confidence_regions.hg38.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta.fai",
        "roh_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38-blacklist.v2.bed",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1-v45-79372c0.uai"
}
    }
  }
}
