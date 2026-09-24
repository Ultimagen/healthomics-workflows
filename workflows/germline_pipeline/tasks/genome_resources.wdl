version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File? calling_interval_list_without_artefacts
  File? cnv_normalization_cohort
  File efficient_dv_target_intervals
  File exome_intervals
  File par_regions
  File? ploidy_exclude_regions
  File ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File? roh_blacklist
  File? sv_blacklist
  File ua_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "b37": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/wgs_calling_regions.v1.nobl.interval_list",
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg19/germline_CNV_cohort/v3.0/solaris2_60samples_cohort_v3.0.b37.ReadsCount.rds",
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/wgs_calling_regions.v1.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg19/annotation_intervals/exome.twist.hg19.sort.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/par_regions.b37.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
        "roh_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg19/hg19-blacklist.v2.bed",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg19/hg19-blacklist.v2.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/UA/b19-v45-79372c0.uai"
},
      "hg38": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/sv/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/germline_CNV_cohort/v3.0/solaris2_60samples_cohort_v3.0.hg38.ReadsCount.rds",
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/wgs_calling_regions.hg38.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed",
        "ploidy_exclude_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/ploidy/chrY_low_confidence_regions.hg38.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "roh_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38-blacklist.v2.bed",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
},
      "hg38_nist_v3_with_decoy": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "cnv_normalization_cohort": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/germline_CNV_cohort/v3.0/solaris2_60samples_cohort_v3.0.GIABv3_hs38d1.ReadsCount.rds",
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
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1-v45-79372c0.uai"
},
      "hg38_no_alt": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/wgs_calling_regions.hg38_no_alt.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed",
        "ploidy_exclude_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/ploidy/chrY_low_confidence_regions.hg38.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/empty_file",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        "roh_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38-blacklist.v2.bed",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/hg38_no_alt_v45.uai"
},
      "hg38_taps": {
        "efficient_dv_target_intervals": "s3://ultimagen-workflow-resources-us-east-1/deepvariant/regions/wgs_calling_regions.hg38_Lambda_pUC19.interval_list",
        "exome_intervals": "s3://ultimagen-workflow-resources-us-east-1/hg38/annotation_intervals/exome.twist.bed",
        "par_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/par_regions.hg38.bed",
        "ploidy_exclude_regions": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/ploidy/chrY_low_confidence_regions.hg38.bed",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.uai"
}
    }
  }
}
