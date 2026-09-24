version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File calling_interval_list_without_artefacts
  File? ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File sv_blacklist
  File ua_index
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "b37": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/wgs_calling_regions.v1.nobl.interval_list",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/hg19/hg19-blacklist.v2.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/UA/b19-v45-79372c0.uai"
},
      "hg38": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/sv/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
},
      "hg38_nist_v3_with_decoy": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1.fasta.fai",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3_hs38d1/GRCh38_GIABv3_hs38d1-v45-79372c0.uai"
},
      "hg38_no_alt": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/empty_file",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/sv/gridss/ENCFF356LFX.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/hg38_no_alt/hg38_no_alt_v45.uai"
},
      "mm39": {
        "calling_interval_list_without_artefacts": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.no_chrY.interval_list",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa.fai",
        "sv_blacklist": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/mm39.excluderanges.bed",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa.uai"
}
    }
  }
}
