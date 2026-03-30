version 1.0

# Auto-generated genome resources file
# Do not edit manually - regenerate using: wdls make-genome-resources target=<target>

struct GenomeResources {
  File? ref_alt
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File? trim_align_sort_coverage_intervals
  File? ua_index
  File? ua_meth_index_c2t
  File? ua_meth_index_g2a
}

workflow GenomeResourcesWorkflow {
  output {
    Map[String, GenomeResources] resources = {
      "b37": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg19.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/UA/b19-v45-79372c0.uai"
},
      "b37_ancient_dna": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg19/Homo_sapiens_assembly19_1000genomes_decoy_with_variants_and_pathogens/hs37d5_rsrs_pathogen.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg19/Homo_sapiens_assembly19_1000genomes_decoy_with_variants_and_pathogens/hs37d5_rsrs_pathogen.fa",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/Homo_sapiens_assembly19_1000genomes_decoy_with_variants_and_pathogens/hs37d5_rsrs_pathogen.fa.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg19/Homo_sapiens_assembly19_1000genomes_decoy_with_variants_and_pathogens/hs37d5_rsrs_pathogen-v45-79372c0.uai"
},
      "hg38": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
},
      "hg38_nist_v3": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3/GRCh38_GIABv3.fasta.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3/GRCh38_GIABv3.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3/GRCh38_GIABv3.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3/GRCh38_GIABv3.fasta.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v3/hg38_giab_v3-v45-79372c0.uai"
},
      "hg38_rna_seq": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/rna-seq/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/rna-seq/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/rna-seq/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai"
},
      "hg38_taps": {
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.fai",
        "trim_align_sort_coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.uai",
        "ua_meth_index_c2t": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.uai.c2t",
        "ua_meth_index_g2a": "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/251015/hg38_Lambda_pUC19.fasta.uai.g2a"
},
      "mm10": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/mm10/mm10.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/mm10/mm10.fa",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/mm10/mm10.fa.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/mm10/mm10-v45.uai"
},
      "mm10_methyl": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/mm10/methyl_seq_ref/GRCm38.primary_assembly.genome.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/mm10/methyl_seq_ref/GRCm38.primary_assembly.genome.fa",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/mm10/methyl_seq_ref/GRCm38.primary_assembly.genome.fa.fai",
        "ua_meth_index_c2t": "s3://ultimagen-workflow-resources-us-east-1/mm10/methyl_seq_ref/GRCm38.primary_assembly.genome.v45.uai.c2t",
        "ua_meth_index_g2a": "s3://ultimagen-workflow-resources-us-east-1/mm10/methyl_seq_ref/GRCm38.primary_assembly.genome.v45.uai.g2a"
},
      "mm39": {
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.dict",
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa.fai",
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/bioinfo-resources/tools/mouse_GRCm39_M31/GRCm39.primary_assembly.genome.fa.uai"
}
    }
  }
}
