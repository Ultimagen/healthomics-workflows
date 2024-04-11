version 1.0

import "structs.wdl" as Structs

task ApplyQualFilters {
  input{
    File input_vcf
    File input_vcf_index
    Int min_variant_quality_hmer_indels
    Int min_variant_quality_exome_hmer_indels
    Int min_variant_quality_non_hmer_indels
    Int min_variant_quality_snps
    String final_vcf_base_name
    File monitoring_script
    String gatk_docker
    String gitc_path
  }
  Int disk_size = ceil(3 * size(input_vcf, "GB") + 1 )
  String output_file = "~{final_vcf_base_name}.annotated.qual_filters.vcf.gz"
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    set -eo pipefail

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms10000m \
    -jar ~{gitc_path}/GATK_ultima.jar VariantFiltration \
      -V ~{input_vcf} \
      -O ~{output_file} \
      -filter "QUAL < ~{min_variant_quality_exome_hmer_indels} and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and vc.hasAttribute('EXOME')" \
      --filter-name LowQualHmerIndelsInExome \
      -filter "QUAL < ~{min_variant_quality_hmer_indels} and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and not vc.hasAttribute('EXOME')" \
      --filter-name LowQualHmerIndels \
      -filter "QUAL < ~{min_variant_quality_non_hmer_indels} and VARIANT_TYPE=='non-h-indel' and not vc.isFiltered()" \
      --filter-name LowQualNonHmerIndels \
      -filter "QUAL < ~{min_variant_quality_snps} and VARIANT_TYPE=='snp' and not vc.isFiltered()" \
      --filter-name LowQualSNVs 

  >>>
  runtime {
    memory: "15 GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: gatk_docker
    noAddress: true
    maxRetries: 1
  }
   output {
    File monitoring_log = "monitoring.log"
    File output_vcf = output_file
    File output_vcf_index = "~{output_file}.tbi"
  }
}

task ApplyAlleleFrequencyRatioFilter {
  input{
    File input_vcf
    Float af_ratio
    String final_vcf_base_name
    File monitoring_script
    String bcftools_docker
  }
  Int disk_size = ceil(3 * size(input_vcf, "GB") + 1 )
  String output_file = "~{final_vcf_base_name}.annotated.filt.vcf.gz"
  String tmp_file = "~{final_vcf_base_name}.tmp.vcf.gz"
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
 
    set -eo pipefail
 
    bcftools filter ~{input_vcf} -e '(VARIANT_TYPE="snp" || VARIANT_TYPE="non-h-indel") && (AD[0:1]/DP)/(BG_AD[0:1]/BG_DP) < ~{af_ratio}' -s "LowAFRatioToBackground" -m "+" -Oz -o  ~{output_file}
    bcftools index -t ~{output_file}
   >>>
  runtime {
    memory: "15 GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: bcftools_docker
    noAddress: true
    maxRetries: 1
  }
   output {
    File monitoring_log = "monitoring.log"
    File output_vcf = output_file
    File output_vcf_index = "~{output_file}.tbi"
  }
}

task RemoveRefCalls {
  input{
    File input_vcf
    String final_vcf_base_name
    File monitoring_script
    String bcftools_docker
  }
  Int disk_size = ceil(3 * size(input_vcf, "GB") + 1 )
  String output_file = "~{final_vcf_base_name}.annotated.filt.vcf.gz"

  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    set -eo pipefail

    bcftools view -i 'FILTER!="RefCall"' ~{input_vcf} -Oz -o ~{output_file}
    bcftools index -t ~{output_file}

  >>>
  runtime {
    memory: "15 GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: bcftools_docker
    noAddress: true
    maxRetries: 1
  }
   output {
    File monitoring_log = "monitoring.log"
    File output_vcf = output_file
    File output_vcf_index = "~{output_file}.tbi"
  }
}

task CalibrateBridgingSnvs { 
    input{
        File input_vcf
        File input_vcf_index
        References references
        String final_vcf_base_name
        File monitoring_script
        String ugvc_docker
    }

    Int disk_size = ceil(3 * size(input_vcf, "GB") + size(references.ref_fasta, "GB") + 1 )
    String output_file = "~{final_vcf_base_name}.fix_bridge_snvs.vcf.gz"
    
    command <<<
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      source ~/.bashrc
      conda activate genomics.py3
      set -eo pipefail
      
      python /VariantCalling/ugvc calibrate_bridging_snvs \
      --vcf ~{input_vcf} \
      --reference ~{references.ref_fasta} \
      --output ~{output_file} 

    >>>
    runtime {
      memory: "15 GB"
      disks: "local-disk " + ceil(disk_size) + " HDD"
      docker: ugvc_docker
      noAddress: true
      maxRetries: 1
    }
    output {
      File monitoring_log = "monitoring.log"
      File output_vcf = output_file
      File output_vcf_index = "~{output_file}.tbi"
  }
}
