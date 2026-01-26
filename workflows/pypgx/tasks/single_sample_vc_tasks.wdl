version 1.0
import "structs.wdl" as Structs


# Call variants on a single sample with HaplotypeCaller to produce a VCF
task HaplotypeCaller {
    input {
        Array[File] input_bam_list
        Array[File] input_bam_index_list
        File interval_list
        File? interval_list_index
        String vcf_basename
        References references
        Int preemptible_tries
        String docker
        Int memory_gb
        String? extra_args
        String? annotation_extra_args
        File monitoring_script
        Boolean no_address
        Boolean make_gvcf
        Boolean make_bamout
        Boolean with_alleles = false
        String gitc_path = "/usr/gitc/"
    }

    parameter_meta {
        input_bam_list: {
            localization_optional: true
        }
        input_bam_index_list: {
            localization_optional: true
        }
        interval_list: {
            localization_optional: true
        }
        interval_list_index: {
            localization_optional: true
        }
    }

    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String output_filename = vcf_basename + output_suffix
    Int VCF_disk_size = 10
    Float ref_size = size(references.ref_fasta, "GB") + size(references.ref_fasta_index, "GB") + size(references.ref_dict, "GB")
    Float additional_disk = 20
    Int disk_size = ceil(VCF_disk_size + ref_size + additional_disk)

    String alleles_str = if with_alleles then ("--alleles " + interval_list) else ""
    command {
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        touch ~{output_filename}.tbi
        touch ~{output_filename}
        touch realigned.bam
        touch realigned.bai

        java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{memory_gb-2}g \
          -jar ~{gitc_path}GATK_ultima.jar \
          HaplotypeCaller \
          -R ~{references.ref_fasta} \
          -O ~{output_filename} \
          -I ~{sep = ' -I ' input_bam_list} \
          --intervals ~{interval_list} \
          --smith-waterman FASTEST_AVAILABLE \
          ~{alleles_str} \
          ~{true="-ERC GVCF" false="" make_gvcf} \
          ~{true="--bamout realigned.bam" false="" make_bamout} \
          ~{extra_args} ~{annotation_extra_args}
    }

    runtime {
        preemptible: preemptible_tries
        memory: "~{memory_gb} GB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        continueOnReturnCode: [0,134,139]
        bootDiskSizeGb: 15
        maxRetries: 1
    }
    output {
        File output_vcf = "~{output_filename}"
        File output_vcf_index = "~{output_filename}.tbi"
        File monitoring_log = "monitoring.log"
        File bamout = "realigned.bam"
        File bamout_index = "realigned.bai"
    }
}


task MoveAnnotationsToGvcf {
 input {
 File monitoring_script
 File filtered_vcf
 File filtered_vcf_index
 File gvcf
 File gvcf_index
 String annotation = "TREE_SCORE"
 }
 String filename = basename(gvcf, ".g.vcf.gz")
 Int disk_size = ceil(size(filtered_vcf, "GB") + size(gvcf, "GB") * 2 + 30)

 command {
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

     bcftools annotate -a ~{filtered_vcf} -c ~{annotation} ~{gvcf} -o ~{filename}.annotated.g.vcf.gz -O z
     tabix -p vcf ~{filename}.annotated.g.vcf.gz
 }
 output {
     File output_gvcf = "~{filename}.annotated.g.vcf.gz"
     File output_gvcf_index = "~{filename}.annotated.g.vcf.gz.tbi"
     File monitoring_log = "monitoring.log"
 }

runtime {
     docker: "us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0"
     memory: "7 GiB"
     disks: "local-disk " + disk_size + " HDD"
	}
}

task ConvertGVCFtoVCF {
  input {
    File monitoring_script
    File input_gvcf
    File input_gvcf_index
    String output_vcf_name
    References references
    Int disk_size = ceil(size(input_gvcf, "GB") + 5)
    Int preemptible_tries
    Boolean no_address
    String docker
    String gitc_path = "/usr/gitc/"
  }
  command {
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms10000m \
    -jar ~{gitc_path}GATK_ultima.jar GenotypeGVCFs \
    -R ~{references.ref_fasta} \
    -V ~{input_gvcf} \
    -O ~{output_vcf_name} \
    -A  StrandBiasBySample
  }
  runtime {
    preemptible: preemptible_tries
    memory: "12 GB"
    cpu: "1"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    noAddress: no_address
    maxRetries: 1
  }
    output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task CreateSECBlacklist {
   input {
      File monitoring_script
      File input_gvcf
      File input_gvcf_index
      File? blacklist_file
      Array[File]? sec_models
      String output_blacklist_path
      Int disk_size = ceil(2 * size(input_gvcf, "GB") +
                size(blacklist_file, "GB") +
                size(select_first([ sec_models, []]), "GB") + 20)
      Int preemptible_tries
      Boolean no_address
      String docker
  }
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

    bedtools intersect -a ~{input_gvcf} -b ~{blacklist_file} -header | uniq > tmp.g.vcf

    correct_systematic_errors \
      --model ~{sep=' --model ' sec_models} \
      --gvcf tmp.g.vcf \
      --relevant_coords ~{blacklist_file} \
      --novel_detection_only \
      --output_file ~{output_blacklist_path}

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "12 GB"
    cpu: "1"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    noAddress: no_address
    maxRetries: 1
  }
    output {
    File monitoring_log = "monitoring.log"
    File output_blacklist = "~{output_blacklist_path}"

  }
}

task PrepareTrainingSet {
  input { 
      File input_vcf
      File input_vcf_index
      String base_file_name 
      File? blacklist_file
      File? base_vcf
      File? base_vcf_index
      File? hcr
      File? ref_fasta
      File? ref_sdf
      String test_chromosome
      Array[String] train_and_test_chromosomes
      Array[String]? custom_annotations
      Int preemptible_tries
      String docker
      File monitoring_script
      Boolean no_address
      String gt_type = "approximate"
  }
  Int disk_size = ceil(2*size(input_vcf, "GB") + size(blacklist_file, "GB") + 20)
  Boolean have_sdf = defined(ref_sdf)

  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail
    
    if [ ~{have_sdf} = true ]
    then
      python -m tarfile -e ~{ref_sdf} _reference.sdf
    fi

    training_prep_pipeline \
              --call_vcf ~{input_vcf} \
              ~{"--blacklist " + blacklist_file} \
              ~{"--base_vcf " + base_vcf} \
              ~{"--hcr " + hcr} \
              ~{"--reference " + ref_fasta} \
              "--reference_sdf _reference.sdf" \
              --contigs_to_read ~{sep=" " train_and_test_chromosomes} \
              --contig_for_test ~{test_chromosome} \
              ~{true="--custom_annotations " false="" defined(custom_annotations)} ~{sep=" --custom_annotations " custom_annotations} \
              --gt_type ~{gt_type} \
              --output_prefix ~{base_file_name}
    >>>
  runtime {
    preemptible: preemptible_tries
    memory: "16 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File monitoring_log = "monitoring.log"
    File output_training_set = "~{base_file_name}.h5"
    File output_test_set = "~{base_file_name}_test.h5"
  }
}

task TrainModel {
  input {
    File monitoring_script
    Array[File] train_data
    Array[File] test_data
    String base_file_name
    String gt_type = "approximate"
    Array[String]? custom_annotations
    Int preemptible_tries
    Int disk_size = ceil(size(train_data, "GB") +
                         size(test_data, "GB") + 20)
    String docker
    Boolean no_address
  }
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

    train_models_pipeline \
            --train_dfs ~{sep =" " train_data} \
            --test_dfs  ~{sep =" " test_data} \
            --output_file_prefix ~{base_file_name}.model \
            --gt_type ~{gt_type} \
            ~{true="--custom_annotations " false="" defined(custom_annotations)} ~{sep=" --custom_annotations " custom_annotations} \
    >>>
  runtime {
    preemptible: preemptible_tries
    memory: "64 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File model_h5 = "~{base_file_name}.model.h5"
    File model_pkl = "~{base_file_name}.model.pkl"
    File monitoring_log = "monitoring.log"
  }
}

task FilterVCF {
  input {
    File monitoring_script
    File input_vcf
    File input_vcf_index
    File? input_model
    Boolean filter_cg_insertions
    Int? decision_threshold
    Boolean? overwrite_quality
    File? ref_fasta
    File? ref_fasta_idx
    File? blacklist_file
    String final_vcf_base_name
    Array[String]? custom_annotations
    Boolean recalibrate_gt
    Int disk_size = 10
    Int preemptible_tries
    String docker
    Boolean no_address
  }

  Boolean should_overwrite_quality = defined(overwrite_quality) && select_first([overwrite_quality])
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

    filter_variants_pipeline --input_file ~{input_vcf} \
                                       ~{"--model_file " + input_model} \
                                       ~{true="--blacklist_cg_insertions" false="" filter_cg_insertions} \
                                       ~{"--blacklist " + blacklist_file} \
                                       ~{"--ref_fasta " + ref_fasta} \
                                       ~{true="--custom_annotations " false="" defined(custom_annotations)}~{sep=" --custom_annotations " custom_annotations} \
                                       ~{true="--recalibrate_genotype --treat_multiallelics" false="" recalibrate_gt} \
                                       ~{"--decision_threshold " + decision_threshold} \
                                       ~{true="--overwrite_qual_tag" false="" should_overwrite_quality} \
                                       --output_file ~{final_vcf_base_name}.filtered.vcf.gz
                                       
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "64 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File monitoring_log = "monitoring.log"
    File output_vcf_filtered = "~{final_vcf_base_name}.filtered.vcf.gz"
    File output_vcf_filtered_index = "~{final_vcf_base_name}.filtered.vcf.gz.tbi"
  }
}

task CompressAndIndexVCF {
  input { 
    File input_vcf
    String base_file_name 
    File monitoring_script
    String docker
    Boolean no_address
    Int preemptible_tries
    Int disk_size = ceil(2*size(input_vcf,"GB")) 
    Int cpus
  }
  command <<< 
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail
    source /opt/conda/etc/profile.d/conda.sh
    conda activate genomics.py3

    bcftools view --threads ~{cpus} -Oz -o ~{base_file_name}.vcf.gz ~{input_vcf}
    bcftools index -t ~{base_file_name}.vcf.gz
  >>>

  output { 
    File output_vcf = "~{base_file_name}.vcf.gz"
    File output_vcf_index = "~{base_file_name}.vcf.gz.tbi"
    File monitoring_log = "monitoring.log"
  }

  runtime {
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }

}


task NormalizeVariants {
  input{
    File vcf_file
    File vcf_file_index
    References references
    Boolean output_differences
    File monitoring_script
    Int preemptible_tries
    Boolean no_address
    String docker
  }
  Int disk_size = ceil(2 * size(vcf_file, "GB") + 10 )
  String output_filename = basename(vcf_file, ".vcf.gz") + ".norm.vcf.gz"
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

    source /opt/conda/etc/profile.d/conda.sh
    conda activate genomics.py3

    bcftools norm -f ~{references.ref_fasta} ~{vcf_file} --threads 2 -O z -o ~{output_filename} > norm_output_stats.txt 2>&1
    cat norm_output_stats.txt
    bcftools index -t ~{output_filename}

    if [ ~{true="true" false="false" output_differences}=true ]
    then
      bcftools isec --complement -p orig_not_norm ~{vcf_file} ~{output_filename}
      mv orig_not_norm/0000.vcf orig_not_norm.vcf
    fi

    touch orig_not_norm.vcf
 
  >>>
  runtime {
    memory: "3 GB"
    preemptible: preemptible_tries
    noAddress: no_address
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }
    output {
    File monitoring_log = "monitoring.log"
    File normalized_vcf = '~{output_filename}'
    File normalized_vcf_index = '~{output_filename}.tbi'
    File norm_output_stats = "norm_output_stats.txt"
    File orig_not_norm = "orig_not_norm.vcf"
  }
}