version 1.0
import "structs.wdl" as Structs


task UGMakeExamples{
  input {
    Array[File] cram_files
    Array[File] cram_index_files
    File interval
    Float total_number_of_shards
    Float overall_calling_regions_length
    Float genome_length
    References references

    # Background sample inputs
    Array[File] background_cram_files
    Array[File] background_cram_index_files

    Int min_base_quality
    Int pileup_min_mapping_quality
    Int min_read_count_snps
    Int min_read_count_hmer_indels
    Int min_read_count_non_hmer_indels
    Float min_fraction_snps
    Float min_fraction_hmer_indels
    Float min_fraction_non_hmer_indels
    Int candidate_min_mapping_quality
    Int max_reads_per_partition
    Int assembly_min_base_quality
    Boolean make_gvcf
    Float p_error = 0
    Boolean prioritize_alt_supporting_reads    
    Array[Int] optimal_coverages
    Boolean cap_at_optimal_coverage
    Int interval_nreads = 10000
    Boolean output_realignment = false
    String? ug_make_examples_extra_args

    String docker
    File monitoring_script
    Float memory
    Int cpu
    Int preemptible_tries
    String cloud_provider
  }
  # Estimate output_size that fits candidate generated parameters (assuming constant image size)
  #   More sensitive thresholds (such as used for somatic variant detection) yield more examples (images)
  #   which consume more disk-space, regardless of the input-size.
  Float min_fraction_hmer_indels_trimmed = if (min_fraction_hmer_indels < 0.12) then min_fraction_hmer_indels else 0.12
  Float min_fraction_non_hmer_indels_trimmed = if (min_fraction_non_hmer_indels < 0.12) then min_fraction_non_hmer_indels else 0.12
  Float min_fraction_snps_trimmed  = if (min_fraction_snps < 0.12) then min_fraction_snps else 0.12
  # We calculate the disk_size additivly on the 3 thresholds we have

  # min_fraction_hmer_indels
  Array[Float] min_fraction_hmer_indels_values = [15786, 7659, 3597, 2169, 1266, 879, 627, 487, 381, 306, 283, 243]
  Float min_fraction_hmer_indels_value = min_fraction_hmer_indels_values[
          if (min_fraction_hmer_indels_trimmed*100)-1 < 0
          then 0 else floor(min_fraction_hmer_indels_trimmed*100)-1]/1.5

  # min_fraction_hmer_indels
  Array[Float] min_fraction_snps_values = [9531, 5067, 2358, 1776, 1380, 1134, 948, 817, 708, 672, 650, 626]
  Float min_fraction_snps_value = min_fraction_snps_values[
        if (min_fraction_snps_trimmed*100)-1 < 0
        then 0 else floor(min_fraction_snps_trimmed*100)-1]/1.5

  # min_fraction_hmer_indels
  Array[Float] min_fraction_non_hmer_values = [5244, 2829, 1446, 1074, 798, 624, 483, 392, 314, 251, 217, 182]
  Float min_fraction_non_hmer_value = min_fraction_non_hmer_values[
        if (min_fraction_non_hmer_indels_trimmed*100)-1 < 0
        then 0 else floor(min_fraction_non_hmer_indels_trimmed*100)-1]/1.5

  Float expected_genome_wide_output_size = min_fraction_hmer_indels_value +
                                            min_fraction_snps_value +
                                            min_fraction_non_hmer_value
  Float shard_region_length = overall_calling_regions_length / total_number_of_shards
  Float shard_region_fraction_of_genome = shard_region_length / genome_length
  Int realigned_bam_size = if output_realignment then
                           7 * ceil(size(cram_files, "GB") / total_number_of_shards) else 0
  Int expected_output_size = ceil(expected_genome_wide_output_size * shard_region_fraction_of_genome + realigned_bam_size)
  Int inputs_size =  ceil(size(cram_files, "GB") / total_number_of_shards + size(references.ref_fasta, "GB"))

  Float c_i = 1.3  # inputs safety factor
  Float c_o = 1.3 # outputs safety factor
  Int disk_size = ceil(c_i * inputs_size + c_o * expected_output_size)

  String output_prefix = basename(interval, ".interval_list")
  String gvcf_string = if make_gvcf then ("--gvcf  --p-error " + p_error) else ""
  Boolean defined_background = length(background_cram_files) > 0


  parameter_meta {
      cram_files: {
          localization_optional: true
      }
      background_cram_files: {
          localization_optional: true
      }
  }

  command <<<
    set -xeo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &

    if [[ "~{cloud_provider}" != "aws" ]]; then
      gatk --java-options "-Xms2G" PrintReads \
          -I ~{sep=' -I ' cram_files} \
          -O /dev/stdout \
          -L ~{interval} \
          -R ~{references.ref_fasta} | \
      samtools view -C --reference ~{references.ref_fasta} -o input.cram
      samtools index input.cram
      input=input.cram
      input_index=input.cram.crai

      background=''
      background_index=''
      if [[ "~{defined_background}" == "true" ]]; then
        gatk --java-options "-Xms2G" PrintReads \
            -I  ~{sep=' -I ' background_cram_files} \
            -O /dev/stdout \
            -L ~{interval} \
            -R ~{references.ref_fasta} | \
        samtools view -C --reference ~{references.ref_fasta} -o background.cram
        samtools index background.cram
        background=background.cram
        background_index=background.cram.crai
      fi

    else
      input=~{sep=',' cram_files}
      input_index=~{sep=',' cram_index_files}
      background=~{sep=',' background_cram_files}
      background_index=~{sep=',' background_cram_index_files}
    fi

    echo 'Input files are:'
    echo $input
    echo $input_index
    echo $background
    echo $background_index

    input_string="~{true='$input;$background' false='$input' defined_background}"
    input_index_string="~{true='$input_index;$background_index' false='$input_index' defined_background}"

    echo 'Input strings are:'
    echo $input_string
    echo $input_index_string

    # Create empty files such that outputs will still be valid even if cloud is aws
    touch input.cram
    touch input.cram.crai
    touch background.cram
    touch background.cram.crai

    cat ~{interval} | grep -v @ | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3}' >> interval.bed

    echo 'interval.bed content:'
    cat interval.bed

    tool \
      --input "$input_string" \
      --cram-index "$input_index_string" \
      --output ~{output_prefix} \
      --reference ~{references.ref_fasta} \
      --bed interval.bed \
      --min-base-quality ~{min_base_quality} \
      --min-mapq ~{pileup_min_mapping_quality} \
      --cgp-min-count-snps ~{min_read_count_snps} \
      --cgp-min-count-hmer-indels ~{min_read_count_hmer_indels} \
      --cgp-min-count-non-hmer-indels ~{min_read_count_non_hmer_indels} \
      --cgp-min-fraction-snps ~{min_fraction_snps} \
      --cgp-min-fraction-hmer-indels ~{min_fraction_hmer_indels} \
      --cgp-min-fraction-non-hmer-indels ~{min_fraction_non_hmer_indels} \
      --cgp-min-mapping-quality ~{candidate_min_mapping_quality} \
      --max-reads-per-region ~{max_reads_per_partition} \
      --assembly-min-base-quality ~{assembly_min_base_quality} \
      --gzip-output \
      ~{true="--no-realigned-sam" false="" !output_realignment} \
      --interval-nreads ~{interval_nreads} \
      ~{true="--somatic" false="" defined_background} \
      ~{gvcf_string} \
      --optimal-coverages "~{sep=";" optimal_coverages}" \
      ~{true="--cap-at-optimal-coverage " false="" cap_at_optimal_coverage} \
      ~{true="--prioritize-alt-supporting-reads " false="" prioritize_alt_supporting_reads} \
      --cycle-examples-min 100000 \
      ~{ug_make_examples_extra_args}

    #TODO - play with --interval-nreads, in order to control memory
    #--progress for debug

    if [ ~{output_realignment} == "true" ]
    then
      ls -lh ~{output_prefix}_hap_out.sam
      samtools sort ~{output_prefix}_hap_out.sam | samtools view -C --reference ~{references.ref_fasta} -o ~{output_prefix}_realign.cram
      samtools index ~{output_prefix}_realign.cram
      ls -lh ~{output_prefix}_realign.cram
    fi
    touch "~{output_prefix}_realign.cram"
    touch "~{output_prefix}_realign.cram.crai"

    touch "~{output_prefix}.gvcf.tfrecord.gz"

  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    preemptible: preemptible_tries
  }
  output {
    File monitoring_log = "monitoring.log"
    File input_cram = "input.cram"
    File input_cram_index = "input.cram.crai"
    File background_cram  = "background.cram"
    File background_cram_index = "background.cram.crai"
    File realigned_cram = "~{output_prefix}_realign.cram"
    File realigned_cram_index = "~{output_prefix}_realign.cram.crai"
    Array[File] output_examples = glob("~{output_prefix}.tfrecord*.gz")
    Array[File] gvcf_records = glob("~{output_prefix}.gvcf.tfrecord*.gz")
    Array[File] output_jsons = glob("~{output_prefix}*.json")
    Array[File] output_gatk = glob("~{output_prefix}*.gatk")
  }
}


task UGCallVariants{
  input{
    Array[File] examples
    File model_onnx
    File? model_serialized
    String docker
    Int call_variants_uncompr_buf_size_gb
    String gpu_type
    Int num_gpus
    Int num_cpus
    Int num_threads
    File monitoring_script
    Int disk_size = ceil(size(examples, 'GB') + 10)
    Int? call_variants_extra_mem
  }
  Int num_examples = length(examples)
  Int extra_mem = select_first([call_variants_extra_mem, 8])
  Int mem = num_threads * call_variants_uncompr_buf_size_gb + extra_mem
  command <<<
    set -eo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &

    nvidia-smi --query-gpu=timestamp,name,driver_version,temperature.gpu,utilization.gpu,utilization.memory,memory.used --format=csv -l 60 -f  ./nvidia-smi.log & #todo change from old

    # Rename the serialized model to fit the onnx file (allow various versions of serialized model files with different file names)
    if [ "~{defined(model_serialized)}" == "true" ] && [ "~{model_serialized}" != "~{model_onnx}.serialized" ]
    then
      echo 'Renaming serialized model to fit onnx file'
      mv ~{model_serialized} ~{model_onnx}.serialized
    fi

    printf "%b\n" "[RT classification]" \
      "onnxFileName = ~{model_onnx}" \
      "useSerializedModel = 1" \
      "trtWorkspaceSizeMB = 2000" \
      "numInferTreadsPerGpu = 2" \
      "useGPUs = ~{num_gpus}" \
      "gpuid = 0\n" \
      "[debug]" \
      "logFileFolder = .\n" \
      "[general]" \
      "tfrecord = 1" \
      "compressed = 1" \
      "outputInOneFile = 0" \
      "numUncomprThreads = ~{num_threads}" \
      "uncomprBufSizeGB = ~{call_variants_uncompr_buf_size_gb}" \
      "outputFileName = call_variants" \
      "numConversionThreads = 2" \
      "numExampleFiles = ~{num_examples}\n" > params.ini

    cat ~{write_lines(examples)} | awk '{print "exampleFile" NR " = "$0}' >> params.ini

    call_variants --param params.ini --fp16

  >>>
  runtime {
    memory: "~{mem} GB"
    cpu: "~{num_cpus}"
    disks: "local-disk " + disk_size + " LOCAL"
    docker: docker
    gpuType: "nvidia-tesla-p100"
    gpuCount: num_gpus
    acceleratorType : gpu_type
    acceleratorCount : num_gpus
  }
output {
    File monitoring_log = "monitoring.log"
    File nvidia_smi_log = "nvidia-smi.log"
    File params = "params.ini"
    Array[File] log = glob('call_variants*.log')
    Array[File] output_records = glob('call_variants*.gz')
    # File output_model_serialized = "~{model_onnx}.serialized" # uncomment to output the serilized model (need also to uncomment output_model_serialized in efficient_dv.wdl outputs)
  }
}

task UGPostProcessing{
  input{
    Array[File] called_records
    File ref
    File ref_index
    String docker
    String output_prefix
    File exome_intervals
    Array[File]? annotation_intervals
    File dbsnp
    File dbsnp_index
    String flow_order
    Int qual_filter
    Array[File]? gvcf_records
    Boolean make_gvcf
    Int min_variant_quality_hmer_indels
    Int min_variant_quality_exome_hmer_indels
    Int min_variant_quality_non_hmer_indels
    Int min_variant_quality_snps
    Boolean show_bg_fields

    File monitoring_script

    Int disk_size = ceil(48 * size(called_records, "GB") +
                         size(ref, "GB") +
                         size(dbsnp, "GB") +
                         (if make_gvcf then size(select_first([gvcf_records]), "GB") else 0) +
                         4 + (if make_gvcf then 5 else 0))

  }
  String gvcf_args = if make_gvcf then "--gvcf_outfile ~{output_prefix}.g.vcf.gz --nonvariant_site_tfrecord_path @gvcf_records.txt" else ""
  Array[File] gvcf_records_not_opt = select_first([gvcf_records, []])
  Array[File] empty_array_of_files = []
  Array[File] annotation_intervals_or_empty = select_first([annotation_intervals, empty_array_of_files])
  Array[File] exome_and_annotations = flatten([[exome_intervals], annotation_intervals_or_empty])
  command <<<
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      set -eo pipefail

      printf "%b\n" "LowQualInExome" \
        "QUAL < ~{min_variant_quality_exome_hmer_indels} and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and vc.hasAttribute('EXOME')" \
        "LowQual" \
        "QUAL < ~{min_variant_quality_hmer_indels} and VARIANT_TYPE=='h-indel' and not vc.isFiltered() and not vc.hasAttribute('EXOME')" \
        "LowQual" \
        "QUAL < ~{min_variant_quality_non_hmer_indels} and VARIANT_TYPE=='non-h-indel' and not vc.isFiltered()" \
        "LowQual" \
        "QUAL < ~{min_variant_quality_snps} and VARIANT_TYPE=='snp' and not vc.isFiltered()" \
        "LargeDeletion" \
        "REFLEN > 220 and vc.isFiltered()" \
        > filters.txt

      if [ ~{make_gvcf} == "true" ]
      then
        cp ~{write_lines(gvcf_records_not_opt)} gvcf_records.txt
      fi

      ug_postproc \
        --infile ~{sep="," called_records} \
        --ref ~{ref} \
        --outfile ~{output_prefix}.vcf.gz \
        ~{gvcf_args} \
        --consider_strand_bias \
        --flow_order ~{flow_order} \
        --annotate \
        --bed_annotation_files ~{sep="," exome_and_annotations} \
        --qual_filter ~{qual_filter} \
        --filter \
        --filters_file filters.txt \
        --dbsnp ~{dbsnp} \
        ~{if show_bg_fields then "--consider_bg_fields" else ""}

      bcftools index -t ~{output_prefix}.vcf.gz
      if [ ~{make_gvcf} == "true" ]
      then
        bcftools index -t ~{output_prefix}.g.vcf.gz
      fi

      touch ~{output_prefix}.g.vcf.gz
      touch ~{output_prefix}.g.vcf.gz.tbi

      bcftools view -f PASS -O z ~{output_prefix}.vcf.gz -o ~{output_prefix}.pass.vcf.gz
      bcftools index -t ~{output_prefix}.pass.vcf.gz

      # Run QC reprort
      echo 'Running QC report...'
      python /opt/deepvariant/qc_report/run_no_gt_report.py \
        --input_file ~{output_prefix}.pass.vcf.gz \
        --reference ~{ref} \
        --output_prefix ~{output_prefix} \
        --output_metrics_h5

  >>>
  runtime {
    memory: "8 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }
  output {
    File monitoring_log = "monitoring.log"
    File vcf_file  = '~{output_prefix}.vcf.gz'
    File vcf_index = '~{output_prefix}.vcf.gz.tbi'
    File gvcf_file = '~{output_prefix}.g.vcf.gz'
    File gvcf_file_index = '~{output_prefix}.g.vcf.gz.tbi'
    File qc_h5     = '~{output_prefix}.h5'
    File qc_report = '~{output_prefix}_report.html'
    File qc_metrics_h5 = '~{output_prefix}_metrics.h5'
  }
}
