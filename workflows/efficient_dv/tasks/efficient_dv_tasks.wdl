version 1.0
import "structs.wdl" as Structs
task UGMakeExamples{
  input{
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

    File? germline_vcf

    Int min_base_quality
    Int pileup_min_mapping_quality
    Int min_read_count_snps
    Int min_read_count_hmer_indels
    Int min_read_count_non_hmer_indels
    Float min_fraction_snps
    Float min_fraction_hmer_indels
    Float min_fraction_non_hmer_indels
    Float min_fraction_single_strand_non_snps
    Int? min_hmer_plus_one_candidate
    Int candidate_min_mapping_quality
    Int max_reads_per_partition
    Int assembly_min_base_quality
    Boolean make_gvcf
    Float p_error = 0
    Boolean prioritize_alt_supporting_reads    
    Array[Int] optimal_coverages
    Boolean cap_at_optimal_coverage
    Boolean is_somatic
    Boolean output_realignment = false
    Boolean single_strand_filter = false
    Boolean keep_duplicates = true
    Boolean add_ins_size_channel = true
    String? extra_args
    Boolean log_progress = false
    Boolean count_candidates_with_dvtools = false
    Boolean normalize_strand_bias = false
    Array[Float]? strand_bias_normalization_thresholds

    String docker
    File monitoring_script
    Float memory
    Int cpu
    Int preemptible_tries
    String cloud_provider
    Boolean no_address = true
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
  Boolean defined_background = length(background_cram_files) > 0

  Array[Float] strand_bias_normalization_thresholds_defined = select_first([strand_bias_normalization_thresholds, [0.5,0.5]])

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
          -R ~{references.ref_fasta} |\
      samtools view -C -T ~{references.ref_fasta} -o input.cram --output-fmt-option embed_ref=1 -

      samtools index input.cram -@ ~{cpu}
      input=input.cram
      input_index=input.cram.crai

      background=''
      background_index=''
      if [[ "~{defined_background}" == "true" ]]; then
        gatk --java-options "-Xms2G" PrintReads \
            -I  ~{sep=' -I ' background_cram_files} \
            -O /dev/stdout \
            -L ~{interval} \
            -R ~{references.ref_fasta} |\
        samtools view -C -T ~{references.ref_fasta} -o background.cram --output-fmt-option embed_ref=1 -

        samtools index background.cram -@ ~{cpu}
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

    # Process interval file
    grep -v @ ~{interval} | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3}' >> interval.bed

    # split the interval file into parts for parallel processing (number of parts as number of cpus)
    echo 'splitting interval into ~{cpu} parts'
    total_lines=$(wc -l < interval.bed)
    lines_per_file=$((total_lines / ~{cpu}))
    remainder=$((total_lines % ~{cpu}))
    if [ $remainder -ne 0 ]; then
      lines_per_file=$((lines_per_file + 1))
    fi

    split -l $lines_per_file -d --additional-suffix .bed interval.bed interval_ #should create files: interval_00.bed, interval_01.bed, etc.

    echo 'Splited interval.bed content:'
    for file in interval_*.bed; do
      echo "$file content:"
      cat "$file"
      echo -n "interval size is: "
      awk '{sum += $3 - $2} END {print sum}' "$file"
    done
    
    #shellcheck disable=SC2034
    st_bias_norm_thresholds_string="--strand-bias-threshold-to-normalize ~{sep=',' strand_bias_normalization_thresholds_defined}"

    # parallel processing: run the tool in different process, each on a different interval part
    pids=()
    for interval_part in $(find . -name "interval_*.bed" | sort); do
      part_number=$(echo "$interval_part" | grep -o -E '[0-9]+') #extract the part number from the file name
      tool \
        --input "$input_string" \
        --cram-index "$input_index_string" \
        --output "~{output_prefix}_$part_number" \
        --reference ~{references.ref_fasta} \
        --bed "$interval_part" \
        --min-base-quality ~{min_base_quality} \
        --min-mapq ~{pileup_min_mapping_quality} \
        --cgp-min-count-snps ~{min_read_count_snps} \
        --cgp-min-count-hmer-indels ~{min_read_count_hmer_indels} \
        --cgp-min-count-non-hmer-indels ~{min_read_count_non_hmer_indels} \
        --cgp-min-fraction-snps ~{min_fraction_snps} \
        --cgp-min-fraction-hmer-indels ~{min_fraction_hmer_indels} \
        --cgp-min-fraction-non-hmer-indels ~{min_fraction_non_hmer_indels} \
        --cgp-min-fraction-single-strand-non-snps ~{min_fraction_single_strand_non_snps} \
        --cgp-min-mapping-quality ~{candidate_min_mapping_quality} \
        --cgp-min-hmer-plus-one-candidate ~{if defined(min_hmer_plus_one_candidate) then min_hmer_plus_one_candidate else 7} \
        --max-reads-per-region ~{max_reads_per_partition} \
        --assembly-min-base-quality ~{assembly_min_base_quality} \
        ~{true="--realigned-sam" false="" output_realignment} \
        ~{true="--somatic" false="" is_somatic} \
        ~{if make_gvcf then "--gvcf  --p-error ~{p_error}" else ""} \
        --optimal-coverages "~{sep=";" optimal_coverages}" \
        ~{true="--cap-at-optimal-coverage " false="" cap_at_optimal_coverage} \
        ~{true="--prioritize-alt-supporting-reads " false="" prioritize_alt_supporting_reads} \
        ~{true="--normalize-strand-bias" false="" normalize_strand_bias} \
        ~{true="$st_bias_norm_thresholds_string" false="" normalize_strand_bias } \
        --cycle-examples-min 100000 \
        --prefix-logging-with "${part_number}>> " \
        ~{true="--single-strand-filter" false="" single_strand_filter} \
        ~{true="--keep-duplicates" false="" keep_duplicates} \
        ~{true="--add-ins-size-channel" false="" add_ins_size_channel} \
        ~{extra_args} \
        ~{true="--progress" false="" log_progress} \
        ~{if defined(germline_vcf) then "--region-haplotypes-vcf ~{germline_vcf}" else ""} \
         &
        
      # Save the PID of the process
      pids+=($!)
    done

  # Wait for the process running the tool to finish (don't wait for the process running the monitor log)
  # if one process is failed, kill all the other processes and exit with error
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      echo "ERROR occurred. View error using: *** for prefix in \"00>>\" \"01>>\"; do cat log | grep \"^\$prefix\" | tail; done; *** Killing all other processes and exiting."
      for other_pid in "${pids[@]}"; do
        if [ "$other_pid" != "$pid" ]; then
          kill "$other_pid"
        fi
      done
      exit 1
    fi
  done

  if [ ~{output_realignment} == "true" ]
  then
      # concat the output ~{output_prefix}_hap_out.sam files into one (using samtools)
      samtools cat -@ ~{cpu} -o ~{output_prefix}_hap_out.sam ~{output_prefix}_*_hap_out.sam | \
      samtools sort -@ ~{cpu} ~{output_prefix}_hap_out.sam | \
      samtools view -C -@ ~{cpu} --reference ~{references.ref_fasta} -o ~{output_prefix}_realign.cram
      samtools index ~{output_prefix}_realign.cram -@ ~{cpu}
      ls -lh ~{output_prefix}_realign.cram
  fi
  touch "~{output_prefix}_realign.cram"
  touch "~{output_prefix}_realign.cram.crai"

  ls -lh -- *tfrecord*

  if [ ~{count_candidates_with_dvtools} == "true" ]
  then
    for tfrec in $(find . -name '*tfrecord*' | sort); do
      dvtools --infile "$tfrec" --filetype dv --op vcf --outfile "$tfrec.vcf"
      echo "$tfrec: $(wc -l < "$tfrec.vcf")"
      rm "$tfrec.vcf"
    done
  fi

  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    preemptible: preemptible_tries
    noAddress: no_address
  }
  output {
    File monitoring_log = "monitoring.log"
    File input_cram = "input.cram"
    File input_cram_index = "input.cram.crai"
    File background_cram  = "background.cram"
    File background_cram_index = "background.cram.crai"
    File realigned_cram = "~{output_prefix}_realign.cram"
    File realigned_cram_index = "~{output_prefix}_realign.cram.crai"
    Array[File] output_examples = glob("~{output_prefix}*[0-9].tfrecord.gz")
    Array[File?] gvcf_records = glob("~{output_prefix}*.gvcf.tfrecord.gz")
    Array[File] output_jsons = glob("~{output_prefix}*.json")
    Array[File] output_gatk = glob("~{output_prefix}*.gatk")
  }
}


task UGCallVariants{
  input{
    Array[File] examples
    File model_onnx
    File? model_serialized
    Boolean is_somatic
    String docker
    Int call_variants_uncompr_buf_size_gb
    String gpu_type
    Int num_gpus
    Int num_cpus
    Int num_threads
    File monitoring_script
    Int disk_size = ceil(1.05*size(examples, 'GB') + 10)
    Int? call_variants_extra_mem
    Int? optimization_level
    Boolean no_address = true
  }
  Int num_examples = length(examples)
  Int extra_mem = select_first([call_variants_extra_mem, 8])
  Int builder_optimization_level  = select_first([optimization_level, if is_somatic then 5 else 1])
  Int mem = num_threads * call_variants_uncompr_buf_size_gb + extra_mem
  String onnx_base_name = basename(model_onnx)
  command <<<
    set -eo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &

    nvidia-smi --query-gpu=timestamp,name,driver_version,temperature.gpu,utilization.gpu,utilization.memory,memory.used --format=csv -l 60 -f  ./nvidia-smi.log & #todo change from old

    cp ~{model_onnx} ~{onnx_base_name}
    # Rename the serialized model to fit the onnx file (allow various versions of serialized model files with different file names)        
    if [ "~{defined(model_serialized)}" == "true" ] && [ "~{model_serialized}" != "~{model_onnx}.serialized" ]
    then
      echo 'Renaming serialized model to fit onnx file'
      cp ~{model_serialized} ~{onnx_base_name}.serialized
    fi

    printf "%b\n" "[RT classification]" \
      "onnxFileName = ~{onnx_base_name}" \
      "builderOptimizationLevel = ~{builder_optimization_level}" \
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

    awk '{print "exampleFile" NR " = "$0}' ~{write_lines(examples)} >> params.ini

    call_variants --param params.ini --fp16

    num_candidates_val=$(grep -oP 'total batch size \K\d+(?= vectors)' call_variants*.log)
    echo "$num_candidates_val" > "num_candidates_${num_candidates_val}"
    echo "$num_candidates_val" > nc.txt

  >>>
  runtime {
    memory: "~{mem} GB"
    cpu: "~{num_cpus}"
    disks: "local-disk " + disk_size + " LOCAL"
    docker: docker
    gpuType: "nvidia-tesla-p100"
    gpuCount: num_gpus
    acceleratorType : gpu_type #!UnknownRuntimeKey
    acceleratorCount : num_gpus #!UnknownRuntimeKey
    noAddress: no_address
  }
output {
    File monitoring_log = "monitoring.log"
    File nvidia_smi_log = "nvidia-smi.log"
    File params = "params.ini"
    Array[File] log = glob('call_variants*.log')
    Array[File] output_records = glob('call_variants*.gz')
    Array[File] num_candidates = glob("num_candidates_*")
    Int num_candidates_as_int = read_int("nc.txt")
    File output_model_serialized = "~{onnx_base_name}.serialized"
  }
}

task UGPostProcessing{
  input{
    Array[File] called_records
    Array[File] cram_files
    Array[File] cram_index_files

    # Background sample inputs
    Array[File] background_cram_files
    Array[File] background_cram_index_files

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
    Boolean recalibrate_vaf
    Boolean is_somatic
    Int min_variant_quality_hmer_indels
    Int min_variant_quality_exome_hmer_indels
    Int min_variant_quality_non_hmer_indels
    Int min_variant_quality_snps
    Boolean show_bg_fields
    String extra_args = ""

    File monitoring_script

    Int disk_size = ceil(48 * size(called_records, "GB") +
                         size(ref, "GB") +
                         size(dbsnp, "GB") +
                         (if make_gvcf then size(select_first([gvcf_records]), "GB") else 0) +
                         4 + (if make_gvcf then 5 else 0)) + 
                         ceil(size(background_cram_files, "GB")) + 
                         ceil(size(cram_files, "GB")) +
                         ceil(size(cram_index_files, "GB")) +
                         ceil(size(background_cram_index_files, "GB"))
    Boolean no_address = true
  }

  Int indel_threshold_for_recalibration = 30000

  String gvcf_args = if make_gvcf then "--gvcf_outfile ~{output_prefix}.g.vcf.gz --nonvariant_site_tfrecord_path @gvcf_records.txt --hcr_bed_file ~{output_prefix}.hcr.bed " else ""
  Array[File] gvcf_records_not_opt = select_first([gvcf_records, []])
  Array[File] empty_array_of_files = []
  Array[File] annotation_intervals_or_empty = select_first([annotation_intervals, empty_array_of_files])
  Array[File] exome_and_annotations = flatten([[exome_intervals], annotation_intervals_or_empty])
  Boolean defined_background = length(background_cram_files) > 0

  command <<<
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      set -xeo pipefail

      cp ~{write_lines(called_records)} called_records.txt

      echo 'Defining filters...'
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

      # shellcheck disable=SC2034
      foreground=~{sep=',' cram_files} 
      # shellcheck disable=SC2034
      background=~{sep=',' background_cram_files} 
      cram_string="~{true='$foreground;$background' false='$foreground' defined_background}"

      echo "Calculating approximate INDEL variant count"
      ug_postproc \
          --infile @called_records.txt \
          --ref ~{ref} \
          --outfile "~{output_prefix}.vcf.gz" \
          --qual_filter ~{qual_filter} \
          --count_indels --group_variants false |& tee indel.count.log

      indel_count=$( grep indel_count indel.count.log | cut -d " " -f 4 )
      echo "Approximate INDEL variant count: $indel_count"
      if [ "$indel_count" -gt ~{indel_threshold_for_recalibration} ]
      then
        echo "INDEL variant count is too high, skipping post-processing"
        recalibration_string=""
      else 
        echo "INDEL variant count is low, recalibrating VAF"
        # shellcheck disable=SC2034
        recalibration_string="--fix_allele_coverage --fix_allele_indels_only --fix_allele_crams $cram_string"
      fi
     

      echo 'Running UG post-processing...'
      ug_postproc \
        --infile @called_records.txt \
        --ref ~{ref} \
        --outfile "~{output_prefix}.vcf.gz" \
        ~{gvcf_args} \
        --consider_strand_bias \
        --flow_order ~{flow_order} \
        --annotate \
        --bed_annotation_files ~{sep="," exome_and_annotations} \
        --qual_filter ~{qual_filter} \
        --filter \
        --filters_file filters.txt \
        --dbsnp ~{dbsnp} \
        ~{if show_bg_fields then "--consider_bg_fields" else ""} \
        ~{if recalibrate_vaf then '$recalibration_string' else ""} \
        ~{if is_somatic then "--ignore_multi_allelic_cvos" else ""} \
        ~{extra_args}

      bcftools index -t "~{output_prefix}.vcf.gz"
      if [ ~{make_gvcf} == "true" ]
      then
        bcftools index -t "~{output_prefix}.g.vcf.gz"
      fi

      touch "~{output_prefix}.g.vcf.gz"
      touch "~{output_prefix}.g.vcf.gz.tbi"
      touch "~{output_prefix}.hcr.bed"

    echo 'Saving header IDs to a file...'
    export header_file=header.hdr
    export id_file=header.ID
    export all_ids=all_ids.txt

    for f in ~{sep=" " exome_and_annotations}
    do
      head -1 $f > $header_file
      sed 's/[<>]/ /' "$header_file" | sed 's/[=]/ /' | sed 's/[=]/ /' | sed 's/[,]/ /' | awk '{print $3}' > $id_file
      cat $id_file >> $all_ids
    done

  >>>
  runtime {
    memory: "8 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File monitoring_log = "monitoring.log"
    File vcf_file  = '~{output_prefix}.vcf.gz'
    File vcf_index = '~{output_prefix}.vcf.gz.tbi'
    File gvcf_file = '~{output_prefix}.g.vcf.gz'
    File gvcf_file_index = '~{output_prefix}.g.vcf.gz.tbi'
    File gvcf_hcr = '~{output_prefix}.hcr.bed'
    Array[String] interval_annotation_names = read_lines('all_ids.txt')
  }
}


task QCReport{
  input{
    File input_vcf
    File input_vcf_index
    File? callable_bed
    String output_prefix
    File ref
    String docker
    File monitoring_script

    Int disk_size = ceil(2 * size(input_vcf, "GB") +
                         size(ref, "GB") + 4)
    Boolean no_address = true
  }

  command <<<
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      set -eo pipefail

      echo 'Filtering PASS variants...'
      bcftools view -f PASS -O z ~{input_vcf} -o ~{output_prefix}.pass.vcf.gz
      bcftools index -t ~{output_prefix}.pass.vcf.gz

      echo 'Running QC reprort...'
      python /opt/deepvariant/qc_report/run_no_gt_report.py \
        --input_file "~{output_prefix}.pass.vcf.gz" \
        --reference ~{ref} \
        ~{"--hcr_bed " + callable_bed} \
        --output_prefix "~{output_prefix}" \
        --output_metrics_h5

  >>>
  runtime {
    memory: "8 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: no_address
  }
  output {
    File monitoring_log = "monitoring.log"
    File qc_h5     = '~{output_prefix}.h5'
    File qc_report = '~{output_prefix}_report.html'
    File qc_metrics_h5 = '~{output_prefix}_metrics.h5'
  }
}
