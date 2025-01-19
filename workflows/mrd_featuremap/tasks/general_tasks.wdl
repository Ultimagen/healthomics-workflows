version 1.0
# LICENSE
#   Copyright 2022 Ultima Genomics
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# DESCRIPTION
#   Tasks used in multiple workflows
# CHANGELOG
#

import "structs.wdl" as Structs

task ExtractSampleNameFlowOrder{
    input{
        File input_bam
        File monitoring_script
        Int preemptible_tries
        String docker
        References references
        Boolean no_address
        String? cloud_provider_override = "gcp"
    }
    String cloud_provider = select_first([cloud_provider_override, 'gcp'])

    parameter_meta {
        input_bam: {
            localization_optional: true
        }
    }

    command <<<
        set -e
        set -o pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &


        gatk GetSampleName  \
            -I ~{input_bam} \
            -R ~{references.ref_fasta} \
            -O sample_name.txt

        if [[ "~{cloud_provider}" == "aws" ]]
        then
            gatk ViewSam -I ~{input_bam} \
                    --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                    | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/FO:/) {print substr($i,4,4)}}}' | sed '1!d' \
                    > flow_order.txt

            gatk ViewSam -I ~{input_bam} \
                    --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                    | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/BC:/) {print substr($i,4)}}}' | sed '1!d' \
                    > barcode.txt

            gatk ViewSam -I ~{input_bam} \
                    --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                    | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/ID:/) {print substr($i,4)}}}' | sed '1!d' \
                    > id.txt
        else
            gsutil cat ~{input_bam} | gatk ViewSam -I /dev/stdin \
                --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/FO:/) {print substr($i,4,4)}}}' | sed '1!d' \
                > flow_order.txt

            gsutil cat ~{input_bam} | gatk ViewSam -I /dev/stdin \
                    --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                    | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/BC:/) {print substr($i,4)}}}' | sed '1!d' \
                    > barcode.txt

            gsutil cat ~{input_bam} | gatk ViewSam -I /dev/stdin \
                    --HEADER_ONLY true --ALIGNMENT_STATUS All --PF_STATUS All \
                    | grep "^@RG" | awk '{for (i=1;i<=NF;i++){if ($i ~/ID:/) {print substr($i,4)}}}' | sed '1!d' \
                    > id.txt
        fi
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        preemptible: preemptible_tries
        noAddress: no_address
        docker: docker
    }

    output {
        String sample_name = read_string("sample_name.txt")
        String flow_order = read_string("flow_order.txt")
        String barcode_seq = read_string("barcode.txt")
        String readgroup_id = read_string("id.txt")
        File sample_name_file = "sample_name.txt"
        File flow_order_file = "flow_order.txt"
        File barcode_seq_file = "barcode.txt"
        File readgroup_id_file = "id.txt"
        File monitoring_log = "monitoring.log"
    }
}

# only works with simple maps (two columns)
task MakeStringMap {
    input {
        Array[String] keys
        Array[String] values
        String docker
    }

    String results_path = "results.tsv"
    command <<<
        cat ~{write_tsv(transpose([keys,values]))} > ~{results_path}
    >>>
    runtime {
        docker: docker
    }
    output {
        Map[String, String] map = read_map(results_path)
    }
}


# Aggregate picard metrics and coverage and variant calling reports into one h5 and convert h5 to json
task AggregateMetricsAndConvertToJson {
    input {
        File monitoring_script
        Array[File]? picard_files
        String base_file_name
        File? coverage_all_h5
        File? short_report_h5
        File? extended_report_h5
        File? no_gt_report_h5
        Float? contamination_stdout
        String docker
        Int disk_size = 20
        Int preemptible_tries
        Boolean no_address
    }

    command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    collect_existing_metrics \
        ~{true="--metric_files " false="" defined(picard_files)} ~{ sep=' ' picard_files} \
        ~{"--coverage_h5 " + coverage_all_h5} \
        ~{"--short_report_h5 " + short_report_h5} \
        ~{"--extended_report_h5 " + extended_report_h5} \
        ~{"--no_gt_report_h5 " + no_gt_report_h5} \
        ~{"--contamination_stdout " + contamination_stdout} \
        --output_h5 ~{base_file_name}.aggregated_metrics.h5
    convert_h5_to_json \
            --root_element "metrics" \
            --ignored_h5_key_substring histogram \
            --input_h5 ~{base_file_name}.aggregated_metrics.h5 \
            --output_json ~{base_file_name}.aggregated_metrics.json
    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_h5 = "~{base_file_name}.aggregated_metrics.h5"
        File aggregated_metrics_json = "~{base_file_name}.aggregated_metrics.json"
    }
}

task ConvertSorterStatsToH5 {
    input {
        File monitoring_script
        Int? file_name_suffix

        # runtime arguments
        Int disk_size = 20
        Int preemptible_tries
        String docker
        Boolean no_address

        # sorter_to_h5 args
        File input_csv_file
        File input_json_file
    }
    
    String input_basename = basename(input_csv_file, ".csv")
    String suffix = if defined(file_name_suffix) then ".~{file_name_suffix}" else ""
    String aggregated_metrics_h5_file = "~{input_basename}~{suffix}.aggregated_metrics.h5"
    String aggregated_metrics_json_file = "~{input_basename}~{suffix}.aggregated_metrics.json"
    
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        export HDF5_USE_FILE_LOCKING=FALSE #to prevent unable to lock file, errno = 11, error message = 'Resource temporarily unavailable' errors

        sorter_to_h5 \
            --input_csv_file "~{input_csv_file}" \
            --input_json_file "~{input_json_file}" \
            --output_h5_file "~{aggregated_metrics_h5_file}"

        convert_h5_to_json \
            --root_element "metrics" \
            --ignored_h5_key_substring histogram \
            --input_h5 "~{aggregated_metrics_h5_file}" \
            --output_json "~{aggregated_metrics_json_file}"

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_h5 = "~{aggregated_metrics_h5_file}"
        File aggregated_metrics_json = "~{aggregated_metrics_json_file}"
    }
}

task ScatterIntervalList {
    input {
        File monitoring_script
        File interval_list
        Int scatter_count
        Int break_bands_at_multiples_of
        String docker
        String gitc_path
        Boolean no_address
        String dummy_input_for_call_caching  # !UnusedDeclaration
    }
    command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    echo ~{dummy_input_for_call_caching}
    set -eo pipefail
    mkdir out
    java -Xms4g -jar ~{gitc_path}picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i+1).rjust(len(str(len(intervals))),"0") + filename)
      os.rename(interval, newName)
    with open("interval_count.txt", "w") as fh:
        fh.write(str(len(intervals)))
    CODE
    >>>
    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int('interval_count.txt')
        File monitoring_log = "monitoring.log"
    }
    runtime {
        memory: "6 GB"
        docker: docker
        noAddress: no_address
    }
}

task IntervalListOfGenome {
  input  {
    File ref_fai
    File ref_dict
    String docker
    Int disk_size
    Int preemptible_tries
    Boolean no_address
    File monitoring_script
  }
  command <<<
    set -e
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    cat ~{ref_fai} | grep -v '[ME_-]' | awk '{print $1"\t"1"\t"$2"\t+\t."}' > modifai
    cat ~{ref_dict} modifai > genome.interval_list
  >>>
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    continueOnReturnCode: true
    maxRetries: 1
    noAddress: no_address
  }
  output{
    File monitoring_log = "monitoring.log"
    File interval_list = "genome.interval_list"
  }
}


task IntervalListFromString {
  input  {
    String intervals_string
    File ref_fai
    File ref_dict
    String docker
    Int disk_size
    Int preemptible_tries
    Boolean no_address
    File monitoring_script
  }
  command <<<
    set -e
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    echo "~{intervals_string}" | tr ':-' '\t\t' | sed 's/;/\n/g' | awk '{print $0 "\t+\t. intersection ACGTmer.1"}' > interval.tsv
    cat ~{ref_dict} interval.tsv > intervals.interval_list
  >>>
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    continueOnReturnCode: true
    maxRetries: 1
    noAddress: no_address
  }
  output{
    File monitoring_log = "monitoring.log"
    File interval_list = "intervals.interval_list"
  }
}


task IntervalListTotalLength {
  input {
    File interval_list
    String docker
    Boolean no_address
    File monitoring_script

  }
   command <<<
    set -e
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    grep -v @ ~{interval_list} | awk '{sum+=$3-$2}; END {printf "%0.f\n", sum}' > sum.txt
    cat sum.txt
  >>>

  output {
    Float interval_list_length = read_float('sum.txt')
    File monitoring_log = "monitoring.log"
  }

   runtime {
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + 4 + " HDD"
    docker: docker
    noAddress: no_address
  }
}

task FastaLengthFromIndex {
  input {
    File fasta_index
    String docker
    Boolean no_address
    File monitoring_script

  }
   command <<<
     set -e
     bash ~{monitoring_script} >&2 &
     awk '{sum+=$2} END {printf "%0.f\n", sum}' ~{fasta_index} > sum.txt
  >>>

  output {
    Float fasta_length = read_float('sum.txt')
  }

   runtime {
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + 4 + " HDD"
    docker: docker
    noAddress: no_address
  }
}

task ConcatMetricsJsons {
    input {
        File monitoring_script
        Array[File] jsons
        String base_file_name
        String docker
        Int preemptible_tries
        Boolean no_address
        Int disk_size = 20
    }

    command <<<
    set -eo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    echo ~{sep=',' jsons} > json_files.txt

    python <<CODE
    import json

    json_out = {'metrics': {}}
    with open('json_files.txt') as file_handler:
        json_files = file_handler.read().rstrip().split(',')
    for jf in json_files:
        with open(jf) as file_handler:
            json_content = json.load(file_handler)
        if "metrics" not in json_content:
            raise ValueError(f"Invalid json file {jf} - does not contain a 'metrics' key")
        joint_keys = set(json_content["metrics"]).intersection(set(json_out["metrics"]))
        if len(joint_keys) > 0:
            raise ValueError(f"Encountered duplicate metrics - {joint_keys}")
        json_out['metrics'].update(json_content['metrics'])
    with open('~{base_file_name}.aggregated_metrics.json','w') as file_handler:
        json.dump(json_out, file_handler, indent=2)
    CODE

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "1 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File aggregated_metrics_json = "~{base_file_name}.aggregated_metrics.json"
    }
}

task DownsampleCramBam {
    input {
        String base_file_name
        File input_cram_bam
        Float downsample_frac
        Int? seed 
        String output_format = "cram"
        References references
        File monitoring_script
        String docker
        Int preemptibles
        Int disk_size = ceil((1.1+(downsample_frac*4))*size(input_cram_bam,"GB") + 20)
        Int cpus = 16
        Int memory_gb = 1
    }
    String output_format_flag = if output_format == "cram" then "-C" else "-b"
    String output_format_extension = if output_format == "cram" then ".cram" else ".bam"
    String output_index_format_extension = if output_format == "cram" then ".crai" else ".bai"
    String output_cram_bam_name = base_file_name + output_format_extension
    
    command <<<
        set -eo pipefail
        set -o xtrace

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        if ~{defined(seed)}
        then
            seed_var="~{seed}"
        else
            seed_var=$RANDOM
        fi
        echo "Seed used " $seed_var
        echo $seed_var >> downsampling_seed.txt

        samtools view \
          -s ${seed_var}~{downsample_frac} \
          ~{output_format_flag} \
          -T ~{references.ref_fasta} \
          -o ~{output_cram_bam_name} \
          --threads ~{cpus} \
          ~{input_cram_bam}

        samtools index ~{output_cram_bam_name}
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: "~{cpus}"
        memory: "~{memory_gb} GB"
        preemptible: preemptibles
        docker: docker
        noAddress: false
    }
    output {
        File output_cram_bam = "~{output_cram_bam_name}"
        File output_cram_bam_index = "~{output_cram_bam_name}~{output_index_format_extension}"
        File downsampling_seed = "downsampling_seed.txt"
        File monitoring_log = "monitoring.log"
    }
}

task ToFastq {
    input {
        File input_cram
        String base_file_name
        References references
        File monitoring_script
        String docker
        Boolean no_address
        Int preemptibles
        Int disk_size = ceil(3*size(input_cram,"GB") + 20)
    }
    String output_fq_name = base_file_name + ".fq.gz"
    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        source ~/.bashrc
        conda activate genomics.py3

        samtools fastq \
          --reference ~{references.ref_fasta}\
          -F 0x900 \
          -0 ~{output_fq_name} \
          --threads 2 \
          ~{input_cram}
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 2
        memory: "1 GB"
        preemptible: preemptibles
        docker: docker
        noAddress: no_address
    }
    output {
        File output_fastq = "~{output_fq_name}"
        File monitoring_log = "monitoring.log"
    }
}

task ConcatHtmls {
   input {
        # input for concat
        Array[File] htmls

        # call arguments
        File monitoring_script
        String base_file_name

        # runtime arguments
        Int preemptible_tries
        String docker
        Int disk_size = ceil(2*size(htmls,"GB") + 1)



    }
    command <<<
        set -eo pipefail

        start=$(date +%s)
        source ~/.bashrc
        conda activate genomics.py3

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo ~{sep=',' htmls} > htmls_files.txt
        ls -lstr
        cat htmls_files.txt

        i=1
        i=1; for file in $(cat htmls_files.txt| sed "s/,/ /g"); do cp $file $i.html ;
        cat $i.html >> ~{base_file_name}_aggregated.html; rm $i.html ;i=$((i + 1)) ; done

        ls -lstr
        echo "**************** D O N E ****************"

        end=$(date +%s)
        mins_elapsed=$(( ($end - $start) / 60))
        secs_elapsed=$(( ($end - $start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"
        ls -ltrsa

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        docker: docker
        cpu: "1"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        noAddress: true
    }

    output {
        File report_html = "~{base_file_name}_aggregated.html"
        File monitoring_log = "monitoring.log"
    }
}

task ZipAndIndexVcf{
    input{
        File input_vcf
        File monitoring_script
        String docker
        Boolean no_address
        Int preemptibles
        Int disk_size = ceil(2*size(input_vcf,"GB")+20)
    }
    String output_vcf = basename(input_vcf) + ".gz"
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        source ~/.bashrc
        conda activate genomics.py3

        bcftools view -O z -o ~{output_vcf} ~{input_vcf}
        bcftools index -t ~{output_vcf}
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 4
        memory: "2 GB"
        preemptible: preemptibles
        docker: docker
        noAddress: no_address
    }
    output {
        File output_vcf_file = "~{output_vcf}"
        File output_vcf_index_file = "~{output_vcf}.tbi"
        File monitoring_log = "monitoring.log"
    }
}


task RenameSampleInBam {
    input{
        File monitoring_script
        File input_cram_bam
        File input_cram_bam_index
        String name
        String docker
        Int disk_size = ceil(2*size(input_cram_bam, "GB"))+20
        Int preemptible_tries
        Int cpus=1
        Boolean no_address
    }
    String detect_input_ending = if sub(input_cram_bam, ".*\\.cram$", "is_cram") == "is_cram" then "is_cram" else "is_bam"
    String extension = if detect_input_ending=="is_cram" then ".cram" else ".bam"
    String extension_index = if detect_input_ending=="is_cram" then ".crai" else ".bai"
    String basename = basename(input_cram_bam, extension)
    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        samtools view -H ~{input_cram_bam} -@ ~{cpus} | sed "s/SM:[^\t]*/SM:~{name}/g" | samtools reheader - ~{input_cram_bam} > ~{basename}~{extension}
        samtools index ~{basename}~{extension}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " LOCAL"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File monitoring_log = "monitoring.log"
        File output_bam = "~{basename}~{extension}"
        File output_bam_index = "~{basename}~{extension}~{extension_index}"
    }
}

task MergeCramFiles {
    input {
            File monitoring_script
            String output_base_name
            String sample_name
            Array[File] crams
            References references
            File? cache_tarball
            String docker
            Boolean no_address
            Int? cpus
            Int preemptible_tries
    }
    String output_cram_name = "~{output_base_name}.cram"
    Int cpus_to_use = select_first([cpus, ceil(size(crams, "GB") / 10)])
    command <<<
            bash ~{monitoring_script} | tee monitoring.log >&2 &
            source ~/.bashrc
            conda activate genomics.py3

            ~{"tar -zxf "+cache_tarball}

            REF_CACHE=cache/%2s/%2s/ REF_PATH='.' samtools merge -@ ~{cpus_to_use} --reference ~{references.ref_fasta} - ~{sep=" " crams} | \
                samtools reheader - --command 'sed "s/SM:[^\t]*/SM:~{sample_name}/g"' | \
                samtools view -h - -@ ~{cpus_to_use} -C --reference ~{references.ref_fasta} -o ~{output_cram_name} --output-fmt-option embed_ref
            samtools index ~{output_cram_name}
    >>>
    output {
            File monitoring_log = "monitoring.log"
            File output_cram = "~{output_cram_name}"
            File output_cram_index = "~{output_cram_name}.crai"
    }
    runtime {
            memory: "8GB"
            disks: "local-disk " + (ceil(size(cache_tarball, "GB") + size(crams, "GB")) * 3 + 10) + " HDD"
            docker: docker
            noAddress: no_address
            cpu: "~{cpus_to_use}"
            preemptible: preemptible_tries

    }
}

task MergeBams {
    input {
        File monitoring_script
        Array[File] inputs
        String output_prefix
        String docker
        Int disk_size = ceil(2*size(inputs, "GB")) + 20
        Int preemptible_tries
        Boolean no_address
    }
    command<<<
        set -eo pipefail
        bash ~{monitoring_script} > monitoring.log &
        samtools merge -c -@ 8 ~{output_prefix}.bam ~{sep=' ' inputs}
        samtools index ~{output_prefix}.bam
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "16 GB"
        cpu: "8"
        disks: "local-disk " + disk_size + " LOCAL"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File monitoring_log = "monitoring.log"
        File output_bam = "~{output_prefix}.bam"
        File output_bam_index = "~{output_prefix}.bam.bai"
    }
}

task MergeVCFs {
  input {
    File monitoring_script
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int disk_size = ceil(2*size(input_vcfs,"GB")+5)
    Int preemptible_tries
    String docker
    String gitc_path = "/usr/gitc/"
    Boolean no_address
  }
  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail
    java -Xms9000m -jar ~{gitc_path}picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
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

task FilterVcfWithBcftools {
    input {
        String docker
        String? base_file_name
        File monitoring_script
        File input_vcf
        String? bcftools_extra_args
        Array[File]? exclude_regions
        Array[File]? include_regions
        Int preemptible_tries = 1
        Int disk_size = ceil((size(input_vcf, "GB") + size(select_first([include_regions, []]), "GB") + size(select_first([exclude_regions, []]), "GB")) * 2 + 10)
        Int memory_gb = 4 + ceil(ceil(size(select_first([include_regions, []]), "GB") + size(select_first([exclude_regions, []]), "GB")) * 0.25)
        Int cpus = 1 + length(select_first([include_regions, []])) + length(select_first([exclude_regions, []]))  # a process for each bcftools command
    }
        String output_base_name = select_first([base_file_name, basename(input_vcf, ".vcf.gz")])
        String output_vcf_filename = output_base_name + ".filtered.vcf.gz"
        Boolean defined_include_regions = length(select_first([include_regions, []])) > 0
        Boolean defined_exclude_regions = length(select_first([exclude_regions, []])) > 0
    
    meta {
        description : "Filter input vcf file using bcftools, with 'bcftools view' args and with genomic regions to include and/or exclude."
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "preemptible_tries"
        ]}
    }
    parameter_meta {
        docker: {
            help: "Docker image to use for this task, must be able to run 'bcftools' with no preceding commands.",
            type: "String",
            category: "input_required"
        }
        base_file_name: {
            help: "Base file name for output files.",
            type: "String", 
            category: "input_optional"
        }
        input_vcf: {
            help: "Input vcf file. Note that an index is not needed.",
            type: "File",
            category: "input_required"
        }
        bcftools_extra_args: {
            help: "Extra arguments to pass to bcftools view. For example, use '-f PASS' to only include variants with a PASS filter.",
            type: "String",
            category: "input_optional"
        }
        exclude_regions: {
            help: "Regions to exclude from the output vcf. Supported formats are bed, bed.gz, vcf, vcf.gz. VCF exclusion is done using bcftools view -T ^regions, by position and not by ref and alt.",
            type: "Array[File]",
            category: "input_optional"
        }
        include_regions: {
            help: "Regions to include in the output vcf. Supported formats are bed, bed.gz. Inclusion is done using 'bcftools view -T'.",
            type: "Array[File]",
            category: "input_optional"
        }
        preemptible_tries: {
            help: "Number of times to retry this task if it fails.",
            type: "Int",
            category: "input_optional"
        }
        disk_size: {
            help: "Size of the local disk to use for this task, in GB. By default it is calculated from the input file sizes.",
            type: "Float",
            category: "input_optional"
        }
        memory_gb: {
            help: "Amount of memory to use for this task, in GB. Default is 4 (GB).",
            type: "Int",
            category: "input_optional"
        }
        cpus: {
            help: "Number of cpus to use for this task. Default is 1.",
            type: "Int",
            category: "input_optional"
        }
        monitoring_script: {
            help: "UG resource monitoring script.",
            type: "File",
            category: "input_required"
        }
    }
    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        bcftools view \
            --threads ~{cpus} \
            ~{bcftools_extra_args} \
            ~{input_vcf} \
            ~{true=" | bcftools view - -T " false="" defined_include_regions}~{sep=" | bcftools view - -T " include_regions} \
            ~{true=" | bcftools view - -T ^" false="" defined_exclude_regions}~{sep=" | bcftools view - -T ^" exclude_regions} \
            -Oz \
            -o ~{output_vcf_filename}
        bcftools index -t ~{output_vcf_filename}
    >>>
    output {
        File monitoring_log = "monitoring.log"
        File output_vcf = "~{output_vcf_filename}"
        File output_vcf_index = "~{output_vcf_filename}.tbi"
    }
    runtime {
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        cpu: "~{cpus}"
        preemptible: preemptible_tries
    }
}

task GetMeanCoverageFromSorterStats {
    input{
      File sorter_json_stats_file
      String docker
      Int preemptible_tries
      File monitoring_script
    }
    String output_file = basename(sorter_json_stats_file, ".json") + ".mean_coverage.txt"
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        sorter_stats_to_mean_coverage \
            --sorter-stats-json "~{sorter_json_stats_file}" \
            --output-file "~{output_file}"

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        docker: docker
        cpu: "1"
    }
    output {
        Int mean_coverage = read_int("~{output_file}")
        File monitoring_log = "monitoring.log"
    }
  }

task CopyFiles {
    input {
        Array[File] input_files
        String docker
    }
    command <<<
    >>>
    runtime {
        docker: docker
        preemptible: 1
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk " +ceil(2*size(input_files,"GB") + 1) + " HDD"
        noAddress: true
    }
    output {
        Array[File] output_files = input_files
    }
}
# creates a bed file from a VCF file and expands if needed
task VcfToIntervalListAndBed {
    input {
        File monitoring_script
        File input_vcf
        Int bedtools_expand
        File reference_dict
        File reference_index
        String base_file_name
        String docker
        Int disk_size = ceil(2*size(input_vcf, "GB")) + 20
        Int preemptible_tries
        Boolean no_address
    }
    command <<< 
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source /opt/conda/etc/profile.d/conda.sh
        conda activate genomics.py3

        picard VcfToIntervalList \
            I=~{input_vcf} \
            O=step1.interval_list 

        picard IntervalListToBed \
            I=step1.interval_list \
            O=step2.bed
        awk 'BEGIN {FS=OFS="\t"} {print $1, $2}' ~{reference_index} > genome_file.txt
        bedtools slop -i step2.bed -g genome_file.txt -b ~{bedtools_expand} | bedtools merge -i - > ~{base_file_name}.bed
        picard BedToIntervalList \
            I=~{base_file_name}.bed \
            O=~{base_file_name}.interval_list \
            SD=~{reference_dict}
    >>>

    runtime {
        memory: "2 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_tries
        noAddress: no_address
    }

    output {
        File interval_list = "~{base_file_name}.interval_list"
        File bed = "~{base_file_name}.bed"
        File monitoring_log = "monitoring.log"
    }
}