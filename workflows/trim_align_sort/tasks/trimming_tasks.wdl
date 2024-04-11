version 1.0
import "structs.wdl"

# Note regarding dummy_input_for_call_caching
# When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your
# own workspace. This will prevent call caching from failing with "Cache Miss (10 failed copy attempts)". Outside of
# Terra this can be left as the default empty String. This dummy input is only needed for tasks that have no inputs
# specific to the sample being run (such as GetBwaVersion which does not take in any sample data).


# TASK DEFINITIONS
task Trimmer {
    # Input:
    # Trimmer can accept either a cram or a bam as input. Cram is better read natively (without samtools piping) and this can effect the performance.
    # For cram input file, if the refernce is not embedded in the cram you need to provide a reference fasta in trimmer_parameters.cram_reference.
    # Output: output format is always ucram (unmapped cram)
    input {
        File monitoring_script
        File input_cram_bam
        TrimmerParameters trimmer_parameters
        String base_file_name
        String docker
        Int disk_size = ceil(3 * size(input_cram_bam, "GB") + 20)
        Int cpus
        Int preemptible_tries
        Boolean no_address
    }
    Int memory_gb = select_first([trimmer_parameters.memory_gb, 8])

    String trimmer_mode = if trimmer_parameters.untrimmed_reads_action == "" then "" else "--~{trimmer_parameters.untrimmed_reads_action}"
    String trimmer_format_flag = if defined(trimmer_parameters.format) then '--format="~{trimmer_parameters.format}"' else ""
    String trimmer_extra_args = if defined(trimmer_parameters.extra_args) then " ~{trimmer_parameters.extra_args}" else ""
    Array[File] pattern_files = select_first([trimmer_parameters.pattern_files, []])
    String local_description_file = select_first([trimmer_parameters.local_formats_description, "/trimmer/formats/formats.json"])

    String output_file_name_prefix = if defined(trimmer_parameters.output_demux_format) then "~{base_file_name}_~{trimmer_parameters.output_demux_format}" else base_file_name
    String output_file_name_suffix = ".trimmed.ucram"
    String output_file_name = "~{output_file_name_prefix}~{output_file_name_suffix}"
    String output_fc_file_name = "~{output_file_name_prefix}.failure_codes.csv"
    String failure_codes_args = if trimmer_mode == "--discard" then "--failure-code-file ~{output_fc_file_name}" else "--failure-code-tag fc --failure-code-file ~{output_fc_file_name}"  # --discard and error codes don't mix

    
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        python <<CODE
        # validating inputs
        if "~{trimmer_mode}" not in ("", "--filter", "--discard"):
            raise ValueError(f'Trimmer mode for reads that do not match all required segments can be either "" (do nothing), "filter" (mark in sam flag) or "discard"\nGot ~{trimmer_mode}')
        CODE

        # determine the formats base path
        formats_base_path=$(dirname ~{local_description_file})
        mkdir -p "$formats_base_path"

        # copy external pattern files to the formats base path
        for pattern_file in ~{sep=" " pattern_files}; do
            cp "$pattern_file" "$formats_base_path"
        done 

        # determine the description file - required because it can be either a File or a String
        # Determine whether to use the file content or the string
        trimmer_description_file=""
        if [[ -f "~{trimmer_parameters.formats_description}" ]]; then
            trimmer_description_file="~{trimmer_parameters.formats_description}"
        else
            trimmer_description_file="~{local_description_file}"
        fi

        # update the command according to the input file type (bam/cram)
        filename="~{input_cram_bam}"
        extension="${filename##*.}"
        
        if [ "$extension" = "cram" ]; then
            trimmer \
            --input=~{input_cram_bam} \
            --description=$trimmer_description_file \
            ~{trimmer_format_flag} \
            --statistics=trimmer_stats.csv \
            --directory="$formats_base_path" \
            --skip-unused-pattern-lists=true \
            ~{trimmer_mode} \
            ~{trimmer_extra_args} \
            ~{failure_codes_args} \
            --progress \
            --vector \
            --nthreads=~{cpus} \
            ~{if defined(trimmer_parameters.cram_reference) then "--reference" else ""} ~{default="" trimmer_parameters.cram_reference} \
            --cram true \
            --output ~{output_file_name}
        else
            samtools view -h ~{input_cram_bam} -@ ~{cpus} | \
            trimmer \
            --description=$trimmer_description_file \
            ~{trimmer_format_flag} \
            --statistics=trimmer_stats.csv \
            --directory="$formats_base_path" \
            --skip-unused-pattern-lists=true \
            ~{trimmer_mode} \
            ~{trimmer_extra_args} \
            ~{failure_codes_args} \
            --progress \
            --vector \
            --nthreads=~{cpus} \
            ~{if defined(trimmer_parameters.cram_reference) then "--reference" else ""} ~{default="" trimmer_parameters.cram_reference} \
            --cram true \
            --output ~{output_file_name}
        fi

    >>>
    runtime {
        cpuPlatform: "Intel Skylake"
        cpu: "~{cpus}"
        preemptible: preemptible_tries
        memory: "~{memory_gb} GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        Array[File] trimmed_ucram_list = glob("*~{output_file_name_suffix}")
        File trimmer_stats = "trimmer_stats.csv"
        File trimmer_failure_codes_csv = "~{output_fc_file_name}"
        File monitoring_log = "monitoring.log"
        Array[File?] histogram = glob("*histogram.csv")  # aim for this to be only one file (determined by the extra_args) to simplify downstream processing
        Array[File?] histogram_extra = glob("*histogram_extra.csv")  # aim for this to be only one file (determined by the extra_args) to simplify downstream processing
    }
}

task TrimmerGenerateFormatsJson {
    input {
        SimpleReadTrimmingParameters read_trimming_parameters
        File monitoring_script
        String docker
        Int preemptible_tries
        Int disk_size = 1
        Boolean no_address
    }
    String output_trimmer_formats_filename = "trimmer_formats.json"
    String min_insert_len = if defined(read_trimming_parameters.min_insert_length) then "--min-insert-length ~{read_trimming_parameters.min_insert_length}" else ""
    String max_insert_len = if defined(read_trimming_parameters.max_insert_length) then "--max-insert-length ~{read_trimming_parameters.max_insert_length}" else ""
    String umi_length_5p = if defined(read_trimming_parameters.umi_length_5p) then "--umi-length-5p ~{read_trimming_parameters.umi_length_5p}" else ""
    String umi_length_3p = if defined(read_trimming_parameters.umi_length_3p) then "--umi-length-3p ~{read_trimming_parameters.umi_length_3p}" else ""

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        /trimmer/generate_format_json.py ~{min_insert_len} ~{max_insert_len} ~{umi_length_5p} ~{umi_length_3p} \
            --adapter-5p ~{read_trimming_parameters.adapter_5p} \
            --adapter-5p-required \
            --adapter-5p-min-overlap ~{read_trimming_parameters.min_overlap_5p} \
            --adapter-5p-max-error-rate ~{read_trimming_parameters.max_error_rate_5p} \
            --adapter-3p ~{read_trimming_parameters.adapter_3p} \
            --adapter-3p-min-overlap ~{read_trimming_parameters.min_overlap_3p} \
            --adapter-3p-max-error-rate ~{read_trimming_parameters.max_error_rate_3p} \
            --output-file ~{output_trimmer_formats_filename}

        cat  ~{output_trimmer_formats_filename}
    >>>
    runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "5 GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    noAddress: no_address
    maxRetries: 1
    }
    output {
        File trimmer_formats = "~{output_trimmer_formats_filename}"
        File monitoring_log = "monitoring.log"
    }
}

task TrimmerAggregateStats {
    input {
        Array[File] trimmer_stat_multiple_files
        String base_file_name

        File monitoring_script
        String docker
        Int preemptible_tries
        Int disk_size = ceil(1)
        Boolean no_address
    }

    String aggregated_trimmer_stats_filename = "~{base_file_name}.aggregated_trimmer_stats.csv"

    command <<<
        set -eo pipefail
        
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        conda activate genomics.py3

        python <<CODE
        import pandas as pd
        from functools import reduce
        
        stats_files = "~{sep="," trimmer_stat_multiple_files}".split(",")
        stats_df = [pd.read_csv(stat).set_index(['read group', 'segment index', 'segment label', 'is required', 'is trimmed']) for stat in stats_files]
        aggregate_stats = reduce(lambda a, b: a.add(b, fill_value=0), stats_df).reset_index()
        aggregate_stats.to_csv("~{aggregated_trimmer_stats_filename}")
        
        CODE
    >>>
    runtime {
        preemptible: preemptible_tries
        cpu: "1"
        memory: "5 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File aggregated_trimmer_stats = "~{aggregated_trimmer_stats_filename}"
        File monitoring_log = "monitoring.log"
    }
}

task CutadaptMarkAdapter {
    input {
        File monitoring_script
        String docker
        String? base_file_name
        Boolean no_address
        File input_ubam
        SimpleReadTrimmingParameters read_trimming_parameters
        Int preemptible_tries
        Int disk_size = ceil(3 * size(input_ubam, "GB") + 20)
        Int memory_gb
        Int cpus
    }

    String env='cutadaptenv'
    Int umi_length_5p = select_first([read_trimming_parameters.umi_length_5p, 0])
    Int umi_length_3p = select_first([read_trimming_parameters.umi_length_3p, 0])
    String output_base_name = select_first([base_file_name, basename(input_ubam, ".bam")])

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        . ~/.bashrc
        conda activate genomics.py3
        export PATH=$PATH:$CONDA_PREFIX/lib/python3.7/site-packages/ugvc/bash/

        . find_adapter_coords.sh ~{input_ubam} ~{env} \
            "~{read_trimming_parameters.adapter_5p}" "~{read_trimming_parameters.adapter_3p}" \
            ~{umi_length_5p} ~{umi_length_3p} \
            ~{read_trimming_parameters.max_error_rate_5p} ~{read_trimming_parameters.max_error_rate_3p} \
            ~{read_trimming_parameters.min_overlap_5p} ~{read_trimming_parameters.min_overlap_3p}
    >>>

    output {
        File output_ubam = "~{output_base_name}.with_adapter_tags.bam"
        File report = "~{output_base_name}.cutadapt_report.txt"
        File stats_json = "~{output_base_name}.cutadapt.json"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "~{memory_gb} GB"
        docker: docker
        noAddress: no_address
        preemptible: preemptible_tries
        disks: "local-disk " + disk_size + " HDD"
    }
}

task CutadaptTrimAdapter {
    input {
        File monitoring_script
        String docker
        Boolean no_address
        File input_ubam
        SimpleReadTrimmingParameters read_trimming_parameters
        String? base_file_name
        String gitc_path = "/usr/gitc/"
        Int preemptible_tries
        Int disk_size = ceil(3 * size(input_ubam, "GB") + 20)
    }
    String output_base_name = select_first([base_file_name, basename(input_ubam, ".with_adapter_tags.bam")])
    Int min_insert_length = select_first([read_trimming_parameters.min_insert_length, 1])

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        java -Xmx10g -jar ~{gitc_path}GATK_ultima.jar ClipReads --input ~{input_ubam} \
        -O ~{output_base_name}.trimmed.bam \
        -os ~{output_base_name}.trimming.report.txt \
        -CR  HARDCLIP_BASES \
        --clip-adapter \
        --min-read-length-to-output ~{min_insert_length}
    >>>

    output {
        File trimmed_ubam = "~{output_base_name}.trimmed.bam"
        File output_report = "~{output_base_name}.trimming.report.txt"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        memory: "10 GB"
        docker: docker
        noAddress: no_address
        preemptible: preemptible_tries
        disks: "local-disk " + disk_size + " HDD"
    }
}
