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
        Array[File] input_cram_bam_list
        TrimmerParameters parameters
        File? cache_tarball
        String base_file_name
        String docker
        Int disk_size = ceil(3 * size(input_cram_bam_list, "GB") + size(cache_tarball, "GB") + 20)
        Int cpus
        Int preemptible_tries
        Boolean no_address
    }
    Int memory_gb = select_first([parameters.memory_gb, 8])

    String trimmer_mode = if defined(parameters.untrimmed_reads_action) then "--~{parameters.untrimmed_reads_action}" else ""
    String trimmer_format_flag = if defined(parameters.format) then "--format=\"~{parameters.format}\"" else ""
    String trimmer_extra_args = if defined(parameters.extra_args) then " ~{parameters.extra_args}" else ""
    Array[File] pattern_files = select_first([parameters.pattern_files, []])
    String formats_base_path = "trimmer/formats"
    String trimmer_prefix_sep = if defined(parameters.filename_prefix_sep) then "~{parameters.filename_prefix_sep}" else "_"

    String output_file_name_prefix = if defined(parameters.output_demux_format) then "~{base_file_name}~{trimmer_prefix_sep}~{parameters.output_demux_format}" else base_file_name
    String output_file_name_suffix = ".trimmed.ucram"
    String output_file_name = "~{output_file_name_prefix}~{output_file_name_suffix}"
    String trimmer_stats_file = "~{output_file_name_prefix}.trimmer_stats.csv"
    String output_trimmed_failed_file_name = if defined(parameters.output_failed_file_name_suffix) then "~{base_file_name}~{trimmer_prefix_sep}~{parameters.output_failed_file_name_suffix}" else "PLACEHOLDER_THAT_SHOULD_NOT_EXIST"
    String failure_file_args = if defined(parameters.output_failed_file_name_suffix) then "--failure-file=~{output_trimmed_failed_file_name}" else ""
    String output_fc_file_name = "~{output_file_name_prefix}.failure_codes.csv"
    String failure_codes_args = if trimmer_mode == "--discard" then "--failure-code-file ~{output_fc_file_name}" else "--failure-code-tag fc --failure-code-file ~{output_fc_file_name}"  # --discard and error codes don't mix
    String failure_read_group_args = if defined(parameters.failure_read_group) then "--failure-field rg:Z:~{parameters.failure_read_group}" else ""
    String minor_read_group_args = if defined(parameters.minor_read_group) then "--output-field rg:Z:~{parameters.minor_read_group}" else ""
    File? cram_reference = parameters.cram_reference


    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir -p ~{formats_base_path}

        # copy external pattern files to the formats base path
        for pattern_file in ~{sep=" " pattern_files}; do
            cp "$pattern_file" ~{formats_base_path}
        done

        # update the command according to the input file type (bam/cram)
        filename="~{input_cram_bam_list[0]}"
        extension="${filename##*.}"

        # this next part has code duplication, the parameters are the same for both cram and bam, but the input is different, make sure to keep them in sync
        if [ "$extension" = "cram" ]; then
            trimmer \
            --input=~{sep=" --input=" input_cram_bam_list} \
            --description=~{parameters.formats_description} \
            ~{trimmer_format_flag} \
            --statistics=~{trimmer_stats_file} \
            --directory=~{formats_base_path} \
            --skip-unused-pattern-lists=true \
            ~{trimmer_mode} \
            ~{failure_codes_args} \
            ~{failure_read_group_args} \
            ~{minor_read_group_args} \
            --progress \
            --vector \
            --nthreads=~{cpus} \
            ~{failure_file_args} \
            ~{"--reference=" + cram_reference} \
            --cram true \
            --output ~{output_file_name} \
            ~{trimmer_extra_args}
        else  # bam extension
            echo "~{sep='\n'input_cram_bam_list}" > bam_list.txt
            ~{"tar -zxf "+cache_tarball}

            export REF_CACHE=cache/%2s/%2s/
            export REF_PATH='.'

            samtools cat -b bam_list.txt | \
            samtools view -h -@ 2 - | \
            trimmer \
            --description=~{parameters.formats_description} \
            ~{trimmer_format_flag} \
            --statistics=~{trimmer_stats_file} \
            --directory=~{formats_base_path} \
            --skip-unused-pattern-lists=true \
            ~{trimmer_mode} \
            ~{failure_codes_args} \
            ~{failure_read_group_args} \
            ~{minor_read_group_args} \
            --progress \
            --vector \
            --nthreads=~{cpus} \
            ~{failure_file_args} \
            ~{"--reference=" + cram_reference} \
            --cram true \
            --output ~{output_file_name} \
            ~{trimmer_extra_args}
        fi

        
        OUT_CRAM_SUFFIX=~{output_file_name_suffix} 
        # If remove_small_files is set to true then remove trimmed ucrams under 500M in file size
        
        OUT_CRAM_SUFFIX=~{output_file_name_suffix} 
        # If remove_small_files is set to true then remove trimmed ucrams under 500M in file size
        if [ ~{parameters.remove_small_files} = true ]; then
            find . -type f -size -500M -name "*${OUT_CRAM_SUFFIX}" | xargs -I {} rm {}
        fi

        # To make sure an output cram was created, we check for any files ending in the specified suffix
        # The -type f ensures we only match regular files.
        # -maxdepth 1 ensures we only look in the current directory (optional).
        find . -maxdepth 1 -type f -name "*${OUT_CRAM_SUFFIX}"
        if ! find . -maxdepth 1 -type f -name "*${OUT_CRAM_SUFFIX}" | grep -q .; then
            echo "ERROR: No file found ending with '${OUT_CRAM_SUFFIX}', this indicates that no output cram was created, likely because all the reads did not match on Trimmer." >&2
            ls -ltr           
            exit 1
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
        File trimmer_stats = "~{trimmer_stats_file}"
        File trimmer_failure_codes_csv = "~{output_fc_file_name}"
        File? output_trimmed_failed_file = "~{output_trimmed_failed_file_name}"
        File monitoring_log = "monitoring.log"
        Array[File?] histogram = glob("*histogram.csv")  # aim for this to be only one file (determined by the extra_args) to simplify downstream processing
        Array[File?] histogram_extra = glob("*histogram_extra.csv")  # aim for this to be only one file (determined by the extra_args) to simplify downstream processing
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
        set -xeuo pipefail
        
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
        set -xeuo pipefail
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
        set -xeuo pipefail
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
