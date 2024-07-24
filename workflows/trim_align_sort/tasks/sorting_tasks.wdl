version 1.0
import "structs.wdl"

# TASK DEFINITIONS
task Sorter {
    input {
        File input_file
        String base_file_name
        File reference_fasta
        Int? mapq_override
        SorterParams sorter_params

        File monitoring_script

        Int cpu
        Int preemptible_tries
        String docker
    }

    Int disk_size = ceil(4.5*size(input_file, "GB") + 3*size(reference_fasta, "GB") + 20 )
    Int local_ssd_size = 375
    # add one local ssd disk (375 GB) when spare disk is < 50 GB
    # the actual size that will be required to google will be a multiple of 375
    Int rounded_disk_size = ceil(disk_size/local_ssd_size) * local_ssd_size
    Int total_disk_size = if (rounded_disk_size - disk_size) < 50 then disk_size + local_ssd_size else disk_size

    Int local_ssd_threshold = 3000
    Int mapped_bam_size_local_ssd = if total_disk_size < local_ssd_threshold then total_disk_size else 9000

    Int memory_gb = select_first([sorter_params.memory_gb, 64])
    Int preemptible_tries_final = if (size(input_file, "GB") < 250) then preemptible_tries else 0

    String demux_output_path = "demux_output/"
    String output_path_demux = select_first([sorter_params.demux_output_path,"{runID}-{outputGroup}/{runID}-{outputGroup}"])
    String output_group_demux = select_first([sorter_params.demux_output_group,base_file_name])

    String timestamp = "" # use this to override sorter output timestemp
    String run_id_override = "output" # use this to override sorter output runid (used only for the output folder structure)

    String sorter_out_dir = "sorter_output"
    String sorter_output_path = "~{sorter_out_dir}/~{run_id_override}-~{timestamp}"
    String demux_align_flag = if defined(sorter_params.demux_align) then "--align=~{sorter_params.demux_align}" else "--align=true"
    Int mapq = select_first([mapq_override,0])

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Print instance available memory
        free -g -h -t

        # Open coverage intervals if given
        file_name=~{default="" sorter_params.coverage_intervals}
        if [[ -n $file_name && -f $file_name ]]; then
            echo "Unzipping coverage intervals"
            tar xvzf ~{sorter_params.coverage_intervals} -C .
            tsv_file=$(find . -name "*.tsv")
            echo "Coverage intervals file: $tsv_file"
            coverage_intervals_flag="--intervals=$tsv_file"
            echo "Coverage intervals flag: $coverage_intervals_flag"
        else
            echo "WARNING: No coverage intervals file given, respective statistics will not be calculated"
        fi
        mkdir -p ~{sorter_out_dir}

        samtools view -h ~{input_file} -@ ~{cpu} ~{"-q " + mapq} -T ~{reference_fasta}|\
        demux \
            --input=- \
            --output-dir=~{demux_output_path} \
            --runid=~{base_file_name} \
            --nthreads ~{cpu} \
            --progress \
            --reference ~{reference_fasta} \
            --mark-duplicates=~{sorter_params.mark_duplicates} \
            --output-group=~{output_group_demux} \
            --output-path=~{output_path_demux} \
            ~{"--umi=" + sorter_params.umi_tag} \
            ~{demux_align_flag} \
            $coverage_intervals_flag \
            ~{default="" sorter_params.demux_extra_args}

        sorter \
            --runid=~{run_id_override} \
            --input-dir=~{demux_output_path} \
            --output-dir=~{sorter_out_dir} \
            --nthreads ~{cpu} \
            --progress \
            --timestamp=~{timestamp} \
            ~{default="" sorter_params.sort_extra_args}

        ls -R ~{sorter_output_path}/
        q_table.py ~{sorter_output_path}/

    >>>
    runtime {
        cpuPlatform: "Intel Skylake"
        cpu: "~{cpu}"
        preemptible: preemptible_tries_final
        memory: "~{memory_gb} GB"
        disks: "local-disk " + ceil(mapped_bam_size_local_ssd) + " LOCAL"
        docker: docker
        noAddress: false
    }
    output {
        File monitoring_log = "monitoring.log"
        File sorted_cram = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*.cram")[0]
        File sorted_cram_index = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*.cram.crai")[0]
        File sorter_stats_csv = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*.csv")[0]
        File sorter_stats_json = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*.json")[0]
        File sorter_out_bedgraph_mapq0 = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*_0.bedGraph.gz")[0]
        File sorter_out_bedgraph_mapq1 = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*_1.bedGraph.gz")[0]
        File sorter_out_bedgraph_extents = glob("~{sorter_output_path}/~{base_file_name}*/~{base_file_name}*_0.bedGraph.gz.extents.tsv")[0]
    }
}

task Demux {
    input {
        File input_file
        String base_file_name
        String? output_dir
        String? demux_extra_args
        String? samtools_extra_args
        String? output_format # cram/fastq/sam . According to demux, default is cram
        File? reference_fasta # when input is CRAM, adding a reference file imporve performance

        Int cpu
        File monitoring_script
        Int preemptible_tries
        String docker
    }
    Int disk_size = ceil(3 * size(input_file, "GB") + 20)
    Int local_ssd_size = 375
    # add one local ssd disk (375 GB) when spare disk is < 50 GB
    # the actual size that will be required to google will be a multiple of 375
    Int rounded_disk_size = ceil(disk_size/local_ssd_size) * local_ssd_size
    Int total_disk_size = if (rounded_disk_size - disk_size) < 50 then disk_size + local_ssd_size else disk_size

    Int local_ssd_threshold = 3000
    Int local_ssd_size_ask = if total_disk_size < local_ssd_threshold then total_disk_size else 9000

    String output_dir_override = select_first([output_dir, "demux_output/"])
    Boolean is_output_fastq = if (defined(output_format) && select_first([output_format]) == "fastq") then true else false
    String reference_flag = if defined(reference_fasta) then "--reference=~{reference_fasta}" else ""
    String reference_flag_for_fastq = if is_output_fastq then "" else reference_flag # in case output_format is fastq, don't use reference

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 & 

        # If output_format is fastq, discard QC failed reads (mostly ones that failed in Trimmer) 0x200, secondary alignments 0x100, and supplementary alignment 0x800
        samtools view -h ~{input_file} -@ ~{cpu} ~{"-T " + reference_fasta} ~{true="-F 0x100 -F 0x200 -F 0x800" false="" is_output_fastq} ~{default="" samtools_extra_args} | \
        demux \
            --input=- \
            --output-dir=~{output_dir_override} \
            --output-group=~{base_file_name} \
            --output-path={runID}-{outputGroup}/{runID}-{outputGroup} \
            --runid= \
            --nthreads=~{cpu} \
            --progress \
            ~{reference_flag_for_fastq} \
            ~{"--output-format=" + output_format} \
            ~{true="--align=false" false="" is_output_fastq} \
            ~{default="" demux_extra_args}

        ls -R ~{output_dir_override}
    >>>
    runtime {
        cpuPlatform: "Intel Skylake"
        cpu: "~{cpu}"
        preemptible: preemptible_tries
        memory: "16 GB"
        disks: "local-disk " + ceil(local_ssd_size_ask) + " LOCAL"
        docker: docker
        noAddress: false
    }
    output {
        File? output_fastq = "~{output_dir_override}/~{base_file_name}.fastq.gz"
        File? output_cram = "~{output_dir_override}/~{base_file_name}.cram"
        File? output_sam = "~{output_dir_override}/~{base_file_name}.sam"
        File monitoring_log = "monitoring.log"
    }
}