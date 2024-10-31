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

    Int mapq = select_first([mapq_override,0])
    String demux_output_path = "demux_output/"
    String align_flag = if defined(sorter_params.aligned) then "--align=~{sorter_params.aligned}" else "--align=true"
    String output_group = select_first([sorter_params.output_group, base_file_name])
    String output_path = select_first([sorter_params.output_path, "{outputGroup}/{outputGroup}"])

    String timestamp = "000" # use this to override sorter output timestemp (used only for the output folder structure)
    String sorter_out_dir = "sorter_output"
    String sorter_output_path = "~{sorter_out_dir}/~{base_file_name}-~{timestamp}" 
    
    # Sorter dir strucutre:
    # <sorter_out_dir>/
    #       <run_id>_<timestamp>/
    #                   <output_path>
    #
    # - sorter_out_dir: output directory for sorter (pass to sorter as --output-dir).
    # - run_id: run id is override with base_file_name (pass to sorter as --runid).
    # - timestamp: timestamp for the sorter output (pass to sorter as --timestamp).
    # - output_path: DEMUX pararmter, determines the filesystem path for a user group (pass to DEMUX as --output-path). 
    #               Read sorter documentation for more information.

    String output_unmatched_crm = "~{sorter_output_path}/~{base_file_name}_unmatched.cram"

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Print instance available memory
        echo "[DEBUG] Available memory:"
        free -g -h -t

        # Open coverage intervals if given
        file_name="~{default="" sorter_params.coverage_intervals}"
        coverage_intervals_flag=""
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
            --output-group=~{output_group} \
            --output-path=~{output_path} \
            ~{"--umi=" + sorter_params.umi_tag} \
            ~{align_flag} \
            $coverage_intervals_flag \
            ~{default="" sorter_params.demux_extra_args}

        sorter \
            --runid=~{base_file_name} \
            --input-dir=~{demux_output_path} \
            --output-dir=~{sorter_out_dir} \
            --nthreads ~{cpu} \
            --progress \
            --timestamp=~{timestamp} \
            ~{default="" sorter_params.sort_extra_args}

        ls -R ~{sorter_output_path}/
        q_table.py ~{sorter_output_path}/

        echo "Moving files around for better organization"

        # Move fastq files to fastq directory (to separate them from the sample.fastq files)
        mkdir -p fastq
        for file in $(ls ~{sorter_output_path}/**/*.fastq.gz | grep -v '_sample'); do
            echo "Moving $file to fastq/"
            mv $file fastq/
        done
        
        # Move the output cram file to seperate it from sub-sample or unmatched crams.
        mkdir -p cram
        for file in $(ls ~{sorter_output_path}/**/*.cram | grep -v '_sample' | grep -v 'unmatched'); do
            echo "Moving $file to cram/ folder"
            mv $file cram/
        done

        # Move the output cram.crai file to seperate it from sub-sample or unmatched files.
        for file in $(ls ~{sorter_output_path}/**/*.cram.crai | grep -v '_sample' | grep -v 'unmatched'); do
            echo "Moving $file to cram/ folder"
            mv $file cram/
        done

        # Rename the unmatched cram file (if exists, we assume there is only one file)
        for file in $(ls ~{sorter_output_path}/**/*unmatched*.cram); do
            echo "Moving $file to ~{output_unmatched_crm}"
            mv $file ~{output_unmatched_crm}
        done

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
        Array[File?] sorted_cram = glob("cram/*.cram")
        Array[File?] sorted_cram_index = glob("cram/*.cram.crai")
        File? unmatched_cram = "~{output_unmatched_crm}"
        Array[File] sorter_stats_csv = glob("~{sorter_output_path}/**/*.csv")
        Array[File] sorter_stats_json = glob("~{sorter_output_path}/**/*.json")
        Array[File?] sorter_out_bedgraph_mapq0 = glob("~{sorter_output_path}/**/*_0.bedGraph.gz")
        Array[File?] sorter_out_bedgraph_mapq1 = glob("~{sorter_output_path}/**/*_1.bedGraph.gz")
        Array[File?] sorter_out_bedgraph_extents = glob("~{sorter_output_path}/**/*_0.bedGraph.gz.extents.tsv")
        Array[File?] fastq_files = glob("fastq/*") # in case output format is fastq (like in single cell)
        Array[File?] sub_sampled_output = glob("~{sorter_output_path}/**/*_sample.*")
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
        set -xeuo pipefail
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