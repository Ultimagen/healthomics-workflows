version 1.0
import "structs.wdl"

# TASK DEFINITIONS
task Demux {
    input {
        Array[File] input_cram_bam_list
        File? cache_tarball
        String base_file_name
        File reference_fasta
        Int? mapq_override
        SorterParams sorter_params

        File monitoring_script

        Int cpu
        Int preemptible_tries
        String docker
    }

    Int disk_size = ceil(2.5*size(input_cram_bam_list, "GB") + 3*size(reference_fasta, "GB") + 20 )
    Int local_ssd_size = 375
    # add one local ssd disk (375 GB) when spare disk is < 50 GB
    # the actual size that will be required to google will be a multiple of 375
    Int rounded_disk_size = ceil(disk_size/local_ssd_size) * local_ssd_size
    Int total_disk_size = if (rounded_disk_size - disk_size) < 50 then disk_size + local_ssd_size else disk_size

    Int local_ssd_threshold = 3000
    Int mapped_bam_size_local_ssd = if total_disk_size < local_ssd_threshold then total_disk_size else 9000

    Int cpu_samtools = ceil(cpu / 5)
    Int cpu_demux_tmp = cpu - cpu_samtools
    Int cpu_demux = if cpu_demux_tmp < 1 then 1 else cpu_demux_tmp

    Int memory_gb = select_first([sorter_params.demux_memory_gb, 2*cpu])
    Int preemptible_tries_final = if (size(input_cram_bam_list, "GB") < 250) then preemptible_tries else 0

    Int mapq = select_first([mapq_override,0])
    String demux_output_path = "demux_output/"
    String align_flag = if defined(sorter_params.aligned) then "--align=~{sorter_params.aligned}" else "--align=true"
    String output_group = select_first([sorter_params.output_group, base_file_name])
    String output_path = select_first([sorter_params.output_path, "{outputGroup}/{outputGroup}"])
    String reference_fasta_base = basename(reference_fasta)
    

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

        mkdir reference_folder/
        ln -s ~{reference_fasta} reference_folder/~{reference_fasta_base}


        echo "~{sep='\n'input_cram_bam_list}" > bam_list.txt
        ~{"tar -zxf "+cache_tarball}
        
        # for compatibility with the old image where ua was in /ua/ua and not in PATH
        export REF_CACHE=cache/%2s/%2s/ 
        export REF_PATH='.' 
        
        samtools cat -b bam_list.txt | \
        samtools view -h -@ ~{cpu_samtools} ~{"-q " + mapq} -T reference_folder/~{reference_fasta_base} - | \
        demux \
            --input=- \
            --output-dir=~{demux_output_path} \
            --runid=~{base_file_name} \
            --nthreads ~{cpu_demux} \
            --progress \
            --reference reference_folder/~{reference_fasta_base} \
            --mark-duplicates=~{sorter_params.mark_duplicates} \
            --output-group=~{output_group} \
            --output-path=~{output_path} \
            ~{"--umi=" + sorter_params.umi_tag} \
            ~{align_flag} \
            $coverage_intervals_flag \
            ~{default="" sorter_params.demux_extra_args}

        ls -R ~{demux_output_path}/

        echo "Extracting required memory for Sorter:"
        awk -F, 'NR > 1 {for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+:[0-9]+$/) {split($i, parts, ":"); if (parts[1] > max) max = parts[1]}} END {print max}' ~{demux_output_path}/*region-counters.csv | tee max_region_size.txt

    >>>
    runtime {
        cpuPlatform: "Intel Skylake"
        cpu: "~{cpu}"
        preemptible: preemptible_tries_final
        memory: "~{memory_gb} GiB"
        disks: "local-disk " + ceil(mapped_bam_size_local_ssd) + " LOCAL"
        docker: docker
        noAddress: false
    }
    output {
        File monitoring_log = "monitoring.log"
        Int max_region_size = read_int("max_region_size.txt")
        Array[File] demux_output = glob("~{demux_output_path}/*.*")
    }
}

task Sorter {
    input {
        Array[File] demux_output
        Int max_region_size
        String base_file_name
        File reference_fasta
        SorterParams sorter_params

        File monitoring_script

        Int preemptible_tries
        String docker
    }

    Int disk_size = ceil(4.5*size(demux_output, "GB") + 3*size(reference_fasta, "GB") + 20 )
    Int local_ssd_size = 375
    # add one local ssd disk (375 GB) when spare disk is < 50 GB
    # the actual size that will be required to google will be a multiple of 375
    Int rounded_disk_size = ceil(disk_size/local_ssd_size) * local_ssd_size
    Int total_disk_size = if (rounded_disk_size - disk_size) < 50 then disk_size + local_ssd_size else disk_size

    Int local_ssd_threshold = 3000
    Int mapped_bam_size_local_ssd = if total_disk_size < local_ssd_threshold then total_disk_size else 9000

    Float memory_required = if defined(sorter_params.umi_tag) then max_region_size / 0.4 else max_region_size / 0.8
    Int minimum_memory_required_gb = ceil(3 * memory_required / 1073741824) + 2  # 1073741824 is 1GB in bytes, need at least 3x the minimum region
    Int memory_gb_ = select_first([sorter_params.memory_gb, 64])
    Int memory_gb = if memory_gb_ < minimum_memory_required_gb then minimum_memory_required_gb else memory_gb_
    Int cpu = ceil(memory_gb / 3) + 1  # 3GB per CPU
    Int preemptible_tries_final = if (size(demux_output, "GB") < 250) then preemptible_tries else 0

    String timestamp = "000" # use this to override sorter output timestemp (used only for the output folder structure)
    String sorter_out_dir = "sorter_output"
    String sorter_output_path = "~{sorter_out_dir}/~{base_file_name}-~{timestamp}" 
    String reference_fasta_base = basename(reference_fasta)
    
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

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Print instance available memory
        echo "[DEBUG] Available memory:"
        free -g -h -t

        mkdir -p ~{sorter_out_dir}

        mkdir reference_folder/
        ln -s ~{reference_fasta} reference_folder/~{reference_fasta_base}

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

        sorter \
            --runid=~{base_file_name} \
            --input-dir=$(dirname ~{demux_output[0]}) \
            --output-dir=~{sorter_out_dir} \
            --nthreads ~{cpu} \
            --progress \
            --timestamp=~{timestamp} \
            ~{default="" sorter_params.sort_extra_args}

        ls -R ~{sorter_output_path}/

        # Create sorter stats summary
        q_table.py ~{sorter_output_path}/

        echo "Moving files around for better organization"

        output_found=0  # we will make sure an output file, cram or fastq, was created

        # Move fastq files to fastq directory (to separate them from the sample.fastq files)
        mkdir -p fastq
        for file in $(ls ~{sorter_output_path}/**/*.fastq.gz 2>/dev/null | grep -v '_sample'); do
            echo "Moving $file to fastq/"
            mv "$file" fastq/
            output_found=1
        done
        
        # Move the output cram file to seperate it from sub-sample or unmatched crams.
        mkdir -p cram
        for file in $(ls ~{sorter_output_path}/**/*.cram 2>/dev/null | grep -v '_sample' | grep -v '_unmatched.cram'); do
            echo "Moving $file to cram/"
            mv "$file" cram/
            output_found=1
        done

        # Move the output cram.crai file to seperate it from sub-sample or unmatched files.
        for file in $(ls ~{sorter_output_path}/**/*.cram.crai 2>/dev/null | grep -v '_sample' | grep -v '_unmatched.cram'); do
            mv "$file" cram/
            echo "Moving $file to cram/"
        done

        # Rename the unmatched cram file (if exists, we assume there is only one file)
        mkdir -p unmatched
        unmatched_cram_found=0
        for file in $(ls ~{sorter_output_path}/**/*_unmatched.* 2>/dev/null); do
            echo "Moving $file to unmatched/"
            mv "$file" unmatched/
            unmatched_cram_found=1
        done

        # If output_found is still 0, then no main output files were found
        if [ "$output_found" -eq 0 ]; then
            echo "ERROR: No main output file found. Both 'fastq/' and 'cram/' are effectively empty."
            if [ "$unmatched_cram_found" -eq 1 ]; then
                echo "An _unmatched.cram file was found, it is possible all the reads were not matched."
            fi
            exit 1
        fi

    >>>
    runtime {
        cpuPlatform: "Intel Skylake"
        cpu: "~{cpu}"
        preemptible: preemptible_tries_final
        memory: "~{memory_gb} GiB"
        disks: "local-disk " + ceil(mapped_bam_size_local_ssd) + " LOCAL"
        docker: docker
        noAddress: false
    }
    output {
        File monitoring_log = "monitoring.log"
        Array[File?] sorted_cram = glob("cram/*.cram")
        Array[File?] sorted_cram_index = glob("cram/*.cram.crai")
        Array[File] sorter_stats_csv = glob("~{sorter_output_path}/**/*.csv")
        Array[File] sorter_stats_json = glob("~{sorter_output_path}/**/*.json")
        Array[File?] unmatched_cram = glob("unmatched/*_unmatched.cram")
        Array[File?] unmatched_sorter_stats_csv = glob("unmatched/*_unmatched.csv")
        Array[File?] unmatched_sorter_stats_json = glob("unmatched/*_unmatched.json")
        Array[File?] sorter_out_bedgraph_mapq0 = glob("~{sorter_output_path}/**/*_0.bedGraph.gz")
        Array[File?] sorter_out_bedgraph_mapq1 = glob("~{sorter_output_path}/**/*_1.bedGraph.gz")
        Array[File?] sorter_out_bedgraph_extents = glob("~{sorter_output_path}/**/*_0.bedGraph.gz.extents.tsv")
        Array[File?] fastq_files = glob("fastq/*") # in case output format is fastq (like in single cell)
        Array[File?] sub_sampled_output = glob("~{sorter_output_path}/**/*_sample.*")
    }
}

task ConvertToFastq {
    input {
        File input_file
        String base_file_name
        String? output_dir
        String? demux_extra_args
        String? samtools_extra_args
        String output_format = "fastq"  # cram/sam/fastq supported, recommended not to change in this task
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
        memory: "16 GiB"
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