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
#
# CHANGELOG in reverse chronological order
#   - Include tasks for running demux (to create fastq file with Ultima headers) and 
#     CreateSyntheticPairedEnd (to convert demux fastq output into splited fastq files).
#   - Add new task to run STAR solo (based on bash and not python).


import "structs.wdl" as Structs

task CreateSyntheticPairedEnd {
    input {
        File input_fastq
        String base_file_name
        String? barcode_fastq_file_suffix
        String? insert_fastq_file_suffix
        String? barcode_fastq_header_suffix
        String? insert_fastq_header_suffix

        Int cpu
        File monitoring_script
        Int disk_size = ceil(3 * size(input_fastq, "GB") + 20)
        Int preemptible_tries
        String docker
    }

    String barcode_fastq_file = "~{base_file_name}" + if defined(barcode_fastq_file_suffix) then "~{barcode_fastq_file_suffix}" else "_barcode.fastq.gz"
    String insert_fastq_file = "~{base_file_name}" + if defined(insert_fastq_file_suffix) then "~{insert_fastq_file_suffix}" else "_insert.fastq.gz"
    String barcode_fastq_suffix = if defined(barcode_fastq_header_suffix) then " ~{barcode_fastq_header_suffix}" else ""
    String insert_fastq_suffix = if defined(insert_fastq_header_suffix) then " ~{insert_fastq_header_suffix}" else ""

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 & 
    
        zcat ~{input_fastq} | \
        tee >( awk '{if (NR % 4 == 1) {print substr($0,0,length($0)-length($15))"~{insert_fastq_suffix}"} else {print}}' FS=: |\
             pigz > ~{insert_fastq_file} ) | \
            awk 'NR % 4 == 1 {print substr($0,0,length($0)-length($15))"~{barcode_fastq_suffix}""\n"$15"\n+";\
             for (i = 0; i < length($15); i++) {printf "I"}; printf "\n"}' FS=: |\
             pigz > ~{barcode_fastq_file}

    >>>
    runtime {
        preemptible: preemptible_tries
        cpu: "~{cpu}"
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: false
    }
    output {
        File output_barcodes_fastq  = "~{barcode_fastq_file}"
        File output_insert_fastq    = "~{insert_fastq_file}"
        File monitoring_log         = "monitoring.log"
    }
}

# new task to run STAR solo (based on bash and not python)
task StarSolo {
  input {
    File insert_fastq
    File barcode_fastq
    String base_file_name
    StarSoloParams star_solo_params

    Int cpu
    File monitoring_script
    Int disk_size = ceil(3*size(insert_fastq, "GB") +3*size(barcode_fastq, "GB") + 3*size(star_solo_params.genome, "GB") + 20 )
    Int preemptible_tries
    Boolean no_address
    String docker
    }
    String genome_dir = "genome_dir"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir ~{genome_dir}
        unzip ~{star_solo_params.genome} -d ~{genome_dir}

        STAR \
            --readFilesIn ~{insert_fastq} ~{barcode_fastq} \
            --readFilesCommand zcat \
            --genomeDir ~{genome_dir} \
            --runThreadN ~{cpu} \
            --soloType CB_UMI_Simple \
            --soloFeatures Gene \
            --soloCellFilter ~{star_solo_params.cell_filter} \
            --soloCBwhitelist ~{star_solo_params.barcode_whitelist} \
            --soloCBstart 1 \
            --soloCBlen ~{star_solo_params.cell_barcode_length} \
            --soloUMIstart ~{star_solo_params.cell_barcode_length+1} \
            --soloUMIlen ~{star_solo_params.umi_length} \
            --soloBarcodeReadLength ~{star_solo_params.cell_barcode_length + star_solo_params.umi_length} \
            --soloStrand ~{star_solo_params.strand} \
            --outSAMattributes ~{sep=" " star_solo_params.out_sam_attributes} \
            ~{"--sjdbGTFfile " + star_solo_params.gtf_override} \
            ~{default="" star_solo_params.extra_args} \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix ~{base_file_name}. 
    
        mkdir -p "~{base_file_name}.Solo.out/Gene/filtered/"
        touch "~{base_file_name}.Solo.out/Gene/filtered/barcodes.tsv"
        touch "~{base_file_name}.Solo.out/Gene/filtered/features.tsv"
        touch "~{base_file_name}.Solo.out/Gene/filtered/matrix.mtx"
    >>>

    runtime {
        preemptible: "~{preemptible_tries}"
        cpu: "~{cpu}"
        memory: "64 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File star_log_file = "~{base_file_name}.Log.final.out"
        File star_log_params_file = "~{base_file_name}.Log.out"
        File barcode_file = '~{base_file_name}.Solo.out/Barcodes.stats'
        File output_bam = "~{base_file_name}.Aligned.out.bam"
        File gene_features_stats = "~{base_file_name}.Solo.out/Gene/Features.stats"
        File gene_summary_csv = "~{base_file_name}.Solo.out/Gene/Summary.csv"
        File gene_umi_per_cell_sorted = "~{base_file_name}.Solo.out/Gene/UMIperCellSorted.txt"
        File gene_filtered_features = "~{base_file_name}.Solo.out/Gene/filtered/features.tsv"
        File gene_filtered_barcodes = "~{base_file_name}.Solo.out/Gene/filtered/barcodes.tsv"
        File gene_filtered_matrix = "~{base_file_name}.Solo.out/Gene/filtered/matrix.mtx"
        File monitoring_log = "monitoring.log"
    }
}

# Combine all the statistics from all single cell workflow steps into a single file
task CombineStatistics {
    input {
        File? trimmer_stats
        File? star_solo_stats_csv
        File? star_stats_csv
        String base_file_name

        Int cpu
        File monitoring_script
        Int? disk_size
        Int preemptible_tries
        Boolean no_address
        String docker
    }
    String output_combined_statistics = "~{base_file_name}_combined_statistics.csv"
    Int disk_size_override = select_first([disk_size,  20])
    
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        conda activate genomics.py3
        
        python <<CODE
        import pandas as pd
        
        columns=["Stat_Type", "Measure", "Value_Type", "Measure_Value"]
        combined_df = pd.DataFrame(columns=columns)

        # Trimmer stats
        if "~{default="None" trimmer_stats}" != "None":
            trimmer_df = pd.read_csv("~{trimmer_stats}")

            combined_df["Measure"] = trimmer_df.loc[trimmer_df["is required"], "segment label"] + " failure rate"
            combined_df["Stat_Type"] = "trimmer"
            combined_df["Value_Type"] = "Percentage"
            combined_df["Measure_Value"] = 100 * trimmer_df.loc[trimmer_df["is required"], "num failures"] / trimmer_df.loc[trimmer_df["is required"], "num input reads"]

            precent_trimmed_read = 100 * trimmer_df["num trimmed reads"][0] / trimmer_df["num input reads"][0]
            tmp_df = pd.DataFrame({
                "Stat_Type": "trimmer",
                "Measure": "percent trimmed read",
                "Value_Type": "Percentage",
                "Measure_Value": [precent_trimmed_read]})
            combined_df = combined_df.append(tmp_df, ignore_index=True)

            tmp_df = pd.DataFrame({
                "Stat_Type": "trimmer",
                "Measure": ["num input reads", "num trimmed reads"],
                "Value_Type": "Value",
                "Measure_Value": [trimmer_df["num input reads"][0], trimmer_df["num trimmed reads"][0]]})
            combined_df = combined_df.append(tmp_df, ignore_index=True)

        # STAR solo stats
        if "~{default="None" star_solo_stats_csv}" != "None":
            star_solo_df = pd.read_csv("~{star_solo_stats_csv}")
            combined_df = combined_df.append(star_solo_df, ignore_index=True)

        # STAR stats
        if "~{default="None" star_stats_csv}" != "None":
            star_df = pd.read_csv("~{star_stats_csv}")

            star_df['metric'] = star_df['metric_type']+ ' - ' + star_df['metric']
            star_df.drop('metric_type', axis=1, inplace=True)
            star_df.rename(columns={'metric': 'Measure', 'value': 'Measure_Value', 'value_type': 'Value_Type'}, inplace=True)
            star_df['Value_Type'] = star_df['Value_Type'].replace({'percentage': 'Percentage', 'datetime': 'Date', 'number': 'Value'})

            star_df['Stat_Type'] =  'STAR'

            combined_df = combined_df.append(star_df, ignore_index=True)

        combined_df.to_csv("~{output_combined_statistics}", index=False)
        
        CODE
    >>>
    runtime {
        preemptible: preemptible_tries
        cpu: "~{cpu}"
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size_override) + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File combined_statistics  = "~{output_combined_statistics}"
        File monitoring_log       = "monitoring.log"
    }
}

task SeparateReads {
    input {
        # Script arguments params
        File single_cell_fastq_or_cram
        String base_file_name
        Int cdna_trimming_length
        Int num_last_bases_to_mask_umi
        Int num_first_cdna_bases_to_clip
        Int min_cdna_len
        Int quality_cutoff
        Int umi_len
        Int umi_quality_threshold
        Int cbc_len
        String ilmn_suffix = "_S1_L001"
        String text_to_prepend_to_fastq_ids
        Boolean cdna_reversed
        String? extra_args
        Adapters10x split_rna_adapters

        # For using some other references sequences, must be taken from gs://gcp-public-data--broad-references (or https://github.com/broadinstitute/gatk/tree/master/src/test/resources/large)
        References references

        String docker
        # Technical issues params
        RuntimeParams separate_reads_runtime_params
        File monitoring_script
    }

    String detect_input_ending = if sub(single_cell_fastq_or_cram, ".*\\.cram$", "is_cram") == "is_cram" then "is_cram" else "is_bam"
    String utility_command = 
        if (detect_input_ending == "is_cram") then
            "SingleCellPipelineTool"
        else
            "SingleCellPipelineFastqTool"
    
    String reverse_complement_read2 = if(cdna_reversed) then "true" else "false"
    # Set some command parameters with default values
    String reference_param =
        if (utility_command == "SingleCellPipelineTool") then
            "--reference " + references.ref_fasta
        else # utility_command is "SingleCellPipelineFastqTool"
            ""
    String extra_args_string = select_first([extra_args, ""])
    # Set some runtime parameters with values given from template or use default values
    String gitc_path = select_first([separate_reads_runtime_params.gitc_path_override, "/usr/gitc/"])
    Int preemptibles = select_first([separate_reads_runtime_params.preemptible_tries_override, 1])
    Boolean no_address = select_first([separate_reads_runtime_params.no_address_override, true ])
    String disk_type = select_first([separate_reads_runtime_params.disk_type_override, "HDD"])
    Int memory_gb = select_first([separate_reads_runtime_params.memory_gb_override, 15])
    Int cpu_num = select_first([separate_reads_runtime_params.cpu_num_override, 1])
    Int calc_disk_size = ceil(size(single_cell_fastq_or_cram, "GB") * 2)
    Int disk_size_min_gb = ceil(select_first([separate_reads_runtime_params.disk_size_gb_override, 50]))
    Int disk_size = if calc_disk_size > disk_size_min_gb then calc_disk_size else disk_size_min_gb
    Int max_retries = select_first([separate_reads_runtime_params.max_retries_override, 1])

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        start=$(date +%s)
        echo "export GATK_STACKTRACE_ON_USER_EXCEPTION=true"
        export GATK_STACKTRACE_ON_USER_EXCEPTION=true

        echo "calling docker ~{docker} with memory ~{memory_gb}GB, \
        disk ~{disk_size} ~{disk_type} and cpus num ~{cpu_num} no address ~{no_address}"


        echo "java -cp ~{gitc_path}GATK_ultima.jar:~{gitc_path}single_cell_pipeline.jar -Xmx5g \
        org.broadinstitute.hellbender.Main ~{utility_command} --input ~{single_cell_fastq_or_cram} \
        ~{reference_param} --multiproc true --base-filename ~{base_file_name}~{ilmn_suffix} \
        --adapter-5p-override ~{split_rna_adapters.adapter_5p},~{split_rna_adapters.adapter_min_error_rate_5p},~{split_rna_adapters.adapter_min_overlap_5p} \
        --adapter-3p-override ~{split_rna_adapters.adapter_3p},~{split_rna_adapters.adapter_min_error_rate_3p},~{split_rna_adapters.adapter_min_overlap_3p} \
        --adapter-middle-override ~{split_rna_adapters.adapter_middle},~{split_rna_adapters.adapter_min_error_rate_middle},~{split_rna_adapters.adapter_min_overlap_middle} \
        --cdna-trimming-length ~{cdna_trimming_length} --compressed-output true \
        --cbc-umi-mask-last-bytes ~{num_last_bases_to_mask_umi} --cdna-first-bases-to-clip ~{num_first_cdna_bases_to_clip} \
        --min-cdna-length ~{min_cdna_len} --quality-cutoff ~{quality_cutoff} \
        --umi-length ~{umi_len} --umi-quality-threshold ~{umi_quality_threshold} \
        --cbc-length ~{cbc_len} --reverse-complement-read2 ~{reverse_complement_read2} \
        --output-id-prefix "~{text_to_prepend_to_fastq_ids}" ~{extra_args_string}"


        java -cp ~{gitc_path}GATK_ultima.jar:~{gitc_path}single_cell_pipeline.jar -Xmx5g \
        org.broadinstitute.hellbender.Main ~{utility_command} --input ~{single_cell_fastq_or_cram} \
        ~{reference_param} --multiproc true --base-filename ~{base_file_name}~{ilmn_suffix} \
        --adapter-5p-override ~{split_rna_adapters.adapter_5p},~{split_rna_adapters.adapter_min_error_rate_5p},~{split_rna_adapters.adapter_min_overlap_5p} \
        --adapter-3p-override ~{split_rna_adapters.adapter_3p},~{split_rna_adapters.adapter_min_error_rate_3p},~{split_rna_adapters.adapter_min_overlap_3p} \
        --adapter-middle-override ~{split_rna_adapters.adapter_middle},~{split_rna_adapters.adapter_min_error_rate_middle},~{split_rna_adapters.adapter_min_overlap_middle} \
        --cdna-trimming-length ~{cdna_trimming_length} --compressed-output true \
        --cbc-umi-mask-last-bytes ~{num_last_bases_to_mask_umi} --cdna-first-bases-to-clip ~{num_first_cdna_bases_to_clip} \
        --min-cdna-length ~{min_cdna_len} --quality-cutoff ~{quality_cutoff} \
        --umi-length ~{umi_len} --umi-quality-threshold ~{umi_quality_threshold} \
        --cbc-length ~{cbc_len} --reverse-complement-read2 ~{reverse_complement_read2} \
        --output-id-prefix "~{text_to_prepend_to_fastq_ids}" ~{extra_args_string}

        mv ~{base_file_name}~{ilmn_suffix}_report.json ~{base_file_name}_report.json

        end=$(date +%s)
        mins_elapsed=$(( ($end - $start) / 60))
        secs_elapsed=$(( ($end - $start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"
    >>>

    runtime {
        docker: docker
        preemptible: preemptibles
        cpu: cpu_num
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size + " " + disk_type
        maxRetries: max_retries
        noAddress: no_address
    }

    output {
        File monitoring_log = "monitoring.log"
        File cbcumi_fastq = "~{base_file_name}~{ilmn_suffix}_R1_001.fastq.gz"
        File cdna_fastq = "~{base_file_name}~{ilmn_suffix}_R2_001.fastq.gz"
        File report_file = "~{base_file_name}_report.json"
    }
}

task STARSolo {
    input {
        # STARSolo params
        File cbcumi_fastq_input
        File cdna_fastq_input

        String base_file_name
        File genome
        File whitelist

        String library_direction

        Int umi_len
        Int cbc_len
        Boolean cdna_reversed

        Int? starsolo_clip3p_nbases_override

        Float? starsolo_match_nmin_override
        String? starsolo_align_ends_type_override
        String? starsolo_multi_mapper_override
        String? starsolo_solo_umi_filtering_override
        String? starsolo_extra_args

        StarsoloScores starsolo_scores
        StarsoloBamParams starsolo_bam

        # If empty, uses the gtf used when creating STAR db:
        File? starsolo_gtf

        String docker
        # Technical issues params
        RuntimeParams starsolo_runtime_params

        File monitoring_script
    }

    # TODO: Maybe this can be part of starsolo docker?  # override with defaul value
    Int cbc_umi_len = umi_len + cbc_len

    Float starsolo_match_nmin = select_first([starsolo_match_nmin_override, 0.66])
    Int starsolo_clip3p_nbases = select_first([starsolo_clip3p_nbases_override, 0])

    Float starsolo_score_del_open = select_first([starsolo_scores.starsolo_score_del_open_override, -2])
    Float starsolo_score_ins_open = select_first([starsolo_scores.starsolo_score_ins_open_override, -2])
    Float starsolo_score_del_base = select_first([starsolo_scores.starsolo_score_del_base_override, -2])
    Float starsolo_score_ins_base = select_first([starsolo_scores.starsolo_score_ins_base_override, -2])

    String  starsolo_align_ends_type = select_first([starsolo_align_ends_type_override, 'Local'])
    String  starsolo_multi_mapper = select_first([starsolo_multi_mapper_override, 'Unique'])
    String  starsolo_solo_umi_filtering = select_first([starsolo_solo_umi_filtering_override, '-'])

    Boolean starsolo_save_bam = starsolo_bam.save_bam
    Boolean starsolo_sort_bam = if (defined(starsolo_bam)) then
                                select_first([starsolo_bam.sort_bam_override, false])
                                else false
    Boolean starsolo_include_unmapped_in_bam = if (defined(starsolo_bam)) then
                                               select_first([starsolo_bam.include_unmapped_override, false])
                                               else false
    String starsolo_limit_bam_sort_ram = if (defined(starsolo_bam)) then
                                           select_first([starsolo_bam.limit_sort_ram_override, ""])
                                           else ""
    Boolean starsolo_transcript_bam = if (defined(starsolo_bam)) then
                                      select_first([starsolo_bam.transcript_override, false])
                                      else false
    String starsolo_extra_args_string = select_first([starsolo_extra_args, ""])

    # Set some parameters with values given in template or use default values
    Int preemptibles = select_first([starsolo_runtime_params.preemptible_tries_override, 1])
    Boolean no_address = select_first([starsolo_runtime_params.no_address_override, true ])
    String disk_type = select_first([starsolo_runtime_params.disk_type_override, "HDD"])
    # TODO create formulas for cpus and memory based on input size
    Int cpu_num = select_first([starsolo_runtime_params.cpu_num_override, 2])
    Int disk_size_minimum = ceil(select_first([starsolo_runtime_params.disk_size_gb_override, 1]))
    # factor 3 because output bam is ~2 as large as the input fastq files, and the genome is compressed
    Int disk_size_formula = ceil(size(cbcumi_fastq_input, "GB") + size(cdna_fastq_input, "GB") + size(genome, "GB")) * 6 + 10
    Int disk_size =  if disk_size_minimum > disk_size_formula then disk_size_minimum else disk_size_formula
    Int max_retries = select_first([starsolo_runtime_params.max_retries_override, 1])
    Int memory_gb = select_first([starsolo_runtime_params.memory_gb_override, 120])

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        start=$(date +%s)

        echo "calling docker ~{docker} with memory ~{memory_gb} GiB, disk ~{disk_size} ~{disk_type} and cpus num ~{cpu_num} no address ~{no_address}"

        export TMPDIR=tmp
        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        mkdir result

        python <<CODE
        import os
        from subprocess import check_call

        call_args = ['STAR', '--soloType', 'CB_UMI_Simple', '--genomeDir', 'genome_ref', '--runThreadN', '~{cpu_num}', '--soloFeatures', 'Gene', '--soloCellFilter','EmptyDrops_CR']
        call_args.extend(['--soloBarcodeReadLength','~{cbc_umi_len}'])

        # mate 1 is the cdna, mate 2 is the cbcumi.  When giving two values, there was an error that there is just one mate
        call_args.extend(['--clip3pNbases','~{starsolo_clip3p_nbases}'])

        # This value is 0.66 by default (as the starsolo defaults)
        call_args.extend(['--outFilterMatchNminOverLread', '~{starsolo_match_nmin}', '--outFilterScoreMinOverLread', '~{starsolo_match_nmin}'])

        # These two options shouold be uncommented once starsolo version is updated:
        #call_args.extend(['--soloMultiMappers', '~{starsolo_multi_mapper}'])
        #call_args.extend(['--soloUMIfiltering', '~{starsolo_solo_umi_filtering}'])
        call_args.extend(['--scoreDelOpen', '~{starsolo_score_del_open}', '--scoreDelBase', '~{starsolo_score_del_base}'])
        call_args.extend(['--scoreInsOpen', '~{starsolo_score_ins_open}', '--scoreInsBase', '~{starsolo_score_ins_base}'])
        call_args.extend(['--alignEndsType', '~{starsolo_align_ends_type}'])

        if '~{starsolo_save_bam}' == 'true':
            if '~{starsolo_sort_bam}' == 'true':
                call_args.extend(['--outSAMtype', 'BAM', 'SortedByCoordinate'])
                call_args.extend(['--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'MD', 'jM', 'jI', 'XS', 'MC', 'ch', 'CR', 'UR', 'CY', 'UY', 'GX', 'GN', 'CB', 'UB'])
            else:
                # If not sorting, cannot include CB, UB fields
                call_args.extend(['--outSAMtype', 'BAM', 'Unsorted'])
                call_args.extend(['--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'MD', 'jM', 'jI', 'XS', 'MC', 'ch', 'CR', 'UR', 'CY', 'UY', 'GX', 'GN'])
            if '~{starsolo_include_unmapped_in_bam}' == 'true':
                call_args.extend(['--outSAMunmapped', 'Within'])
            if len('~{starsolo_limit_bam_sort_ram}') > 0:
                call_args.extend(['--limitBAMsortRAM', '~{starsolo_limit_bam_sort_ram}'])
            if '~{starsolo_transcript_bam}' == 'true':
                call_args.extend(['--quantMode', 'TranscriptomeSAM'])
        else:
            call_args.extend(['--outSAMtype', 'None'])

        call_args.extend(['--soloCBwhitelist', '~{whitelist}', '--soloCBstart', '1', '--soloCBlen', '~{cbc_len}', '--soloUMIstart', '~{cbc_len+1}', '--soloUMIlen', '~{umi_len}'])

        if '~{library_direction}' is 'three_prime':
            if '~{cdna_reversed}' == 'true':
               call_args.extend(['--soloStrand', 'Forward'])
            else:
               call_args.extend(['--soloStrand', 'Reverse'])
        elif '~{library_direction}' is 'five_prime':
            if '~{cdna_reversed}' == 'true':
               call_args.extend(['--soloStrand', 'Reverse'])
            else:
               call_args.extend(['--soloStrand', 'Forward'])



        call_args.extend(['--readFilesCommand', 'zcat'])
        call_args.extend(['--readFilesIn', '~{cdna_fastq_input}',  '~{cbcumi_fastq_input}'])
        call_args.extend(['--outFileNamePrefix', 'result/~{base_file_name}_'])

        if '~{starsolo_gtf}':
            call_args.extend(['--sjdbGTFfile', '~{starsolo_gtf}'])

        if len('~{starsolo_extra_args_string}') > 0:
            call_args.append('~{starsolo_extra_args_string}')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        end=$(date +%s)
        mins_elapsed=$(( ($end - $start) / 60))
        secs_elapsed=$(( ($end - $start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"
    >>>

    runtime {
        docker: docker
        preemptible: preemptibles
        cpu: cpu_num
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size + " " + "~{disk_type}"
        maxRetries: max_retries
        noAddress: no_address
    }

    output {
        File monitoring_log = "monitoring.log"
        File star_log_file = "result/~{base_file_name}_Log.final.out"
        File star_log_params_file = "result/~{base_file_name}_Log.out"
        File barcode_file = 'result/~{base_file_name}_Solo.out/Barcodes.stats'
        File gene_features_stats = "result/~{base_file_name}_Solo.out/Gene/Features.stats"
        File gene_summary_csv = "result/~{base_file_name}_Solo.out/Gene/Summary.csv"
        File gene_umi_per_cell_sorted = "result/~{base_file_name}_Solo.out/Gene/UMIperCellSorted.txt"
        File gene_filtered_features = "result/~{base_file_name}_Solo.out/Gene/filtered/features.tsv"
        File gene_filtered_barcodes = "result/~{base_file_name}_Solo.out/Gene/filtered/barcodes.tsv"
        File gene_filtered_matrix = "result/~{base_file_name}_Solo.out/Gene/filtered/matrix.mtx"
        # Array[File] gene_files = ['result/~{base_file_name}_Solo.out/Gene/Features.stats', 'result/~{base_file_name}_Solo.out/Gene/Summary.csv', 'result/~{base_file_name}_Solo.out/Gene/UMIperCellSorted.txt']
        # Array[File] gene_count_files = ['result/~{base_file_name}_Solo.out/Gene/filtered/features.tsv', 'result/~{base_file_name}_Solo.out/Gene/filtered/barcodes.tsv', 'result/~{base_file_name}_Solo.out/Gene/filtered/matrix.mtx']
        Array[File?] output_bams = glob("result/*.bam")
    }
}

task GatherStatistics {
    input {
        File gene_features_stats
        File gene_summary_csv
        File gene_umi_per_cell_sorted
        File star_log_file
        File barcode_file

        File? separate_reads_statistics
        String base_file_name
        String? library_direction
        Int umi_len
        File monitoring_script

        String docker
        # Technical issues params
        RuntimeParams statistics_runtime_params

    }

    Int preemptibles = select_first([statistics_runtime_params.preemptible_tries_override, 1])
    Boolean no_address = select_first([statistics_runtime_params.no_address_override, true ])
    String disk_type = select_first([statistics_runtime_params.disk_type_override, "HDD"])
    Int cpu_num = select_first([statistics_runtime_params.cpu_num_override, 2])
    Int disk_size = ceil(select_first([statistics_runtime_params.disk_size_gb_override, 1]))
    Int max_retries = select_first([statistics_runtime_params.max_retries_override, 1])

    command {
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        start=$(date +%s)

        echo "calling docker ~{docker} with memory 1 MiB, disk ~{disk_size} ~{disk_type} and cpus num ~{cpu_num} no address ~{no_address}"

        source ~/.bashrc
        conda activate genomics.py3

        echo "~{base_file_name} ~{umi_len} ~{default="" library_direction} ~{barcode_file} ~{star_log_file}"
        echo "Params" > ~{base_file_name}.stats.txt
        echo "sample: ~{base_file_name}"
        echo "UMI-length: ~{umi_len}" >> ~{base_file_name}.stats.txt
        echo "Library-direction: ~{default="" library_direction}" >> ~{base_file_name}.stats.txt

        echo "mkdir output_stats"
        mkdir output_stats
        cp  ~{gene_features_stats} ~{gene_summary_csv} ~{gene_umi_per_cell_sorted} output_stats/
        cp ~{star_log_file} ~{barcode_file} ~{base_file_name}.stats.txt output_stats/
        #tar -zcf ~{base_file_name}.tgz output_stats

        python <<CODE
        import json
        import pandas as pd
        import numpy as np
        from collections import defaultdict

        # Putting stats into single csv file
        print("read_csv stats txt")
        param_df = pd.read_csv("~{base_file_name}.stats.txt", header = None, sep = '\s+', skiprows = 1)
        param_df.columns = ['Measure', 'Value']
        param_df['Measure'] = param_df['Measure'].str.replace(':$', '')
        param_df['Stat_Type'] = 'Params'

        print("read_csv barcode file")
        barcode_df = pd.read_csv("~{barcode_file}", header = None, sep = '\s+')
        barcode_df.columns = ['Measure', 'Value']
        barcode_df['Stat_Type'] = 'Barcode'

        print("read_csv gene_files")
        feature_df = pd.read_csv("~{gene_features_stats}", header = None, sep = '\s+')
        feature_df.columns = ['Measure', 'Value']
        feature_df['Stat_Type'] = 'Feature'

        summary_df = pd.read_csv("~{gene_summary_csv}", header = None, sep = ',')
        summary_df.columns = ['Measure', 'Value']
        is_percentage = summary_df['Value'].astype('float') <= 1
        summary_df['Percentage'] = np.nan
        summary_df.loc[is_percentage, 'Percentage'] = summary_df.loc[is_percentage, 'Value'].multiply(100).astype(str) + '%'
        summary_df.loc[is_percentage, 'Value'] = np.nan
        summary_df['Stat_Type'] = 'Summary'

        print("read_csv star_log_file")
        alignment_df = pd.read_csv("~{star_log_file}", header = None, sep = '\t')
        alignment_df.columns = ['Measure', 'Value']
        alignment_df = alignment_df[alignment_df['Measure'].str.endswith("|")]
        alignment_df['Measure'] = alignment_df['Measure'].str.replace(' \|\s*', '')
        has_percentage = alignment_df['Value'].str.endswith('%')
        alignment_df['Percentage'] = np.nan
        alignment_df.loc[has_percentage, 'Percentage'] = alignment_df.loc[has_percentage, 'Value']
        alignment_df.loc[has_percentage, 'Value'] = np.nan
        alignment_df['Stat_Type'] = 'Alignment'

        this_json_dict = defaultdict(lambda:0)
        adapter_df = pd.DataFrame()
        if "~{default="None" separate_reads_statistics}" != "None":
            with open("~{separate_reads_statistics}", "r") as json_file:
                data = json.load(json_file)
                this_json_dict['Total read pairs processed'] += data["readsIn"]
                this_json_dict['Reads written (passing filters)'] += int(data["readsOut"].split()[0])
                this_json_dict['Reads filtered'] +=  int(data["filtered"].split()[0])
                this_json_dict['Reads with 5p adapter'] +=  int(data["adapter5p"].split()[0])
                this_json_dict['Reads with 3p adapter'] += int(data["adapter3p"].split()[0])
                this_json_dict['Reads with middle sequence'] += int(data["adapterMiddle"].split()[0])
                this_json_dict['Total base pairs processed'] += data["bpIn"]
                this_json_dict['Total base pairs output'] += int(data["bpOut"].split()[0])
                this_json_dict['Reads too short after quality trimming'] += int(data["trimmedTooShort"].split()[0])
                this_json_dict['Discarded base pairs, quality trimming'] += int(data["bpCutoff"].split()[0])
                this_json_dict['Reads with low umi quality'] += int(data["umiQualityDropped"].split()[0])
                this_json_dict['Read 2 too short'] += int(data["read2TooShortDropped"].split()[0])
                this_json_dict['Read 1 too short'] += int(data["read1TooShortDropped"].split()[0])

            adapter_df = pd.DataFrame.from_dict(this_json_dict, orient = 'index', columns = ['Value'])
            adapter_df['Stat_Type'] = 'Adapter removal and prefiltering'
            adapter_df['Percentage'] = np.nan

            # Calc percentage metrics
            read_columns_for_percentage = ['Reads written (passing filters)', 'Reads too short after quality trimming', 'Reads with low umi quality', 'Read 1 too short', 'Read 2 too short']
            adapter_df.loc[read_columns_for_percentage, 'Percentage'] = \
                100*adapter_df['Value'][read_columns_for_percentage] / adapter_df.loc['Total read pairs processed', 'Value']

            read_columns_for_percentage_of_filtered = ['Reads with 5p adapter', 'Reads with middle sequence', 'Reads with 3p adapter']
            adapter_df.loc[read_columns_for_percentage_of_filtered, 'Percentage'] = \
                100*adapter_df['Value'][read_columns_for_percentage_of_filtered] / adapter_df.loc['Reads filtered', 'Value']

            base_columns_for_percentage = ['Total base pairs output']
            adapter_df.loc[base_columns_for_percentage, 'Percentage'] = \
                100*adapter_df['Value'][base_columns_for_percentage] / adapter_df.loc['Total base pairs processed', 'Value']

            adapter_df.index.names = ['Measure']
            adapter_df = adapter_df.reset_index()
            adapter_df['Stat_Type'] = 'Adapter removal and prefiltering'

        # TODO: add cbcum adapter stats
        print("create combined_df")
        adapter_cbcumi_df = pd.DataFrame()
        combined_df = pd.concat([param_df, adapter_df, adapter_cbcumi_df, alignment_df, feature_df, barcode_df, summary_df]).reset_index()
        combined_df['Measure'] = combined_df['Measure'].str.strip()
        combined_df['index_orig'] = combined_df.index
        combined_df = combined_df[['index_orig', 'Stat_Type', 'Measure', 'Value', 'Percentage']]
        combined_df = combined_df.melt(id_vars=['index_orig', 'Stat_Type', 'Measure'], var_name='Value_Type', value_name='Measure_Value')
        combined_df = combined_df.dropna()
        looks_like_date = combined_df['Measure'].str.endswith(" on")
        combined_df.loc[looks_like_date, 'Value_Type'] = 'Date'
        combined_df = combined_df.sort_values(by=['index_orig'], ascending =True)
        combined_df.loc[combined_df["Measure_Value"].str.endswith("%", na=False), "Measure_Value"] = combined_df.loc[combined_df["Measure_Value"].str.endswith("%", na=False), "Measure_Value"].str[:-1].astype(float)
        print("to_csv stats summary")
        combined_df.to_csv("~{base_file_name}.stat_summary.csv", index = False)

        CODE

        end=$(date +%s)
        mins_elapsed=$(( ($end - $start) / 60))
        secs_elapsed=$(( ($end - $start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"
    }

    runtime {
        docker: docker
        memory: "1 MiB"
        disks: "local-disk " + ceil(disk_size) + " " + disk_type
        cpu: cpu_num
        preemptible: preemptibles
        noAddress: no_address
        maxRetries: max_retries
    }
    output {
        File stats_csv="~{base_file_name}.stat_summary.csv"
        File monitoring_log="monitoring.log"
    }
}

# create html report
task CreateReport {
   input {
       # call arguments
       String base_file_name
       # python notebook args
       File input_csv_file
       String docker
       # Technical issues params
       RuntimeParams report_runtime_params
       File monitoring_script
   }
   String notebook_file_in = "/VariantCalling/ugvc/reports/single_cell_qc_report.ipynb"
   Int preemptibles = select_first([report_runtime_params.preemptible_tries_override, 1])
   Boolean no_address = select_first([report_runtime_params.no_address_override, true ])
   String disk_type = select_first([report_runtime_params.disk_type_override, "HDD"])
   Int memory_gb = select_first([report_runtime_params.memory_gb_override, 4])
   Int cpu_num = select_first([report_runtime_params.cpu_num_override, 1])
   Int disk_size = ceil(select_first([report_runtime_params.disk_size_gb_override, 4]) + size(input_csv_file, "GB"))
   Int max_retries = select_first([report_runtime_params.max_retries_override, 1])

   command <<<
       set -eo pipefail
       start=$(date +%s)
       source ~/.bashrc
       conda activate genomics.py3

       bash ~{monitoring_script} | tee monitoring.log >&2 &

       cp ~{notebook_file_in} .
       cp ~{input_csv_file} input_for_html_report.csv

       basename_notebook=$(basename ~{notebook_file_in} ".ipynb")

       echo "jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute ${basename_notebook}.ipynb"
       jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute ${basename_notebook}.ipynb

       echo "jupyter nbconvert --to html ${basename_notebook}.nbconvert.ipynb --template classic --no-input --output ~{base_file_name}.single_cell_qc_report.html"
       jupyter nbconvert --to html ${basename_notebook}.nbconvert.ipynb --template classic --no-input --output ~{base_file_name}.single_cell_qc_report.html

       echo "**************** D O N E ****************"

       end=$(date +%s)
       mins_elapsed=$(( ($end - $start) / 60))
       secs_elapsed=$(( ($end - $start) % 60 ))
       if [ $secs_elapsed -lt 10 ]; then
           secs_elapsed=0$secs_elapsed
       fi
       echo "Total run time: $mins_elapsed:$secs_elapsed"

   >>>
   runtime {
       preemptible: preemptibles
       memory: "~{memory_gb} GB"
       disks: "local-disk " + disk_size + " " + disk_type
       docker: docker
       cpu: cpu_num
       noAddress: no_address
       maxRetries: max_retries
   }
   output {
       File monitoring_log = "monitoring.log"
       File report_html = "~{base_file_name}.single_cell_qc_report.html"
   }
}
