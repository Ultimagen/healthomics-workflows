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

task FindInsertBarcodeFastq {
    input {
        Array[File] input_fastq_list
        Array[File] sub_sumple_fastq_list
        Array[File] sorter_csv_stats_list
        Array[File] sorter_json_stats_list
        String base_file_name
        String insert_rg
        String barcode_rg
        String? additional_rg

        File monitoring_script

        String docker
        Boolean no_address
    }

    Int disk_size = round(1.5*size(input_fastq_list, "GB") + 3*size(sub_sumple_fastq_list, "GB") + 3*size(sorter_csv_stats_list, "GB") + 3*size(sorter_json_stats_list, "GB") + 20)

    #String output_insert_fastq = "~{base_file_name}_~{insert_rg}.fastq.gz"
    #String output_barcode_fastq = "~{base_file_name}_~{barcode_rg}.fastq.gz"
    String additional_rg_override = "~{if defined(additional_rg) then additional_rg else ""}"
    String output_sub_sample_fastq = "~{base_file_name}_~{insert_rg}_sample.fastq.gz"
    String output_sorter_stats_csv = "~{base_file_name}_~{insert_rg}.csv"
    String output_sorter_stats_json = "~{base_file_name}_~{insert_rg}.json"

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Find the relevant barcode/insert files from the list of input fastq
        python <<CODE
        import sys
        import shutil

        input_fastq_list = "~{sep=',' input_fastq_list}".split(',')
        insert_fastq = None
        barcode_fastq = None
        additional_fastq = None
        # Writing a -1 into the additional fastq index file in case the additional read group is not defined
        with open("additional_fastq_index.txt", "w") as f:
                f.write(str(-1))

        for fastq_idx, fastq in enumerate(input_fastq_list):
            if "~{insert_rg}" in fastq:
                insert_fastq = fastq
                with open("insert_fastq_index.txt", "w") as f:
                    f.write(str(fastq_idx))
            elif "~{barcode_rg}" in fastq:
                barcode_fastq = fastq
                with open("barcode_fastq_index.txt", "w") as f:
                    f.write(str(fastq_idx))
            elif "~{additional_rg_override}" != "" and "~{additional_rg_override}" in fastq:
                additional_fastq = fastq
                with open("additional_fastq_index.txt", "w") as f:
                    f.write(str(fastq_idx))

        sub_sumple_fastq_list = "~{sep=',' sub_sumple_fastq_list}".split(",")
        for fastq in sub_sumple_fastq_list:
            if "~{insert_rg}" in fastq:
                sub_sample_fastq = fastq

        sorter_csv_stats_list = "~{sep=',' sorter_csv_stats_list}".split(",")
        for csv in sorter_csv_stats_list:
            if "~{insert_rg}" in csv:
                sorter_stats_csv = csv
        
        sorter_json_stats_list = "~{sep=',' sorter_json_stats_list}".split(",")
        for json_file in sorter_json_stats_list:
            if "~{insert_rg}" in json_file:
                sorter_stats_json = json_file

        try:
            print(f"{insert_fastq}\n{barcode_fastq}\n{sub_sample_fastq}\n{sorter_stats_csv}\n{sorter_stats_json}")
        except NameError as e:
            print(f"ERROR: couldn't find one of the input files. Missing file error: {e}", file=sys.stderr)
            print(f"ERROR: Given files are:{input_fastq_list},{sub_sumple_fastq_list},{sorter_csv_stats_list},{sorter_stats_json}", file=sys.stderr)
            raise e

        # copy files to new location
        print(f"Copying files to new location")
        shutil.copy(sub_sample_fastq, "~{output_sub_sample_fastq}")
        shutil.copy(sorter_stats_csv, "~{output_sorter_stats_csv}")
        shutil.copy(sorter_stats_json, "~{output_sorter_stats_json}")

        print("[DEBUG] Insert fastq:",insert_fastq)
        print("[DEBUG] Barcode fastq:", barcode_fastq)
        if additional_fastq is not None:
            print("[DEBUG] Additional fastq:", additional_fastq)
        print("[DEBUG] Additional fastq: ", additional_fastq)
        print("[DEBUG] Sub sample fastq: ~{output_sub_sample_fastq}")
        print("[DEBUG] Sorter stats csv: ~{output_sorter_stats_csv}")
        print("[DEBUG] Sorter stats json: ~{output_sorter_stats_json}")
        CODE
    >>>

    runtime{
        cpu: 1
        memory: "4 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        Int insert_fastq_index            = read_int("insert_fastq_index.txt")
        Int barcode_fastq_index           = read_int("barcode_fastq_index.txt")
        Int additional_fastq_index = read_int("additional_fastq_index.txt")
        File insert_sub_sample_fastq = "~{output_sub_sample_fastq}"
        File insert_sorter_stats_csv = "~{output_sorter_stats_csv}"
        File insert_sorter_stats_json = "~{output_sorter_stats_json}"
    }
}

task SingleCellQc {
    input{
        File trimmer_stats
        File trimmer_failure_codes
        File sorter_stats_csv
        File sorter_stats_json
        File star_stats
        File star_reads_per_gene
        File insert_sub_sample_fastq
        String base_file_name
        SingleCellQcThresholds qc_thresholds
        String star_db

        File monitoring_script

        Int memory_gb
        Int preemptible_tries
        Boolean no_address
        String docker

    }

    Int disk_size = round(3*size(trimmer_stats, "GB") + 3*size(trimmer_failure_codes, "GB") + 3*size(sorter_stats_csv, "GB") + 3*size(star_stats, "GB") + 3*size(star_reads_per_gene, "GB") + 3*size(insert_sub_sample_fastq, "GB") + 20)
    
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo "[DEBUG $(date)] Running single cell qc script"
        mkdir -p sc_qc

        single_cell_qc \
            --sample-name ~{base_file_name} \
            --trimmer-stats ~{trimmer_stats} \
            --trimmer-failure-codes ~{trimmer_failure_codes} \
            --sorter-stats ~{sorter_stats_csv} \
            --sorter-stats-json ~{sorter_stats_json} \
            --star-stats ~{star_stats} \
            --star-reads-per-gene ~{star_reads_per_gene} \
            --insert ~{insert_sub_sample_fastq} \
            --output-path sc_qc \
            --pass-trim-rate ~{qc_thresholds.pass_trim_rate} \
            --read-length ~{qc_thresholds.read_length} \
            --fraction-below-read-length ~{qc_thresholds.fraction_below_read_length} \
            --percent-aligned ~{qc_thresholds.percent_aligned} \
            --star-db ~{star_db}

      convert_h5_to_json \
            --root_element "metrics" \
            --ignored_h5_key_substring histogram \
            --input_h5 sc_qc/~{base_file_name}.scRNA.applicationQC.h5 \
            --output_json sc_qc/~{base_file_name}.scRNA.applicationQC.json
    >>>

    runtime{
        preemptible: preemptible_tries
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }
    output{
        File report         = "sc_qc/~{base_file_name}.scRNA.applicationQC.html"
        File h5             = "sc_qc/~{base_file_name}.scRNA.applicationQC.h5"
        File json           = "sc_qc/~{base_file_name}.scRNA.applicationQC.json"
        File monitoring_log = "monitoring.log"
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

        echo "~{base_file_name} ~{umi_len} ~{default="" library_direction} ~{barcode_file} ~{star_log_file}"
        echo "Params" > ~{base_file_name}.stats.txt
        echo "sample: ~{base_file_name}"
        echo "UMI-length: ~{umi_len}" >> ~{base_file_name}.stats.txt
        echo "Library-direction: ~{default="" library_direction}" >> ~{base_file_name}.stats.txt

        #echo "mkdir output_stats"
        #mkdir output_stats
        #cp  ~{gene_features_stats} ~{gene_summary_csv} ~{gene_umi_per_cell_sorted} output_stats/
        #cp ~{star_log_file} ~{barcode_file} ~{base_file_name}.stats.txt output_stats/
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

