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
#   QC pipeline tasks
# CHANGELOG
#

import "structs.wdl"

task CheckContamination {
    input {
        File input_bam
        File input_bam_index
        File contamination_sites_ud
        File contamination_sites_bed
        File contamination_sites_mu
        References references
        String output_prefix
        String docker
        String exec_path
        Int preemptible_tries
        File monitoring_script
        Boolean disable_sanity_check = false

        # since the task and the docker are very specific, we set the sizes in the task
        Int disk_size_gb = ceil(if ceil((size(input_bam, "GB")) +
                            (size(references.ref_fasta, "GB") +
                            size(references.ref_fasta_index, "GB") +
                            size(references.ref_dict, "GB")) + 80) > 510 then ceil((size(input_bam, "GB")) +
                            (size(references.ref_fasta, "GB") +
                            size(references.ref_fasta_index, "GB") +
                            size(references.ref_dict, "GB")) + 80) else 510)
        Int cpu = 1
        Int memory_mb = 2000
        Int max_retries = 1
  }

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        ~{exec_path}/VerifyBamID \
        --Verbose \
        --NumPC 4 \
        --Output ~{output_prefix} \
        --BamFile ~{input_bam} \
        --Reference ~{references.ref_fasta} \
        --UDPath ~{contamination_sites_ud} \
        --MeanPath ~{contamination_sites_mu} \
        --BedPath ~{contamination_sites_bed} \
        --adjust-MQ 0 \
        ~{true="--DisableSanityCheck " false="" disable_sanity_check} \
        1>/dev/null

        # used to read from the selfSM file and calculate contamination, which gets printed out
        python3 <<CODE
        import csv
        import sys
        with open('~{output_prefix}.selfSM') as selfSM:
            reader = csv.DictReader(selfSM, delimiter='\t')
            i = 0
            for row in reader:
                if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
                    # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
                    # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
                    # vcf and bam.
                    sys.stderr.write("Found zero likelihoods. A zero value for the likelihoods implies no data. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
                    sys.exit(1)
                with open('contamination.txt', 'w') as f:
                    print(float(row["FREEMIX"]), file=f)

                with open('coverage.txt', 'w') as f2:
                    print(float(row["AVG_DP"]), file = f2)

                i = i + 1
                # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
                # and the results are not reliable.
                if i != 1:
                    sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
                    sys.exit(2)
        CODE
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File monitoring_log = "monitoring.log"
        File selfSM = "~{output_prefix}.selfSM"
        Float contamination = read_float("contamination.txt")
        Float contamination_coverage = read_float("coverage.txt")
    }
}



task CollectDuplicateMetrics {
    input {
        File input_bam
        References references

        File monitoring_script
        String metrics_filename
        String docker
        Boolean no_address
        Int preemptible_tries
        Int disk_size = ceil(size(input_bam, "GB") + size(references.ref_fasta, "GB") + 20)
        String gitc_path = "/usr/gitc"
    }

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        java -Xms3000m -jar ~{gitc_path}/picard.jar CollectDuplicateMetrics \
        -I ~{input_bam} \
        -R ~{references.ref_fasta} \
        -M ~{metrics_filename}
    >>>

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 1
        memory: "4 GB"
        preemptible: preemptible_tries
        docker: docker
        noAddress: no_address
    }

    output {
        File monitoring_log = "monitoring.log"
        File duplicate_metrics = "~{metrics_filename}"
    }
}


# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
    input {
        File monitoring_script
        File input_bam
        References references
        String metrics_filename
        Int preemptible_tries
        String docker
        Boolean no_address
        Int disk_size = ceil(size(input_bam, "GB") + size(references.ref_fasta, "GB") + 20)
        String gitc_path = "/usr/gitc"
    }

    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms3000m -jar ~{gitc_path}/picard.jar \
            CollectQualityYieldMetrics \
            INPUT=~{input_bam} \
            R=~{references.ref_fasta} \
            OQ=true \
            FLOW_MODE=true \
            OUTPUT=~{metrics_filename}
    >>>

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "4 GB"
        preemptible: preemptible_tries
        docker: docker
        noAddress: no_address
    }
    output {
        File metrics = "~{metrics_filename}"
        File monitoring_log = "monitoring.log"
    }
}


# Note these tasks will break if the read lengths in the bam are greater than 250.
task CollectWgsMetrics {
    input {
        File monitoring_script
        File input_bam
        File input_bam_index
        String metrics_filename
        References references
        File wgs_coverage_interval_list
        Int? read_length
        Int preemptible_tries
        String docker
        Boolean no_address
        Int disk_size_gb = if ceil((size(input_bam, "GB")) +
                                  size(references.ref_fasta, "GB") +
                                  size(references.ref_fasta_index, "GB") +
                                  size(references.ref_dict, "GB") + 160) > 510 then ceil((size(input_bam, "GB")) +
                                  size(references.ref_fasta, "GB") +
                                  size(references.ref_fasta_index, "GB") +
                                  size(references.ref_dict, "GB") + 160) else 510

        String gitc_path = "/usr/gitc"

    }
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms8000m -jar ~{gitc_path}/picard.jar \
            CollectWgsMetrics \
            INPUT=~{input_bam} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE=~{references.ref_fasta} \
            INCLUDE_BQ_HISTOGRAM=true \
            INTERVALS=~{wgs_coverage_interval_list} \
            OUTPUT=~{metrics_filename} \
            USE_FAST_ALGORITHM=false \
            COUNT_UNPAIRED=true \
            COVERAGE_CAP=12500 \
            READ_LENGTH=~{default=250 read_length}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "10 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File metrics = "~{metrics_filename}"
        File monitoring_log = "monitoring.log"
    }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
    input {
        File monitoring_script
        File input_bam
        File input_bam_index
        String metrics_filename
        References references
        File wgs_coverage_interval_list
        Int? read_length
        Int preemptible_tries
        String docker
        Boolean no_address
        Int disk_size_gb = if ceil((size(input_bam, "GB")) +
                                  size(references.ref_fasta, "GB") +
                                  size(references.ref_fasta_index, "GB") +
                                  size(references.ref_dict, "GB") + 160) > 510 then ceil((size(input_bam, "GB")) +
                                  size(references.ref_fasta, "GB") +
                                  size(references.ref_fasta_index, "GB") +
                                  size(references.ref_dict, "GB") + 160) else 510

        String gitc_path = "/usr/gitc"

    }
    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms8000m -jar ~{gitc_path}/picard.jar \
            CollectRawWgsMetrics \
            INPUT=~{input_bam} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE=~{references.ref_fasta} \
            INCLUDE_BQ_HISTOGRAM=true \
            INTERVALS=~{wgs_coverage_interval_list} \
            OUTPUT=~{metrics_filename} \
            USE_FAST_ALGORITHM=false \
            COUNT_UNPAIRED=true \
            COVERAGE_CAP=12500 \
            READ_LENGTH=~{default="250" ""+read_length}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "10 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File metrics = "~{metrics_filename}"
        File monitoring_log = "monitoring.log"
    }
}

# Add flag for enabling adding extra argument for methylation data (bisulfite)
# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
    input {
        File monitoring_script
        File input_bam
        File input_bam_index
        String output_bam_prefix
        References references
        Int preemptible_tries
        String docker
        String ug_adapter
        Boolean no_address
        Boolean is_methylation_flag

        String gitc_path = "/usr/gitc"
        Int disk_size = if ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) > 510 then ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) else 510

    }
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms5000m -jar ~{gitc_path}/picard.jar \
            CollectMultipleMetrics \
            INPUT=~{input_bam} \
            REFERENCE_SEQUENCE=~{references.ref_fasta} \
            OUTPUT=~{output_bam_prefix} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            EXTRA_ARGUMENT="CollectAlignmentSummaryMetrics::ADAPTER_SEQUENCE=~{ug_adapter}" \
            EXTRA_ARGUMENT="CollectAlignmentSummaryMetrics::IS_BISULFITE_SEQUENCED"=~{is_methylation_flag} \
            PROGRAM="CollectGcBiasMetrics" \
            EXTRA_ARGUMENT="CollectGcBiasMetrics::IS_BISULFITE_SEQUENCED"=~{is_methylation_flag} \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY"
    >>>
    runtime {
        memory: "7 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        docker: docker
        noAddress: no_address
    }
    output {
        File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
        File? alignment_summary_pdf = "~{output_bam_prefix}.read_length_histogram.pdf"
        File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
        File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
        File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
        File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
        File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
        File monitoring_log = "monitoring.log"
    }
}

task CollectIntervalCoverages {
    input  {
        File input_cram_bam
        File input_cram_bam_index
        References references
        Int min_mapq
        Array[String]? region
        String docker
        File monitoring_script
        Int preemptible_tries
        # SSD is requested in batches of 375 GB
        Int disk_size = ceil((size(input_cram_bam, "GB")+120)/375) * 375
    }
    Array[String] region_str = select_first([region,[]])
    String region_argument = if defined(region) then '-r' else ''
    
    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        set -eo pipefail

        # shellcheck source=/dev/null
        source ~/.bashrc
        conda activate genomics.py3

        start=$(date +%s)
        COVERAGE_ANALYSIS="coverage_analysis.py collect_coverage"

        OUTPUT=coverage
        mkdir $OUTPUT
        echo "Collecting $OUTPUT"
        $COVERAGE_ANALYSIS \
            -i ~{input_cram_bam} \
            -o $OUTPUT \
            -Q ~{min_mapq} \
            --reference ~{references.ref_fasta} \
            ~{region_argument} ~{sep=" " region_str} 


        end=$(date +%s)
        mins_elapsed=$(( (end - start) / 60))
        secs_elapsed=$(( (end - start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
          secs_elapsed=0$secs_elapsed
        fi
        echo "Run time: $mins_elapsed:$secs_elapsed"
    >>>

    runtime {
        preemptible: preemptible_tries
        cpu: 4
        memory: "8 GB"
        disks: "local-disk " + disk_size + " LOCAL"
        docker: docker
        noAddress: false # no_address=true caused problems attaching SSD
    }

    output{
        Array[File] coverage_bw_w1 = glob("coverage/*w1.*depth.bw")
        Array[File] coverage_bw = glob("coverage/*.depth.bw")
        File monitoring_log = "monitoring.log"
    }
}


task CollectIntervalCoverageStats {
    input  {
        File input_cram_bam
        File input_cram_bam_index
        File coverage_intervals
        References references
        File reference_gaps
        File centromeres
        Int min_mapq
        String docker
        File monitoring_script
        Int preemptible_tries
        # SSD is requested in batches of 375 GB
        Int coverage_stats_disk = ceil((size(input_cram_bam, "GB")+120)/375) * 375
    }

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        set -eo pipefail

        # shellcheck source=/dev/null
        source ~/.bashrc
        conda activate genomics.py3

        download_start=$(date +%s)

        download_dest=$(dirname "$(realpath ~{coverage_intervals})")
        top_level_dir=$(echo "$download_dest" | awk -F/ '{print FS $2}')

        cloud_src="${download_dest/$top_level_dir/gs:\/}"

        # shellcheck disable=SC2001
        src=$(echo "$cloud_src" | sed 's~\/*$~/~')

        echo "Download from $src to $download_dest"
        gsutil -m cp -r "$src*" "$download_dest"

        download_end=$(date +%s)
        mins_elapsed=$(( (download_end - download_start) / 60))
        secs_elapsed=$(( (download_end - download_start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
          secs_elapsed=0$secs_elapsed
        fi
        echo "Download time: $mins_elapsed:$secs_elapsed"

        start=$(date +%s)
        COVERAGE_ANALYSIS="coverage_analysis.py "
        mkdir /data  # kind of a patch because of the localization issue

        OUTPUT=coverage
        mkdir $OUTPUT
        echo "Collecting $OUTPUT"
        $COVERAGE_ANALYSIS full_analysis \
            -i ~{input_cram_bam} \
            --coverage_intervals ~{coverage_intervals} \
            -o $OUTPUT \
            -Q ~{min_mapq} \
            --reference ~{references.ref_fasta} \
            --reference-gaps ~{reference_gaps} \
            --centromeres ~{centromeres}

        end=$(date +%s)
        mins_elapsed=$(( (end - start) / 60))
        secs_elapsed=$(( (end - start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
          secs_elapsed=0$secs_elapsed
        fi
        echo "Run time: $mins_elapsed:$secs_elapsed"
    >>>

    runtime {
        preemptible: preemptible_tries
        cpu: "16"
        memory: "16 GB"
        disks: "local-disk " + coverage_stats_disk + " LOCAL"
        docker: docker
        noAddress: false # no_address=true caused problems attaching SSD
    }

    output{
        Array[File] coverage_bw_w1 = glob("coverage/*w1.*depth.bw")
        Array[File] coverage_bw = glob("coverage/*.depth.bw")
        Array[File] coverage_parquet = glob("coverage/*.depth.parquet")
        File coverage_stats = glob("coverage/*.coverage_stats.q0.Q~{min_mapq}.l0.h5")[0]
        Array[File] coverage_boxplot = glob("coverage/*coverage_boxplot.q0.Q~{min_mapq}.l0.png")
        Array[File] coverage_profileplots = glob("coverage/*coverage_profile.w*.q0.Q~{min_mapq}.l0.png")
        File monitoring_log = "monitoring.log"
    }
}

task CheckPreValidation {
    input {
        File duplication_metrics
        File chimerism_metrics
        Float max_duplication_in_reasonable_sample
        Float max_chimerism_in_reasonable_sample
        Int preemptible_tries = 3
    }

    command <<<
        set -o pipefail
        set -e

        grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
        grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

        python3 <<CODE

        import csv
        with open('duplication.csv') as dupfile:
          reader = csv.DictReader(dupfile, delimiter='\t')
          for row in reader:
            with open("duplication_value.txt","w") as file:
              file.write(row['PERCENT_DUPLICATION'])
              file.close()

        with open('chimerism.csv') as chimfile:
          reader = csv.DictReader(chimfile, delimiter='\t')
          for row in reader:
            with open("chimerism_value.txt","w") as file:
              file.write(row['PCT_CHIMERAS'])
              file.close()

        CODE

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
        preemptible: preemptible_tries
        memory: "2 GiB"
    }
    output {
        Float duplication_rate = read_float("duplication_value.txt")
        Float chimerism_rate = read_float("chimerism_value.txt")
        Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
    }
}

task FastQC {
  input {
     Array[File] input_fastq
     String? adapter_5p
     File? limits
     String docker
     # Technical issues params
     RuntimeParams fastqc_runtime_params
     File monitoring_script
  }
     Int memory_gb = select_first([fastqc_runtime_params.memory_gb_override, 8])
     String disk_type = select_first([fastqc_runtime_params.disk_type_override, "HDD"])
     Int disk_size = ceil(select_first([fastqc_runtime_params.disk_size_gb_override, 2 * size(input_fastq, "GB")]))
     Boolean no_address = select_first([fastqc_runtime_params.no_address_override, false])
     Int cpu = select_first([fastqc_runtime_params.cpu_num_override, length(input_fastq)])

  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    mkdir fastqc_out
    fastqc -v

    if ~{defined(adapter_5p)} ; then
        echo adapter_5p$'\t'~{select_first([adapter_5p, ''])} > adapter.txt
        fastqc ~{sep=" " input_fastq} -a adapter.txt -o fastqc_out -t ~{cpu} ~{"--limits " + limits}
    else
        fastqc ~{sep=" " input_fastq} -o fastqc_out -t ~{cpu} ~{"--limits " + limits}
    fi

  >>>
  output {
    Array[File] reports_html = glob("fastqc_out/*.html")
    File monitoring_log="monitoring.log"
  }
  runtime
  {
    docker: docker
    cpu: "~{cpu}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + disk_size + " " + disk_type
    noAddress: no_address
  }
}


task CreateReportSingleSampleQC {
   input {
        # call arguments
        File monitoring_script
        String base_file_name

        # runtime arguments
        Float disk_size
        Int preemptible_tries
        String docker

        # python notebook args
        File input_h5_file
        String notebook_file_in= "/VariantCalling/ugvc/reports/single_sample_qc_create_html_report.ipynb"
        String top_metrics_file = "/VariantCalling/ugvc/reports/top_metrics_for_tbl.csv"

    }
    command <<<
        set -eo pipefail
        start=$(date +%s)

        # shellcheck source=/dev/null
        source ~/.bashrc
        conda activate genomics.py3
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        basename_notebook=$(basename ~{notebook_file_in} ".ipynb")

        papermill ~{notebook_file_in} "${basename_notebook}.papermill.ipynb" \
            -p top_metrics_file ~{top_metrics_file} \
            -p input_h5_file ~{input_h5_file} \
            -p input_base_file_name ~{base_file_name}

        jupyter nbconvert --to html "${basename_notebook}.papermill.ipynb" --template classic --no-input --output ~{base_file_name}.html

        echo "**************** D O N E ****************"

        end=$(date +%s)
        mins_elapsed=$(( (end - start) / 60))
        secs_elapsed=$(( (end - start) % 60 ))
        if [ $secs_elapsed -lt 10 ]; then
            secs_elapsed=0$secs_elapsed
        fi
        echo "Total run time: $mins_elapsed:$secs_elapsed"


    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "3 GB"
        docker: docker
        cpu: "1"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        noAddress: true
        maxRetries: preemptible_tries
    }

    output {
        File monitoring_log = "monitoring.log"
        File report_html = "~{base_file_name}.html"
    }
}


task CollectHsMetrics {
    input {
        File monitoring_script
        File input_bam
        File input_bam_index
        File target_exome_interval_list
        File bait_exome_interval_list
        String output_bam_prefix
        References references
        Int preemptible_tries
        String docker
        Boolean no_address
        String gitc_path = "/usr/gitc"
        Int disk_size = if ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) > 510 then ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) else 510

    }
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms4000m -jar ~{gitc_path}/picard.jar \
            CollectHsMetrics \
            I=~{input_bam} \
            R=~{references.ref_fasta} \
            O=~{output_bam_prefix}.HSmetrics \
            BAIT_INTERVALS=~{bait_exome_interval_list} \
            TARGET_INTERVALS=~{target_exome_interval_list}


        ls
    >>>
    runtime {
        memory: "10 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        docker: docker
        noAddress: no_address
    }
    output {
        File hybrid_selection_metrics = "~{output_bam_prefix}.HSmetrics"
        File monitoring_log = "monitoring.log"
    }
}


task TargetsBedcov {
    input {
        File monitoring_script
        File input_bam
        File input_bam_index
        File target_interval_list
        String base_file_name
        References references
        String? samtools_bedcov_extra_args
        Int preemptible_tries
        String docker
        Boolean no_address
        String gitc_path = "/usr/gitc"
        Int disk_size = if ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) > 510 then ceil((size(input_bam, "GB")) +
                               size(references.ref_fasta, "GB") +
                               size(references.ref_fasta_index, "GB") +
                               size(references.ref_dict, "GB") + 160) else 510

    }
    String target_intervals_basename = basename(target_interval_list, ".interval_list")
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms3000m -jar ~{gitc_path}/picard.jar \
            IntervalListToBed \
            I=~{target_interval_list} \
            O=~{target_intervals_basename}.targets.bed

        samtools bedcov ~{if defined(samtools_bedcov_extra_args) then samtools_bedcov_extra_args else ""} --reference ~{references.ref_fasta} ~{target_intervals_basename}.targets.bed ~{input_bam} | gzip > ~{base_file_name}.~{target_intervals_basename}.bedcov.gz
        ls
    >>>
    runtime {
        memory: "10 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        docker: docker
        noAddress: no_address
    }
    output {
        File targets_bedcov = "~{base_file_name}.~{target_intervals_basename}.bedcov.gz"
        File monitoring_log = "monitoring.log"
    }
}
