version 1.0

import "structs.wdl" as Structs


task ppmSeqQC {
  input {
    String adapter_version
    File trimmer_histogram_csv
    File? trimmer_histogram_extra_csv
    File trimmer_failure_codes_csv
    File sorter_stats_csv
    File? sorter_stats_json
    String base_file_name
    String? ppmSeq_analysis_extra_args
    File monitoring_script
    Float memory_gb = 8
    Int cpu = 1
    Int preemptible_tries = 1
    String docker
    Float disk_size=10
  }
  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail

    ppmSeq_qc_analysis \
      --adapter-version ~{adapter_version} \
      --trimmer-histogram-csv ~{trimmer_histogram_csv} \
      ~{true="--trimmer-histogram-extra-csv " false="" defined(trimmer_histogram_extra_csv)}~{trimmer_histogram_extra_csv} \
      --trimmer-failure-codes-csv "~{trimmer_failure_codes_csv}" \
      --sorter-stats-csv ~{sorter_stats_csv} \
      --sorter-stats-json ~{sorter_stats_json} \
      --output-path output/ \
      --output-basename ~{base_file_name} \
      ~{default="" ppmSeq_analysis_extra_args}
  >>>
  runtime {
    preemptible: preemptible_tries
    cpu: "~{cpu}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    noAddress: false
  }
  output {
    File monitoring_log = "monitoring.log"
    File aggregated_metrics_h5 = "output/~{base_file_name}.ppmSeq.applicationQC.h5"
    File aggregated_metrics_json = "output/~{base_file_name}.ppmSeq.applicationQC.json"
    File report_html = "output/~{base_file_name}.ppmSeq.applicationQC.html"
  }
}