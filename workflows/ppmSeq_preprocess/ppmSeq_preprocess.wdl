version 1.0

# LICENSE
#   Copyright 2023 Ultima Genomics
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
# This wdl generates a QC report based on preprocessing of ppmSeq sequencing data
# 
# CHANGELOG
# 1.12 updated ppmSeq version names 
# 1.12 updated workflow name
# 1.10.1 Added LA_v7 adapter version


import "trim_align_sort.wdl" as TrimAlignSortSubWF
import "tasks/ppmSeq_preprocess.wdl" as ppmSeqTasks
import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks

workflow ppmSeqPreprocess {
  input {
    # Workflow args
    String pipeline_version = "1.27.3" # !UnusedDeclaration

    # Data inputs
    Array[File] input_cram_bam_list
    Array[File] ref_fastas_cram
    String base_file_name
    String adapter_version
    String? ppmSeq_analysis_extra_args  # extra args for python ugvc ppmSeq_analysis

    # References
    References references
    
    # trimmer parameters
    TrimmerParameters trimmer_parameters

    # alignment parameters
    UaParameters ua_parameters
    
    # sorter parameters (for marking duplicates)
    SorterParams sorter_params

    TrimAlignSortSteps steps

    # general parameters
    Boolean no_address = true
    Int preemptible_tries = 1
    Int cpu

    # Used for running on other clouds (aws)
    File? monitoring_script_input
    String? cloud_provider_override

    Boolean create_md5_checksum_outputs = false

    # winval validations
    # base_file_name
    #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)

    # adapter version 
    #@wv adapter_version in {"v1", "legacy_v5", "legacy_v5_start", "legacy_v5_end", "dmbl"}

    # references
    #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
    #@wv suffix(references['ref_dict']) == '.dict'
    #@wv suffix(references['ref_fasta_index']) == '.fai'
    #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']

    ## Trimmer checks
    #@wv 'trim' in steps and steps['trim'] -> defined(trimmer_parameters)
    #@wv defined(trimmer_parameters) and not 'formats_description' in trimmer_parameters -> 'format' in trimmer_parameters

    ## UA
    #@wv 'align' in steps and steps['align'] -> defined(ua_parameters)
    #@wv 'align' in steps and steps['align'] and 'ua_index' in ua_parameters -> suffix(ua_parameters['ua_index']) == '.uai'
    #@wv 'align' in steps and steps['align'] -> suffix(ua_parameters['ref_alt']) == '.alt'

    #@wv 'align' in steps and steps['align'] 
    #@wv 'trim' in steps and steps['trim'] 
    #@wv 'sort' in steps and steps['sort'] 
  }

  meta {
    description: "The ppmSeq Preprocess pipeline is designed to process untrimmed ppmSeq sequencing data through a complete preprocessing workflow including adapter trimming, alignment to reference genome, sorting, duplicate marking, and QC report generation. This analysis is generally done on the UG sequencer, this pipeline is intended for cases where it did not happen or was improperly configured. For more details on ppmSeq, see https://www.ultimagenomics.com/products/ppmseq-tm/ and https://www.biorxiv.org/content/10.1101/2025.08.11.669689v1. \n\nThe following input templates are available for different input data: \n\n1) `ppmSeq_preprocess_template-ppmSeq.json` | Use this template for ppmSeq data. The input CRAM file should NOT be trimmed. If it contains the ppmSeq tags (e.g. st, et), it was trimmed. \n\n2) `ppmSeq_preprocess_template-ppmSeq_legacy_v5.json` | Use this template for LEGACY v5 ppmSeq data. This is an older version of the ppmSeq adapters, generally not available since 2024. The input CRAM file should NOT be trimmed. If it contains the ppmSeq tags (e.g. as, ts), it was trimmed. \n\n3) `ppmSeq_preprocess_template-ppmSeq_post_native_adapter_trimming.json` | Use this template only for the case where the UG native adapters were trimmed, but not the ppmSeq adapters and loop. This generally happens if the application_type is configured to be 'native' instead of 'ppmSeq', can be verified by the presence of an 'a3' tag in some reads but the absence of ppmSeq tags (e.g. st, et). "
    author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "cloud_provider_override",
            "no_address",
            "preemptible_tries",
            "monitoring_script_input",
            "cpu",
            "ppmSeqQC.cpu",
            "ppmSeqQC.memory_gb",
            "ppmSeqQC.preemptible_tries",
            "ppmSeqQC.disk_size",
            "steps",
            "Globals.glob",
            "trimmer_histogram",
            "trimmer_histogram_extra",
            "trimmer_stats",
            "trimmer_failure_codes_csv",
            "Demux.mapq_override",
            "MergeMd5sToJson.output_json"
        ]}
  }

  parameter_meta {
    base_file_name: {
      help: "Base file name for output files.",
      type: "String", 
      category: "input_required"
    }
    input_cram_bam_list: {
      help: "Input CRAM or BAM file list",
      type: "Array[File]",
      category: "input_required"
    }
    ref_fastas_cram: {
      help: "Reference fasta file for cache tarball",
      type: "Array[File]",
      category: "ref_required"
    }
    adapter_version: {
      help: "ppmSeq adapter version",
      type: "String",
      category: "input_required"
    }
    ppmSeq_analysis_extra_args: {
      help: "Extra arguments for ppmSeq analysis",
      type: "String",
      category: "input_optional"
    }
    references: {
      help: "Reference files",
      type: "References",
      category: "ref_required"
    }
    trimmer_parameters: {
      help: "Trimmer parameters",
      type: "TrimmerParameters",
      category: "param_required"
    }
    ua_parameters: {
      help: "UA alignment parameters",
      type: "UaParameters",
      category: "param_required"
    }
    sorter_params: {
      help: "Sorter parameters",
      type: "SorterParams",
      category: "param_required"
    }
    create_md5_checksum_outputs: {
       help: "Create md5 checksum for requested output files",
       type: "Boolean",
       category: "input_optional"
    }
    trimmer_histogram_csv_out: {
      help: "Trimmer histogram csv output",
      type: "File",
      category: "output"
    }
    trimmer_histogram_csv_extra_out: {
      help: "Trimmer histogram csv extra output",
      type: "File",
      category: "output"
    }
    trimmer_failure_codes: {
      help: "Trimmer failure codes csv",
      type: "File",
      category: "output"
    }
    sorter_stats_csv: {
      help: "Sorter stats csv output",
      type: "File",
      category: "output"
    }
    sorter_stats_json: {
      help: "Sorter stats json output",
      type: "File",
      category: "output"
    }
    output_cram_bam: {
      help: "Output CRAM or BAM file, trimmed aligned and sorted",
      type: "File",
      category: "output"
    }
    output_cram_bam_index: {
      help: "Output CRAM or BAM index file",
      type: "File",
      category: "output"
    }
    unmatched_cram: {
      help: "Unmatched cram file output from sorter (if defined).",
      type: "File",
      category: "output"
    }
    unmatched_sorter_stats_csv: {
      help: "Unmatched output cram files sorter stats in csv format (if defined).",
      type: "File",
      category: "output"
    }
    unmatched_sorter_stats_json: {
        help: "Unmatched output cram files sorter stats in json format (if defined).",
        type: "File",
        category: "output"
    }
    bedgraph_mapq0: {
      help: "Bedgraph mapq0 file.",
      type: "File",
      category: "output"
    }
    bedgraph_mapq1: {
      help: "Bedgraph mapq1 file.",
      type: "File",
      category: "output"
    }
    report_html: {
      help: "ppmSeq QC report html",
      type: "File",
      category: "output"
    }
    application_qc_h5: {
      help: "ppmSeq QC aggregated metrics h5",
      type: "File",
      category: "output"
    }
    aggregated_metrics_json: {
      help: "ppmSeq QC aggregated metrics json",
      type: "File",
      category: "output"
    }
    md5_checksums_json: {
        help: "json file that will contain md5 checksums for requested output files",
        type: "File",
        category: "output"
    }

  }

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers

  call TrimAlignSortSubWF.TrimAlignSort as TrimAlignSort {
    input:
        input_cram_bam_list = input_cram_bam_list,
        base_file_name =      base_file_name,
        steps =               steps,
        references =          references,
        ref_fastas_cram =     ref_fastas_cram,
        trimmer_parameters =  trimmer_parameters,
        aligner =             "ua",
        ua_parameters =       ua_parameters,
        sorter_params =       sorter_params,
        no_address =          no_address,
        preemptible_tries =   preemptible_tries,
        cpu =                 cpu,
        monitoring_script_input = monitoring_script_input
  }
  File trimmer_histogram_csv_out = select_first([select_first([TrimAlignSort.trimmer_histogram])[0]])  # convert Array[File?]? to File
  Array[File?] trimmer_histogram_csv_extra_out_arr = select_first([TrimAlignSort.trimmer_histogram_extra])  # convert Array[File?]? to Array[File]
  if (defined(trimmer_histogram_csv_extra_out_arr) && (length(select_all(trimmer_histogram_csv_extra_out_arr)) > 0)){
    File trimmer_histogram_csv_extra_out = select_all(trimmer_histogram_csv_extra_out_arr)[0]  # convert Array[File?]? to File
  }
  File trimmer_failure_codes = select_first([TrimAlignSort.trimmer_failure_codes_csv])
  File sorter_stats_csv_out = select_first([TrimAlignSort.sorter_stats_csv])
  File sorter_stats_json_out = select_first([TrimAlignSort.sorter_stats_json])

  call ppmSeqTasks.ppmSeqQC as ppmSeqQC {
    input:
        adapter_version =                     adapter_version,
        trimmer_histogram_csv =               trimmer_histogram_csv_out,
        trimmer_histogram_extra_csv =         trimmer_histogram_csv_extra_out,
        trimmer_failure_codes_csv =           trimmer_failure_codes,
        sorter_stats_csv =                    sorter_stats_csv_out,
        sorter_stats_json =                   sorter_stats_json_out,
        base_file_name =                      base_file_name,
        ppmSeq_analysis_extra_args = ppmSeq_analysis_extra_args,
        docker =                              global.ugbio_ppmseq_docker,
        monitoring_script =                   global.monitoring_script, # !FileCoercion
  }

    File output_cram_bam_                 = select_first([TrimAlignSort.output_cram_bam])
    File output_cram_bam_index_           = select_first([TrimAlignSort.output_cram_bam_index])
    File report_html_                     = ppmSeqQC.report_html
    if (create_md5_checksum_outputs) {

        Array[File] output_files = select_all(flatten([
                                                      select_first([[output_cram_bam_], []]),
                                                      select_first([[output_cram_bam_index_], []]),
                                                      select_first([[report_html_], []]),
                                                      ]))

        scatter (file in output_files) {
            call UGGeneralTasks.ComputeMd5 as compute_md5 {
                input:
                    input_file = file,
                    docker = global.ubuntu_docker,
            }
        }

        call UGGeneralTasks.MergeMd5sToJson {
            input:
                md5_files = compute_md5.checksum,
                docker = global.ugbio_core_docker
        }
    }
output {
        File trimmer_histogram          = trimmer_histogram_csv_out
        File? trimmer_histogram_extra   = trimmer_histogram_csv_extra_out
        File trimmer_stats              = select_first([TrimAlignSort.trimmer_stats])
        File trimmer_failure_codes_csv  = trimmer_failure_codes

        File output_cram_bam                = output_cram_bam_
        File output_cram_bam_index          = output_cram_bam_index_
        File sorter_stats_csv               = sorter_stats_csv_out
        File sorter_stats_json              = sorter_stats_json_out
        File? unmatched_cram                = TrimAlignSort.unmatched_cram
        File? unmatched_sorter_stats_csv    = TrimAlignSort.unmatched_sorter_stats_csv
        File? unmatched_sorter_stats_json   = TrimAlignSort.unmatched_sorter_stats_json
        File? bedgraph_mapq0                = TrimAlignSort.bedgraph_mapq0
        File? bedgraph_mapq1                = TrimAlignSort.bedgraph_mapq1

        File report_html                = report_html_
        File application_qc_h5          = ppmSeqQC.aggregated_metrics_h5
        File aggregated_metrics_json    = ppmSeqQC.aggregated_metrics_json

        File? md5_checksums_json = MergeMd5sToJson.md5_json
    }
}