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


workflow ppmSeqPreprocess {
  input {
    # Workflow args
    String pipeline_version = "1.15.6" # !UnusedDeclaration

    # Data inputs
    Array[File] input_cram_bam_list
    Array[File] ref_fastas_cram
    String base_file_name
    String adapter_version
    String? ppmSeq_analysis_extra_args  # extra args for python ugvc ppmSeq_analysis
    String? sample_name  # for when merging cram files, if not provided base_file_name is used

    # References
    References references
    
    # trimmer parameters
    TrimmerParameters trimmer_parameters

    # alignment parameters
    UaReferences ua_parameters
    
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
    description: "This workflow takes untrimmed ppmSeq sequencing data, trims, aligns and sorts, and generates a QC report"
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
            "sort_cram",
            "sort_cram_index",
            "sort_stats_csv",
            "sort_stats_json"
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
    sample_name: {
      help: "Sample name used when merging cram files, if not provided base_file_name is used",
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
      type: "UaReferences",
      category: "param_required"
    }
    sorter_params: {
      help: "Sorter parameters",
      type: "SorterParams",
      category: "param_required"
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
    sort_stats_csv_out: {
      help: "Sorter stats csv output",
      type: "File",
      category: "output"
    }
    sort_stats_json_out: {
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
    report_html: {
      help: "ppmSeq QC report html",
      type: "File",
      category: "output"
    }
    aggregated_metrics_h5: {
      help: "ppmSeq QC aggregated metrics h5",
      type: "File",
      category: "output"
    }
    aggregated_metrics_json: {
      help: "ppmSeq QC aggregated metrics json",
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
        sample_name =         select_first([sample_name, base_file_name]),
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
  File sort_stats_csv_out = select_first([select_first([TrimAlignSort.sort_stats_csv])[0]]) # convert Array[File?]? to File
  File sort_stats_json_out = select_first([select_first([TrimAlignSort.sort_stats_json])[0]]) # convert Array[File?]? to File
  File output_cram_bam = TrimAlignSort.output_cram_bam
  File output_cram_bam_index = select_first([TrimAlignSort.output_cram_bam_index])

  call ppmSeqTasks.ppmSeqQC as ppmSeqQC {
    input:
        adapter_version =                     adapter_version,
        trimmer_histogram_csv =               trimmer_histogram_csv_out,
        trimmer_histogram_extra_csv =         trimmer_histogram_csv_extra_out,
        trimmer_failure_codes_csv =           trimmer_failure_codes,
        sorter_stats_csv =                    sort_stats_csv_out,
        sorter_stats_json =                   sort_stats_json_out,
        base_file_name =                      base_file_name,
        ppmSeq_analysis_extra_args = ppmSeq_analysis_extra_args,
        docker =                              global.ugbio_ppmseq_docker,
        monitoring_script =                   global.monitoring_script, # !FileCoercion
  }

output {
        File trimmer_histogram          = trimmer_histogram_csv_out
        File? trimmer_histogram_extra   = trimmer_histogram_csv_extra_out
        File trimmer_stats              = select_first([TrimAlignSort.trimmer_stats])
        File trimmer_failure_codes_csv  = trimmer_failure_codes

        File sort_cram                  = output_cram_bam
        File sort_cram_index            = output_cram_bam_index
        File sort_stats_csv             = sort_stats_csv_out
        File sort_stats_json            = sort_stats_json_out

        File report_html                = ppmSeqQC.report_html 
        File aggregated_metrics_h5      = ppmSeqQC.aggregated_metrics_h5
        File aggregated_metrics_json    = ppmSeqQC.aggregated_metrics_json
    }
}