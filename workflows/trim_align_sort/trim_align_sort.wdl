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
#   This workflow allows to perform trimming, alignemnt and sorting to mark duplicates, togheter or as seperate steps.
# CHANGELOG
# 1.10.3 Works now when the input is realigned to a different reference genome
# 1.9.2 Added STAR use-case and migrated to omics
# 1.9.0 Aggregates sorter metrics in QC report
# 1.6.0 More trimmer formats [BIOIN-1105]
# 1.5.0 Initial commit

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/alignment_tasks.wdl" as UGAlignment
import "tasks/trimming_tasks.wdl" as TrimmingTasks
import "ua_align.wdl" as UaAlignWorkflow
import "ua_meth_align.wdl" as UaMethAlignWorkflow
import "star_align_gene_count.wdl" as StarAlignWorkflow
import "tasks/sorting_tasks.wdl" as SortTasks
import "tasks/qc_tasks.wdl" as QCTasks

workflow TrimAlignSort {
    input {
        String pipeline_version = "1.11" # !UnusedDeclaration
        Array[File] input_cram_bam_list
        Array[File] ref_fastas_cram
        String base_file_name
        TrimAlignSortSteps steps
        References references

        # merged cram parameters
        String? sample_name

        # trimmer parameters
        TrimmerParameters? trimmer_parameters

        # alignment parameters
        String? aligner # ua, ua-meth, star
        UaReferences? ua_parameters
        UaMethReferences? ua_meth_parameters

        ## STAR param
        File? star_genome
        StarGenomeGenerateParams? star_genome_generate_params
        String? star_align_extra_args
        File? star_align_gtf_override

        # sorter parameters (for marking duplicates)
        SorterParams? sorter_params

        # fastqc parameters        
        File? fastqc_limits

        # general parameters
        Boolean no_address = true
        Int preemptible_tries = 1
        Int cpu

        # Used for running on other clouds (aws)
        File? monitoring_script_input
        String dummy_input_for_call_caching = ""

    
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv len(input_cram_bam_list) > 0
        #@wv suffix(input_cram_bam_list) <= {".bam", ".cram"}
        #@wv 'trim' in steps or 'align' in steps or 'sort' in steps or 'fastqc' in steps -> (steps['trim'] or steps['align'] or steps['sort'] or steps['fastqc'])
        #@wv defined(references) -> len(references) == 3
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa','.fna'}
        #@wv suffix(references['ref_dict']) == '.dict'
        #@wv suffix(references['ref_fasta_index']) == '.fai'
        #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
        #@wv prefix(references['ref_dict']) == prefix(references['ref_fasta'])
        
        ## Trimmer checks
        #@wv 'trim' in steps and steps['trim'] -> defined(trimmer_parameters)
        #@wv defined(trimmer_parameters) and not 'formats_description' in trimmer_parameters -> 'format' in trimmer_parameters
        ##@wv 'formats_description' in trimmer_parameters -> not('local_formats_description' in trimmer_parameters)
        ##@wv 'local_formats_description' in trimmer_parameters -> not('formats_description' in trimmer_parameters)

        ## Align checks
        #@wv 'align' in steps and steps['align'] -> defined(aligner) and aligner in {"ua", "ua-meth", "star"}
        ## UA
        #@wv 'align' in steps and steps['align'] and aligner == "ua" -> defined(ua_parameters)
        #@wv 'align' in steps and steps['align'] and aligner == "ua" and 'ua_index' in ua_parameters -> suffix(ua_parameters['ua_index']) == '.uai'
        #@wv 'align' in steps and steps['align'] and aligner == "ua" -> suffix(ua_parameters['ref_alt']) == '.alt'
        ## UA-meth
        #@wv 'align' in steps and steps['align'] and aligner == "ua-meth" and defined(ua_meth_parameters) -> suffix(ua_meth_parameters['index_g2a']) == '.g2a'
        #@wv 'align' in steps and steps['align'] and aligner == "ua-meth" and defined(ua_meth_parameters) -> suffix(ua_meth_parameters['index_c2t']) == '.c2t'
        ## STAR
        #@wv 'align' in steps and steps['align'] and aligner == "star" -> defined(star_genome) or defined(star_genome_generate_params)
        #@wv 'align' in steps and steps['align'] and aligner == "star" and defined(star_genome) -> suffix(star_genome) == '.zip'
        #@wv 'align' in steps and steps['align'] and aligner == "star" and defined(star_genome_generate_params) -> suffix(star_genome_generate_params['fasta_files']) <= {'.fasta','.fa'}
        #@wv 'align' in steps and steps['align'] and aligner == "star" and defined(star_genome_generate_params) -> suffix(star_genome_generate_params['gtf_file']) == '.gtf'

        ## Sort checks
        #@wv 'sort' in steps and steps['sort'] -> defined(sorter_params)

        ## Fastqc checks - Fastqc cannot run alone
        #@wv 'fastqc' in steps and steps['fastqc'] -> (steps['trim'] or steps['align'] or steps['sort'])

    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    
    Boolean trim = select_first([steps.trim, false])
    Boolean align = select_first([steps.align, false])
    Boolean sort = select_first([steps.sort, false])
    Boolean fastqc = select_first([steps.fastqc, false])

    call UGAlignment.CreateReferenceCache {
        input:
            references = ref_fastas_cram,
            cache_populate_script = global.ref_cache_script, #!StringCoercion
            preemptible_tries = preemptible_tries,
            docker = global.perl_docker,
            dummy_input_for_call_caching = dummy_input_for_call_caching
    }

    if (length( input_cram_bam_list) > 1) {
        call UGGeneralTasks.MergeCramFiles {
            input: 
                output_base_name    = base_file_name,
                references          = references,
                cache_tarball       = CreateReferenceCache.cache_tarball,
                crams               = input_cram_bam_list,
                sample_name         = select_first([sample_name, "sm1"]),
                docker              = global.ug_vc_docker,
                monitoring_script   = monitoring_script,
                no_address          = no_address,
                cpus                = cpu,
                preemptible_tries   = preemptible_tries
        }
    }

    File single_cram_bam = select_first([MergeCramFiles.output_cram,  input_cram_bam_list[0]])

    if (trim) {
        call TrimmingTasks.Trimmer {
            input:
                input_cram_bam      = single_cram_bam,
                trimmer_parameters  = select_first([trimmer_parameters]),
                base_file_name      = base_file_name,
                docker              = global.trimmer_docker,
                preemptible_tries   = preemptible_tries,
                monitoring_script   = monitoring_script,  # !FileCoercion
                no_address          = no_address,
                cpus                = cpu
        }
        File trimmer_output_ucram = Trimmer.trimmed_ucram_list[0]
    }

    if (align) {
        File input_for_alignment = select_first([trimmer_output_ucram, single_cram_bam])
        String aligner_override  = select_first([aligner])

        if (aligner_override == "ua" ){
            call UaAlignWorkflow.UAAlignment {
                input:
                    input_files             = [input_for_alignment],
                    base_file_name          = base_file_name,
                    cache_tarball           = CreateReferenceCache.cache_tarball,
                    ua_parameters           = select_first([ua_parameters]),
                    references              = references,
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
            }
        }

        if (aligner_override == "ua-meth" ){
            call UaMethAlignWorkflow.UAMethAlignment {
                input:
                    input_files             = [input_for_alignment],
                    base_file_name          = base_file_name,
                    ua_meth_parameters      = select_first([ua_meth_parameters]),
                    references              = references,
                    cache_tarball           = CreateReferenceCache.cache_tarball,
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
            }
        }

        if (aligner_override == "star" ){
            call StarAlignWorkflow.StarAlignment {
                input:
                    genome                  = star_genome,
                    genome_generate_params  = star_genome_generate_params,
                    star_align_extra_args   = star_align_extra_args,
                    star_align_gtf_override = star_align_gtf_override,
                    input_bams              = [input_for_alignment],
                    base_file_name          = base_file_name + ".star.aln",
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
                    cpu                     = cpu
            }
        }
        File align_output = select_first([UAAlignment.ua_output_bam, UAMethAlignment.ua_output_bam ,StarAlignment.output_bam])
    }

    if (sort) {
        File input_for_sort = select_first([align_output, single_cram_bam])

        call SortTasks.Sorter {
            input:
                input_file         = input_for_sort,
                base_file_name     = base_file_name,
                reference_fasta    = references.ref_fasta,
                sorter_params      = select_first([sorter_params]),
                monitoring_script  = monitoring_script,  # !FileCoercion
                docker             = global.sorter_docker,
                preemptible_tries  = preemptible_tries,
                cpu                = cpu
        }
    }


    if (fastqc) {
        File input_for_fastqc = select_first([Sorter.sorted_cram, align_output, trimmer_output_ucram])
        call SortTasks.Demux as ConvertToFastq{
            input:
                input_file          = input_for_fastqc,
                base_file_name      = base_file_name,
                output_format       = "fastq",
                reference_fasta     = references.ref_fasta,
                docker              = global.sorter_docker,
                preemptible_tries   = preemptible_tries,
                monitoring_script   = monitoring_script,  # !FileCoercion
                cpu                 = cpu
        }

    RuntimeParams fastqc_runtime_params = {}
        call QCTasks.FastQC {
            input:
                input_fastq             = [select_first([ConvertToFastq.output_fastq])],
                limits                  = fastqc_limits,
                fastqc_runtime_params   = fastqc_runtime_params,
                docker                  = global.fastqc_docker,
                monitoring_script       = monitoring_script #!FileCoercion
        }
    }

    if (sort) {
        call UGGeneralTasks.ConvertSorterStatsToH5 as ConvertSorterStatsToH5 {
            input:
                monitoring_script = monitoring_script,
                base_file_name    = base_file_name,
                preemptible_tries = preemptible_tries,
                docker            = global.ug_vc_docker,
                no_address        = no_address,
                input_csv_file    = select_first([Sorter.sorter_stats_csv]),
                input_json_file   = select_first([Sorter.sorter_stats_json])
        }

        call QCTasks.CreateReportSingleSampleQC as CreateReport {
            input:
                monitoring_script = monitoring_script,
                base_file_name    = base_file_name,
                preemptible_tries = preemptible_tries,
                docker            = global.ug_vc_docker,
                disk_size         = 4,
                input_h5_file = ConvertSorterStatsToH5.aggregated_metrics_h5
        }
    }

    # select the output of the last stage that ran as the workflow output, the rest remain intermediate outputs
    File output_cram_bam_ = select_first([Sorter.sorted_cram, align_output, trimmer_output_ucram])
    
    output {
        # main output
        File output_cram_bam            = output_cram_bam_
        File? output_cram_bam_index     = Sorter.sorted_cram_index

        # single sample qc outputs
        File? aggregated_metrics_h5      = ConvertSorterStatsToH5.aggregated_metrics_h5
        File? aggregated_metrics_json    = ConvertSorterStatsToH5.aggregated_metrics_json
        File? report_html                = CreateReport.report_html

        # trim outputs
        File? trimmer_stats             = Trimmer.trimmer_stats
        File? trimmer_failure_codes_csv = Trimmer.trimmer_failure_codes_csv
        Array[File?]? trimmer_histogram = Trimmer.histogram
        Array[File?]? trimmer_histogram_extra = Trimmer.histogram_extra

        ## STAR outputs
        File? align_star_reads_per_gene_file = StarAlignment.reads_per_gene_file
        File? align_star_stats          = StarAlignment.star_stats

        # sorter outputs
        File? sort_stats_csv            = Sorter.sorter_stats_csv
        File? sort_stats_json           = Sorter.sorter_stats_json
        File? bedgraph_mapq0            = Sorter.sorter_out_bedgraph_mapq0
        File? bedgraph_mapq1            = Sorter.sorter_out_bedgraph_mapq1

        # fastqc outputs
        Array[File]? fastqc_reports     = FastQC.reports_html
    }
}
