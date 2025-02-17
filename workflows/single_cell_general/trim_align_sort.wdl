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
        String pipeline_version = "1.17.2" # !UnusedDeclaration
        Array[File] input_cram_bam_list
        Array[File] ref_fastas_cram
        String base_file_name
        TrimAlignSortSteps steps
        References references

        # trimmer parameters
        TrimmerParameters? trimmer_parameters

        # alignment parameters
        String? aligner # ua, ua-meth, star
        UaParameters? ua_parameters
        UaMethReferences? ua_meth_parameters

        ## STAR param
        File? star_genome
        StarGenomeGenerateParams? star_genome_generate_params
        String? star_align_extra_args
        File? star_align_gtf_override

        # sorter parameters (for marking duplicates)
        SorterParams? sorter_params

        # general parameters
        Boolean no_address = true
        Int preemptible_tries = 1
        Int cpu

        # Used for running on other clouds (aws)
        File? monitoring_script_input
        String dummy_input_for_call_caching = ""


        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name or 'sample' in base_file_name)
        #@wv not("_sample" in base_file_name or "_unmatched" in base_file_name)
        #@wv not("_sample" in base_file_name or "_unmatched" in base_file_name)
        #@wv len(input_cram_bam_list) > 0
        #@wv suffix(input_cram_bam_list) <= {".bam", ".cram"}
        #@wv 'trim' in steps or 'align' in steps or 'sort' in steps -> (steps['trim'] or steps['align'] or steps['sort'])
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
        #@wv defined(trimmer_parameters) and ('failure_read_group_args' in trimmer_parameters) -> not('output_trimmed_failed_file_name' in trimmer_parameters)

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
        #@wv 'sort' in steps and steps['sort'] <-> defined(sorter_params)

        ## Integration checks
        #@wv 'trim' in steps and steps['trim'] and 'align' in steps and not(steps['align']) -> not('sort' in steps and steps['sort'] and defined(sorter_params) and 'aligned' in sorter_params and sorter_params['aligned'])

    }


    meta {
        description : "Pipeline for trimming, aligning and sorting Ultima data in a fast, cost-effective, and easier-to-maintain way, similiar to the way it is done on-tool. It can be used to run the entire pipeline or just a subset of the steps."
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "no_address",
            "preemptible_tries",
            "dummy_input_for_call_caching",
            "monitoring_script_input",
            "Globals.glob",
            "CreateReferenceCache.disk_size",
            "Trimmer.disk_size",
            "UAAlignment.Globals.glob",
            "UAAlignment.BuildUaIndex.disk_size",
            "UAAlignment.AlignWithUA.v_aware_vcf",
            "UAAlignment.AlignWithUA.disk_size",
            "UAAlignment.AlignWithUA.cpu",
            "UAMethAlignment.UaMethIntensiveMode",
            "UAMethAlignment.Globals.glob",
            "UAMethAlignment.BuildUaMethIndex.disk_size",
            "UAMethAlignment.AlignWithUAMeth.disk_size",
            "StarAlignment.monitoring_script_input",
            "StarAlignment.Globals.glob",
            "StarAlignment.StarGenomeGenerate.disk_size",
            "StarAlignment.StarAlign.disk_size",
            "StarAlignment.StarAlignStats.disk_size",
            "Sorter.mapq_override",
            "ConvertSorterStatsToH5.disk_size",
            "ConvertSorterStatsToH5.metric_mapping_file",
            "CreateReport.notebook_file_in",
            "CreateReport.top_metrics_file",
            "Demux.mapq_override"
        ]}
    }

    parameter_meta {
        input_cram_bam_list: {
            help: "List of input cram or bam files to be processed.",
            type: "String",
            category: "input_required"
        }
        ref_fastas_cram: {
            help: "List of references for CreateReferenceCache task.",
            type: "String",
            category: "input_required"
        }
        base_file_name: {
            help: "Base name for the output files.",
            type: "String",
            category: "input_required"
        }
        steps: {
            help: "Steps to be executed in the pipeline.Options are: trim, align, sort",
            type: "String",
            category: "input_required"
        }
        references: {
            help: "References for merging inputs into one file, alignment, and sorting.",
            type: "String",
            category: "input_required"
        }
        trimmer_parameters: {
            help: "Parameters for the trimmer task. Mandatory if trim step is selected.",
            type: "String",
            category: "input_optional"
        }
        aligner: {
            help: "Aligner to be used. Options are: ua, ua-meth, star. Mandatory if align step is selected.",
            type: "String",
            category: "input_optional"
        }
        ua_parameters: {
            help: "Parameters for the UA aligner. Mandatory if aligner is ua.",
            type: "String",
            category: "input_optional"
        }
        ua_meth_parameters: {
            help: "Parameters for the UA meth aligner. Mandatory if aligner is ua-meth.",
            type: "String",
            category: "input_optional"
        }
        star_genome: {
            help: "Star genome file. If aligner is star, supllay either star genome file or generate new genome index.",
            type: "String",
            category: "input_optional"
        }
        star_genome_generate_params: {
            help: "Parameters for generating the star genome. Mandatory if aligner is star and not given star genome file.",
            type: "String",
            category: "input_optional"
        }
        star_align_extra_args: {
            help: "Extra arguments for the STAR aligner.",
            type: "String",
            category: "input_optional"
        }
        star_align_gtf_override: {
            help: "GTF file to be used for STAR aligner.",
            type: "String",
            category: "input_optional"
        }
        sorter_params: {
            help: "Parameters for the sorter task. Mandatory if sort step is selected.",
            type: "String",
            category: "input_optional"
        }
        cpu: {
            help: "Number of cpus to be used for the tasks.",
            type: "Int",
            category: "input_required"
        }
        output_cram_bam: {
            help: "Output file after the pipeline is executed.",
            type: "File",
            category: "output"
        }
        output_cram_bam_index: {
            help: "Index file for the output cram file.",
            type: "File",
            category: "output"
        }
        aggregated_metrics_h5: {
            help: "Aggregated metrics in h5 format.",
            type: "File",
            category: "output"
        }
        aggregated_metrics_json: {
            help: "Aggregated metrics in json format.",
            type: "File",
            category: "output"
        }
        report_html: {
            help: "Report in html format.",
            type: "File",
            category: "output"
        }
        trimmer_stats: {
            help: "Trimmer stats file.",
            type: "File",
            category: "output"
        }
        trimmer_failure_codes_csv: {
            help: "Trimmer failure codes in csv format.",
            type: "File",
            category: "output"
        }
        trimmer_histogram: {
            help: "Trimmer histogram files.",
            type: "Array[File]",
            category: "output"
        }
        trimmer_histogram_extra: {
            help: "Trimmer histogram extra files.",
            type: "Array[File]",
            category: "output"
        }
        align_star_reads_per_gene_file: {
            help: "STAR reads per gene file.",
            type: "File",
            category: "output"
        }
        align_star_stats: {
            help: "STAR stats file.",
            type: "File",
            category: "output"
        }
        sorter_stats_csv: {
            help: "Sorter stats in csv format.",
            type: "File",
            category: "output"
        }
        sorter_stats_json: {
            help: "Sorter stats in json format.",
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
        fastq_file: {
            help: "Sorter output in Fastq format.",
            type: "File",
            category: "output"
        }
        fastq_file_list: {
            help: "Sorter output in Fastq format (allow multiple files for applications like single-cell).",
            type: "Array[File]",
            category: "output"
        }
        sub_sampled_output: {
            help: "Sub-sampling files as part of sorter output.",
            type: "Array[File]",
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
        ua_stats_jsons: {
            help: "Statistics file in json format from UA.",
            type: "Array[File]",
            category: "output"
        }
        output_cram_bam_list: {
            help: "Output file list after the pipeline is executed in multiple outputs mode",
            type: "Array[File]",
            category: "output"
        }
        output_cram_bam_index_list: {
            help: "Index files for the output cram files (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        sorter_stats_csv_list: {
            help: "Sorter stats in csv format (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        sorter_stats_json_list: {
            help: "Sorter stats in json format (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        bedgraph_mapq0_list: {
            help: "Bedgraph mapq0 files (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        bedgraph_mapq1_list: {
            help: "Bedgraph mapq1 files (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        unmatched_cram_list: {
            help: "Unmatched cram file output from sorter (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        unmatched_sorter_stats_csv_list: {
            help: "Unmatched output cram files sorter stats in csv format (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        unmatched_sorter_stats_json_list: {
            help: "Unmatched output cram files sorter stats in json format (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        aggregated_metrics_h5_list: {
            help: "Aggregated metrics in h5 format (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        aggregated_metrics_json_list: {
            help: "Aggregated metrics in json format (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        report_html_list: {
            help: "Reports in html format (if defined, when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    Boolean trim = select_first([steps.trim, false])
    Boolean align = select_first([steps.align, false])
    Boolean sort = select_first([steps.sort, false])

    call UGAlignment.CreateReferenceCache {
        input:
            references = ref_fastas_cram,
            cache_populate_script = global.ref_cache_script, #!StringCoercion
            preemptible_tries = preemptible_tries,
            docker = global.perl_docker,
            dummy_input_for_call_caching = dummy_input_for_call_caching
    }

    if (trim) {
        call TrimmingTasks.Trimmer {
            input:
                input_cram_bam_list = input_cram_bam_list,
                cache_tarball       = CreateReferenceCache.cache_tarball,
                parameters  = select_first([trimmer_parameters]),
                base_file_name      = base_file_name,
                docker              = global.trimmer_docker,
                preemptible_tries   = preemptible_tries,
                monitoring_script   = monitoring_script,  # !FileCoercion
                no_address          = no_address,
                cpus                = cpu
        }
        # trimmer output should be one file (trimming tasks allow multiple outputs in case trimmer is 
        # doing demultiplexing, but this is not a valid option for this workflow)
        Array[File] trimmer_output_ucram_list = Trimmer.trimmed_ucram_list
    }

    if (align) {
        Array[File] input_for_alignment_list = select_first([trimmer_output_ucram_list, input_cram_bam_list])
        String aligner_override  = select_first([aligner])

        if (aligner_override == "ua" ){
            call UaAlignWorkflow.UAAlignment {
                input:
                    input_files             = input_for_alignment_list,
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
                    input_files             = input_for_alignment_list,
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
                    input_bams_or_fastqs    = input_for_alignment_list,
                    base_file_name          = base_file_name + ".star.aln",
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
                    cpu                     = cpu
            }
        }
        Array[File] align_output_list = [select_first([UAAlignment.ua_output_bam, UAMethAlignment.ua_output_bam ,StarAlignment.output_bam])]
    }

    Array[File] input_for_sort_list = select_first([align_output_list, trimmer_output_ucram_list, input_cram_bam_list])
    if (!sort) {
        if (length(input_for_sort_list) == 1) {
            File output_cram_bam_without_sorter_ = select_first(input_for_sort_list)
        }
        if (length(input_for_sort_list) > 1) {
            Array[File] output_cram_bam_list_without_sorter_ = select_all(select_first([input_for_sort_list]))
        }
    }
    if (sort) {
        call SortTasks.Demux {
            input:
                input_cram_bam_list= input_for_sort_list,
                cache_tarball      = CreateReferenceCache.cache_tarball,
                base_file_name     = base_file_name,
                reference_fasta    = references.ref_fasta,
                sorter_params      = select_first([sorter_params]),
                monitoring_script  = monitoring_script,  # !FileCoercion
                docker             = global.sorter_docker,
                preemptible_tries  = preemptible_tries,
                cpu                = cpu
        }

        call SortTasks.Sorter {
            input:
                demux_output       = Demux.demux_output,
                max_region_size    = Demux.max_region_size,
                base_file_name     = base_file_name,
                reference_fasta    = references.ref_fasta,
                sorter_params      = select_first([sorter_params]),
                monitoring_script  = monitoring_script,  # !FileCoercion
                docker             = global.sorter_docker,
                preemptible_tries  = preemptible_tries,
        }

        # Set workflow output variables - a single file is there is just one, a list if there is more
        if (length(Sorter.sorted_cram) == 1) {
            File output_cram_bam_with_sorter_ = select_first(Sorter.sorted_cram)
        }
        if (length(Sorter.sorted_cram) > 1) {
            Array[File] output_cram_bam_list_with_sorter_ = select_all(select_first([Sorter.sorted_cram]))
        }

        if (length(Sorter.fastq_files) == 1) {
            File fastq_file_ = select_first(Sorter.fastq_files)
        }
        if (length(Sorter.fastq_files) > 1) {
            Array[File] fastq_file_list_ = select_all(select_first([Sorter.fastq_files]))
        }

        if (length(Sorter.sorted_cram_index) == 1) {
            File output_cram_bam_index_ = select_first(Sorter.sorted_cram_index)
        }
        if (length(Sorter.sorted_cram_index) > 1) {
            Array[File] output_cram_bam_index_list_ = select_all(select_first([Sorter.sorted_cram_index]))
        }

        if (length(Sorter.sorter_out_bedgraph_mapq0) == 1) {
            File bedgraph_mapq0_ = select_first(Sorter.sorter_out_bedgraph_mapq0)
        }
        if (length(Sorter.sorter_out_bedgraph_mapq0) > 1) {
            Array[File] bedgraph_mapq0_list_ = select_all(select_first([Sorter.sorter_out_bedgraph_mapq0]))
        }

        if (length(Sorter.sorter_out_bedgraph_mapq1) == 1) {
            File bedgraph_mapq1_ = select_first(Sorter.sorter_out_bedgraph_mapq1)
        }
        if (length(Sorter.sorter_out_bedgraph_mapq1) > 1) {
            Array[File] bedgraph_mapq1_list_ = select_all(select_first([Sorter.sorter_out_bedgraph_mapq1]))
        }

        if (length(Sorter.sorter_stats_json) == 1) {
            File sorter_stats_json_ = select_first(Sorter.sorter_stats_json)
        }
        if (length(Sorter.sorter_stats_json) > 1) {
            Array[File] sorter_stats_json_list_ = select_all(select_first([Sorter.sorter_stats_json]))
        }

        if (length(Sorter.sorter_stats_csv) == 1) {
            File sorter_stats_csv_ = select_first(Sorter.sorter_stats_csv)
        }
        if (length(Sorter.sorter_stats_csv) > 1) {
            Array[File] sorter_stats_csv_list_ = select_all(select_first([Sorter.sorter_stats_csv]))
        }

        if (length(Sorter.unmatched_cram) == 1) {
            File unmatched_cram_ = select_first(Sorter.unmatched_cram)
        }
        if (length(Sorter.unmatched_cram) > 1) {
            Array[File?] unmatched_cram_list_ = Sorter.unmatched_cram
        }

        if (length(Sorter.unmatched_sorter_stats_csv) == 1) {
            File unmatched_sorter_stats_csv_ = select_first(Sorter.unmatched_sorter_stats_csv)
        }
        if (length(Sorter.unmatched_sorter_stats_csv) > 1) {
            Array[File?] unmatched_sorter_stats_csv_list_ = Sorter.unmatched_sorter_stats_csv
        }

        if (length(Sorter.unmatched_sorter_stats_json) == 1) {
            File unmatched_sorter_stats_json_ = select_first(Sorter.unmatched_sorter_stats_json)
        }
        if (length(Sorter.unmatched_sorter_stats_json) > 1) {
            Array[File?] unmatched_sorter_stats_json_list_ = Sorter.unmatched_sorter_stats_json
        }

        Boolean create_qc_reports = length(Sorter.sorter_stats_csv) > 0
        if (create_qc_reports) {
            # If there are any output csvs, convert the sorter stats to h5 and create a report for each.
            scatter (i in range(length(Sorter.sorter_stats_csv))) {
                call UGGeneralTasks.ConvertSorterStatsToH5 {
                    input:
                        monitoring_script = monitoring_script,
                        file_name_suffix  = i,
                        preemptible_tries = preemptible_tries,
                        docker            = global.ugbio_core_docker,
                        no_address        = no_address,
                        input_csv_file    = Sorter.sorter_stats_csv[i],
                        input_json_file   = Sorter.sorter_stats_json[i]
                }

                call QCTasks.CreateReportSingleSampleQC as CreateReport {
                    input:
                        monitoring_script = monitoring_script,
                        base_file_name    = base_file_name,
                        preemptible_tries = preemptible_tries,
                        docker            = global.ugbio_core_docker,
                        disk_size         = 4,
                        input_h5_file     = ConvertSorterStatsToH5.aggregated_metrics_h5
                }
            }
            
            # Set workflow output variables - a single file is there is just one, a list if there is more
            if (length(ConvertSorterStatsToH5.aggregated_metrics_h5) == 1) {
                File aggregated_metrics_h5_ = select_first(ConvertSorterStatsToH5.aggregated_metrics_h5)
            }
            if (length(ConvertSorterStatsToH5.aggregated_metrics_h5) > 1) {
                Array[File] aggregated_metrics_h5_list_ = ConvertSorterStatsToH5.aggregated_metrics_h5
            }

            if (length(ConvertSorterStatsToH5.aggregated_metrics_json) == 1) {
                File aggregated_metrics_json_ = select_first(ConvertSorterStatsToH5.aggregated_metrics_json)
            }
            if (length(ConvertSorterStatsToH5.aggregated_metrics_json) > 1) {
                Array[File] aggregated_metrics_json_list_ = ConvertSorterStatsToH5.aggregated_metrics_json
            }

            if (length(CreateReport.report_html) == 1) {
                File report_html_ = select_first(CreateReport.report_html)
            }
            if (length(CreateReport.report_html) > 1) {
                Array[File] report_html_list_ = CreateReport.report_html
            }
        }
    }

    # define final output cram and index files
    if (defined(output_cram_bam_without_sorter_) || defined(output_cram_bam_with_sorter_)) {
        File output_cram_bam_ = select_first([output_cram_bam_with_sorter_, output_cram_bam_without_sorter_])
    }
    if (defined(output_cram_bam_list_without_sorter_) || defined(output_cram_bam_list_with_sorter_)) {
        Array[File] output_cram_bam_list_ = select_first([output_cram_bam_list_with_sorter_, output_cram_bam_list_without_sorter_])
    }

    output {
        # single output mode outputs
        File? output_cram_bam               = output_cram_bam_
        File? output_cram_bam_index         = output_cram_bam_index_
        File? fastq_file                    = fastq_file_
        File? sorter_stats_csv              = sorter_stats_csv_
        File? sorter_stats_json             = sorter_stats_json_
        File? unmatched_cram                = unmatched_cram_
        File? unmatched_sorter_stats_csv    = unmatched_sorter_stats_csv_
        File? unmatched_sorter_stats_json   = unmatched_sorter_stats_json_
        File? bedgraph_mapq0                = bedgraph_mapq0_
        File? bedgraph_mapq1                = bedgraph_mapq1_
        # single sample qc outputs
        File? aggregated_metrics_h5         = aggregated_metrics_h5_
        File? aggregated_metrics_json       = aggregated_metrics_json_
        File? report_html                   = report_html_

        # FastQ output mode
        Array[File?]? fastq_file_list       = fastq_file_list_
        Array[File?]? sub_sampled_output = Sorter.sub_sampled_output

        # multiple output cram mode outputs
        Array[File?]? output_cram_bam_list       = output_cram_bam_list_
        Array[File?]? output_cram_bam_index_list = output_cram_bam_index_list_
        Array[File?]? sorter_stats_csv_list      = sorter_stats_csv_list_
        Array[File?]? sorter_stats_json_list     = sorter_stats_json_list_
        Array[File?]? unmatched_cram_list        = unmatched_cram_list_
        Array[File?]? unmatched_sorter_stats_csv_list = unmatched_sorter_stats_csv_list_
        Array[File?]? unmatched_sorter_stats_json_list = unmatched_sorter_stats_json_list_
        Array[File?]? bedgraph_mapq0_list        = bedgraph_mapq0_list_
        Array[File?]? bedgraph_mapq1_list        = bedgraph_mapq1_list_
#        # single sample qc outputs
        Array[File?]? aggregated_metrics_h5_list     = aggregated_metrics_h5_list_
        Array[File?]? aggregated_metrics_json_list   = aggregated_metrics_json_list_

        Array[File?]? report_html_list                = report_html_list_

        # trim outputs
        File? trimmer_stats                     = Trimmer.trimmer_stats
        File? trimmer_failure_codes_csv         = Trimmer.trimmer_failure_codes_csv
        Array[File?]? trimmer_histogram         = Trimmer.histogram
        Array[File?]? trimmer_histogram_extra   = Trimmer.histogram_extra

        ## UA outputs
        Array[File]? ua_stats_jsons             = UAAlignment.ua_output_json

        ## STAR outputs
        File? align_star_reads_per_gene_file    = StarAlignment.reads_per_gene_file
        File? align_star_stats                  = StarAlignment.star_stats
    }
}
