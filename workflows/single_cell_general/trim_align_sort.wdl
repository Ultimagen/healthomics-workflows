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
        String pipeline_version = "1.15.5" # !UnusedDeclaration
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

        # general parameters
        Boolean no_address = true
        Int preemptible_tries = 1
        Int cpu

        # Used for running on other clouds (aws)
        File? monitoring_script_input
        String dummy_input_for_call_caching = ""

    
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name or 'sample' in base_file_name)
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
        #@wv 'sort' in steps and steps['sort'] -> defined(sorter_params))

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
            "CreateReport.top_metrics_file"
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
        sample_name: {
            help: "Sample name for the merged cram file.",
            type: "String", 
            category: "input_optional"
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
        sort_stats_csv: {
            help: "Sorter stats in csv format.",
            type: "File", 
            category: "output"
        }
        sort_stats_json: {
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
        fastq_files: {
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
                parameters  = select_first([trimmer_parameters]),
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
        File input_for_sort = select_first([align_output, trimmer_output_ucram,  single_cram_bam])

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
        if (length(Sorter.sorter_stats_csv) == 1) {
            call UGGeneralTasks.ConvertSorterStatsToH5 as ConvertSorterStatsToH5 {
                input:
                    monitoring_script = monitoring_script,
                    base_file_name    = base_file_name,
                    preemptible_tries = preemptible_tries,
                    docker            = global.ug_vc_docker,
                    no_address        = no_address,
                    input_csv_file    = Sorter.sorter_stats_csv[0],
                    input_json_file   = Sorter.sorter_stats_json[0]
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

        if (length(Sorter.sorted_cram) > 0) {
            File sorted_cram = select_first(Sorter.sorted_cram)
        }
        if (length(Sorter.sorted_cram_index) > 0) {
            File output_cram_bam_index_ = select_first(Sorter.sorted_cram_index)
        }
        if (length(Sorter.sorter_out_bedgraph_mapq0) > 0) {
            File bedgraph_mapq0_ = select_first(Sorter.sorter_out_bedgraph_mapq0)
        }
        if (length(Sorter.sorter_out_bedgraph_mapq1) > 0) {
            File bedgraph_mapq1_ = select_first(Sorter.sorter_out_bedgraph_mapq1)
        }
    }

    # select the output of the last stage that ran as the workflow output, the rest remain intermediate outputs.
    # Note about sorter output: sometimes there is no cram output (only fastq), so the output is selected from the previous stages
    File output_cram_bam_ = select_first([sorted_cram, align_output, trimmer_output_ucram, single_cram_bam])


    output {
        # main output
        File output_cram_bam            = output_cram_bam_
        File? output_cram_bam_index     = output_cram_bam_index_

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
        Array[File?]? sort_stats_csv     = Sorter.sorter_stats_csv
        Array[File?]? sort_stats_json    = Sorter.sorter_stats_json
        File? bedgraph_mapq0             = bedgraph_mapq0_
        File? bedgraph_mapq1             = bedgraph_mapq1_
        Array[File?]? fastq_files        = Sorter.fastq_files
        Array[File?]? sub_sampled_output = Sorter.sub_sampled_output
        File? unmatched_cram             = Sorter.unmatched_cram
    }
}
