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
# LICENSE for using 10x barcodes:
#   Copyright (c) 2020 10x Genomics
#             
#              Permission is hereby granted, free of charge, to any person obtaining a copy
#              of this software and associated documentation files (the "Software"), to deal
#              in the Software without restriction, including without limitation the rights
#              to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#              copies of the Software, and to permit persons to whom the Software is
#              furnished to do so, subject to the following conditions:
#             
#              1. The above copyright notice and this permission notice shall be included in all
#              copies or substantial portions of the Software.
#             
#              2. The above rights granted in the Software may be exercised only in connection 
#              with a 10x Genomics Product, rightfully purchased from 10x Genomics or an 
#              authorized reseller, or data generated using such a 10x Genomics Product. A 
#              10X Genomics Product means, collectively, 10x Genomics branded instruments, 
#              reagents, consumables, kits, and labware used in accordance with 10X Genomics
#              Product Terms and Conditions of Sale or, if applicable, any written contract 
#              between you and 10x Genomics. The rights granted may also be exercised in 
#              connection with other products when doing so is an integral part of an experiment 
#              where the data is generated primarily using a 10x Genomics Product.
#           
#              3. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#              IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#              FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#              AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#              LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#              OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#              SOFTWARE.
#
#
#
#
# DESCRIPTION
#   This workflow processes Ultima single-cell data and creates simulated paired-end reads.
# CHANGELOG
# 1.13.0 [BIOIN-1663] Use new versions of Trimmer and Sorter and run costum single cell qc script
# 1.4.0 [BIOIN-885] Renamed from single_cell_10x.wdl StarSolo based downstream analysis, updated trimmer, optimizations

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/single_cell_tasks.wdl" as SingleCellTasks
import "star_solo.wdl" as StarSoloWdl
import "star_align_gene_count.wdl" as StarAlignWorkflow
import "trim_align_sort.wdl" as TrimAlignSortSubWF

workflow SingleCell {
    input {
        String pipeline_version = "1.14.1" # !UnusedDeclaration

        File input_file
        String base_file_name

        # Trimming and sorting parameters
        TrimAlignSortSteps steps
        Array[File] ref_fastas_cram
        # References
        References references
        # trimmer parameters
        TrimmerParameters trimmer_parameters
        # sorter parameters
        SorterParams sorter_params

        String insert_rg
        String barcode_rg
        File star_genome

        SingleCellQcThresholds qc_thresholds

        # general parameters
        Boolean no_address = true
        Int preemptible_tries = 1
        Int cpu
        File? monitoring_script_input
        
        # STAR and STAR solo parameters
        String? downstream_analysis
        StarSoloParams? star_solo_params
        StarGenomeGenerateParams? genome_generate_params
        String? star_align_extra_args
        File? star_align_gtf_override

        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name or 'sample' in base_file_name)
        #@wv suffix(input_file) in {".bam", ".cram"}
        #@wv defined(trimmer_parameters) and not 'formats_description' in trimmer_parameters -> 'format' in trimmer_parameters

        #@wv 'trim' in steps or 'align' in steps or 'sort' in steps -> (steps['trim'] or steps['align'] or steps['sort'])
        #@wv defined(references) -> len(references) == 3
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa','.fna'}
        #@wv suffix(references['ref_dict']) == '.dict'
        #@wv suffix(references['ref_fasta_index']) == '.fai'
        #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
        #@wv prefix(references['ref_dict']) == prefix(references['ref_fasta'])
        
        # Trimmer checks
        #@wv 'trim' in steps and steps['trim'] -> defined(trimmer_parameters)
        #@wv defined(trimmer_parameters) and not 'formats_description' in trimmer_parameters -> 'format' in trimmer_parameters
        ##@wv 'formats_description' in trimmer_parameters -> not('local_formats_description' in trimmer_parameters)
        ##@wv 'local_formats_description' in trimmer_parameters -> not('formats_description' in trimmer_parameters)
        # validate that the user updated trimmer extra_args with the correct read_group of the sample
        #@wv not("<!READ GROUP!>" in trimmer_parameters['extra_args'])

        # Sort checks
        #@wv 'sort' in steps and steps['sort'] -> defined(sorter_params))

        # STAR solo validations
        #@wv defined(downstream_analysis) and downstream_analysis == "star_solo" -> defined(star_solo_params)
        #@wv defined(downstream_analysis) and downstream_analysis == "star_solo" -> defined(star_solo_params['genome']) or defined(genome_generate_params)
        #@wv defined(downstream_analysis) and downstream_analysis == "star_solo" and defined(star_solo_params['genome']) -> suffix(star_solo_params['genome']) == '.zip'
        #@wv defined(downstream_analysis) and downstream_analysis == "star_solo" -> star_solo_params['strand'] == 'Reverse' or star_solo_params['strand'] == 'Forward'
        
        # STAR genome generation validations
        #@wv defined(downstream_analysis) and (downstream_analysis == "star_solo" or downstream_analysis == "star") and defined(genome_generate_params) -> suffix(genome_generate_params['fasta_files']) <= {'.fasta','.fa'}
        #@wv defined(downstream_analysis) and (downstream_analysis == "star_solo" or downstream_analysis == "star") and defined(genome_generate_params) -> suffix(genome_generate_params['gtf_file']) == '.gtf'
        
        # STAR validations
        #@wv defined(downstream_analysis) and downstream_analysis == "star" -> defined(star_genome) or defined(genome_generate_params)
        #@wv defined(downstream_analysis) and downstream_analysis == "star" and defined(star_genome) -> suffix(star_genome) == '.zip'

    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers
    
    meta {
        description : "Create simulated paired end fastq reads from Ultima single-ended CRAM or BAM.\n\n Intended for reads with a cell barcode, UMI and insert (e.g., 10x, Parse Biosciences, Fluent).\n\n In addition, runs fastqc on the insert and, optionally, runs FastQC and STAR or STARsolo."
        author: "Ultima Genomics"
        WDL_AID: { 
            exclude: [
                "pipeline_version",
                "monitoring_script_input",
                "StarSoloWorkflow.StarGenomeGenerate.disk_size",
                "StarSoloWorkflow.StarSolo.disk_size",
                "StarSoloWorkflow.StarAlignStats.disk_size",
                "StarSoloWorkflow.GatherStatistics.separate_reads_statistics",
                "StarAlignment.StarGenomeGenerate.disk_size",
                "StarAlignment.StarAlign.disk_size",
                "StarAlignment.StarAlignStats.disk_size",
                "StarSoloWorkflow.Globals.glob",
                "StarAlignment.Globals.glob",
                "Globals.glob",
                "TrimAlignSort.StarAlignment.memory_gb",
                "StarAlignment.memory_gb",
                "SingleCell.StarAlignSubSample.memory_gb"
            ]
        }
    }


    parameter_meta {
        base_file_name: {
            help: "Base file name for output files. The output files will be named <tt>[base_file_name]*.fastq.gz</tt>",
            type: "String", 
            category: "input_required"
        }
        input_file: {
            help: "Input CRAM or BAM file",
            type: "File",
            category: "input_required"
        }
        steps: {
            help: "The steps to run in the workflow (trim+sort)",
            type: "TrimAlignSortSteps",
            category: "input_required"
        }
        trimmer_parameters: {
            help: "Parameters for Trimmer task.  See input template",
            type: "TrimmerParameters",
            category: "optional"
        }
        sorter_params: {
            help: "Parameters for Sorter task.  See input template",
            type: "SorterParams",
            category: "optional"
        }
        insert_rg: {
            help: "Read group name for the insert reads, e.g. S1_L001_R2_001",
            type: "String",
            category: "input_required"
        }
        barcode_rg: {
            help: "Read group name for the barcode reads, e.g. S1_L001_R1_001",
            type: "String",
            category: "input_required"
        }
        references: {
            help: "References for the workflow",
            type: "References",
            category: "input_required"
        }
        ref_fastas_cram: {
            help: "Reference fasta files for the CRAM file",
            type: "Array[File]",
            category: "input_required"
        }
        qc_thresholds: {
            help: "Thresholds for the single cell qc",
            type: "SingleCellQcThresholds",
            category: "input_required"
        }
        downstream_analysis: {
            help: "Can be either star_solo, star, or undefined (default)",
            type: "String",
            category: "optional"
        }
        star_solo_params: {
            help: "Parameters for running the StarSolo task.  See input template",
            type: "StarSoloParams",
            category: "optional"
        }
        genome_generate_params: {
            help: "Parameters that can be passed to the StarSolo or StarAlignment task, for formatting a STAR genome.  See input template",
            type: "StarGenomeGenerateParams",
            category: "optional"
        }
        star_genome: {
            help: "The genome to use in the StarAlignment task.  See input template",
            type: "File",
            category: "reference_optional"
        }
        star_align_extra_args: {
            help: "Extra parameters to pass to the StarAlignment task.  See input template",
            type: "String",
            category: "param_optional"
        }
        star_align_gtf_override: {
            help: "The gtf to use in the StarAlignment task.  See input template",
            type: "File",
            category: "optional"
        }
        no_address: {
            type: "Boolean",
            help: "Should the instances used be without external IP address. Allows for more parallelization, but not supported with Dockerhub dockers Default: true",
            category: "param_required"
        }
        preemptible_tries: {
            type: "Int",
            help: "Number of preemptible tries",
            category: "param_required"
        }
        cpu: {
            type: "Int",
            help: "Number of CPUs to use (for Trimmer, conversion to fastq, and alignment tasks)", 
            category: "param_required"
        }
        output_barcodes_fastq: {
            type: "File",
            help: "The fastq with the barcodes portion of the read", 
            category: "output"
        }
        output_insert_fastq: {
            type: "File",
            help: "The fastq with the insert portion of the read", 
            category: "output"
        }
        report_html: {
            type: "File",
            help: "The report from the single cell qc", 
            category: "output"
        }
        aggregated_metrics_h5: {
            type: "File",
            help: "The h5 store from the single cell qc", 
            category: "output"
        }
        trimmer_stats_output: {
            type: "File",
            help: "Trimmer output statistics", 
            category: "output"
        }
        trimmer_failure_codes_csv: {
            type: "File",
            help: "Trimmer failure codes csv", 
            category: "output"
        }
        trimmer_histogram: {
            type: "Array[File]",
            help: "Trimmer histograms", 
            category: "output"
        }
        trimmer_histogram_extra: {
            type: "Array[File]",
            help: "Trimmer extra histograms", 
            category: "output"
        }
        sort_stats_csv : {
            type: "Array[File]",
            help: "Sorter statistics csv", 
            category: "output"
        }
        sort_stats_json : {
            type: "Array[File]",
            help: "Sorter statistics json", 
            category: "output"
        }
        star_solo_outputs: {
            type: "StarSoloOutputs",
            help: "The outputs from the StarSolo workflow.  Created only if STARsolo is run", 
            category: "output_optional"
        }
        star_bam: {
            type: "File",
            help: "The Star-aligned bam.  Created only if STAR is run.", 
            category: "output"
        }
        star_reads_per_gene_file: {
            type: "File",
            help: "The Star reads per gene file", 
            category: "output"
        }
        star_stats: {
            type: "File",
            help: "The Star output statistics",
            category: "output"
        }
        unmatched_cram: {
            help: "Unmatched cram file output from sorter",
            type: "File", 
            category: "output"
        }
    }


    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call TrimAlignSortSubWF.TrimAlignSort {
        input:
            input_cram_bam_list     = [input_file],
            base_file_name          = base_file_name,
            steps                   = steps,
            references              = references,
            ref_fastas_cram         = ref_fastas_cram,
            sample_name             = base_file_name,
            trimmer_parameters      = trimmer_parameters,
            sorter_params           = sorter_params,
            no_address              = no_address,
            preemptible_tries       = preemptible_tries,
            cpu                     = cpu,
            monitoring_script_input = monitoring_script
    }

    call SingleCellTasks.FindInsertBarcodeFastq {
        input:
            input_fastq_list        = select_all(select_first([TrimAlignSort.fastq_files])),
            sub_sumple_fastq_list   = select_all(select_first([TrimAlignSort.sub_sampled_output])),
            sorter_csv_stats_list   = select_all(select_first([TrimAlignSort.sort_stats_csv])),
            base_file_name          = base_file_name,
            insert_rg               = insert_rg,
            barcode_rg              = barcode_rg,
            no_address              = no_address,
            monitoring_script       = monitoring_script,
            docker                  = global.single_cell_qc_docker
    }

    String sub_sample_star_align_extra_args = "--outSAMunmapped Within --chimOutType WithinBAM SoftClip --clip3pNbases 0 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --scoreDelOpen -2 --scoreDelBase -2 --scoreInsOpen -2--scoreInsBase -2 --alignEndsType Local --outSAMmapqUnique 60"
    call StarAlignWorkflow.StarAlignment as StarAlignSubSample{
        input:
            genome                  = star_genome,
            input_bams              = [FindInsertBarcodeFastq.insert_sub_sample_fastq],
            base_file_name          = base_file_name,
            star_align_extra_args   = sub_sample_star_align_extra_args,
            preemptible_tries       = preemptible_tries,
            no_address              = no_address,
            cpu                     = cpu,
            monitoring_script_input = monitoring_script
    }

    call SingleCellTasks.SingleCellQc{
        input:
            trimmer_stats               = select_first([TrimAlignSort.trimmer_stats]),
            trimmer_histogram           = select_first([select_first([TrimAlignSort.trimmer_histogram])[0]]),
            trimmer_failure_codes       = select_first([TrimAlignSort.trimmer_failure_codes_csv]),
            sorter_stats_csv            = FindInsertBarcodeFastq.insert_sorter_stats_csv,
            star_stats                  = StarAlignSubSample.raw_star_log_file,
            star_reads_per_gene         = StarAlignSubSample.reads_per_gene_file,
            insert_sub_sample_fastq     = FindInsertBarcodeFastq.insert_sub_sample_fastq,
            base_file_name              = base_file_name,
            qc_thresholds               = qc_thresholds,
            monitoring_script           = monitoring_script,
            memory_gb                   = 16,
            preemptible_tries           = preemptible_tries,
            no_address                  = no_address,
            docker                      = global.single_cell_qc_docker
    }

    if (defined(downstream_analysis)) {
        String analysis_type = select_first([downstream_analysis])
        if(analysis_type == "star_solo"){
            call StarSoloWdl.StarSoloWorkflow {
                input:
                    insert_fastq            = FindInsertBarcodeFastq.insert_fastq,
                    barcode_fastq           = FindInsertBarcodeFastq.barcode_fastq,
                    base_file_name          = base_file_name,
                    star_solo_params        = select_first([star_solo_params]),
                    genome_generate_params  = genome_generate_params,
                    no_address              = no_address,
                    preemptible_tries       = preemptible_tries,
                    cpu                     = cpu,
                    monitoring_script_input = monitoring_script
            }
        }

        if (analysis_type == "star") {
            call StarAlignWorkflow.StarAlignment {
                input:
                    genome                  = star_genome,
                    genome_generate_params  = genome_generate_params,
                    star_align_extra_args   = star_align_extra_args,
                    star_align_gtf_override = star_align_gtf_override,
                    input_bams              = [TrimAlignSort.output_cram_bam],
                    base_file_name          = base_file_name + ".star.aln",
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
                    cpu                     = cpu,
                    monitoring_script_input = monitoring_script
            }
        }
    }

    output {
        #TrimAlignSort outputs
        File? trimmer_stats_output          = TrimAlignSort.trimmer_stats
        File? trimmer_failure_codes_csv     = TrimAlignSort.trimmer_failure_codes_csv
        Array[File?]? trimmer_histogram     = TrimAlignSort.trimmer_histogram
        Array[File?]? trimmer_histogram_extra = TrimAlignSort.trimmer_histogram_extra
        Array[File?]? sort_stats_csv        = TrimAlignSort.sort_stats_csv
        Array[File?]? sort_stats_json       = TrimAlignSort.sort_stats_json
        File? unmatched_cram                = TrimAlignSort.unmatched_cram

        # Fastq outputs
        File output_barcodes_fastq          = FindInsertBarcodeFastq.barcode_fastq
        File output_insert_fastq            = FindInsertBarcodeFastq.insert_fastq

        # SingleCellQc outputs
        File report_html                    = SingleCellQc.report
        File aggregated_metrics_h5          = SingleCellQc.h5

        # Star solo outputs
        StarSoloOutputs? star_solo_outputs = StarSoloWorkflow.outputs

        # STAR outputs
        File? star_bam                    = StarAlignment.output_bam
        File? star_reads_per_gene_file    = StarAlignment.reads_per_gene_file
        File? star_stats                  = StarAlignment.star_stats
    }
}
