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
# 1.4.0 [BIOIN-885] Renamed from single_cell_10x.wdl StarSolo based downstream analysis, updated trimmer, optimizations

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/trimming_tasks.wdl" as TrimmingTasks
import "tasks/qc_tasks.wdl" as QCTasks
import "tasks/single_cell_tasks.wdl" as SingleCellTasks
import "star_solo.wdl" as StarSoloWdl
import "star_align_gene_count.wdl" as StarAlignWorkflow
import "tasks/sorting_tasks.wdl" as SortTasks

workflow SingleCell {
    # Trimming step is optional. If trimmer_parameters are not provided, the input file is used as is.
    input {
        File input_file
        String base_file_name
        TrimmerParameters? trimmer_parameters
        File? trimmer_stats # if trimmer run outside of this workflow, you can still combine the stats in the final report

        String? demux_extra_args
        String? barcode_fastq_file_suffix
        String? insert_fastq_file_suffix
        String? barcode_fastq_header_suffix
        String? insert_fastq_header_suffix

        String? fastqc_adapter
        File? fastqc_limits
        String pipeline_version = "1.13.2" # !UnusedDeclaration
        String? downstream_analysis
        # STAR and STAR solo parameters
        StarSoloParams? star_solo_params
        StarGenomeGenerateParams? genome_generate_params
        File? star_genome
        String? star_align_extra_args
        File? star_align_gtf_override

        Boolean no_address = true
        Int preemptible_tries = 1
        Int cpu

        # Used for running on other clouds (aws)
        File? monitoring_script_input

        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv suffix(input_file) in {".bam", ".cram"}
        #@wv defined(trimmer_parameters) and not 'formats_description' in trimmer_parameters -> 'format' in trimmer_parameters

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
                "Trimmer.disk_size",
                "ConvertToFastq.output_dir",
                "ConvertToFastq.samtools_extra_args",
                "ConvertToFastq.reference_fasta",
                "CreateSyntheticPairedEnd.disk_size",
                "StarSoloWorkflow.StarGenomeGenerate.disk_size",
                "StarSoloWorkflow.StarSolo.disk_size",
                "StarSoloWorkflow.StarAlignStats.disk_size",
                "StarSoloWorkflow.GatherStatistics.separate_reads_statistics",
                "StarAlignment.StarGenomeGenerate.disk_size",
                "StarAlignment.StarAlign.disk_size",
                "StarAlignment.StarAlignStats.disk_size",
                "CombineStatistics.disk_size",
                "StarSoloWorkflow.Globals.glob",
                "StarAlignment.Globals.glob",
                "Globals.glob"
            ]}
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
        trimmer_parameters: {
            help: "Parameters for Trimmer task.  See input template",
            type: "TrimmerParameters",
            category: "optional"
        }
        trimmer_stats: {
            help: "If trimmer was run outside of this workflow, the stats can still be combined in the final report",
            type: "File",
            category: "optional"
        }
        demux_extra_args: {
            help: "Extra parameters to pass to SortTasks.Demux , when converting to fastq",
            type: "String",
            category: "param_optional"
        }
        barcode_fastq_file_suffix: {
            help: "Suffix to add to the name of the file with the barcode reads (in case downstream software has name requirements)",
            type: "String",
            category: "optional"
        }
        insert_fastq_file_suffix: {
            help: "Suffix to add to the name of the file with the insert reads (in case downstream software has name requirements)",
            type: "String",
            category: "optional"
        }
        barcode_fastq_header_suffix: {
            help: "Suffix to add to the end of the fastq headers for the barcode reads (in case downstream software has fastq header requirements)",
            type: "String",
            category: "optional"
        }
        insert_fastq_header_suffix: {
            help: "Suffix to add to the end of the fastq headers for the insert reads (in case downstream software has fastq header requirements)",
            type: "String",
            category: "optional"
        }
        fastqc_adapter: {
            help: "Adapter that can be passed to fastqc with the --adapters options",
            type: "String",
            category: "optional"
        }
        fastqc_limits: {
            help: "Adapter that can be passed to fastqc with the --limits option",
            type: "File",
            category: "optional"
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
        combined_statistics: {
            type: "File",
            help: "A csv with the trimming and alignment statistics", 
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
        fastqc_reports: {
            type: "Array[File]",
            help: "Fastqc output report for insert", 
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
    }


    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    if (defined(trimmer_parameters)){
        call TrimmingTasks.Trimmer {
            input:
                input_cram_bam      = input_file,
                parameters  = select_first([trimmer_parameters]),
                base_file_name      = base_file_name,
                docker              = global.trimmer_docker,
                preemptible_tries   = preemptible_tries,
                monitoring_script   = monitoring_script,  # !FileCoercion
                no_address          = no_address,
                cpus                = cpu
        }
        File trimmer_output = Trimmer.trimmed_ucram_list[0]
    }

    File trimmed_file = select_first([trimmer_output, input_file])

    call SortTasks.Demux as ConvertToFastq {
        input:
            input_file          = trimmed_file,
            base_file_name      = base_file_name,
            demux_extra_args    = demux_extra_args,
            output_format       = "fastq",
            docker              = global.sorter_docker,
            preemptible_tries   = preemptible_tries,
            monitoring_script   = monitoring_script,  # !FileCoercion
            cpu                 = cpu
    }

    call SingleCellTasks.CreateSyntheticPairedEnd {
        input:
            input_fastq             = select_first([ConvertToFastq.output_fastq]),
            base_file_name          = base_file_name,
            barcode_fastq_file_suffix = barcode_fastq_file_suffix,
            insert_fastq_file_suffix  = insert_fastq_file_suffix,
            barcode_fastq_header_suffix = barcode_fastq_header_suffix,
            insert_fastq_header_suffix = insert_fastq_header_suffix,
            preemptible_tries       = preemptible_tries,
            monitoring_script       = monitoring_script,  # !FileCoercion
            cpu                     = cpu,
            docker                  = global.pigz_docker
    }

    RuntimeParams fastqc_runtime_params = {}
    call QCTasks.FastQC {
        input:
            input_fastq             = [CreateSyntheticPairedEnd.output_insert_fastq],
            adapter_5p              = fastqc_adapter,
            limits                  = fastqc_limits,
            fastqc_runtime_params   = fastqc_runtime_params,
            docker                  = global.fastqc_docker,
            monitoring_script       = monitoring_script #!FileCoercion
    }

    if (defined(downstream_analysis)) {
        String analysis_type = select_first([downstream_analysis])
        if(analysis_type == "star_solo"){
            call StarSoloWdl.StarSoloWorkflow {
                input:
                    insert_fastq            = CreateSyntheticPairedEnd.output_insert_fastq,
                    barcode_fastq           = CreateSyntheticPairedEnd.output_barcodes_fastq,
                    base_file_name          = base_file_name,
                    star_solo_params        = select_first([star_solo_params]),
                    genome_generate_params  = genome_generate_params,
                    no_address              = no_address,
                    preemptible_tries       = preemptible_tries,
                    cpu                     = cpu,
                    monitoring_script_input = monitoring_script_input
            }
            File star_solo_stats_csv = StarSoloWorkflow.outputs.gathered_star_stats_csv
        }

        if (analysis_type == "star") {
            call StarAlignWorkflow.StarAlignment {
                input:
                    genome                  = star_genome,
                    genome_generate_params  = genome_generate_params,
                    star_align_extra_args   = star_align_extra_args,
                    star_align_gtf_override = star_align_gtf_override,
                    input_bams              = [trimmed_file],
                    base_file_name          = base_file_name + ".star.aln",
                    preemptible_tries       = preemptible_tries,
                    no_address              = no_address,
                    cpu                     = cpu,
                    monitoring_script_input = monitoring_script_input
            }
        }
    }

    call SingleCellTasks.CombineStatistics {
        input:
            trimmer_stats           = select_first([trimmer_stats,Trimmer.trimmer_stats, ""]),
            star_solo_stats_csv     = star_solo_stats_csv,
            star_stats_csv          = StarAlignment.star_stats,
            base_file_name          = base_file_name,
            docker                  = global.ug_vc_docker,
            monitoring_script       = monitoring_script, #!FileCoercion
            no_address              = no_address,
            preemptible_tries       = preemptible_tries,
            cpu                     = 2
    }

    output {
        File? trimmer_stats_output          = Trimmer.trimmer_stats
        File? trimmer_failure_codes_csv     = Trimmer.trimmer_failure_codes_csv
        Array[File?]? trimmer_histogram     = Trimmer.histogram
        Array[File?]? trimmer_histogram_extra = Trimmer.histogram_extra
        File output_barcodes_fastq          = CreateSyntheticPairedEnd.output_barcodes_fastq
        File output_insert_fastq            = CreateSyntheticPairedEnd.output_insert_fastq
        Array[File] fastqc_reports          = FastQC.reports_html
        File combined_statistics            = CombineStatistics.combined_statistics

        # Star solo outputs
        StarSoloOutputs? star_solo_outputs = StarSoloWorkflow.outputs

        # STAR outputs
        File? star_bam                    = StarAlignment.output_bam
        File? star_reads_per_gene_file    = StarAlignment.reads_per_gene_file
        File? star_stats                  = StarAlignment.star_stats
    }
}
