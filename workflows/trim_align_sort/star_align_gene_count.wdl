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
#   This workflow aligns reads to a reference genome using STAR.
# CHANGELOG
import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/alignment_tasks.wdl" as UGRealign


workflow StarAlignment {
    input {
        File? genome
        StarGenomeGenerateParams? genome_generate_params

        Array[File] input_bams
        String base_file_name
        String? star_align_extra_args
        File? star_align_gtf_override

        Boolean no_address
        Int preemptible_tries
        Int cpu
        Int? memory_gb

        # Used for running on other clouds (aws)
        File? monitoring_script_input

        #@wv not(defined(genome)) -> defined(genome_generate_params)
        #@wv not(defined(genome_generate_params)) -> defined(genome)
        #@wv defined(genome) -> suffix(genome) == '.zip'
        #@wv defined(genome_generate_params) -> suffix(genome_generate_params['fasta_files']) <= {'.fasta','.fa'}
        #@wv defined(genome_generate_params) -> suffix(genome_generate_params['gtf_file']) == '.gtf'
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv len(input_bams) >= 1
        #@wv suffix(input_bams) <= {".bam", ".cram", ".ucram"}
        

    }
    meta {
        description: "Aligns reads to a reference genome using STAR aligner."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: [
                "no_address",
                "preemptible_tries",
                "monitoring_script_input",
                "Globals.glob",
                "StarGenomeGenerate.disk_size",
                "StarAlign.disk_size",
                "StarAlignStats.disk_size",
                "memory_gb"
            ]
        }
    }

    parameter_meta {
        genome: {
            help: "The reference genome to align reads to. If not provided, the genome will be generated using the genome_generate_params.",
            type: "File",
            category: "input_optional"
        }
        genome_generate_params: {
            help: "Parameters for generating the reference genome.",
            type: "StarGenomeGenerateParams",
            category: "input_optional"
        }
        input_bams: {
            help: "The input BAM files to align.",
            type: "Array[File]",
            category: "input_required"
        }
        base_file_name: {
            help: "The base file name for the output files.",
            type: "String",
            category: "input_required"
        }
        star_align_extra_args: {
            help: "Extra arguments to pass to STAR aligner.",
            type: "String",
            category: "input_optional"
        }
        star_align_gtf_override: {
            help: "Override the GTF file used for STAR alignment.",
            type: "File",
            category: "input_optional"
        }
        cpu: {
            help: "The number of CPUs to use.",
            type: "Int",
            category: "input_required"
        }
        genome_zip_output: {
            help: "The output genome zip file.",
            type: "File",
            category: "output"
        }
        output_bam: {
            help: "The aligned output BAM file.",
            type: "File",
            category: "output"
        }
        reads_per_gene_file: {
            help: "The reads per gene file.",
            type: "File",
            category: "output"
        }
        star_stats: {
            help: "The STAR alignment stats file parsed into a csv format.",
            type: "File",
            category: "output"
        }
        raw_star_log_file: {
            help: "The raw STAR log file.",
            type: "File",
            category: "output"
        }
    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    if (!defined(genome)) {

        StarGenomeGenerateParams genome_generate_params_override = select_first([genome_generate_params])

        call UGRealign.StarGenomeGenerate {
            input:
                fasta_files     =   genome_generate_params_override.fasta_files,
                gtf_file        =   genome_generate_params_override.gtf_file,
                output_basename =   genome_generate_params_override.output_basename,
                extra_args      =   genome_generate_params_override.extra_args,
                docker          =   global.star_docker,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script,  # !FileCoercion
                no_address      =   no_address,
                cpu             =   cpu
        }
    }

    File genome_override = select_first([StarGenomeGenerate.genome_zip, genome])

    call UGRealign.StarAlign {
        input:
            input_bams      =   input_bams,
            genome          =   genome_override,
            base_file_name  =   base_file_name,
            gtf_override    =   star_align_gtf_override,
            extra_args      =   star_align_extra_args,
            docker          =   global.star_docker,
            preemptible_tries = preemptible_tries,
            monitoring_script = monitoring_script,  # !FileCoercion
            no_address      =   no_address,
            cpu             =   cpu,
            memory_gb       =   memory_gb
    }

    call UGRealign.StarAlignStats {
        input:
            star_log_file   =   StarAlign.star_log_file,
            base_file_name  =   base_file_name,
            docker          =   global.ug_vc_docker,
            preemptible_tries = preemptible_tries,
            monitoring_script = monitoring_script,  # !FileCoercion
            no_address      =   no_address
    }

    output {
        File genome_zip_output = genome_override
        File output_bam = StarAlign.output_bam
        File reads_per_gene_file = StarAlign.reads_per_gene_file
        File star_stats = StarAlignStats.star_stats
        File raw_star_log_file = StarAlign.star_log_file
    }
}