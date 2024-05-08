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
#   This workflow use STAR solo to analyze 10x libraries.
# CHANGELOG
import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/alignment_tasks.wdl" as UGRealign
import "tasks/single_cell_tasks.wdl" as SingleCellTasks


workflow StarSoloWorkflow {
    input {
        File barcode_fastq
        File insert_fastq
        String base_file_name
        
        StarSoloParams star_solo_params
        StarGenomeGenerateParams? genome_generate_params

        Boolean no_address
        Int preemptible_tries
        Int cpu

        # Used for running on other clouds (aws)
        File? monitoring_script_input
    
        #@wv not(defined(star_solo_params['genome'])) -> defined(genome_generate_params)
        #@wv not(defined(genome_generate_params)) -> defined(star_solo_params['genome'])
        #@wv defined(star_solo_params['genome']) -> suffix(star_solo_params['genome']) == '.zip'
        #@wv defined(genome_generate_params) -> suffix(genome_generate_params['fasta_files']) <= {'.fasta','.fa'}
        #@wv defined(genome_generate_params) -> suffix(genome_generate_params['gtf_file']) == '.gtf'
        #@wv star_solo_params['strand'] == 'Reverse' or star_solo_params['strand'] == 'Forward'
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv suffix(barcode_fastq) == '.fastq.gz'
        #@wv suffix(insert_fastq) == '.fastq.gz'
        

    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    
    if (!defined(star_solo_params.genome)) {

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

    File genome_override = select_first([StarGenomeGenerate.genome_zip, star_solo_params.genome])

    call SingleCellTasks.StarSolo {
        input:
            insert_fastq        = insert_fastq,
            barcode_fastq       = barcode_fastq,
            base_file_name      = base_file_name,
            star_solo_params    = star_solo_params,
            docker              = global.star_docker,
            preemptible_tries   = preemptible_tries,
            monitoring_script   = monitoring_script,  # !FileCoercion
            no_address          = no_address,
            cpu                 = cpu
    }

    call UGRealign.StarAlignStats {
        input:
            star_log_file   =   StarSolo.star_log_file,
            base_file_name  =   base_file_name,
            docker          =   global.ug_vc_docker,
            preemptible_tries = preemptible_tries,
            monitoring_script = monitoring_script,  # !FileCoercion
            no_address      =   no_address
    }

    RuntimeParams statistics_runtime_params = {}
    call SingleCellTasks.GatherStatistics as GatherStatistics {
        input:
            gene_features_stats        = StarSolo.gene_features_stats,
            gene_summary_csv           = StarSolo.gene_summary_csv,
            gene_umi_per_cell_sorted   = StarSolo.gene_umi_per_cell_sorted,
            star_log_file               = StarSolo.star_log_file,
            barcode_file                = StarSolo.barcode_file,
            base_file_name              = base_file_name,
            library_direction           = star_solo_params.library_direction,
            umi_len                     = star_solo_params.umi_length,
            statistics_runtime_params   = statistics_runtime_params,
            docker                      = global.ug_vc_docker,
            monitoring_script           = monitoring_script #!FileCoercion
    }

    StarSoloOutputs star_solo_outputs = {
            "genome_zip_output"           : genome_override,
            "output_bam"                  : StarSolo.output_bam,
            "gene_features_stats"         : StarSolo.gene_features_stats,
            "gene_summary_csv"            : StarSolo.gene_summary_csv,
            "gene_umi_per_cell_sorted"    : StarSolo.gene_umi_per_cell_sorted,
            "gene_filtered_features"      : StarSolo.gene_filtered_features,
            "gene_filtered_barcodes"      : StarSolo.gene_filtered_barcodes,
            "gene_filtered_matrix"        : StarSolo.gene_filtered_matrix,
            "star_log_file"               : StarSolo.star_log_file,
            "star_log_params_file"        : StarSolo.star_log_params_file,
            "barcode_file"                : StarSolo.barcode_file,
            "star_stats"                  : StarAlignStats.star_stats,
            "gathered_star_stats_csv"     : GatherStatistics.stats_csv
        }
    output {
        StarSoloOutputs outputs = star_solo_outputs
    }
}
