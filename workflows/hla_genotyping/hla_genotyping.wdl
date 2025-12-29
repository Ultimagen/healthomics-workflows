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
#   HLA Genotyping workflow
#   This workflow is used to run HLA genotyping on a given CRAM or BAM file
#   The workflow will extract the sample name from the input file and run HLA-LA on the input file


import 'tasks/structs.wdl'
import 'tasks/globals.wdl' as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks

workflow HLAGenotyping {
input{
    String pipeline_version = "1.26.1" # !UnusedDeclaration
    String base_file_name

    File input_cram_bam
    File input_cram_bam_index
    File graphs_files_tar

    References references

    Int? preemptible_tries_override
    Boolean? no_address_override
    File? monitoring_script_input
    # Used for running on other clouds (aws)
    String? cloud_provider_override

    # winval validations
    #@wv suffix(input_cram_bam) in {".bam", ".cram"}
    #@wv suffix(input_cram_bam_index) in {".bai", ".crai"}
    #@wv defined(references) -> len(references) == 3
    #@wv defined(references) -> defined(references['ref_fasta'])
    #@wv defined(references) -> defined(references['ref_fasta_index'])
    #@wv prefix(input_cram_bam_index) == input_cram_bam
    }
    meta {
        description: "HLA Genotyping"
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "no_address",
            "preemptible_tries_override",
            "no_address_override",
            "cloud_provider_override",
            "monitoring_script_input",
            "Globals.glob"
    ]}
    }
    parameter_meta {
        base_file_name: {
        help: "Base file name for the output files (to be used as the prefix)",
        type: "string",
        category: "required"
        }
        input_cram_bam: {
            type: "File",
            help: "Input CRAM or BAM file for annalysing HLA genotyping",
            category: "required"
        }
        input_cram_bam_index: {
            type: "File",
            help: "Input CRAM or BAM index file for annalysing HLA genotyping",
            category: "required"
        }
        graphs_files_tar: {
            type: "File",
            help: "HLA-LA graphs files tar",
            category: "optional"
        }
        references: {
            type: "References",
            help: "Reference files: fasta, dict and fai, recommended value set in the template",
            category: "required"
        }
        cloud_provider_override: {
            type: "String",
            help: "Cloud provider to use for the workflow. Currently supported: aws, gcp default: gcp",
            category: "optional"
        }
        output_hla: {
            type: "File",
            help: "HLA genotyping output file",
            category: "output"
        }

    }
    Int preemptibles = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])

    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call UGGeneralTasks.ExtractSampleNameFlowOrder as ExtractSampleName {
            input:
                input_bam = input_cram_bam,
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                docker = global.broad_gatk_docker,
                references = references,
                no_address = no_address,
                cloud_provider_override = cloud_provider_override
    }

    call HLALAGenotyping {
        input:
            base_file_name = base_file_name,
            input_cram = input_cram_bam,
            reference = references,
            input_cram_index = input_cram_bam_index,
            sample_name = ExtractSampleName.sample_name,
            docker = global.hla_la_docker,
            graphs_files_tar = graphs_files_tar,
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script,
            no_address = no_address

    }
    output {
        File output_hla = HLALAGenotyping.output_hla
    }
}

task HLALAGenotyping {
    input {
        String base_file_name
        File monitoring_script
        File input_cram
        File input_cram_index
        File graphs_files_tar
        References reference
        String sample_name
        Int preemptible_tries
        String docker
        Boolean no_address
    }
    Int disk_size = ceil(size(input_cram,"GB") + 10*size(graphs_files_tar,"GB") + size(reference.ref_fasta,"GB")) + 40
    command {
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        # Extract the base directory and create necessary directories
        base_dir=$(echo ~{graphs_files_tar} | cut -d'/' -f2)

        mkdir -p graphs
        tar -xzvf ~{graphs_files_tar} -C graphs

        mkdir -p working

        /usr/local/bin/HLA-LA/src/HLA-LA.pl \
        --BAM ~{input_cram} \
        --workingDir working/ \
        --customGraphDir graphs/ \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ~{sample_name} \
        --maxThreads 7 \
        --samtools_T ~{reference.ref_fasta} \
        --longReads ultimagen

        mv working/~{sample_name}/hla/R1_bestguess_G.txt R1_bestguess_G_~{base_file_name}.txt
    }
    runtime {
        preemptible: preemptible_tries
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File output_hla = "R1_bestguess_G_~{base_file_name}.txt"
    }
}