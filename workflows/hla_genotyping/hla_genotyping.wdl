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
#   This workflow is used to run HLA genotyping on a given CRAM or BAM file (by T1K or HLA-LA). 
#   T1K is the recommended tool


import 'tasks/structs.wdl'
import 'tasks/globals.wdl' as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/alignment_tasks.wdl" as UGAlignment
import "tasks/genome_resources.wdl" as GenomeResourcesLib

workflow HLAGenotyping {
input{
    String pipeline_version = "1.33.0" # !UnusedDeclaration
    String base_file_name

    File input_cram_bam
    File input_cram_bam_index
    File? graphs_files_tar

    String hla_genotyping_tool = "T1K"  # Options: "HLA-LA", "T1K"
    File? t1k_index_tar

    # Genome type selector
    String reference_genome = "hg38"

    Int? preemptible_tries_override
    Boolean? no_address_override
    File? monitoring_script_input
    # Used for running on other clouds (aws)
    String? cloud_provider_override

    # winval validations
    #@wv suffix(input_cram_bam) in {".bam", ".cram"}
    #@wv suffix(input_cram_bam_index) in {".bai", ".crai"}
    #@wv prefix(input_cram_bam_index) == input_cram_bam
    #@wv reference_genome in {"hg38", "hg38_nist_v3_with_decoy"}
    #@wv hla_genotyping_tool in {"HLA-LA", "T1K"}
    #@wv hla_genotyping_tool == "HLA-LA" -> defined(graphs_files_tar)
    #@wv hla_genotyping_tool == "T1K" -> defined(t1k_index_tar)
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
            "Glob.glob",
            "CreateReferenceCache.disk_size",
            "CreateReferenceCache.cache_populate_script_path"

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
            help: "HLA-LA graphs files tar (required if using HLA-LA)",
            category: "optional"
        }
        hla_genotyping_tool: {
            type: "String",
            help: "HLA genotyping tool to use. Options: 'HLA-LA' or 'T1K'",
            category: "required"
        }
        t1k_index_tar: {
            type: "File",
            help: "T1K index tar.gz containing hlaidx/ and kiridx/ directories with all index files (_seq.fa and _coord.fa). Required if using T1K.",
            category: "optional"
        }
        reference_genome: {
            type: "String",
            help: "Genome type selector (hg38 or hg38_nist_v3_with_decoy). Determines which reference files to use.",
            category: "required"
        }
        cloud_provider_override: {
            type: "String",
            help: "Cloud provider to use for the workflow. Currently supported: aws, gcp default: gcp",
            category: "input_optional"
        }
        monitoring_script_input: {
            type: "File",
            help: "Monitoring script override for AWS HealthOmics workflow templates multi-region support",
            category: "input_optional"
        }
        output_hla: {
            type: "File",
            help: "HLA genotyping output file",
            category: "output"
        }
        output_kir: {
            type: "File",
            help: "KIR genotyping output file",
            category: "output"
        }

    }
    Int preemptibles = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    # Get genome resources based on reference_genome
    call GenomeResourcesLib.GenomeResourcesWorkflow as GenomeResourcesWorkflow

    References references = object {
        ref_fasta: GenomeResourcesWorkflow.resources[reference_genome].ref_fasta,
        ref_fasta_index: GenomeResourcesWorkflow.resources[reference_genome].ref_fasta_index,
        ref_dict: GenomeResourcesWorkflow.resources[reference_genome].ref_dict
    }

    call UGAlignment.CreateReferenceCache {
        input:
            references = [references.ref_fasta],
            preemptible_tries = preemptibles,
            docker = global.ugbio_core_docker,
            dummy_input_for_call_caching = ""
    }

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

    # HLA-LA path
    if (hla_genotyping_tool == "HLA-LA") {
        call HLALAGenotyping {
            input:
                base_file_name = base_file_name,
                input_cram = input_cram_bam,
                reference = references,
                input_cram_index = input_cram_bam_index,
                sample_name = ExtractSampleName.sample_name,
                docker = global.hla_la_docker,
                graphs_files_tar = select_first([graphs_files_tar]),
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }
    }

    # T1K path
    if (hla_genotyping_tool == "T1K") {
        call T1KHLAGenotyping {
            input:
                base_file_name = base_file_name,
                input_cram = input_cram_bam,
                reference = references,
                input_cram_index = input_cram_bam_index,
                t1k_index_tar = select_first([t1k_index_tar]),
                cache_tarball = CreateReferenceCache.cache_tarball,
                docker = global.t1k_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }

        call T1KKIRGenotyping {
                input:
                    base_file_name = base_file_name,
                    input_cram = input_cram_bam,
                    reference = references,
                    input_cram_index = input_cram_bam_index,
                    t1k_index_tar = select_first([t1k_index_tar]),
                    cache_tarball = CreateReferenceCache.cache_tarball,
                    docker = global.t1k_docker,
                    preemptible_tries = preemptibles,
                    monitoring_script = monitoring_script,
                    no_address = no_address
        }
    }

    output {
        File output_hla = select_first([HLALAGenotyping.output_hla, T1KHLAGenotyping.output_hla])
        File? output_kir = T1KKIRGenotyping.output_kir
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

task T1KHLAGenotyping {
    input {
        String base_file_name
        File monitoring_script
        File input_cram
        File input_cram_index
        File t1k_index_tar
        References reference
        File cache_tarball
        Int preemptible_tries
        String docker
        Boolean no_address
    }

    Int disk_size = ceil(size(input_cram, "GB") + 10*size(t1k_index_tar, "GB") + size(reference.ref_fasta, "GB")) + 40

    command <<<
        set -exo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        ~{"tar -zxf "+cache_tarball}
        export REF_CACHE=cache/%2s/%2s/ 
        export REF_PATH='.' 

        # Extract T1K index files (contains both hlaidx/ and kiridx/)
        tar -xzf ~{t1k_index_tar}

        # Run T1K HLA genotyping
        perl /usr/local/bin/run-t1k \
            -b ~{input_cram} \
            -f hlaidx/_dna_seq.fa \
            -c hlaidx/_dna_coord.fa \
            --skipPostAnalysis \
            --preset hla-wgs \
            -t 4 \
            -o ~{base_file_name}_hla \
            --od ./

        # Outputs: {base_file_name}_hla_genotype.tsv, {base_file_name}_hla_allele.tsv
    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "16 GiB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File monitoring_log = "monitoring.log"
        File output_hla = "~{base_file_name}_hla_genotype.tsv"
        File output_hla_allele = "~{base_file_name}_hla_allele.tsv"
    }
}

task T1KKIRGenotyping {
    input {
        String base_file_name
        File monitoring_script
        File input_cram
        File input_cram_index
        File t1k_index_tar
        References reference
        File cache_tarball
        Int preemptible_tries
        String docker
        Boolean no_address
    }

    Int disk_size = ceil(size(input_cram, "GB") + 10*size(t1k_index_tar, "GB") + size(reference.ref_fasta, "GB")) + 40

    command <<<
        set -exo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        ~{"tar -zxf "+cache_tarball}
        export REF_CACHE=cache/%2s/%2s/ 
        export REF_PATH='.' 

        # Extract T1K index files (contains both hlaidx/ and kiridx/)
        tar -xzf ~{t1k_index_tar}

        # Run T1K KIR genotyping
        perl /usr/local/bin/run-t1k \
            -b ~{input_cram} \
            -f kiridx/_dna_seq.fa \
            -c kiridx/_dna_coord.fa \
            --preset kir-wgs \
            --skipPostAnalysis \
            -t 4 \
            -o ~{base_file_name}_kir \
            --od ./

        # Outputs: {base_file_name}_kir_genotype.tsv, {base_file_name}_kir_allele.tsv
    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "16 GiB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File monitoring_log = "monitoring.log"
        File output_kir = "~{base_file_name}_kir_genotype.tsv"
        File output_kir_allele = "~{base_file_name}_kir_allele.tsv"
    }
}