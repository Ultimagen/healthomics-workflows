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
# Runs the Copy Number Variation calling using cnvpytor for Ultima Genomics data.

# CHANGELOG in reverse chronological order

import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks

workflow SingleSampleCNVpytorCalling {

    input {
        String pipeline_version = "1.27.3" # !UnusedDeclaration

        String base_file_name
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        Array[String] ref_seq_names
        Array[Int] window_lengths
        
        Boolean? no_address_override
        Int? preemptible_tries_override
        File? monitoring_script_input

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv prefix(input_bam_file_index) == input_bam_file
        #@wv suffix(input_bam_file) in {".bam", ".cram"}

        #@wv reference_genome == prefix(reference_genome_index)
        #@wv suffix(reference_genome) in {'.fasta', '.fa', '.fna'}
        #@wv suffix(reference_genome_index) == '.fai'

        #@wv len(ref_seq_names) > 0
        #@wv len(window_lengths) > 0

    }

    meta {
        description: "Runs single sample germline CNV calling workflow based on [CNVpytor](https://github.com/abyzovlab/CNVpytor).\n The workflow takes as input a BAM/CRAM file and reference genome fasta file and generates CNV calls in TSV and VCF formats.\n The calls are generated using multiple window sizes to improve robustness of the calls."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "no_address_override",
                "preemptible_tries_override",
                "Glob.glob"
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        input_bam_file: {
            help: "Input sample BAM/CRAM file.",
            type: "File",
            category: "input_required"
         }
        input_bam_file_index: {
            help:"Input sample BAI/CRAI index file",
            type: "File",
            category: "input_required"
       }
        reference_genome: {
            help: "Genome fasta file associated with the CRAM file",
            type: "File",
            category: "ref_required"
         }
        reference_genome_index: {
            help : "Fai index of the fasta file",
            type: "File",
            category: "ref_required"
        }
        ref_seq_names: {
            help : "Chromosome names for which coverage will be calculated",
            type: "Array[String]",
            category: "param_required"
        }
        window_lengths: {
            help: "Window lengths on which the read counts will be aggregated",
            type: "Array[Int]",
            category: "param_required"
        }
        cnvpytor_cnv_calls_tsv: {
            help: "CNVpytor CNV calls in tsv format",
            type: "File",
            category: "output"
        }
        cnvpytor_cnv_calls_vcf: {
            help: "CNVpytor CNV calls in VCF format",
            type: "File",
            category: "output"
        }
        cnvpytor_cnv_calls_vcf_index: {
            help: "Index file for the CNVpytor CNV calls VCF",
            type: "File",
            category: "output"  
        }
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
     
    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers   #!FileCoercion
    
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script]) #!FileCoercion

    scatter (window_size in window_lengths) {
        call RunCNVpytor {
            input:
                sample_name = base_file_name,
                input_bam = input_bam_file,
                input_bam_index = input_bam_file_index,
                reference_fasta = reference_genome,
                reference_fasta_index = reference_genome_index,
                window_size = window_size,
                chr_list = ref_seq_names,
                docker = global.ugbio_cnv_docker,
                monitoring_script = monitoring_script,
                no_address = no_address,
                preemptible_tries = preemptible_tries
        }
    }

    call CombineCNVVcfs { 
        input:
            base_file_name = base_file_name,
            cnv_vcfs = RunCNVpytor.cnvpytor_cnv_calls_vcf,
            reference_fasta_index = reference_genome_index,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    call UGGeneralTasks.ConcatFiles as CombineTsvs {
        input:
            out_file_name = base_file_name + ".cnvpytor.cnv_calls.tsv",
            files = RunCNVpytor.cnvpytor_cnv_calls_tsv,
            docker = global.ubuntu_docker,
    }

    output {
        File cnvpytor_cnv_calls_tsv = CombineTsvs.out_merged_file
        File cnvpytor_cnv_calls_vcf = CombineCNVVcfs.output_vcf
        File cnvpytor_cnv_calls_vcf_index = CombineCNVVcfs.output_vcf_index
    }
}

task RunCNVpytor {
    input {
        String sample_name
        File input_bam
        File input_bam_index
        File reference_fasta
        File reference_fasta_index
        Int window_size
        Array[String] chr_list
        
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float input_bam_file_size = size(input_bam, "GB")
    Float reference_fasta_file_size = size(reference_fasta, "GB")
    Float additional_disk = 100
    Int disk_size = ceil(2* input_bam_file_size + reference_fasta_file_size + additional_disk )
    
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        cnvpytor -root ~{sample_name}.pytor \
            -rd ~{input_bam} \
            -chrom ~{sep=' ' chr_list} \
            -T ~{reference_fasta}
        
        cnvpytor -root ~{sample_name}.pytor \
            -his ~{window_size}

        cnvpytor -root ~{sample_name}.pytor \
            -partition ~{window_size}
        
        cnvpytor -root ~{sample_name}.pytor \
            -call ~{window_size} > ~{sample_name}.pytor.bin~{window_size}.CNVs.1based.tsv
        
        cnvpytor -root ~{sample_name}.pytor  -view ~{window_size} <<EOF
                set print_filename ~{sample_name}.~{window_size}.CNV.vcf
                print calls
                quit 
        EOF
        
        bcftools view -Oz -o ~{sample_name}.~{window_size}.CNV.vcf.gz ~{sample_name}.~{window_size}.CNV.vcf
        bcftools index -tf ~{sample_name}.~{window_size}.CNV.vcf.gz
        #making cnvpytor output coordinates 0-based
        awk 'BEGIN {OFS="\t"}
        {
            split($2, coords, ":");
            split(coords[2], range, "-");
            start = range[1] - 1;
            $2 = coords[1] ":" start "-" range[2];
            print
        }' ~{sample_name}.pytor.bin~{window_size}.CNVs.1based.tsv > ~{sample_name}.pytor.bin~{window_size}.CNVs.tsv
        
        # The following command is for debugging purposes only. it lists down the cnvpytor project content.
        cnvpytor -root ~{sample_name}.pytor -ls

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "16 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 8
    }
    output {
        File cnvpytor_cnv_calls_tsv = "~{sample_name}.pytor.bin~{window_size}.CNVs.tsv"
        File cnvpytor_cnv_calls_vcf = "~{sample_name}.~{window_size}.CNV.vcf.gz"
        File cnvpytor_cnv_calls_vcf_index = "~{sample_name}.~{window_size}.CNV.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task CombineCNVVcfs {
    input {
        String base_file_name
        Array[File] cnv_vcfs
        File reference_fasta_index
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        combine_cnmops_cnvpytor_cnv_calls concat \
            --cnvpytor_vcf ~{sep=" " cnv_vcfs} \
            --output_vcf ~{base_file_name}.cnvpytor.cnv_calls.vcf.gz \
            --fasta_index ~{reference_fasta_index}
    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "4 GB"
        disks: "local-disk 4 HDD"
        docker: docker
        noAddress: no_address
        cpu: 1
    }
    output {
        File output_vcf = "~{base_file_name}.cnvpytor.cnv_calls.vcf.gz"
        File output_vcf_index = "~{base_file_name}.cnvpytor.cnv_calls.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }    
}
