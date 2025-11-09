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

workflow SingleSampleCNVpytorCalling {

    input {
        String pipeline_version = "1.23.0" # !UnusedDeclaration

        String base_file_name
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        Array[String] ref_seq_names
        Int window_length
        Int mapq
        
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

        #@wv window_length > 0
        #@wv len(ref_seq_names) > 0

    }

    meta {
        description: "Runs single sample germline CNV calling workflow based on [CNVpytor](https://github.com/abyzovlab/CNVpytor)\n</b>"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "no_address_override",
                "preemptible_tries_override",
                "Globals.glob"
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
        window_length: {
            help: "Window length on which the read counts will be aggregated",
            type: "Int",
            category: "param_required"
        }
        mapq: {
            help: "Minimum mapping quality for read to be included in the analysis",
            type: "Int",
            category: "param_required"
        }
        cnvpytor_cnv_calls_tsv: {
            help: "CNVpytor CNV calls in tsv format",
            type: "File",
            category: "output"
        }
        
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
     
    call Globals.Globals as Globals
      GlobalVariables global = Globals.global_dockers

    
    call RunCNVpytor {
        input:
        sample_name = base_file_name,
        input_bam = input_bam_file,
        input_bam_index = input_bam_file_index,
        reference_fasta = reference_genome,
        reference_fasta_index = reference_genome_index,
        window_size = window_length,
        chr_list = ref_seq_names,
        mapq = mapq,
        docker = global.ugbio_cnv_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }

    output {
        File cnvpytor_cnv_calls_tsv = RunCNVpytor.cnvpytor_cnv_calls_tsv
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
        Int mapq
        
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
        File monitoring_log = "monitoring.log"
    }
    
}
