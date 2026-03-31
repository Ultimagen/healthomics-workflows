version 1.0
# LICENSE
#   Copyright 2026 Ultima Genomics
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

# DESCRIPTION
# SexDetermination: Determines biological sex (XX/XY) from CRAM file using chrY coverage ratio
#
# Method:
#   - Calculates mean coverage genome-wide
#   - Calculates mean coverage on chrY non-PAR region
#   - Computes ratio (chrY_coverage / genome_coverage)
#   - Threshold-based classification:
#       * ratio > 1.0 → XY (male)
#       * ratio < 0.1 → XX (female)
#       * otherwise → ambiguous

import "tasks/structs.wdl"
import "tasks/globals.wdl" as Globals

workflow SexDetermination {
    input {
        String pipeline_version = "1.29.1" # !UnusedDeclaration
        String base_file_name

        # Input CRAM file
        File input_cram
        File input_cram_index

        # Y chromosome non-PAR BED file for coverage calculation
        File chry_nonpar_bed

        # Reference genome (optional, for CRAM decompression)
        File? ref_fasta
        File? ref_fasta_index

        # Execution parameters
        Int memory_gb = 8
        Int cpus = 2
        Int? preemptible_tries_override
        Boolean? no_address_override
        File? monitoring_script_input
        String? cloud_provider_override

        # winval validations
        #@wv suffix(input_cram) == '.cram'
        #@wv suffix(input_cram_index) == '.crai'
        #@wv prefix(input_cram_index) == input_cram
        #@wv suffix(chry_nonpar_bed) == '.bed'
    }

    meta {
        description: "SexDetermination: Determine biological sex from CRAM using chrY coverage"
        author: "Ultima Genomics"
        version: "1.0"
        WDL_AID: {exclude: [
            "pipeline_version",
            "preemptible_tries_override",
            "no_address_override",
            "cloud_provider_override",
            "monitoring_script_input",
            "Global.glob"
            ]
        }
    }

    parameter_meta {
        base_file_name: {
            help: "Base name for output files",
            category: "input_required"
        }
        input_cram: {
            help: "Input CRAM file with aligned reads",
            type: "File",
            category: "input_required"
        }
        input_cram_index: {
            help: "CRAM index file (.crai)",
            type: "File",
            category: "input_required"
        }
        chry_nonpar_bed: {
            help: "BED file defining Y chromosome non-PAR regions for coverage calculation",
            type: "File",
            category: "input_required"
        }
        ref_fasta: {
            help: "Reference genome FASTA (required for CRAM decompression if not embedded)",
            type: "File",
            category: "input_optional"
        }
        ref_fasta_index: {
            help: "Reference genome FASTA index",
            type: "File",
            category: "input_optional"
        }
        memory_gb: {
            help: "Memory in GB for sex determination task",
            type: "Int",
            category: "input_optional"
        }
        cpus: {
            help: "Number of CPUs for sex determination task",
            type: "Int",
            category: "input_optional"
        }
        sex: {
            help: "Determined sex: XY (male), XX (female), or ambiguous",
            category: "output",
            type: "String"
        }
        sex_metrics: {
            help: "Sex determination metrics file with coverage values and ratio",
            category: "output",
            type: "File"
        }
    }

    Boolean no_address = true
    Int preemptibles = 1

    call Globals.Globals as Global
    GlobalVariables global = Global.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script]) #!FileCoercion

    call DetermineSex {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            chry_nonpar_bed = chry_nonpar_bed,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            base_file_name = base_file_name,
            memory_gb = memory_gb,
            cpus = cpus,
            docker = global.ugbio_core_docker,
            monitoring_script = monitoring_script
    }

    output {
        String sex = DetermineSex.sex
        File sex_metrics = DetermineSex.metrics_file
    }
}

task DetermineSex {
    input {
        File input_cram
        File input_cram_index
        File chry_nonpar_bed
        File? ref_fasta
        File? ref_fasta_index
        String base_file_name
        Int memory_gb
        Int cpus
        String docker
        File monitoring_script
    }

    String metrics_file_name = base_file_name + ".sex_determination.txt"

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        # Build reference flag if provided
        ref_flag=""
        if [ -n "~{ref_fasta}" ] && [ "~{ref_fasta}" != "" ]; then
            ref_flag="--reference ~{ref_fasta}"
        fi

        # Calculate mean coverage genome-wide
        echo "Calculating genome-wide coverage..." >&2
        wgs_mean=$(samtools depth -a $ref_flag ~{input_cram} | \
            awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
        echo "Genome-wide mean coverage: $wgs_mean" >&2

        # Calculate mean coverage on chrY non-PAR
        echo "Calculating chrY non-PAR coverage..." >&2
        y_mean=$(samtools depth -a $ref_flag -b ~{chry_nonpar_bed} ~{input_cram} | \
            awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
        echo "ChrY non-PAR mean coverage: $y_mean" >&2

        # Calculate ratio
        ratio=$(awk -v y="$y_mean" -v w="$wgs_mean" \
            'BEGIN {if (w>0) print y/w; else print 0}')
        echo "Ratio (chrY/genome): $ratio" >&2

        # Determine sex based on thresholds
        sex=$(awk -v r="$ratio" '
            BEGIN {
                if (r > 1) print "XY";
                else if (r < 0.1) print "XX";
                else print "ambiguous";
            }')

        echo "Determined sex: $sex" >&2

        # Write metrics file
        cat > ~{metrics_file_name} <<EOF
sample	genome_mean_coverage	chry_nonpar_mean_coverage	ratio	sex
~{base_file_name}	$wgs_mean	$y_mean	$ratio	$sex
EOF

        # Output sex to stdout for WDL to capture
        echo "$sex" > sex.txt
    >>>

    output {
        String sex = read_string("sex.txt")
        File metrics_file = metrics_file_name
    }

    runtime {
        docker: docker
        memory: memory_gb + " GB"
        cpu: cpus
        disks: "local-disk 100 HDD"
    }
}
