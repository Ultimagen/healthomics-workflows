version 1.0
# LICENSE
#   Copyright 2025 Ultima Genomics
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
# Alignment-based STR Genotype Caller WDL Workflow
# This workflow runs the alignment genotype caller on a CRAM file
# and outputs detailed alignment results, summary, and genotype calls.
#
# OUTPUT FILES:
# - detailed_csv: Per-read alignment details with scores (CSV format for analysis)
# - summary_csv: Per-locus summary statistics aggregated from detailed results (CSV)
# - genotypes_bed: Final genotype calls in BED format (for genome browser visualization)
# - genotypes_vcf: Final genotype calls in VCF format with per-allele support counts

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl"

workflow STRGenotyper {
    input {
        # Required inputs
        String base_file_name
        File cram_file
        File cram_index
        File variant_catalog
        
        # Reference files using the standard References struct
        References references
        
        # Optional parameters
        Int ref_padding = 500
        Int min_repeat = 1
        Int max_repeat = 100
        Float min_score_ratio = 0.85
        Int spanning_flank_bases = 5
        
        # Resource configuration
        Int threads = 8
        Int? memory_gb_override
        
        # Runtime parameters
        Int preemptible_tries = 1
        Boolean no_address = true
        File? monitoring_script_input
    }
    
    # Calculate disk size from input files automatically
    Int disk_gb = ceil(size(cram_file, "GB") + size(references.ref_fasta, "GB") + size(variant_catalog, "GB")) + 10
    
    # Default memory is 16GB, can be overridden
    Int memory_gb = select_first([memory_gb_override, 16])
    
    meta {
        description: "Alignment-based STR genotype caller using Smith-Waterman alignment"
        author: "Ultima Genomics"
        version: "1.0"
        WDL_AID: {
            exclude: [
                "GlobalsCall.glob",
                "no_address",
                "preemptible_tries",
                "monitoring_script_input"
            ]
        }
    }
    
    parameter_meta {
        base_file_name: {
            help: "Prefix for name of all output files",
            type: "String",
            category: "input_required"
        }
        cram_file: {
            help: "Input CRAM file for STR genotyping",
            type: "File",
            category: "input_required"
        }
        cram_index: {
            help: "CRAM index file (.crai)",
            type: "File",
            category: "input_required"
        }
        variant_catalog: {
            help: "JSON file containing STR variant catalog with locus definitions",
            type: "File",
            category: "input_required"
        }
        references: {
            help: "Reference genome files (fasta, fasta.fai, dict) as References struct",
            type: "References",
            category: "input_required"
        }
        ref_padding: {
            help: "Number of bases to extend around the STR repeat region when building auxiliary references for alignment. Larger values provide more flanking sequence context for accurate alignment.",
            type: "Int",
            category: "param_optional"
        }
        min_repeat: {
            help: "Minimum number of repeat units to include in auxiliary reference sequences",
            type: "Int",
            category: "param_optional"
        }
        max_repeat: {
            help: "Maximum number of repeat units to include in auxiliary reference sequences",
            type: "Int",
            category: "param_optional"
        }
        min_score_ratio: {
            help: "Minimum ratio of alignment score to the theoretical maximum score (read_length * match_score). Alignments below this threshold are filtered out. Range: 0.0-1.0, where 1.0 requires perfect alignment.",
            type: "Float",
            category: "param_optional"
        }
        spanning_flank_bases: {
            help: "Minimum number of bases that must align on each side of the STR repeat region for a read to be considered 'spanning' the locus",
            type: "Int",
            category: "param_optional"
        }
        threads: {
            help: "Number of threads for parallel processing",
            type: "Int",
            category: "param_optional"
        }
        memory_gb_override: {
            help: "Optional memory allocation override in GB (default: 16)",
            type: "Int",
            category: "param_advanced"
        }
        preemptible_tries: {
            help: "Number of preemptible tries before running on non-preemptible",
            type: "Int",
            category: "param_advanced"
        }
        detailed_csv: {
            help: "Detailed per-read alignment results in CSV format, containing alignment scores, repeat counts, and read metadata for each alignment",
            type: "File",
            category: "output"
        }
        summary_csv: {
            help: "Per-locus summary statistics in CSV format, aggregating alignment results across all reads for each STR locus",
            type: "File",
            category: "output"
        }
        genotypes_bed: {
            help: "Final genotype calls in BED format for visualization in genome browsers (IGV, UCSC). Contains chromosome, start, end, and genotype information",
            type: "File",
            category: "output"
        }
        genotypes_vcf: {
            help: "Final genotype calls in compressed VCF format with per-allele support counts (ADSP, ADFL). Compatible with standard VCF tools.",
            type: "File",
            category: "output"
        }
        genotypes_vcf_index: {
            help: "Tabix index for the genotypes VCF file",
            type: "File",
            category: "output"
        }
    }

    call Globals.Globals as GlobalsCall
    GlobalVariables global = GlobalsCall.global_dockers

    String monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    
    call GenotypeSTR {
        input:
            base_file_name = base_file_name,
            cram_file = cram_file,
            cram_index = cram_index,
            variant_catalog = variant_catalog,
            reference_fasta = references.ref_fasta,
            reference_fasta_index = references.ref_fasta_index,
            ref_padding = ref_padding,
            min_repeat_count = min_repeat,
            max_repeat_count = max_repeat,
            min_score_ratio = min_score_ratio,
            spanning_flank_bases = spanning_flank_bases,
            threads = threads,
            memory_gb = memory_gb,
            disk_gb = disk_gb,
            docker = global.str_genotyper_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptible_tries,
            no_address = no_address
    }
    
    output {
        File detailed_csv = GenotypeSTR.detailed_csv
        File summary_csv = GenotypeSTR.summary_csv
        File genotypes_bed = GenotypeSTR.genotypes_bed
        File genotypes_vcf = GenotypeSTR.genotypes_vcf
        File genotypes_vcf_index = GenotypeSTR.genotypes_vcf_index
    }
}

task GenotypeSTR {
    input {
        String base_file_name
        File cram_file
        File cram_index
        File variant_catalog
        File reference_fasta
        File reference_fasta_index
        
        Int ref_padding
        Int min_repeat_count
        Int max_repeat_count
        Float min_score_ratio
        Int spanning_flank_bases
        
        Int threads
        Int memory_gb
        Int disk_gb
        String docker
        String monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    
    command <<<
        set -exo pipefail
        
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        
        # Create output directory
        mkdir -p ./results
        
        # Run the alignment genotype caller with CLI arguments
        python -m alignment_str_len_caller.main \
            --cram-file ~{cram_file} \
            --cram-index ~{cram_index} \
            --reference ~{reference_fasta} \
            --reference-index ~{reference_fasta_index} \
            --variant-catalog ~{variant_catalog} \
            --ref-padding ~{ref_padding} \
            --min-repeat ~{min_repeat_count} \
            --max-repeat ~{max_repeat_count} \
            --min-score-ratio ~{min_score_ratio} \
            --spanning-flank-bases ~{spanning_flank_bases} \
            --threads ~{threads} \
            --output-dir ./results \
            --output-prefix ~{base_file_name}
    >>>
    
    output {
        File monitoring_log = "monitoring.log"
        File detailed_csv = "results/~{base_file_name}_detailed.csv"
        File summary_csv = "results/~{base_file_name}_summary.csv"
        File genotypes_bed = "results/~{base_file_name}_genotypes.bed"
        File genotypes_vcf = "results/~{base_file_name}_genotypes.vcf.gz"
        File genotypes_vcf_index = "results/~{base_file_name}_genotypes.vcf.gz.tbi"
    }
    
    runtime {
        docker: docker
        cpu: threads
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible_tries
        noAddress: no_address
        maxRetries: 1
    }
    
    meta {
        description: "Run alignment-based STR genotype caller on a single CRAM file"
    }
}
