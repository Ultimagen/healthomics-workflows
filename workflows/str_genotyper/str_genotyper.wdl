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
import "tasks/genome_resources.wdl" as GenomeResourcesLib

workflow STRGenotyper {
    input {
        String pipeline_version = "1.35.1" # !UnusedDeclaration
        # Required inputs
        String base_file_name
        File cram_file
        File cram_index

        # Genome resources
        String reference_genome = "hg38"
        File variant_catalog  # Optional override for custom catalog
        
        # Optional parameters
        Int ref_padding = 500
        Int min_repeat = 1
        Int max_repeat = 100
        Float min_score_ratio = 0.85
        Int spanning_flank_bases = 10
        Int min_mapping_quality = 1
        
        # Resource configuration
        Int threads = 2
        Int? memory_gb_override
        
        # Output options (set to false to skip CSV outputs for large catalogs)
        Boolean output_detailed_csv = true
        Boolean output_summary_csv = true
        
        # Genotype calling options
        Boolean haploid = false
        Boolean report_micro_alleles = false
        Float micro_allele_consensus_ratio = 0.8
        Int micro_allele_min_reads = 10

        # Runtime parameters
        Int preemptible_tries = 1
        Boolean no_address = true
        File? monitoring_script_input
    }

    meta {
        description: "Alignment-based STR genotype caller using Smith-Waterman alignment"
        author: "Ultima Genomics"
        version: "1.0"
        WDL_AID: {
            exclude: [
                "pipeline_version",
                "GlobalsCall.glob",
                "GenomeResourcesCall",
                "no_address",
                "preemptible_tries",
                "monitoring_script_input"
            ]
        }
    }

    #@wv reference_genome in {"hg38", "hg38_nist_v3_with_decoy"}
    
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
        reference_genome: {
            help: "Reference genome name (supported: 'hg38', 'hg38_nist_v3_with_decoy'). Automatically loads genome-specific reference files.",
            type: "input_required",
            category: "input_required"
        }
        variant_catalog: {
            help: "Variant catalog (json). Example: https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/variant_catalogs/variant_catalog_with_offtargets.GRCh38.json",
            type: "File",
            category: "input_optional"
        }
        ref_padding: {
            help: "Number of bases to extend around the STR repeat region when building auxiliary references for alignment. Larger values provide more flanking sequence context for accurate alignment.",
            type: "Int",
            category: "input_optional"
        }
        min_repeat: {
            help: "Minimum number of repeat units to include in auxiliary reference sequences",
            type: "Int",
            category: "input_optional"
        }
        max_repeat: {
            help: "Maximum number of repeat units to include in auxiliary reference sequences",
            type: "Int",
            category: "input_optional"
        }
        min_score_ratio: {
            help: "Minimum ratio of alignment score to the theoretical maximum score (read_length * match_score). Alignments below this threshold are filtered out. Range: 0.0-1.0, where 1.0 requires perfect alignment.",
            type: "Float",
            category: "input_optional"
        }
        spanning_flank_bases: {
            help: "Minimum number of bases that must align on each side of the STR repeat region for a read to be considered 'spanning' the locus",
            type: "Int",
            category: "input_optional"
        }
        min_mapping_quality: {
            help: "Minimum mapping quality for reads to be included in analysis",
            type: "Int",
            category: "input_optional"
        }
        threads: {
            help: "Number of threads for parallel processing",
            type: "Int",
            category: "input_advanced"
        }
        memory_gb_override : {
            help: "Optional memory allocation override in GB (default: 4)",
            type: "Int",
            category: "input_advanced"
        }
        output_detailed_csv: {
            help: "Whether to output detailed per-read CSV file. Set to false for large catalogs to reduce I/O.",
            type: "Boolean",
            category: "input_advanced"
        }
        output_summary_csv: {
            help: "Whether to output summary per-locus CSV file. Set to false for large catalogs to reduce I/O.",
            type: "Boolean",
            category: "input_optional"
        }
        haploid: {
            help: "Enable haploid mode: report single allele instead of diploid pairs. Use for X/Y chromosomes in males or haploid organisms.",
            type: "Boolean",
            category: "input_optional"
        }
        report_micro_alleles: {
            help: "Report micro-alleles (e.g. 15.3) for haploid loci when a partial-repeat insertion is present. REPCN stays the integer floor; the micro-allele decimal appears in the genotype (GT) and a new RCMA VCF/BED field. No-op for diploid loci. Default: off.",
            type: "Boolean",
            category: "input_optional"
        }
        micro_allele_consensus_ratio: {
            help: "Minimum fraction of supporting spanning reads required to report a micro-allele decimal. Only used when report_micro_alleles is true. Range 0.0-1.0.",
            type: "Float",
            category: "input_optional"
        }
        micro_allele_min_reads: {
            help: "Minimum number of spanning reads supporting the consensus tract length required to report a micro-allele. Guards against low-coverage indel artifacts. Only used when report_micro_alleles is true. Default: 10.",
            type: "Int",
            category: "input_optional"
        }
        preemptible_tries: {
            help: "Number of preemptible tries before running on non-preemptible",
            type: "Int",
            category: "input_advanced"
        }
        monitoring_script_input: {
            help: "Monitoring script override for AWS HealthOmics workflow templates multi-region support",
            type: "File",
            category: "input_advanced"
        }
        detailed_csv_files: {
            help: "Detailed per-read alignment results in CSV format, containing alignment scores, repeat counts, and read metadata for each alignment. Empty array if output_detailed_csv=false.",
            type: "Array[File]",
            category: "output"
        }
        summary_csv_files: {
            help: "Per-locus summary statistics in CSV format, aggregating alignment results across all reads for each STR locus. Empty array if output_summary_csv=false.",
            type: "Array[File]",
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

    call GenomeResourcesLib.GenomeResourcesWorkflow as GenomeResourcesCall

    References references = object {
        ref_fasta: GenomeResourcesCall.resources[reference_genome].ref_fasta,
        ref_fasta_index: GenomeResourcesCall.resources[reference_genome].ref_fasta_index,
        ref_dict: GenomeResourcesCall.resources[reference_genome].ref_dict
    }



    # Default memory is 4GB, can be overridden
    Int memory_gb = select_first([memory_gb_override, 4])
    
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
            min_mapping_quality = min_mapping_quality,
            output_detailed_csv = output_detailed_csv,
            output_summary_csv = output_summary_csv,
            haploid = haploid,
            report_micro_alleles = report_micro_alleles,
            micro_allele_consensus_ratio = micro_allele_consensus_ratio,
            micro_allele_min_reads = micro_allele_min_reads,
            threads = threads,
            memory_gb = memory_gb,
            docker = global.str_genotyper_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptible_tries,
            no_address = no_address
    }
    
    output {
        Array[File] detailed_csv_files = GenotypeSTR.detailed_csv_files
        Array[File] summary_csv_files = GenotypeSTR.summary_csv_files
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
        Int min_mapping_quality
        
        Boolean output_detailed_csv
        Boolean output_summary_csv
        Boolean haploid
        Boolean report_micro_alleles
        Float micro_allele_consensus_ratio
        Int micro_allele_min_reads

        Int threads
        Int memory_gb
        String docker
        String monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    # Calculate disk size from input files automatically
    Int disk_gb = ceil(size(cram_file, "GB") + size(reference_fasta, "GB") + size(variant_catalog, "GB")) + 10

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
            --min-mapping-quality ~{min_mapping_quality} \
            --threads ~{threads} \
            --output-dir ./results \
            --output-prefix ~{base_file_name} \
            ~{if output_detailed_csv then "" else "--skip-detailed-csv"} \
            ~{if output_summary_csv then "" else "--skip-summary-csv"} \
            ~{if haploid then "--haploid" else ""} \
            ~{if report_micro_alleles then "--report-micro-alleles --micro-allele-consensus-ratio " + micro_allele_consensus_ratio + " --micro-allele-min-reads " + micro_allele_min_reads else ""}
    >>>
    
    output {
        File monitoring_log = "monitoring.log"
        Array[File] detailed_csv_files = glob("results/*_detailed.csv")
        Array[File] summary_csv_files = glob("results/*_summary.csv")
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
