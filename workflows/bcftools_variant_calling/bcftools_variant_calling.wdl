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
# BcftoolsVariantCalling: A variant calling pipeline that processes CRAM files using bcftools
# to generate gVCF files with optional quality filtering.
#
# Main features:
#   - Automatic sex determination from CRAM (optional, enabled by default)
#   - Parallel processing of autosomes split into configurable number of shards
#   - Sex-specific handling of X and Y chromosomes with proper ploidy
#   - Optional filtering by depth and genotype quality

import "tasks/structs.wdl"
import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/genome_resources.wdl" as GenomeResourcesLib
import "sex_determination.wdl" as SexDeterminationWorkflow

workflow BcftoolsVariantCalling {
    input {
        String pipeline_version = "1.29.1" # !UnusedDeclaration
        String base_file_name

        # Main input: aligned reads in CRAM format
        File input_cram
        File input_cram_index

        # Sample metadata
        String? sex # Optional: XX or XY. If not provided, will be auto-determined when auto_determine_sex=true
        Boolean auto_determine_sex = true # Automatically determine sex from CRAM if sex not provided

        # Genomic regions for variant calling
        File autosomes_bed
        File chrx_par_bed
        File chrx_nonpar_bed
        File chry_nonpar_bed

        # Genome type for reference selection
        String reference_genome = "hg38"

        # Filtering options
        Boolean run_filtering = true

        # Execution parameters
        Int num_autosome_shards = 10
        Int memory_gb = 4
        Int cpus = 2
        Int? preemptible_tries_override
        Boolean? no_address_override
        File? monitoring_script_input
        String? cloud_provider_override

        # winval validations
        #@wv not defined(sex) or sex in {'XX', 'XY'}
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv suffix(input_cram) == '.cram'
        #@wv suffix(input_cram_index) == '.crai'
        #@wv prefix(input_cram_index) == input_cram
        #@wv suffix(autosomes_bed) == '.bed'
        #@wv suffix(chrx_par_bed) == '.bed'
        #@wv suffix(chrx_nonpar_bed) == '.bed'
        #@wv suffix(chry_nonpar_bed) == '.bed'
        #@wv reference_genome in {"hg38", "b37"}
        #@wv num_autosome_shards > 0
    }

    meta {
        description: "BcftoolsVariantCalling: bcftools-based variant calling for CRAM files"
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
            help: "Base name for all output files",
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
        sex: {
            help: "Sample sex: XX (female) or XY (male). If not provided and auto_determine_sex=true, will be determined automatically from CRAM",
            type: "String",
            category: "input_optional"
        }
        auto_determine_sex: {
            help: "Automatically determine sex from CRAM file if sex parameter not provided (default: true)",
            type: "Boolean",
            category: "input_optional"
        }
        autosomes_bed: {
            help: "BED file with autosomal targets",
            type: "File",
            category: "input_required"
        }
        chrx_par_bed: {
            help: "BED file with X chromosome PAR (pseudoautosomal) regions",
            type: "File",
            category: "input_required"
        }
        chrx_nonpar_bed: {
            help: "BED file with X chromosome non-PAR regions",
            type: "File",
            category: "input_required"
        }
        chry_nonpar_bed: {
            help: "BED file with Y chromosome non-PAR regions (used for XY samples)",
            type: "File",
            category: "input_required"
        }
        reference_genome: {
            help: "Genome type for reference selection (hg38, b37). Reference files are automatically selected based on this value.",
            type: "String",
            category: "input_required"
        }
        run_filtering: {
            help: "Whether to generate filtered gVCF using depth and GQ thresholds",
            type: "Boolean",
            category: "input_optional"
        }
        num_autosome_shards: {
            help: "Number of shards to split autosomes for parallel processing",
            type: "Int",
            category: "input_optional"
        }
        memory_gb: {
            help: "Memory in GB for variant calling tasks",
            type: "Int",
            category: "input_optional"
        }
        cpus: {
            help: "Number of CPUs for variant calling tasks",
            type: "Int",
            category: "input_optional"
        }
        final_gvcf: {
            help: "Final merged gVCF file",
            category: "output",
            type: "File"
        }
        final_gvcf_index: {
            help: "Index for final gVCF",
            category: "output",
            type: "File"
        }
        filtered_gvcf: {
            help: "Filtered gVCF with depth and GQ filtering applied",
            category: "output",
            type: "File"
        }
        filtered_gvcf_index: {
            help: "Index for filtered gVCF",
            category: "output",
            type: "File"
        }
        determined_sex: {
            help: "Final sex used for variant calling (provided or auto-determined): XX or XY",
            category: "output",
            type: "String"
        }
        sex_metrics: {
            help: "Sex determination metrics file (only present if sex was auto-determined)",
            category: "output",
            type: "File"
        }
    }

    Boolean no_address = true
    Int preemptibles = 1

    call Globals.Globals as Global
    GlobalVariables global = Global.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script]) #!FileCoercion

    # Get genome resources based on reference_genome
    call GenomeResourcesLib.GenomeResourcesWorkflow as GenomeResources

    # Construct References struct from genome resources
    References references = object {
        ref_fasta: GenomeResources.resources[reference_genome].ref_fasta,
        ref_fasta_index: GenomeResources.resources[reference_genome].ref_fasta_index,
        ref_dict: GenomeResources.resources[reference_genome].ref_dict
    }

    # Auto-determine sex if not provided
    if (!defined(sex) && auto_determine_sex) {
        call SexDeterminationWorkflow.DetermineSex {
            input:
                input_cram = input_cram,
                input_cram_index = input_cram_index,
                chry_nonpar_bed = chry_nonpar_bed,
                ref_fasta = references.ref_fasta,
                ref_fasta_index = references.ref_fasta_index,
                base_file_name = base_file_name,
                memory_gb = memory_gb,
                cpus = cpus,
                docker = global.ugbio_core_docker,
                monitoring_script = monitoring_script
        }
    }

    # Use provided sex or determined sex
    String final_sex = select_first([sex, DetermineSex.sex])

    # Split autosomes into shards for parallel processing
    call SplitAutosomesIntoBed {
        input:
            autosomes_bed = autosomes_bed,
            num_shards = num_autosome_shards,
            docker = global.ubuntu_docker,
            monitoring_script = monitoring_script
    }

    # Run variant calling for each autosome shard in parallel
    scatter (shard_num in range(num_autosome_shards)) {
        Int shard_number = shard_num + 1
        String shard_name = "autosomes.part" + (if shard_number < 10 then "0" else "") + shard_number
        File shard_bed = SplitAutosomesIntoBed.shard_beds[shard_num]

        call BcftoolsCallVariants {
            input:
                input_cram = input_cram,
                input_cram_index = input_cram_index,
                regions_bed = shard_bed,
                sample_name = base_file_name,
                part_name = shard_name,
                ploidy = 2,
                ref_fasta = references.ref_fasta,
                ref_fasta_index = references.ref_fasta_index,
                memory_gb = memory_gb,
                cpus = cpus,
                docker = global.bcftools_docker,
                monitoring_script = monitoring_script
        }
    }

    # Call variants for X PAR (diploid)
    call BcftoolsCallVariants as CallChrxPar {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            regions_bed = chrx_par_bed,
            sample_name = base_file_name,
            part_name = "chrX_PAR",
            ploidy = 2,
            ref_fasta = references.ref_fasta,
            ref_fasta_index = references.ref_fasta_index,
            memory_gb = memory_gb,
            cpus = cpus,
            docker = global.bcftools_docker,
            monitoring_script = monitoring_script
    }

    # X non-PAR region (ploidy depends on sex)
    call BcftoolsCallVariants as CallChrxNonpar {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            regions_bed = chrx_nonpar_bed,
            sample_name = base_file_name,
            part_name = "chrX_nonPAR",
            ploidy = if final_sex == "XX" then 2 else 1,
            ref_fasta = references.ref_fasta,
            ref_fasta_index = references.ref_fasta_index,
            memory_gb = memory_gb,
            cpus = cpus,
            docker = global.bcftools_docker,
            monitoring_script = monitoring_script
    }

    # Y non-PAR region (only for males)
    if (final_sex == "XY") {
        call BcftoolsCallVariants as CallChryNonparMale {
            input:
                input_cram = input_cram,
                input_cram_index = input_cram_index,
                regions_bed = chry_nonpar_bed,
                sample_name = base_file_name,
                part_name = "chrY_nonPAR",
                ploidy = 1,
                ref_fasta = references.ref_fasta,
                ref_fasta_index = references.ref_fasta_index,
                memory_gb = memory_gb,
                cpus = cpus,
                docker = global.bcftools_docker,
                monitoring_script = monitoring_script
        }
    }

    # Merge autosome shards
    call BcftoolsConcatGvcf {
        input:
            gvcf_files = BcftoolsCallVariants.gvcf,
            gvcf_files_index = BcftoolsCallVariants.gvcf_index,
            output_vcf = base_file_name + ".autosomes.g.vcf.gz",
            docker = global.bcftools_docker,
            monitoring_script = monitoring_script
    }

    # Build list of files to concatenate based on sex
    Array[File] base_gvcf_files = [
        BcftoolsConcatGvcf.merged_gvcf,
        CallChrxPar.gvcf,
        CallChrxNonpar.gvcf
    ]
    Array[File] base_gvcf_index_files = [
        BcftoolsConcatGvcf.merged_gvcf_index,
        CallChrxPar.gvcf_index,
        CallChrxNonpar.gvcf_index
    ]

    # For males, add Y chromosome; for females, skip it
    Array[File] final_gvcf_files = if final_sex == "XY" then flatten([base_gvcf_files, [select_first([CallChryNonparMale.gvcf])]]) else base_gvcf_files
    Array[File] final_gvcf_index_files = if final_sex == "XY" then flatten([base_gvcf_index_files, [select_first([CallChryNonparMale.gvcf_index])]]) else base_gvcf_index_files

    # Create final gVCF with all regions in order
    call BcftoolsConcatGvcf as ConcatFinalGvcf {
        input:
            gvcf_files = final_gvcf_files,
            gvcf_files_index = final_gvcf_index_files,
            output_vcf = base_file_name + ".final.g.vcf.gz",
            docker = global.bcftools_docker,
            monitoring_script = monitoring_script
    }

    # Optional: Filter gVCF by depth and GQ
    if (run_filtering) {
        call BcftoolsFilterGvcf {
            input:
                input_gvcf = ConcatFinalGvcf.merged_gvcf,
                input_gvcf_index = ConcatFinalGvcf.merged_gvcf_index,
                output_gvcf = base_file_name + ".filtered.g.vcf.gz",
                docker = global.bcftools_docker,
                monitoring_script = monitoring_script
        }
    }

    output {
        File final_gvcf = ConcatFinalGvcf.merged_gvcf
        File final_gvcf_index = ConcatFinalGvcf.merged_gvcf_index
        File? filtered_gvcf = BcftoolsFilterGvcf.filtered_gvcf
        File? filtered_gvcf_index = BcftoolsFilterGvcf.filtered_gvcf_index
        String determined_sex = final_sex
        File? sex_metrics = DetermineSex.metrics_file
    }
}

# Task to split autosomes into N shards
task SplitAutosomesIntoBed {
    input {
        File autosomes_bed
        Int num_shards
        String docker
        File monitoring_script
    }

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        # Count total lines and calculate per-shard count
        total=$(wc -l < ~{autosomes_bed})
        per=$((($total + ~{num_shards} - 1) / ~{num_shards}))

        # Split into N parts
        awk -v per="$per" '
        {
          f=sprintf("autosomes.part%02d.bed", int((NR-1)/per)+1)
          print >> f
        }' ~{autosomes_bed}

        # List all shard files
        ls -1 autosomes.part*.bed | sort -V
    >>>

    output {
        Array[File] shard_beds = glob("autosomes.part*.bed")
    }

    runtime {
        docker: docker
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
    }
}

# Task: Run bcftools mpileup + call for a genomic region
task BcftoolsCallVariants {
    input {
        File input_cram
        File input_cram_index
        File regions_bed
        String sample_name
        String part_name
        Int ploidy
        File ref_fasta
        File ref_fasta_index
        Int memory_gb = 4
        Int cpus = 2
        String docker
        File monitoring_script
    }

    String output_vcf = sample_name + "." + part_name + ".g.vcf.gz"

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        # Create a sample name file to force consistent sample naming
        echo "~{sample_name}" > sample_name.txt

        bcftools mpileup \
            --ignore-RG \
            ~{input_cram} \
            --fasta-ref ~{ref_fasta} \
            -R ~{regions_bed} \
            -a FORMAT/DP,FORMAT/AD \
            -q 20 -Q 20 \
            -Ou | \
        bcftools call \
            -m \
            -A \
            -a FORMAT/GQ \
            --ploidy ~{ploidy} \
            -Ou | \
        bcftools reheader --samples sample_name.txt -o - | \
        bcftools view -O z -o ~{output_vcf}

        bcftools index -t ~{output_vcf}
    >>>

    output {
        File gvcf = output_vcf
        File gvcf_index = output_vcf + ".tbi"
    }

    runtime {
        docker: docker
        memory: memory_gb + " GB"
        cpu: cpus
        disks: "local-disk 100 HDD"
    }
}

# Task: Concatenate gVCF files
task BcftoolsConcatGvcf {
    input {
        Array[File] gvcf_files
        Array[File] gvcf_files_index
        String output_vcf
        String docker
        File monitoring_script
    }

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        bcftools concat -a \
            ~{sep=' ' gvcf_files} \
            -Oz -o ~{output_vcf}

        bcftools index -t ~{output_vcf}
    >>>

    output {
        File merged_gvcf = output_vcf
        File merged_gvcf_index = output_vcf + ".tbi"
    }

    runtime {
        docker: docker
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 HDD"
    }
}

# Task: Filter gVCF by depth and GQ
task BcftoolsFilterGvcf {
    input {
        File input_gvcf
        File input_gvcf_index
        String output_gvcf
        String docker
        File monitoring_script
    }

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        bcftools filter \
            -s LowDepth -e 'FORMAT/DP<8' -m + \
            ~{input_gvcf} -Ou | \
        bcftools filter \
            -s LowGQ -e 'FORMAT/GQ<20' -m + \
            -Oz -o ~{output_gvcf}

        bcftools index -t ~{output_gvcf}
    >>>

    output {
        File filtered_gvcf = output_gvcf
        File filtered_gvcf_index = output_gvcf + ".tbi"
    }

    runtime {
        docker: docker
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 50 HDD"
    }
}

# Task: Create empty gVCF (for Y chromosome in females)
task CreateEmptyGvcf {
    input {
        String output_vcf
        String docker
        File monitoring_script
    }

    command <<<
        set -euo pipefail
        bash ~{monitoring_script} &

        echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | \
            bgzip -c > ~{output_vcf}

        tabix -p vcf ~{output_vcf}
    >>>

    output {
        File empty_gvcf = output_vcf
        File empty_gvcf_index = output_vcf + ".tbi"
    }

    runtime {
        docker: docker
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
    }
}

