version 1.0

# Apache License 2.0
# Haplotype Sampling Workflow
# 
# This workflow performs haplotype sampling from a CRAM file using VG tools.
# It creates pangenome indices, samples haplotypes based on k-mer content,
# extracts sequences, and aligns to a reference genome.
#
# Based on vgteam/vg_wdl haplotype_sampling workflow and Jupyter notebook commands.
#
# CHANGELOG:
# 1.0.0 - Initial implementation
import "tasks/globals.wdl" as Globals

workflow HaplotypeSampling {
    input {
        # Required inputs
        Array[File] input_cram_bam_list        # Input CRAM files
        File cram_reference_fasta          # Reference FASTA used for CRAM encoding
        File? cram_reference_fasta_index   # Optional .fai index for CRAM reference
        File gbz_file                      # Pangenome GBZ index file
        File hapl_file                     # Pre-computed haplotype index file (.hapl)
        File alignment_reference_fasta     # Reference FASTA for final minimap2 alignment
        File alignment_reference_fasta_index # .fai index for the alignment reference FASTA
        String sample_name                 # Sample name for output files


        # Optional parameters
        Int kmer_length = 29               # K-mer length for KMC counting
        Int num_haplotypes = 32            # Number of haplotypes to sample
        Boolean include_reference = true   # Include reference in sampled haplotypes
        Boolean diploid_sampling = true    # Use diploid sampling strategy
        Int window_size = 50000            # Sliding window size for seqkit
        Int step_size = 50000              # Sliding window step size
        String minimap2_preset = "asm5"    # Minimap2 preset for alignment

        String pipeline_version = "1.29.2"   #!UnusedDeclaration
        # Resource parameters
        Int index_cores = 32
        Int index_mem_gb = 64
        Int map_cores = 30
        Int map_mem_gb = 64

    }

    meta {
        description: "## Haplotype Sampling Workflow\n\nThis workflow performs haplotype sampling from a CRAM file using VG pangenome tools. The workflow:\n1. Converts CRAM to FASTQ\n2. Counts k-mers using KMC\n3. Samples haplotypes based on k-mer content\n4. Extracts FASTA sequences from sampled graph\n5. Creates sliding windows for alignment\n6. Aligns to reference using minimap2\n\nNote: Pre-computed index files (GBZ, hapl) are required as inputs."
        author: "Ultima Genomics"
        WDL_AID: {
          exclude: ["pipeline_version","index_cores",
          "index_mem_gb",
            "map_cores",
            "map_mem_gb",
            "Glob.glob",
            "ConvertCramToFastq.mem_gb",
            "ConvertCramToFastq.disk_size_gb",
            "KmerCountingKMC.disk_size_gb",
            "SampleHaplotypes.disk_size_gb",
            "ExtractPathsFasta.cores",
            "ExtractPathsFasta.mem_gb",
            "ExtractPathsFasta.disk_size_gb",
            "SlidingWindows.cores",
            "SlidingWindows.mem_gb",
            "SlidingWindows.disk_size_gb",
            "Minimap2Align.disk_size_gb"
          ]}

    }

    parameter_meta {
        input_cram_bam_list: {
            help: "Input CRAM files containing sequencing reads",
            type: "Array[File]",
            category: "param_required"
        }
        cram_reference_fasta: {
            help: "Reference FASTA file used for CRAM encoding/decoding",
            type: "File",
            category: "param_required"
        }
        cram_reference_fasta_index: {
            help: ".fai index file for the CRAM reference FASTA",
            type: "File",
            category: "param_required"
        }
        gbz_file: {
            help: "Pangenome GBZ index file (e.g., hprc-v1.1-mc-grch38.gbz)",
            type: "File",
            category: "param_required"
        }
        alignment_reference_fasta: {
            help: "Reference FASTA for final minimap2 alignment (e.g., hg38.primary.fa)",
            type: "File",
            category: "param_required"
        }
        alignment_reference_fasta_index: {
            help: ".fai index file for the alignment reference FASTA",
            type: "File",
            category: "param_required"
        }
        hapl_file: {
            help: "Pre-computed haplotype index file (.hapl) corresponding to the GBZ",
            type: "File",
            category: "param_required"
        }

        sample_name: {
            help: "Sample name used for naming output files",
            type: "String",
            category: "param_required"
        }
        kmer_length: {
            help: "K-mer length for KMC counting (default: 29)",
            type: "Int",
            category: "param_optional"
        }
        num_haplotypes: {
            help: "Number of haplotypes to sample (default: 32)",
            type: "Int",
            category: "param_optional"
        }
        include_reference: {
            help: "Include reference in sampled haplotypes (default: true)",
            type: "Boolean",
            category: "param_optional"
        }
        diploid_sampling: {
            help: "Use diploid sampling strategy (default: true)",
            type: "Boolean",
            category: "param_optional"
        }
        window_size: {
            help: "Sliding window size for seqkit (default: 50000)",
            type: "Int",
            category: "param_optional"
        }
        step_size: {
            help: "Sliding window step size for seqkit (default: 50000)",
            type: "Int",
            category: "param_optional"
        }
        minimap2_preset: {
            help: "Minimap2 preset for alignment (default: asm5)",
            type: "String",
            category: "param_optional"
        }
        output_cram: {
            help: "Final CRAM file containing aligned haplotype sequences with embedded reference",
            type: "File",
            category: "output"
        }
        output_cram_index: {
            help: "Index file (.crai) for the output CRAM file",
            type: "File",
            category: "output"
        }
    }

    # Default docker images
    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers
    String vg_docker = global.giraffe_docker
    String samtools_docker = global.giraffe_docker
    String kmc_docker = global.giraffe_docker
    String seqkit_docker = global.giraffe_docker
    String minimap2_docker = global.giraffe_docker
    String monitoring_script = global.monitoring_script
    # Step 1: Convert CRAM to FASTQ
    call ConvertCramToFastq {
        input:
            input_cram_bam_list = input_cram_bam_list,
            reference_fasta = cram_reference_fasta,
            reference_fasta_index = cram_reference_fasta_index,
            sample_name = sample_name,
            cores = index_cores,
            samtools_docker = samtools_docker,
            monitoring_script = monitoring_script #!FileCoercion
    }

    # Step 2: K-mer counting with KMC
    call KmerCountingKMC {
        input:
            fastq_file = ConvertCramToFastq.fastq_file,
            sample_name = sample_name,
            kmer_length = kmer_length,
            mem_gb = 128,
            cores = 16,
            kmc_docker = kmc_docker,
            monitoring_script = monitoring_script   #!FileCoercion
    }

    # Step 3: Sample haplotypes
    call SampleHaplotypes {
        input:
            gbz_file = gbz_file,
            hapl_file = hapl_file,
            kff_file = KmerCountingKMC.kff_file,
            sample_name = sample_name,
            num_haplotypes = num_haplotypes,
            include_reference = include_reference,
            diploid_sampling = diploid_sampling,
            cores = 16,
            mem_gb = index_mem_gb,
            vg_docker = vg_docker,
            monitoring_script = monitoring_script   #!FileCoercion
    }

    # Step 4: Extract FASTA paths from sampled GBZ
    call ExtractPathsFasta {
        input:
            gbz_file = SampleHaplotypes.sampled_gbz,
            sample_name = sample_name + ".sampled",
            vg_docker = vg_docker,
            monitoring_script = monitoring_script    #!FileCoercion
    }

    # Step 5: Create sliding windows
    call SlidingWindows {
        input:
            input_fasta = ExtractPathsFasta.fasta_file,
            sample_name = sample_name + ".sampled",
            window_size = window_size,
            step_size = step_size,
            seqkit_docker = seqkit_docker,
            monitoring_script = monitoring_script    #!FileCoercion
    }

    # Step 6: Align to reference with minimap2 and convert to CRAM
    call Minimap2Align {
        input:
            query_fasta = SlidingWindows.windowed_fasta,
            reference_fasta = alignment_reference_fasta,
            reference_fasta_index = alignment_reference_fasta_index,
            sample_name = sample_name,
            preset = minimap2_preset,
            cores = map_cores,
            mem_gb = map_mem_gb,
            minimap2_docker = minimap2_docker,
            monitoring_script = monitoring_script    #!FileCoercion
    }

    output {
        # Final alignment (CRAM with embedded reference)
        File output_cram = Minimap2Align.output_cram
        File output_cram_index = Minimap2Align.output_cram_index
    }
}

# ============================================================================
# TASK DEFINITIONS
# ============================================================================

task ConvertCramToFastq {
    input {
        Array[File] input_cram_bam_list
        File reference_fasta
        File? reference_fasta_index
        File monitoring_script
        String sample_name
        Int cores = 32
        Int mem_gb = 16
        Int disk_size_gb = 200
        String samtools_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} > monitoring.log &

        output_fastq=~{sample_name}.fastq
        mkdir -p fastq_parts
        i=0

        for cram in ~{sep=' ' input_cram_bam_list}; do
            part_fastq="fastq_parts/part_${i}.fastq"
            samtools fastq \
                -0 "$part_fastq" \
                -@ ~{cores} \
                --reference ~{reference_fasta} \
                "$cram"
            i=$((i+1))
        done

        cat fastq_parts/*.fastq > "$output_fastq"
        rm -rf fastq_parts
    >>>

    output {
        File fastq_file = "~{sample_name}.fastq"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: samtools_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}

task KmerCountingKMC {
    input {
        File fastq_file
        File monitoring_script
        String sample_name
        Int kmer_length = 29
        Int mem_gb = 128
        Int cores = 16
        Int disk_size_gb = 200
        String kmc_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Create temp directory for KMC
        mkdir -p kmc_tmp

        # Decompress fastq if gzipped (KMC needs uncompressed for kff output)
        INPUT_FILE=~{fastq_file}

        # Run KMC
        kmc \
            -k~{kmer_length} \
            -m~{mem_gb} \
            -okff \
            -t~{cores} \
            -hp \
            $INPUT_FILE \
            ~{sample_name} \
            kmc_tmp
    >>>

    output {
        File kff_file = "~{sample_name}.kff"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: kmc_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}

task SampleHaplotypes {
    input {
        File gbz_file
        File hapl_file
        File kff_file
        File monitoring_script
        String sample_name
        Int num_haplotypes = 32
        Boolean include_reference = true
        Boolean diploid_sampling = true
        Int cores = 16
        Int mem_gb = 64
        Int disk_size_gb = 100
        String vg_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        vg haplotypes \
            -v 2 \
            -t ~{cores} \
            ~{if include_reference then "--include-reference" else ""} \
            ~{if diploid_sampling then "--diploid-sampling" else ""} \
            --num-haplotypes ~{num_haplotypes} \
            -i ~{hapl_file} \
            -k ~{kff_file} \
            -g ~{sample_name}.sampled.gbz \
            ~{gbz_file}
    >>>

    output {
        File sampled_gbz = "~{sample_name}.sampled.gbz"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: vg_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}

task ExtractPathsFasta {
    input {
        File gbz_file
        File monitoring_script
        String sample_name
        Int cores = 4
        Int mem_gb = 32
        Int disk_size_gb = 100
        String vg_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} > monitoring.log &

        vg paths \
            -F \
            -x ~{gbz_file} \
            > ~{sample_name}.fasta
    >>>

    output {
        File fasta_file = "~{sample_name}.fasta"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: vg_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}

task SlidingWindows {
    input {
        File input_fasta
        File monitoring_script
        String sample_name
        Int window_size = 50000
        Int step_size = 50000
        Int cores = 4
        Int mem_gb = 16
        Int disk_size_gb = 50
        String seqkit_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        seqkit sliding \
            -W ~{window_size} \
            -s ~{step_size} \
            ~{input_fasta} \
            -o ~{sample_name}.windowed.fasta
    >>>

    output {
        File windowed_fasta = "~{sample_name}.windowed.fasta"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: seqkit_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}

task Minimap2Align {
    input {
        File query_fasta
        File reference_fasta
        File reference_fasta_index
        File monitoring_script
        String sample_name
        String preset = "asm5"
        Int cores = 30
        Int mem_gb = 64
        Int disk_size_gb = 100
        String minimap2_docker
    }

    command <<<
        set -euxo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Run minimap2 and pipe to samtools for sorting and CRAM conversion
        # Using embedded reference for self-contained CRAM file
        minimap2 \
            -x ~{preset} \
            -t ~{cores} \
            -a \
            ~{reference_fasta} \
            ~{query_fasta} \
        | samtools sort \
            -@ ~{cores} \
            -O CRAM \
            --reference ~{reference_fasta} \
            --output-fmt-option embed_ref=1 \
            -o ~{sample_name}.haplotypes.cram \
            -

        # Create index for the CRAM file
        samtools index ~{sample_name}.haplotypes.cram
    >>>

    output {
        File output_cram = "~{sample_name}.haplotypes.cram"
        File output_cram_index = "~{sample_name}.haplotypes.cram.crai"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: minimap2_docker
        cpu: cores
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
    }
}
