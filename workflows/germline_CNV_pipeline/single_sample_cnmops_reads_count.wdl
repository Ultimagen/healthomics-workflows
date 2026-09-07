version 1.0

import "tasks/globals.wdl" as Globals
import "tasks/cnv_calling_tasks.wdl" as CnvTasks
workflow SingleSampleCnmopsReadsCount{

    input{
        String pipeline_version = "1.35.1" # !UnusedDeclaration

        File input_bam_file
        File input_bam_file_index
        File reference_genome      #ref-genome+idx to enable cram as input file
        File reference_genome_index
    
        Int mapq
        # Either genome_windows or ref_seq_names + window_length must be given
        File? genome_windows
        Array[String]? ref_seq_names
        Int? window_length
        String base_file_name
        Boolean? save_hdf_override
        Boolean? no_address_override
        Int? preemptible_tries_override

       # Used for running on other clouds (aws)
        File? monitoring_script_input
    }
    meta {
        description: "Runs single sample germline CNV calling workflow based on [cn.mops](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)\n\nThe pipeline uses a given cohort's coverage profile for normalization.\n\nThe pipeline can recieve one of the following options as input:\n\n&nbsp;&nbsp;1. Input CRAM/BAM file. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_template.json\n\n&nbsp;&nbsp;2. A rds file which stores a GenomicRanges object with coverage collected in the same windows as the given cohort. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_skip_reads_count_template.json\n\n&nbsp;&nbsp;3. A BedGraph holding the coverage per location. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_input_bedGraph_template.json\n\nThe pipeline calls CNVs for the given sample and filters them by length (>10,000b) and overlap with UG-CNV-LCR.\n\n<b>When Running in AWS HealthOmics this pipeline should run with [dynamic storage](https://docs.omics.ai/products/workbench/engines/parameters/aws-healthomics#storage_type-dynamic-or-static)</b>"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "Glob.glob",
                "monitoring_script_input"
                ]}
    }

    parameter_meta {
                base_file_name: {
            help: "Base name for the output file, if sample_name not provided - will also be the sample name in the VCF",
            type: "String",
            category: "input_required"
        }

        input_bam_file: {
            help: "Input sample bam/cram file",
            type: "File",
            category: "input_required"
        }
        input_bam_file_index: {
            help: "Input sample bai/crai index file",
            type: "File",
            category: "input_required"
        }
        reference_genome: {
            help: "Genome fasta file associated with the CRAM file",
            type: "File",
            category: "input_required"
        }
        reference_genome_index: {
            help: "Index of the fasta file associated with the CRAM file",
            type: "File",
            category: "input_required"
        }
        mapq: {
            help: "Reads mapping-quality cutoff for reads count calculation",
            type: "Int",
            category: "input_required"
        }
        genome_windows: {
            help: "Bed file with the windows in which reads counts will be calculated (e.g. the cohort's windows). Mutually exclusive with ref_seq_names/window_length",
            type: "File?",
            category: "input_optional"
        }
        ref_seq_names: {
            help: "Chromosome names for which reads counts will be calculated. Mutually exclusive with genome_windows",
            type: "Array[String]?",
            category: "input_optional"
        }
        window_length: {
            help: "Window lenght for which reads counts will be calculated for. Mutually exclusive with genome_windows",
            type: "Int?",
            category: "input_optional"
        }
        parallel: {
            help: "Number of cpus to use, default is set in the template",
            type: "Int",
            category: "input_advanced"
        }
        sample_name: {
            help: "Sample name to be used as the prefix for the output files.",
            type: "String",
            category: "input_advanced"
        }
        save_hdf_override: {
            help: "Whether to save sample reads counts in hdf5 format. (additionally to RDS format)",
            type: "Boolean?",
            category: "input_advanced"
        }
        preemptible_tries_override: {
            help: "Number of tries for preemptible instances. Default is 1.",
            type: "Int?",
            category: "input_advanced"
        }
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    Boolean save_hdf = select_first([save_hdf_override , false])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call CnvTasks.CnmopsGetReadCountsFromBam {
      input:
        input_bam_file = input_bam_file,
        input_bai_file =  input_bam_file_index,
        reference_genome=reference_genome,
        reference_genome_index=reference_genome_index,
        mapq = mapq,
        genome_windows = genome_windows,
        ref_seq_names = ref_seq_names,
        window_length = window_length,
        base_file_name = base_file_name,
        save_hdf = save_hdf,
        docker = global.ugbio_cnv_docker,
        preemptible_tries = preemptible_tries,
        monitoring_script = monitoring_script,
        no_address = no_address
    }
    # this is done to fix the caching issue
    if (save_hdf){
        File out_reads_count_hdf5_maybe = CnmopsGetReadCountsFromBam.out_reads_count_hdf5
    }
    output {
       File out_reads_count = CnmopsGetReadCountsFromBam.out_reads_count
       String out_sample_name = CnmopsGetReadCountsFromBam.out_sample_name
       File? out_reads_count_hdf5 = out_reads_count_hdf5_maybe
    }
}


