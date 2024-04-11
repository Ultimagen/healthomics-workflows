version 1.0

import "tasks/globals.wdl" as Globals
import "tasks/cnv_calling_tasks.wdl" as CnvTasks
workflow SingleSampleCnmopsReadsCount{

    input{
        String pipeline_version = "1.11" # !UnusedDeclaration

        File input_bam_file
        File input_bam_file_index
        File reference_genome      #ref-genome+idx to enable cram as input file
        File reference_genome_index
    
        Int mapq
        Array[String] ref_seq_names
        Int window_length
        String base_file_name
        Boolean? save_hdf_override
        Boolean? no_address_override
        Int? preemptible_tries_override

       # Used for running on other clouds (aws)
        File? monitoring_script_input
    }

    meta {
        description: "## single sample reads count workflow."
    }
    parameter_meta {
        pipeline_version: "Pipeline version"
        input_bam_file: "Input sample bam/cram file"
        input_bam_file_index: "Input sample bai/crai index file"
        reference_genome: "Genome fasta file associated with the CRAM file"
        reference_genome_index: "Index of the fasta file associated with the CRAM file"
        mapq: "Reads mapping-quality cutoff for reads count calculation"
        ref_seq_names: "Chromosome names for which reads counts will be calculated"
        window_length: "Window lenght for which reads counts will be calculated for"
        parallel: "Number of cpus"
        sample_name: "Sample name"
        save_hdf_override: "(OPTIONAL) Whether to save sample reads counts/cohort including sample/cnmops output data in hdf5 format. (additionally to RDS format)"
        no_address_override: "(OPTIONAL) no_address_override"
        preemptible_tries_override: "(OPTIONAL) number of preemptible tries"
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    Boolean save_hdf = select_first([save_hdf_override , false])

    call Globals.Globals as Globals
      GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call CnvTasks.CnmopsGetReadCountsFromBam {
      input:
        input_bam_file = input_bam_file,
        input_bai_file =  input_bam_file_index,
        reference_genome=reference_genome,
        reference_genome_index=reference_genome_index,
        mapq = mapq,
        ref_seq_names = ref_seq_names,
        window_length = window_length,
        base_file_name = base_file_name,
        save_hdf = save_hdf,
        docker = global.ug_vc_docker,
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
       File? out_reads_count_hdf5 = out_reads_count_hdf5_maybe
    }
}


