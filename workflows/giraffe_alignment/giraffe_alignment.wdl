version 1.0

import "tasks/alignment_tasks.wdl" as UGRealign
import "haplotype_sampling.wdl" as HSampling
import "tasks/structs.wdl" as Structs  #!UnusedImport
import "tasks/globals.wdl" as Globals
import "tasks/sorting_tasks.wdl" as SortTasks

workflow GiraffeAlignment {
    input {
        String pipeline_version = "1.29.1" # !UnusedDeclaration

        Array[File] input_cram_list
        String base_file_name

        Int? preemptible_tries

        # File cache_populate_script
        Array[File] ref_files_for_tarball

        References references
        GiraffeReferences giraffe_parameters

        SorterParams sorter_params

        Boolean output_haplotypes_cram = true

        Int reads_per_split = 10000000

        Boolean? no_address_override

        String dummy_input_for_call_caching = ""

        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv len(input_cram_list) >= 1
        #@wv suffix(input_cram_list) <= {".bam", ".cram", ".ubam", ".ucram"}

        #@wv defined(references) -> len(references) == 3
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa','.fna'}
        #@wv suffix(references['ref_dict']) == '.dict'
        #@wv suffix(references['ref_fasta_index']) == '.fai'
        #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
        #@wv prefix(references['ref_dict']) == prefix(references['ref_fasta'])

        #@wv suffix(giraffe_parameters['ref_gbz']) == '.gbz'
        #@wv suffix(giraffe_parameters['ref_dist']) == '.dist'
        #@wv suffix(giraffe_parameters['ref_min']) == '.min'
        #@wv output_haplotypes_cram -> defined(giraffe_parameters['ref_gbz_for_haplotypes'])
        #@wv output_haplotypes_cram -> defined(giraffe_parameters['ref_hapl'])
        #@wv output_haplotypes_cram -> defined(giraffe_parameters['alignment_reference_fasta_for_haplotypes'])
    }

    meta {
        description: "Single-sample alignment + sort + mark-duplicates using vg giraffe only. This is a simplified variant of AlignCramOrBam with all non-giraffe aligner options removed."
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "no_address_override",
            "dummy_input_for_call_caching",
            # "cache_populate_script",
            "Glob.glob",
            "CreateReferenceCache.disk_size",
            "CreateReferenceCache.cache_populate_script_path",
            "ConvertToUbam.disk_size",
            "SplitInputBam.disk_size",
            "SplitInputCram.disk_size",
            "SamToFastqAndGiraffeAndMba.disk_size",
            "SamToFastqAndGiraffeAndMba.threads",
            "Demux.mapq_override",
            "Demux.coverage_intervals",
            "Sorter.coverage_intervals"
        ]}
    }

    parameter_meta {
        input_cram_list: {
            help: "List of input CRAM/BAM files to align using vg giraffe. All files in the list are processed together.",
            type: "Array[File]",
            category: "input_required"
        }
        base_file_name: {
            help: "Base name for all output files. Should not contain spaces, #, or comma characters.",
            type: "String",
            category: "input_required"
        }
        preemptible_tries: {
            help: "Number of preemptible retries for scatter tasks.",
            type: "Int",
            category: "param_optional"
        }
        ref_files_for_tarball: {
            help: "List of reference FASTA files used to build the reference cache.",
            type: "Array[File]",
            category: "input_required"
        }
        references: {
            help: "Linear reference genome (FASTA, FAI, DICT) used for read processing and downstream sorting/markdup.",
            type: "References",
            category: "ref_required"
        }
        giraffe_parameters: {
            help: "Giraffe graph reference bundle (GBZ/DIST/MIN/ZIPCODES/PATH_LIST) plus linear reference (FASTA/FAI/DICT).",
            type: "GiraffeReferences",
            category: "ref_required"
        }
        sorter_params: {
            help: "Parameters for sorting and marking duplicates.",
            type: "SorterParams",
            category: "input_required"
        }
        output_haplotypes_cram: {
            help: "When true, run haplotype sampling and output a haplotype CRAM.",
            type: "Boolean",
            category: "param_optional"
        }
        reads_per_split: {
            help: "Approximate number of reads per split chunk for parallel processing.",
            type: "Int",
            category: "param_optional"
        }
        no_address_override: {
            help: "Override noAddress runtime setting for tasks (defaults to true when not provided).",
            type: "Boolean",
            category: "param_optional"
        }
        dummy_input_for_call_caching: {
            help: "Dummy input for call caching stability across environments.",
            type: "String",
            category: "param_optional"
        }
        output_cram: {
            help: "Giraffe-aligned cram file.",
            type: "File",
            category: "output"
        }
        output_cram_index: {
            help: "Index file for the Giraffe-aligned cram file.",
            type: "File",
            category: "output"
        }
        sorter_stats_csv: {
            help: "Sorter stats in csv format.",
            type: "File",
            category: "output"
        }
        sorter_stats_json: {
            help: "Sorter stats in json format.",
            type: "File",
            category: "output"
        }
        ua_stats_jsons: {
            help: "Statistics file in json format from UA.",
            type: "Array[File]",
            category: "output"
        }
        output_cram_bam_list: {
            help: "Output file list after the pipeline is executed in multiple outputs mode",
            type: "Array[File]",
            category: "output"
        }
        output_cram_bam_index_list: {
            help: "Index files for the output cram files (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        sorter_stats_csv_list: {
            help: "Sorter stats in csv format (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        sorter_stats_json_list: {
            help: "Sorter stats in json format (when running in multiple outputs mode).",
            type: "Array[File]",
            category: "output"
        }
        haplotypes_cram: {
            help: "Haplotype-sampled CRAM output (when enabled).",
            type: "File",
            category: "output"
        }
        haplotypes_cram_index: {
            help: "Index for haplotype-sampled CRAM output (when enabled).",
            type: "File",
            category: "output"
        }
    }

    Int preemptibles = select_first([preemptible_tries, 1])
    Boolean no_address = select_first([no_address_override, true])

    # String detect_input_ending_from_file = if sub(input_cram_list[0], ".*\\.cram$", "is_cram") == "is_cram" then "is_cram" else "is_bam" # !StringCoercion
    # String detect_input_ending = select_first([override_input_ending, detect_input_ending_from_file])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    scatter(input_cram in input_cram_list) {
        call UGRealign.SplitCram as SplitInputCram {
            input:
                monitoring_script = global.monitoring_script, # !FileCoercion
                input_cram_bam = input_cram,
                base_file_name = base_file_name,
                reads_per_file = reads_per_split,
                docker = global.crammer_docker,
                preemptible_tries = preemptibles,
                no_address = false
        }
    }
    Array[File] split_cram = flatten(SplitInputCram.split_outputs)

    call UGRealign.CreateReferenceCache {
        input:
            references = ref_files_for_tarball,
            preemptible_tries = preemptibles,
            docker = global.ugbio_core_docker,
            dummy_input_for_call_caching = dummy_input_for_call_caching, 
    }

    scatter(split_chunk in split_cram) {
        call UGRealign.ConvertCramOrBamToUBam as ConvertToUbam {
            input:
                monitoring_script = global.monitoring_script, # !FileCoercion
                input_file = split_chunk,
                cache_tarball = CreateReferenceCache.cache_tarball,
                base_file_name = base_file_name,
                preemptible_tries = preemptibles,
                docker = global.ug_gatk_picard_docker,
                no_address = no_address
        }

        String current_name = basename(basename(basename(ConvertToUbam.unmapped_bam, ".u.bam"), ".bam"), ".ucram")

        call UGRealign.SamToFastqAndGiraffeAndMba as SamToFastqAndGiraffeAndMba {
            input:
                input_bam = ConvertToUbam.unmapped_bam,
                output_bam_basename = current_name,
                references = references,
                giraffe_references = giraffe_parameters,
                preemptible_tries = preemptibles,
                docker = global.giraffe_docker,
                monitoring_script = global.monitoring_script, # !FileCoercion
                no_address = no_address
        }
    }

    call SortTasks.Demux {
        input:
            input_cram_bam_list= SamToFastqAndGiraffeAndMba.output_bam,
            cache_tarball      = CreateReferenceCache.cache_tarball,
            base_file_name     = base_file_name,
            reference_fasta    = references.ref_fasta,
            sorter_params      = sorter_params,
            monitoring_script  = global.monitoring_script,  # !FileCoercion
            docker             = global.sorter_docker,
            preemptible_tries  = preemptibles,
            cpu_input         = 32,
    }

    call SortTasks.Sorter {
        input:
            demux_output       = Demux.demux_output,
            max_region_size    = Demux.max_region_size,
            base_file_name     = base_file_name,
            reference_fasta    = references.ref_fasta,
            sorter_params      = sorter_params,
            monitoring_script  = global.monitoring_script,  # !FileCoercion
            docker             = global.sorter_docker,
            preemptible_tries  = preemptibles,
    }

    # Sorter always outputs lists of files, even if there is only one file,
    # so we need to adjust the outputs if there is a single file
    if (length(Sorter.sorted_cram) == 1) {
        File output_cram_bam_ = select_first(Sorter.sorted_cram)
    }
    if (length(Sorter.sorted_cram) > 1) {
        Array[File] output_cram_bam_list_ = select_all(Sorter.sorted_cram)
    }

    if (length(Sorter.sorted_cram_index) == 1) {
        File output_cram_bam_index_ = select_first(Sorter.sorted_cram_index)
    }
    if (length(Sorter.sorted_cram_index) > 1) {
        Array[File] output_cram_bam_index_list_ = select_all(Sorter.sorted_cram_index)
    }


    if (length(Sorter.sorter_stats_json) == 1) {
        File sorter_stats_json_ = Sorter.sorter_stats_json[0]
    }
    if (length(Sorter.sorter_stats_json) > 1) {
        Array[File] sorter_stats_json_list_ = Sorter.sorter_stats_json
    }

    if (length(Sorter.sorter_stats_csv) == 1) {
        File sorter_stats_csv_ = Sorter.sorter_stats_csv[0]
    }
    if (length(Sorter.sorter_stats_csv) > 1) {
        Array[File] sorter_stats_csv_list_ = Sorter.sorter_stats_csv
    }


    if (output_haplotypes_cram) {
        call HSampling.HaplotypeSampling as HaplotypeSampling {
            input:
                input_cram_bam_list = input_cram_list,
                cram_reference_fasta = references.ref_fasta,
                cram_reference_fasta_index = references.ref_fasta_index,
                gbz_file = select_first([giraffe_parameters.ref_gbz_for_haplotypes]),
                hapl_file = select_first([giraffe_parameters.ref_hapl]),
                alignment_reference_fasta =   select_first([giraffe_parameters.alignment_reference_fasta_for_haplotypes]),
                alignment_reference_fasta_index = select_first([giraffe_parameters.alignment_reference_fasta_index_for_haplotypes]),
                sample_name = base_file_name
        }
    }

    output {
        # Outputs when there is only a single output cram
        File? output_cram             = output_cram_bam_
        File? output_cram_index         = output_cram_bam_index_
        File? sorter_stats_csv              = sorter_stats_csv_
        File? sorter_stats_json             = sorter_stats_json_

        # Outputs when there are multiple output crams
        Array[File?]? output_cram_bam_list       = output_cram_bam_list_
        Array[File?]? output_cram_bam_index_list = output_cram_bam_index_list_
        Array[File?]? sorter_stats_csv_list      = sorter_stats_csv_list_
        Array[File?]? sorter_stats_json_list     = sorter_stats_json_list_

        File? haplotypes_cram = HaplotypeSampling.output_cram
        File? haplotypes_cram_index = HaplotypeSampling.output_cram_index
    }
}
