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
# Runs the Copy Number Variation calling pipeline for Ultima Genomics data. Includes the following main steps:
# 1. Reads count calculation .
# 2. Merge sample reads count data to a given cohort reads count matrix.
# 3. CNV calling using cn.mops algorithm
# 4. Filtering sample's CNV calls

# CHANGELOG in reverse chronological order
# 1.11.0 Added an option to call CNVs that are less than duplications (mosaic?)
# 1.9.0 Support bedGraph input format as external pre-calculated sample's coverage.
# 1.6.0 Normalization now genome-wide rather than per chromosome (allows correct calling of chrY, chrX)
# 1.5.0 Faster

import "single_sample_cnmops_reads_count.wdl" as ReadsCount
import "tasks/cnv_calling_tasks.wdl" as CnvTasks
import "tasks/globals.wdl" as Globals

workflow SingleSampleCnmopsCNVCalling {

    input {
        String pipeline_version = "1.11" # !UnusedDeclaration

        String base_file_name

        File? input_bam_file
        File? input_bam_file_index
        File reference_genome      #ref-genome+idx to enable cram as input file
        File reference_genome_index
        Int mapq
        Array[String] ref_seq_names
        Int window_length
        Int parallel

        #extenal reads count
        File? input_sample_reads_count
        File? bed_graph
        File? genome_windows

        File cohort_reads_count_matrix
        File? merged_cohort_ploidy_file
        String? chrX_name
        String? chrY_name
        Boolean? cap_coverage_override
        Int min_width_value = 2

        Int min_cnv_length = 10000
        Float intersection_cutoff = 0.5
        File cnv_lcr_file
        Boolean? enable_moderate_amplifications_override

        Boolean? save_hdf_override
        Boolean? save_csv_override
        Boolean? no_address_override
        Int? preemptible_tries_override

        # Used for running on other clouds (aws)
        File? monitoring_script_input

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv defined(input_bam_file) -> (prefix(input_bam_file_index) == input_bam_file)
        #@wv defined(input_bam_file) -> (suffix(input_bam_file) in {".bam", ".cram"})

        #@wv reference_genome == prefix(reference_genome_index)
        #@wv suffix(reference_genome) in {'.fasta', '.fa', '.fna'}
        #@wv suffix(reference_genome_index) == '.fai'

        #@wv defined(bed_graph) -> defined(genome_windows)
        #@wv defined(bed_graph) -> not(defined(input_sample_reads_count))
        #@wv defined(input_sample_reads_count) -> not(defined(bed_graph))

        #@wv intersection_cutoff <= 1 and intersection_cutoff >= 0
        #@wv min_width_value > 0
        #@wv len(ref_seq_names) >= 0
        #@wv min_cnv_length >= window_length

    }

    meta {
        description: "Runs single sample germline CNV calling workflow based on \<a href=\"https://bioconductor.org/packages/release/bioc/html/cn.mops.html\"\>cn.mops</a>\n\nThe pipeline uses a given cohort's coverage profile for normalization.\n\nThe pipeline can recieve one of the following options as input:\n\n&nbsp;&nbsp;1. Input CRAM/BAM file. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_template.json\n\n&nbsp;&nbsp;2. A rds file which stores a GenomicRanges object with coverage collected in the same windows as the given cohort. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_skip_reads_count_template.json\n\n&nbsp;&nbsp;3. A BedGraph holding the coverage per location. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_input_bedGraph_template.json\n\nThe pipeline calls CNVs for the given sample and filters them by length (>10,000b) and overlap with UG-CNV-LCR."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "SingleSampleReadsCount.monitoring_script_input",
                "no_address_override",
                "Globals.glob",
                "SingleSampleReadsCount.Globals.glob"
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        input_bam_file: {
            help: "Input sample BAM/CRAM file. one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set",
            type: "File",
            category: "input_optional"
         }
        input_bam_file_index: {
            help:"Input sample BAI/CRAI index file",
            type: "File",
            category: "input_optional"
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
        mapq: {
            help : "Reads mapping-quality cutoff for coverage aggregation, recommended value set in the template",
            type: "Int",
            category: "param_required"
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
        parallel: {
            help: "Number of cpus for cn.mops run",
            type: "Int",
            category: "param_advanced"
        }
        input_sample_reads_count: {
            help: "Inputs sample windowed coverage stored as GenomicRanges object in rds file. can be calculated using cn.mops::getReadCountsFromBAM R function.  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set",
            type: "File",
            category: "input_optional"
        }
        bed_graph: {
            help: "Previously calculated input bedGraph holding the coverage per base (outputs with the sequencing data).  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set",
            type: "File",
            category: "input_optional"
        }
        genome_windows: {
            help: "Bed file of the genome binned to equal sized windows similar to the cohort_reads_count_matrix. if bed_graph input is set, this file must be given. ",
            type: "File",
            category: "input_optional"
        }
        cohort_reads_count_matrix: {
            help : "GenomicRanges object of the cohort reads count matrix in rds file format. default cohort can be found in the template. can be created by cn.mops::getReadCountsFromBAM R function ",
            type: "File",
            category: "input_required"
         }
        merged_cohort_ploidy_file: {
            help : "Cohort ploidy file indicating 1 for male and 2 for female, per sample. The number of lines should be the same as the number of samples in cohort + current_sample. if not given, defaults to 2 for all samples.",
            type: "File",
            category: "input_optional"
        }
        chrX_name: {
            help: "The name of the female sex chromosome in the genome. default is: chrX",
            type: "String",
            category: "param_optional"
        }
        chrY_name: {
           help: "The name of the male sex chromosome in the genome. default is: chrY",
           type: "String",
           category: "param_optional"
       }
        min_width_value: {
            help: "Minimum of consecutive windows with a significant signal to consider for CNV reporting. Default is: 2",
            type: "Int",
            category: "param_advanced"
        }
        min_cnv_length: {
            help: "Minimum length for reporting CNV. Default is: 10,000",
            type: "Int",
            category: "param_required"
        }
        intersection_cutoff: {
            help: "Intersection cutoff with UG-CNV-LCR regions to filter out CNV calls. Default is:  0.5",
            type: "Float",
            category: "param_required"
        }
        cnv_lcr_file: {
            help: "UG-CNV-LCR bed file",
            type: "File",
            category: "param_required"
        }
        save_hdf_override: {
            help: "Whether to save sample reads counts/cohort including sample/cnmops output data in hdf5 format (additionally to RDS format). Default is: False.",
            type: "Boolean",
            category: "param_optional"
        }
        save_csv_override: {
            help: "Whether to save sample reads counts/cohort including sample/cnmops output data in csv format (additionally to RDS format). Default is: False.",
            type: "Boolean",
            category: "param_optional"
        }
        preemptible_tries_override: {
            help: "Number of preemptible tries,default is: 1",
            type: "Int",
            category: "param_optional"
        }
        out_sample_cnvs_bed:{
            help: "Bed file with sample's called CNVs",
            type: "File",
            category: "output"
        }
        out_sample_cnvs_filtered_bed:{
            help: "Bed file with CNVs filtered by length and overlap with low confidence regions",
            type: "File",
            category: "output"
        }
        out_sample_reads_count : {
            help: "GenomicRanges object of the sample's reads count in rds file format",
            type: "File",
            category: "output"
        }
        out_sample_reads_count_hdf5:{
            help: "GenomicRanges object of the sample's reads count in hdf5 file format",
            type: "File",
            category: "output_optional"
        }
        enable_moderate_amplifications_override:
        {
            help: "whether to call moderate amplifications (Fold-Change>1.5 & < 2 will be tagged as CN2.5) Default is: False",
            type: "Boolean",
            category: "param_optional"
        }
        cap_coverage_override:
        {
            help: "whether to cap extremely high average coverage windows to 2*cohort's average coverage quantile 99.9% value",
            type: "Boolean",
            category: "param_optional"
        }
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    Boolean skip_reads_count = (defined(input_sample_reads_count) || defined(bed_graph))
    Boolean run_convert_bedGraph_to_Granges = defined(bed_graph)
    Boolean save_hdf = select_first([save_hdf_override , false])
    Boolean save_csv = select_first([save_csv_override , false])
    Boolean enable_moderate_amplifications = select_first([enable_moderate_amplifications_override, false])
    Boolean cap_coverage = select_first([cap_coverage_override, false])

    call Globals.Globals as Globals
      GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    if(skip_reads_count == false) {
        File input_bam = select_first([input_bam_file])
        File input_bai = select_first([input_bam_file_index])
        call ReadsCount.SingleSampleCnmopsReadsCount as SingleSampleReadsCount {
                input:
                input_bam_file = input_bam,
                input_bam_file_index = input_bai,
                reference_genome = reference_genome,      #ref-genome+idx to enable cram as input file
                reference_genome_index = reference_genome_index,
                mapq = mapq,
                ref_seq_names = ref_seq_names,
                window_length = window_length,
                base_file_name = base_file_name,
                save_hdf_override = save_hdf,
                no_address_override = no_address_override,
                preemptible_tries_override = preemptible_tries_override
        }
    }
    if(run_convert_bedGraph_to_Granges)
    {
        File input_bed_graph = select_first([bed_graph])
        File input_genome_windows = select_first([genome_windows])
        call CnvTasks.ConvertBedGraphToGranges as ConvertBedGraphToGranges{
        input:
            sample_name = base_file_name,
            input_bed_graph = input_bed_graph,
            genome_windows = input_genome_windows,
            genome_file = reference_genome_index,
            docker = global.ug_vc_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
        }
    }

    File sample_reads_count_file = select_first([input_sample_reads_count,ConvertBedGraphToGranges.out_RC_Granges, SingleSampleReadsCount.out_reads_count])

    call CnvTasks.AddCountsToCohortMatrix {
        input:
        sample_reads_count = sample_reads_count_file,
        cohort_reads_count_matrix = cohort_reads_count_matrix,
        docker = global.ug_vc_docker,
        save_hdf = save_hdf,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }

    call CnvTasks.RunCnmops {
        input:
        merged_cohort_reads_count_matrix = AddCountsToCohortMatrix.merged_cohort_reads_count_matrix,
        min_width_value = min_width_value,
        ploidy = merged_cohort_ploidy_file,
        chrX_name = chrX_name,
        chrY_name = chrY_name,
        cap_coverage = cap_coverage,
        docker = global.ug_vc_docker,
        save_hdf = save_hdf,
        save_csv = save_csv,
        moderate_amplificiations = enable_moderate_amplifications,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries,
        parallel = parallel
    }

    Array[String] sample_names = [base_file_name]
    call CnvTasks.FilterSampleCnvs {
        input:
        cohort_cnvs_csv = RunCnmops.cohort_cnvs_csv,
        sample_names = sample_names,
        min_cnv_length = min_cnv_length,
        intersection_cutoff = intersection_cutoff,
        cnv_lcr_file = cnv_lcr_file,
        docker = global.ug_vc_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }
    output {
        File out_sample_reads_count = sample_reads_count_file
        File? out_sample_reads_count_hdf5 = SingleSampleReadsCount.out_reads_count_hdf5
        Array[File] out_sample_cnvs_bed = FilterSampleCnvs.sample_cnvs_bed
        Array[File] out_sample_cnvs_filtered_bed = FilterSampleCnvs.sample_cnvs_filtered_bed
    }
}

