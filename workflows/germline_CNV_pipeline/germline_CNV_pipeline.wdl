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
# Runs the Copy Number Variation calling pipeline for Ultima Genomics data. Includes the following main steps:
# 1. Runs cn.mops .
# 2. Runs cnvpytor.
# 3. Runs results combination of the two tools.


# CHANGELOG in reverse chronological order
# 1.26.0 - Updated annotations and filtering model, quality significantly improved
# 1.24.0 - Removed filtering on UG-CNV-LCR, updated filtering model

import "single_sample_cnmops_CNV_calling.wdl" as SingleSampleCnmopsCNVCalling
import "single_sample_CNVpytor_calling.wdl" as SingleSampleCNVpytorCalling
import "combine_germline_CNV_calls.wdl" as CombineGermlineCNVCalls
import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks

workflow GermlineCNVPipeline {

    input {
        String pipeline_version = "1.27.2" # !UnusedDeclaration

        String base_file_name
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        Array[String] ref_seq_names
        File? ug_cnv_lcr_file

        Boolean skip_filtering
        File? filtering_model
        Int? filtering_model_decision_threshold
        #cnmops params
        Int? cnmops_mapq_override
        Int? cnmops_window_length_override
        Int? cnmops_parallel_override
        Array[File] bed_graph
        File genome_windows
        File cohort_reads_count_matrix
        File ploidy_file
        Int? cnmops_min_width_value_override
        Int? cnmops_min_cnv_length_override
        Float? cnmops_intersection_cutoff_override
        Boolean? disable_mod_cnv

        #cnvpytor params
        Array[Int]? cnvpytor_window_length_override
        
        Int cushion_size
        
        Boolean? skip_figure_generation
        Boolean? no_address_override
        Int? preemptible_tries_override
        File? monitoring_script_input
        Boolean create_md5_checksum_outputs = false

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv prefix(input_bam_file_index) == input_bam_file
        #@wv suffix(input_bam_file) in {".bam", ".cram"}
        #@wv suffix(input_bam_file_index) in {".bai", ".crai"}
        #@wv reference_genome == prefix(reference_genome_index)
        #@wv suffix(reference_genome) in {'.fasta', '.fa', '.fna'}
        #@wv suffix(reference_genome_index) == '.fai'
        #@wv suffix(filtering_model) == '.pkl'
    }

    meta {
        description: "Runs: <br>1. single sample germline CNV calling workflow based on [cn.mops](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)<br>2. cnvpytor workflow<br> 3. combines results, verifies them using split reads and jump alignments<br>4. Applies ML model to estimate quality of the CNV</b>"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "SingleSampleReadsCount.monitoring_script_input",
                "no_address_override",
                "preemptible_tries_override",
                "Glob.glob",
                "SingleSampleReadsCount.Globals.glob",
                "MergeMd5sToJson.output_json", 
                'CombineCNVCalls.FilterVCF.ref_fasta',
                'CombineCNVCalls.FilterVCF.ref_fasta_idx',
                'CombineCNVCalls.FilterVCF.blacklist_file',
                'CombineCNVCalls.FilterVCF.custom_annotations',
                'CombineCNVCalls.FilterVCF.disk_size'
        ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        input_bam_file:{
            help: "Input sample BAM/CRAM file",
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
        ug_cnv_lcr_file: {
            help: "UG-CNV-LCR bed file",
            type: "File",
            category: "input_optional"
        }
        filtering_model: {
            help: "CNV filtering model, default in template, calls are not filtered if not provided",
            type: "File",
            category: "input_optional"
        }
        filtering_model_decision_threshold: {
            help: "Decision threshold for the filtering model, default is set in template. Lower- less stringent, Higher- more stringent",
            type: "Int",
            category: "param_optional"
        }
        skip_filtering: {
            help: "Whether to skip CNV filtering step, default is False",
            type: "Boolean",
            category: "param_required"
        }
        cnmops_mapq_override: {
            help : "Reads mapping-quality cutoff for coverage aggregation used in cn.mops, default value is 1",
            type: "Int",
            category: "param_advanced"
        }
        cnmops_window_length_override: {
            help: "Window length on which the read counts will be aggregated, default value is 500",
            type: "Int",
            category: "param_advanced"
        }
        window_length_override: {
            help: "Window length on which the read counts will be aggregated, default value is 500",
            type: "Int",
            category: "param_advanced"
        }
        cnmops_parallel_override: {
            help: "Number of cpus for cn.mops run. Default value is 4",
            type: "Int",
            category: "param_advanced"
        }
        bed_graph: {
            help: "Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data).",
            type: "Array[File]",
            category: "input_required"
        }
        genome_windows: {
            help: "Bed file of the genome binned to equal sized windows similar to the cohort_reads_count_matrix.",
            type: "File",
            category: "input_required"
        }
        cohort_reads_count_matrix: {
            help : "GenomicRanges object of the cohort reads count matrix in rds file format. default cohort can be found in the template.",
            type: "File",
            category: "input_required"
        }
        ploidy_file: {
            help : "X chromosome ploidy of the cohort and the additional sample. Each sample is represented on a number on a separate row. Ploidy of the default cohort can be found in the template. The last row corresponds to the sample being called",
            type: "File",
            category: "input_required"
        }
        cnmops_min_width_value_override: {
            help: "Minimum of consecutive windows with a significant signal to consider for CNV reporting. Default is: 2",
            type: "Int",
            category: "param_advanced"
        }
        cnmops_min_cnv_length_override: {
            help: "Minimum length for reporting CNV. Default is: 0",
            type: "Int",
            category: "param_advanced"
        }
        cnmops_intersection_cutoff_override: {
            help: "Intersection cutoff with UG-CNV-LCR regions to filter out CNV calls. Default is:  0.5",
            type: "Float",
            category: "param_advanced"
        }
        disable_mod_cnv:
        {
            help: "whether to call moderate cnvs (Fold-Change~1.5 will be tagged as CN2.5 and Fold-Change~0.7 will be tagged as CN1.5). Default is: True",
            type: "Boolean",
            category: "param_advanced"
        }
        cnvpytor_window_length_override: {
            help: "Window length on which the read counts will be aggregated, default value is [500,2500]",
            type: "Array[Int]",
            category: "param_advanced"
        }
        skip_figure_generation: {
            help: "Skip CNV calls figure generation. please set to True if reference genome is not hg38. Default is: False",
            type: "Boolean",
            category: "param_optional"
        }
        create_md5_checksum_outputs: {
           help: "Create md5 checksum for requested output files",
           type: "Boolean",
           category: "input_optional"
        }
        cnmops_cnv_calls_bed: {
            help: "CNMOPS CNV calls in bed format",
            type: "File",
            category: "output"
        }
        cnmops_cnv_calls_vcf: {
            help: "CNMOPS CNV calls in VCF format",
            type: "File",
            category: "output"
        }
        cnmops_cnv_calls_vcf_index: {
            help: "Index file for the CNMOPS CNV calls VCF",
            type: "File",
            category: "output"
        }
        cnvpytor_cnv_calls_bed: {
            help: "CNVpytor CNV calls in bed format",
            type: "File",
            category: "output"
        }
        cnvpytor_cnv_calls_vcf: {
            help: "CNVpytor CNV calls in VCF format",
            type: "File",
            category: "output"
        }
        cnvpytor_cnv_calls_vcf_index: {
            help: "Index file for the CNVpytor CNV calls VCF",
            type: "File",
            category: "output"  
        }
        cushion_size: {
            help: "Cushion size around CNV breakpoints for split-read analysis and jump alignment analysis",
            type: "Int",
            category: "param_required"
        }
        combined_cnv_calls_bed: {
            help: "Final (combined) CNV calls in bed format",
            type: "File",
            category: "output"
        }

        combined_cnv_calls_bed_vcf: {
            help: "Combined CNV calls in vcf format",
            type: "File",
            category: "output"
        }
        combined_cnv_calls_bed_vcf_index: {
            help: "Index of the combined CNV calls in vcf format",
            type: "File",
            category: "output"
        }
        combine_read_evidence: {
            help: "BAM file with read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        combine_read_evidence_index: {
            help: "Index file for the BAM with read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        combine_read_scores_csv: {
            help: "CSV file with jalign scores for each read",
            type: "File",
            category: "output"
        }
        md5_checksums_json: {
            help: "json file that will contain md5 checksums for requested output files",
            type: "File",
            category: "output"
        }
    }
    Int cnmops_mapq = select_first([cnmops_mapq_override, 1])
    Int cnmops_window_length = select_first([cnmops_window_length_override, 500])
    Int cnmops_parallel = select_first([cnmops_parallel_override, 4])
    Int cnmops_min_width_value = select_first([cnmops_min_width_value_override, 2])
    Int cnmops_min_cnv_length = select_first([cnmops_min_cnv_length_override, 0])
    Float cnmops_intersection_cutoff = select_first([cnmops_intersection_cutoff_override, 0.5])
    Boolean enable_mod_cnv = select_first([disable_mod_cnv, true])
    Array[Int] cnvpytor_window_lengths = select_first([cnvpytor_window_length_override, [500,2500]])
    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    Boolean skip_figure_generation_value = select_first([skip_figure_generation, false])
    
    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script]) #!FileCoercion

    call SingleSampleCnmopsCNVCalling.SingleSampleCnmopsCNVCalling as CnmopsCNVCalling{
        input:
            base_file_name = base_file_name,
            reference_genome = reference_genome,
            reference_genome_index = reference_genome_index,
            mapq = cnmops_mapq,
            ref_seq_names = ref_seq_names,
            window_length = cnmops_window_length,
            parallel = cnmops_parallel,
            bed_graph = bed_graph,
            genome_windows = genome_windows,
            cohort_reads_count_matrix = cohort_reads_count_matrix,
            ploidy_file = ploidy_file,
            min_width_value = cnmops_min_width_value,
            min_cnv_length = cnmops_min_cnv_length,
            intersection_cutoff = cnmops_intersection_cutoff,
            cnv_lcr_file = ug_cnv_lcr_file,
            enable_mod_cnv_override = enable_mod_cnv,
            skip_figure_generation = skip_figure_generation_value,
            preemptible_tries_override = preemptible_tries,
            no_address_override = no_address,
            monitoring_script_input = monitoring_script
    }

    call SingleSampleCNVpytorCalling.SingleSampleCNVpytorCalling as CnvpytorCNVCalling{
        input:
        base_file_name = base_file_name,
        input_bam_file = input_bam_file,
        input_bam_file_index = input_bam_file_index,
        reference_genome = reference_genome,
        reference_genome_index = reference_genome_index,
        ref_seq_names = ref_seq_names,
        window_lengths = cnvpytor_window_lengths,
        preemptible_tries_override = preemptible_tries,
        no_address_override = no_address,
        monitoring_script_input = monitoring_script

    }
    
    call CombineGermlineCNVCalls.CombineGermlineCNVCalls as CombineCNVCalls{
        input:
            base_file_name = base_file_name,
            cnmops_cnvs_vcf = CnmopsCNVCalling.out_sample_cnvs_vcf,
            cnmops_cnvs_vcf_index = CnmopsCNVCalling.out_sample_cnvs_vcf_index,
            cnvpytor_cnvs_vcf = CnvpytorCNVCalling.cnvpytor_cnv_calls_vcf,
            cnvpytor_cnvs_vcf_index = CnvpytorCNVCalling.cnvpytor_cnv_calls_vcf_index,            
            input_bam_file = input_bam_file,
            input_bam_file_index = input_bam_file_index,
            cushion_size = cushion_size,
            reference_genome = reference_genome,
            reference_genome_index = reference_genome_index,
            cnv_lcr_file = ug_cnv_lcr_file,
            skip_filtering = skip_filtering,
            filtering_model = filtering_model,
            filtering_model_decision_threshold = filtering_model_decision_threshold,
            monitoring_script_input = monitoring_script,
            preemptible_tries_override = preemptible_tries,
            no_address_override = no_address
    }

    File combined_cnv_calls_bed_vcf_ = CombineCNVCalls.out_sample_cnvs_vcf
    File combined_cnv_calls_bed_vcf_index_ = CombineCNVCalls.out_sample_cnvs_vcf_index
    if (create_md5_checksum_outputs) {
        Array[File] output_files = [combined_cnv_calls_bed_vcf_, combined_cnv_calls_bed_vcf_index_]

        scatter (file in output_files) {
            call UGGeneralTasks.ComputeMd5 as compute_md5 {
                input:
                    input_file = file,
                    docker = global.ubuntu_docker,
            }
        }

        call UGGeneralTasks.MergeMd5sToJson {
            input:
                md5_files = compute_md5.checksum,
                docker = global.ugbio_core_docker
        }
    }
    output {
        File cnmops_cnv_calls_bed = CnmopsCNVCalling.out_sample_cnvs_bed
        File cnmops_cnv_calls_vcf = CnmopsCNVCalling.out_sample_cnvs_vcf
        File cnmops_cnv_calls_vcf_index = CnmopsCNVCalling.out_sample_cnvs_vcf_index
        File cnvpytor_cnv_calls_bed = CnvpytorCNVCalling.cnvpytor_cnv_calls_tsv
        File cnvpytor_cnv_calls_vcf = CnvpytorCNVCalling.cnvpytor_cnv_calls_vcf
        File cnvpytor_cnv_calls_vcf_index = CnvpytorCNVCalling.cnvpytor_cnv_calls_vcf_index
        File combined_cnv_calls_bed = CombineCNVCalls.out_sample_cnvs_bed
        File combined_cnv_calls_bed_vcf = combined_cnv_calls_bed_vcf_
        File combined_cnv_calls_bed_vcf_index = combined_cnv_calls_bed_vcf_index_
        File combine_read_evidence = CombineCNVCalls.read_evidence
        File combine_read_evidence_index = CombineCNVCalls.read_evidence_index
        File combine_read_scores_csv = CombineCNVCalls.read_scores_csv
        File? md5_checksums_json = MergeMd5sToJson.md5_json
    }
}

