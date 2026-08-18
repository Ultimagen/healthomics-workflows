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
# 1.35.0 Removed UG-CNV-LCR  annotation - no longer needed
#        The intervals and the chromosome list are taken from the cohort reads count matrix, ref_seq_names was removed
#        The reference is selected with reference_genome (hg38, b37, hg38_nist_v3_with_decoy) instead of a fasta+fai pair
#        The cohort reads count matrix comes from the genome resources and can be overridden with cohort_reads_count_matrix_override
#        The sex chromosome names are derived from reference_genome so the ploidy file takes effect on b37
# 1.11.0 Added an option to call CNVs that are less than duplications (mosaic?)
# 1.9.0 Support bedGraph input format as external pre-calculated sample's coverage.
# 1.6.0 Normalization now genome-wide rather than per chromosome (allows correct calling of chrY, chrX)
# 1.5.0 Faster

import "single_sample_cnmops_reads_count.wdl" as ReadsCount
import "tasks/cnv_calling_tasks.wdl" as CnvTasks
import "tasks/globals.wdl" as Globals
import "tasks/genome_resources.wdl" as GenomeResourcesLib

workflow SingleSampleCnmopsCNVCalling {

    input {
        String pipeline_version = "1.35.0" # !UnusedDeclaration

        String base_file_name
        String? sample_name
        File? input_bam_file
        File? input_bam_file_index
        String reference_genome = "hg38"
        Int mapq
        Int window_length
        Int parallel

        #extenal reads count
        File? input_sample_reads_count
        Array[File]? bed_graph

        File? cohort_reads_count_matrix_override
        File? ploidy_file
        String? chrX_name_override
        String? chrY_name_override
        Boolean? cap_coverage_override
        Int min_width_value = 2

        Int min_cnv_length = 10000
        Boolean? enable_mod_cnv_override

        Boolean? skip_figure_generation
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

        #@wv reference_genome in {"hg38", "b37", "hg38_nist_v3_with_decoy"}

        #@wv defined(bed_graph) -> not(defined(input_sample_reads_count))
        #@wv defined(input_sample_reads_count) -> not(defined(bed_graph))

        #@wv min_width_value > 0

    }

    meta {
        description: "Runs single sample germline CNV calling workflow based on [cn.mops](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)\n\nThe pipeline uses a given cohort's coverage profile for normalization.\n\nThe pipeline can receive one of the following options as input:\n\n&nbsp;&nbsp;1. Input CRAM/BAM file. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_template.json\n\n&nbsp;&nbsp;2. A rds file which stores a GenomicRanges object with coverage collected in the same windows as the given cohort. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_skip_reads_count_template.json\n\n&nbsp;&nbsp;3. A BedGraph holding the coverage per location. Corresponding template: Input_templates/single_sample_cnmops_CNV_calling_input_bedGraph_template.json\n\nThe pipeline calls CNVs for the given sample and filters them by length (>10,000b).\n\n<b>When Running in AWS HealthOmics this pipeline should run with [dynamic storage](https://docs.omics.ai/products/workbench/engines/parameters/aws-healthomics#storage_type-dynamic-or-static)</b>"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "SingleSampleReadsCount.monitoring_script_input",
                "no_address_override",
                "Glob.glob",
                'ProcessCnmopsCnvs.intersection_cutoff',
                'ProcessCnmopsCnvs.cnv_lcr_file'
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Base name for the output file, if sample_name not provided - will also be the sample name in the VCF",
            type: "String",
            category: "input_required"
        }
        sample_name: {
            help: "Sample name for the output VCF. if not provided, base_file_name will be used as sample name in the VCF",
            type: "String",
            category: "input_optional"
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
            help: "Genome type selector. Supported values are hg38, b37 and hg38_nist_v3_with_decoy. The reference files and the cn.mops cohort reads count matrix are taken from the genome resources accordingly",
            type: "String",
            category: "ref_required"
         }
        mapq: {
            help : "Reads mapping-quality cutoff for coverage aggregation, recommended value set in the template",
            type: "Int",
            category: "param_required"
       }
        window_length: {
            help: "Window length on which the read counts will be aggregated. The cohort reads count matrix is rebinned to this window length (must be a multiple of the cohort's window length)",
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
            help: "Previously calculated input bedGraph files holding the coverage per base (outputs with the sequencing data).  one of the `input_bam_file`, `input_sample_reads_count` or `bed_graph` must be set",
            type: "File",
            category: "input_optional"
        }
        cohort_reads_count_matrix_override: {
            help : "GenomicRanges object of the cohort reads count matrix in rds file format. By default the cohort matching reference_genome is taken from the genome resources. Can be created by cn.mops::getReadCountsFromBAM R function ",
            type: "File",
            category: "input_optional"
         }
        # TODO(BIOIN-3020): bundle the ploidy file with the cohort in genome_resources so a future
        # genome-specific ploidy file is selected automatically. Kept as a plain File for now because
        # all three v3.0 cohorts share a byte-identical ploidy file (same 60-sample batch order).
        ploidy_file: {
            help : "X chromosome ploidy in the cohort. 1 for male and 2 for female, per sample. The number of lines should be the same as the number of samples in cohort + current_sample. if not given, defaults to 2 for all samples. Genome independent.",
            type: "File",
            category: "input_optional"
        }
        chrX_name_override: {
            help: "The name of the female sex chromosome in the cohort reads count matrix. By default derived from reference_genome ('X' for b37, 'chrX' otherwise)",
            type: "String",
            category: "param_optional"
        }
        chrY_name_override: {
           help: "The name of the male sex chromosome in the cohort reads count matrix. By default derived from reference_genome ('Y' for b37, 'chrY' otherwise)",
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
        skip_figure_generation: {
            help: "Whether to skip figure generation. Default is: False",
            type: "Boolean",
            category: "param_optional"
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
        out_sample_cnvs_vcf:{
            help: "VCF file with sample's called CNVs",
            type: "File",
            category: "output"
        }
        out_sample_cnvs_vcf_index:{
            help: "Index file for the VCF file with sample's called CNVs",
            type: "File",
            category: "output"
        }
        out_sample_cnvs_filtered_bed:{
            help: "Bed file with CNVs filtered by length",
            type: "File",
            category: "output"
        }
        out_sample_norm_coverage_bed: {
            help: "Normalized coverage bed file for the sample",
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
        out_sample_merged_bedGraph:{
            help: "Merged bedGraph file of the sample's coverage",
            type: "File",
            category: "output_optional"
        }
        out_coverage_plot_files:{
            help: "List of coverage figures (file per sample)",
            type: "Array[File]",
            category: "output"
        }
        out_dup_del_plot_files:{
            help: "List of duplication/deletion figures (file per sample)",
            type: "Array[File]",
            category: "output"
        }
        out_copy_number_plot_files:{
            help: "List of copy number figures (file per sample)",
            type: "Array[File]",
            category: "output"
        }
        sample_norm_read_counts_bed:{
            help: "Bed file with normalized read counts of the sample",
            type: "File",
            category: "output"
        }

        enable_mod_cnv_override:
        {
            help: "whether to call moderate cnvs (Fold-Change~1.5 will be tagged as CN2.5 and Fold-Change~0.7 will be tagged as CN1.5). Default is: False",
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
    Boolean enable_mod_cnv = select_first([enable_mod_cnv_override, false])
    Boolean cap_coverage = select_first([cap_coverage_override, false])
    Boolean skip_figure_generation_value = select_first([skip_figure_generation, false])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])    #!FileCoercion

    # Get genome resources based on reference_genome
    call GenomeResourcesLib.GenomeResourcesWorkflow as GenomeResources

    File ref_fasta = GenomeResources.resources[reference_genome].ref_fasta
    File ref_fasta_index = GenomeResources.resources[reference_genome].ref_fasta_index

    File cohort_reads_count_matrix = select_first([cohort_reads_count_matrix_override,
        GenomeResources.resources[reference_genome].cnv_normalization_cohort])

    # The sex chromosome names must match the contig naming of the cohort reads count matrix,
    # otherwise the ploidy correction in cn.mops is silently skipped.
    # TODO(BIOIN-3021): move these into genome_resources once it supports String/Array[String]
    # values (today it only holds cloud-path File values, so the derivation stays inline here).
    String chrX_name = select_first([chrX_name_override, if reference_genome == "b37" then "X" else "chrX"])
    String chrY_name = select_first([chrY_name_override, if reference_genome == "b37" then "Y" else "chrY"])

    # Always rebin cohort to match workflow window_length
    # R script handles no-op case when cohort already at correct window size
    call CnvTasks.RebinCohortReadsCount {
        input:
            cohort_reads_count_matrix = cohort_reads_count_matrix,
            new_window_length = window_length,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    # Extract the genome windows and the chromosome names from the rebinned cohort matrix.
    # The cohort defines the intervals in which CNVs are called, so no interval input is needed.
    call CnvTasks.ExtractGenomeWindows as ExtractGenomeWindows {
        input:
            cohort_reads_count_matrix = RebinCohortReadsCount.rebinned_cohort_reads_count_matrix,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    if(skip_reads_count == false) {
        File input_bam = select_first([input_bam_file])
        File input_bai = select_first([input_bam_file_index])
        call ReadsCount.SingleSampleCnmopsReadsCount as SingleSampleReadsCount {
                input:
                input_bam_file = input_bam,
                input_bam_file_index = input_bai,
                reference_genome = ref_fasta,      #ref-genome+idx to enable cram as input file
                reference_genome_index = ref_fasta_index,
                mapq = mapq,
                genome_windows = ExtractGenomeWindows.genome_windows,
                base_file_name = base_file_name,
                save_hdf_override = save_hdf,
                no_address_override = no_address_override,
                preemptible_tries_override = preemptible_tries_override
        }
    }

    String sample_name_defined = select_first([sample_name, SingleSampleReadsCount.out_sample_name,base_file_name])

    if(run_convert_bedGraph_to_Granges)
    {
        Array[File] input_bed_graph = select_first([bed_graph])

        call CnvTasks.ConvertBedGraphToGranges as ConvertBedGraphToGranges{
        input:
            sample_name = sample_name_defined,
            input_bed_graph = input_bed_graph,
            genome_windows = ExtractGenomeWindows.genome_windows,
            genome_file = ref_fasta_index,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
        }
    }

    File sample_reads_count_file = select_first([input_sample_reads_count,ConvertBedGraphToGranges.out_RC_Granges, SingleSampleReadsCount.out_reads_count])

    call CnvTasks.AddCountsToCohortMatrix {
        input:
        sample_reads_count = sample_reads_count_file,
        cohort_reads_count_matrix = RebinCohortReadsCount.rebinned_cohort_reads_count_matrix,
        docker = global.ugbio_cnv_docker,
        save_hdf = save_hdf,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }

    call CnvTasks.RunCnmops {
        input:
            merged_cohort_reads_count_matrix = AddCountsToCohortMatrix.merged_cohort_reads_count_matrix,
            min_width_value = min_width_value,
            ploidy = ploidy_file,
            chrX_name = chrX_name,
            chrY_name = chrY_name,
            cap_coverage = cap_coverage,
            docker = global.ugbio_cnv_docker,
            save_hdf = save_hdf,
            save_csv = save_csv,
            mod_cnv = enable_mod_cnv,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries,
            parallel = parallel
    }

    Array[String] sample_names = [sample_name_defined]

    call CnvTasks.ExtractNormalizedReadCount{
        input: 
            cohort_reads_count_norm = RunCnmops.cohort_reads_count_norm,
            sample_names = sample_names,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    call CnvTasks.ProcessCnmopsCnvs {
        input:
            cohort_cnvs_csv = RunCnmops.cohort_cnvs_csv,
            sample_names = sample_names,
            min_cnv_length = min_cnv_length,
            sample_norm_coverage_file = ExtractNormalizedReadCount.sample_reads_count_bed,
            cohort_norm_avg_coverage_file = ExtractNormalizedReadCount.cohort_reads_count_bed,
            skip_figure_generation = skip_figure_generation_value,
            ref_genome_file = ref_fasta_index,
            germline_coverage_rds = sample_reads_count_file,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }



    scatter (vcf_file in zip(sample_names, ProcessCnmopsCnvs.sample_cnvs_vcf)) {
        call CnvTasks.AddCIPOS { 
            input:
                input_vcf = vcf_file.right,
                base_file_name = base_file_name,
                window_size = window_length,
                docker = global.ugbio_cnv_docker,
                monitoring_script = monitoring_script,
                no_address = no_address,
                preemptible_tries = preemptible_tries
        }

        call CnvTasks.CnvVcfToBed {
            input:
                input_cnv_vcf = AddCIPOS.output_vcf,
                base_file_name = vcf_file.left,
                docker = global.ugbio_cnv_docker,
                monitoring_script = monitoring_script,
                no_address = no_address,
                preemptible_tries = preemptible_tries
        }
    }

    output {
        File out_sample_reads_count = sample_reads_count_file
        File? out_sample_merged_bedGraph = ConvertBedGraphToGranges.merged_bedGraph
        File? out_sample_reads_count_hdf5 = SingleSampleReadsCount.out_reads_count_hdf5
        File out_sample_cnvs_vcf = AddCIPOS.output_vcf[0]
        File out_sample_cnvs_vcf_index = AddCIPOS.output_vcf_index[0]
        File out_sample_cnvs_bed = CnvVcfToBed.output_cnv_bed[0]
        File out_sample_norm_coverage_bed = ExtractNormalizedReadCount.sample_reads_count_bed[0]
        Array[File] out_coverage_plot_files = ProcessCnmopsCnvs.coverage_plot
        Array[File] out_dup_del_plot_files = ProcessCnmopsCnvs.dup_del_plot
        Array[File] out_copy_number_plot_files = ProcessCnmopsCnvs.copy_number_plot
    }
}

