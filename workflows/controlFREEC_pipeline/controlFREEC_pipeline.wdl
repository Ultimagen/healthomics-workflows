version 1.0
# LICENSE
#   Copyright 2023 Ultima Genomics
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
# Runs Somatic Copy Number Variation calling pipeline for Tumor-Normal pair using Ultima Genomics data.
# Includes the following main steps:
# 1. Reads count calculation for tumor and normal samples seperatley.
# 2. mPileup generation for tumor and normal samples seperatley.
# 3. CNV calling using controlFREEC algorithm
# 4. Filter CNV calls

# CHANGELOG in reverse chronological order
# 1.11.0 Better parallelization in coverage collection
# 1.9.2 Aded GEM mappability file
import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/qc_tasks.wdl" as UGQCTasks

workflow SomaticCNVCallingControlFREEC{
    input{
        String pipeline_version = "1.15.2" # !UnusedDeclaration
        String base_file_name

        # input bam files need to be supplied even if coverage and pileup are supplied externally.
        Array[File]+ input_tumor_cram_bam_file
        Array[File]+ input_tumor_cram_bam_file_index
        Array[File]+ input_normal_cram_bam_file
        Array[File]+ input_normal_cram_bam_file_index

        References references

        # Scatter interval list args
        File interval_list
        Int num_shards
        Int scatter_intervals_break
        String dummy_input_for_call_caching = ""

        #Create mpileup args
        File snp_file
        File snp_file_index

        #external mpileup
        File? normal_mpileup_override
        File? tumor_mpileup_override

        #coverage collection args
        Int? mapq_override
        File genome_windows
        Array[String]? collect_coverage_region

        #external coverage cpn format
        File? normal_coverage_cpn
        File? tumor_coverage_cpn
        #external coverage - sorter format
        Array[File]? normal_sorter_zipped_bed_graph
        Array[File]? tumor_sorter_zipped_bed_graph

        #controlFREEC args
        Boolean high_sensitivity_mode
        Boolean naive_normalization
        Boolean bed_graph_output
        Boolean contamination_adjustment
        Float? contamination_fraction
        String input_format
        Int window
        File? chrLenFile_override
        Int degree
        Int? ploidy
        Int? max_threads_override
        String? sex
        File? gem_mappability_file
        Boolean? no_address_override
        Int? preemptible_tries_override

        #Filter CNV calls
        Float? CNV_gain_cutoff_override
        Float? CNV_loss_cutoff_override
        Int min_cnv_length
        Float intersection_cutoff
        File cnv_lcr_file

        # Used for running on other clouds (aws)
        File? monitoring_script_input
        String? cloud_provider_override

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)

        #@wv len(input_tumor_cram_bam_file) > 0
        #@wv len(input_tumor_cram_bam_file) == len(input_tumor_cram_bam_file_index)
        #@wv suffix(input_tumor_cram_bam_file) <= {".bam", ".cram"}
        #@wv suffix(input_tumor_cram_bam_file_index) <= {".bai", ".crai"}
        #@wv prefix(input_tumor_cram_bam_file_index) == input_tumor_cram_bam_file
        
        #@wv len(input_normal_cram_bam_file) > 0
        #@wv len(input_normal_cram_bam_file) == len(input_normal_cram_bam_file_index)
        #@wv suffix(input_normal_cram_bam_file) <= {".bam", ".cram"}
        #@wv suffix(input_normal_cram_bam_file_index) <= {".bai", ".crai"}
        #@wv prefix(input_normal_cram_bam_file_index) == input_normal_cram_bam_file
        
        #@wv references['ref_fasta'] == prefix(references['ref_fasta_index'])
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa', '.fna'}
        #@wv suffix(references['ref_fasta_index']) == '.fai'

        #@wv defined(normal_sorter_zipped_bed_graph) <-> defined(tumor_sorter_zipped_bed_graph)
        #@wv defined(normal_coverage_cpn) <-> defined(tumor_coverage_cpn)
        #@wv defined(normal_sorter_zipped_bed_graph) -> not(defined(normal_coverage_cpn))
        #@wv defined(normal_coverage_cpn) -> not(defined(normal_sorter_zipped_bed_graph))

        #@wv defined(normal_mpileup_override) <-> defined(tumor_mpileup_override)

        #@wv mapq_override >= 0
        #@wv window > 0
        #@wv degree > 0
        #@wv defined(ploidy) -> ploidy > 0
        #@wv defined(sex) -> sex in {'XX','XY'}

        #@wv min_cnv_length >  0
        #@wv intersection_cutoff <1 and intersection_cutoff >0

        #@wv defined(CNV_gain_cutoff_override) -> CNV_gain_cutoff_override > 1
        #@wv defined(CNV_loss_cutoff_override) -> CNV_loss_cutoff_override < 1
    
    }
    Boolean run_createMpileup = !(defined(normal_mpileup_override))
    Boolean run_bedgraph_to_cpn = !(defined(normal_coverage_cpn))
    Boolean run_collect_coverage = length(select_all([normal_coverage_cpn,normal_sorter_zipped_bed_graph])) == 0

    Int mapq = select_first([mapq_override,1])
    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    File chrLenFile = select_first([chrLenFile_override,references.ref_fasta_index])
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    String input_tumor_cram_bam_file_name = input_tumor_cram_bam_file[0]
    String input_normal_cram_bam_file_name = input_normal_cram_bam_file[0]
    Int maxThreads = select_first([max_threads_override,8])
    Float CNV_gain_cutoff = select_first([CNV_gain_cutoff_override,1.03])
    Float CNV_loss_cutoff = select_first([CNV_gain_cutoff_override,0.97])
    
    meta {
        description: "Runs single sample somatic CNV calling workflow based on \<a href=\"https://boevalab.inf.ethz.ch/FREEC/\"\>ControlFREEC</a>.\n\nCNVs are called based on both coverage and allele frequencies in the tumor and the matched germline sample.\n\nThe pipeline will gather coverage and allele frequencies, run controlFREEC and filter called CNVs by length and low confidence regions.\n\ncoverage will be calculted based on the input cram/bam. Alternativley, it can recieve coverage input as one of:bedGraph, cpn formats.\n\nAllele frequencies will be calculated based on the input cram/bam and a given vcf file to specify locations. Alternativley, it can recieve precalculated frequencies as mpileup format.\n\nPipeline has an option to run in High-Sensitivity-Mode which can be used for low tumor purity samples. in this case segmentation results will be outputed and filtered by their average fold change.\n\nfor High-Sensitivity-Mode fold changes for gain and loss calls can be defined by the user. (default cutoff values are [gain,loss]=[1.03,0.97])\n\nThe pipeline outputs: \n\n&nbsp;&nbsp;- calculated coverage for tumor and normal samples\n\n&nbsp;&nbsp;- calculated mpileup for tumor and normal samples\n\n&nbsp;&nbsp;- called CNVs + filtered called CNVs\n\n&nbsp;&nbsp;- controlFREEC run-summary\n\n&nbsp;&nbsp;-coverage plot that shows normalized (log scale) coverage along the genome for the germline and tumor samples.\n\n&nbsp;&nbsp;-duplications and deletions figure - showing gains and losses along the genome.\n\n&nbsp;&nbsp;-copy-number figure  shows the copy number along the genome.\n\n"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "cloud_provider_override",
                "no_address_override",
                "dummy_input_for_call_caching",
                "tumor_createMpileup_scatter.max_depth",
                "tumor_createMpileup_scatter.min_BaseQ",
                "normal_createMpileup_scatter.max_depth",
                "normal_createMpileup_scatter.min_BaseQ",
                "TumorCollectIntervalCoverages.disk_size",
                "NormalCollectIntervalCoverages.disk_size",
                "Globals.glob",
                "Sentieon.Globals.glob",
                "AnnotateVCF.Globals.glob",
                "SingleSampleQC.Globals.glob",
                "VariantCallingEvaluation.Globals.glob"
        ]}
    }
    parameter_meta {
        base_file_name: {
            help:"Base file name used for some output files",
            type:"String",
            category: "input_required"
        }
        input_tumor_cram_bam_file: {
            help:"Input tumor BAM/CRAM files",
            type: "File",
            category: "input_required"
        }
        input_tumor_cram_bam_file_index: {
            help:"Input tumor BAI/CRAI index files",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam_file: {
            help:"Input normal BAM/CRAM files",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam_file_index: {
            help:"Input normal BAI/CRAI index files",
            type: "File",
            category: "input_required"
        }
        references: {
        help: "Struct of reference objects holding reference fasta with corresponding fai and dict files",
        type: "References",
        category: "ref_required"
        }
        interval_list: {
            help:"Interval list defining the regions to gather allele frequencies on, recommended value set in the template",
            type: "File",
            category: "param_required"
        }
        num_shards: {
            help:"Shards for scatter-interval-list",
            type: "Int",
            category: "param_advanced"
        }
        scatter_intervals_break: {
            help:"Breaks for scatter interval list, see template for the default",
            type: "Int",
            category: "param_advanced"
        }
        snp_file: {
            help: "Vcf file holding locations of the common variants to calculate pileup statistics on",
            type: "File",
            category: "input_required"
        }
        snp_file_index : {
            help: "Vcf.tbi index file for snp_file",
            type: "File",
            category: "input_required"
        }
        normal_mpileup_override: {
            help: "Pre-calculated mpileup for normal sample",
            type: "File",
            category: "input_optional"
        }
        tumor_mpileup_override: {
            help:"Pre-calculated mpileup for tumor sample",
            type: "File",
            category: "input_optional"
        }
        mapq_override: {
            help:"Reads mapping-quality cutoff for coverage calculation. Default is 1",
            type: "Int",
            category: "param_optional"
        }
        genome_windows: {
            help:"Bed file of the genome binned to equal sized windows",
            type: "File",
            category: "param_required"
        }
        collect_coverage_region: {
            help:"Genomic region to limit the CNV calling to",
            type: "Array[String]",
            category: "param_optional"
        }
        normal_coverage_cpn: {
            help:"Pre-calculated binned coverage for the normal sample in the format needed for controlFREEC (cpn)",
            type: "File",
            category: "input_optional"
        }
        tumor_coverage_cpn: {
            help:"Pre-calculated binned coverage for the tumor sample in the format needed for controlFREEC (cpn)",
            type: "File",
            category: "input_optional"
        }
        normal_sorter_zipped_bed_graph: {
            help:"Pre-calculated bedGraph files containing per-base coverage for the normal sample",
            type: "File",
            category: "input_optional"
        }
        tumor_sorter_zipped_bed_graph: {
            help: "Pre-calculated bedGraph files containing per-base coverage for the tumor sample",
            type: "File",
            category: "input_optional"
        }
        bed_graph_output: {
            help: "Whether to add bed_graph_output to controlFREEC outputs, recommended value set in the template",
            type: "Boolean",
            category: "param_advanced"
        }
        contamination_adjustment: {
            help: "Whether to run controlFREEC with contamination_adjustment option, recommended value set in the template",
            type: "Boolean",
            category: "param_advanced"
        }
        contamination_fraction : {
            help: "a priori known value of tumor sample contamination by normal cells. Default: contamination=0",
            type: "Float",
            category: "param_optional"
        }
        input_format: {
            help: "controlFREEC input format, recommended value set in the template",
            type: "String",
            category: "param_advanced"
        }
        window: {
            help: "The size of the window over which the coverage is aggregated",
            type: "Int",
            category: "param_required"
        }
        chrLenFile_override: {
            help: "Chromosome lengths file focusing controlFREEC regions. file is expected to be tab-delimited where the first column indicates the chromosome name and the second column indicates the chromosome length. By default, the reference.fai file will be used.",
            type: "File",
            category: "param_optional"
        }
        degree: {
            help: "Degree of polynomial for GC normalization. Default is 3",
            type: "Int",
            category: "param_advanced"
        }
        ploidy: {
            help: "Average sample ploidy (if known)",
            type: "Int",
            category: "param_optional"
        }
        max_threads_override: {
            help: "maximal threads for controlFREEC. Default is 8",
            type: "Int",
            category: "param_optional"
        }
        sex: {
            help: "Sample's sex value, should be 'XX' or 'XY'",
            type: "String",
            category: "param_optional"
        }
        gem_mappability_file :{
            help: "Gem file holding mappablity biased regions. ",
            type: "File",
            category: "param_optional"
        }
        CNV_gain_cutoff_override: {
            help: "Gain cutoff for CNV filtering. Default is 1.03",
            type: "Float",
            category: "param_optional"
        }
        CNV_loss_cutoff_override: {
            help: "Loss cutoff for CNV filtering. Default is 0.97",
            type: "Float",
            category: "param_optional"
        }
        high_sensitivity_mode : {
            help: "Whether to run controlFREEC in high sensitivity mode",
            type: "Boolean",
            category: "param_advanced"
        }
        naive_normalization : {
            help: "Whether to run controlFREEC with naive normalization mode (Rather than using FREEC Polynomial fitting)",
            type: "Boolean",
            category: "param_advanced"
        }
        min_cnv_length : {
            help: "Minimum length for reported CNVs. Default is 10,000",
            type: "Int",
            category: "param_required"
        }
        intersection_cutoff: {
            help: "Intersection cutoff with UG-CNV-LCR regions to filter out CNV calls. Default is  0.5",
            type: "Float",
            category: "param_required"
        }
        cnv_lcr_file: {
            help: "UG-CNV-LCR bed file",
            type: "File",
            category: "param_required"
        }
        preemptible_tries_override: {
            help: "Number of preemptible tries",
            type: "Int",
            category: "param_optional"
        }

        tumor_segments: {
            help: "controlFREEC segmentation for tumor sample",
            type: "File",
            category: "output"
        }
        tumor_CNVs : {
            help: "controlFREEC predicted copy number alterations for tumor sample",
            type: "File",
            category: "output"
        }
        controlFREEC_info : {
            help: "controlFREEC run summary",
            type: "File",
            category: "output"
        }
        tumor_ratio_bedgraph: {
            help: "controlFREEC ratios in BedGraph format",
            type: "File",
            category: "output"
        }
        tumor_CNVs_annotated_bed_file : {
            help: "Called CNVs for tumor sample with LCR and LENGTH annotations",
            type: "File",
            category: "output"
        }
        tumor_CNVs_filtered_bed_file : {
            help: "Filtered called CNVs for tumor sample",
            type: "File",
            category: "output"
        }
        tumor_mpileup: {
            help: "mpileup file for tumor sample",
            type: "File",
            category: "output"
        }
        normal_mpileup: {
            help: "mpileup file for normal sample",
            type: "File",
            category: "output"
        }
        tumor_coverage : {
            help: "Coverage file for tumor sample",
            type: "File",
            category: "output"
        }
        normal_coverage : {
            help: "Coverage file for normal sample",
            type: "File",
            category: "output"
        }
        coverage_plot : {
            help: "Coverage plot that shows normalized (log scale) coverage along the genome for the germline and tumor samples",
            type: "File",
            category: "output"
        }
        dup_del_plot : {
            help: "Duplications and deletions figure - showing gains and losses along the genome",
            type: "File",
            category: "output"
        }
        copy_number_plot : {
            help: "Copy-number figure  shows the copy number along the genome",
            type: "File",
            category: "output"
        }
        neutral_AF_plot : {
            help: "Neutral allele frequency plot",
            type: "File",
            category: "output"
        }
        neutral_AF_bed : {
            help: "Neutral allele frequency bed file",
            type: "File",
            category: "output"
        }
    }

    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    if(run_createMpileup){
        call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList{
            input:
              interval_list = interval_list,
              scatter_count = num_shards,
              break_bands_at_multiples_of = scatter_intervals_break,
              dummy_input_for_call_caching = dummy_input_for_call_caching,
              docker = global.gitc_docker,
              gitc_path = global.gitc_jar_path,
              no_address = true,
              monitoring_script = monitoring_script
        }

        scatter (interval in ScatterIntervalList.out){
             call CreateMpileup as tumor_createMpileup_scatter {
                input:
                    input_bam_files = input_tumor_cram_bam_file,
                    input_bam_files_index = input_tumor_cram_bam_file_index,
                    reference_fasta = references.ref_fasta,
                    reference_fai = references.ref_fasta_index,
                    reference_dict = references.ref_dict,
                    min_MapQ = mapq,
                    snp_file = snp_file,
                    snp_file_index = snp_file_index,
                    interval = interval,
                    docker = global.ug_vc_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script,
                    cloud_provider = cloud_provider_override
                }

                call CreateMpileup as normal_createMpileup_scatter{
                input:
                    input_bam_files = input_normal_cram_bam_file,
                    input_bam_files_index = input_normal_cram_bam_file_index,
                    reference_fasta = references.ref_fasta,
                    reference_fai = references.ref_fasta_index,
                    reference_dict = references.ref_dict,
                    min_MapQ = mapq,
                    snp_file = snp_file,
                    snp_file_index = snp_file_index,
                    interval = interval,
                    docker = global.ug_vc_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script,
                    cloud_provider = cloud_provider_override
                }
            }

            Array[File] tumor_pileup_files = tumor_createMpileup_scatter.out_pileup
            Array[File] normal_pileup_files = normal_createMpileup_scatter.out_pileup

            call ConcatFiles as tumor_ConcatMpileupFiles{
                input:
                    files = tumor_pileup_files,
                    out_file_name = basename(input_tumor_cram_bam_file[0])+"_minipileup.pileup",
                    docker = global.ubuntu_docker
            }
            call ConcatFiles as normal_ConcatMpileupFiles{
                input:
                    files = normal_pileup_files,
                    out_file_name = basename(input_normal_cram_bam_file[0])+"_minipileup.pileup",
                    docker = global.ubuntu_docker
            }
    }
    File tumor_pileup = select_first([tumor_ConcatMpileupFiles.out_merged_file,tumor_mpileup_override])
    File normal_pileup = select_first([normal_ConcatMpileupFiles.out_merged_file,normal_mpileup_override])

    if(run_collect_coverage){
        ## Tumor coverage collection
        scatter(tumor_bam_file in zip(input_tumor_cram_bam_file, input_tumor_cram_bam_file_index)) {
            call UGQCTasks.CollectIntervalCoverages as TumorCollectIntervalCoverages{
                input:
                    input_cram_bam = tumor_bam_file.left,
                    input_cram_bam_index = tumor_bam_file.right,
                    references = references,
                    min_mapq = mapq,
                    region = collect_coverage_region,
                    preemptible_tries = preemptible_tries,
                    docker = global.ug_vc_docker,
                    monitoring_script = monitoring_script
                }

            scatter(bw_file in TumorCollectIntervalCoverages.coverage_bw) {
                call BigWigToBedGraph as TumorBigWigToBedGraph{
                    input:
                        bigwig_file = bw_file,
                        genome_windows = genome_windows,
                        genome_file = references.ref_fasta_index,
                        docker = global.ug_vc_docker,
                        no_address = no_address,
                        preemptible_tries = preemptible_tries,
                        monitoring_script = monitoring_script
                }
            }
            Array[File] tumor_bedgraph_files = TumorBigWigToBedGraph.coverage_file
            
            call ConcatFilesBwBedgraphFiles as TumorConcatFilesBwBedgraphFiles{
                input:
                    files = tumor_bedgraph_files,
                    out_file_name = basename(tumor_bedgraph_files[0])+".bedGraph",
                    genome_file = chrLenFile,
                    monitoring_script = monitoring_script,
                    docker = global.ug_vc_docker
            }
        }
        Array[File] t_cov_merged_file = TumorConcatFilesBwBedgraphFiles.out_merged_file

        ## Normal coverage collection
        scatter(normal_bam_file in zip(input_normal_cram_bam_file, input_normal_cram_bam_file_index)) {
            call UGQCTasks.CollectIntervalCoverages as NormalCollectIntervalCoverages{
                input:
                    input_cram_bam = normal_bam_file.left,
                    input_cram_bam_index = normal_bam_file.right,
                    references = references,
                    min_mapq = mapq,
                    region = collect_coverage_region,
                    preemptible_tries = preemptible_tries,
                    docker = global.ug_vc_docker,
                    monitoring_script = monitoring_script
            }
            scatter(bw_file in NormalCollectIntervalCoverages.coverage_bw) {
                call BigWigToBedGraph as NormalBigWigToBedGraph{
                    input:
                        bigwig_file = bw_file,
                        genome_windows = genome_windows,
                        genome_file = references.ref_fasta_index,
                        docker = global.ug_vc_docker,
                        no_address = no_address,
                        preemptible_tries = preemptible_tries,
                        monitoring_script = monitoring_script
                }
            }
            Array[File] normal_bedgraph_files = NormalBigWigToBedGraph.coverage_file

            call ConcatFilesBwBedgraphFiles as NormalConcatFilesBwBedgraphFiles{
                input:
                    files = normal_bedgraph_files,
                    out_file_name = basename(normal_bedgraph_files[0])+".bedGraph",
                    genome_file = chrLenFile,
                    monitoring_script = monitoring_script,
                    docker = global.ug_vc_docker
            }
        }
        Array[File] n_cov_merged_file = NormalConcatFilesBwBedgraphFiles.out_merged_file
    }
    #parsing sorter input
    if(run_bedgraph_to_cpn)
    {
        Array[File] normal_bedgraph=select_first([normal_sorter_zipped_bed_graph,n_cov_merged_file])
        Array[File] tumor_bedgraph=select_first([tumor_sorter_zipped_bed_graph,t_cov_merged_file])
        call BedGraphToCpn as normal_BedGraphToCpn {
            input:
                bedgraph_files = normal_bedgraph,
                genome_windows = genome_windows,
                genome_file = chrLenFile,
                docker = global.ug_vc_docker,
                no_address = no_address,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script
        }
        
        call BedGraphToCpn as tumor_BedGraphToCpn {
            input:
                bedgraph_files = tumor_bedgraph,
                genome_windows = genome_windows,
                genome_file = chrLenFile,
                docker = global.ug_vc_docker,
                no_address = no_address,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script
        }
    }
    File tumor_cpn = select_first([tumor_BedGraphToCpn.coverage_file,tumor_coverage_cpn])
    File normal_cpn = select_first([normal_BedGraphToCpn.coverage_file,normal_coverage_cpn])
    
    call runControlFREEC{
            input:
        sample_name =  base_file_name,
        tumor_bam_file = input_tumor_cram_bam_file_name,
        tumor_mpileup =  tumor_pileup,
        tumor_cpn = tumor_cpn,
        normal_bam_file = input_normal_cram_bam_file_name,
        normal_mpileup = normal_pileup,
        normal_cpn = normal_cpn,
        reference_fasta = references.ref_fasta,
        chrLenFile = chrLenFile,
        input_format = input_format,
        snp_file = snp_file,
        snp_file_index = snp_file_index,
        high_sensitivity_mode = high_sensitivity_mode,
        naive_normalization = naive_normalization,
        bed_graph_output = bed_graph_output,
        contamination_adjustment = contamination_adjustment,
        contamination_fraction = contamination_fraction,
        ploidy = ploidy,
        gem_mappability_file = gem_mappability_file,
        maxThreads = maxThreads,
        sex = sex,
        window =  window,
        degree = degree,
        docker = global.ug_control_freec_docker,
        no_address = no_address,
        preemptible_tries = preemptible_tries,
        monitoring_script = monitoring_script
        }

    call FilterControlFREECCnvs {
        input:
        sample_name = base_file_name,
        CNV_calls = runControlFREEC.tumor_CNVs,
        segments_file = runControlFREEC.segments_file,
        gain_cutoff = CNV_gain_cutoff,
        loss_cutoff = CNV_loss_cutoff,
        high_sensitivity_mode = high_sensitivity_mode,
        min_cnv_length = min_cnv_length,
        intersection_cutoff = intersection_cutoff,
        cnv_lcr_file = cnv_lcr_file,
        tumor_coverage_cpn=tumor_cpn,
        normal_coverage_cpn=normal_cpn,
        ploidy = runControlFREEC.ploidy_value,
        tumor_mpileup = tumor_pileup,
        docker = global.ugbio_cnv_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
        }
    
    output {
         File tumor_mpileup = tumor_pileup
         File normal_mpileup = normal_pileup
         File tumor_coverage = tumor_cpn
         File normal_coverage = normal_cpn
         File tumor_segments = runControlFREEC.segments_file
         File controlFREEC_info = runControlFREEC.controlFREEC_info
         File tumor_ratio_bedgraph = runControlFREEC.tumor_ratio_bedgraph
         File tumor_CNVs_annotated_bed_file = FilterControlFREECCnvs.sample_cnvs_bed
         File tumor_CNVs_filtered_bed_file =  FilterControlFREECCnvs.sample_cnvs_filtered_bed
         File coverage_plot = FilterControlFREECCnvs.coverage_plot
         File dup_del_plot = FilterControlFREECCnvs.dup_del_plot
         File copy_number_plot = FilterControlFREECCnvs.copy_number_plot
         File neutral_AF_plot = FilterControlFREECCnvs.neutral_AF_plot
         File neutral_AF_bed = FilterControlFREECCnvs.neutral_AF_bed
    }
}
task CreateMpileup {
    input{

        Array[File] input_bam_files
        Array[File] input_bam_files_index
        File reference_fasta
        File reference_fai
        File reference_dict
        Int max_depth = 8000
        Int min_MapQ
        Int min_BaseQ = 0
        File snp_file
        File snp_file_index
        File interval
        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
        String? cloud_provider
    }
    Int disk_size = ceil(size(reference_fasta,"GB") + size(snp_file,"GB") + 50)
    String base_input_name = basename(input_bam_files[0])
    Boolean is_aws = defined(cloud_provider)


    parameter_meta {
      input_bam_files: {
          localization_optional: true
      }
    }

command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        set -eo pipefail
        conda activate genomics.py3

    
    # if ~{is_aws}
    # then 
    #     inputs_string="~{sep=' ' input_bam_files}"
    # else
    echo "DEBUG start PrintReads $(date)"
    # download only the region of interval
    gatk --java-options "-Xms1G" PrintReads \
        -I ~{sep=' -I ' input_bam_files} \
        -O input.bam \
        -L ~{interval} \
        -R ~{reference_fasta}
    inputs_string="input.bam"
    echo "DEBUG end PrintReads $(date)"
    # fi

    echo "inputs_string:"
    echo $inputs_string
    echo "DEBUG start intersect snp_file with current interval $(date)"
    # intersect snp_file with current interval
    gatk IntervalListToBed -I ~{interval} -O interval.bed
    bcftools view -R interval.bed -O z -o out.vcf.gz ~{snp_file}
    tabix out.vcf.gz
    echo "DEBUG end intersect snp_file with current interval $(date)"

    echo "DEBUG start mpileup $(date)"
    #caclulate mpileup for current interval
    samtools mpileup -f ~{reference_fasta} \
    -d ~{max_depth} \
    -Q ~{min_BaseQ} \
    -q ~{min_MapQ} \
    -l out.vcf.gz \
    $inputs_string \
    >  ~{base_input_name}_minipileup.pileup
    echo "DEBUG end mpileup $(date)"

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu:1
    }
    output {
        File out_pileup = "~{base_input_name}_minipileup.pileup"
        File monitoring_log = "monitoring.log"
    }
}

task ConcatFiles{
    input{
        Array[File] files
        String out_file_name
        String docker
    }
    Int disk_size = ceil(2 * size(files,"GB") + 2)
    command <<<
        cat ~{sep=" " files} > ~{out_file_name}
    >>>

    runtime {
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        cpu:1
    }
    output{
        File out_merged_file = "~{out_file_name}"
    }
}


task ConcatFilesBwBedgraphFiles{
    input{
        Array[File] files
        File genome_file
        String out_file_name
        String docker
        File monitoring_script
    }
    Int disk_size = ceil(2 * size(files,"GB") + 2)
    
    command {
bash ~{monitoring_script} | tee monitoring.log >&2 &
source ~/.bashrc
set -eo pipefail
conda activate genomics.py3
echo ~{out_file_name}
echo -e "~{sep="\n" files}"
echo -e "~{sep="\n" files}" > files.txt
cat ~{genome_file} | cut -f1
cat ~{genome_file} | cut -f1 > chroms_order.txt
python3 <<CODE
import pandas as pd
import os
df_files = pd.read_csv("files.txt",header=None)
df_files.columns = ["filename"]
df_files["chrom"] = df_files["filename"].apply(lambda x: os.path.basename(x).split(".")[1])
df_chrom_order=pd.read_csv("chroms_order.txt",header=None)
df_chrom_order.columns = ["chrom"]
df_sorted_files = df_chrom_order.merge(df_files, how='inner', on='chrom')
cmd = "cat " + " ".join(df_sorted_files["filename"].tolist()) + " > ~{out_file_name}"
os.system(cmd)
CODE
}
    output{
        File out_merged_file = "~{out_file_name}"
    }
    runtime {
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        cpu:1
    }
    
}


task BigWigToBedGraph{
    input{
        File bigwig_file
        File genome_windows
        File genome_file
        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
    }
    Int disk_size = ceil(size(bigwig_file,"GB") + size(genome_windows,"GB") +  5)
    String base_file_name = basename(bigwig_file, ".bw")
    command<<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        set -eo pipefail
        conda activate ucsc

        bigWigToBedGraph ~{bigwig_file} ~{base_file_name}.bedGraph
        echo ~{bigwig_file} 
        echo ~{base_file_name}.bedGraph
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File coverage_file = "~{base_file_name}.bedGraph"
        File monitoring_log = "monitoring.log"
    }
}

task BedGraphToCpn{
    input{
        Array[File] bedgraph_files
        File genome_windows
        File genome_file

        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
    }
    Int disk_size = ceil((2*size(bedgraph_files,"GB")) + (2*8*size(bedgraph_files,"GB")) + size(genome_windows,"GB") + 5)
    String base_file_name = basename(bedgraph_files[0])
    Boolean multiple_bedgraph_files = length(bedgraph_files) > 1

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        set -eo pipefail
        conda activate genomics.py3

        ##unzip in case of zipped bedgraph
        
        COUNTER=0
        for bedgraph in ~{sep=" " bedgraph_files}; 
        do 
            echo "in_bedgraph:"
            echo $bedgraph
            file_basename=$(basename $bedgraph)
            echo "file_basename:"
            echo $file_basename

            if [[ $bedgraph =~ \.gz$ ]]; 
            then 
                gzip -d -c $bedgraph > $file_basename.bedgraph; 
            else
                cp $bedgraph $file_basename.bedgraph;
            fi

            ## fetch only relevant chromosomes from bedgraph
            awk '{print $1"\t0\t"$2}' ~{genome_file} > ~{genome_file}.bed
            bedtools intersect -a $file_basename.bedgraph -b ~{genome_file}.bed -wa > $file_basename.relevant_chrs.bedgraph
            
            bedtools map -g ~{genome_file} -a ~{genome_windows} -b $file_basename.relevant_chrs.bedgraph -c 4 -o mean | \
            awk '{if($4=="."){print $1"\t"$2"\t"$3"\t"0}else{print $1"\t"$2"\t"$3"\t"$4}}' | \
            grep -v "chrY" > $file_basename.$COUNTER.bedgraph.mean;
            echo $file_basename.$COUNTER.bedgraph.mean
            COUNTER=$[$COUNTER +1]
        done
        
        echo "finished unzipping and converting bedgraph to bedgraph.mean"
        printf "Number of files in array: %d\n" $i
        ls -lrta

        if ~{multiple_bedgraph_files}
        then    
            bedtools unionbedg -i *.bedgraph.mean | \
            awk -v OFS="\t" -F "\t" 'NR>0{sum=0; for(i=4; i<=NF; i++) sum += $i; NF++; $NF=sum } 1' | \
            awk '{print $1"\t"$2"\t"$3"\t"$NF}' > ~{base_file_name}.union.bedgraph
        else
            cp ~{base_file_name}.0.bedgraph.mean ~{base_file_name}.union.bedgraph
        fi            
        echo "finished unionbedg"
        ls -lrta

        cat ~{base_file_name}.union.bedgraph | \
        sed 's/^chr//' > \
        ~{base_file_name}.cpn
        echo "finished converting to cpn"
        ls -lrta
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "4 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File bedgraph_union = "~{base_file_name}.union.bedgraph"
        File coverage_file = "~{base_file_name}.cpn"
        File monitoring_log = "monitoring.log"
    }
}

task runControlFREEC{
        input{
            String sample_name
            String tumor_bam_file
            File tumor_mpileup
            File tumor_cpn
            String normal_bam_file
            File normal_mpileup
            File normal_cpn

            File reference_fasta
            File chrLenFile

            File snp_file
            File snp_file_index

            Boolean high_sensitivity_mode
            String input_format
            Boolean bed_graph_output
            Boolean naive_normalization
            Boolean contamination_adjustment
            Float? contamination_fraction
            Int degree
            Int? ploidy
            File? gem_mappability_file
            Int? maxThreads
            String? sex
            Int window
            String docker
            Boolean no_address
            Int preemptible_tries
            File monitoring_script
        }
    
    Int disk_size = ceil(size(tumor_mpileup,"GB") + size(tumor_cpn,"GB") +
                         size(normal_mpileup,"GB") + size(normal_cpn,"GB") +
                         2*size(reference_fasta,"GB") + size(snp_file,"GB") + size(snp_file_index,"GB") +
                         size(gem_mappability_file,"GB")+ 5)

    String base_input_tumor = basename(tumor_bam_file)
    String base_input_normal = basename(normal_bam_file)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        touch ~{sample_name}.config.txt
        touch ~{base_input_normal}_BAF.txt
        touch ~{base_input_tumor}_segments.txt
        touch ~{base_input_tumor}_BAF.txt
        touch ~{base_input_tumor}_CNVs
        touch ~{base_input_tumor}_info.txt
        touch ~{base_input_tumor}_normal_CNVs
        touch ~{base_input_tumor}_normal_ratio.BedGraph
        touch ~{base_input_tumor}_normal_ratio.txt
        touch ~{base_input_tumor}_ratio.BedGraph
        touch ~{base_input_tumor}_ratio.txt
        touch GC_profile.~{window}bp.cnp

        #split reference to file per chromosome
        mkdir chrFiles_dir
        cd chrFiles_dir
        faidx -x ~{reference_fasta}
        cd ../

        if ~{high_sensitivity_mode}
        then
            ploidy_contam_args="--ploidy 4 --contamination 0.1"
        else
            ploidy_contam_args=~{"--ploidy " + ploidy} ~{" --contamination " + contamination_fraction}
        fi

        #create controlFREEC config file
        python /generate_controlFREEC_config.py \
            --sample_name ~{sample_name} \
            ~{true="--BedGraphOutput TRUE" false='--BedGraphOutput FALSE' bed_graph_output} \
            ~{true="--NaiveNormalization TRUE" false='--NaiveNormalization FALSE' naive_normalization} \
            --chrLenFile ~{chrLenFile}\
            ~{true="--contaminationAdjustment TRUE" false='--contaminationAdjustment FALSE' contamination_adjustment} \
            ~{"--maxThreads " + maxThreads} \
            --window ~{window} \
            --chrFiles chrFiles_dir \
            --degree ~{degree} \
            ~{"--sex " + sex} \
            ~{"--gemMappabilityFile " + gem_mappability_file} \
            --sample_mateFile ~{tumor_bam_file} \
            --sample_mateCopyNumberFile ~{tumor_cpn} \
            --sample_miniPileupFile ~{tumor_mpileup} \
            --sample_inputFormat ~{input_format} \
            --sample_mateOrientation 0 \
            --control_mateFile ~{normal_bam_file} \
            --control_mateCopyNumberFile ~{normal_cpn} \
            --control_miniPileupFile ~{normal_mpileup} \
            --control_inputFormat ~{input_format} \
            --control_mateOrientation 0 \
            --baf_makePileup ~{snp_file} \
            --baf_fastaFile ~{reference_fasta} \
            --baf_SNPfile ~{snp_file} \
            $ploidy_contam_args 

        #run controlFREEC
        /freec -conf ~{sample_name}.config.txt
        
        cat ~{base_input_tumor}_info.txt | grep "Output_Ploidy" | cut -f2 > ploidy_value.txt
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "16 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: maxThreads
    }
    output {
        File sample_config_file = "~{sample_name}.config.txt"
        File normal_BAF = "~{base_input_normal}_BAF.txt"
        File segments_file = "~{base_input_tumor}_segments.txt"
        File tumor_BAF = "~{base_input_tumor}_BAF.txt"
        File tumor_CNVs = "~{base_input_tumor}_CNVs"
        File controlFREEC_info = "~{base_input_tumor}_info.txt"
        File normal_CNVs = "~{base_input_tumor}_normal_CNVs"
        File noraml_ratio_bedgraph = "~{base_input_tumor}_normal_ratio.BedGraph"
        File normal_ratio =  "~{base_input_tumor}_normal_ratio.txt"
        File tumor_ratio_bedgraph =  "~{base_input_tumor}_ratio.BedGraph"
        File tumor_ratio = "~{base_input_tumor}_ratio.txt"
        Int ploidy_value =  read_int("ploidy_value.txt")
        File monitoring_log = "monitoring.log"
    }
}

task FilterControlFREECCnvs {
    input {
        String sample_name
        File segments_file
        File CNV_calls
        Float gain_cutoff
        Float loss_cutoff
        Int min_cnv_length
        Float intersection_cutoff
        File cnv_lcr_file
        File tumor_coverage_cpn
        File normal_coverage_cpn
        Int ploidy
        Boolean high_sensitivity_mode
        File tumor_mpileup
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float CNV_calls_file_size = size(CNV_calls, "GB")
    Float segments_file_size = size(segments_file, "GB")
    Float cnv_lcr_file_size = size(cnv_lcr_file, "GB")
    Float additional_disk = 5
    Int disk_size = ceil(CNV_calls_file_size + segments_file_size + cnv_lcr_file_size + additional_disk)
    String out_annotated_cnv_file = basename(segments_file)+"_annotated.txt"
    String out_cnvs_file = basename(segments_file)+"_CNVs.bed"
    String out_filtered_cnv_file = basename(segments_file)+"_CNVs.filter.bed"
    String mpileup_basename = basename(tumor_mpileup)    
    
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &
        if ~{high_sensitivity_mode}
        then
            #annotate segments to CNVs
            annotate_FREEC_segments \
                --input_segments_file ~{segments_file}\
                --gain_cutoff ~{gain_cutoff} \
                --loss_cutoff ~{loss_cutoff}
        else
            #convert to bedfile
            cat ~{CNV_calls} | sed 's/^/chr/' | grep -v "neutral" | cut -f1-4 > ~{out_cnvs_file}
        fi
              
        #annotate cnvs bed file
        filter_sample_cnvs \
                --input_bed_file ~{out_cnvs_file} \
                --intersection_cutoff ~{intersection_cutoff} \
                --cnv_lcr_file ~{cnv_lcr_file} \
                --min_cnv_length ~{min_cnv_length}
        
        cat ~{normal_coverage_cpn} | awk '{print $1"\t"$2"\t"$2+999"\t"$NF}' | sed 's/^/chr/' > normal_coverage.cpn.bed
        cat ~{tumor_coverage_cpn} | awk '{print $1"\t"$2"\t"$2+999"\t"$NF}' | sed 's/^/chr/' > tumor_coverage.cpn.bed
        
        if ~{high_sensitivity_mode}
        then 
            cat ~{out_filtered_cnv_file} | awk '$4>~{gain_cutoff}' > ~{sample_name}.cnvs.filter.DUP.bed
            cat ~{out_filtered_cnv_file} | awk '$4<~{loss_cutoff}' > ~{sample_name}.cnvs.filter.DEL.bed
        else
            cat  ~{out_filtered_cnv_file} | sed 's/CN//' | awk '$4>~{ploidy}' > ~{sample_name}.cnvs.filter.DUP.bed
            cat  ~{out_filtered_cnv_file} | sed 's/CN//' | awk '$4<~{ploidy}' > ~{sample_name}.cnvs.filter.DEL.bed
        fi

        mkdir CNV_figures
        touch CNV_figures/~{sample_name}.CNV.coverage.jpeg
        touch CNV_figures/~{sample_name}.dup_del.calls.jpeg
        touch CNV_figures/~{sample_name}.CNV.calls.jpeg


        plot_cnv_results \
            --germline_coverage normal_coverage.cpn.bed \
            --tumor_coverage tumor_coverage.cpn.bed \
            --duplication_cnv_calls ~{sample_name}.cnvs.filter.DUP.bed \
            --deletion_cnv_calls ~{sample_name}.cnvs.filter.DEL.bed \
            --sample_name ~{sample_name} \
            --out_directory CNV_figures
        
        plot_FREEC_neutral_AF \
            --mpileup ~{tumor_mpileup} \
            --cnvs_file ~{out_filtered_cnv_file} \
            --sample_name ~{sample_name} \
            --out_directory CNV_figures

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File sample_cnvs_bed = "~{out_cnvs_file}"
        File? sample_annotate_segments = "~{out_annotated_cnv_file}"
        File sample_cnvs_filtered_bed = "~{out_filtered_cnv_file}"
        File coverage_plot = "CNV_figures/~{sample_name}.CNV.coverage.jpeg"
        File dup_del_plot = "CNV_figures/~{sample_name}.dup_del.calls.jpeg"
        File copy_number_plot = "CNV_figures/~{sample_name}.CNV.calls.jpeg"
        File neutral_AF_plot = "CNV_figures/~{mpileup_basename}.freq.SNP.neutral.hist.jpeg"
        File neutral_AF_bed = "CNV_figures/~{mpileup_basename}.freq.SNP.neutral.bed"
        File monitoring_log = "monitoring.log"
    }
    
}