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
        String pipeline_version = "1.11" # !UnusedDeclaration
        String base_file_name

        # input bam files need to be supplied even if coverage and pileup are supplied externally.
        File input_tumor_cram_bam_file
        File input_tumor_cram_bam_file_index
        File input_normal_cram_bam_file
        File input_normal_cram_bam_file_index

        References references

        # Scatter interval list args
        File interval_list
        Int num_shards
        Int scatter_intervals_break
        String dummy_input_for_call_caching = ""

        #Create mpileup args
        File SNPfile
        File SNPfile_index

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
        File? normal_sorter_zipped_bed_graph
        File? tumor_sorter_zipped_bed_graph

        #controlFREEC args
        Boolean BedGraphOutput
        Boolean contaminationAdjustment
        String inputFormat
        Int window
        File? chrLenFile_override
        Int degree
        Int? ploidy
        Int? maxThreads_override
        String? sex
        File? gemMappabilityFile
        Boolean? no_address_override
        Int? preemptible_tries_override

        #Filter CNV calls
        Int min_cnv_length
        Float intersection_cutoff
        File cnv_lcr_file

        # Used for running on other clouds (aws)
        File? monitoring_script_input

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)

        #@wv defined(input_tumor_cram_bam_file) -> (prefix(input_tumor_cram_bam_file_index) == input_tumor_cram_bam_file)
        #@wv defined(input_normal_cram_bam_file) -> (prefix(input_normal_cram_bam_file_index) == input_normal_cram_bam_file)
        #@wv defined(input_tumor_cram_bam_file) -> (suffix(input_tumor_cram_bam_file) in {".bam", ".cram"})
        #@wv defined(input_normal_cram_bam_file) -> (suffix(input_normal_cram_bam_file) in {".bam", ".cram"})

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
    }
    Boolean run_createMpileup = !(defined(normal_mpileup_override))
    Boolean run_bedgraph_to_cpn = !(defined(normal_coverage_cpn))
    Boolean run_collect_coverage = length(select_all([normal_coverage_cpn,normal_sorter_zipped_bed_graph])) == 0

    Int mapq = select_first([mapq_override,1])
    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    File chrLenFile = select_first([chrLenFile_override,references.ref_fasta_index])
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    String input_tumor_cram_bam_file_name = input_tumor_cram_bam_file
    String input_normal_cram_bam_file_name = input_normal_cram_bam_file
    Int maxThreads = select_first([maxThreads_override,8])

    meta {
        description: "Runs single sample somatic CNV calling workflow based on \<a href=\"https://boevalab.inf.ethz.ch/FREEC/\"\>ControlFREEC</a>.\n\nCNVs are called based on both coverage and allele frequencies in the tumor and the matched germline sample.\n\nThe pipeline will gather coverage and allele frequencies, run controlFREEC and filter called CNVs by length and low confidence regions.\n\ncoverage will be calculted based on the input cram/bam. Alternativley, it can recieve coverage input as one of:bedGraph, cpn formats.\n\nAllele frequencies will be calculated based on the input cram/bam and a given vcf file to specify locations. Alternativley, it can recieve precalculated frequencies as mpileup format.\n\nThe pipeline outputs: \n\n&nbsp;&nbsp;- calculated coverage for tumor and normal samples\n\n&nbsp;&nbsp;- calculated mpileup for tumor and normal samples\n\n&nbsp;&nbsp;- called CNVs + filtered called CNVs\n\n&nbsp;&nbsp;- controlFREEC run-summary\n\n"
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
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
            help:"Input tumor BAM/CRAM file",
            type: "File",
            category: "input_required"
        }
        input_tumor_cram_bam_file_index: {
            help:"Input tumor BAI/CRAI index file",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam_file: {
            help:"Input normal BAM/CRAM file",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam_file_index: {
            help:"Input normal BAI/CRAI index file",
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
        SNPfile: {
            help: "Vcf file holding locations of the common variants to calculate pileup statistics on",
            type: "File",
            category: "input_required"
        }
        SNPfile_index : {
            help: "Vcf.tbi index file for SNPfile",
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
            help:"Pre-calculated bedGraph file containing per-base coverage for the normal sample",
            type: "File",
            category: "input_optional"
        }
        tumor_sorter_zipped_bed_graph: {
            help: "Pre-calculated bedGraph file containing per-base coverage for the tumor sample",
            type: "File",
            category: "input_optional"
        }
        BedGraphOutput: {
            help: "Whether to add BedGraphOutput to controlFREEC outputs, recommended value set in the template",
            type: "Boolean",
            category: "param_advanced"
        }
        contaminationAdjustment: {
            help: "Whether to run controlFREEC with contaminationAdjustment option, recommended value set in the template",
            type: "Boolean",
            category: "param_advanced"
        }
        inputFormat: {
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
        maxThreads_override: {
            help: "maximal threads for controlFREEC. Default is 8",
            type: "Int",
            category: "param_optional"
        }
        sex: {
            help: "Sample's sex value, should be 'XX' or 'XY'",
            type: "String",
            category: "param_optional"
        }
        gemMappabilityFile :{
            help: "Gem file holding mappablity biased regions. ",
            type: "File",
            category: "param_optional"
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
                    input_bam_file = input_tumor_cram_bam_file,
                    input_bam_file_index = input_tumor_cram_bam_file_index,
                    reference_fasta = references.ref_fasta,
                    reference_fai = references.ref_fasta_index,
                    reference_dict = references.ref_dict,
                    min_MapQ = mapq,
                    SNPfile = SNPfile,
                    SNPfile_index = SNPfile_index,
                    interval = interval,
                    docker = global.ug_vc_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script
                }

                call CreateMpileup as normal_createMpileup_scatter{
                input:
                    input_bam_file = input_normal_cram_bam_file,
                    input_bam_file_index = input_normal_cram_bam_file_index,
                    reference_fasta = references.ref_fasta,
                    reference_fai = references.ref_fasta_index,
                    reference_dict = references.ref_dict,
                    min_MapQ = mapq,
                    SNPfile = SNPfile,
                    SNPfile_index = SNPfile_index,
                    interval = interval,
                    docker = global.ug_vc_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script
                }
            }

            Array[File] tumor_pileup_files = tumor_createMpileup_scatter.out_pileup
            Array[File] normal_pileup_files = normal_createMpileup_scatter.out_pileup

            call ConcatFiles as tumor_ConcatMpileupFiles{
                input:
                    files = tumor_pileup_files,
                    out_file_name = basename(input_tumor_cram_bam_file)+"_minipileup.pileup",
                    docker = global.ubuntu_docker
            }
            call ConcatFiles as normal_ConcatMpileupFiles{
                input:
                    files = normal_pileup_files,
                    out_file_name = basename(input_normal_cram_bam_file)+"_minipileup.pileup",
                    docker = global.ubuntu_docker
            }
    }
    File tumor_pileup = select_first([tumor_ConcatMpileupFiles.out_merged_file,tumor_mpileup_override])
    File normal_pileup = select_first([normal_ConcatMpileupFiles.out_merged_file,normal_mpileup_override])

    if(run_collect_coverage){
        ## Tumor coverage collection
        call UGQCTasks.CollectIntervalCoverages as TumorCollectIntervalCoverages{
              input:
                input_cram_bam = input_tumor_cram_bam_file,
                input_cram_bam_index = input_tumor_cram_bam_file_index,
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
                out_file_name = basename(input_tumor_cram_bam_file)+".bedGraph",
                genome_file = chrLenFile,
                monitoring_script = monitoring_script,
                docker = global.ug_vc_docker
        }

        ## Normal coverage collection
        call UGQCTasks.CollectIntervalCoverages as NormalCollectIntervalCoverages{
              input:
                input_cram_bam = input_normal_cram_bam_file,
                input_cram_bam_index = input_normal_cram_bam_file_index,
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
                out_file_name = basename(input_normal_cram_bam_file)+".bedGraph",
                genome_file = chrLenFile,
                monitoring_script = monitoring_script,
                docker = global.ug_vc_docker
        }
    }
    #parsing sorter input
    if(run_bedgraph_to_cpn)
    {
        File normal_bedgraph=select_first([normal_sorter_zipped_bed_graph,NormalConcatFilesBwBedgraphFiles.out_merged_file])
        File tumor_bedgraph=select_first([tumor_sorter_zipped_bed_graph,TumorConcatFilesBwBedgraphFiles.out_merged_file])
        call BedGraphToCpn as normal_BedGraphToCpn {
            input:
                bedgraph = normal_bedgraph,
                genome_windows = genome_windows,
                genome_file = chrLenFile,
                docker = global.ug_vc_docker,
                no_address = no_address,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script
        }
        
        call BedGraphToCpn as tumor_BedGraphToCpn {
            input:
                bedgraph = tumor_bedgraph,
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
        inputFormat = inputFormat,
        SNPfile = SNPfile,
        SNPfile_index = SNPfile_index,
        BedGraphOutput = BedGraphOutput,
        contaminationAdjustment = contaminationAdjustment,
        ploidy = ploidy,
        gemMappabilityFile = gemMappabilityFile,
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
        min_cnv_length = min_cnv_length,
        intersection_cutoff = intersection_cutoff,
        cnv_lcr_file = cnv_lcr_file,
        docker = global.ug_vc_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
        }

    output {
         File tumor_mpileup = tumor_pileup
         File normal_mpileup = normal_pileup
         File tumor_coverage = tumor_cpn
         File normal_coverage = normal_cpn
         File tumor_CNVs = runControlFREEC.tumor_CNVs
         File controlFREEC_info = runControlFREEC.controlFREEC_info
         File tumor_ratio_bedgraph = runControlFREEC.tumor_ratio_bedgraph
         File tumor_CNVs_annotated_bed_file = FilterControlFREECCnvs.sample_cnvs_bed
         File tumor_CNVs_filtered_bed_file =  FilterControlFREECCnvs.sample_cnvs_filtered_bed
    }
}
task CreateMpileup {
    input{

        File input_bam_file
        File input_bam_file_index
        File reference_fasta
        File reference_fai
        File reference_dict
        Int max_depth = 8000
        Int min_MapQ
        Int min_BaseQ = 0
        File SNPfile
        File SNPfile_index
        File interval
        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
    }
    Int disk_size = ceil(size(reference_fasta,"GB") + size(SNPfile,"GB") + 50)
    String base_input_name = basename(input_bam_file)

    parameter_meta {
      input_bam_file: {
          localization_optional: true
      }
    }

command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        set -eo pipefail
        conda activate genomics.py3

    echo "DEBUG start PrintReads $(date)"
    # download only the region of interval
    gatk --java-options "-Xms1G" PrintReads \
        -I ~{input_bam_file} \
        -O input.bam \
        -L ~{interval} \
        -R ~{reference_fasta}
    echo "DEBUG end PrintReads $(date)"

    echo "DEBUG start intersect SNPfile with current interval $(date)"
    # intersect SNPfile with current interval
    gatk IntervalListToBed -I ~{interval} -O interval.bed
    bcftools view -R interval.bed -O z -o out.vcf.gz ~{SNPfile}
    tabix out.vcf.gz
    echo "DEBUG end intersect SNPfile with current interval $(date)"

    echo "DEBUG start mpileup $(date)"
    #caclulate mpileup for current interval
    samtools mpileup -f ~{reference_fasta} \
    -d ~{max_depth} \
    -Q ~{min_BaseQ} \
    -q ~{min_MapQ} \
    -l out.vcf.gz \
    input.bam \
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
echo -e "~{sep="\n" files}" > files.txt
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
        File bedgraph
        File genome_windows
        File genome_file

        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
    }
    Int disk_size = ceil(8*size(bedgraph,"GB") + size(genome_windows,"GB") + 5)
    String base_file_name = basename(bedgraph)
    
    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        set -eo pipefail
        conda activate genomics.py3

        ##unzip in case of zipped bedgraph
        if [[ ~{bedgraph} =~ \.gz$ ]]; 
        then 
            gzip -d -c ~{bedgraph} > in.bedgraph; 
        else
            cp ~{bedgraph} in.bedgraph;
        fi

        bedtools map -g ~{genome_file} -a ~{genome_windows} -b in.bedgraph -c 4 -o mean | \
        awk '{if($4=="."){print $1"\t"$2"\t"0}else{print $1"\t"$2"\t"$4}}' | \
        grep -v "chrY" | \
        sed 's/^chr//' > \
        ~{base_file_name}.cpn
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

            File SNPfile
            File SNPfile_index

            String inputFormat
            Boolean BedGraphOutput
            Boolean contaminationAdjustment
            Int degree
            Int? ploidy
            File? gemMappabilityFile
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
                         2*size(reference_fasta,"GB") + size(SNPfile,"GB") + size(SNPfile_index,"GB") +
                         size(gemMappabilityFile,"GB")+ 5)

    String base_input_tumor = basename(tumor_bam_file)
    String base_input_normal = basename(normal_bam_file)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        touch ~{sample_name}.config.txt
        touch ~{base_input_normal}_BAF.txt
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

        #create controlFREEC config file
        python /generate_controlFREEC_config.py \
            --sample_name ~{sample_name} \
            ~{true="--BedGraphOutput TRUE" false='--BedGraphOutput FALSE' BedGraphOutput} \
            --chrLenFile ~{chrLenFile}\
            ~{true="--contaminationAdjustment TRUE" false='--contaminationAdjustment FALSE' contaminationAdjustment} \
            ~{"--maxThreads " + maxThreads} \
            --window ~{window} \
            --chrFiles chrFiles_dir \
            --degree ~{degree} \
            ~{"--sex " + sex} \
            ~{"--ploidy " + ploidy} \
            ~{"--gemMappabilityFile " + gemMappabilityFile} \
            --sample_mateFile ~{tumor_bam_file} \
            --sample_mateCopyNumberFile ~{tumor_cpn} \
            --sample_miniPileupFile ~{tumor_mpileup} \
            --sample_inputFormat ~{inputFormat} \
            --sample_mateOrientation 0 \
            --control_mateFile ~{normal_bam_file} \
            --control_mateCopyNumberFile ~{normal_cpn} \
            --control_miniPileupFile ~{normal_mpileup} \
            --control_inputFormat ~{inputFormat} \
            --control_mateOrientation 0 \
            --baf_makePileup ~{SNPfile} \
            --baf_fastaFile ~{reference_fasta} \
            --baf_SNPfile ~{SNPfile}

        #run controlFREEC
        /freec -conf ~{sample_name}.config.txt
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
        File tumor_BAF = "~{base_input_tumor}_BAF.txt"
        File tumor_CNVs = "~{base_input_tumor}_CNVs"
        File controlFREEC_info = "~{base_input_tumor}_info.txt"
        File normal_CNVs = "~{base_input_tumor}_normal_CNVs"
        File noraml_ratio_bedgraph = "~{base_input_tumor}_normal_ratio.BedGraph"
        File normal_ratio =  "~{base_input_tumor}_normal_ratio.txt"
        File tumor_ratio_bedgraph =  "~{base_input_tumor}_ratio.BedGraph"
        File tumor_ratio = "~{base_input_tumor}_ratio.txt"
        File monitoring_log = "monitoring.log"
    }
}

task FilterControlFREECCnvs {
    input {
        String sample_name
        File CNV_calls
        Int min_cnv_length
        Float intersection_cutoff
        File cnv_lcr_file
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float CNV_calls_file_size = size(CNV_calls, "GB")
    Float cnv_lcr_file_size = size(cnv_lcr_file, "GB")
    Float additional_disk = 5
    Int disk_size = ceil(CNV_calls_file_size + cnv_lcr_file_size + additional_disk)

    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        conda activate genomics.py3

        #convert to bedfile
        cat ~{CNV_calls} | sed 's/^/chr/' | grep -v "neutral" | cut -f1-4 > ~{sample_name}.cnvs.bed

        #annotate cnvs bed file
        python /VariantCalling/ugvc filter_sample_cnvs \
                --input_bed_file ~{sample_name}.cnvs.bed \
                --intersection_cutoff ~{intersection_cutoff} \
                --cnv_lcr_file ~{cnv_lcr_file} \
                --min_cnv_length ~{min_cnv_length}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 1
    }
    output {
        File sample_cnvs_bed = "~{sample_name}.cnvs.annotate.bed"
        File sample_cnvs_filtered_bed = "~{sample_name}.cnvs.filter.bed"
        File monitoring_log = "monitoring.log"
    }
}