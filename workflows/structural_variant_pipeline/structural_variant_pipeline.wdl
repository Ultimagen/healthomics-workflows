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
#

# DESCRIPTION
# Runs Structural variant pipeline
# This pipeline supports germline and somatic modes
# The input of that pipeline is cram files and the output is vcf file
# The steps of the pipeline are as following:
# create an assembly file out of the cram files
# run UA alingnment on that
# fix the UA alignment which are secondarily mapped to decoy or with low mapq
# Run gridss.IdentifyVariants and gridss.AnnotateVariants
# Run R script / GRIPSS for filtering and linkage the variants
#
# CHANGELOG
# 1.11.0 Pipeline optimized for the AWS
# 1.10.4 Improved support for duplication and insertion identification
# 1.10.3 Supports tumor only somatic calling, bug fixes in the assembly
# 1.10.1 Initial revision

import 'tasks/structs.wdl'
import 'tasks/globals.wdl' as Globals
import 'tasks/alignment_tasks.wdl' as AlignmentTasks
import "tasks/general_tasks.wdl" as UGGeneralTasks

workflow SVPipeline {
    input {
        # Workflow args
        String pipeline_version = "1.27.2" # !UnusedDeclaration

        String base_file_name
        Array[File] input_germline_crams = []
        Array[File] input_germline_crams_indexes = []
        Array[File] input_tumor_crams = []
        Array[File] input_tumor_crams_indexes = []
        References references
        UaParameters ua_references
        GiraffeReferences? giraffe_parameters
        File wgs_calling_interval_list
        Int min_base
        Int min_mapq
        Int max_reads_per_partition
        Int max_reads_per_working_area
        Int? max_num_haps
        Int realign_mapq
        String? min_indel_sc_size_to_include
        String? min_mismatch_count_to_include
        Int homopolymer_length
        String config_file_string
        File? blacklist_bed
        Boolean is_somatic
        String reference_name
        Boolean run_ua
        Boolean run_giraffe
        String? prefilter_query
        
        String? gridss_metrics_interval
        # gripss optional args
        File? pon_sgl_file
        File? pon_sv_file
        File? repeat_mask_file
        File? known_hotspot_file
        Int? min_normal_coverage
        String? exclude_filters
        Boolean symbolic_vcf_format

        Int num_shards
        Boolean no_address = true
        Int? preemptible_tries_override
        Int? create_assembly_memory_override
        Int? rematching_memory_override
        Int? annotate_variants_cpu_override
        Int? annotate_variants_memory_override
        Int? convert_vcf_format_memory_override
        Int? germline_link_variants_memory_override
        # Used for running on other clouds (aws)
        String? cloud_provider_override
        Int scatter_intervals_break # Maximal resolution for scattering intervals
        String dummy_input_for_call_caching = ""
        File? monitoring_script_input
        Boolean create_md5_checksum_outputs = false

        # Winval validations
        #@wv min_base >= 0
        #@wv min_mapq >= 0
        #@wv max_num_haps >= 0
        #@wv num_shards > 0
        #@wv realign_mapq >= 0
        #@wv max_reads_per_partition > 0
        #@wv scatter_intervals_break > 0
        #@wv reference_name in {"38","19"}
        #@wv run_giraffe -> defined(giraffe_parameters)
        # tumor + germline
        #@wv is_somatic -> defined(input_tumor_crams) and len(input_tumor_crams)>0
        #@wv is_somatic -> defined(input_tumor_crams_indexes) and len(input_tumor_crams_indexes)>0
        # germline only
        #@wv not is_somatic -> len(input_germline_crams)>0
        #@wv not is_somatic -> len(input_germline_crams_indexes)>0
        #@wv not is_somatic -> len(input_tumor_crams)==0
        #@wv not is_somatic -> len(input_tumor_crams_indexes)==0
        #@wv prefix(input_germline_crams_indexes) == input_germline_crams
        #@wv prefix(input_tumor_crams_indexes) == input_tumor_crams

    }
    meta {
        description : "Runs Structural variant pipeline\nThis pipeline supports germline and somatic modes\nThe input of that pipeline is cram files and the output is vcf file\nThe steps of the pipeline are as following:\n-Create an assembly file out of the cram files\n-Run UA alingnment on that\n-Fix the UA alignment which are secondarily mapped to decoy or with low mapq\n-Run gridss.IdentifyVariants and gridss.AnnotateVariants\n-Run R script / GRIPSS for filtering and linkage the variants\n\n<b>When Running in AWS HealthOmics this pipeline should run with [dynamic storage](https://docs.omics.ai/products/workbench/engines/parameters/aws-healthomics#storage_type-dynamic-or-static)</b>"
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "no_address",
            "preemptible_tries_override",
            "dummy_input_for_call_caching",
            "monitoring_script_input",
            "CreateAssembly.no_address",
            "CreateAssembly.memory_override",
            "MergeBams.disk_size",
            "AlignWithUA.disk_size",
            "AlignWithUA.cpu",
            "AlignWithUA.cache_tarball",
            "IdentifyVariants.no_adress",
            "MergeVCFs.disk_size",
            "MergeVCFs.gitc_path",
            "AlignWithUA.v_aware_vcf",
            "IdentifyVariants.input_crams",
            "IdentifyVariants.input_crams_indexes",
            "Glob.glob",
            "Sentieon.Globals.glob",
            "AnnotateVCF.Globals.glob",
            "SingleSampleQC.Globals.glob",
            "VariantCallingEvaluation.Globals.glob",
            "aws_AnnotateVariants.interval",
            "AnnotateVariants.germline_metrics",
            "AnnotateVariants.tumor_metrics",
            "AnnotateVariants.assembly_metrics",
            "AnnotateVariants.cpu_override",
            "AnnotateVariants.memory_override",
            "AlignWithUA.memory_gb",
            "ConvertToUbam.cache_tarball",
            "ConvertToUbam.disk_size",
            "AlignWithGiraffe.in_prefix_to_strip",
            "SortGiraffeAlignment.disk_size",
            "SortGiraffeAlignment.gitc_path",
            "IndexGiraffeAlignment.disk_size",
            "MergeMd5sToJson.output_json"
            ]}
    }
    parameter_meta {
        base_file_name: {
        help: "Base file name for the output files (to be used as the prefix)",
        type: "string",
        category: "required"
        }
        input_germline_crams: {
        type: "File",
        help: "Input CRAM file for the germline or matched normal sample; optinal for supporting somatic calling tumor only, default []",
        category: "optional"
        }
        input_germline_crams_indexes: {
        type: "File",
        help: "Input CRAM index for the germline or matched normal sample; optinal for supporting somatic calling tumor only",
        category: "optional"
        }
        input_tumor_crams: {
            type: "File",
            help: "Input CRAM file for the tumor (in case of matched T/N calling)",
            category: "optional"
        }
        input_tumor_crams_indexes: {
            type: "File",
            help: "Input CRAM index for the tumor (in case of matched T/N calling)",
            category: "optional"
        }
        references: {
            type: "References",
            help: "Reference files: fasta, dict and fai, recommended value set in the template",
            category: "required"
        }
        ua_references: {
            type: "UaParameters",
            help: "UAReference files: ua_index, ref_alt, v_aware_alignment_flag and ua_extra_args, recommended value set in the template",
            category: "required"
        }
        giraffe_parameters: {
            type: "GiraffeReferences",
            help: "vg giraffe index files to improve haplotype interpretation using population graphs",
            category: "optional"
        }
        wgs_calling_interval_list: {
            type: "File",
            help: "interval list defining the region to perform variant calling on, recommended value set in the template",
            category: "required"
        }
        min_base: {
            type: "Int",
            help: "Assembly parameter: Minimum base quality for using in DeBruijn graph construction. Default value in template",
            category:"required"
        }
        min_mapq: {
            type: "Int",
            help: "Assembly parameter: Minimum mapping quality. Default value in template",
            category:"required"
        }
        min_indel_sc_size_to_include: {
            type: "String",
            help: "Assembly parameter: Minimum size of an indel and soft-clipping in the read to include the read in the assembly. ;-separated between samples",
            category:"optional"
        }
        min_mismatch_count_to_include: {
            type: "String",
            help: "Assembly parameter: Minimal number of counts to require to include the read in the assembly. ;-separated between samples",
            category:"optional"

        }
        max_num_haps: {
            type: "Int",
            help: "Assembly parameter: Maximum number of haplotypes showing an evidence of SV to report",
            category:"required"
        }
        realign_mapq: {
            type: "Int",
            help: "Realignment parameter: Below this value we skip realignment on the supplementary alignment",
            category:"required"
        }
        homopolymer_length: {
            type: "Int",
            help: "Realignment parameter: do realignment on homopolymeres longer than this value",
            category:"required"
        }
        max_reads_per_partition: {
            type: "Int",
            help: "Assembly parameter: Maximal number of reads that are stored in memory when analyzing an active region",
            category:"advanced"
        }
        max_reads_per_working_area: {
            type: "Int",
            help: "Rematching parameter: Maximal number of reads that are stored in memory when rematching reads to haplotypes (similar to max_reads_per_partition in assembly)",
            category:"advanced"
        }
        prefilter_query: {
            type: "String",
            help: "Expression (in bcftools view format) to filter the variants before annotation",
            category:"optional"
        }
        gridss_metrics_interval: {
            type: "String",
            help: "Interval for collecting gridss metrics",
            category:"optional"
        }
        config_file_string: {
            type: "String",
            help: "Gridss config file content",
            category:"advanced"
        }
        blacklist_bed: {
            type: "File",
            help: "Gridss blacklist file",
            category:"optional"
        }
        is_somatic: {
            type: "Boolean",
            help: "run in somatic mode or in germline mode",
            category:"required"
        }

        reference_name: {
            type: "Boolean",
            help: "Can be 38 or 19",
            category:"required"
        }
        run_ua: {
            type: "Boolean",
            help: "Whether to run UA realignment on the output of the assembly (helps resolving some deletions) or not",
            category:"required"
        }
        run_giraffe: {
            type: "Boolean",
            help: "Whether to run Giraffe haplotype aware alignment or not",
            category:"optional"
        }

        pon_sgl_file: {
            type: "File",
            help: "gripss paramter: Panel of normals for single end breakend (partially resolved) calls. Note that the default value is in template",
            category:"optional"
        }
        pon_sv_file: {
            type: "File",
            help: "gripss paramter: panel of normals for breakpoint (fully resolved) calls. Note that the default value is in template",
            category:"optional"
        }
        repeat_mask_file: {
            type: "File",
            help: "gripss paramter: Repeat mask file. Note that the default value is in template",
            category:"optional"
        }
        known_hotspot_file: {
            type: "File",
            help: "gripss paramter: Known locations that are hot spot for SVs (see https://github.com/hartwigmedical/hmftools/tree/master/linx), filtered less stringently",
            category:"optional"
        }
        min_normal_coverage: {
            type: "Int",
            help: "gripss paramter: Minimum coverage in the normal sample to determine somatic status. Default value:8",
            category:"optional"
        }
        exclude_filters: {
                type: "String",
                help: "gripss paramter: Exclude filters from the output vcf, separated by ;",
                category: "optional"
        }
        symbolic_vcf_format: {
            type: "Boolean",
            help: "Whether to convert the output vcf to the region format or not, default True",
            category:"optional"
        }
        cloud_provider_override: {
            type: "String",
            help: "Cloud provider to use for the workflow. Currently supported: aws, gcp default: gcp",
            category: "optional"
        }
        num_shards: {
            type: "Int",
            help: "Relevant for scatter tasks, which are CreateAssembly and gridss.AnnotateVariants",
            category:"required"
        }
        scatter_intervals_break: {
            type: "Int",
            help: "Maximal resolution for scattering intervals",
            category: "advanced"
        }
        create_assembly_memory_override: {
            type: "Int",
            help: "memory override for create_assembly task",
            category: "advanced"
        }
        rematching_memory_override: {
            type: "Int",
            help: "memory override for rematching task",
            category: "advanced"
        }
        annotate_variants_memory_override: {
            type: "Int",
            help: "memory override for annotate_variants task",
            category: "advanced"
        }
        convert_vcf_format_memory_override: {
            type: "Int",
            help: "memory override for convert_vcf_format task",
            category: "advanced"
        }
        germline_link_variants_memory_override: {
            type: "Int",
            help: "memory override for germline_link_variants task",
            category: "advanced"
        }
        annotate_variants_cpu_override: {
            type: "Int",
            help: "cpu override for annotate_variants task",
            category: "advanced"
        }
        create_md5_checksum_outputs: {
           help: "Create md5 checksum for requested output files",
           type: "Boolean",
           category: "input_optional"
        }
        output_vcf: {
            type: "File",
            help: "Final VCF",
            category: "output"
        }
        output_vcf_index: {
            type: "File",
            help: "Final VCF index",
            category: "output"
        }
        assembly: {
            type: "File",
            help: "Raw assembly - before the realignment",
            category: "output"
        }
        assembly_index: {
            type: "File",
            help: "Raw assembly - before the realignment - index",
            category: "output"
        }

        realigned_assembly: {
            type: "File",
            help: "Assembly output after UA realingment",
            category: "output"
        }
        realigned_assembly_index: {
            type: "File",
            help: "Assembly output index after UA realingment",
            category: "output"
        }
        giraffe_realigned_assembly: {
            type: "File",
            help: "Giraffe realigned assembly output",
            category: "output"
        }
        giraffe_realigned_assembly_index: {
            type: "File",
            help: "Giraffe realigned assembly output index",
            category: "output"
        }
        converted_vcf: {
            type: "File",
            help: "Final VCF file in the region (non-breakend) format",
            category: "output"
        }
        converted_vcf_index: {
            type: "File",
            help: "Final VCF index file in the region (non-breakend) format",
            category: "output"
        }
        annotated_vcf_out: {
            type: "File",
            help: "Annotated VCF file, before GRIPSS or GermlineLinkVariants",
            category: "output"
        }
        annotated_vcf_index_out: {
            type: "File",
            help: "Annotated VCF index file",
            category: "output"
        }
        md5_checksums_json: {
            help: "json file that will contain md5 checksums for requested output files",
            type: "File",
            category: "output"
        }
    }
    Int preemptibles = select_first([preemptible_tries_override, 1])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false

    if(length(input_germline_crams)>0){
        call UGGeneralTasks.ExtractSampleNameFlowOrder as ExtractSampleNameGermline {
            input:
                input_bam = input_germline_crams[0],
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                docker = global.broad_gatk_docker,
                references = references,
                no_address = no_address,
                cloud_provider_override = cloud_provider_override
        }
    }

    if (is_somatic){
       call UGGeneralTasks.ExtractSampleNameFlowOrder as ExtractSampleNameTumor {
            input:
                input_bam = input_tumor_crams[0],
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                docker = global.broad_gatk_docker,
                references = references,
                no_address = no_address,
                cloud_provider_override = cloud_provider_override
        }
    }

    call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList {
        input:
            interval_list = wgs_calling_interval_list,
            scatter_count = num_shards,
            break_bands_at_multiples_of = scatter_intervals_break,
            dummy_input_for_call_caching = dummy_input_for_call_caching,
            docker = global.broad_gatk_docker,
            no_address = no_address,
            monitoring_script = monitoring_script
    }

    scatter (interval in ScatterIntervalList.out) {
        call CreateAssembly {
            input:
                monitoring_script = monitoring_script,
                input_germline_crams = input_germline_crams,
                input_germline_crams_indexes = input_germline_crams_indexes,
                input_tumor_crams = input_tumor_crams,
                input_tumor_crams_indexes = input_tumor_crams_indexes,
                output_prefix = base_file_name +"_assembly",
                make_examples_docker = global.ug_make_examples_docker,
                references = references,
                interval = interval,
                min_base = min_base,
                min_mapq = min_mapq,
                max_num_haps = max_num_haps,
                max_reads_per_partition = max_reads_per_partition,
                min_sc_indel_size = min_indel_sc_size_to_include,
                min_mismatch_count = min_mismatch_count_to_include, 
                is_somatic = is_somatic,
                number_of_shards = ScatterIntervalList.interval_count,
                cloud_provider_override = cloud_provider_override,
                no_address = no_address,
                preemptible_tries = preemptibles,
                memory_override = create_assembly_memory_override
        }
    }

    call UGGeneralTasks.MergeBams {
        input:
            inputs = CreateAssembly.assembly,
            output_prefix = base_file_name + '_merged_assembly',
            docker = global.broad_gatk_docker,
            monitoring_script=monitoring_script,
            preemptible_tries = preemptibles,
            no_address = no_address
    }

    if(run_ua){
      call AlignmentTasks.AlignWithUA as AlignWithUA{
            input:
                input_bams = CreateAssembly.assembly,
                output_bam_basename = base_file_name +"_assembly_file_ua_aligned",
                ua_index = select_first([ua_references.ua_index]),
                ref_alt = ua_references.ref_alt,
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                ua_docker = global.ua_docker,
                extra_args = ua_references.ua_extra_args,
                use_v_aware_alignment =  ua_references.v_aware_alignment_flag,
                no_address = no_address
        }

    }

    if(run_giraffe){
        call AlignmentTasks.ConvertCramOrBamToUBam as ConvertToUbam {
            input:
                monitoring_script = global.monitoring_script, # !FileCoercion
                input_file = MergeBams.output_bam,
                base_file_name = base_file_name,
                preemptible_tries = preemptibles,
                docker = global.ug_gatk_picard_docker,
                no_address = no_address
        }

        call AlignmentTasks.SamToFastqAndGiraffeAndMba as AlignWithGiraffe {
            input:
                input_bam = ConvertToUbam.unmapped_bam,
                output_bam_basename = base_file_name,
                giraffe_references = select_first([giraffe_parameters]),
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                no_address = no_address,
                docker = global.giraffe_docker
        }
    }

    call ChooseBestAlignment {
        input:
            inputs = select_all([MergeBams.output_bam, AlignWithUA.ua_output_bam, AlignWithGiraffe.output_bam]),
            tie_breaker_idx = 2, # index of giraffe alignment (if exists), might need smarter way to change
            outptut_prefix = base_file_name +"_assembly_realigned",
            realign_mapq = realign_mapq,
            gridss_docker = global.gridss_docker,
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script,
            no_address = no_address
    }

    call LongHomopolymersAlignment {
        input:
            input_bam = ChooseBestAlignment.realigned_assembly,
            input_bam_index = ChooseBestAlignment.realigned_assembly_index,
            output_bam_basename = base_file_name +"_assembly_realigned_long_homopolymers_aligned",
            references = references,
            homopolymer_length = homopolymer_length,
            gridss_docker = global.gridss_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptibles,
            no_address = no_address

    }

    # Find for each read the best haplotype that it matches to. 
    scatter (interval in ScatterIntervalList.out) {
        # First we do this in chunks
        call MatchReadsToHaplotypes {
            input: 
                germline_reads = input_germline_crams,
                germline_reads_indexes = input_germline_crams_indexes,
                tumor_reads = input_tumor_crams,
                tumor_reads_indexes = input_tumor_crams_indexes,
                assembly_bam = LongHomopolymersAlignment.realigned_assembly,
                assembly_bam_index = LongHomopolymersAlignment.realigned_assembly_index,
                base_file_name = base_file_name,
                is_somatic = is_somatic,
                min_indel_sc_size_to_include = min_indel_sc_size_to_include,
                min_mismatch_count_to_include = min_mismatch_count_to_include,
                min_mapq = min_mapq,
                max_reads_per_working_area = max_reads_per_working_area,
                number_of_shards = ScatterIntervalList.interval_count,
                rematching_interval = interval,
                reference = references,
                docker = global.rematching_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address,
                cloud_provider_override = cloud_provider_override,
                memory_override = rematching_memory_override
        }
    }

    # combine the chunks, choosing best scoring assignment for read per chunk and sorting by haplotype
    call MergeMatchReadsFiles {
        input:
            matched_reads_files = MatchReadsToHaplotypes.read_haplotype_assignment,
            base_file_name = base_file_name,
            docker = global.ubuntu_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptibles,
            no_address = no_address 
    }

    # add the information into  the assembly BAM
    call MergeMatchedReadsIntoAssembly {
        input:
            assembly_bam = LongHomopolymersAlignment.realigned_assembly,
            assembly_bam_index = LongHomopolymersAlignment.realigned_assembly_index,
            matched_reads = MergeMatchReadsFiles.read_haplotype_assignment,
            base_file_name = base_file_name,
            reference = references,
            docker = global.rematching_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptibles,
            no_address = no_address 
    }

    File assembly_file = MergeMatchedReadsIntoAssembly.updated_assembly
    File assembly_file_index = MergeMatchedReadsIntoAssembly.updated_assembly_index

    call IdentifyVariants {
        input:
            sample_names = select_all([ExtractSampleNameTumor.sample_name, ExtractSampleNameGermline.sample_name]),
            references = references,
            assembly = assembly_file,
            assembly_index = assembly_file_index,
            config_file_string = config_file_string,
            output_vcf_prefix = base_file_name,
            blacklist_bed = blacklist_bed,
            gridss_docker = global.gridss_docker,
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script,
            no_address = no_address
    }

    if (defined(prefilter_query)) {
        call PreFilterCandidates{
            input:
                input_vcf = IdentifyVariants.identify_vcf,
                filter_string = select_first([prefilter_query]),
                output_vcf_prefix = base_file_name,
                docker = global.gridss_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }
    }

    if(is_aws)
    {

            if (length(input_germline_crams)>0) {
                call CollectGridssMetrics as germline_CollectGridssMetrics {
                    input:
                        input_cram = input_germline_crams[0],
                        input_cram_index = input_germline_crams_indexes[0],
                        references = references,
                        interval = gridss_metrics_interval,
                        gridss_docker = global.gridss_docker,
                        monitoring_script = monitoring_script,
                        preemptible_tries = preemptibles,
                        no_address = no_address
                }
            }
            if(length(input_tumor_crams)>0) {
                call CollectGridssMetrics as tumor_CollectGridssMetrics {
                    input:
                        input_cram = input_tumor_crams[0],
                        input_cram_index = input_tumor_crams_indexes[0],
                        references = references,
                        interval = gridss_metrics_interval,
                        gridss_docker = global.gridss_docker,
                        monitoring_script = monitoring_script,
                        preemptible_tries = preemptibles,
                        no_address = no_address
                }
            }

            call CollectGridssMetrics as assembly_CollectGridssMetrics {
            input:
                input_cram = assembly_file,
                input_cram_index = assembly_file_index,
                references = references,
                interval = gridss_metrics_interval,
                gridss_docker = global.gridss_docker,
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                no_address = no_address
            }

            call AnnotateVariants as aws_AnnotateVariants {
            input:
                input_vcf = select_first([PreFilterCandidates.filtered_vcf, IdentifyVariants.identify_vcf]),
                input_vcf_index = select_first([PreFilterCandidates.filtered_vcf_index, IdentifyVariants.identify_vcf_index]),
                input_germline_crams = input_germline_crams,
                input_germline_crams_indexes = input_germline_crams_indexes,
                input_tumor_crams = input_tumor_crams,
                input_tumor_crams_indexes = input_tumor_crams_indexes,
                references =references,
                assembly = assembly_file,
                assembly_index = assembly_file_index,
                config_file_string = config_file_string,
                output_vcf_prefix = base_file_name,
                blacklist_bed = blacklist_bed,
                is_somatic = is_somatic,
                number_of_shards = 1,
                germline_metrics = germline_CollectGridssMetrics.out_metrics,
                tumor_metrics = tumor_CollectGridssMetrics.out_metrics,
                assembly_metrics = assembly_CollectGridssMetrics.out_metrics,
                gridss_docker = global.gridss_docker,
                cloud_provider_override = cloud_provider_override,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address,
                cpu_override = annotate_variants_cpu_override,
                memory_override = annotate_variants_memory_override,
            }
    }
    if(!is_aws){
        scatter (interval in ScatterIntervalList.out){
            call AnnotateVariants {
            input:
                input_vcf = select_first([PreFilterCandidates.filtered_vcf, IdentifyVariants.identify_vcf]),
                input_vcf_index = select_first([PreFilterCandidates.filtered_vcf_index, IdentifyVariants.identify_vcf_index]),
                input_germline_crams = input_germline_crams,
                input_germline_crams_indexes = input_germline_crams_indexes,
                input_tumor_crams = input_tumor_crams,
                input_tumor_crams_indexes = input_tumor_crams_indexes,
                references =references,
                interval = interval,
                assembly = assembly_file,
                assembly_index = assembly_file_index,
                config_file_string = config_file_string,
                output_vcf_prefix = base_file_name,
                blacklist_bed = blacklist_bed,
                is_somatic = is_somatic,
                number_of_shards = ScatterIntervalList.interval_count,
                gridss_docker = global.gridss_docker,
                cloud_provider_override = cloud_provider_override,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address,
                memory_override = annotate_variants_memory_override
            }
    }

        call UGGeneralTasks.MergeVCFs as MergeVCFs {
            input:
            input_vcfs = AnnotateVariants.annotated_vcf,
            input_vcfs_indexes = AnnotateVariants.annotated_vcf,
            output_vcf_name = base_file_name + ".ann.vcf.gz",
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script,
            docker = global.broad_gatk_docker,
            no_address = no_address
        }
    }

    File annotated_vcf = select_first([aws_AnnotateVariants.annotated_vcf, MergeVCFs.output_vcf])
    File annotated_vcf_index = select_first([aws_AnnotateVariants.annotated_vcf_index, MergeVCFs.output_vcf_index])

    if(!is_somatic) {
        call GermlineLinkVariants {
            input:
                input_vcf = annotated_vcf,
                input_vcf_index = annotated_vcf_index,
                references = references,
                output_vcf_prefix = base_file_name,
                gridss_docker = global.gridss_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address,
                memory_override = germline_link_variants_memory_override,
        }
    }
    if(is_somatic){
        call SomaticGripss {
            input:
                input_vcf = annotated_vcf,
                input_vcf_index = annotated_vcf_index,
                tumor_sample = select_first([ExtractSampleNameTumor.sample_name]),
                germline_sample = ExtractSampleNameGermline.sample_name,
                references = references,
                reference_name = reference_name,
                output_vcf_prefix = base_file_name,
                pon_sgl_file = pon_sgl_file,
                pon_sv_file = pon_sv_file,
                known_hotspot_file = known_hotspot_file,
                repeat_mask_file = repeat_mask_file,
                min_normal_coverage = min_normal_coverage,
                exclude_filters = exclude_filters,
                gripss_docker = global.gripss_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }

    }

    if(symbolic_vcf_format) {
        call ConvertVcfFormat {
            input:
                input_vcf = select_first([GermlineLinkVariants.linked_vcf, SomaticGripss.gripss_vcf]),
                input_vcf_index = select_first([GermlineLinkVariants.linked_vcf_index, SomaticGripss.gripss_vcf_index]),
                output_vcf_prefix = base_file_name,
                gridss_docker = global.gridss_docker,
                references = references,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address,
                memory_override = convert_vcf_format_memory_override
        }
    }

    File output_vcf_ = select_first([GermlineLinkVariants.linked_vcf, SomaticGripss.gripss_vcf])
    File output_vcf_index_ = select_first([GermlineLinkVariants.linked_vcf_index, SomaticGripss.gripss_vcf_index])
    
    if (create_md5_checksum_outputs) {
        Array[File] output_files = [output_vcf_, output_vcf_index_]
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
        File annotated_vcf_out = annotated_vcf
        File annotated_vcf_index_out = annotated_vcf_index
        File output_vcf = output_vcf_
        File output_vcf_index = output_vcf_index_
        File assembly = MergeBams.output_bam
        File assembly_index = MergeBams.output_bam_index
        File realigned_assembly = assembly_file
        File realigned_assembly_index = assembly_file_index
        File? converted_vcf = ConvertVcfFormat.output_vcf
        File? converted_vcf_index = ConvertVcfFormat.output_vcf_index

        File? md5_checksums_json = MergeMd5sToJson.md5_json
    }
}

task CreateAssembly {
    input {
        Array[File] input_germline_crams
        Array[File] input_germline_crams_indexes
        Array[File] input_tumor_crams
        Array[File] input_tumor_crams_indexes
        Boolean is_somatic
        String output_prefix
        String make_examples_docker
        References references
        File interval
        Int min_base
        Int min_mapq
        String? min_sc_indel_size
        String? min_mismatch_count
        Int? max_num_haps
        Int max_reads_per_partition
        Int number_of_shards
        String? cloud_provider_override
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
        Int? memory_override
    }
    Int disk_size = 3*ceil(size(input_germline_crams,"GB")/number_of_shards +
                    (size(input_tumor_crams, "GB")/number_of_shards)  +
                       size(references.ref_fasta,"GB") + 10)
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false
    Boolean defined_germline = if(length(input_germline_crams)>0) then true else false
    Int mem = select_first([memory_override,4])

    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # convert interval list to bed file
        gatk IntervalListToBed -I ~{interval} -O interval.bed

        if ~{is_aws}
        then
            tumor=~{sep=',' input_tumor_crams}
            tumor_index=~{sep=',' input_tumor_crams_indexes}
            germline=~{sep=',' input_germline_crams}
            germline_index=~{sep=',' input_germline_crams_indexes}

            echo 'Input files are:'
                echo $tumor
                echo $tumor_index
                echo $germline
                echo $germline_index

            if ~{is_somatic}
            then
                input_string="~{true='$tumor;$germline' false='$tumor' defined_germline}"
                input_index_string="~{true='$tumor_index;$germline_index' false='$tumor_index' defined_germline}"
            else
                input_string=$germline
                input_index_string=$germline_index
            fi

            echo 'Input strings are:'
            echo $input_string
            echo $input_index_string

            tool \
            --input "$input_string" \
            --cram-index "$input_index_string" \
            ~{if is_somatic then "--somatic" else ""} \
            --output ~{output_prefix} \
            --ref ~{references.ref_fasta} \
            --bed interval.bed \
            --min-base ~{min_base} \
            --min-mapq ~{min_mapq} \
            --max-reads-per-region ~{max_reads_per_partition} \
            ~{"--min-feature-length '" + min_sc_indel_size + "'"} \
            ~{"--min-mismatch-count '" + min_mismatch_count + "'"} \
            ~{"--max-num-haps " + max_num_haps} \
            --prog \
            --sv \
            --realigned-sam


        else
            if ~{defined_germline}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' input_germline_crams} \
                -O /dev/stdout  \
                -L ~{interval} \
                -R ~{references.ref_fasta} | \
                samtools view -C -T ~{references.ref_fasta} -o input_germline.cram --output-fmt-option embed_ref=1 -

                samtools index input_germline.cram
            fi

            if ~{is_somatic}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' input_tumor_crams} \
                -O /dev/stdout \
                -L ~{interval} \
                -R ~{references.ref_fasta} | \
                samtools view -C -T ~{references.ref_fasta} -o input_tumor.cram --output-fmt-option embed_ref=1 -

                samtools index input_tumor.cram
            fi

            tool \
            ~{if is_somatic then (if defined_germline then "--input input_tumor.cram\\;input_germline.cram" else "--input input_tumor.cram") else "--input input_germline.cram"} \
            ~{if is_somatic then (if defined_germline then "--cram-index input_tumor.cram.crai\\;input_germline.cram.crai" else "--cram-index input_tumor.cram.crai") else "--cram-index input_germline.cram.crai"} \
            ~{if is_somatic then "--somatic" else ""} \
            --output ~{output_prefix} \
            --ref ~{references.ref_fasta} \
            --bed interval.bed \
            --min-base ~{min_base} \
            --min-mapq ~{min_mapq} \
            --max-reads-per-region ~{max_reads_per_partition} \
            ~{"--min-feature-length '" + min_sc_indel_size + "'"} \
            ~{"--min-mismatch-count '" + min_mismatch_count + "'"} \
            ~{"--max-num-haps " + max_num_haps} \
            --prog \
            --interval-target-nreads 10000 \
            --sv \
            --realigned-sam 

        fi

        samtools view -bS ~{output_prefix}_hap_out.sam > ~{output_prefix}_hap_out.bam
        samtools sort ~{output_prefix}_hap_out.bam -o ~{output_prefix}_hap_out_sorted.bam
        samtools index ~{output_prefix}_hap_out_sorted.bam

    >>>
    parameter_meta {
      input_tumor_crams: {
          localization_optional: true
      }
      input_germline_crams: {
          localization_optional: true
      }
    }
    runtime {
        memory: mem + " GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: make_examples_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File assembly = "~{output_prefix}_hap_out_sorted.bam"
        File assembly_index = "~{output_prefix}_hap_out_sorted.bam.bai"
    }
}


task ChooseBestAlignment {
    input {
        Array[File] inputs
        Int tie_breaker_idx
        Int realign_mapq
        String outptut_prefix
        String gridss_docker
        Int preemptible_tries
        File monitoring_script
        Boolean no_address
    }
    Int disk_size = ceil(3*size(inputs,"GB"))+20
    Int cpus = 8
    command <<<
        set -x 
        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # go over the list of inputs and samtools sort -n them (sort by query name)
        # this is required by the choose_best_haplotype_realignment.py script
        alignment_source_flag=""
        for i in ~{sep=" " inputs}
        do
            output_file=${i%.bam}_sorted.bam
            samtools sort -@~{cpus} -n ${i} -o ${output_file}
            alignment_source_flag+=" --alignment_source ${output_file}"
        done

        python3 /opt/gridss/choose_best_haplotype_realignment.py \
            ${alignment_source_flag} \
            --output ~{outptut_prefix}_unsorted.bam \
            --min_mapping_quality ~{realign_mapq} \
            --tie_breaker_idx ~{tie_breaker_idx} \
            --overwrite_mapq 2,3,30
        samtools sort -o ~{outptut_prefix}.bam ~{outptut_prefix}_unsorted.bam
        samtools index ~{outptut_prefix}.bam
    >>>
    runtime {
        memory: "16 GiB"
        cpu: cpus
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File realigned_assembly = "~{outptut_prefix}.bam"
        File realigned_assembly_index = "~{outptut_prefix}.bam.bai"
    }
}

task LongHomopolymersAlignment {
        input {
            File input_bam
            File input_bam_index
            String output_bam_basename
            String gridss_docker
            References references
            Int homopolymer_length
            File monitoring_script
            Int preemptible_tries
            Boolean no_address
        }
        Int disk_size = ceil(4*size(input_bam,"GB") + size(references.ref_fasta,"GB")) + 10
        command <<<
            set -xeo pipefail

            bash ~{monitoring_script} | tee monitoring.log >&2 &

            python3 /opt/gridss/align_long_homopolymers.py \
                --input ~{input_bam} \
                --output ~{output_bam_basename}.bam \
                --reference ~{references.ref_fasta} \
                --homopolymer_length ~{homopolymer_length}

        >>>
        runtime {
            memory: "16 GiB"
            cpu: 8
            disks: "local-disk " + disk_size + " HDD"
            docker: gridss_docker
            preemptible: preemptible_tries
            noAddress: no_address
        }
        output {
            File monitoring_log = "monitoring.log"
            File realigned_assembly = "~{output_bam_basename}.bam"
            File realigned_assembly_index = "~{output_bam_basename}.bam.bai"
        }
}

task MatchReadsToHaplotypes {
    input {
        References reference
        File assembly_bam
        File assembly_bam_index
        Array[File] tumor_reads
        Array[File] tumor_reads_indexes
        Array[File] germline_reads
        Array[File] germline_reads_indexes
        Boolean is_somatic
        Int max_reads_per_working_area
        String? min_indel_sc_size_to_include
        String? min_mismatch_count_to_include
        Int min_mapq
        File rematching_interval
        File monitoring_script
        String base_file_name
        Int number_of_shards
        String docker
        String? cloud_provider_override
        Int preemptible_tries
        Int? memory_override
        Boolean no_address
    }

    Int disk_size = 3*ceil(size(germline_reads,"GB")/number_of_shards +
                    (size(tumor_reads, "GB")/number_of_shards)  +
                       size(reference.ref_fasta,"GB") + 10)
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false
    Boolean defined_germline = if(length(germline_reads)>0) then true else false
    Int mem = select_first([memory_override,4])

    command <<<
        set -o pipefail
        set -ex
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # convert interval list to bed file
        gatk IntervalListToBed -I ~{rematching_interval} -O interval.bed

        if ~{is_aws}
        then
            tumor=~{sep=',' tumor_reads}
            germline=~{sep=',' germline_reads}

            if ~{is_somatic}
            then
                input_string="~{true='$tumor;$germline' false='$tumor' defined_germline}"
            else
                input_string=$germline
            fi

            echo 'Input strings are:'
            echo $input_string

            sv_rematch -b interval.bed \
            -j 2 -t -a -v \
            -l '~{min_indel_sc_size_to_include}' \
            -s '~{min_mismatch_count_to_include}' \
            -m ~{min_mapq} \
            -L ~{max_reads_per_working_area} \
            $input_string \
            ~{assembly_bam} \
            ~{reference.ref_fasta} \
            ~{base_file_name}_rematched_hap.txt

        else
            if ~{defined_germline}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' germline_reads} \
                -O /dev/stdout \
                -L ~{rematching_interval} \
                -R ~{reference.ref_fasta} | \
                samtools view -C -T ~{reference.ref_fasta} -o input_germline.cram --output-fmt-option embed_ref=1 -

                samtools index input_germline.cram
            fi

            if ~{is_somatic}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' tumor_reads} \
                -O /dev/stdout \
                -L ~{rematching_interval} \
                -R ~{reference.ref_fasta} | \
                samtools view -C -T ~{reference.ref_fasta} -o input_tumor.cram --output-fmt-option embed_ref=1 -

                samtools index input_tumor.cram
            fi

            sv_rematch -b interval.bed \
            -j 2 -t -a -v \
            ~{"-l '" + min_indel_sc_size_to_include + "'"} \
            ~{"-s '" + min_mismatch_count_to_include + "'"} \
            -m ~{min_mapq} \
            -L ~{max_reads_per_working_area} \
            ~{if is_somatic then (if defined_germline then "input_tumor.cram\\;input_germline.cram" else "input_tumor.cram") else "input_germline.cram"} \
            ~{assembly_bam} \
            ~{reference.ref_fasta} \
            ~{base_file_name}_rematched_hap.txt
        fi

    >>>

    parameter_meta {
      tumor_reads: {
          localization_optional: true
      }
      germline_reads: {
          localization_optional: true
      }
    }
    runtime {
        memory: mem + " GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File read_haplotype_assignment = "~{base_file_name}_rematched_hap.txt"
    }
}

task MergeMatchReadsFiles { 
    input {
        Array[File] matched_reads_files
        String base_file_name
        String docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    
    Int disk_size = ceil(4*size(matched_reads_files,"GB")) + 10

    command <<<
        set -o pipefail
        set -ex
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        
        # sort by read name, choose the best score, then sort by haplotype
        LC_ALL=C sort -T "." -m -k1,1 -k3,3rn ~{sep=' ' matched_reads_files} | awk '$1 != prev { print; prev = $1 }'  | sort -T "." -k2,2 > ~{base_file_name}.all_matches_by_hap.txt
    >>>

    runtime {
        memory: "4 GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_tries
        noAddress: no_address   
    }    
    output {
        File monitoring_log = "monitoring.log"
        File read_haplotype_assignment = "~{base_file_name}.all_matches_by_hap.txt" 
    }
}

task MergeMatchedReadsIntoAssembly { 
    input {
        File matched_reads
        File assembly_bam
        File assembly_bam_index
        References reference
        String base_file_name
        File monitoring_script
        String docker
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = ceil(2*size(matched_reads,"GB") + 4*size(assembly_bam,"GB")) + 10
    command <<< 
        set -o pipefail
        set -ex
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        samtools sort -N ~{assembly_bam} -o assembly_name_sorted.bam
        sv_rematch -v -S -M ~{matched_reads} \
            assembly_name_sorted.bam \
            ~{reference.ref_fasta} \
            ~{base_file_name}_assembly.support.name_sorted.bam

        samtools sort ~{base_file_name}_assembly.support.name_sorted.bam -o ~{base_file_name}_assembly.support.bam
        samtools index ~{base_file_name}_assembly.support.bam
    >>> 
    runtime {
        memory: "16 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_tries
        noAddress: no_address
    }   
    output {
        File monitoring_log = "monitoring.log"
        File updated_assembly = "~{base_file_name}_assembly.support.bam"
        File updated_assembly_index = "~{base_file_name}_assembly.support.bam.bai"
    }
}

task IdentifyVariants {
    input {
        Array[File]? input_crams
        Array[File]? input_crams_indexes
        Array[String] sample_names
        References references
        File assembly
        File assembly_index
        String output_vcf_prefix
        File? blacklist_bed
        String config_file_string
        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = ceil((if defined(input_crams) then size(select_first([input_crams]), "GB") else 0) + 2*size(assembly,"GB") + size(references.ref_fasta,"GB")) + 20
    Int mem = 16
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        python3 <<CODE
        import re
        text = "~{config_file_string}"
        pattern = r'(\w+(?:\.\w+)+\s*=\s*(?:\d+(?:\.\d+)?|true|false|\w+))'
        split_text = re.findall(pattern, text)

        # Join the split strings with a newline character
        output_text = "\n".join(split_text)

        with open("gridss.config", "w") as file:
            file.write(output_text)
        CODE

        java -Xmx~{mem-2}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
        SAMPLE_NAMES=~{sep=' SAMPLE_NAMES=' sample_names} \
        R=~{references.ref_fasta} \
        O=~{output_vcf_prefix}.vcf \
        ~{"BLACKLIST= " + blacklist_bed} \
        ASSEMBLY=~{assembly} \
        C=gridss.config

        bcftools view ~{output_vcf_prefix}.vcf -Oz -o ~{output_vcf_prefix}.vcf.gz
        bcftools index -t ~{output_vcf_prefix}.vcf.gz

    >>>
    runtime {
        memory: mem + " GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File identify_vcf = "~{output_vcf_prefix}.vcf.gz"
        File identify_vcf_index = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

task PreFilterCandidates {
    input {
        File input_vcf
        String filter_string
        String output_vcf_prefix
        String docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = ceil(2*size(input_vcf,"GB")) + 3
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &
        bcftools view -Oz -i "~{filter_string}" ~{input_vcf} -o ~{output_vcf_prefix}.vcf.gz
        bcftools index -t ~{output_vcf_prefix}.vcf.gz
    >>>
    runtime{
        memory: "4 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File filtered_vcf = "~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_index = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

task CollectGridssMetrics{
    input{
        File input_cram
        File input_cram_index
        References references
        String? interval

        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }

    Int disk_size = ceil(size(input_cram,"GB") + size(input_cram_index,"GB") + size(references.ref_fasta,"GB")) + 20
    String output_file = input_cram + ".gridss.working/" + basename(input_cram)
    String output_folder = input_cram + ".gridss.working/"
    Int mem = 8
    Int java_mem = mem - 2
    Int cpu = 2

    command{
set -xeo pipefail

bash ~{monitoring_script} | tee monitoring.log >&2 &
mkdir -p ~{input_cram}.gridss.working
java -Xmx~{java_mem}g -Xms~{java_mem}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.analysis.CollectGridssMetrics \
INPUT=~{input_cram} \
OUTPUT=~{output_file} \
THRESHOLD_COVERAGE=25000 \
FILE_EXTENSION=null \
GRIDSS_PROGRAM=CollectCigarMetrics \
GRIDSS_PROGRAM=CollectMapqMetrics \
GRIDSS_PROGRAM=CollectTagMetrics \
GRIDSS_PROGRAM=CollectIdsvMetrics \
GRIDSS_PROGRAM=ReportThresholdCoverage \
PROGRAM=null \
PROGRAM=CollectInsertSizeMetrics \
COMPRESSION_LEVEL=5 \
CREATE_INDEX=false \
CREATE_MD5_FILE=false \
GA4GH_CLIENT_SECRETS=client_secrets.json \
MAX_RECORDS_IN_RAM=500000 \
QUIET=false \
REFERENCE_SEQUENCE=~{references.ref_fasta} \
VALIDATION_STRINGENCY=STRICT \
VERBOSITY=INFO \
~{"INTERVALS=" + interval}
    }
    parameter_meta {
      input_cram: {
          localization_optional: true
      }
    }
    runtime {
        memory: mem + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        Array[File] out_metrics = [output_file + ".cigar_metrics", output_file + ".mapq_metrics",output_file + ".coverage.blacklist.bed", output_file + ".tag_metrics", output_file + ".idsv_metrics"]
    }
}

task AnnotateVariants {
    input {
        File input_vcf
        File input_vcf_index
        Array[File] input_germline_crams
        Array[File] input_germline_crams_indexes
        Array[File] input_tumor_crams
        Array[File] input_tumor_crams_indexes
        References references
        File assembly
        File assembly_index
        File? interval
        String output_vcf_prefix
        String config_file_string
        File? blacklist_bed
        Boolean is_somatic
        Int number_of_shards
        Array[File]? germline_metrics
        Array[File]? tumor_metrics
        Array[File]? assembly_metrics
        String gridss_docker
        String? cloud_provider_override
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
        Int? cpu_override
        Int? memory_override
    }
    Int disk_size = 3*ceil(
            size(input_germline_crams,"GB")/number_of_shards +
            size(input_tumor_crams, "GB")/number_of_shards +
            size(assembly,"GB")/number_of_shards + size(references.ref_fasta,"GB")) + 10
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false
    Boolean defined_germline = if(length(input_germline_crams)>0) then true else false
    Int cpu = select_first([ cpu_override, if(is_aws) then 2 else 1])
    Int mem = select_first([ memory_override, if(!is_somatic) then 8 else 128 ])
    # variables for metrics files reordering
    Boolean defined_germline_metrics = if(defined(germline_metrics)) then true else false
    Boolean defined_tumor_metrics = if(defined(tumor_metrics)) then true else false
    Boolean defined_assembly_metrics = if(defined(assembly_metrics)) then true else false

    command <<<
        set -o pipefail
        set -o xtrace
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        python3 <<CODE
        import re
        text = "~{config_file_string}"
        pattern = r'(\w+(?:\.\w+)+\s*=\s*(?:\d+(?:\.\d+)?|true|false|\w+))'
        split_text = re.findall(pattern, text)

        # Join the split strings with a newline character
        output_text = "\n".join(split_text)

        with open("gridss.config", "w") as file:
            file.write(output_text)
        CODE

        ### Order metrics files to be in the same order as the input cram files ###
        if ~{defined_germline_metrics}
        then
            input_germline_crams_string="~{sep=',' input_germline_crams}"
            first_cram=$(echo $input_germline_crams_string | awk -F "," '{print $1}')
            germline_working_folder=$first_cram".gridss.working/"
            mkdir -p $germline_working_folder
            for germline_metric_file in ~{sep=' ' germline_metrics}; do
                cp $germline_metric_file $germline_working_folder
            done
            ls -lrta $germline_working_folder
        fi

        if ~{defined_tumor_metrics}
        then
            input_tumor_crams_string="~{sep=',' input_tumor_crams}"
            first_cram=$(echo $input_tumor_crams_string | awk -F "," '{print $1}')
            tumor_working_folder=$first_cram".gridss.working/"
            mkdir -p $tumor_working_folder
            for tumor_metric_file in ~{sep=' ' tumor_metrics}; do
                cp $tumor_metric_file $tumor_working_folder
            done
            ls -lrta $tumor_working_folder
        fi

        if ~{defined_assembly_metrics}
        then
            assembly_working_folder="~{assembly}.gridss.working/"
            mkdir -p $assembly_working_folder
            for assembly_metric_file in ~{sep=' ' assembly_metrics}; do
                cp $assembly_metric_file $assembly_working_folder
            done
            ls -lrta $assembly_working_folder
        fi
        ### Done order metrics files ###

        if ~{is_aws}
        then
            if ~{is_somatic}
            then
                input_tumor_crams_string="I=~{sep=' I=' input_tumor_crams}"
            else
                input_tumor_crams_string=""
            fi
            if ~{defined_germline}
            then
                input_germline_crams_string="I=~{sep=' I=' input_germline_crams}"
            else
                input_germline_crams_string=""
            fi

            echo $input_tumor_crams_string
            echo $input_germline_crams_string


            java -Xmx~{mem-4}g -Xms~{mem/2}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \
                INPUT_VCF=~{input_vcf} \
                R=~{references.ref_fasta} \
                $input_tumor_crams_string \
                $input_germline_crams_string \
                ~{"BLACKLIST= " + blacklist_bed} \
                C=gridss.config \
                ASSEMBLY=~{assembly} \
                OUTPUT_VCF=~{output_vcf_prefix}.ann.vcf \
                THREADS=1

        else
            # convert interval list to bed file
            gatk  \
            IntervalListToBed -I ~{interval} -O interval.bed

            # take only the relevant part from the vcf
            bcftools view -R interval.bed ~{input_vcf} -o output_region.vcf
            input_vcf=output_region.vcf

            # download only the relevant part of the crams
            if ~{is_somatic}
            then
                gatk  \
                    PrintReads \
                -I ~{sep=' -I ' input_tumor_crams} \
                -O /dev/stdout \
                -L ~{interval} \
                -R ~{references.ref_fasta} | \
                samtools view -C -T ~{references.ref_fasta} -o input_tumor.cram --output-fmt-option embed_ref=1 -

                samtools index input_tumor.cram
            fi

            if ~{defined_germline}
            then
                gatk  \
                    PrintReads \
                -I ~{sep=' -I ' input_germline_crams} \
                -O /dev/stdout \
                -L ~{interval} \
                -R ~{references.ref_fasta} | \
                samtools view -C -T ~{references.ref_fasta} -o input_germline.cram --output-fmt-option embed_ref=1 -

                samtools index input_germline.cram
            fi

        gatk \
            PrintReads \
          -I ~{assembly} \
          -O /dev/stdout \
          -L ~{interval} \
          -R ~{references.ref_fasta} | \
            samtools view -C -T ~{references.ref_fasta} -o assembly_partial.cram --output-fmt-option embed_ref=1 -

        samtools index assembly_partial.cram

            java -Xmx~{mem-4}g -Xms~{mem/2}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \
            INPUT_VCF=$input_vcf \
            R=~{references.ref_fasta} \
            ~{if is_somatic then "I=input_tumor.cram" else ""} \
            ~{if defined_germline then "I=input_germline.cram" else ""} \
            ~{"BLACKLIST= " + blacklist_bed} \
            C=gridss.config \
            ASSEMBLY=assembly_partial.cram \
            OUTPUT_VCF=~{output_vcf_prefix}.ann.vcf \
            THREADS=1
        fi
    >>>
    parameter_meta {
      input_tumor_crams: {
          localization_optional: true
      }
      input_germline_crams: {
          localization_optional: true
      }
      assembly: {
          localization_optional: true
      }
    }
    runtime {
        memory: mem + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File annotated_vcf = "~{output_vcf_prefix}.ann.vcf"
        File annotated_vcf_index = "~{output_vcf_prefix}.ann.vcf.idx"
    }
}

task GermlineLinkVariants {
    input {
        File input_vcf
        File input_vcf_index
        References references
        String output_vcf_prefix
        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
        Int? memory_override
    }
    Int disk_size = 4*ceil(size(input_vcf,"GB") + size(references.ref_fasta,"GB")) + 10
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        Rscript /opt/gridss/link_breakpoints.R \
            --input ~{input_vcf} \
            --fulloutput ~{output_vcf_prefix}_linked.vcf \
            --ref ~{references.ref_fasta} \
            --scriptdir /opt/gridss/
    >>>
    runtime {
        memory: select_first([memory_override, 64]) + " GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File linked_vcf = "~{output_vcf_prefix}_linked.vcf.bgz"
        File linked_vcf_index = "~{output_vcf_prefix}_linked.vcf.bgz.tbi"
    }
}

task SomaticGripss {
    input {
        File input_vcf
        File input_vcf_index
        String tumor_sample
        String? germline_sample
        String reference_name
        String output_vcf_prefix
        References references
        File? pon_sgl_file
        File? pon_sv_file
        File? repeat_mask_file
        File? known_hotspot_file
        Int? min_normal_coverage
        String? exclude_filters
        String gripss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address

    }
    Int disk_size = ceil(2*size(input_vcf,"GB") + size(references.ref_fasta,"GB")) + 10
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir gripss_output
        java -XX:MaxRAMPercentage=80.0 -jar /opt/wtsi-cgp/java/gripss.jar \
        -vcf ~{input_vcf} \
        -sample ~{tumor_sample} \
        ~{"-reference " + germline_sample} \
        -ref_genome_version ~{reference_name} \
        -ref_genome ~{references.ref_fasta} \
        ~{"-pon_sgl_file " + pon_sgl_file} \
        ~{"-pon_sv_file " + pon_sv_file} \
        ~{"-known_hotspot_file " + known_hotspot_file} \
        ~{"-repeat_mask_file " + repeat_mask_file} \
        ~{"-min_normal_coverage " + min_normal_coverage} \
        ~{"-exclude_filters " + exclude_filters} \
        -output_dir gripss_output

        mv gripss_output/~{tumor_sample}.gripss.vcf.gz  ~{output_vcf_prefix}.gripss.vcf.gz
        mv gripss_output/~{tumor_sample}.gripss.vcf.gz.tbi  ~{output_vcf_prefix}.gripss.vcf.gz.tbi
    >>>

    runtime {
        memory: "64 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: gripss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File gripss_vcf = "~{output_vcf_prefix}.gripss.vcf.gz"
        File gripss_vcf_index = "~{output_vcf_prefix}.gripss.vcf.gz.tbi"
    }
}

task ConvertVcfFormat {
    input {
        File input_vcf
        File input_vcf_index
        String output_vcf_prefix
        References references
        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
        Int? memory_override
    }

    Int cpu = 8
    Int mem = select_first([memory_override, 64])

    Int disk_size = ceil(2*size(input_vcf,"GB") + size(references.ref_fasta,"GB")) + 10
    command <<<
        set -xeo pipefail
        
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        Rscript /opt/gridss/convert_vcf_format.R \
            --input_vcf ~{input_vcf} \
            --output_vcf ~{output_vcf_prefix}.vcf \
            --reference ~{references.ref_fasta} \
            --n_jobs ~{cpu}
        bcftools view -Oz -o ~{output_vcf_prefix}.vcf.gz ~{output_vcf_prefix}.vcf.bgz
        bcftools index -t ~{output_vcf_prefix}.vcf.gz

    >>>
    runtime {
        memory: mem + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File output_vcf = "~{output_vcf_prefix}.vcf.gz"
        File output_vcf_index = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}