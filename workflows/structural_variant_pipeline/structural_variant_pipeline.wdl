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
        String pipeline_version = "1.12.0" # !UnusedDeclaration

        String base_file_name
        Array[File] input_germline_crams = []
        Array[File] input_germline_crams_indexes = []
        Array[File] input_tumor_crams = []
        Array[File] input_tumor_crams_indexes = []
        References references
        UaReferences ua_references
        File wgs_calling_interval_list
        Int min_base
        Int min_mapq
        Int? max_num_haps
        Int realign_mapq
        String? min_indel_sc_size_to_include
        Int homopolymer_length
        String config_file_string
        File? blacklist_bed
        Boolean is_somatic
        String reference_name
        Boolean run_ua
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
        # Used for running on other clouds (aws)
        String? cloud_provider_override
        Int scatter_intervals_break # Maximal resolution for scattering intervals
        String dummy_input_for_call_caching = ""
        File? monitoring_script_input

        # Winval validations
        #@wv min_base >= 0
        #@wv min_mapq >= 0
        #@wv max_num_haps >= 0
        #@wv num_shards > 0
        #@wv realign_mapq >= 0
        #@wv scatter_intervals_break > 0
        #@wv reference_name in {"38","19"}
        # tumor + germline
        #@wv is_somatic -> defined(input_tumor_crams) and len(input_tumor_crams)>0
        #@wv is_somatic -> defined(input_tumor_crams_indexes) and len(input_tumor_crams_indexes)>0
        # germline only
        #@wv not is_somatic -> defined(input_germline_crams) and len(input_germline_crams)>0
        #@wv not is_somatic -> defined(input_germline_crams_indexes) and len(input_germline_crams_indexes)>0
        #@wv not is_somatic -> not defined(input_tumor_crams) or len(input_tumor_crams)==0
        #@wv not is_somatic -> not defined(input_tumor_crams_indexes) or len(input_tumor_crams_indexes)==0

    }
    meta {
        description : "Runs Structural variant pipeline\nThis pipeline supports germline and somatic modes\nThe input of that pipeline is cram files and the output is vcf file\nThe steps of the pipeline are as following:\n-Create an assembly file out of the cram files\n-Run UA alingnment on that\n-Fix the UA alignment which are secondarily mapped to decoy or with low mapq\n-Run gridss.IdentifyVariants and gridss.AnnotateVariants\n-Run R script / GRIPSS for filtering and linkage the variants"
        author: "Ultima Genomics"
        WDL_AID: { exclude: [
            "pipeline_version",
            "no_address",
            "preemptible_tries_override",
            "dummy_input_for_call_caching",
            "monitoring_script_input",
            "CreateAssembly.no_address",
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
            "Globals.glob",
            "Sentieon.Globals.glob",
            "AnnotateVCF.Globals.glob",
            "SingleSampleQC.Globals.glob",
            "VariantCallingEvaluation.Globals.glob",
            "aws_AnnotateVariants.interval",
            "AnnotateVariants.germline_metrics",
            "AnnotateVariants.tumor_metrics",
            "AnnotateVariants.assembly_metrics",
            "AnnotateVariants.germline_metrics_folder",
            "AnnotateVariants.tumor_metrics_folder",
            "AnnotateVariants.assembly_metrics_folder"
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
        help: "Input CRAM file for the germline or matched normal sample; optinal for supporting somatic calling tumor only",
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
        type: "UaReferences",
        help: "UAReference files: ua_index, ref_alt, v_aware_alignment_flag and ua_extra_args, recommended value set in the template",
        category: "required"
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
            help: "Assembly parameter: Minimum size of an indel and soft-clipping in the read to include the read in the assembly.",
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
        output_vcf: {
            type: "File",
            help: "Final VCF file",
            category: "output"
        }
        output_vcf_index: {
            type: "File",
            help: "Final VCF index file",
            category: "output"
        }
        assembly: {
            type: "File",
            help: "Assembly output before UA realingment",
            category: "output"
        }
        assembly_index: {
            type: "File",
            help: "Assembly output index before UA realingment",
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
    }
    Int preemptibles = select_first([preemptible_tries_override, 1])

    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

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
            docker = global.gitc_docker,
            gitc_path = global.gitc_jar_path,
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
                assembly_docker = global.assembly_docker,
                references = references,
                interval = interval,
                min_base = min_base,
                min_mapq = min_mapq,
                max_num_haps = max_num_haps,
                min_sc_indel_size = min_indel_sc_size_to_include,
                is_somatic = is_somatic,
                number_of_shards = ScatterIntervalList.interval_count,
                cloud_provider_override = cloud_provider_override,
                no_address = no_address,
                preemptible_tries = preemptibles
        }
    }

    call UGGeneralTasks.MergeBams {
        input:
            inputs = CreateAssembly.assembly,
            output_prefix = base_file_name + '_merged_assembly',
            docker = global.gitc_docker,
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

        call RevertLowMAPQUASecondaryAlignment {
            input:
                before_UA = MergeBams.output_bam,
                before_UA_index = MergeBams.output_bam_index,
                after_UA = AlignWithUA.ua_output_bam,
                outptut_prefix = base_file_name +"_assembly_ua_realigned",
                realign_mapq = realign_mapq,
                gridss_docker = global.gridss_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }

        call LongHomopolymersAlignmnet {
            input:
                input_bam = RevertLowMAPQUASecondaryAlignment.realigned_assembly,
                input_bam_index = RevertLowMAPQUASecondaryAlignment.realigned_assembly_index,
                output_bam_basename = base_file_name +"_assembly_ua_realigned_long_homopolymers_aligned",
                references = references,
                homopolymer_length = homopolymer_length,
                gridss_docker = global.gridss_docker,
                monitoring_script = monitoring_script,
                preemptible_tries = preemptibles,
                no_address = no_address

        }
    }
    File assembly_file = select_first([LongHomopolymersAlignmnet.realigned_assembly, MergeBams.output_bam])
    File assembly_file_index = select_first([LongHomopolymersAlignmnet.realigned_assembly_index, MergeBams.output_bam_index])

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
                germline_metrics_folder = germline_CollectGridssMetrics.out_metrics_folder,
                tumor_metrics = tumor_CollectGridssMetrics.out_metrics,
                tumor_metrics_folder = tumor_CollectGridssMetrics.out_metrics_folder,
                assembly_metrics = assembly_CollectGridssMetrics.out_metrics,
                assembly_metrics_folder = assembly_CollectGridssMetrics.out_metrics_folder,
                gridss_docker = global.gridss_docker,
                cloud_provider_override = cloud_provider_override,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
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
                no_address = no_address
            }
    }

        call UGGeneralTasks.MergeVCFs as MergeVCFs {
            input:
            input_vcfs = AnnotateVariants.annotated_vcf,
            input_vcfs_indexes = AnnotateVariants.annotated_vcf,
            output_vcf_name = base_file_name + ".ann.vcf.gz",
            preemptible_tries = preemptibles,
            monitoring_script = monitoring_script,
            docker = global.gitc_docker,
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
                reference_name = reference_name,
                output_vcf_prefix = base_file_name,
                gridss_docker = global.gridss_docker,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
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
                reference_name = reference_name,
                preemptible_tries = preemptibles,
                monitoring_script = monitoring_script,
                no_address = no_address
        }
    }
    output {
        File output_vcf = select_first([GermlineLinkVariants.linked_vcf, SomaticGripss.gripss_vcf])
        File output_vcf_index = select_first([GermlineLinkVariants.linked_vcf_index, SomaticGripss.gripss_vcf_index])
        File assembly = MergeBams.output_bam
        File assembly_index = MergeBams.output_bam_index
        File? realigned_assembly = LongHomopolymersAlignmnet.realigned_assembly
        File? realigned_assembly_index = LongHomopolymersAlignmnet.realigned_assembly_index
        File? converted_vcf = ConvertVcfFormat.output_vcf
        File? converted_vcf_index = ConvertVcfFormat.output_vcf_index
    }
}

task CreateAssembly {
    input {
        Array[File]? input_germline_crams
        Array[File]? input_germline_crams_indexes
        Array[File]? input_tumor_crams
        Array[File]? input_tumor_crams_indexes
        Boolean is_somatic
        String output_prefix
        String assembly_docker
        References references
        File interval
        Int min_base
        Int min_mapq
        String? min_sc_indel_size
        Int? max_num_haps
        Int number_of_shards
        String? cloud_provider_override
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = 3*ceil((if defined(input_germline_crams) then size(select_first([input_germline_crams]),"GB")/number_of_shards else 0) +
                    (if is_somatic then size(select_first([input_tumor_crams]), "GB")/number_of_shards else 0)  +
                       size(references.ref_fasta,"GB")) + 10
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false
    Boolean defined_germline = if(defined(input_germline_crams)) then true else false
    
    command <<<
        set -o pipefail
        set -e
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
            ~{"--min-feature-length '" + min_sc_indel_size + "'"} \
            ~{"--max-num-haps " + max_num_haps} \
            --prog \
            --interval-nreads 10000 \
            --sv
        
        else        
            if ~{defined(input_germline_crams)}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' input_germline_crams} \
                -O input_germline.cram \
                -L ~{interval} \
                -R ~{references.ref_fasta}

                samtools index input_germline.cram
            fi

            if ~{is_somatic}
            then
                gatk --java-options "-Xms2G" PrintReads \
                -I ~{sep=' -I ' input_tumor_crams} \
                -O input_tumor.cram \
                -L ~{interval} \
                -R ~{references.ref_fasta}

                samtools index input_tumor.cram
            fi

            tool \
            ~{if is_somatic then (if defined(input_germline_crams) then "--input input_tumor.cram\\;input_germline.cram" else "--input input_tumor.cram") else "--input input_germline.cram"} \
            ~{if is_somatic then (if defined(input_germline_crams) then "--cram-index input_tumor.cram.crai\\;input_germline.cram.crai" else "--cram-index input_tumor.cram.crai") else "--cram-index input_germline.cram.crai"} \
            ~{if is_somatic then "--somatic" else ""} \
            --output ~{output_prefix} \
            --ref ~{references.ref_fasta} \
            --bed interval.bed \
            --min-base ~{min_base} \
            --min-mapq ~{min_mapq} \
            ~{"--min-feature-length '" + min_sc_indel_size + "'"} \
            ~{"--max-num-haps " + max_num_haps} \
            --prog \
            --interval-nreads 10000 \
            --sv

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
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        docker: assembly_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File assembly = "~{output_prefix}_hap_out_sorted.bam"
        File assembly_index = "~{output_prefix}_hap_out_sorted.bam.bai"
    }
}


task RevertLowMAPQUASecondaryAlignment {
    input {
        File before_UA
        File before_UA_index
        File after_UA
        Int realign_mapq
        String outptut_prefix
        String gridss_docker
        Int preemptible_tries
        File monitoring_script
        Boolean no_address
    }
    Int disk_size = ceil(size(before_UA,"GB") + 3*size(after_UA,"GB"))+20
    command <<<
        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        samtools sort ~{after_UA} -o ~{after_UA}_sorted.bam
        samtools index ~{after_UA}_sorted.bam
        # just to save disk space
        rm ~{after_UA}

        python3 /opt/gridss/revert_sup_low_mapq_ua_alignment.py \
            --before ~{before_UA} \
            --after ~{after_UA}_sorted.bam \
            --output ~{outptut_prefix}_unsorted.bam \
            --min_mapping_quality ~{realign_mapq}
        samtools sort -o ~{outptut_prefix}.bam ~{outptut_prefix}_unsorted.bam
        samtools index ~{outptut_prefix}.bam
    >>>
    runtime {
        memory: "4 GB"
        cpu: 1
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

task LongHomopolymersAlignmnet {
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
            set -o pipefail
            set -e
            bash ~{monitoring_script} | tee monitoring.log >&2 &

            python3 /opt/gridss/align_long_homopolymers.py \
                --input ~{input_bam} \
                --output ~{output_bam_basename}.bam \
                --reference ~{references.ref_fasta} \
                --homopolymer_length ~{homopolymer_length}

        >>>
        runtime {
            memory: "8 GB"
            cpu: 10
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
    command <<<
        set -o pipefail
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

        java -Xmx10g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
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
        memory: "32 GB"
        cpu: 1
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
        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        bcftools view -Oz -i "~{filter_string}" ~{input_vcf} -o ~{output_vcf_prefix}.vcf.gz
        bcftools index -t ~{output_vcf_prefix}.vcf.gz
    >>>
    runtime{
        memory: "4 GB"
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
    Int mem = 16
    Int java_mem = 14
    Int cpu = 4
    
    command{
set -o pipefail
set -e
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
        memory: mem + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        Array[File] out_metrics = [output_file + ".cigar_metrics", output_file + ".mapq_metrics",output_file + ".coverage.blacklist.bed", output_file + ".tag_metrics", output_file + ".idsv_metrics"]
        String out_metrics_folder = "~{output_folder}"
    }
}

task AnnotateVariants {
    input {
        File input_vcf
        File input_vcf_index
        Array[File]? input_germline_crams
        Array[File]? input_germline_crams_indexes
        Array[File]? input_tumor_crams
        Array[File]? input_tumor_crams_indexes
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
        String? germline_metrics_folder
        Array[File]? tumor_metrics
        String? tumor_metrics_folder
        Array[File]? assembly_metrics
        String? assembly_metrics_folder
        String gridss_docker
        String? cloud_provider_override
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = 3*ceil(
            (if defined(input_germline_crams) then size(select_first([input_germline_crams]),"GB")/number_of_shards else 0) +
            (if is_somatic then size(select_first([input_tumor_crams]), "GB")/number_of_shards else 0) +
            size(assembly,"GB")/number_of_shards + size(references.ref_fasta,"GB")) + 10
    Boolean is_aws = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then true else false
    Boolean defined_germline = if(defined(input_germline_crams)) then true else false
    Int cpu = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then 16 else 1
    Int mem = if(defined(cloud_provider_override) && select_first([cloud_provider_override]) == "aws") then 100 else 32
    Int aws_gridss_cpu = 14
    
    #variables for metrics files reordering
    Boolean defined_germline_metrics = if(defined(germline_metrics_folder)) then true else false
    Boolean defined_tumor_metrics = if(defined(tumor_metrics_folder)) then true else false
    Boolean defined_assembly_metrics = if(defined(assembly_metrics_folder)) then true else false
            
    command <<<
        set -o pipefail
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
            cp -r ~{germline_metrics_folder}* $germline_working_folder
            ls -lrta $germline_working_folder
        fi

        if ~{defined_tumor_metrics}
        then
            input_tumor_crams_string="~{sep=',' input_tumor_crams}"
            first_cram=$(echo $input_tumor_crams_string | awk -F "," '{print $1}')
            tumor_working_folder=$first_cram".gridss.working/"
            mkdir -p $tumor_working_folder
            cp -r ~{tumor_metrics_folder}* $tumor_working_folder
            ls -lrta $tumor_working_folder
        fi

        if ~{defined_assembly_metrics}
        then
            assembly_working_folder="~{assembly}.gridss.working/"
            mkdir -p $assembly_working_folder
            cp -r ~{assembly_metrics_folder}* $assembly_working_folder
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


            java -Xmx~{mem}g -Xms~{mem}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \
                INPUT_VCF=~{input_vcf} \
                R=~{references.ref_fasta} \
                $input_tumor_crams_string \
                $input_germline_crams_string \
                ~{"BLACKLIST= " + blacklist_bed} \
                C=gridss.config \
                ASSEMBLY=~{assembly} \
                OUTPUT_VCF=~{output_vcf_prefix}.ann.vcf \
                THREADS=~{aws_gridss_cpu}
            
        else
            # convert interval list to bed file
            java -Xms4g -jar /opt/gatk/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar  \
            IntervalListToBed -I ~{interval} -O interval.bed

            # take only the relevant part from the vcf
            bcftools view -R interval.bed ~{input_vcf} -o output_region.vcf
            input_vcf=output_region.vcf

            # download only the relevant part of the crams
            if ~{is_somatic}
            then
                java -Xms4g -jar /opt/gatk/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar  \
                    PrintReads \
                -I ~{sep=' -I ' input_tumor_crams} \
                -O input_tumor.cram \
                -L ~{interval} \
                -R ~{references.ref_fasta}

                samtools index input_tumor.cram
            fi

            if ~{defined(input_germline_crams)}
            then
                java -Xms4g -jar /opt/gatk/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar  \
                    PrintReads \
                -I ~{sep=' -I ' input_germline_crams} \
                -O input_germline.cram \
                -L ~{interval} \
                -R ~{references.ref_fasta}

                samtools index input_germline.cram
            fi

        java -Xms4g -jar /opt/gatk/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar  \
            PrintReads \
          -I ~{assembly} \
          -O assembly_partial.cram \
          -L ~{interval} \
          -R ~{references.ref_fasta}

        samtools index assembly_partial.cram

            java -Xmx~{mem}g -Xms~{mem}g -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \
            INPUT_VCF=$input_vcf \
            R=~{references.ref_fasta} \
            ~{if is_somatic then "I=input_tumor.cram" else ""} \
            ~{if defined(input_germline_crams) then "I=input_germline.cram" else ""} \
            ~{"BLACKLIST= " + blacklist_bed} \
            C=gridss.config \
            ASSEMBLY=assembly_partial.cram \
            OUTPUT_VCF=~{output_vcf_prefix}.ann.vcf \
            THREADS=~{cpu}
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
        memory: mem + " GB"
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
        String reference_name
        String output_vcf_prefix
        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = 4*ceil(size(input_vcf,"GB")) + 10
    command <<<
        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        Rscript /opt/gridss/link_breakpoints.R \
            --input ~{input_vcf} \
            --fulloutput ~{output_vcf_prefix}_linked.vcf \
            --ref BSgenome.Hsapiens.UCSC.hg~{reference_name} \
            --scriptdir /opt/gridss/
    >>>
    runtime {
        memory: "64 GB"
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

        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir gripss_output
        java -jar /opt/wtsi-cgp/java/gripss.jar \
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
        memory: "64 GB"
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
        String reference_name
        String gridss_docker
        File monitoring_script
        Int preemptible_tries
        Boolean no_address
    }
    Int disk_size = ceil(2*size(input_vcf,"GB")) + 10
    command <<<
        set -o pipefail
        set -e
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        Rscript /opt/gridss/convert_vcf_format.R \
            --input_vcf ~{input_vcf} \
            --output_vcf ~{output_vcf_prefix}.vcf \
            --reference BSgenome.Hsapiens.UCSC.hg~{reference_name} \
            --n_jobs 8

    >>>
    runtime {
        memory: "64 GB"
        cpu: 8
        disks: "local-disk " + disk_size + " HDD"
        docker: gridss_docker
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File monitoring_log = "monitoring.log"
        File output_vcf = "~{output_vcf_prefix}.vcf.bgz"
        File output_vcf_index = "~{output_vcf_prefix}.vcf.bgz.tbi"
    }
}


