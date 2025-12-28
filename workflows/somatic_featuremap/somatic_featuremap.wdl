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
# 
# given tumor and normal pileup-featuremap vcf files, the pipeline creates a merged somatic-pileup-featurmap vcf file.
# the merged file includes all original information per tumor/normal sample.
# additional information is added to the merged somatic-pileup-featurmap vcf :
# 1. tandem repeat information: proximity of variant to closest tandem repeat and its details.
# 2. ref/non-ref counts: in positions around the variant loci for tumor and normal samples.
# Eventually, given the whole tumor-normal information a pre-trained xgb model is used to infer probability of a candidate to be a true variant.
# The pipeline outputs the merged somatic-pileup-featurmap vcf including the added information and the corresponding probability.

# CHANGELOG in reverse chronological order

import "tasks/globals.wdl" as Globals
import "tasks/mrd.wdl" as UGMrdTasks
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/pileup_tasks.wdl" as UGPileupTasks

workflow SomaticFeaturemap {

    input {
        String pipeline_version = "1.26.0" # !UnusedDeclaration

        String base_file_name
        File? tumor_featuremap_vcf
        File? tumor_featuremap_vcf_index
        File? normal_featuremap_vcf
        File? normal_featuremap_vcf_index
        File ref_tandem_repeats

        File? somatic_featuremap_given
        File? somatic_featuremap_given_index
        
        Boolean? keep_non_pass_tumor_candidates_override
        Array[File]+? input_tumor_cram_bam
        Array[File]+? input_tumor_cram_bam_index
        Array[File]+? input_normal_cram_bam
        Array[File]+? input_normal_cram_bam_index
        Boolean? run_mpileup_override

        File? sfm_xgb_model

        Boolean? no_address_override
        Int? preemptible_tries_override
        File? monitoring_script_input

        Int pad_size=2
        Int min_mapq_for_samtools_mpileup=0
        File? interval_list
        Int? num_shards
        Int? scatter_intervals_break=100
        String dummy_input_for_call_caching = ""

        References? references

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv defined(tumor_featuremap_vcf) -> prefix(tumor_featuremap_vcf_index) == tumor_featuremap_vcf
        #@wv defined(normal_featuremap_vcf) -> prefix(normal_featuremap_vcf_index) == normal_featuremap_vcf
        #@wv defined(tumor_featuremap_vcf_index) -> suffix(tumor_featuremap_vcf_index) in {".tbi", ".csi"}
        #@wv defined(normal_featuremap_vcf_index) -> suffix(normal_featuremap_vcf_index) in {".tbi", ".csi"}

    }

    meta {
        description: "The SomaticFeaturemap workflow generates a merged somatic featuremap VCF file from tumor and normal featuremap VCF files. This workflow is designed to combine information from tumor-normal pairs and enrich the data with additional features for somatic variant analysis."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "no_address_override",
                "preemptible_tries_override",
                "Globals.glob",
                'PadSomaticFeaturemapVcf.disk_size',
                'PadSomaticFeaturemapVcf.memory_gb',
                'PadSomaticFeaturemapVcf.cpus',
                'tumor_createMpileup_scatter.max_depth',
                'tumor_createMpileup_scatter.min_BaseQ',
                'tumor_createMpileup_scatter.snp_file_index',
                'tumor_createMpileup_scatter.cloud_provider',
                'normal_createMpileup_scatter.max_depth',
                'normal_createMpileup_scatter.min_BaseQ',
                'normal_createMpileup_scatter.snp_file_index',
                'normal_createMpileup_scatter.cloud_provider',
                'ConcatScoredSfmVcfs.disk_size',
                'ConcatScoredSfmVcfs.preemptible_tries'
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        tumor_featuremap_vcf: {
            help: "Input tumor featuremap VCF file.",
            type: "File",
            category: "input_required"
        }
        tumor_featuremap_vcf_index: {
            help: "Input tumor featuremap VCF index file.",
            type: "File",
            category: "input_required"
        }
        normal_featuremap_vcf: {
            help: "Input normal featuremap VCF file.",
            type: "File",
            category: "input_required"
        }
        normal_featuremap_vcf_index: {
            help: "Input normal featuremap VCF index file.",
            type: "File",
            category: "input_required"
        }
        somatic_featuremap_given: {
            help: "If provided, use this somatic featuremap file instead of generating a new one.",
            type: "File?",
            category: "input_optional"
        }
        somatic_featuremap_given_index: {
            help: "If provided, use this somatic featuremap index file instead of generating a new one.",
            type: "File?",
            category: "input_optional"
        }
        ref_tandem_repeats: {
            help: "Reference tandem repeats file.",
            type: "File",
            category: "input_required"
        }
        keep_non_pass_tumor_candidates_override: {
            help: "If true, keep non-PASS tumor candidates in the somatic featuremap.",
            type: "Boolean",
            default: false,
            category: "input_optional"
        }
        no_address_override: {
            help: "If true, do not use the noAddress runtime attribute.",
            type: "Boolean",
            default: true,
            category: "input_optional"
        }
        preemptible_tries_override: {
            help: "Number of preemptible retries.",
            type: "Int",
            default: 1,
            category: "input_optional"
        }
        monitoring_script_input: {
            help: "Path to the monitoring script.",
            type: "File",
            category: "input_optional"
        }
        pad_size: {
            help: "Size of padding to add to the positions of variants defined as tumor-PASS in somatic-pileup-featuremap.",
            type: "Int",
            default: 2,
            category: "input_optional"
        }
        min_mapq_for_samtools_mpileup: {
            help: "Minimum mapping quality for samtools mpileup.",
            type: "Int",
            default: 0,
            category: "input_optional"
        }
        interval_list: {
            help: "Input interval list file.",
            type: "File",
            category: "input_required"
        }
        num_shards: {
            help: "Number of shards for splitting the interval list.",
            type: "Int",
            category: "input_optional"
        }
        scatter_intervals_break: {
            help: "Break scatter intervals at multiples of this value.",
            type: "Int",
            default: 100,
            category: "input_optional"
        }
        dummy_input_for_call_caching: {
            help: "Dummy input for call caching.",
            type: "File",
            category: "input_optional"
        }
        input_tumor_cram_bam: {
            help: "Input tumor CRAM/BAM file.",
            type: "File",
            category: "input_required"
        }
        input_tumor_cram_bam_index: {
            help: "Input tumor CRAM/BAM index file.",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam: {
            help: "Input normal CRAM/BAM file.",
            type: "File",
            category: "input_required"
        }
        input_normal_cram_bam_index: {
            help: "Input normal CRAM/BAM index file.",
            type: "File",
            category: "input_required"
        }
        run_mpileup_override: {
            help: "If true, run the mpileup step.",
            type: "Boolean",
            default: true,
            category: "input_optional"
        }
        sfm_xgb_model: {
            help: "XGBoost model file for somatic featuremap pileup scoring.",
            type: "File?",
            category: "input_optional"
        }
        references: {
            help: "Input reference genome file, index file and dict.",
            type: "Struct",
            category: "input_required"
        }
        somatic_featuremap: {
            help: "Output somatic featuremap pileup file.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_index: {
            help: "Output somatic featuremap pileup index file.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_tumor_pass: {
            help: "Output somatic featuremap pileup file with tumor-only PASS variants.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_tumor_pass_index: {
            help: "Output somatic featuremap pileup index file with tumor-only PASS variants.",
            type: "File",
            category: "output"
        }
        tumor_mpileup_vcf: {
            help: "Output tumor mpileup VCF file.",
            type: "File",
            category: "output"
        }
        tumor_mpileup_vcf_index: {
            help: "Output tumor mpileup VCF index file.",
            type: "File",
            category: "output"
        }
        normal_mpileup_vcf: {
            help: "Output normal mpileup VCF file.",
            type: "File",
            category: "output"
        }
        normal_mpileup_vcf_index: {
            help: "Output normal mpileup VCF index file.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_mpileup_vcf: {
            help: "Output somatic featuremap pileup VCF with integrated mpileup info.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_mpileup_vcf_index: {
            help: "somatic_featuremap_mpileup_vcf index file.",
            type: "File",
            category: "output"
        }
        tumor_mpileup: {
            help: "Output tumor mpileup file.",
            type: "File",
            category: "output"
        }
        normal_mpileup: {
            help: "Output normal mpileup file.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_tumor_pass_mpileup_xgb_vcf: {
            help: "Output somatic featuremap pileup VCF tumor-only PASS variants with integrated mpileup info and XGB scores.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_tumor_pass_mpileup_xgb_vcf_index: {
            help: "somatic_featuremap_tumor_pass_mpileup_xgb_vcf index file.",
            type: "File",
            category: "output"
        }
        final_out_sfm_vcf: {
            help: "Final output somatic featuremap pileup VCF file.",
            type: "File",
            category: "output"
        }
        final_out_sfm_vcf_index: {
            help: "Final output somatic featuremap pileup VCF index file.",
            type: "File",
            category: "output"
        }
        somatic_featuremap_padded_bed_file: {
            help: "Output padded BED file for somatic featuremap variants.",
            type: "File",
            category: "output"
        }
    }

    Boolean keep_non_pass_tumor_candidates = select_first([keep_non_pass_tumor_candidates_override, false])
    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])
    Boolean run_mpileup = select_first([run_mpileup_override, true])

    call Globals.Globals as Globals
      GlobalVariables global = Globals.global_dockers

    if (!defined(somatic_featuremap_given)){
        call CreateSomaticFeaturemap {
            input:
            sample_name = base_file_name,
            tumor_featuremap_vcf = select_first([tumor_featuremap_vcf]),
            tumor_featuremap_vcf_index = select_first([tumor_featuremap_vcf_index]),
            normal_featuremap_vcf = select_first([normal_featuremap_vcf]),
            normal_featuremap_vcf_index = select_first([normal_featuremap_vcf_index]),
            keep_non_pass_tumor_candidates = keep_non_pass_tumor_candidates,
            docker = global.ugbio_featuremap_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
        }
    }
    
    File somatic_featuremap_vcf = select_first([somatic_featuremap_given, CreateSomaticFeaturemap.somatic_featuremap])
    File somatic_featuremap_vcf_index = select_first([somatic_featuremap_given_index, CreateSomaticFeaturemap.somatic_featuremap_index])
    References references_variable = select_first([references])
    
    if (run_mpileup)
    {    
        # 2. Make a padded regions file for tumor-pass variants
        call UGMrdTasks.PadVcf as PadSomaticFeaturemapVcf {
            input:
            input_vcf = somatic_featuremap_vcf,
            ref_fai = select_first([references_variable.ref_fasta_index]),
            pad_size = pad_size,
            docker = global.ugbio_core_docker,
            preemptible_tries = preemptible_tries,
            monitoring_script = monitoring_script
        }
        File somatic_featuremap_padded_bed = PadSomaticFeaturemapVcf.padded_bed
        # 3. Run samtools mpileup: Tumor and Normal separately
        call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList{
            input:
                interval_list = select_first([interval_list]),
                scatter_count = select_first([num_shards]),
                break_bands_at_multiples_of = select_first([scatter_intervals_break]),
                dummy_input_for_call_caching = select_first([dummy_input_for_call_caching]),
                docker = global.broad_gatk_docker,
                no_address = true,
                monitoring_script = monitoring_script
        }

        scatter (interval in ScatterIntervalList.out){
                call UGPileupTasks.CreateMpileup as tumor_createMpileup_scatter {
                input:
                    input_bam_files = select_first([input_tumor_cram_bam]),
                    input_bam_files_index = select_first([input_tumor_cram_bam_index]),
                    reference_fasta = references_variable.ref_fasta,
                    reference_fai = references_variable.ref_fasta_index,
                    reference_dict = references_variable.ref_dict,
                    min_MapQ = min_mapq_for_samtools_mpileup,
                    snp_file = somatic_featuremap_padded_bed,
                    interval = interval,
                    docker = global.ugbio_freec_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script
                }

                call UGPileupTasks.CreateMpileup as normal_createMpileup_scatter{
                input:
                    input_bam_files = select_first([input_normal_cram_bam]),
                    input_bam_files_index = select_first([input_normal_cram_bam_index]),
                    reference_fasta = select_first([references_variable.ref_fasta]),
                    reference_fai = select_first([references_variable.ref_fasta_index]),
                    reference_dict = select_first([references_variable.ref_dict]),
                    min_MapQ = select_first([min_mapq_for_samtools_mpileup]),
                    snp_file = somatic_featuremap_padded_bed,
                    interval = interval,
                    docker = global.ugbio_freec_docker,
                    no_address = no_address,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script
                }
            }

        Array[File] tumor_pileup_files = tumor_createMpileup_scatter.out_pileup
        Array[File] normal_pileup_files = normal_createMpileup_scatter.out_pileup


        call UGGeneralTasks.ConcatFiles as tumor_ConcatMpileupFiles{
            input:
                files = tumor_pileup_files,
                out_file_name = basename(select_first([input_tumor_cram_bam])[0])+"_minipileup.pileup",
                docker = global.ubuntu_docker
        }
        call UGGeneralTasks.ConcatFiles as normal_ConcatMpileupFiles{
            input:
                files = normal_pileup_files,
                out_file_name = basename(select_first([input_normal_cram_bam])[0])+"_minipileup.pileup",
                docker = global.ubuntu_docker
        }
        File tumor_pileup = tumor_ConcatMpileupFiles.out_merged_file
        File normal_pileup = normal_ConcatMpileupFiles.out_merged_file

        call MpileupIntegrationToSFM {
            input:
                somatic_featuremap_vcf = somatic_featuremap_vcf,
                somatic_featuremap_vcf_index = somatic_featuremap_vcf_index,
                tumor_mpileup = tumor_pileup,
                normal_mpileup = normal_pileup,
                distance_start_to_center = pad_size,
                monitoring_script = monitoring_script,
                no_address = no_address,
                preemptible_tries = preemptible_tries,
                docker = global.ugbio_featuremap_docker
        }
    }
    File somatic_featuremap_tumor_pass_mpileup_vcf = select_first([somatic_featuremap_given, MpileupIntegrationToSFM.out_sfm_vcf])
    File somatic_featuremap_tumor_pass_mpileup_vcf_index = select_first([somatic_featuremap_given_index, MpileupIntegrationToSFM.out_sfm_vcf_index])
        
    
    call UGGeneralTasks.ScatterIntervalList as ScoreSomaticVariants_ScatterIntervalList{
            input:
                interval_list = select_first([interval_list]),
                scatter_count = select_first([num_shards]),
                break_bands_at_multiples_of = select_first([scatter_intervals_break]),
                dummy_input_for_call_caching = select_first([dummy_input_for_call_caching]),
                docker = global.broad_gatk_docker,
                no_address = true,
                monitoring_script = monitoring_script
    }
    scatter (interval in ScoreSomaticVariants_ScatterIntervalList.out)
    {
        call ScoreSomaticVariants {
            input:
                input_sfm_vcf = somatic_featuremap_tumor_pass_mpileup_vcf,
                input_sfm_vcf_index = somatic_featuremap_tumor_pass_mpileup_vcf_index,
                xgb_model = sfm_xgb_model,
                interval_list = interval,
                ref_tandem_repeats = ref_tandem_repeats,
                reference_fai = references_variable.ref_fasta_index,
                docker = global.ugbio_featuremap_docker,
                no_address = no_address,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script
        }
    }

    Array[File] scored_sfm_vcfs = ScoreSomaticVariants.out_sfm_vcf
    Array[File] scored_sfm_vcfs_indexes = ScoreSomaticVariants.out_sfm_vcf_index

    ### also maybe add ConcatParquetFiles task here for future use
    call UGGeneralTasks.ConcatVcfs as ConcatScoredSfmVcfs {
        input:
            input_vcfs = scored_sfm_vcfs,
            input_vcfs_indexes = scored_sfm_vcfs_indexes,
            output_vcf_name = base_file_name + ".somatic_featuremap.tr_info.xgb_proba.vcf.gz",
            docker = global.broad_gatk_docker,
            preemptible_tries = preemptible_tries,
            no_address = no_address,
            monitoring_script = monitoring_script
    }
    
    # optional outputs
    # if somatic_featuremap_given is provided, then the outputs from CreateSomaticFeaturemap are not available
    # if run_mpileup is false, then the outputs from MpileupIntegrationToSFM are not available

    if (!defined(somatic_featuremap_given))
    {
        File out_somatic_featuremap = select_first([CreateSomaticFeaturemap.somatic_featuremap])
        File out_somatic_featuremap_index = select_first([CreateSomaticFeaturemap.somatic_featuremap_index])
    }
    if (run_mpileup)
    {
        File out_tumor_mpileup = select_first([tumor_ConcatMpileupFiles.out_merged_file])
        File out_normal_mpileup = select_first([normal_ConcatMpileupFiles.out_merged_file])
        File out_somatic_featuremap_mpileup_vcf = select_first([MpileupIntegrationToSFM.out_sfm_vcf])
        File out_somatic_featuremap_mpileup_vcf_index = select_first([MpileupIntegrationToSFM.out_sfm_vcf_index])
        File out_somatic_featuremap_padded_bed = select_first([PadSomaticFeaturemapVcf.padded_bed])
    }
    
    File final_sfm_vcf = select_first([ConcatScoredSfmVcfs.output_vcf,MpileupIntegrationToSFM.out_sfm_vcf,CreateSomaticFeaturemap.somatic_featuremap])
    File final_sfm_vcf_index = select_first([ConcatScoredSfmVcfs.output_vcf_index,MpileupIntegrationToSFM.out_sfm_vcf_index,CreateSomaticFeaturemap.somatic_featuremap_index])

    output {
        File? somatic_featuremap = out_somatic_featuremap
        File? somatic_featuremap_index = out_somatic_featuremap_index
        File? tumor_mpileup = out_tumor_mpileup
        File? normal_mpileup = out_normal_mpileup
        File? somatic_featuremap_mpileup_vcf = out_somatic_featuremap_mpileup_vcf
        File? somatic_featuremap_mpileup_vcf_index = out_somatic_featuremap_mpileup_vcf_index
        File final_out_sfm_vcf = final_sfm_vcf
        File final_out_sfm_vcf_index = final_sfm_vcf_index
        File? somatic_featuremap_padded_bed_file = out_somatic_featuremap_padded_bed
    }
}

task CreateSomaticFeaturemap {
    input {
        String sample_name
        File tumor_featuremap_vcf
        File tumor_featuremap_vcf_index
        File normal_featuremap_vcf
        File normal_featuremap_vcf_index
        Boolean keep_non_pass_tumor_candidates
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float tumor_featuremap_vcf_size = size(tumor_featuremap_vcf, "GB")
    Float normal_featuremap_vcf_size = size(normal_featuremap_vcf, "GB")
    Float additional_disk = 100
    Int disk_size = ceil(2 * tumor_featuremap_vcf_size + normal_featuremap_vcf_size + additional_disk)
    
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        
        create_somatic_featuremap \
            --tumor_vcf  ~{tumor_featuremap_vcf} \
            --normal_vcf ~{normal_featuremap_vcf} \
            --sample_name ~{sample_name} \
            --out_directory . \
            ~{true="--keep-non-pass-tumor-candidates" false='' keep_non_pass_tumor_candidates}

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 4
    }
    output {
        File somatic_featuremap = "~{sample_name}.tumor_normal.merged.vcf.gz"
        File somatic_featuremap_index = "~{sample_name}.tumor_normal.merged.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}


task MpileupIntegrationToSFM {
    input{
        File somatic_featuremap_vcf
        File somatic_featuremap_vcf_index
        File tumor_mpileup
        File normal_mpileup
        Int distance_start_to_center
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
        String docker
    }
    String out_sfm_vcf_filename = basename(somatic_featuremap_vcf, ".vcf.gz") + ".mpileup.vcf.gz"

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        integrate_mpileup_to_sfm \
            --sfm_vcf ~{somatic_featuremap_vcf} \
            --tumor_mpileup  ~{tumor_mpileup} \
            --normal_mpileup  ~{normal_mpileup} \
            --distance_start_to_center ~{distance_start_to_center} \
            --out_directory .

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "4 GB"
        disks: "local-disk 100 HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File out_sfm_vcf = "~{out_sfm_vcf_filename}"
        File out_sfm_vcf_index = "~{out_sfm_vcf_filename}.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task ScoreSomaticVariants {
    input{
        File input_sfm_vcf
        File input_sfm_vcf_index
        File? xgb_model
        File interval_list
        File ref_tandem_repeats
        File reference_fai
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int cpus = 4
    String out_sfm_vcf_filename = basename(input_sfm_vcf, ".vcf.gz") + ".tr_info.xgb_proba.vcf.gz"

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        #convert interval list to bed
        # Check if interval list is gzipped and handle accordingly
        if [[ ~{interval_list} == *.gz ]]; then
            gunzip -c ~{interval_list} > interval_list.txt
        else
            cp ~{interval_list} interval_list.txt
        fi
        cat interval_list.txt | grep -v "^@" | awk '{print $1"\t"$2"\t"$3}' > interval_list.bed

        somatic_featuremap_fields_transformation \
            -sfm ~{input_sfm_vcf} \
            -o ~{out_sfm_vcf_filename} \
            -i interval_list.bed \
            -ref_tr ~{ref_tandem_repeats} \
            -g ~{reference_fai} \
            ~{"-xgb_model " + xgb_model}
            
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk 100 HDD"
        docker: docker
        noAddress: no_address
        cpu: cpus
    }
    output {
        File out_sfm_vcf = "~{out_sfm_vcf_filename}"
        File out_sfm_vcf_index = "~{out_sfm_vcf_filename}.tbi"

        File monitoring_log = "monitoring.log"
    }
}   