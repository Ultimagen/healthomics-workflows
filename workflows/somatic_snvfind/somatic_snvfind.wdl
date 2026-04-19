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
# V2 workflow with two modes:
# 1. Full pipeline mode: Takes tumor and normal CRAM files directly, runs snvfind on both together per region,
#    runs classification per region, and merges final VCFs.
# 2. Classifier-only mode: Takes a pre-existing somatic featuremap VCF, runs classification per region, and merges.


# CHANGELOG in reverse chronological order
# V3: Rename from SomaticFeaturemap to SomaticSNVfind to better reflect that this workflow runs the somatic snvfind tool.
# V2.1: Added classifier-only mode for running classification on pre-existing somatic featuremap VCF
# V2: Complete rewrite - takes CRAMs directly, scatters upfront, runs snvfind with models, classifies, merges

import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/single_read_snv_tasks.wdl" as SRSNVTasks
import "tasks/structs.wdl" as Structs

workflow SomaticSNVfind {

    input {
        String pipeline_version = "1.29.2" # !UnusedDeclaration

        String base_file_name
        References references
        File interval_list
        Int num_shards
        File tandem_repeats_bed
        File xgb_model
        Float xgb_proba_threshold = 0.6
        Int? scatter_intervals_break = 10000000

        # Full pipeline mode inputs (optional when using classifier-only mode)
        Array[File]? tumor_crams
        Array[File]? tumor_cram_index_list
        Array[File]? normal_crams
        Array[File]? normal_cram_index_list
        SingleReadSNVModel? tumor_single_read_snv_model
        SingleReadSNVModel? normal_single_read_snv_model
        FeatureMapParams? featuremap_params
        Int? override_memory_gb_CreateFeatureMap

        # annotation files: only needed when generating a featuremap (not in classifier-only mode)
        FeaturemapAnnotationFiles? annotation_files

        # Classifier-only mode inputs (provide pre-existing featuremap VCF)
        File? input_somatic_snvfind_vcf
        File? input_somatic_snvfind_vcf_index

        Boolean? no_address_override
        Int? preemptible_tries_override
        File? monitoring_script_input

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv suffix(xgb_model) == ".json"
        #@wv suffix(tandem_repeats_bed) == ".bed"
        
        # Either classifier-only mode OR full pipeline mode inputs must be provided
        #@wv defined(input_somatic_snvfind_vcf) or (defined(tumor_crams) and defined(normal_crams))
        
        #  Scoring only mode validations
        #@wv defined(input_somatic_snvfind_vcf) <-> defined(input_somatic_snvfind_vcf_index)
        #@wv defined(input_somatic_snvfind_vcf) -> suffix(input_somatic_snvfind_vcf_index) == ".tbi"
        #@wv defined(input_somatic_snvfind_vcf) -> prefix(input_somatic_snvfind_vcf_index) == input_somatic_snvfind_vcf
        
        # Full pipeline mode validations (only when not in scoring only mode)
        # When full pipeline mode is active, all required inputs must be provided
        #@wv (defined(tumor_crams) and defined(normal_crams)) -> (defined(tumor_single_read_snv_model) and defined(normal_single_read_snv_model) and defined(featuremap_params) and defined(tumor_cram_index_list) and defined(normal_cram_index_list))
        #@wv defined(tumor_crams) -> (suffix(tumor_crams) <= {".bam", ".cram"})
        #@wv defined(normal_crams) -> (suffix(normal_crams) <= {".bam", ".cram"})
        #@wv defined(tumor_cram_index_list) -> (suffix(tumor_cram_index_list) <= {".crai", ".bai"})
        #@wv defined(normal_cram_index_list) -> (suffix(normal_cram_index_list) <= {".crai", ".bai"})
        #@wv (defined(tumor_crams) and defined(tumor_cram_index_list)) -> len(tumor_crams) == len(tumor_cram_index_list)
        #@wv (defined(normal_crams) and defined(normal_cram_index_list)) -> len(normal_crams) == len(normal_cram_index_list)

    }

    meta {
        description: "The SomaticSNVfind workflow performs somatic variant calling, producing an annotated somatic FeatureMap VCF file containing information on a tumor and a normal sample. It supports two execution modes: (1) Full pipeline mode, which takes matched tumor and normal CRAM files, divides the genome into regions, runs snvfind jointly on the tumor–normal pair per region using model metadata, applies somatic classification per region, and merges the results into a final VCF. (2) Classifier-only mode, which takes a pre-existing somatic FeatureMap VCF and applies the classification step only."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "no_address_override",
                "preemptible_tries_override",
                "monitoring_script_input",
                "Globals.glob",
                "ScatterIntervalList.dummy_input_for_call_caching",
                "ConcatVcfs.disk_size",
                "ConcatVcfs.preemptible_tries",
                "CreateFeatureMap.random_sample_trinuc_freq_",
                "CreateFeatureMap.cpus",
                "SomaticSNVfindClassifier.memory_gb",
                "SomaticSNVfindClassifier.cpus",
                'CreateFeatureMap.mean_coverage',
                'CreateFeatureMap.annotation_files',
                'CreateFeatureMap.max_coverage_factor'
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        references: {
            help: "Input reference genome file, index file and dict.",
            type: "Struct",
            category: "input_required"
        }
        interval_list: {
            help: "Input interval list file for scattering.",
            type: "File",
            category: "input_required"
        }
        num_shards: {
            help: "Number of shards for splitting the interval list.",
            type: "Int",
            category: "input_required"
        }
        tandem_repeats_bed: {
            help: "Reference tandem repeats file.",
            type: "File",
            category: "input_required"
        }
        xgb_model: {
            help: "XGBoost model file for somatic featuremap classification.",
            type: "File",
            category: "input_required"
        }
        xgb_proba_threshold: {
            help: "XGBoost probability threshold for somatic featuremap classification.",
            type: "Float",
            default: 0.6,
            category: "input_optional"
        }
        scatter_intervals_break: {
            help: "Break scatter intervals at multiples of this value.",
            type: "Int",
            default: 10000000,
            category: "input_optional"
        }
        tumor_crams: {
            help: "Array of input tumor CRAM files. All files must have the same sample name. Required for full pipeline mode, not needed for classifier-only mode.",
            type: "Array[File]",
            category: "input_optional"
        }
        tumor_cram_index_list: {
            help: "Array of input tumor CRAM index files, matching the order of tumor_crams. Required for full pipeline mode.",
            type: "Array[File]",
            category: "input_optional"
        }
        normal_crams: {
            help: "Array of input normal CRAM files. All files must have the same sample name. Required for full pipeline mode, not needed for classifier-only mode.",
            type: "Array[File]",
            category: "input_optional"
        }
        normal_cram_index_list: {
            help: "Array of input normal CRAM index files, matching the order of normal_crams. Required for full pipeline mode.",
            type: "Array[File]",
            category: "input_optional"
        }
        tumor_single_read_snv_model: {
            help: "Tumor SingleReadSNVModel struct containing model_metadata and model_fold_files. Required for full pipeline mode.",
            type: "Struct",
            category: "input_optional"
        }
        normal_single_read_snv_model: {
            help: "Normal SingleReadSNVModel struct containing model_metadata and model_fold_files. Required for full pipeline mode.",
            type: "Struct",
            category: "input_optional"
        }
        featuremap_params: {
            help: "FeatureMap parameters for snvfind. Required for full pipeline mode.",
            type: "Struct",
            category: "input_optional"
        }
        annotation_files: {
            type: "FeaturemapAnnotationFiles",
            help: "Annotation files for featuremap generation: dbSNP, gnomAD, and UG High Confidence Regions with their indices. Required for full pipeline mode.",
            category: "input_optional"
        }
        override_memory_gb_CreateFeatureMap: {
            type: "Int",
            help: "Override memory in GB for the CreateFeatureMap task, default: 2 (GiB). If an out of memory error occurs in the CreateFeatureMap task, try increasing this value, e.g. double it.",
            category: "input_optional"
        }
        input_somatic_snvfind_vcf: {
            help: "Pre-existing somatic featuremap VCF file. When provided, runs classifier-only mode (skips featuremap generation).",
            type: "File",
            category: "input_optional"
        }
        input_somatic_snvfind_vcf_index: {
            help: "Index file for pre-existing somatic featuremap VCF. Required when input_somatic_snvfind_vcf is provided.",
            type: "File",
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
            help: "Monitoring script override for AWS HealthOmics workflow templates multi-region support",
            type: "File",
            category: "input_optional"
        }
        somatic_snvfind_vcf: {
            help: "Output final somatic snvfind VCF file.",
            type: "File",
            category: "output"
        }
        somatic_snvfind_vcf_index: {
            help: "Output final somatic snvfind VCF index file.",
            type: "File",
            category: "output"
        }
        aggregated_parquet: {
            help: "Array of aggregated parquet files, one per genomic region/shard, primarily intended for debugging and troubleshooting. Each parquet file contains structured data from the somatic featuremap classifier including features, classification scores, and results in a columnar format. These files facilitate debugging by providing detailed intermediate results that can be inspected, analyzed, and visualized more easily than VCF format. Useful for investigating classification behavior, feature distributions, and pipeline issues.",
            type: "Array[File]",
            category: "output"
        }
    }

    Boolean no_address = select_first([no_address_override, true])
    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    # Determine if running in classifier-only mode
    Boolean classifier_only_mode = defined(input_somatic_snvfind_vcf)

    # Validate that tumor and normal samples have different sample names (only in full pipeline mode)
    if (!classifier_only_mode) {
        call ValidateTumorNormalSampleNames {
            input:
                tumor_crams = select_first([tumor_crams]),
                normal_crams = select_first([normal_crams]),
                references = references,
                docker = global.broad_gatk_docker,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script,
                no_address = true
        }
    }

    # Scatter genome into regions
    call UGGeneralTasks.ScatterIntervalList as ScatterIntervalList {
        input:
            interval_list = interval_list,
            scatter_count = num_shards,
            break_bands_at_multiples_of = select_first([scatter_intervals_break, 10000000]),
            dummy_input_for_call_caching = "",
            docker = global.broad_gatk_docker,
            no_address = true,
            monitoring_script = monitoring_script,
            convert_to_bed = true
    }

    # Process each region: CreateFeatureMap (full pipeline mode) + SomaticSNVfindClassifier
    scatter (bed_file in select_first([ScatterIntervalList.out_bed, []])) {
        # Full pipeline mode: create featuremap from CRAMs
        if (!classifier_only_mode) {
            FeatureMapParams fm_params = select_first([featuremap_params])
            # Create modified FeatureMapParams with bed_file set
            FeatureMapParams modified_featuremap_params = object {
                min_mapq : fm_params.min_mapq,
                padding_size : fm_params.padding_size,
                score_limit : fm_params.score_limit,
                max_score_to_emit : fm_params.max_score_to_emit,
                min_score_to_emit : fm_params.min_score_to_emit,
                exclude_nan_scores : fm_params.exclude_nan_scores,
                include_dup_reads : fm_params.include_dup_reads,
                keep_supplementary : fm_params.keep_supplementary,
                surrounding_quality_size : fm_params.surrounding_quality_size,
                reference_context_size : fm_params.reference_context_size,
                cram_tags_to_copy : fm_params.cram_tags_to_copy,
                attributes_prefix : fm_params.attributes_prefix,
                bed_file : bed_file,
                pileup_window_width : fm_params.pileup_window_width,
                somatic_filter_mode : fm_params.somatic_filter_mode,
                generate_random_sample : fm_params.generate_random_sample
            }
            
            call SRSNVTasks.CreateFeatureMap {
                input:
                    # do we need the mean_coverage in this case? this is needed for model training, which is not done here.
                    input_cram_bam_list = flatten([select_first([tumor_crams]), select_first([normal_crams])]),
                    input_cram_bam_index_list = flatten([select_first([tumor_cram_index_list]), select_first([normal_cram_index_list])]),
                    references = references,
                    base_file_name = base_file_name + ".shard_" + basename(bed_file, ".bed"),
                    total_aligned_bases = "0",  # Not used when generate_random_sample=false
                    random_sample_size = 0,     # Not used when generate_random_sample=false
                    featuremap_params = modified_featuremap_params,
                    annotation_files = select_first([annotation_files]),
                    model_files = [select_first([tumor_single_read_snv_model]), select_first([normal_single_read_snv_model])],
                    docker = global.featuremap_docker,
                    preemptible_tries = preemptible_tries,
                    monitoring_script = monitoring_script,
                    memory_gb = select_first([override_memory_gb_CreateFeatureMap, 32]),
                    cpus = 2,
            }
        }

        # In classifier-only mode, define the interval for region filtering
        if (classifier_only_mode) {
            File classifier_bed_file = bed_file
        }

        # Determine which VCF to use for classifier
        File shard_vcf = select_first([CreateFeatureMap.featuremap, input_somatic_snvfind_vcf])
        File shard_vcf_index = select_first([CreateFeatureMap.featuremap_index, input_somatic_snvfind_vcf_index])

        call SomaticSNVfindClassifier {
            input:
                somatic_snvfind_vcf = shard_vcf,
                somatic_snvfind_vcf_index = shard_vcf_index,
                ref_fasta_index = references.ref_fasta_index,
                tandem_repeats_bed = tandem_repeats_bed,
                xgb_model = xgb_model,
                xgb_proba_threshold = xgb_proba_threshold,
                regions_bed_file =classifier_bed_file,
                docker = global.ugbio_featuremap_docker,
                preemptible_tries = preemptible_tries,
                monitoring_script = monitoring_script,
                no_address = no_address
        }
    }

    # Merge scattered VCF outputs
    call UGGeneralTasks.ConcatVcfs as ConcatVcfs {
        input:
            input_vcfs = SomaticSNVfindClassifier.classified_vcf,
            input_vcfs_indexes = SomaticSNVfindClassifier.classified_vcf_index,
            output_vcf_name = base_file_name + ".somatic_featuremap.vcf.gz",
            docker = global.broad_gatk_docker,
            preemptible_tries = preemptible_tries,
            no_address = no_address,
            monitoring_script = monitoring_script,
            disk_size = ceil(6*size(SomaticSNVfindClassifier.classified_vcf,"GB")+5)
    }

    output {
        File somatic_snvfind_vcf = ConcatVcfs.output_vcf
        File somatic_snvfind_vcf_index = ConcatVcfs.output_vcf_index
        Array[File] aggregated_parquet = SomaticSNVfindClassifier.aggregated_parquet
    }
}

task SomaticSNVfindClassifier {
    input {
        File somatic_snvfind_vcf
        File somatic_snvfind_vcf_index
        File ref_fasta_index
        File tandem_repeats_bed
        File xgb_model
        Float xgb_proba_threshold
        File? regions_bed_file
        String docker
        Int preemptible_tries
        File monitoring_script
        Boolean no_address
        Int memory_gb = 4
        Int cpus = 2
    }

    Float vcf_size = size(somatic_snvfind_vcf, "GB")
    Int disk_size = ceil(vcf_size * 2 + 20)

    String out_vcf = basename(somatic_snvfind_vcf, ".vcf.gz") + ".classified.vcf.gz"
    String out_parquet = basename(somatic_snvfind_vcf, ".vcf.gz") + ".aggregated.parquet"

    command <<<
        set -xeuo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        somatic_snvfind_classifier \
            --somatic-snvfind-vcf ~{somatic_snvfind_vcf} \
            --output-vcf ~{out_vcf} \
            --genome-index-file ~{ref_fasta_index} \
            --tandem-repeats-bed ~{tandem_repeats_bed} \
            --xgb-model-json ~{xgb_model} \
            --n-threads ~{cpus} \
            ~{if defined(regions_bed_file) then "--regions-bed-file ~{regions_bed_file}" else ""} \
            --xgb-proba-thresh ~{xgb_proba_threshold} \
            --output-parquet ~{out_parquet} \
            --verbose
    >>>

    runtime {
        preemptible: preemptible_tries
        docker: docker
        cpu: cpus
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        noAddress: no_address
    }

    output {
        File classified_vcf = "~{out_vcf}"
        File classified_vcf_index = "~{out_vcf}.tbi"
        File monitoring_log = "monitoring.log"
        File aggregated_parquet = "~{out_parquet}"
    }
}

task ValidateTumorNormalSampleNames {
    input {
        Array[File] tumor_crams
        Array[File] normal_crams
        References references
        String docker
        Int preemptible_tries
        File monitoring_script
        Boolean no_address
    }

    parameter_meta {
        tumor_crams: {
            localization_optional: true
        }
        normal_crams: {
            localization_optional: true
        }
    }

    command <<<
        set -ex
        set -o pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Extract sample names from all tumor CRAM files
        > tumor_sample_names.tmp.txt
        for cram in ~{sep=' ' tumor_crams}; do
            gatk GetSampleName \
                -I "$cram" \
                -R ~{references.ref_fasta} \
                -O sample_name.tmp.txt
            cat sample_name.tmp.txt >> tumor_sample_names.tmp.txt
        done
        sort tumor_sample_names.tmp.txt | uniq > tumor_sample_names.txt

        # Extract sample names from all normal CRAM files
        > normal_sample_names.tmp.txt
        for cram in ~{sep=' ' normal_crams}; do
            gatk GetSampleName \
                -I "$cram" \
                -R ~{references.ref_fasta} \
                -O sample_name.tmp.txt
            cat sample_name.tmp.txt >> normal_sample_names.tmp.txt
        done
        sort normal_sample_names.tmp.txt | uniq > normal_sample_names.txt

        # Verify all tumor CRAMs have the same sample name
        if [[ $(wc -l < tumor_sample_names.txt) -ne 1 ]]
        then
            echo "Error: tumor CRAM files have different sample names:" >&2
            cat tumor_sample_names.txt >&2
            exit 1
        fi

        # Verify all normal CRAMs have the same sample name
        if [[ $(wc -l < normal_sample_names.txt) -ne 1 ]]
        then
            echo "Error: normal CRAM files have different sample names:" >&2
            cat normal_sample_names.txt >&2
            exit 1
        fi

        # Read sample names
        tumor_sample=$(cat tumor_sample_names.txt)
        normal_sample=$(cat normal_sample_names.txt)

        # Validate that sample names are different
        if [[ "$tumor_sample" == "$normal_sample" ]]
        then
            echo "Error: Tumor and normal samples have the same sample name: $tumor_sample" >&2
            echo "Tumor sample name: $tumor_sample" >&2
            echo "Normal sample name: $normal_sample" >&2
            exit 1
        fi

        # Write final sample names for output
        echo "$tumor_sample" > tumor_sample_name.txt
        echo "$normal_sample" > normal_sample_name.txt

        echo "Validation passed: Tumor sample ($tumor_sample) and normal sample ($normal_sample) have different names" >&2
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        preemptible: preemptible_tries
        noAddress: no_address
        docker: docker
    }

    output {
        File monitoring_log = "monitoring.log"
        String tumor_sample_name = read_string("tumor_sample_name.txt")
        String normal_sample_name = read_string("normal_sample_name.txt")
    }
}