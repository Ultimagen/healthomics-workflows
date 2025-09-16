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
#   Structures definitions used mainly in ug_variant_calling_pipeline.
# CHANGELOG in reverse chronological order
#
struct ContaminationSites {
  String contamination_sites_path
  File contamination_sites_vcf
  File contamination_sites_vcf_index
}

struct References {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
}

# BWA alignment
struct AlignmentReferences {
  References references
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
}

struct GiraffeReferences {
  References references
  File ref_gbz
  File ref_dist
  File ref_min
  File ref_path_list
}

# BWA-METH alignment
struct BwaMethReferences {
  File bwa_meth_ref_idx
}

# UA alignment
struct UaParameters {
  File? ua_index
  File? ref_alt
  String ua_extra_args
  Boolean v_aware_alignment_flag
  Int? cpus
  Int? memory_gb
}

# UA-METH alignment
struct UaMethParameters {
  File index_c2t
  File index_g2a
  Int? cpus
}

struct VariantCallingSettings {
  File wgs_calling_interval_list
  Int break_bands_at_multiples_of
  Int haplotype_scatter_count
}

struct VcfPostProcessing {
  Array[File] annotation_intervals
  File? filtering_model_no_gt
  File af_only_gnomad
  File af_only_gnomad_index
  Boolean filter_cg_insertions
  File? filtering_blocklist_file
  File? training_blocklist_file
  Int? exome_weight
  String? exome_weight_annotation
  File? interval_list_override
  File runs_file
  String? filtering_model_with_gt_name_override
  Float max_duplication_in_reasonable_sample
  Float max_chimerism_in_reasonable_sample
  File ref_dbsnp
  File ref_dbsnp_index
  File wgs_coverage_interval_list
}

struct SimpleReadTrimmingParameters {
  String trimming_tool  # allowed values cutadapt or trimmer
  String adapter_5p
  String adapter_3p
  # the 5p and 3p adapters are mandatory because the python script in trimmer requires it, can be changed if needed
  Int? umi_length_3p
  Int? umi_length_5p
  Float max_error_rate_5p
  Float max_error_rate_3p
  Int min_overlap_5p
  Int min_overlap_3p
  Int? min_insert_length
  Int? max_insert_length
  # in practice only trimmer supports it - current implementation of cutadapt can only discard through ClipReads and
  # it's not the behavior we want

  String untrimmed_reads_action
  # in practice only trimmer supports it - either "" (do nothing), "filter" (mark in sam flag) or "discard", cutadapt
  # only applies a tag to trimmed reads - XF:Z:A
}

struct TrimmerParameters {
  File? formats_description         # optional formtas.json file. If not provided, the default trimmer formats will be used (https://github.com/Ultimagen/trimmer/blob/master/formats/formats.json)
  String? local_formats_description # path to description file stored in the docker, default is /trimmer/formats/formats.json
  String? untrimmed_reads_action    # either "" (do nothing), "filter" (mark in sam flag) or "discard"
  String? format                    # format name to be used, as defined in the formats.json file
  String? failure_read_group        # name of field to be used to mark failed reads, generally "unmatched" (will show as "rg:Z:unmatched" in the cram)
  String? minor_read_group          # name of the minor read group to be used in the output file
  String? extra_args                # extra arguments to be passed to trimmer
  Array[File]? pattern_files        # If the formats description uses pattern files (specified with file: or enc: prefix in format file), they must be provided also here.
                                    # Note that built-in pattern files from the Trimmer repo do not need to be provided here.
  File? cram_reference              # Use this in case the reference is not embedded in the input cram (or if you want to read the cram with a different reference)
  Int? memory_gb                    # Ovverride the default memory (in GB) used by trimmer
  String? output_demux_format       # If provided, trimmer will perform demultiplexing and the output file name will be in this format. E.g. "output-2%"
  String? filename_prefix_sep       # default is set to "_", but this can be changed to "-" as in the case of ancient DNA trimmed output
  String? output_failed_file_name_suffix # Setting the name for the failed reads file. Default is "failed.cram", if this parameter is given will be <base_file_name>_<output_failed_file_name_suffix>
  Boolean? remove_small_files       # Allows to remove trimmed files below 1Gb in size
  Boolean? add_run_id_as_ri_tag     # If true, will extract the run id from the input file and add it as a tag to the output reads.
}

struct ReferenceDbSnp {
    File ref_dbsnp
    File ref_dbsnp_index
}

struct FeatureMapParams { 
  # snvfind parameters
  Int? min_mapq                     # -q minimum mapping quality. defaults to 20
  Int? padding_size                 # -p padding size. defaults to 5
  Int? score_limit                  # -L score limit
  Int? max_score_to_emit            # -X max score to emit
  Int? min_score_to_emit            # -N min score to emit
  Boolean? exclude_nan_scores       # -n exclude nan scores
  Boolean? include_dup_reads        # -d include dup reads
  Boolean? keep_supplementary       # -k keep supplementary alignments
  Int? surrounding_quality_size     # -Q surrounding median and mean quality size. defaults to pad
  Int? reference_context_size       # -r reference context size, defualts to 3
  Array[String]? cram_tags_to_copy  # -c list of attributes to copy from sam to vcf
  String? attributes_prefix         # -C prefix for copied attributes
  File? bed_file                    # -b bed file containing ranges to process
}

struct SingleReadSNVParams {
    Int tp_train_set_size
    Int fp_train_set_size
    Float tp_train_set_size_sampling_overhead
    Int random_seed
    Int num_CV_folds
    Int min_coverage_filter
    Float max_coverage_factor
    Float max_vaf_for_fp  # Maximum VAF for false positive filtering
    Array[String] pre_filters  # Pre-filter configuration in format "name=X:field=Y:op=Z:value=W:type=T"
}

struct MrdAnalysisParams {
  String signature_filter_query
  String read_filter_query
  String? tumor_sample # Optional, used to specify the tumor sample name in singature vcf
}

struct StarsoloBamParams {
  Boolean save_bam
  Boolean? sort_bam_override
  String? limit_sort_ram_override
  Boolean? include_unmapped_override
  Boolean? transcript_override
}

struct StarSoloParams {
  # STARsolo params
  Array[String] out_sam_attributes # 'NH', 'HI', 'AS', 'nM', 'MD', 'jM', 'jI', 'XS', 'MC', 'ch', 'CR', 'UR', 'CY', 'UY', 'GX', 'GN'
  String strand                    # 'Reverse' or 'Forward'
  Int cell_barcode_length
  Int umi_length
  File barcode_whitelist
  String cell_filter              # EmptyDrops_CR
  String? library_direction       # should be "five_prime" or "three_prime"
  # general STAR params
  File? genome
  File? gtf_override
  String? extra_args
}

struct StarSoloOutputs {
  File? genome_zip_output
  File output_bam
  File gene_features_stats
  File gene_summary_csv
  File gene_umi_per_cell_sorted
  File gene_filtered_features
  File gene_filtered_barcodes
  File gene_filtered_matrix
  File star_log_file
  File star_log_params_file
  File barcode_file
  File star_stats
  File gathered_star_stats_csv
}

struct RuntimeParams {
  Int? memory_gb_override
  String? disk_type_override
  Float? disk_size_gb_override
  Int? cpu_num_override
  Boolean? no_address_override
  Int? preemptible_tries_override
  Int? max_retries_override
  String? gitc_path_override
}

struct Adapters10x {
  String adapter_5p
  String adapter_3p
  String adapter_middle
  Float adapter_min_error_rate_5p
  Int adapter_min_overlap_5p
  Float adapter_min_error_rate_middle
  Int adapter_min_overlap_middle
  Float adapter_min_error_rate_3p
  Int adapter_min_overlap_3p
}

struct StarsoloScores {
  Float? starsolo_score_del_open_override
  Float? starsolo_score_ins_open_override
  Float? starsolo_score_del_base_override
  Float? starsolo_score_ins_base_override
}

struct SentieonTNScopeParameters {
  String read_filter
  String tnscope_extra_args
  String bcftools_annotate_remove
  String bcftools_view_include
}

struct StarGenomeGenerateParams {
  Array[File] fasta_files
  File gtf_file
  String output_basename
  String? extra_args
}

struct TrimAlignSortSteps {
  Boolean? trim
  Boolean? align
  Boolean? sort
}

struct SorterParams {
  Boolean mark_duplicates
  String? umi_tag           # multiple tags should be separated by comma e.g. "u5,u3"
  Boolean? aligned          # demux arg to mentioned if the data aligned. The default is true.
  String? output_group      # Define a custom read-group e.g. "majorRG-minorRG"  (instead of the default "majorRG"). majorRG is the value of the RG tag of each read. See sorter documentation for more details.
  String? output_path       # Define the output path for a custom read-group. Default is: {outputGroup}/{outputGroup} !NOTE! the path must include a subfolder
  Float? downsample_frac    # Downsample fraction (0.0-1.0) to be used in Demux
  Int? downsample_seed      # Downsample seed to be used in Demux
  Int? mark_duplicates_ends_read_uncertainty   # Number of bases of uncertainty in read ends position to use when marking duplicates
  Boolean? mark_duplicates_flow_use_clipped_location  # If true, use the clipped location of the read to mark duplicates, otherwise add the softclip length to the alignment end position
  Boolean? mark_duplicates_flow_q_is_known_end  # If true, the ends in quality trimmed reads are treated as known when marking duplicates. Otherwise, the ends are treated as unknown so any end position is matched.
  String? demux_extra_args
  String? sort_extra_args
  Int? memory_gb            # Override the default memory (in GB) used by sorter
  Int? demux_memory_gb            # Override the default memory (in GB) used by demux
  Int? demux_cpu           # Override the default cpu used by demux
  File? coverage_intervals  # tar.gz file with the coverage intervals tsv pointing to the relevant coverage intervals files
  File? single_cell_cbc_classifier  # single cell classifier model (json)
}

struct SingleCellQcThresholds {
  Int pass_trim_rate
  Int read_length
  Int fraction_below_read_length
  Int percent_aligned
}
