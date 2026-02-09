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
# Runs substitution error analysis for Ultima Genomics data. Includes three parts:
# 1. Substitution error analysis, including FeatureMap generation
# 2. Optional generation of a new MRD somatic signature using the somatic_mutations_pipeline.wdl
# 3. Intersection of signature from #2 and any other input signatures with FeatureMap from #1 and prepartion of data
#    for MRD analysis

# CHANGELOG in reverse chronological order
# 1.9.0 Removed the option to generate new signature
# 1.7.0 Now using SingleReadSNV pipeline instead of SubstitutionAnalysis, added a GenerateControlSignaturesFromDatabase task
# 1.3.1 Initial release

import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/mrd.wdl" as UGMrdTasks
import "tasks/structs.wdl" as Structs  # !UnusedImport
import "tasks/globals.wdl" as Globals

workflow MRDFeatureMap {
    input {
        String pipeline_version = "1.26.1" # !UnusedDeclaration
        String base_file_name
        # Outputs from single_read_snv.wdl (cfDNA sample)
        File cfdna_featuremap
        File cfdna_featuremap_index
        File featuremap_df_file
        File srsnv_metadata_json
        # for coverage collection with gatk DepthOfCoverage (cfDNA sample)
        File cfdna_cram_bam
        File cfdna_cram_bam_index
        Int mapping_quality_threshold
        # Somatic signature files, matched from the same patient as the cfDNA sample, or control signatures
        Array[File]? external_matched_signatures
        Array[File]? external_control_signatures
        # filter signatures
        String? bcftools_extra_args
        Array[File] include_regions
        Array[File] exclude_regions
        # diluent germline vcfs
        Array[File]? diluent_germline_vcfs
        # final-analysis-level filters
        MrdAnalysisParams mrd_analysis_params
        # for generating db control signatures
        File? snv_database
        Int? n_synthetic_signatures
        # option to increase memory for the coverage extraction task
        Int? memory_extract_coverage_override

        References references
        
        Int? preemptible_tries
        Boolean? no_address_override
        # Used for running on other clouds (aws)
        String? cloud_provider_override
        File? monitoring_script_input

        Boolean create_md5_checksum_outputs = false

        # When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your
        # own workspace. This will prevent call caching from failing with "Cache Miss (10 failed copy attempts)". Outside of
        # Terra this can be left as the default empty String. This dummy input is only needed for tasks that have no inputs
        # specific to the sample being run (such as GetBwaVersion which does not take in any sample data).
        String dummy_input_for_call_caching = ""

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)

        # signatures
        #@wv defined(external_matched_signatures) or defined(external_control_signatures)
        #@wv defined(external_matched_signatures) -> len(external_matched_signatures) > 0
        #@wv defined(external_control_signatures) -> len(external_control_signatures) > 0
        #@wv defined(snv_database) <-> defined(n_synthetic_signatures)
        #@wv defined(snv_database) -> suffix(snv_database) == '.vcf' or suffix(prefix(snv_database)) == '.vcf'

        # references
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
        #@wv suffix(references['ref_dict']) == '.dict'
        #@wv suffix(references['ref_fasta_index']) == '.fai'
        #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']

        # input suffixes
        #@wv prefix(basename(cfdna_cram_bam_index)) == basename(cfdna_cram_bam)
        #@wv suffix(cfdna_cram_bam) in {'.bam', '.cram'}
        #@wv suffix(cfdna_cram_bam_index) in {'.bai', '.crai'}
        #@wv suffix(include_regions) <= {'.bed', '.gz'}
        #@wv suffix(exclude_regions) <= {'.bed', '.gz', '.vcf', '.vcf.gz'}
        #@wv suffix(snv_database) in {'.gz'}
        #@wv suffix(prefix(snv_database)) in {'.vcf', '.vcf.gz'}
        #@wv suffix(featuremap_df_file) in {'.parquet'}
        ##@wv defined(diluent_germline_vcfs) and len(diluent_germline_vcfs)>0 -> suffix(diluent_germline_vcfs) <= {'.vcf', '.vcf.gz'}
    }

  meta {
      description : "The UG pipeline for tumor informed MRD measures the ctDNA variant allele fraction in cfDNA from the presence of tumor-specific SNVs. The input data is generally 3 aligned cram files:\n\n- cfDNA (plasma)\n\n- Tumor tissue (FFPE / FF)\n\n- Normal tissue (buffy coat / PBMCs)\n\n\n\nIt is possible to provide the cfDNA cram file only with an existing somatic vcf file.\n\n\n\nThe analysis of this MRD data is composed of three parts:\n\n1. Tumor signature mutation calling, where the tumor and normal tissues are used for finding the tumor somatic mutations signature with somatic variant calling (by default UG Somatic DeepVariant, though these can be provided from other callers).\n\n2. Single Read SNV pipeline, where all the SNV candidates compared to the reference genome are extracted from the cfDNA cram file to a FeatureMap vcf, annotated and assigned a quality score (SNVQ).\n\n3. Intersection and MRD data analysis, where the FeatureMap and signature are intersected and filtered, then reads supporting the tumor mutations are counted and a ctDNA variant allele fraction is measured. Control signatures can be added to estimate the background noise, e.g. from other cohort patients, and in addition control signatures are generated from a somatic mutation database. \n\n\n\nThis pipeline describes step #3, the intersection and MRD data analysis, once #1 and #2 are completed. The following input templates are available for different kinds of input data: 1) `mrd_featuremap_template-Matched-signature-with-cohort-with-quality-filtering.json` | Use this template to run MRD using a matched signature, including cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and control signatures (suitable for EfficientDV output). 2) `mrd_featuremap_template-Matched-signature-without-cohort-with-quality-filtering.json` | Use this template to run MRD using a matched signature, without cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and signature (suitable for EfficientDV output). 3) `mrd_featuremap_template-Matched-signature-with-cohort-without-quality-filtering.json` | Use this template to run MRD using a matched signature, with cohort controls (non-matched mutation signature vcf). Quality filtering is not applied to matched and signature (suitable for vcf files without a QUAL field). 4) `mrd_featuremap_template-Healthy-without-matched-with-cohort-with-quality-filtering.json` | Use this template to run MRD on a healthy control plasma without a matched signature, with cohort controls (non-matched mutation signature vcf). Quality filtering is applied to matched and signature (suitable for EfficientDV output)."
      author: "Ultima Genomics"
      WDL_AID: { exclude: [
          "pipeline_version",
          "cloud_provider_override",
          "no_address_override",
          "preemptible_tries",
          "dummy_input_for_call_caching",
          "monitoring_script_input",
          "Globals.glob",
          "Sentieon.Globals.glob",
          "AnnotateVCF.Globals.glob",
          "SingleSampleQC.Globals.glob",
          "VariantCallingEvaluation.Globals.glob",
          "MergeVcfsIntoBed.disk_size",
          "PadDiluentVcf.disk_size",
          "PadDiluentVcf.memory_gb",
          "PadDiluentVcf.cpus",
          "MergeMd5sToJson.output_json"
      ]}
  }    

  parameter_meta {
      base_file_name: {
          help: "Base file name for output files.",
          type: "String", 
          category: "input_required"
      }
      cfdna_featuremap: {
          help: "FeatureMap file generated from the cfDNA CRAM file by the SingleReadSNV pipeline",
          type: "File",
          category: "input_required"
      }
      cfdna_featuremap_index: {
          help: "Respective index",
          type: "File",
          category: "input_required"
      }
      featuremap_df_file: {
          help: "Parquet file of the FeatureMap model data generated by the SingleReadSNV pipeline",
          type: "File",
          category: "input_required"
      }
      srsnv_metadata_json: {
          help: "Metadata json file output by the SingleReadSNV pipeline",
          type: "File",
          category: "input_required"
      }
      cfdna_cram_bam: {
          help: "CRAM or BAM file of cfDNA sample, must be equal to the input used for the SingleReadSNV pipeline that generated cfdna_featuremap",
          type: "File",
          category: "input_required"
      }
      cfdna_cram_bam_index: {
          help: "Respective index",
          type: "File",
          category: "input_required"
      }
      mapping_quality_threshold: {
          help: "Mapping quality threshold for reads to be included in the coverage analysis, default 0 as srsnv mapq filtering is included in srsnv_metadata_json",
          type: "Int",
          category: "input_required"
      }
      external_matched_signatures: {
          help: "Optional signatures matched to the patient from whom the cfDNA sample was taken, to be used as matched signatures in the MRD analysis, leave blank if running a healthy control donor",
          type: "Array[File]",
          category: "input_optional"
      }
      external_control_signatures: {
          help: "Optional signatures from other individuals to be used as controls in the MRD analysis",
          type: "Array[File]",
          category: "input_optional"
      }
      bcftools_extra_args: {
          help: "Extra argument given after the 'bcftools view' command that filters the matched and control signatures. E.g. '-f PASS' to filter only PASS variants.",
          type: "String",
          category: "input_required"
      }
      include_regions: {
          help: "Genomic regions to include in the analysis, bed files are accepted",
          type: "Array[File]",
          category: "ref_required"
      }
      exclude_regions: {
          help: "Genomic regions to exclude from the analysis, bed and vcf[.gz] files are accepted",
          type: "Array[File]",
          category: "ref_optional"
      }
      diluent_germline_vcfs: {
          help: "Optional germline VCF files from diluent samples to pad and exclude from analysis",
          type: "Array[File]",
          category: "input_optional"
      }
      mrd_analysis_params: {
          help: "Parameters for the MRD analysis",
          type: "MrdAnalysisParams",
          category: "input_required"
      }
      snv_database: {
          help: "Somatic mutation database to use for generating control signatures",
          type: "File",
          category: "input_optional"
      }
      n_synthetic_signatures: {
          help: "Number of synthetic signatures to generate from the somatic mutation database. Set to 0 to disable generation of control signatures from the database.",
          type: "Int",
          category: "input_optional"
      }
      memory_extract_coverage_override: {
            help: "Memory in GB to use for the coverage extraction task",
            type: "Int",
            category: "input_optional"
      }
      create_md5_checksum_outputs: {
           help: "Create md5 checksum for requested output files",
           type: "Boolean",
           category: "input_optional"
      }
      monitoring_script_input: {
          help: "Monitoring script override for AWS HealthOmics workflow templates multi-region support",
          type: "File",
          category: "input_optional"
      }
      references: {
          help: "Reference files: fasta, dict and fai, recommended value set in the template",
          type: "References",
          category: "ref_required"
      }
      features_dataframe: {
          help: "Parquet file of the FeatureMap dataframe after intersection with the matched and control signatures",
          type: "File",
          category: "output"
      }
      signatures_dataframe: {
          help: "Parquet file of the matched and control signatures dataframe",
          type: "File",
          category: "output"
      }
      mrd_analysis_notebook: {
          help: "Jupyter notebook of the MRD analysis",
          type: "File",
          category: "output"
      }
      report_html: {
          help: "HTML report of the MRD analysis",
          type: "File",
          category: "output"
      }
      ctdna_vaf_h5: {
          help: "HDF5 file of the ctDNA VAF and other results of the MRD analysis",
          type: "File",
          category: "output"
      }
      intersected_featuremaps_parquet: {
          help: "Parquet file of the intersected FeatureMap",
          type: "Array[File]",
          category: "output"
      }
      intersected_featuremaps: {
          help: "VCF file of the intersected FeatureMap",
          type: "Array[File]",
          category: "output"
      }
      intersected_featuremaps_indices: {
          help: "Respective indices",
          type: "Array[File]",
          category: "output"
      }
      control_signatures_vcf: {
          help: "VCF file of the filtered control signatures",
          type: "Array[File]",
          category: "output"
      }
      matched_signatures_vcf: {
          help: "VCF file of the filtered matched signatures",
          type: "Array[File]",
          category: "output"
      }
      db_signatures_vcf: {
          help: "VCF file of the filtered db control signatures",
          type: "Array[File]",
          category: "output"
      }
      coverage_bed: {
        help: "Coverage bed file",
        type: "File",
        category: "output"
      }
      coverage_bed_index: {
        help: "Respective index",
        type: "File",
        category: "output"
      }
      md5_checksums_json: {
        help: "json file that will contain md5 checksums for requested output files",
        type: "File",
        category: "output"
    }
  }

  Int preemptibles = select_first([preemptible_tries, 1])

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers

  Boolean defined_external_matched_signatures = defined(external_matched_signatures)
  Boolean defined_external_control_signatures = defined(external_control_signatures)
  Boolean defined_somatic_snv_database = defined(snv_database) && (select_first([n_synthetic_signatures]) > 0)
  Boolean defined_diluent_germline_vcfs = defined(diluent_germline_vcfs)
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

  # Part 1 - Filter signatures
  # Process diluent germline VCFs if provided
  if (defined_diluent_germline_vcfs) {
    Array[File] diluent_vcfs = select_first([diluent_germline_vcfs])
    scatter (vcf_path in diluent_vcfs) {
      call UGMrdTasks.PadVcf as PadDiluentVcf {
        input:
          input_vcf = vcf_path,
          ref_fai = references.ref_fasta_index,
          pad_size = 2,
          docker = global.ugbio_core_docker,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script
      }
    }
  }

  if (defined_external_matched_signatures) {
    Array[File] all_matched_signatures = select_first([external_matched_signatures,])
    Array[Array[File]] matched_exclude_regions_array = select_all([exclude_regions, PadDiluentVcf.padded_bed])
    Array[File] matched_exclude_regions = flatten(matched_exclude_regions_array)
    # Filter the matched signatures over the genomic regions bed file + apply filters
    scatter (i in range(length(all_matched_signatures))) {
      call UGGeneralTasks.FilterVcfWithBcftools as FilterMatched {
        input:
          input_vcf = all_matched_signatures[i],
          docker = global.bcftools_docker,
          bcftools_extra_args = bcftools_extra_args,  
          exclude_regions = matched_exclude_regions,
          include_regions = include_regions,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script, #!FileCoercion
      }
    }
  }  
  
  # for the external and db controls, exclude the regions from the matched signatures
  Array[Array[File]] control_exclude_regions_array = select_all([exclude_regions, external_matched_signatures, PadDiluentVcf.padded_bed])
  Array[File] control_exclude_regions = flatten(control_exclude_regions_array) 
  if (defined_external_control_signatures) {
    Array[File] external_control_signatures_array = select_first([external_control_signatures,])
    scatter (i in range(length(external_control_signatures_array))) {
      call UGGeneralTasks.FilterVcfWithBcftools as FilterControlSignatures {
        input:
          input_vcf = external_control_signatures_array[i],
          docker = global.bcftools_docker,
          bcftools_extra_args = bcftools_extra_args,  
          exclude_regions = control_exclude_regions,
          include_regions = include_regions,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script, #!FileCoercion
      }
    }
  }

  if (defined_somatic_snv_database) {
    call UGGeneralTasks.FilterVcfWithBcftools as FilterDb {
      input:
        input_vcf = select_first([snv_database]),
        docker = global.bcftools_docker,
        exclude_regions = control_exclude_regions,
        include_regions = include_regions,
        preemptible_tries = preemptibles,
        monitoring_script = monitoring_script, #!FileCoercion
    }

    # use the first matched signature as the reference for the db signatures, unless it is not given, then use the first control signature
    Array[File] filtered_signature_files = select_first([FilterMatched.output_vcf, FilterControlSignatures.output_vcf])
    File filtered_signature_file = filtered_signature_files[0]
    call UGMrdTasks.GenerateControlSignaturesFromDatabase as GenerateControlSignaturesFromDatabase {
      input:
        signature_file = filtered_signature_file,
        snv_database = FilterDb.output_vcf,
        n_synthetic_signatures = select_first([n_synthetic_signatures]),
        ref_fasta = references.ref_fasta,
        ref_fasta_index = references.ref_fasta_index,
        ref_dict = references.ref_dict,
        docker = global.ugbio_mrd_docker,
        disk_size = 2 * (size(snv_database, "GB") + size(FilterDb.output_vcf, "GB")) + 10,
        memory_gb = 8,
        cpus = 2,
        monitoring_script = monitoring_script #!FileCoercion
    }
  }

  Array[File]? filtered_db_control_signatures = GenerateControlSignaturesFromDatabase.db_signatures
  Array[File]? filtered_external_control_signatures = FilterControlSignatures.output_vcf
  Array[File]? filtered_matched_control_signatures = FilterMatched.output_vcf

  # Part 2 - Collect coverage over signatures
  Array[File] all_vcf_files = flatten(select_all([filtered_matched_control_signatures, filtered_external_control_signatures, filtered_db_control_signatures]))
  Int memory_extract_coverage = select_first([memory_extract_coverage_override, 8])

  call UGMrdTasks.MergeVcfsIntoBed as MergeVcfsIntoBed {
    input:
      vcf_files = all_vcf_files,
      docker = global.ugbio_core_docker,
      disk_size = ceil(size(all_vcf_files, "GB") * 2 + 20),
      memory_gb = memory_extract_coverage,
      cpus = 2,
      preemptibles = preemptibles,
      monitoring_script = monitoring_script #!FileCoercion
  }

  call UGMrdTasks.ExtractCoverageOverVcfFiles as ExtractCoverageOverVcfFiles{
    input:
      merged_loci_bed = MergeVcfsIntoBed.merged_loci_bed,
      input_cram_bam = cfdna_cram_bam,
      input_cram_bam_index = cfdna_cram_bam_index,
      base_file_name = base_file_name,
      mapping_quality_threshold = mapping_quality_threshold,
      references = references,
      docker = global.mosdepth_docker,
      memory_gb = memory_extract_coverage,
      cpus = 2,
      preemptibles = preemptibles,
      monitoring_script = monitoring_script  #!FileCoercion
  }

  # Part 3 - Intersect FeatureMap with signatures
  Float featuremap_size = size(cfdna_featuremap, "GB")
  Int tmp_cpus_FeatureMapIntersectWithSignatures = round(length(all_vcf_files) / 2)
  Int cpus_FeatureMapIntersectWithSignatures = if tmp_cpus_FeatureMapIntersectWithSignatures < 2 then 2 else tmp_cpus_FeatureMapIntersectWithSignatures
  Int tmp_memory_FeatureMapIntersectWithSignatures = cpus_FeatureMapIntersectWithSignatures * 2
  Int memory_FeatureMapIntersectWithSignatures = if tmp_memory_FeatureMapIntersectWithSignatures < 4 then 4 else tmp_memory_FeatureMapIntersectWithSignatures
  call UGMrdTasks.FeatureMapIntersectWithSignatures as FeatureMapIntersectWithSignatures{
    input:
      featuremap = cfdna_featuremap,
      featuremap_index = cfdna_featuremap_index,
      matched_signatures = FilterMatched.output_vcf,
      matched_signatures_indexes = FilterMatched.output_vcf_index,
      control_signatures = FilterControlSignatures.output_vcf,
      control_signatures_indexes = FilterControlSignatures.output_vcf_index,
      db_signatures = GenerateControlSignaturesFromDatabase.db_signatures,
      matched_signatures_indices = FilterMatched.output_vcf_index,
      control_signatures_indices = FilterControlSignatures.output_vcf_index,
      db_signatures_indices = GenerateControlSignaturesFromDatabase.db_signatures_indices,
      docker = global.ugbio_mrd_docker,
      disk_size = 2 * featuremap_size + 30,
      memory_gb = memory_FeatureMapIntersectWithSignatures,
      cpus = cpus_FeatureMapIntersectWithSignatures,
      monitoring_script = monitoring_script  #!FileCoercion
  }

  # Part 4 - Integrate all the processed data in the MRD data analysis
  call UGMrdTasks.MrdDataAnalysis as MrdDataAnalysis{
    input:
      intersected_featuremaps_parquet = FeatureMapIntersectWithSignatures.intersected_featuremaps_parquet,
      matched_signatures_vcf = filtered_matched_control_signatures,
      control_signatures_vcf = filtered_external_control_signatures,
      db_signatures_vcf = filtered_db_control_signatures,
      coverage_bed = ExtractCoverageOverVcfFiles.coverage_bed,
      mrd_analysis_params = mrd_analysis_params,
      basename = base_file_name,
      featuremap_df_file = featuremap_df_file,
      srsnv_metadata_json = srsnv_metadata_json,
      docker = global.ugbio_mrd_docker,
      disk_size = 3 * featuremap_size + 30,
      memory_gb = memory_extract_coverage,
      cpus = 4,
      monitoring_script = monitoring_script  #!FileCoercion
  }

    File features_dataframe_ = MrdDataAnalysis.features
    File signatures_dataframe_ = MrdDataAnalysis.signatures
    File report_html_ = MrdDataAnalysis.mrd_analysis_html
    File ctdna_vaf_h5_ = MrdDataAnalysis.ctdna_vaf_h5

    if (create_md5_checksum_outputs) {

        Array[File] output_files = select_all(flatten([
                                                      select_first([[features_dataframe_], []]),
                                                      select_first([[signatures_dataframe_], []]),
                                                      select_first([[report_html_], []]),
                                                      select_first([[ctdna_vaf_h5_], []]),
                                                      ]))

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
    # MRD analysis results
    File features_dataframe = features_dataframe_
    File signatures_dataframe = signatures_dataframe_
    File report_html = report_html_
    File ctdna_vaf_h5 = ctdna_vaf_h5_
    # Intersected featuremaps
    Array[File] intersected_featuremaps_parquet = FeatureMapIntersectWithSignatures.intersected_featuremaps_parquet
    Array[File] intersected_featuremaps = FeatureMapIntersectWithSignatures.intersected_featuremaps
    Array[File] intersected_featuremaps_indices = FeatureMapIntersectWithSignatures.intersected_featuremaps_indices
    # filtered signatures
    Array[File]? control_signatures_vcf = filtered_external_control_signatures
    Array[File]? matched_signatures_vcf = filtered_matched_control_signatures
    Array[File]? db_signatures_vcf = filtered_db_control_signatures
    # Coverage stats
    File coverage_bed = ExtractCoverageOverVcfFiles.coverage_bed
    File coverage_bed_index = ExtractCoverageOverVcfFiles.coverage_bed_index

    File? md5_checksums_json = MergeMd5sToJson.md5_json
  }
}
