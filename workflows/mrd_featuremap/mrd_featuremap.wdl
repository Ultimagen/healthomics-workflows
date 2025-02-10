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
        String pipeline_version = "1.17.1" # !UnusedDeclaration
        String base_file_name
        # Outputs from single_read_snv.wdl (cfDNA sample)
        File cfdna_featuremap
        File cfdna_featuremap_index
        File featuremap_df_file
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
        #@wv prefix(cfdna_cram_bam_index) == cfdna_cram_bam
        #@wv suffix(cfdna_cram_bam) in {'.bam', '.cram'}
        #@wv suffix(cfdna_cram_bam_index) in {'.bai', '.crai'}
        #@wv suffix(include_regions) <= {'.bed', '.gz'}
        #@wv suffix(exclude_regions) <= {'.bed', '.gz', '.vcf'}
        #@wv suffix(snv_database) in {'.gz'}
        #@wv suffix(prefix(snv_database)) in {'.vcf'}
        #@wv suffix(featuremap_df_file) in {'.parquet'}
    }

  meta {
      description : "The UG pipeline for tumor informed MRD measures the tumor fraction in cfDNA from the presence of tumor-specific SNVs. The input data is generally 3 aligned cram files:\n\n- cfDNA (plasma)\n\n- Tumor tissue (FFPE / FF)\n\n- Normal tissue (buffy coat / PBMCs)\n\n\n\nIt is possible to provide the cfDNA cram file only with an existing somatic vcf file.\n\n\n\nThe analysis of this MRD data is composed of three parts:\n\n1. Tumor signature mutation calling, where the tumor and normal tissues are used for finding the tumor somatic mutations signature with somatic variant calling (by default UG Somatic DeepVariant, though these can be provided from other callers).\n\n2. Single Read SNV pipeline, where all the SNV candidates compared to the reference genome are extracted from the cfDNA cram file to a FeatureMap vcf, annotated and assigned a quality score (SNVQ).\n\n3. Intersection and MRD data analysis, where the FeatureMap and signature are intersected and filtered, then reads supporting the tumor mutations are counted and a tumor fraction is measured. Control signatures can be added to estimate the background noise, e.g. from other cohort patients, and in addition control signatures are generated from a somatic mutation database. \n\n\n\nThis pipeline describes step #3, the intersection and MRD data analysis, once #1 and #2 are completed."
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
          "VariantCallingEvaluation.Globals.glob"
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
          help: "Mapping quality threshold for reads to be included in the coverage analysis, corresponding to the value used to filter the srsnv output",
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
      tumor_fraction_h5: {
          help: "HDF5 file of the tumor fraction and other results of the MRD analysis",
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
      coverage_per_locus_csv: {
          help: "CSV file of the coverage per locus across all signatures",
          type: "File",
          category: "output"
      }
      sample_summary: {
          help: "Summary file of the coverage per sample",
          type: "File",
          category: "output"
      }
      sample_statistics: {
          help: "Statistics file of the coverage per sample",
          type: "File",
          category: "output"
      }
      sample_interval_summary: {
          help: "Summary file of the coverage per interval",
          type: "File",
          category: "output"
      }
      sample_interval_statistics: {
          help: "Statistics file of the coverage per interval",
          type: "File",
          category: "output"
      }
      sample_cumulative_coverage_proportions: {
          help: "Cumulative coverage proportions file of the coverage per sample",
          type: "File",
          category: "output"
      }
      sample_cumulative_coverage_counts: {
          help: "Cumulative coverage counts file of the coverage per sample",
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
  File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

  # Part 1 - Filter signatures
  if (defined_external_matched_signatures) {
    Array[File] all_matched_signatures = select_first([external_matched_signatures,])
    # Filter the matched signatures over the genomic regions bed file + apply filters
    scatter (i in range(length(all_matched_signatures))) {
      call UGGeneralTasks.FilterVcfWithBcftools as FilterMatched {
        input:
          input_vcf = all_matched_signatures[i],
          docker = global.bcftools_docker,
          bcftools_extra_args = bcftools_extra_args,  
          exclude_regions = exclude_regions,
          include_regions = include_regions,
          preemptible_tries = preemptibles,
          monitoring_script = monitoring_script, #!FileCoercion
      }
    }
  }  
  
  # for the external and db controls, exclude the regions from the matched signatures
  Array[Array[File]] control_exclude_regions_array = select_all([exclude_regions, external_matched_signatures])
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
  call UGMrdTasks.ExtractCoverageOverVcfFiles as ExtractCoverageOverVcfFiles{
    input:
      vcf_files = all_vcf_files,
      input_cram_bam = cfdna_cram_bam,
      input_cram_bam_index = cfdna_cram_bam_index,
      base_file_name = base_file_name,
      mapping_quality_threshold = mapping_quality_threshold,
      references = references,
      docker = global.ugbio_mrd_docker,
      memory_gb = memory_extract_coverage,
      cpus = 2,
      preemptibles = preemptibles,
      monitoring_script = monitoring_script  #!FileCoercion
  }

  # Part 3 - Intersect FeatureMap with signatures
  Float featuremap_size = size(cfdna_featuremap, "GB")
  call UGMrdTasks.FeatureMapIntersectWithSignatures as FeatureMapIntersectWithSignatures{
    input:
      featuremap = cfdna_featuremap,
      featuremap_index = cfdna_featuremap_index,
      matched_signatures = FilterMatched.output_vcf,
      control_signatures = FilterControlSignatures.output_vcf,
      db_signatures = GenerateControlSignaturesFromDatabase.db_signatures,
      docker = global.ugbio_mrd_docker,
      disk_size = 2 * featuremap_size + 30,
      memory_gb = 8,
      cpus = 2,
      monitoring_script = monitoring_script  #!FileCoercion
  }

  # Part 4 - Integrate all the processed data in the MRD data analysis
  call UGMrdTasks.MrdDataAnalysis as MrdDataAnalysis{
    input:
      intersected_featuremaps_parquet = FeatureMapIntersectWithSignatures.intersected_featuremaps_parquet,
      matched_signatures_vcf = filtered_matched_control_signatures,
      control_signatures_vcf = filtered_external_control_signatures,
      db_signatures_vcf = filtered_db_control_signatures,
      coverage_csv = ExtractCoverageOverVcfFiles.coverage_per_locus_csv,
      mrd_analysis_params = mrd_analysis_params,
      basename = base_file_name,
      featuremap_df_file = featuremap_df_file,
      docker = global.ugbio_mrd_docker,
      disk_size = 3 * featuremap_size + 30,
      memory_gb = memory_extract_coverage,
      cpus = 4,
      monitoring_script = monitoring_script  #!FileCoercion
  }

  output {
    # MRD analysis results
    File features_dataframe = MrdDataAnalysis.features
    File signatures_dataframe = MrdDataAnalysis.signatures
    File report_html = MrdDataAnalysis.mrd_analysis_html
    File tumor_fraction_h5 = MrdDataAnalysis.tumor_fraction_h5
    # Intersected featuremaps
    Array[File] intersected_featuremaps_parquet = FeatureMapIntersectWithSignatures.intersected_featuremaps_parquet
    Array[File] intersected_featuremaps = FeatureMapIntersectWithSignatures.intersected_featuremaps
    Array[File] intersected_featuremaps_indices = FeatureMapIntersectWithSignatures.intersected_featuremaps_indices
    # filtered signatures
    Array[File]? control_signatures_vcf = filtered_external_control_signatures
    Array[File]? matched_signatures_vcf = filtered_matched_control_signatures
    Array[File]? db_signatures_vcf = filtered_db_control_signatures
    # Coverage stats
    File coverage_per_locus_csv = ExtractCoverageOverVcfFiles.coverage_per_locus_csv
    File sample_summary = ExtractCoverageOverVcfFiles.sample_summary
    File sample_statistics = ExtractCoverageOverVcfFiles.sample_statistics
    File sample_interval_summary = ExtractCoverageOverVcfFiles.sample_interval_summary
    File sample_interval_statistics = ExtractCoverageOverVcfFiles.sample_interval_statistics
    File sample_cumulative_coverage_proportions = ExtractCoverageOverVcfFiles.sample_cumulative_coverage_proportions
    File sample_cumulative_coverage_counts = ExtractCoverageOverVcfFiles.sample_cumulative_coverage_counts    
  }
}
