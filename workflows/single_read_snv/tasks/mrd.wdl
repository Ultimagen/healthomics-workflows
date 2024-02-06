version 1.0

import "structs.wdl" as Structs

task MrdDataAnalysis {
  input {
    Array[File] intersected_featuremaps_parquet
    Array[File]? matched_signatures_vcf
    Array[File]? control_signatures_vcf
    Array[File]? db_signatures_vcf
    File coverage_csv
    MrdAnalysisParams mrd_analysis_params
    String basename
    File srsnv_test_X
    File srsnv_test_y
    File srsnv_qual_test
    String docker
    Float disk_size
    Int memory_gb
    Int cpus
    File monitoring_script
  }
  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    echo "********** Aggregating data to neat dataframes **********"
    python /VariantCalling/ugvc prepare_data_from_mrd_pipeline \
      --intersected-featuremaps ~{sep=" " intersected_featuremaps_parquet} \
      --coverage-csv ~{coverage_csv} \
      ~{true="--matched-signatures-vcf " false="" defined(matched_signatures_vcf)}~{sep=" " matched_signatures_vcf} \
      ~{true="--control-signatures-vcf " false="" defined(control_signatures_vcf)}~{sep=" " control_signatures_vcf} \
      ~{true="--db-control-signatures-vcf " false="" defined(db_signatures_vcf)}~{sep=" " db_signatures_vcf} \
      --output-dir "$PWD" \
      --output-basename "~{basename}" \


    echo "********** Creating MRD analysis notebook **********"
    papermill /VariantCalling/ugvc/reports/mrd_automatic_data_analysis.ipynb ~{basename}.mrd_data_analysis.ipynb \
      -p features_file_parquet "~{basename}.features.parquet" \
      -p signatures_file_parquet "~{basename}.signatures.parquet" \
      -p signature_filter_query "~{mrd_analysis_params.signature_filter_query}" \
      -p read_filter_query "~{mrd_analysis_params.read_filter_query}" \
      -p X_test_file "~{srsnv_test_X}" \
      -p y_test_file "~{srsnv_test_y}" \
      -p qual_test_file "~{srsnv_qual_test}" \
      -p output_dir "$PWD" \
      -p basename "~{basename}" \
      -k python3

    echo "********** Creating html version **********"
    jupyter nbconvert ~{basename}.mrd_data_analysis.ipynb --output ~{basename}.mrd_data_analysis.html \
      --to html --template classic --no-input
  >>>
  runtime {
    preemptible: 0
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File monitoring_log = "monitoring.log"
    File features = "~{basename}.features.parquet"
    File signatures = "~{basename}.signatures.parquet"
    File mrd_analysis_notebook = "~{basename}.mrd_data_analysis.ipynb"
    File mrd_analysis_html = "~{basename}.mrd_data_analysis.html"
    File tumor_fraction_h5 = "~{basename}.tumor_fraction.h5"
  }
}

task SelectMRDQualitySNVs {
  input {
    String docker
    String? base_file_name
    File monitoring_script
    File input_vcf
    String? jexl_variant_selectors_for_mrd_snvs
    String output_suffix = ".mrd_quality_snvs.vcf.gz"
    Array[File]? exclude_regions
    Array[File]? exclude_regions_index
    Array[File]? include_regions
    Array[File]? include_regions_index
    Int preemptible_tries
    Float disk_size
    Int memory_gb
    Int cpus
  }
  String base_file_name_ = select_first([base_file_name, basename(input_vcf)])
  String output_file_name = base_file_name_ + output_suffix
  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    bcftools index -t ~{input_vcf}
    gatk \
        SelectVariants -V ~{input_vcf} \
        -O ~{output_file_name} \
        -select "~{jexl_variant_selectors_for_mrd_snvs}" \
        ~{true="-L " false="" defined(include_regions)}~{sep=" -L " include_regions} \
        ~{true="-XL " false="" defined(exclude_regions)}~{sep=" -XL " exclude_regions} \
        --interval-set-rule INTERSECTION \
        --create-output-variant-index \
        --java-options "-Xmx~{memory_gb-2}G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker
    cpu: cpus
	  memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    maxRetries: 1
  }

  output {
    File monitoring_log = "monitoring.log"
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
  }
}

task GenerateControlSignaturesFromDatabase {
  input {
    File signature_file  # number of variants per motif will be used for db sampling
    File snv_database  # unfiltered SNV database file
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int n_synthetic_signatures
    String docker
    Float disk_size
    Int memory_gb
    Int cpus
    File monitoring_script
  }
  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    echo "********** Generating control signatures from database **********"
    python /VariantCalling/ugvc generate_synthetic_signatures \
      --signature_vcf ~{signature_file} \
      --db_vcf ~{snv_database} \
      --n_synthetic_signatures ~{n_synthetic_signatures} \
      --ref_fasta ~{ref_fasta} \
      --output_dir ./

    echo "********** Converting control signatures to dataframe **********"
    find syn*.vcf.gz | \
      xargs -P~{cpus} -I% sh -c "python /VariantCalling/ugvc featuremap_to_dataframe -i %"
    echo "********** DONE **********"

  >>>
  runtime {
    preemptible: 0
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output {
    File monitoring_log = "monitoring.log"
    Array[File] db_signatures = glob("syn*.vcf.gz")
    Array[File] db_signatures_indices = glob("syn*.vcf.gz.tbi")
    Array[File] db_signatures_parquet = glob("*.parquet")
  }
}


task FeatureMapIntersectWithSignatures {
  input {
    File featuremap
    File featuremap_index
    Array[File]? matched_signatures
    Array[File]? control_signatures
    Array[File]? db_signatures
    String docker
    Float disk_size
    Int memory_gb
    Int cpus
    File monitoring_script
  }
  Boolean is_defined_matched_signatures = defined(matched_signatures)
  Boolean is_defined_control_signatures = defined(control_signatures)
  Boolean is_defined_db_signatures = defined(db_signatures)
  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    . ~/.bashrc
    conda activate genomics.py3


    # run intersections
    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_matched_signatures} ]]
    then
        echo "******** Processing matched signatures ********"
        echo "******** 1/2 Run intersections ********"
        echo ~{sep=" " matched_signatures} | tr " " $"\n" | sort | uniq | \
          xargs -P~{cpus} -I% sh -c \
          "python /VariantCalling/ugvc intersect_featuremap_with_signature \
          --featuremap ~{featuremap} \
          --signature % \
          --signature_type matched"
      echo "******** 2/2 Converting to dataframes ********"
      # the next command is not multiprocessed because it uses pyfaidx which is not multiprocessing-safe
      find *.matched.intersection.vcf.gz | \
        xargs -P1 -I% sh -c "python /VariantCalling/ugvc featuremap_to_dataframe -i %"
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_control_signatures} ]]
    then
      echo "******** Processing control signatures ********"
      echo "******** 1/2 Run intersections ********"
      echo ~{sep=" " control_signatures} | tr " " $"\n" | sort | uniq | \
        xargs -P~{cpus} -I% sh -c \
        "python /VariantCalling/ugvc intersect_featuremap_with_signature \
        --featuremap ~{featuremap} \
        --signature % \
        --signature_type control"
      echo "******** 2/2 Converting to dataframes ********"
      find *.control.intersection.vcf.gz | \
        xargs -P1 -I% sh -c "python /VariantCalling/ugvc featuremap_to_dataframe -i %"
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_db_signatures} ]]
    then
      echo "******** Processing db control signatures ********"
      echo "******** 1/2 Run intersections ********"
      echo ~{sep=" " db_signatures} | tr " " $"\n" | sort | uniq | \
        xargs -P~{cpus} -I% sh -c \
        "python /VariantCalling/ugvc intersect_featuremap_with_signature \
        --featuremap ~{featuremap} \
        --signature % \
        --signature_type db_control"
      echo "******** 2/2 Converting to dataframes ********"
      find *.db_control.intersection.vcf.gz | \
        xargs -P1 -I% sh -c "python /VariantCalling/ugvc featuremap_to_dataframe -i %"
    fi

    echo "******** DONE ********"
  >>>
  runtime {
    preemptible: 0
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File monitoring_log = "monitoring.log"
    Array[File] intersected_featuremaps = glob("*.intersection.vcf.gz")
    Array[File] intersected_featuremaps_indices = glob("*.intersection.vcf.gz.tbi")
    Array[File] intersected_featuremaps_parquet = glob("*.parquet")
  }
}


task TrainSnvQualityRecalibrationModel {
  input {
    String basename
    File sorter_json_stats_file
    File hom_snv_featuremap
    File hom_snv_featuremap_index
    File singleton_snv_featuremap
    File singleton_snv_featuremap_index
    Int train_set_size
    Int test_set_size
    File hom_snv_regions_bed
    File single_substitution_regions_bed
    Array[String] numerical_features
    Array[String] categorical_features
    Array[String]? balanced_sampling_info_fields
    String? pre_filter
    Int random_seed
    String? balanced_strand_adapter_version
    References references
    String flow_order
    File monitoring_script
    String docker
    Int preemptible_tries
    Int cpus = 1
    Float memory_gb = 16
    Float disk_size = ceil(size(hom_snv_featuremap, "GB") + size(singleton_snv_featuremap, "GB")+ size(references.ref_fasta, "GB") + 30)
  }
  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    echo "***************************** Running Train Snv Quality Recalibration Model *****************************"
    start_task=$(date +%s)

    python /VariantCalling/ugvc srsnv_training \
    --hom_snv_featuremap ~{hom_snv_featuremap} \
    --single_substitution_featuremap ~{singleton_snv_featuremap} \
    --flow_order "~{flow_order}" \
    --reference_fasta ~{references.ref_fasta}  \
    --cram_stats_file ~{sorter_json_stats_file} \
    --train_set_size ~{train_set_size} \
    --test_set_size ~{test_set_size} \
    ~{true="--balanced_sampling_info_fields " false="" defined(balanced_sampling_info_fields)}~{sep=" " balanced_sampling_info_fields} \
    ~{true="--balanced_strand_adapter_version " false="" defined(balanced_strand_adapter_version)}~{balanced_strand_adapter_version} \
    --hom_snv_regions ~{hom_snv_regions_bed} \
    --single_sub_regions ~{single_substitution_regions_bed} \
    --numerical_features ~{sep=" " numerical_features} \
    --categorical_features ~{sep=" " categorical_features} \
    ~{true="--pre_filter " false="" defined(pre_filter)}"~{pre_filter}" \
    --random_seed ~{random_seed} \
    --output "$PWD" --basename ~{basename}
    #TODO add LoD filters as an optional input (--lod_filters lod_filters.json)

    end=$(date +%s)
    mins_elapsed=$(( (end - start_task) / 60 ))
    secs_elapsed=$(( (end - start_task) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"

  >>>
  runtime {
    preemptible: "~{preemptible_tries}"
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File monitoring_log = "monitoring.log"
    File model_file = "~{basename}.model.joblib"
    File X_test_file = "~{basename}.X_test.parquet"
    File y_test_file = "~{basename}.y_test.parquet"
    File qual_test_file = "~{basename}.qual_test.parquet"    
    File X_train_file = "~{basename}.X_train.parquet"
    File y_train_file = "~{basename}.y_train.parquet"
    File params_file = "~{basename}.params.json"
    File test_set_mrd_simulation_dataframe = "~{basename}.test.df_mrd_simulation.parquet"
    File train_set_mrd_simulation_dataframe = "~{basename}.train.df_mrd_simulation.parquet"
    File test_set_statistics_h5 = "~{basename}.test.statistics.h5"
    File train_set_statistics_h5 = "~{basename}.train.statistics.h5"
    File test_set_statistics_json = "~{basename}.test.statistics.json"
    File train_set_statistics_json = "~{basename}.train.statistics.json"
    File train_report_file_notebook = "~{basename}.train_report.ipynb"
    File train_report_file_html = "~{basename}.train_report.html"
    File test_report_file_notebook = "~{basename}.test_report.ipynb"
    File test_report_file_html = "~{basename}.test_report.html"
    Array[File] training_report_plot_png = glob("~{basename}*.png")
  }
}


task InferenceSnvQualityRecalibrationModel {
  input {
    File model_file
    File params_file
    File featuremap
    File featuremap_index
    File X_train_file
    File test_set_mrd_simulation_dataframe_file
    String output_file

    File monitoring_script
    String docker
    Int preemptible_tries
    Int cpus = 8
    Float memory_gb = 8
    Float disk_size = ceil(3 * size(featuremap, "GB") + 30)
  }

  command <<<
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    set -eo pipefail
    . ~/.bashrc
    conda activate genomics.py3

    echo "***************************** Running Inference Snv Quality Recalibration Model *****************************"
    start_task=$(date +%s)

    python /VariantCalling/ugvc srsnv_inference \
    --featuremap_path "~{featuremap}" \
    --params_path "~{params_file}" \
    --model_path "~{model_file}" \
    --X_train_path "~{X_train_file}" \
    --output_path "~{output_file}" \
    --test_set_mrd_simulation_dataframe_file "~{test_set_mrd_simulation_dataframe_file}" \
    --process_number ~{cpus}

    end=$(date +%s)
    mins_elapsed=$(( (end - start_task) / 60 ))
    secs_elapsed=$(( (end - start_task) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"
  >>>
  runtime {
    preemptible: "~{preemptible_tries}"
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output {
    File monitoring_log = "monitoring.log"
    File featuremap_with_ml_qual = "~{output_file}"
    File featuremap_with_ml_qual_index = "~{output_file}.tbi"
  }
}


task CreateHomSnvFeatureMap {
  input {
      File featuremap
      File featuremap_index
      File? sorter_json_stats_file
      String base_file_name
      Float min_af
      Int min_coverage
      Float disk_size = ceil(size(featuremap, "GB") * 2 + 20)
      Float memory_gb = 2
      Int cpus = 1
      String docker
      Int preemptibles
      File monitoring_script
    }

  String hom_snv_vcf = "~{base_file_name}.hom_snv.vcf.gz"

  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    source ~/.bashrc
    conda activate genomics.py3

    python /VariantCalling/ugvc create_hom_snv_featuremap \
      --featuremap ~{featuremap} \
      ~{true="--sorter_stats_json " false="" defined(sorter_json_stats_file)} ~{sorter_json_stats_file} \
      --hom_snv_featuremap ~{hom_snv_vcf} \
      --requested_min_coverage ~{min_coverage} \
      --min_af ~{min_af} 

  >>>

  runtime {
    preemptible: "~{preemptibles}"
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File monitoring_log = "monitoring.log"
    File hom_snv_featuremap = "~{hom_snv_vcf}"
    File hom_snv_featuremap_index = "~{hom_snv_vcf}.tbi"
  }

}
 

task BedIntersectAndExclude {
    input {
      Array[File] include_regions
      Array[File]? exclude_regions
      String output_basename
      String docker
      Float disk_size = 5 + 10 * ceil(size(include_regions, "GB") + if defined(exclude_regions) then size(select_first([exclude_regions]), "GB") else 0)
      Float memory_gb = 2
      Int cpus = 1
      Int preemptibles
      File monitoring_script
    }
    Boolean exclude_regions_defined = length(select_first([exclude_regions, []])) > 0
    command <<<
      set -eo pipefail
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      source ~/.bashrc
      conda activate genomics.py3

      python /VariantCalling/ugvc intersect_bed_regions \
        --include-regions ~{sep=" " include_regions} \
        ~{true="--exclude-regions " false="" exclude_regions_defined}~{sep=" " exclude_regions} \
        --output-bed "~{output_basename}.bed"

    >>>

    output {
      File merged_bed = "~{output_basename}.bed"
    }

    runtime {
      preemptible: "~{preemptibles}"
      cpu: "~{cpus}"
      memory: "~{memory_gb} GB"
      disks: "local-disk " + ceil(disk_size) + " HDD"
      docker: docker
    }
}


task ExtractCoverageOverVcfFiles {
    # Task: ExtractCoverageOverVcfFiles
    # Description:
    #     This task extracts coverage metrics from a given CRAM file over specified loci provided in VCFs.
    #     The loci from VCFs are combined into a BED file, merged, and then processed using GATK's DepthOfCoverage.
    #     The main output is a CSV file containing per-locus coverage metrics.
    #     The CRAM file is not localized.
    
    # Inputs:
    #     vcf_files: Array of VCF files containing loci for coverage extraction.
    #     input_cram_bam: Input CRAM/BAM file containing read alignments.
    #     input_cram_bam_index: Index of respective CRAM/BAM file.
    #     base_file_name: Base string used to name output files.
    #     mapping_quality_threshold: Minimum mapping quality threshold for reads to be considered.
    #     references: Reference related files - fasta, index, and dictionary.
    #     docker: Docker image to use for task execution.
    #     memory_gb: Amount of memory to allocate for the task.
    #     cpus: Number of CPU cores to allocate for the task.
    #     preemptibles: Number of preemption retries.
    #     monitoring_script: Path to a script to monitor task execution.

    # Outputs:
    #     coverage_per_locus: CSV file containing per-locus coverage metrics.
    #     coverage_summary: Summary of coverage metrics.
    #     coverage_statistics: Statistics related to coverage.
    #     interval_summary: Coverage metrics summarized over intervals.
    #     interval_statistics: Statistics related to interval coverage.
    input {
      Array[File] vcf_files
      File input_cram_bam
      File input_cram_bam_index
      String base_file_name
      Int mapping_quality_threshold = 60  
      References references
      String docker
      Float memory_gb = 2
      Int cpus = 1
      Int preemptibles
      File monitoring_script
    }

    parameter_meta {
      input_cram_bam: {
        localization_optional: true
      }
    }

    Int vcf_files_size = ceil(size(vcf_files, "GB"))
    Int reference_size = ceil(size(references.ref_fasta, "GB"))
    Int disk_size = ceil((2*vcf_files_size) + reference_size) + 30  # Sum of VCF and reference sizes, plus 10GB overhead
    String output_format = "CSV"

    command <<<
      set -eo pipefail
      bash ~{monitoring_script} | tee monitoring.log >&2 &
      source ~/.bashrc
      conda activate genomics.py3

      echo "Combining all the VCF loci into one BED file..."

      for vcf in ~{sep=" " vcf_files}; do
          zcat $vcf | grep -v "^#" | awk '{print $1"\t"($2-1)"\t"$2}' >> combined_loci.bed
      done

      echo "Sorting and merging the combined BED..."
      sort -k1,1 -k2,2n combined_loci.bed | bedtools merge > merged_loci.bed

      echo "Extracting coverage from CRAM for the specified loci..."
      gatk \
        DepthOfCoverage \
        -R ~{references.ref_fasta} \
        -O ~{base_file_name}.coverage \
        -I ~{input_cram_bam} \
        --intervals merged_loci.bed \
        --read-filter MappingQualityReadFilter --minimum-mapping-quality ~{mapping_quality_threshold} \
        --output-format "~{output_format}"

      # The default output does not contain a prefix
      mv "~{base_file_name}.coverage" "~{base_file_name}.coverage.csv"

      echo "Editing csv output..."
      # The default output format is chr1:1000 and we convert it to chr1,1000 
      sed -i -e 's/Locus/Chrom,Pos/' -e 's/:/,/g' "~{base_file_name}.coverage.csv"

      echo "Coverage extraction completed."
    >>>

    output {
      File coverage_per_locus_csv = "~{base_file_name}.coverage.csv"
      File sample_summary = "~{base_file_name}.coverage.sample_summary"
      File sample_statistics = "~{base_file_name}.coverage.sample_statistics"
      File sample_interval_summary = "~{base_file_name}.coverage.sample_interval_summary"
      File sample_interval_statistics = "~{base_file_name}.coverage.sample_interval_statistics"
      File sample_cumulative_coverage_proportions = "~{base_file_name}.coverage.sample_cumulative_coverage_proportions"
      File sample_cumulative_coverage_counts = "~{base_file_name}.coverage.sample_cumulative_coverage_counts"
    }

    runtime {
      preemptible: "~{preemptibles}"
      cpu: "~{cpus}"
      memory: "~{memory_gb} GB"
      disks: "local-disk " + ceil(disk_size) + " HDD"
      docker: docker
    }
}
