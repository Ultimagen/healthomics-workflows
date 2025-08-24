version 1.0

import "structs.wdl" as Structs

task MrdDataAnalysis {
  input {
    Array[File] intersected_featuremaps_parquet
    Array[File]? matched_signatures_vcf
    Array[File]? control_signatures_vcf
    Array[File]? db_signatures_vcf
    File coverage_bed
    MrdAnalysisParams mrd_analysis_params
    String basename
    File featuremap_df_file
    String docker
    Float disk_size
    Int memory_gb
    Int cpus
    File monitoring_script
  }
  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    generate_report \
      --intersected-featuremaps ~{sep=" " intersected_featuremaps_parquet} \
      --coverage-bed ~{coverage_bed} \
      ~{true="--matched-signatures-vcf " false="" defined(matched_signatures_vcf)}~{sep=" " matched_signatures_vcf} \
      ~{true="--control-signatures-vcf " false="" defined(control_signatures_vcf)}~{sep=" " control_signatures_vcf} \
      ~{true="--db-control-signatures-vcf " false="" defined(db_signatures_vcf)}~{sep=" " db_signatures_vcf} \
      --output-dir "$PWD" \
      --output-basename "~{basename}" \
      --signature-filter-query "~{mrd_analysis_params.signature_filter_query}" \
      --read-filter-query "~{mrd_analysis_params.read_filter_query}" \
      --featuremap-file "~{featuremap_df_file}"

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
    set -xeuo pipefail
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
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    echo "********** Generating control signatures from database **********"
    generate_synthetic_signatures \
      --signature_vcf ~{signature_file} \
      --db_vcf ~{snv_database} \
      --n_synthetic_signatures ~{n_synthetic_signatures} \
      --ref_fasta ~{ref_fasta} \
      --output_dir ./

    echo "********** Converting control signatures to dataframe **********"
    find syn*.vcf.gz | \
      xargs -P~{cpus} -I% sh -c "featuremap_to_dataframe -i %"
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
  Int tmp_cpus_featuremap_to_dataframe = round(memory_gb / 4)
  Int cpus_featuremap_to_dataframe = if tmp_cpus_featuremap_to_dataframe < 1 then 1 else tmp_cpus_featuremap_to_dataframe
  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &


    # run intersections
    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_matched_signatures} ]]
    then
        echo "******** Processing matched signatures ********"
        echo "******** 1/2 Run intersections ********"
        echo ~{sep=" " matched_signatures} | tr " " $"\n" | sort | uniq | \
          xargs -P~{cpus} -I% sh -c \
          "intersect_featuremap_with_signature \
          --featuremap ~{featuremap} \
          --signature % \
          --signature_type matched"
      echo "******** 2/2 Converting to dataframes ********"
      # the next command is not multiprocessed because it uses pyfaidx which is not multiprocessing-safe
      find *.matched.intersection.vcf.gz | \
        xargs -P~{cpus_featuremap_to_dataframe} -I% sh -c "featuremap_to_dataframe -i %"
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_control_signatures} ]]
    then
      echo "******** Processing control signatures ********"
      echo "******** 1/2 Run intersections ********"
      echo ~{sep=" " control_signatures} | tr " " $"\n" | sort | uniq | \
        xargs -P~{cpus} -I% sh -c \
        "intersect_featuremap_with_signature \
        --featuremap ~{featuremap} \
        --signature % \
        --signature_type control"
      echo "******** 2/2 Converting to dataframes ********"
      find *.control.intersection.vcf.gz | \
        xargs -P~{cpus_featuremap_to_dataframe} -I% sh -c "featuremap_to_dataframe -i %"
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_db_signatures} ]]
    then
      echo "******** Processing db control signatures ********"
      echo "******** 1/2 Run intersections ********"
      echo ~{sep=" " db_signatures} | tr " " $"\n" | sort | uniq | \
        xargs -P~{cpus} -I% sh -c \
        "intersect_featuremap_with_signature \
        --featuremap ~{featuremap} \
        --signature % \
        --signature_type db_control"
      echo "******** 2/2 Converting to dataframes ********"
      find *.db_control.intersection.vcf.gz | \
        xargs -P~{cpus_featuremap_to_dataframe} -I% sh -c "featuremap_to_dataframe -i %"
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
    SingleReadSNVParams single_read_snv_params
    Map[String, Array[String]] categorical_features
    File hom_snv_regions_bed
    File single_substitution_regions_bed
    References references
    String flow_order
    Boolean raise_exceptions_in_report
    File monitoring_script
    String docker
    String? pipeline_version
    Int preemptible_tries
    Int cpus = 1
    Float memory_gb = 16
    Float disk_size = ceil(size(hom_snv_featuremap, "GB") + size(singleton_snv_featuremap, "GB")+ size(references.ref_fasta, "GB") + 30)
  }
  Boolean random_split = single_read_snv_params.split_folds_by=="random"
  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    start_task=$(date +%s)
    echo "***************************** Create parameters json file *****************************"
    

    # Create a json file with parameters for the model
    python <<CODE
    import json
    
    with open("~{write_json(single_read_snv_params)}", "r") as f1, open("~{write_json(categorical_features)}", "r") as f2, open("single_read_snv_params.json", "w") as f3: 
      data1 = json.load(f1)
      data2 = json.load(f2)
      if isinstance(data2, list):  # this is how write_json works in Cromwell - "[{'left': 'ref', 'right': ['A', 'C', 'G', 'T']}, ...]"
        data1["categorical_features"] = {entry["left"]: entry["right"] for entry in data2}
      else:  # this is how write_json works in Omics - "{'ref': ['A', 'C', 'G', 'T'], ...}"
        data1["categorical_features"] = data2
      if len("~{pipeline_version}") > 0:
        data1["pipeline_version"] = "~{pipeline_version}"
      data1["docker"] = "~{docker}"
      ppmSeq_adapter_version = "~{single_read_snv_params.ppmSeq_adapter_version}"
      if len(ppmSeq_adapter_version) > 0:
        data1["ppmSeq_adapter_version"] = ppmSeq_adapter_version
      json.dump(data1, f3, indent=4)
    CODE
    echo "DEBUG - single_read_snv_params.json file:"
    cat single_read_snv_params.json

    echo "***************************** Running Train Snv Quality Recalibration Model *****************************"
    srsnv_training \
    --hom_snv_featuremap ~{hom_snv_featuremap} \
    --single_substitution_featuremap ~{singleton_snv_featuremap} \
    --dataset_params_json_path single_read_snv_params.json \
    --flow_order "~{flow_order}" \
    --reference_fasta "~{references.ref_fasta}" \
    --reference_dict "~{references.ref_dict}" \
    --cram_stats_file "~{sorter_json_stats_file}" \
    --hom_snv_regions "~{hom_snv_regions_bed}" \
    --single_sub_regions "~{single_substitution_regions_bed}" \
    ~{true="--raise_exceptions_in_report" false="" raise_exceptions_in_report} \
    --output "$PWD" \
    --basename "~{basename}"
    
    end=$(date +%s)
    mins_elapsed=$(( (end - start_task) / 60 ))
    secs_elapsed=$(( (end - start_task) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"

    ls -ltr

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
    File featuremap_df_file = "~{basename}.featuremap_df.parquet"
    File params_file = "~{basename}.params.json"
    File test_set_mrd_simulation_dataframe = "~{basename}.test.df_mrd_simulation.parquet"
    File test_set_statistics_h5 = "~{basename}.single_read_snv.applicationQC.h5"
    File test_set_statistics_json = "~{basename}.test.statistics.json"
    File test_report_file_html = "~{basename}.test_report.html"
  }
}


task InferenceSnvQualityRecalibrationModel {
  input {
    File model_file
    File featuremap
    File featuremap_index
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

    echo "***************************** Running Inference Snv Quality Recalibration Model *****************************"
    start_task=$(date +%s)

    srsnv_inference \
    --featuremap_path "~{featuremap}" \
    --model_joblib_path "~{model_file}" \
    --output_path "~{output_file}" \
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
    File featuremap_with_qual = "~{output_file}"
    File featuremap_with_qual_index = "~{output_file}.tbi"
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

    create_hom_snv_featuremap \
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

      intersect_bed_regions \
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

task MergeVcfsIntoBed {
    input {
      Array[File] vcf_files
      String docker
      Float disk_size
      Float memory_gb = 2
      Int cpus = 1
      Int preemptibles
      File monitoring_script
    }
    command <<<
      set -xeo pipefail
      bash ~{monitoring_script} | tee monitoring.log >&2 &

      echo "Combining all the VCF loci into one BED file..."

      for vcf in ~{sep=" " vcf_files}; do
          zcat $vcf | grep -v "^#" | awk '{print $1"\t"($2-1)"\t"$2}' >> combined_loci.bed
      done

      echo "Sorting and merging the combined BED..."
      sort -k1,1 -k2,2n combined_loci.bed | bedtools merge > merged_loci.bed

    >>>

    output {
      File merged_loci_bed = "merged_loci.bed"
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
    #     This task extracts coverage metrics from a given CRAM file over specified loci provided in a bed file.
    #     Coverage is collected with mosdepth.
    #     The main output is a bed file containing per-locus coverage metrics.
    
    # Inputs:
    #     merged_loci_bed: A bed file containing loci for coverage extraction.
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
    #     coverage_bed: A bed file containing coverage metrics for the specified loci.
    #     coverage_bed_index: Index of the coverage bed file.
    input {
      File merged_loci_bed
      File input_cram_bam
      File input_cram_bam_index
      String base_file_name
      Int mapping_quality_threshold = 60  
      References references
      String docker
      Int memory_gb
      Int cpus = 1
      Int preemptibles
      File monitoring_script
    }

    Int merged_loci_bed_size = ceil(size(merged_loci_bed, "GB"))
    Int reference_size = ceil(size(references.ref_fasta, "GB"))
    Int input_cram_bam_size = ceil(size(input_cram_bam, "GB"))
    Int disk_size = ceil((2*merged_loci_bed_size) + reference_size + input_cram_bam_size) + 30  # Bed and reference sizes, plus 10GB overhead

    command <<<
      set -xeo pipefail
      bash ~{monitoring_script} | tee monitoring.log >&2 &

      echo "Extracting coverage from CRAM for the specified loci..."
      mosdepth --by ~{merged_loci_bed} -f ~{references.ref_fasta} -Q ~{mapping_quality_threshold} --fast-mode \
      ~{base_file_name} ~{input_cram_bam}

      echo "Coverage extraction completed."
    >>>

    output {
      File coverage_bed = "~{base_file_name}.regions.bed.gz"
      File coverage_bed_index = "~{base_file_name}.regions.bed.gz.csi"
    }

    runtime {
      preemptible: "~{preemptibles}"
      cpu: "~{cpus}"
      memory: "~{memory_gb} GB"
      disks: "local-disk " + ceil(disk_size) + " HDD"
      docker: docker
    }
}

task PileupFeatureMapInterval {
  input {
    File featuremap
    File featuremap_index
    String output_vcf_file_name
    File interval_list
    Int disk_size
    Int memory_gb
    Int cpus
    Int preemptibles
    Int? min_qual
    String? sample_name
    String? qual_agg_func
    File monitoring_script
    String docker
  }

  command <<<
    set -eo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    echo "********** Pileups featuremap into a single locus per entry **********"
    pileup_featuremap \
      --featuremap ~{featuremap} \
      --interval_list ~{interval_list} \
      --output_vcf "~{output_vcf_file_name}.vcf.gz" \
      ~{true="--min_qual " false="" defined(min_qual)}~{min_qual} \
      ~{true="--sample_name " false="" defined(sample_name)}~{sample_name} \
      ~{true="--qual_agg_func " false="" defined(qual_agg_func)}~{qual_agg_func}\
  >>>

  runtime {
    preemptible: "~{preemptibles}"
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }

  output {
    File monitoring_log = "monitoring.log"
    File pileup_file = "~{output_vcf_file_name}.vcf.gz"
    File pileup_file_index = "~{output_vcf_file_name}.vcf.gz.tbi"
  }
}

task AddAggregatedVariablesAndXgbScoreToPileupFeaturemap
{
  input {
    File featuremap_pileup_vcf
    File featuremap_pileup_vcf_index
    String output_basename
    String? filter_string
    File interval_list
    File? model 
    Int memory_gb
    Int cpus
    Int preemptibles    
    File monitoring_script
    String docker
  }
  Int disk_size = ceil(4*size(featuremap_pileup_vcf,"GB"))
  Boolean is_defined_filter_string = defined(filter_string)
  command <<<
  set -eo pipefail
  bash ~{monitoring_script} | tee monitoring.log >&2 &
  
  add_aggregate_params_and_xgb_score_to_pileup_featuremap \
    -f ~{featuremap_pileup_vcf} \
    ~{true="-filter_string " false="" defined(filter_string)}~{filter_string} \
    -o ~{output_basename} \
    -i ~{interval_list} \
    -m ~{model}
  
  >>>
  runtime {
    preemptible: "~{preemptibles}"
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }

  output {
    File monitoring_log = "monitoring.log"
    File pileup_xgb_file = "~{output_basename}.vcf.gz"
    File pileup_xgb_file_index = "~{output_basename}.vcf.gz.tbi"
  }    
}

task PadVcf {
  input {
    File input_vcf
    String docker
    Int pad_size 
    Int preemptible_tries
    File monitoring_script
    Int disk_size = ceil(3 * size(input_vcf, "GB") + 10)
    Int memory_gb = 2
    Int cpus = 1
  }
  String basename_vcf = basename(input_vcf, ".vcf.gz")
  String basename_vcf2 = basename(basename_vcf, ".vcf")
  
  command <<<
    set -xeo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    
    echo "Extract variants to bed format using bcftools query"
    # Extract CHROM, POS0, END, REF, ALT and calculate the maximum length
    bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' ~{input_vcf} | awk -F'\t' 'BEGIN{OFS="\t"}{
        s=$2;
        m=length($4);
        n=split($5,a,",");
        for(i=1;i<=n;i++) if(length(a[i])>m) m=length(a[i]);
        print $1, s, s+m
    }' | gzip > variants.bed.gz
    
    echo "Create a genome file for bedtools slop (chromosome sizes)"
    # We'll extract this from the VCF header
    bcftools view -h ~{input_vcf} | grep "^##contig" | \
      sed 's/##contig=<ID=//;s/,length=/\t/;s/>//' > genome.txt
    head genome.txt
    
    echo "Pad the bed file using bedtools slop"
    zcat variants.bed.gz | bedtools slop -i stdin -g genome.txt -b ~{pad_size} | gzip > ~{basename_vcf2}.padded.bed.gz
    
    echo "Padded bed file created successfully"
  >>>
  
  runtime {
    preemptible: preemptible_tries
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }
  
  output {
    File monitoring_log = "monitoring.log"
    File padded_bed = "~{basename_vcf2}.padded.bed.gz"
  }
}
