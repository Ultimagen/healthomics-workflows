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
    File srsnv_metadata_json
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
      ~{true="--tumor-sample " false="" defined(mrd_analysis_params.tumor_sample)}~{mrd_analysis_params.tumor_sample} \
      --featuremap-file "~{featuremap_df_file}" \
      --srsnv-metadata-json "~{srsnv_metadata_json}"

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
  }
}

task FeatureMapIntersectWithSignatures {
  input {
    File featuremap
    File featuremap_index
    Array[File]? matched_signatures
    Array[File]? matched_signatures_indexes
    Array[File]? control_signatures
    Array[File]? control_signatures_indexes
    Array[File]? db_signatures
    Array[File]? matched_signatures_indices
    Array[File]? control_signatures_indices
    Array[File]? db_signatures_indices
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
        echo "******** Run intersections ********"
        for signature in ~{sep=" " matched_signatures}
        do
          featuremap_base=$(basename ~{featuremap})
          signature_base=$(basename "$signature")
          signature_type="matched"
          output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
          bcftools isec -n=2 -w1 ~{featuremap} "$signature" -Oz -o "$output_vcf" && bcftools index -t "$output_vcf"
        done
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_control_signatures} ]]
    then
        echo "******** Processing control signatures ********"
        echo "******** Run intersections ********"
        for signature in ~{sep=" " control_signatures}
        do
          featuremap_base=$(basename ~{featuremap})
          signature_base=$(basename "$signature")
          signature_type="control"
          output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
          bcftools isec -n=2 -w1 ~{featuremap} "$signature" -Oz -o "$output_vcf" && bcftools index -t "$output_vcf"
        done
    fi

    if [[ -z ~{default='"skip"' true='""' false='"skip"' is_defined_db_signatures} ]]
    then
        echo "******** Processing db control signatures ********"
        echo "******** Run intersections ********"
        for signature in ~{sep=" " db_signatures}
        do
          featuremap_base=$(basename ~{featuremap})
          signature_base=$(basename "$signature")
          signature_type="db_control"
          output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.${signature_type}.intersection.vcf.gz"
          bcftools isec -n=2 -w1 ~{featuremap} "$signature" -Oz -o "$output_vcf" && bcftools index -t "$output_vcf"
        done
    fi

    echo "******** 2/2 Converting to dataframes ********"
    find *.intersection.vcf.gz | \
    for signature_vcf in $(cat) 
      do
        # check if vcf is empty
        number_of_lines=$(bcftools view "$signature_vcf" -H | wc -l)
        if [[ $number_of_lines -eq 0 ]]; then
          echo "Skipping empty VCF: $signature_vcf"
          continue
        fi
        output_parquet="${signature_vcf%.vcf.gz}.parquet"
        echo "Converting $signature_vcf to dataframe: $output_parquet"
        featuremap_to_dataframe --in "$signature_vcf" --out "$output_parquet" --jobs ~{cpus_featuremap_to_dataframe} --drop-format AD GT
      done

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
      set -xe
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
      Int mapping_quality_threshold = 0 
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
    File ref_fai
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
    # extract genome file from ref_fasta_index (.fai)
    cut -f1,2 ~{ref_fai} > genome.txt
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
