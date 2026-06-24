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
    File ctdna_vaf_h5 = "~{basename}.ctdna_vaf.h5"
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
    File signature
    File signature_index
    String signature_type
    String docker
    Float disk_size
    Int memory_gb
    Int cpus
    File monitoring_script
  }
  String featuremap_base = sub(basename(featuremap), "\\..*", "")
  String signature_base = sub(basename(signature), "\\..*", "")
  String output_vcf_basename = featuremap_base + "." + signature_base + "." + signature_type + ".intersection"
  command <<<
    set -xeuo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    featuremap_base=$(basename ~{featuremap})
    signature_base=$(basename ~{signature})
    output_vcf="${featuremap_base%%.*}.${signature_base%%.*}.~{signature_type}.intersection.vcf.gz"

    echo "******** making sure files are in the same directory as their index ********"
    workdir="${featuremap_base%%.*}.${signature_base%%.*}"
    sig_dir="${workdir}/signatures"
    fm_dir="${workdir}/featuremap"
    mkdir -p "${sig_dir}" "${fm_dir}"
    ln -s ~{signature} "${sig_dir}/${signature_base}.vcf.gz"
    ln -s ~{signature_index} "${sig_dir}/${signature_base}.vcf.gz.tbi"
    ln -s ~{featuremap} "${fm_dir}/${featuremap_base}.vcf.gz"
    ln -s ~{featuremap_index} "${fm_dir}/${featuremap_base}.vcf.gz.tbi"

    echo "******** Run intersection ********"
    bcftools isec -n=2 -w1 "${fm_dir}/${featuremap_base}.vcf.gz" "${sig_dir}/${signature_base}.vcf.gz" -Oz -o "$output_vcf" --threads ~{cpus} --write-index=tbi

    echo "******** Converting to dataframe ********"
    number_of_lines=$(bcftools view "$output_vcf" -H | wc -l)
    if [[ $number_of_lines -eq 0 ]]; then
      echo "Skipping empty VCF: $output_vcf"
      touch "${output_vcf%.vcf.gz}.parquet"
    else
      featuremap_to_dataframe --in "$output_vcf" --out "${output_vcf%.vcf.gz}.parquet" --jobs ~{cpus} --drop-format AD GT
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
  output {
    File monitoring_log = "monitoring.log"
    File intersected_featuremap = output_vcf_basename + ".vcf.gz"
    File intersected_featuremap_index = output_vcf_basename + ".vcf.gz.tbi"
    File intersected_featuremap_parquet = output_vcf_basename + ".parquet"
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

task FilterSignatureOnExactAltAllele {
  input {
    File signature_vcf
    File signature_vcf_index
    Array[File] exclude_vcfs
    Array[File] exclude_vcf_indices
    String docker
    Int preemptible_tries = 1
    File monitoring_script
    Int disk_size = ceil(2 * size(signature_vcf, "GB") + size(exclude_vcfs, "GB") + 10)
    Int memory_gb = 4
    Int cpus = 2
  }
  parameter_meta {
    signature_vcf: {
      help: "Input signature VCF to filter.",
      type: "File",
      category: "input_required"
    }
    signature_vcf_index: {
      help: "Respective tabix index.",
      type: "File",
      category: "input_required"
    }
    exclude_vcfs: {
      help: "Array of VCFs to exclude. For each VCF, variants matching by locus and exact alt allele will be removed from the signature.",
      type: "Array[File]",
      category: "input_required"
    }
    exclude_vcf_indices: {
      help: "Respective tabix indices for the exclude VCFs (same order).",
      type: "Array[File]",
      category: "input_required"
    }
    docker: {
      help: "Docker image with bcftools.",
      type: "String",
      category: "input_required"
    }
    preemptible_tries: {
      help: "Number of preemption retries.",
      type: "Int",
      category: "input_optional"
    }
    monitoring_script: {
      help: "UG resource monitoring script.",
      type: "File",
      category: "input_required"
    }
    disk_size: {
      help: "Disk size in GB. Default is calculated from input file sizes.",
      type: "Int",
      category: "input_optional"
    }
    memory_gb: {
      help: "Memory in GB. Default is 4.",
      type: "Int",
      category: "input_optional"
    }
    cpus: {
      help: "Number of CPUs. Default is 2.",
      type: "Int",
      category: "input_optional"
    }
  }
  String output_basename = basename(signature_vcf, ".vcf.gz") + ".exact_alt_filtered"

  command <<<
    set -xeo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    EXCLUDE_VCFS=(~{sep=" " exclude_vcfs})
    EXCLUDE_IDXS=(~{sep=" " exclude_vcf_indices})

    ln -s ~{signature_vcf} sig.vcf.gz
    ln -s ~{signature_vcf_index} sig.vcf.gz.tbi

    EXCLUDE_ARGS=()
    for i in "${!EXCLUDE_VCFS[@]}"; do
      ln -s "${EXCLUDE_VCFS[$i]}" "exclude_${i}.vcf.gz"
      ln -s "${EXCLUDE_IDXS[$i]}" "exclude_${i}.vcf.gz.tbi"
      EXCLUDE_ARGS+=("exclude_${i}.vcf.gz")
    done

    bcftools isec -C -w1 \
      sig.vcf.gz \
      "${EXCLUDE_ARGS[@]}" \
      --threads ~{cpus} \
      --collapse some \
      -Oz -o ~{output_basename}.vcf.gz

    bcftools index -t ~{output_basename}.vcf.gz
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
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}
