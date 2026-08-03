version 1.0

import "structs.wdl" as Structs

task MrdDataAnalysis {
  input {
    Array[File] intersected_featuremaps_parquet
    File? matched_signature_vcf
    Array[File]? control_signatures_vcf
    Array[File]? db_signatures_vcf
    File coverage_bed
    MrdAnalysisParams mrd_analysis_params
    String basename
    File featuremap_df_file
    File srsnv_metadata_json
    File? filter_funnel_json
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
      ~{true="--matched-signature-vcf " false="" defined(matched_signature_vcf)}~{matched_signature_vcf} \
      ~{true="--control-signatures-vcf " false="" defined(control_signatures_vcf)}~{sep=" " control_signatures_vcf} \
      ~{true="--db-control-signatures-vcf " false="" defined(db_signatures_vcf)}~{sep=" " db_signatures_vcf} \
      --output-dir "$PWD" \
      --output-basename "~{basename}" \
      --signature-filter-query "~{mrd_analysis_params.signature_filter_query}" \
      --read-filter-query "~{mrd_analysis_params.read_filter_query}" \
      ~{true="--tumor-sample " false="" defined(mrd_analysis_params.tumor_sample)}~{mrd_analysis_params.tumor_sample} \
      --featuremap-file "~{featuremap_df_file}" \
      --srsnv-metadata-json "~{srsnv_metadata_json}" \
      ~{true="--filter-funnel-json " false="" defined(filter_funnel_json)}~{filter_funnel_json} \
      ~{true="--thresh-noise-lq-reads " false="" defined(mrd_analysis_params.thresh_noise_lq_reads)}~{mrd_analysis_params.thresh_noise_lq_reads} \
      ~{true="--thresh-multi-read-pvalue " false="" defined(mrd_analysis_params.thresh_multi_read_pvalue)}~{mrd_analysis_params.thresh_multi_read_pvalue} \
      ~{true="--alpha " false="" defined(mrd_analysis_params.mrd_detection_fpr)}~{mrd_analysis_params.mrd_detection_fpr} \
      ~{true="--lod-fpr " false="" defined(mrd_analysis_params.lod_fpr)}~{mrd_analysis_params.lod_fpr} \
      ~{true="--lod-recall " false="" defined(mrd_analysis_params.lod_recall)}~{mrd_analysis_params.lod_recall}

  >>>
  runtime {
    preemptible: 0
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File monitoring_log = "monitoring.log"
    File features = "~{basename}.features.parquet"
    File signatures = "~{basename}.signatures.parquet"
    File mrd_analysis_html = "~{basename}.mrd_analysis_report.html"
    File mrd_qc_html = "~{basename}.mrd_qc_report.html"
    File detection_result_json = "~{basename}.detection_result.json"
    File ctdna_vaf_h5 = "~{basename}.ctdna_vaf.h5"
    # Only produced when a matched signature is present (filter funnels are matched-only).
    File? output_filter_funnel_json = "~{basename}.filter_funnel.json"
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
    cpu: cpus
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
    echo "{\"intersected_reads\": $number_of_lines}" > intersection_funnel.json
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
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output {
    File monitoring_log = "monitoring.log"
    File intersected_featuremap = output_vcf_basename + ".vcf.gz"
    File intersected_featuremap_index = output_vcf_basename + ".vcf.gz.tbi"
    File intersected_featuremap_parquet = output_vcf_basename + ".parquet"
    File intersection_funnel_json = "intersection_funnel.json"
  }
}

task CollectFilterFunnel {
  input {
    Array[File] filter_funnel_jsons
    Array[File?]? exact_alt_funnel_jsons
    Array[File] intersection_funnel_jsons
    String? bcftools_extra_args
    Array[String] include_region_names = []
    Array[String] exclude_region_names = []
    Array[String] exact_alt_region_names = []
    String docker
    Int disk_size = 2
    Int memory_gb = 2
    Int cpus = 1
  }
  parameter_meta {
    filter_funnel_jsons: {
      help: "Per-signature filter funnel JSON files from FilterVcfWithBcftools.",
      type: "Array[File]",
      category: "input_required"
    }
    exact_alt_funnel_jsons: {
      help: "Per-signature exact alt allele filter count JSON files.",
      type: "Array[File?]",
      category: "input_optional"
    }
    intersection_funnel_jsons: {
      help: "Per-signature intersection count JSON files from FeatureMapIntersectWithSignatures.",
      type: "Array[File]",
      category: "input_required"
    }
    bcftools_extra_args: {
      help: "bcftools view extra args used for signature filtering; recorded as the 'After bcftools extra args' funnel step description.",
      type: "String",
      category: "input_optional"
    }
    include_region_names: {
      help: "Basenames of the include-region files; recorded as the 'After include regions' funnel step description.",
      type: "Array[String]",
      category: "input_optional"
    }
    exclude_region_names: {
      help: "Basenames of the exclude-region files; recorded as the 'After exclude regions' funnel step description.",
      type: "Array[String]",
      category: "input_optional"
    }
    exact_alt_region_names: {
      help: "Basenames of the exact-alt-allele exclude VCFs; recorded as the 'After exact alt allele filter' funnel step description.",
      type: "Array[String]",
      category: "input_optional"
    }
    docker: {
      help: "Docker image with Python 3.",
      type: "String",
      category: "input_required"
    }
    disk_size: {
      help: "Disk size in GB.",
      type: "Int",
      category: "input_optional"
    }
    memory_gb: {
      help: "Memory in GB.",
      type: "Int",
      category: "input_optional"
    }
    cpus: {
      help: "Number of CPUs.",
      type: "Int",
      category: "input_optional"
    }
  }
  Array[File] exact_alt_files_resolved = select_all(select_first([exact_alt_funnel_jsons, []]))
  # Pass free-text / list inputs through files to avoid shell-quoting issues
  # (bcftools_extra_args contains quotes, e.g. -i 'QUAL>10').
  File bcftools_extra_args_file = write_lines(select_all([bcftools_extra_args]))
  File include_region_names_file = write_lines(include_region_names)
  File exclude_region_names_file = write_lines(exclude_region_names)
  File exact_alt_region_names_file = write_lines(exact_alt_region_names)
  command <<<
    set -xeuo pipefail
    python3 <<'PYEOF'
import json
import os

funnel_files = "~{sep="," filter_funnel_jsons}".split(",")
exact_alt_files_str = "~{sep="," exact_alt_files_resolved}"
exact_alt_files = exact_alt_files_str.split(",") if exact_alt_files_str else []
intersection_files = "~{sep="," intersection_funnel_jsons}".split(",")

signatures = []
for i, f in enumerate(funnel_files):
    with open(f.strip()) as fh:
        data = json.load(fh)
    if i < len(exact_alt_files) and exact_alt_files[i].strip():
        with open(exact_alt_files[i].strip()) as fh2:
            data.update(json.load(fh2))
    if i < len(intersection_files) and intersection_files[i].strip():
        with open(intersection_files[i].strip()) as fh3:
            data.update(json.load(fh3))
    signatures.append(data)


def read_lines(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def basename_no_ext(name):
    # Strip directory and common bed/vcf(.gz) extensions -> e.g. ug_hcr.bed -> ug_hcr
    base = os.path.basename(name.strip())
    for ext in (".bed.gz", ".vcf.gz", ".bed", ".vcf", ".gz"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


bcftools_args = read_lines("~{bcftools_extra_args_file}")
include_names = [basename_no_ext(n) for n in read_lines("~{include_region_names_file}")]
exclude_names = [basename_no_ext(n) for n in read_lines("~{exclude_region_names_file}")]
exact_alt_names = [basename_no_ext(n) for n in read_lines("~{exact_alt_region_names_file}")]

# Descriptions shared by all signatures, keyed by the funnel step they annotate.
descriptions = {}
if bcftools_args:
    descriptions["After bcftools extra args"] = bcftools_args[0]
if include_names:
    descriptions["After include regions"] = ", ".join(include_names)
if exclude_names:
    descriptions["After exclude regions"] = ", ".join(exclude_names)
if exact_alt_names:
    descriptions["After exact alt allele filter"] = ", ".join(exact_alt_names)

with open("filter_funnel.json", "w") as out:
    json.dump({"signatures": signatures, "descriptions": descriptions}, out, indent=2)
PYEOF
  >>>
  runtime {
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }
  output {
    File collected_funnel_json = "filter_funnel.json"
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
      cpu: cpus
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
      cpu: cpus
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
      cpu: cpus
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
    cpu: cpus
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

    # Emit count after exact alt allele filtering
    AFTER_EXACT_ALT=$(bcftools view -H ~{output_basename}.vcf.gz | wc -l)
    echo "{\"after_exact_alt_allele_filter\": $AFTER_EXACT_ALT}" > exact_alt_funnel.json
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
    File exact_alt_funnel_json = "exact_alt_funnel.json"
  }
}
