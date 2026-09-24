version 1.0
# LICENSE
#   Copyright 2026 Ultima Genomics
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
#   Tasks for genome ploidy estimation and sex-aware VCF post-processing.
#   EstimatePloidyFromVcf: VCF mode — coverage from SNP DP + BAF analysis
#   EstimatePloidyFromCram: CRAM mode — mosdepth coverage only

import "structs.wdl" as Structs

task RunMosdepthSummary {
    input {
        File input_cram
        File input_cram_index
        References references
        String base_file_name
        String docker
        Int preemptibles
        File monitoring_script
        Int memory_gb = 4
        Int cpus = 4
    }

    Int ref_size = ceil(size(references.ref_fasta, "GB"))
    Int cram_size = ceil(size(input_cram, "GB"))
    Int disk_size = ref_size + cram_size + 20

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo "******** Running mosdepth (summary-only) ********"
        mosdepth \
            --threads ~{cpus} \
            --fasta ~{references.ref_fasta} \
            --no-per-base \
            --fast-mode \
            "~{base_file_name}" \
            ~{input_cram}

        echo "******** DONE ********"
    >>>

    output {
        File mosdepth_summary = "~{base_file_name}.mosdepth.summary.txt"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        preemptible: "~{preemptibles}"
        cpu: cpus
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }
}

task EstimatePloidyFromVcf {
    input {
        File input_vcf
        File input_vcf_index
        String sample_id
        Array[String] sex_chromosomes
        File? ploidy_exclude_regions_bed
        String docker
        Int preemptibles
        File monitoring_script
    }

    parameter_meta {
        ploidy_exclude_regions_bed: {
            help: "Optional BED of low-confidence regions (e.g. PAR-adjacent chrY) to exclude from the SNP coverage evidence used for ploidy estimation.",
            type: "File",
            category: "param_optional"
        }
    }

    Int vcf_size = ceil(size(input_vcf, "GB"))
    Int disk_size = vcf_size + 10

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo "******** Estimating ploidy from VCF ********"
        estimate_ploidy \
            --vcf ~{input_vcf} \
            --sample-id ~{sample_id} \
            --sex-chromosomes ~{sep=' ' sex_chromosomes} \
            --het-sample-count 5000 \
            ~{"--exclude-regions-bed " + ploidy_exclude_regions_bed} \
            --output-dir .

        # Extract karyotype string for downstream conditional logic
        if ! grep -m1 "Karyotype:" ~{sample_id}.ploidy_report.txt | awk '{print $NF}' > karyotype.txt; then
            echo "ERROR: Karyotype not found in ploidy report" >&2
            echo "UNDETERMINED" > karyotype.txt
        fi
        echo "******** DONE (karyotype=$(cat karyotype.txt)) ********"
    >>>

    output {
        File ploidy_report = "~{sample_id}.ploidy_report.txt"
        String karyotype = read_string("karyotype.txt")
        File monitoring_log = "monitoring.log"
    }

    runtime {
        preemptible: "~{preemptibles}"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }
}

task EstimatePloidyFromCram {
    input {
        File mosdepth_summary
        String sample_id
        Array[String] sex_chromosomes
        String docker
        Int preemptibles
        File monitoring_script
    }

    Int disk_size = 10

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo "******** Estimating ploidy from mosdepth summary ********"
        estimate_ploidy \
            --mosdepth-summary ~{mosdepth_summary} \
            --sample-id ~{sample_id} \
            --sex-chromosomes ~{sep=' ' sex_chromosomes} \
            --output-dir .

        if ! grep -m1 "Karyotype:" ~{sample_id}.ploidy_report.txt | awk '{print $NF}' > karyotype.txt; then
            echo "ERROR: Karyotype not found in ploidy report" >&2
            echo "UNDETERMINED" > karyotype.txt
        fi
        echo "******** DONE (karyotype=$(cat karyotype.txt)) ********"
    >>>

    output {
        File ploidy_report = "~{sample_id}.ploidy_report.txt"
        String karyotype = read_string("karyotype.txt")
        File monitoring_log = "monitoring.log"
    }

    runtime {
        preemptible: "~{preemptibles}"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }
}

task SubtractPloidyExcludeRegions {
    # Removes low-confidence ploidy regions (e.g. PAR-adjacent chrY) from the calling-regions
    # BED before it is used for CRAM-mode mosdepth coverage extraction.
    input {
        File input_bed
        File exclude_bed
        String docker
        Int preemptibles
        File monitoring_script
    }

    Int disk_size_subtract = ceil(size(input_bed, "GB") + size(exclude_bed, "GB")) + 10

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        echo "******** Subtracting ploidy-exclude regions from calling intervals ********"
        # input_bed may be a Picard interval_list (1-based, '@' header); convert to 0-based BED first
        if grep -q '^@' ~{input_bed}; then
            grep -v '^@' ~{input_bed} | awk 'BEGIN {OFS="\t"} {print $1, $2 - 1, $3}' > calling_intervals.bed
        else
            ln -s ~{input_bed} calling_intervals.bed
        fi
        bedtools subtract -a calling_intervals.bed -b ~{exclude_bed} > ploidy_calling_intervals.bed
        echo "******** DONE ********"
    >>>

    output {
        File output_bed = "ploidy_calling_intervals.bed"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        preemptible: "~{preemptibles}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk " + disk_size_subtract + " HDD"
        docker: docker
    }
}


