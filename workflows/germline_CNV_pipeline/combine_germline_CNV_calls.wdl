version 1.0
# LICENSE
#   Copyright 2025 Ultima Genomics
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
# Combine cn.mops and cnvpytor CNV calls for a single sample.:
# 1. seperate deletions and duplications cnv calls.
# 3. Runs jalign (an internal tool detecting alternative alignments supporting a CNV candidate) for deletion candidates. 
# 4. combines calls from all tools and filters them based on the intersection with UG-CNV-LCR regions.

# CHANGELOG in reverse chronological order

import "tasks/globals.wdl" as Globals

workflow CombineGermlineCNVCalls {

    input {
        String pipeline_version = "1.23.0" # !UnusedDeclaration

        String base_file_name

        File cnmops_cnvs_bed
        File cnvpytor_cnvs_tsv
        Float? cnvpytor_precent_gaps_threshold_override
        Int? distance_threshold_override
        Int? deletions_length_cutoff_override
        Int? jalign_written_cutoff_override
        Int? jalign_min_mismatches_override
        Int? duplication_length_cutoff_for_cnmops_filter_override

        # jalign parameters
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        
        File? cnv_lcr_file
        
        Boolean? no_address_override
        Int? preemptible_tries_override

        # Used for running on other clouds (aws)
        File? monitoring_script_input

        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv prefix(input_bam_file_index) == input_bam_file
        #@wv suffix(input_bam_file) in {".bam", ".cram"}
        #@wv suffix(input_bam_file_index) in {".bai", ".crai"}
        #@wv reference_genome == prefix(reference_genome_index)
        #@wv suffix(reference_genome) in {'.fasta', '.fa', '.fna'}
        #@wv suffix(reference_genome_index) == '.fai'
    }

    meta {
        description: "Combine cn.mops and cnvpytor CNV calls for a single sample."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "no_address_override",
                "Glob.glob"
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        cnmops_cnvs_bed: {
            help: "cn.mops CNV calls in bed format",
            type: "File",
            category: "input_required"
        }
        cnvpytor_cnvs_tsv: {
            help: "cnvpytor CNV calls in bed format",
            type: "File",
            category: "input_required"
        }
        cnvpytor_precent_gaps_threshold_override: {
            help: "Threshold for pN (fraction of reference genome gaps (N's) in call region) for cnvpytor calls. default=0.9",
            type: "Float",
            category: "param_advanced"
        }
        distance_threshold_override: {
            help: "Distance threshold for merging CNV calls. default=1500",
            type: "Int",
            category: "param_advanced"
        }
        deletions_length_cutoff_override: {
            help: "Minimum length of deletions to be considered without jalign support. default=3000",
            type: "Int",
            category: "param_advanced"
        }
        jalign_written_cutoff_override: {
            help: "Minimal number of supporting jaligned reads for deletions. default=1",
            type: "Int",
            category: "param_advanced"
        }
        jalign_min_mismatches_override: {
            help: "Minimum number of mismatches for jalign to consider a read as supporting a CNV candidate. default=1",
            type: "Int",
            category: "param_advanced"
        }
        duplication_length_cutoff_for_cnmops_filter_override: {
            help: "Minimum length of duplications to be considered for cn.mops calls. default=10000",
            type: "Int",
            category: "param_advanced"
        }
        input_bam_file: {
            help: "Input sample BAM/CRAM file.",
            type: "File",
            category: "input_required"
         }
        input_bam_file_index: {
            help:"Input sample BAI/CRAI index file",
            type: "File",
            category: "input_required"
       }
        reference_genome: {
            help: "Genome fasta file associated with the CRAM file",
            type: "File",
            category: "ref_required"
         }
        reference_genome_index: {
            help : "Fai index of the fasta file",
            type: "File",
            category: "ref_required"
        }
        cnv_lcr_file: {
            help: "UG-CNV-LCR bed file",
            type: "File",
            category: "input_optional"
        }
        no_address_override: {
            help: "Whether to use the --no-address flag in docker run commands. Default is: true",
            type: "Boolean",
            category: "param_advanced"
        }
        preemptible_tries_override: {
            help: "Number of preemptible tries,default is: 1",
            type: "Int",
            category: "param_optional"
        }
        monitoring_script_input: {
            help: "Monitoring script for the docker run commands",
            type: "File",
            category: "input_advanced"
        }        
        out_jalign_del_bed:{
            help: "Output file with DEL candidates after jalign",
            type: "File",
            category: "output"
        }
        out_sample_cnvs_vcf:{
            help: "VCF file with sample's called CNVs",
            type: "File",
            category: "output"
        }
        out_sample_cnvs_vcf_index:{
            help: "Index file for the VCF file with sample's called CNVs",
            type: "File",
            category: "output"
        }
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])
    Int distance_threshold = select_first([distance_threshold_override,1500])
    Int deletions_length_cutoff = select_first([deletions_length_cutoff_override,3000])
    Int jalign_written_cutoff = select_first([jalign_written_cutoff_override,1])
    Int duplication_length_cutoff_for_cnmops_filter = select_first([duplication_length_cutoff_for_cnmops_filter_override,10000])
    Float cnvpytor_precent_gaps_threshold = select_first([cnvpytor_precent_gaps_threshold_override, 0])
    Int jalign_min_mismatches = select_first([jalign_min_mismatches_override,1])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call RunJalignForDelCandidates {
        input:
        base_file_name = base_file_name,
        cnmops_cnvs_bed = cnmops_cnvs_bed,
        cnvpytor_cnvs_tsv = cnvpytor_cnvs_tsv,
        cnvpytor_precent_gaps_threshold = cnvpytor_precent_gaps_threshold,
        distance_threshold = distance_threshold,
        
        input_bam_file = input_bam_file,
        input_bam_file_index = input_bam_file_index,
        reference_genome = reference_genome,
        reference_genome_index = reference_genome_index,
        jalign_min_mismatches = jalign_min_mismatches,
        
        docker = global.ug_jalign_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }
    call ProcessCnvCalls  {
        input:
        cnmops_cnvs_bed = cnmops_cnvs_bed,
        cnvpytor_cnvs_tsv = cnvpytor_cnvs_tsv,
        jalign_del_candidates = RunJalignForDelCandidates.out_jalign_del_bed,
        base_file_name = base_file_name,
        distance_threshold = distance_threshold,
        deletions_length_cutoff = deletions_length_cutoff,
        jalign_written_cutoff = jalign_written_cutoff,
        duplication_length_cutoff_for_cnmops_filter = duplication_length_cutoff_for_cnmops_filter,
        cnv_lcr_file = cnv_lcr_file,
        reference_fasta = reference_genome,
        fasta_index = reference_genome_index,
        docker = global.ugbio_cnv_docker,
        monitoring_script = monitoring_script,
        no_address = no_address,
        preemptible_tries = preemptible_tries
    }   
    
    output {
        File out_jalign_del_bed = RunJalignForDelCandidates.out_jalign_del_bed
        File out_sample_cnvs_vcf = ProcessCnvCalls.sample_cnvs_vcf_file
        File out_sample_cnvs_vcf_index = ProcessCnvCalls.sample_cnvs_vcf_index_file
    }
}

task RunJalignForDelCandidates {
    input {
        String base_file_name
        File cnmops_cnvs_bed
        File cnvpytor_cnvs_tsv
        Float cnvpytor_precent_gaps_threshold
        Int distance_threshold
        
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index

        Int jalign_min_mismatches

        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int cpu = 48
    Float input_bam_file_size = size(input_bam_file, "GB")
    Float reference_fasta_file_size = size(reference_genome, "GB")
    Float additional_disk = 100

    Int disk_size = ceil(input_bam_file_size + input_bam_file_size + reference_fasta_file_size + additional_disk)
     
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
                
        echo "Generating DEL candidates"
        #cnmops output format example:
        #chr1    632000  634000  CN0
        #chr1    124740000       124742000       CN0,CN1
        cat ~{cnmops_cnvs_bed} | sed 's/UG-CNV-LCR//g' | sed 's/LEN//g' | sed 's/|//g' | grep -E "CN0|CN1" | \
            bedtools merge -c 4 -o distinct -d ~{distance_threshold} -i - \
            > ~{base_file_name}.cnmops.DEL.merged.bed
        
        #cnvpytor output format example: 
        #chr1    123468001       124437000       deletion,969000
        #chr1    124440001       124511500       deletion,16500,deletion,54000
        cat ~{cnvpytor_cnvs_tsv} | grep "deletion" | awk '$(NF-1)<=~{cnvpytor_precent_gaps_threshold}' | \
            cut -f1-3 | sed 's/:/\t/' | sed 's/-/\t/' | \
             awk '{print $2"\t"$3"\t"$4"\t"$1","$5}' | \
            bedtools merge -c 4 -o distinct -d ~{distance_threshold} -i - \
            > ~{base_file_name}.cnvpytor.DEL.merged.bed

        cat ~{base_file_name}.cnmops.DEL.merged.bed ~{base_file_name}.cnvpytor.DEL.merged.bed | \
            bedtools sort -i - \
            > ~{base_file_name}.cnmops_cnvpytor.DEL.bed

        echo "split DEL candidates for parallel processing"
        mkdir cnmops500mod_cnvpytor500_DEL_split
        cat ~{base_file_name}.cnmops_cnvpytor.DEL.bed | split -l 10 -d --additional-suffix .bed
        mv x*.bed cnmops500mod_cnvpytor500_DEL_split/

        echo "Running jalign for DEL candidates"
        python /jalign/scripts/parallel_run_cnv_realign.py \
            --folder_with_cnv_del_bed_files cnmops500mod_cnvpytor500_DEL_split \
            --input_cram ~{input_bam_file} \
            --ref_fasta ~{reference_genome} \
            --out_folder out_jalign \
            --sample_name ~{base_file_name} \
            --min_mismatches ~{jalign_min_mismatches} \
            --mode DEL \
            --num_jobs ~{cpu} 
            
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "96 GiB"
        disks: "local-disk " + ceil(disk_size) + " LOCAL"
        docker: docker
        noAddress: no_address
        cpu: cpu
    }

    output {
        File out_jalign_del_bed = "out_jalign/DEL_jalign_merged_results/~{base_file_name}.DEL.jalign.bed"
        File out_jalign_del_bam = "out_jalign/DEL_jalign_merged_results/~{base_file_name}.DEL.jalign.bam"
        File out_jalign_del_bam_index = "out_jalign/DEL_jalign_merged_results/~{base_file_name}.DEL.jalign.bam.bai"
    }
}

task ProcessCnvCalls  {
    input {
        File cnmops_cnvs_bed
        File cnvpytor_cnvs_tsv
        File jalign_del_candidates
        String base_file_name
        File? cnv_lcr_file
        File reference_fasta
        File fasta_index
        Int? distance_threshold
        Int? deletions_length_cutoff
        Int? jalign_written_cutoff
        Int? duplication_length_cutoff_for_cnmops_filter
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Float cnmops_cnvs_bed_size = size(cnmops_cnvs_bed, "GB")
    Float cnvpytor_cnvs_tsv_size = size(cnvpytor_cnvs_tsv, "GB")
    Float jalign_del_candidates_size = size(jalign_del_candidates, "GB")
    Float cnv_lcr_file_size = size(cnv_lcr_file, "GB")   
    Float additional_disk = 10
    Int disk_size = ceil(cnmops_cnvs_bed_size + cnvpytor_cnvs_tsv_size + jalign_del_candidates_size + cnv_lcr_file_size + additional_disk)

     
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        combine_cnmops_cnvpytor_cnv_calls \
            --cnmops_cnv_calls ~{cnmops_cnvs_bed} \
            --cnvpytor_cnv_calls ~{cnvpytor_cnvs_tsv} \
            --del_jalign_merged_results ~{jalign_del_candidates} \
            ~{"--deletions_length_cutoff " + deletions_length_cutoff} \
            ~{"--jalign_written_cutoff " + jalign_written_cutoff} \
            ~{"--distance_threshold "+ distance_threshold} \
            ~{"--duplication_length_cutoff_for_cnmops_filter " + duplication_length_cutoff_for_cnmops_filter} \
            ~{"--ug_cnv_lcr " + cnv_lcr_file} \
            --ref_fasta ~{reference_fasta} \
            --fasta_index ~{fasta_index} \
            --out_directory . \
            --sample_name ~{base_file_name}

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "4 GiB"
        disks: "local-disk " + ceil(disk_size) + " LOCAL"
        docker: docker
        noAddress: no_address
        cpu: 2
    }

    output
    {
        File sample_cnvs_vcf_file = "~{base_file_name}.cnv.vcf.gz"
        File sample_cnvs_vcf_index_file = "~{base_file_name}.cnv.vcf.gz.tbi"
    }
}
