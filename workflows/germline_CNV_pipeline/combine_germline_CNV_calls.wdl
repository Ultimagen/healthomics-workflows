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
# 1. Concatenate CNV calls from cn.mops and cnvpytor.
# 2. Annotate with gap information.
# 3. Annotate with split-read information.
# 4. Run jalign (an internal tool detecting alternative alignments supporting a CNV candidate) for CNV candidates.
# 5. Filter (apply ML model that returns confidence to variants - and places it on QUAL/TREE_SCORE)
# 6. Collapse the combined callset.(merge close PASS-filter CNVs)

import "tasks/globals.wdl" as Globals
import "tasks/single_sample_vc_tasks.wdl" as Filtering
import "tasks/cnv_calling_tasks.wdl" as CnvTasks
workflow CombineGermlineCNVCalls {

    input {
        String pipeline_version = "1.29.2" # !UnusedDeclaration

        String base_file_name

        File cnmops_cnvs_vcf
        File cnmops_cnvs_vcf_index
        File cnvpytor_cnvs_vcf
        File cnvpytor_cnvs_vcf_index
        Int cushion_size

        # jalign parameters
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        
        File? filtering_model 
        Int? filtering_model_decision_threshold
        Boolean skip_filtering
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
        #@wv suffix(filtering_model) == ".pkl"
    }

    meta {
        description: "Combine cn.mops and cnvpytor CNV calls for a single sample."
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                "monitoring_script_input",
                "no_address_override",
                "Glob.glob", 
                'FilterVCF.ref_fasta',
                'FilterVCF.ref_fasta_idx',
                'FilterVCF.blacklist_file',
                'FilterVCF.custom_annotations',
                'FilterVCF.disk_size'
                ]}
    }
    parameter_meta {
        base_file_name: {
            help: "Sample name",
            type: "String",
            category: "input_required"
        }
        cnmops_cnvs_vcf: {
            help: "cn.mops CNV calls in VCF format",
            type: "File",
            category: "input_required"
        }
        cnmops_cnvs_vcf_index: {
            help: "Index file for cn.mops CNV calls VCF",
            type: "File",
            category: "input_required"
        }
        cnvpytor_cnvs_vcf: {
            help: "cnvpytor CNV calls in VCF format",
            type: "File",
            category: "input_required"
        }
        cnvpytor_cnvs_vcf_index: {
            help: "Index file for cnvpytor CNV calls VCF",
            type: "File",
            category: "input_required"
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
       cushion_size: {
            help: "Cushion size around CNV breakpoints for split-read analysis and jump alignment analysis",
            type: "Int",
            category: "param_required"
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
        filtering_model: {
            help: "CNV filtering model file, set in template, can be removed, the callset is not filtered if not provided",
            type: "File",
            category: "input_optional"
        }
        filtering_model_decision_threshold: {
            help: "Decision threshold for the filtering model, default is set in template. Lower- less stringent, Higher- more stringent",
            type: "Int",
            category: "param_advanced"
        }
        skip_filtering: {
            help: "Whether to skip the filtering step even if a filtering model is provided. Default is false, useful for training models",
            type: "Boolean",
            category: "param_advanced"
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
        out_sample_cnvs_bed: {
            help: "Final (combined) CNV calls in bed format",
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
        realign_read_evidence:{
            help: "BAM file with read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        realign_read_evidence_index:{
            help: "Index file for the BAM with read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        split_read_evidence:{
            help: "BAM file with split-read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        split_read_evidence_index:{
            help: "Index file for the BAM with split-read evidence supporting combined CNV calls",
            type: "File",
            category: "output"
        }
        read_scores_csv:{
            help: "CSV file with jalign scores for each read",
            type: "File",
            category: "output"  
        }
    }

    Int preemptible_tries = select_first([preemptible_tries_override, 1])
    Boolean no_address = select_first([no_address_override, true ])

    call Globals.Globals as Glob
    GlobalVariables global = Glob.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script]) #!FileCoercion
    call ValidateSameSampleName {
        input:
            vcf1 = cnmops_cnvs_vcf,
            vcf2 = cnvpytor_cnvs_vcf,
            docker = global.bcftools_docker
    }

    call ConcatAndAnnotateCallsets {
        input:
            base_file_name = base_file_name,
            cnmops_vcf = cnmops_cnvs_vcf,
            cnmops_vcf_index = cnmops_cnvs_vcf_index,
            cnvpytor_vcf = cnvpytor_cnvs_vcf,
            cnvpytor_vcf_index = cnvpytor_cnvs_vcf_index,
            reference_genome = reference_genome,
            reference_genome_index = reference_genome_index,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }
    call AnnotateWithSplitReadsInfo {
        input:
            input_cnv_vcf = ConcatAndAnnotateCallsets.combined_cnv_calls_vcf,
            base_file_name = base_file_name,
            input_cram_bam_file = input_bam_file,
            input_cram_bam_file_index = input_bam_file_index,
            reference_genome = reference_genome,
            reference_genome_index = reference_genome_index,
            cushion_size = cushion_size,
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    call RunJalignForCNVCandidates {
        input:
            base_file_name = base_file_name,
            input_vcf = AnnotateWithSplitReadsInfo.annotated_cnv_calls_vcf,
            input_vcf_index = AnnotateWithSplitReadsInfo.annotated_cnv_calls_vcf_index,
            
            input_bam_file = input_bam_file,
            input_bam_file_index = input_bam_file_index,
            reference_genome = reference_genome,
            reference_genome_index = reference_genome_index,            
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptible_tries,
            no_address = no_address
    }

    call RefineCNVBreakpoints {
        input: 
            base_file_name = base_file_name,
            input_vcf = RunJalignForCNVCandidates.output_vcf,
            input_vcf_index = RunJalignForCNVCandidates.output_vcf_index,
            evidence_bam = [ RunJalignForCNVCandidates.output_bam, AnnotateWithSplitReadsInfo.evidence_bam ],
            evidence_bam_index = [ RunJalignForCNVCandidates.output_bam_index, AnnotateWithSplitReadsInfo.evidence_bam_index ],
            docker = global.ugbio_cnv_docker,
            monitoring_script = monitoring_script,
            preemptible_tries = preemptible_tries,
            no_address = no_address
    }

    Array[String] custom_annotations =     [
        "CNV_SOURCE",
        "RoundedCopyNumber",
        "CopyNumber",
        "pytorQ0",
        "pytorP2",
        "pytorRD",
        "pytorP1",
        "pytorP3",
        "CN",
        "GAP_PERCENTAGE",
        "CNV_DUP_READS",
        "CNV_DEL_READS",
        "CNV_DUP_FRAC",
        "CNV_DEL_FRAC",
        "JALIGN_DUP_SUPPORT",
        "JALIGN_DEL_SUPPORT",
        "JALIGN_DUP_SUPPORT_STRONG",
        "JALIGN_DEL_SUPPORT_STRONG",
        "DUP_READS_MEDIAN_INSERT_SIZE",
        "DEL_READS_MEDIAN_INSERT_SIZE",
        "SVTYPE",
        "SVLEN",
        "CIPOS"
    ]
    if (defined(filtering_model) && !skip_filtering) {
        call Filtering.FilterVCF{
            input:
                input_vcf = RefineCNVBreakpoints.output_vcf,
                input_vcf_index = RefineCNVBreakpoints.output_vcf_index,
                input_model = filtering_model,
                filter_cg_insertions = false,
                final_vcf_base_name = base_file_name,
                recalibrate_gt = false,
                custom_annotations = custom_annotations,
                decision_threshold = filtering_model_decision_threshold, 
                overwrite_quality = true,
                monitoring_script = monitoring_script,
                preemptible_tries = preemptible_tries,
                docker = global.ugbio_filtering_docker, 
                no_address = no_address
        }   

        call CollapseCallset{
            input:
                input_vcf = FilterVCF.output_vcf_filtered,
                input_vcf_index = FilterVCF.output_vcf_filtered_index,
                base_file_name = base_file_name,
                docker = global.ugbio_cnv_docker,
                monitoring_script = monitoring_script,
                no_address = no_address,
                preemptible_tries = preemptible_tries
        } 
    }
    call CnvTasks.CnvVcfToBed {
        input:
            input_cnv_vcf = select_first([CollapseCallset.output_vcf,RefineCNVBreakpoints.output_vcf]),
            base_file_name = base_file_name,
            docker = global.bcftools_docker,
            monitoring_script = monitoring_script,
            no_address = no_address,
            preemptible_tries = preemptible_tries
    }

    output {
        File out_sample_cnvs_bed = CnvVcfToBed.output_cnv_bed
        File out_sample_cnvs_vcf = select_first([CollapseCallset.output_vcf, RefineCNVBreakpoints.output_vcf])
        File out_sample_cnvs_vcf_index = select_first([CollapseCallset.output_vcf_index, RefineCNVBreakpoints.output_vcf_index])
        File realign_read_evidence = RunJalignForCNVCandidates.output_bam
        File realign_read_evidence_index = RunJalignForCNVCandidates.output_bam_index
        File split_read_evidence = AnnotateWithSplitReadsInfo.evidence_bam
        File split_read_evidence_index = AnnotateWithSplitReadsInfo.evidence_bam_index
        File read_scores_csv = RunJalignForCNVCandidates.scores_csv
    }
}

task ValidateSameSampleName {
    input {
        File vcf1
        File vcf2
        String docker
    }
    command <<<
        set -xeo pipefail
        sample_name_1=$(bcftools query -l ~{vcf1})
        sample_name_2=$(bcftools query -l ~{vcf2})

        if [ "$sample_name_1" != "$sample_name_2" ]; then
            echo "Error: Sample names do not match: $sample_name_1 != $sample_name_2" >&2
            exit 1
        fi
    >>>
    runtime {
        docker: docker
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 2 HDD"
    }
}

task ConcatAndAnnotateCallsets {
    input {
        String base_file_name
        File cnmops_vcf
        File cnmops_vcf_index
        File cnvpytor_vcf
        File cnvpytor_vcf_index
        File reference_genome
        File reference_genome_index
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        combine_cnmops_cnvpytor_cnv_calls concat \
            --make_ids_unique \
            --cnmops_vcf ~{cnmops_vcf} \
            --cnvpytor_vcf ~{cnvpytor_vcf} \
            --output_vcf ~{base_file_name}.step1.vcf.gz \
            --fasta_index ~{reference_genome_index} 
            
        combine_cnmops_cnvpytor_cnv_calls annotate_gaps \
            --calls_vcf ~{base_file_name}.step1.vcf.gz \
            --output_vcf ~{base_file_name}.combined.vcf.gz \
            --ref_fasta ~{reference_genome}
         
    >>>
    runtime {
        docker: docker
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk 10 HDD"
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File combined_cnv_calls_vcf = "~{base_file_name}.combined.vcf.gz"
        File combined_cnv_calls_vcf_index = "~{base_file_name}.combined.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task AnnotateWithSplitReadsInfo {
    input {
        File input_cnv_vcf
        String base_file_name
        File input_cram_bam_file
        File input_cram_bam_file_index
        File reference_genome
        File reference_genome_index
        Int cushion_size
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int disk_size = ceil(size(input_cram_bam_file, "GB") + size(input_cnv_vcf, "GB") + 10)
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        combine_cnmops_cnvpytor_cnv_calls analyze_breakpoint_reads \
            --bam-file ~{input_cram_bam_file} \
            --vcf-file ~{input_cnv_vcf} \
            --reference-fasta ~{reference_genome} \
            --cushion ~{cushion_size} \
            --output-file ~{base_file_name}.split.annotated.vcf.gz \
            --output-bam ~{base_file_name}.split.annotated.bam
        
        bcftools index -t ~{base_file_name}.split.annotated.vcf.gz
        samtools sort -o ~{base_file_name}.split.annotated.sort.bam ~{base_file_name}.split.annotated.bam
        samtools index ~{base_file_name}.split.annotated.sort.bam
    >>>
    runtime {
        docker: docker
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        noAddress: no_address
    }
    output {
        File annotated_cnv_calls_vcf = "~{base_file_name}.split.annotated.vcf.gz"
        File annotated_cnv_calls_vcf_index = "~{base_file_name}.split.annotated.vcf.gz.tbi"
        File evidence_bam = "~{base_file_name}.split.annotated.sort.bam"
        File evidence_bam_index = "~{base_file_name}.split.annotated.sort.bam.bai"
        File monitoring_log = "monitoring.log"
    }
}

task RunJalignForCNVCandidates {
    input {
        String base_file_name
        File input_vcf
        File input_vcf_index
        
        File input_bam_file
        File input_bam_file_index
        File reference_genome
        File reference_genome_index
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int cpu = 32
    Float input_bam_file_size = size(input_bam_file, "GB")
    Float reference_fasta_file_size = size(reference_genome, "GB")
    Float additional_disk = 10

    Int disk_size = ceil(input_bam_file_size + input_bam_file_size + reference_fasta_file_size + additional_disk)
     
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        run_jalign --tool-path /opt/para_jalign \
                ~{input_bam_file} \
                ~{input_vcf} \
                ~{reference_genome} \
                ~{base_file_name} \
                --max-reads-per-cnv 300 \
                --threads ~{cpu}        
        bcftools index -tf ~{base_file_name}.jalign.vcf.gz
        samtools sort -o ~{base_file_name}.jalign.sort.bam ~{base_file_name}.jalign.bam
        samtools index ~{base_file_name}.jalign.sort.bam
    >>>
    runtime {
        memory: "32 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        preemptible: preemptible_tries
        cpu: cpu
    }

    output {
        File output_vcf = "~{base_file_name}.jalign.vcf.gz"
        File output_vcf_index = "~{base_file_name}.jalign.vcf.gz.tbi"
        File output_bam = "~{base_file_name}.jalign.sort.bam"
        File output_bam_index = "~{base_file_name}.jalign.sort.bam.bai"
        File scores_csv = "~{base_file_name}.jalign.csv"
        File monitoring_log = "monitoring.log"
    }
}

task RefineCNVBreakpoints {
    input{
        String base_file_name
        File input_vcf
        File input_vcf_index
        Array[File] evidence_bam
        Array[File] evidence_bam_index
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int disk_size = ceil(2 * size(input_vcf, "GB") + size(evidence_bam, "GB") + 10)
    command <<< 
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        refine_cnv_breakpoints \
            --input-vcf ~{input_vcf} \
            --output-vcf ~{base_file_name}.refined.tmp.vcf.gz \
            --bam-files ~{sep=" " evidence_bam}
        bcftools sort -Oz -o ~{base_file_name}.refined.vcf.gz ~{base_file_name}.refined.tmp.vcf.gz
        bcftools index -tf ~{base_file_name}.refined.vcf.gz
    >>>

    runtime {
        memory: "4 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        preemptible: preemptible_tries
        cpu: 1
    }

    output {
        File output_vcf = "~{base_file_name}.refined.vcf.gz"
        File output_vcf_index = "~{base_file_name}.refined.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }


}
task CollapseCallset {
    input{ 
        File input_vcf
        File input_vcf_index
        String base_file_name
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        combine_cnmops_cnvpytor_cnv_calls merge_records \
            --input_vcf ~{input_vcf} \
            --output_vcf ./~{base_file_name}.mrg.vcf.gz \
            --distance 0 \
            --enable_smoothing \
            --max_gap_absolute 50000 \
            --gap_scale_fraction 0.05 \
            --cipos_threshold 50
    >>>
    runtime {
        memory: "4 GiB"
        disks: "local-disk " + "4 HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
        preemptible: preemptible_tries
    }

    output {
        File output_vcf = "~{base_file_name}.mrg.vcf.gz"
        File output_vcf_index = "~{base_file_name}.mrg.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }

}
