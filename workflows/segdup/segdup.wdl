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

# DESCRIPTION
# Segmental duplication processing workflow
#
# CHANGELOG

import "tasks/structs.wdl" as structs
import "tasks/general_tasks.wdl" as UGGeneral
import "tasks/globals.wdl" as GlobalsWDL
import "efficient_dv.wdl" as EDV

workflow SegDupAnalysis {
	input {
        String pipeline_version = "1.25.0" # !UnusedDeclaration
        String base_file_name
        File input_cram_bam
        File input_crai_bai
        References references
        File homology_table
        File homology_table_index
        File segdup_regions
        File background_bed
        File cn_model
        Int n_threads
        File model_onnx
        File? model_serialized
        File exome_intervals
        File dbsnp
        File dbsnp_index
        String? cloud_provider_override
        Int preemptible_tries = 3
        Boolean no_address = true
        # winval validations
        #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
        #@wv suffix(input_cram_bam) in {".bam", ".cram"}
        #@wv suffix(input_crai_bai) in {".bai", ".crai"}
        #@wv suffix(references['ref_fasta']) in {'.fasta', '.fa'}
        #@wv suffix(references['ref_dict']) == '.dict'
        #@wv suffix(references['ref_fasta_index']) == '.fai'
        #@wv prefix(references['ref_fasta_index']) == references['ref_fasta']
        #@wv suffix(homology_table) == '.gz'
        #@wv suffix(homology_table_index) == '.tbi'
	}
    meta {
        description: "Processes segmental duplications in the genome by collapsing all copies on a single copy of the segmental duplication" 
        author: "Ultima Genomics"
        WDL_AID: {
            exclude: ["pipeline_version",
                    "Globals.glob",
                    "DV.FilterVCF.ref_fasta",
                    "DV.FilterVCF.ref_fasta_idx",
                    "DV.UGMakeExamples.count_candidates_with_dvtools",
                    "DV.single_strand_filter", 
                    "DV.keep_duplicates", 
                    "DV.add_ins_size_channel",
                    "BedToIntervalList.disk_size"
            ]
        }
    }
    parameter_meta{
        base_file_name: {
            help: "Prefix of the output files",
            type: "String",
            category: "input_required"
        }
        input_cram_bam: {
            help: "Input CRAM/BAM file",
            type: "File",
            category: "input_required"
        }
        input_crai_bai: {
            help: "Input CRAM/BAM index file",
            type: "File",
            category: "input_required"
        }
        references: {
            help: "Reference genome files",
            type: "References",
            category: "input_required"
        }
        homology_table: {
            help: "Segmental duplication table (see parascopy), see template",
            type: "File",
            category: "input_advanced"
        }
        homology_table_index: {
            help: "Segmental duplication table index (see parascopy), see template",
            type: "File",
            category: "input_advanced"
        }
        segdup_regions: {
            help: "Segmental duplication regions. All reads will be remapped to `segdup_regions` and the calling will happen only on these regions, see template (BED file)",
            type: "File",
            category: "input_advanced"
        }
        background_bed: {
            help: "Background regions (non-segmental duplicated) for CNV calling, see template and `parascopy`",
            type: "File",
            category: "input_advanced"
        }
        cn_model: {
            help: "CNV model file from parascopy, see template",
            type: "File",
            category: "input_advanced"
        }
        n_threads: {
            help: "Number of threads to use",
            type: "Int",
            category: "input_required"
        }
        model_onnx: {
            help: "DeepVariant model for variant calling on segmental duplications, see template",
            type: "File",
            category: "input_advanced"
        }
        model_serialized: {
            help: "Serialized model for variant calling",
            type: "File",
            category: "input_advanced"
        }
        exome_intervals: {
            help: "Exome intervals for variant calling (required for deepVariant, otherwise not important)",
            type: "File",
            category: "input_required"
        }
        dbsnp: {
            help: "dbSNP reference file (for annotation)",
            type: "File",
            category: "input_required"
        }
        dbsnp_index: {
            help: "dbSNP reference index file (for annotation)",
            type: "File",
            category: "input_required"
        }
        cloud_provider_override: {
            help: "Cloud provider override (for running on other clouds): gcp or aws",
            type: "String",
            category: "input_optional"
        }
        preemptible_tries: {
            help: "Number of preemptible tries",
            type: "Int",
            category: "input_optional"
        }
        no_address: {
            help: "Start instances with no public IP address",
            type: "Boolean",
            category: "input_advanced"
        }
        remap_bam: {
            help: "Remapped BAM file",
            type: "File",
            category: "output"
        }
        remap_bam_index: {
            help: "Remapped BAM index file",
            type: "File",
            category: "output"
        }
        acnv_calls :{
            help: "CNV calls",
            type: "File",
            category: "output"
        }
        acnv_calls_index: {
            help: "CNV calls index",
            type: "File",
            category: "output"
        }
        small_variants: {
            help: "Small variants (VCF)",
            type: "File",
            category: "output"
        }
        small_variants_idx: {
            help: "Small variants index",
            type: "File",
            category: "output"
        }
    }
    call GlobalsWDL.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    call PoolReads {
        input:
            base_file_name = base_file_name,
            input_cram_bam = input_cram_bam,
            input_crai_bai = input_crai_bai,
            references = references,
            segdup_docker = global.segdup_docker, 
            homology_table = homology_table,
            homology_table_index = homology_table_index, 
            segdup_regions = segdup_regions,
            preemptible_tries = preemptible_tries,
            monitoring_script = global.monitoring_script    #!FileCoercion
    }

    call CallCNV {
        input: 
            input_cram_bam = input_cram_bam, 
            input_crai_bai = input_crai_bai, 
            homology_table = homology_table, 
            homology_table_index = homology_table_index,
            background_regions = background_bed, 
            cn_model = cn_model, 
            references = references, 
            segdup_docker = global.segdup_docker,
            base_file_name = base_file_name,
            n_threads = n_threads, 
            preemptible_tries = preemptible_tries,
            monitoring_script = global.monitoring_script   #!FileCoercion
    }

    call UGGeneral.BedToIntervalList {
        input:
            monitoring_script = global.monitoring_script,  #!FileCoercion
            input_file = segdup_regions,
            reference_dict = references.ref_dict,
            base_file_name = base_file_name, 
            docker = global.broad_gatk_docker,
            preemptible_tries = preemptible_tries,
            no_address = no_address
    }

    call EDV.EfficientDV as DV {
        input:
            base_file_name = base_file_name,
            cram_files = [PoolReads.remap_cram],
            cram_index_files = [PoolReads.remap_cram_index],
            references = references,
            make_gvcf = false,
            is_somatic = false,
            recalibrate_vaf = false,
            num_shards = 3,
            cap_at_optimal_coverage = false,
            optimal_coverages = [70],
            min_fraction_snps = 0.05,
            min_fraction_non_hmer_indels = 0.05,
            min_fraction_hmer_indels = 0.05,
            min_variant_quality_exome_hmer_indels = 5,
            min_fraction_single_strand_non_snps = 0.15,
            normalize_strand_bias = false,
            # Call variants args
            model_onnx = model_onnx,
            model_serialized = model_serialized,
            target_intervals = BedToIntervalList.interval_list, 
            exome_intervals = exome_intervals,
            ref_dbsnp = dbsnp,
            ref_dbsnp_index = dbsnp_index,

            # Used for running on other clouds (aws)
            cloud_provider_override = cloud_provider_override,
            preemptible_tries = preemptible_tries
    }

    call ParascopyCall {
        input: 
            base_file_name = base_file_name, 
            input_cram = input_cram_bam, 
            input_crai = input_crai_bai,
            homology_table = homology_table,
            homology_table_index = homology_table_index, 
            dv_vcf = DV.output_vcf,
            dv_vcf_index = DV.output_vcf_index, 
            cn_model = CallCNV.cnv_results, 
            references  = references,
            monitoring_script = global.monitoring_script,  #!FileCoercion
            segdup_docker = global.segdup_docker,
            n_threads = n_threads,
            preemptible_tries = preemptible_tries,
    }


    output { 
        File remap_bam = PoolReads.remap_cram
        File remap_bam_index = PoolReads.remap_cram_index
        File acnv_calls = CallCNV.acnv_calls
        File acnv_calls_index = CallCNV.acnv_calls_index
        File small_variants = ParascopyCall.small_variants
        File small_variants_idx = ParascopyCall.small_variants_index
    }

}

task PoolReads { 
    input {
        File input_cram_bam
        File input_crai_bai
        File homology_table
        File homology_table_index
        File segdup_regions
        References references
        String segdup_docker
        String base_file_name
        File monitoring_script
        Int preemptible_tries
    }
    Int disk_size = ceil(1.75 * size(input_cram_bam, "GB"))
    command <<<
        set -eo pipefail
        set -x
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        while IFS= read -r line  
        do 
            l1=$( echo "$line" | awk '{print $1":"$2+1"-"$3}' )
            parascopy pool -i ~{input_cram_bam} \
                       -t ~{homology_table} \
                    -f ~{references.ref_fasta} \
                    -o ~{base_file_name}."$l1".remap.bam \
                    -m 0 \
                    --tags_to_reverse t0 tp \
                    --tags_to_retain XA XB \
                    -r "$l1"
            echo "Output file" ~{base_file_name}."$l1".remap.bam
        done < ~{segdup_regions}

        find . -name "~{base_file_name}*:*.bam" > file.lst
        echo "Merging files:"
        cat file.lst
        samtools cat -o ~{base_file_name}.remap.tmp.bam -b file.lst 
        echo "Sorting files:"
        samtools sort --reference ~{references.ref_fasta} -o ~{base_file_name}.remap.cram ~{base_file_name}.remap.tmp.bam --output-fmt-option embed_ref=1
        samtools index ~{base_file_name}.remap.cram

    >>> 

    output {
        File remap_cram = "~{base_file_name}.remap.cram"
        File remap_cram_index = "~{base_file_name}.remap.cram.crai"
    }

    runtime { 
        docker: segdup_docker
        memory: "16 GB"
        preemptible: preemptible_tries
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"    
    }
}

task CallCNV { 
    input {
        File input_cram_bam
        File input_crai_bai
        File homology_table
        File homology_table_index
        File background_regions
        File cn_model
        References references
        String segdup_docker
        String base_file_name
        Int n_threads
        Int preemptible_tries
        File monitoring_script
    }
    Int disk_size = ceil(1.5 * size(input_cram_bam, "GB"))

    command <<< 
        set -eo pipefail
        set -x
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        
        tar --no-same-owner --no-same-permissions -xvf ~{cn_model}
        
        parascopy depth --input ~{input_cram_bam} \
        --bed-regions ~{background_regions} \
        --fasta-ref ~{references.ref_fasta} \
        -o ~{base_file_name}.depth -@~{n_threads}\
         --clipped-perc 20 --unpaired-perc 120

        parascopy "cn-using" model \
        --input ~{input_cram_bam} \
        --fasta-ref ~{references.ref_fasta} \
        -d ~{base_file_name}.depth -t ~{homology_table} \
        -o ~{base_file_name}.cn -@~{n_threads}

        #return the CN results
        tar cvzf ~{base_file_name}.cn.tar.gz ~{base_file_name}.cn
    >>>

    output {
        File acnv_calls = base_file_name + ".cn/res.samples.bed.gz"
        File acnv_calls_index = base_file_name + ".cn/res.samples.bed.gz.tbi"
        File pcnv_calls = base_file_name + ".cn/res.paralog.bed.gz"
        File pcnv_calls_index = base_file_name + ".cn/res.paralog.bed.gz.tbi"
        File cnv_results = base_file_name + ".cn.tar.gz"
    }

    runtime {
        docker: segdup_docker
        memory: "32 GB"
        preemptible: preemptible_tries
        cpu: n_threads
        disks: "local-disk " + disk_size + " HDD"
    }
}

task ParascopyCall {
    input {
        String base_file_name
        File input_cram
        File input_crai
        File homology_table
        File homology_table_index
        File dv_vcf
        File dv_vcf_index
        File cn_model
        References references
        File monitoring_script
        String segdup_docker
        Int n_threads
        Int preemptible_tries
    }

    Int disk_size = ceil(1.5 * size(input_cram, "GB"))

    command <<<
        set -eo pipefail
        set -x
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        
        tar --no-same-owner --no-same-permissions -xvzf ~{cn_model}

        parascopy call \
        -p ~{base_file_name}.cn \
        -i ~{input_cram} \
        --fasta-ref ~{references.ref_fasta} \
        -t ~{homology_table} \
        --freebayes /usr/local/bin/freebayes \
        --precalled-variants ~{dv_vcf} \
        -o ~{base_file_name}.calls -@~{n_threads} 
    >>>

    output {
        File small_variants = "~{base_file_name}.calls/variants.vcf.gz"
        File small_variants_index = "~{base_file_name}.calls/variants.vcf.gz.tbi"
    }
    runtime {
        docker: segdup_docker # Placeholder Docker image
        memory: "32 GB"
        cpu: n_threads
        preemptible: preemptible_tries
        disks: "local-disk " + disk_size + " HDD"
    }
}