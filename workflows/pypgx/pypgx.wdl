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
# Runs PyPGx workflow for pharmacogenomics analysis.
#
# CHANGELOG
# 1.0 - Initial version

import "tasks/general_tasks.wdl" as UGGeneralTasks
import "tasks/structs.wdl" as Structs
import "efficient_dv.wdl" as UGEfficientDV
import "tasks/globals.wdl" as Globals
import "tasks/alignment_tasks.wdl" as UGAlignment

workflow PyPGx {
    input {
        String base_file_name
        File cram_file
        File cram_index_file
        Array[String] gene_symbols
        File? input_vcf_file
        File? input_vcf_index_file

        References references

        # EfficientDV parameters
        File? model_onnx
        File? exome_intervals
        File? ref_dbsnp
        File? ref_dbsnp_index

        Boolean? no_address_override
        Int preemptible_tries = 1

        # Used for running on other clouds (aws)
        String? cloud_provider_override
        File? monitoring_script_input
        Array[File] ref_files_for_tarball

        ##@wv suffix(cram_file) <= {".bam", ".cram"}
        ##@wv suffix(cram_index_file) <= {".bai", ".crai"}

    }

    meta {
        description: "Runs pharmacogenomics analysis on several genes."
        author: "Ultima Genomics"
        version: "1.0"
        WDL_AID: {exclude: ["pipeline_version",
                "cloud_provider_override",
                "monitoring_script_input",
                "no_address_override",
                "DepthOfCoverage.disk_size",
                "Globals.glob",
                "preemptible_tries",
                "CreateReferenceCache.disk_size"
                ]
        }
    }

    parameter_meta {
        base_file_name: {
            help: "Prefix for name of all output files",
            category: "input_required"
        }
        cram_file: {
            help: "Input (bam/cram) file",
            type: "File",
            category: "input_required"
        }
        cram_index_file: {
            help: "Input (bai/crai) index file",
            type: "File",
            category: "input_required"
        }
        gene_symbols: {
            help: "List of gene symbols to analyze",
            type: "Array[String]",
            category: "input_required"
        }
        references: {
            type: "References",
            help: "Reference files: fasta, dict and fai",
            category: "ref_required"
        }
        input_vcf_file: {
            help: "Input VCF file with variants. If not provided, Efficient DV will be run",
            type: "File",
            category: "input_optional"
        }
        input_vcf_index_file: {
            help: "Input VCF index file",
            type: "File",
            category: "input_optional"
        }
        preemptible_tries: {
            help: "Number of preemptible tries",
            category: "param_advanced"
        }
        model_onnx: {
            help: "TensorRT model for calling variants (onnx format)",
            category: "input_optional"
        }
        exome_intervals: {
            help: "A bed file with exome intervals. Used at the post-processing step to annotate the vcf and modify the FILTER of variants in the exome.",
            category: "input_optional"
        }
        ref_dbsnp: {
            help: "DbSNP vcf for the annotation of known variants",
            category: "input_optional"
        }
        ref_dbsnp_index: {
            help: "DbSNP vcf index",
            category: "input_optional"
        }
        ref_files_for_tarball: {
            help: "List of references for CreateReferenceCache task.",
            category: "input_required"
        }
        allele_fraction_profiles: {
            help: "Allele fraction profiles for each gene",
            category: "output"
        }
        alleles: {
            help: "Alleles for each gene",
            category: "output"
        }
        cnv_calls: {
            help: "CNV calls for each gene",
            category: "output"
        }
        consolidated_variants: {
            help: "Consolidated variants for each gene",
            category: "output"
        }
        copy_number_profiles: {
            help: "Copy number profiles for each gene",
            category: "output"
        }
        copy_numbers: {
            help: "Copy numbers for each gene",
            category: "output"
        }
        genotypes: {
            help: "Genotypes for each gene",
            category: "output"
        }
        imported_variants: {
            help: "Imported variants for each gene",
            category: "output"
        }
        phased_variants: {
            help: "Phased variants for each gene",
            category: "output"
        }
        phenotypes: {
            help: "Phenotypes for each gene",
            category: "output"
        }
        read_depths: {
            help: "Read depths for each gene",
            category: "output"
        }
        results: {
            help: "Results for each gene",
            category: "output"
        }
     }


    String cloud_provider = select_first([cloud_provider_override, 'gcp'])
    Boolean no_address = select_first([no_address_override, true ])

    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    File monitoring_script = select_first([monitoring_script_input, global.monitoring_script])

    call UGAlignment.CreateReferenceCache {
        input:
            references = ref_files_for_tarball,
            cache_populate_script = global.ref_cache_script, #!StringCoercion
            preemptible_tries = preemptible_tries,
            docker = global.perl_docker,
            dummy_input_for_call_caching = ""
    }


    if (!defined(input_vcf_file)) {
        call ExtractGeneIntervals {
            input:
                gene_symbols = gene_symbols,
                ref_dict = references.ref_dict,
                monitoring_script = monitoring_script,
                no_address = no_address,
                docker = global.pypgx_docker,
                preemptible_tries = preemptible_tries
        }
        # PyPGx.preemptible_tries
        call UGEfficientDV.EfficientDV as EfficientDV {
        input:
            base_file_name = base_file_name,
            cram_files = [cram_file],
            cram_index_files = [cram_index_file],
            background_cram_files = [],
            background_cram_index_files = [],
            references = select_first([references]),
            make_gvcf = false,
            recalibrate_vaf = false,
            is_somatic = false,
            candidate_min_mapping_quality = 0, # Set to 0 so variants will be called also in regions of high homology such as CYP2D6
            pileup_min_mapping_quality = 0, # Set to 0 so variants will be called also in regions of high homology such as CYP2D6
            keep_duplicates = true,
            target_intervals = ExtractGeneIntervals.interval_list,
            model_onnx = select_first([model_onnx]),
            exome_intervals = select_first([exome_intervals]),
            ref_dbsnp       = select_first([ref_dbsnp]),
            ref_dbsnp_index = select_first([ref_dbsnp_index]),
            min_variant_quality_hmer_indels = 5,
            min_variant_quality_non_hmer_indels = 0,
            min_variant_quality_snps = 0,
            min_variant_quality_exome_hmer_indels = 7,
            optimal_coverages = [ 50 ],
            cap_at_optimal_coverage = false,
            normalize_strand_bias = false,
            num_shards = 20,
            single_strand_filter = false,
            min_fraction_single_strand_non_snps = 0,
            add_ins_size_channel = true,
            cloud_provider_override = cloud_provider,
            monitoring_script_input = monitoring_script,
            preemptible_tries = preemptible_tries
        }
    }

    File vcf_file = select_first([input_vcf_file, EfficientDV.output_vcf])
    File vcf_index_file = select_first([input_vcf_index_file, EfficientDV.output_vcf_index])

    call DepthOfCoverage{
        input:
            base_file_name = base_file_name,
            cram_file = cram_file,
            cram_index_file = cram_index_file,
            cache_tarball = CreateReferenceCache.cache_tarball,
            monitoring_script = monitoring_script,
            no_address_override = no_address,
            preemptible_tries = preemptible_tries,
            docker = global.pypgx_docker
    }

    scatter (gene in gene_symbols){
        call NGSPipeline {
            input:
                gene_symbol = gene,
                vcf_file = vcf_file,
                vcf_index_file = vcf_index_file,
                depth_of_coverage = DepthOfCoverage.depth_of_coverage,
                control_statistics = DepthOfCoverage.control_statistics,
                monitoring_script = monitoring_script,
                no_address_override = no_address,
                preemptible_tries = preemptible_tries,
                docker = global.pypgx_docker
        }
    }

    call ConcatResults{
        input:
            base_file_name = base_file_name,
            result_files = NGSPipeline.results,
            gene_symbols = gene_symbols,
            monitoring_script = monitoring_script,
            no_address_override = no_address,
            preemptible_tries = preemptible_tries,
            docker = global.pypgx_docker
    }

    output {
        Array[File] allele_fraction_profiles = NGSPipeline.allele_fraction_profile
        Array[File] alleles = NGSPipeline.alleles
        Array[File] cnv_calls = NGSPipeline.cnv_calls
        Array[File] consolidated_variants = NGSPipeline.consolidated_variants
        Array[File] copy_number_profiles = NGSPipeline.copy_number_profile
        Array[File] copy_numbers = NGSPipeline.copy_number
        Array[File] genotypes = NGSPipeline.genotypes
        Array[File] imported_variants = NGSPipeline.imported_variants
        Array[File] phased_variants = NGSPipeline.phased_variants
        Array[File] phenotypes = NGSPipeline.phenotypes
        Array[File] read_depths = NGSPipeline.read_depth
        File results = ConcatResults.results
    }
}

# ----- TASKS -----

task ExtractGeneIntervals {
    input {
        Array[String] gene_symbols
        File ref_dict
        File monitoring_script
        Boolean no_address = true
        String docker
        Int preemptible_tries
    }

    command <<<
    set -xeo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &

    pypgx create-regions-bed --assembly GRCh38 --add-chr-prefix --genes ~{sep=" " gene_symbols} > unsorted.bed

    gatk BedToIntervalList -I unsorted.bed -O genes.interval_list --SEQUENCE_DICTIONARY ~{ref_dict}

    >>>

    output {
        File interval_list = "genes.interval_list"
    }

  runtime {
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
    docker: docker
    preemptible: preemptible_tries
    noAddress: no_address
  }
}


task DepthOfCoverage {
    input {
        String base_file_name
        File cram_file
        File cram_index_file
        File cache_tarball
        File monitoring_script
        Int disk_size = ceil(1.1*size(cram_file, 'GB') + size(cache_tarball, 'GB') + 11)
        Boolean no_address_override = true
        String docker
        Int preemptible_tries
    }
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # When the genome reference is not embedded in the cram, it has to be provided (patch)
        ~{"tar -zxf "+cache_tarball}
        export REF_CACHE=cache/%2s/%2s/ 
        export REF_PATH='.' 

        pypgx prepare-depth-of-coverage --assembly GRCh38 ~{base_file_name}.coverage.zip ~{cram_file}
        pypgx compute-control-statistics --assembly GRCh38 VDR ~{base_file_name}_VDR_coverage.zip ~{cram_file}

    >>>

    output {
        File depth_of_coverage = "~{base_file_name}.coverage.zip"
        File control_statistics = "~{base_file_name}_VDR_coverage.zip"
    }

  runtime {
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    preemptible: preemptible_tries
    noAddress: no_address_override
  }
}


task NGSPipeline {
    input {
        String gene_symbol
        File vcf_file
        File vcf_index_file
        File depth_of_coverage
        File control_statistics
        File monitoring_script
        Boolean no_address_override = true
        String docker
        Int preemptible_tries
    }

    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        pypgx run-ngs-pipeline ~{gene_symbol} outputs  --assembly GRCh38 --variants ~{vcf_file} --depth-of-coverage ~{depth_of_coverage} --control-statistics ~{control_statistics}

        # Some files do not get created, so I use touch to create empty files if they do not exist
        touch outputs/imported-variants.zip
        touch outputs/consolidated-variants.zip
        touch outputs/genotypes.zip
        touch outputs/phenotypes.zip
        touch outputs/alleles.zip
        touch outputs/cnv-calls.zip
        touch outputs/copy-number.zip
        touch outputs/read-depth.zip
        touch outputs/phased-variants.zip

        mv outputs/results.zip results.~{gene_symbol}.zip
        mv outputs/imported-variants.zip imported-variants.~{gene_symbol}.zip
        mv outputs/consolidated-variants.zip consolidated-variants.~{gene_symbol}.zip
        mv outputs/genotypes.zip genotypes.~{gene_symbol}.zip
        mv outputs/phenotypes.zip phenotypes.~{gene_symbol}.zip
        mv outputs/alleles.zip alleles.~{gene_symbol}.zip
        mv outputs/cnv-calls.zip cnv-calls.~{gene_symbol}.zip
        mv outputs/copy-number.zip copy-number.~{gene_symbol}.zip
        mv outputs/read-depth.zip read-depth.~{gene_symbol}.zip
        mv outputs/phased-variants.zip phased-variants.~{gene_symbol}.zip

        # Copy the png images in the subdirs only of the subdir copy-number-profile exists else touch an empty file
        if [ -d outputs/copy-number-profile ]; then
            src=$(ls outputs/copy-number-profile/*.png | head -n 1)
            cp "$src" "copy-number-profile.~{gene_symbol}.png"
        else
            touch "copy-number-profile.~{gene_symbol}.png"
        fi

        if [ -d outputs/allele-fraction-profile ]; then
            src=$(ls outputs/allele-fraction-profile/*.png | head -n 1)
            cp "$src" "allele-fraction-profile.~{gene_symbol}.png"
        else
            touch "allele-fraction-profile.~{gene_symbol}.png"
        fi

    >>>

    output {
        File alleles = "alleles.~{gene_symbol}.zip"
        File cnv_calls = "cnv-calls.~{gene_symbol}.zip"
        File consolidated_variants = "consolidated-variants.~{gene_symbol}.zip"
        File copy_number = "copy-number.~{gene_symbol}.zip"
        File genotypes = "genotypes.~{gene_symbol}.zip"
        File imported_variants = "imported-variants.~{gene_symbol}.zip"
        File phased_variants = "phased-variants.~{gene_symbol}.zip"
        File phenotypes = "phenotypes.~{gene_symbol}.zip"
        File read_depth = "read-depth.~{gene_symbol}.zip"
        File results = "results.~{gene_symbol}.zip"
        File allele_fraction_profile = "allele-fraction-profile.~{gene_symbol}.png"
        File copy_number_profile = "copy-number-profile.~{gene_symbol}.png"
    }

  runtime {
    memory: "16 GB"
    cpu: 2
    disks: "local-disk 10 HDD"
    docker: docker
    preemptible: preemptible_tries
    noAddress: no_address_override
  }
}


task ConcatResults {
    input {
        String base_file_name
        Array[File] result_files
        Array[String] gene_symbols
        File monitoring_script
        Boolean no_address_override = true
        String docker
        Int preemptible_tries
    }

    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # Write header line
        echo -e "Symbol\tSampleName\tGenotype\tPhenotype\tHaplotype1\tHaplotype2\tAlternativePhase\tVariantData\tCNV" > ~{base_file_name}.pypgx.tsv

        result_files_path=~{write_lines(result_files)}
        # Concat the contents of the results, with the gene symbol as the first column
        for gene_symbol in ~{sep=' ' gene_symbols}; do
            src=$(grep -w "$gene_symbol" $result_files_path)
            pypgx print-data "$src" | tail -1 | awk -v gs="$gene_symbol" '{print gs"\t"$0}' >> ~{base_file_name}.pypgx.tsv
        done

    >>>

    output {
        File results = "~{base_file_name}.pypgx.tsv"
    }
  runtime {
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 10 HDD"
    docker: docker
    preemptible: preemptible_tries
    noAddress: no_address_override
  }

}