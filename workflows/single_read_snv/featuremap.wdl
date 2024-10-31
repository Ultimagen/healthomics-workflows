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
#

# DESCRIPTION
# Runs FeatureMap generation - aggregates all the substitutions from a bam file

# CHANGELOG in reverse chronological order
# 1.7.0 FeatureMapProcess now adds annotations to FeatureMap (support for balanced strand annotations) + merged dataframe no longer created
# 1.4.0 BIOIN-870 Featuremap flow updates

import "tasks/structs.wdl" as Structs
import "tasks/globals.wdl" as Globals
import "tasks/general_tasks.wdl" as UGGeneralTasks


workflow FeatureMap {
input {
  String pipeline_version = "1.5.0"  # !UnusedDeclaration
  File input_cram_bam
  File input_cram_bam_index
  References references

  File wgs_calling_interval_list
  Int break_bands_at_multiples_of
  FeatureMapParams featuremap_params
  Int scatter_count

  String flow_order
  String? ppmSeq_adapter_version

  String base_file_name
  Int? preemptible_tries
  Int? process_featuremap_memory_gb_override
  # winval validations
  #@wv not(" " in base_file_name or "#" in base_file_name or ',' in base_file_name)
}
  Int preemptibles = select_first([preemptible_tries, 1])
  Int ref_fasta_size = ceil(size(references.ref_fasta, "GB"))
  Float input_cram_bam_size = size(input_cram_bam, "GB")
  Int dynamic_featuremap_disk_size = ceil(input_cram_bam_size) + 50
  Int secure_disk_size_threshold = 512
  Int featuremap_disk_size = if dynamic_featuremap_disk_size > secure_disk_size_threshold then dynamic_featuremap_disk_size else secure_disk_size_threshold
  Int process_featuremap_memory_gb = select_first([process_featuremap_memory_gb_override, 16])

  call Globals.Globals as Globals
  GlobalVariables global = Globals.global_dockers

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call UGGeneralTasks.ScatterIntervalList {
    input:
      interval_list = wgs_calling_interval_list,
      scatter_count = scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      docker = global.gitc_docker,
      gitc_path = global.gitc_jar_path,
      no_address = true,
      dummy_input_for_call_caching = "",
      monitoring_script = global.monitoring_script   #!FileCoercion
  }


  # Create featuremap in parallel over WGS calling intervals, then merge
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate VCF by interval
    Int create_featuremap_memory_gb = 4
    call FeatureMapCreate {
      input:
        input_cram_bam = input_cram_bam,
        input_cram_bam_index = input_cram_bam_index,
        references = references,
        featuremap_params = featuremap_params,
        output_basename = base_file_name,
        interval_list = ScatterIntervalList.out[index],
        memory_gb = create_featuremap_memory_gb,
        disk_size = ceil(2 * featuremap_disk_size / ScatterIntervalList.interval_count) + ceil(size(references.ref_fasta, "GB")) + 2,
        docker = global.ug_gatk_picard_docker,
        gitc_path = global.gitc_jar_path,
        preemptibles = preemptibles,
        monitoring_script = global.monitoring_script  #!FileCoercion
    }
  
  }

  Int merge_featuremap_disk_size_gb = (ceil(size(FeatureMapCreate.featuremap,"GB")) * 7) + 3
  Int merge_featuremap_memory_gb = 4
  Int merge_featuremap_cpus = 4

  call FeatureMapMerge {
    input:
      featuremap_parts = FeatureMapCreate.featuremap,
      featuremap_parts_indices = FeatureMapCreate.featuremap_index,
      output_basename = base_file_name,
      tag = "featuremap_no_annots",
      disk_size = merge_featuremap_disk_size_gb,
      memory_gb = merge_featuremap_memory_gb,
      cpus = merge_featuremap_cpus,
      docker = global.ug_vc_docker,
      preemptibles = preemptibles,
      monitoring_script = global.monitoring_script  #!FileCoercion
  }

  Float featuremap_size_gb = size(FeatureMapMerge.featuremap, "GB")
  Int process_featuremap_disk_size_gb = ceil(featuremap_size_gb * 8) + ref_fasta_size
  Int process_featuremap_cpus = 8

  call FeatureMapAnnotations {
    input:
      featuremap = FeatureMapMerge.featuremap,
      featuremap_index = FeatureMapMerge.featuremap_index,
      references = references,
      featuremap_params = featuremap_params,
      flow_order = flow_order,
      output_basename = base_file_name,
      ppmSeq_adapter_version = ppmSeq_adapter_version,
      docker = global.ug_vc_docker,
      memory_gb = process_featuremap_memory_gb,
      cpus = process_featuremap_cpus,
      disk_size = process_featuremap_disk_size_gb,
      preemptibles = preemptibles,
      monitoring_script = global.monitoring_script  #!FileCoercion
  }

  call FeatureMapSingleSubstitutions {
    input:
      featuremap = FeatureMapAnnotations.annotated_featuremap,
      featuremap_index = FeatureMapAnnotations.annotated_featuremap_index,
      output_basename = base_file_name,
      docker = global.vcflite_docker,
      memory_gb = 4,
      cpus = 2,
      preemptibles = preemptibles,
      monitoring_script = global.monitoring_script  #!FileCoercion
  }

  output {
    File featuremap = FeatureMapAnnotations.annotated_featuremap
    File featuremap_index = FeatureMapAnnotations.annotated_featuremap_index
    File featuremap_single_substitutions = FeatureMapSingleSubstitutions.featuremap_single_substitutions
    File featuremap_single_substitutions_index = FeatureMapSingleSubstitutions.featuremap_single_substitutions_index
  }
}

task FeatureMapCreate {
  input {
    File input_cram_bam
    File input_cram_bam_index
    References references
    FeatureMapParams featuremap_params
    String output_basename
    File interval_list
    Int memory_gb
    Int disk_size
    String docker
    String gitc_path
    Int preemptibles
    File monitoring_script
  }

  parameter_meta {
    input_cram_bam: {
      localization_optional: true
    }
  }
  String interval_list_basename = basename(interval_list, ".interval_list")
  command <<<
    set -xeo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    source ~/.bashrc
    start=$(date +%s)

    echo "***************************** Running FeatureMap *****************************"
    start_task=$(date +%s)
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{memory_gb-2}g -jar ~{gitc_path}GATK_ultima.jar  \
      FlowFeatureMapper \
      -I "~{input_cram_bam}" \
      -O tmp.vcf.gz \
      -R "~{references.ref_fasta}"  \
      --intervals "~{interval_list}" \
      --snv-identical-bases ~{featuremap_params.snv_identical_bases} \
      --snv-identical-bases-after ~{featuremap_params.snv_identical_bases_after} \
      --min-score ~{featuremap_params.min_score} \
      --limit-score ~{featuremap_params.limit_score} \
      --read-filter MappingQualityReadFilter --minimum-mapping-quality ~{featuremap_params.min_mapq} \
      --threaded-walker --threaded-writer \
      ~{featuremap_params.extra_args} 

    echo "***************************** Filtering on region *****************************"
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms~{memory_gb-2}g -jar ~{gitc_path}GATK_ultima.jar  \
      VariantFiltration \
      -V tmp.vcf.gz \
      -O "~{output_basename}.~{interval_list_basename}.vcf.gz" \
      -R "~{references.ref_fasta}"  \
      --intervals "~{interval_list}"
    echo "***************************** Done *****************************"

    end=$(date +%s)
    mins_elapsed=$(( ($end - $start_task) / 60))
    secs_elapsed=$(( ($end - $start_task) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"
    echo "***************************** Finished Running FeatureMap *****************************"
    ls -ltr

    end=$(date +%s)
    mins_elapsed=$(( ($end - $start) / 60))
    secs_elapsed=$(( ($end - $start) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"
  >>>
  runtime {
    preemptible: preemptibles
    cpu: 2
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    maxRetries: 0
  }
  output{
    File monitoring_log = "monitoring.log"
    File featuremap = "~{output_basename}.~{interval_list_basename}.vcf.gz"
    File featuremap_index = "~{output_basename}.~{interval_list_basename}.vcf.gz.tbi"
  }
}

task FeatureMapAnnotations {
  input  {
    File featuremap
    File featuremap_index
    FeatureMapParams featuremap_params
    References references
    String flow_order
    String output_basename
    String? ppmSeq_adapter_version
    String docker
    Int memory_gb
    Int cpus
    Int disk_size
    Int preemptibles
    File monitoring_script
  }
  
  String annotated_featuremap_vcf = basename(featuremap, ".vcf.gz") + ".annotated.vcf.gz"

  command <<<
    set -xeo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    source ~/.bashrc
    conda activate genomics.py3

    echo "***************************** Add additional features to featuremap *****************************"
    start=$(date +%s)
    start_task=$(date +%s)
    python /VariantCalling/ugvc annotate_featuremap \
      -i ~{featuremap} \
      -o ~{annotated_featuremap_vcf} \
      --ref_fasta ~{references.ref_fasta} \
      ~{true="--ppmSeq_adapter_version " false="" defined(ppmSeq_adapter_version)}~{default="" ppmSeq_adapter_version} \
      --flow_order ~{flow_order} \
      --motif_length_to_annotate ~{featuremap_params.motif_length_to_annotate} \
      --max_hmer_length ~{featuremap_params.max_hmer_length} \
      --process_number ~{cpus}

    end=$(date +%s)
    mins_elapsed=$(( ($end - $start_task) / 60))
    secs_elapsed=$(( ($end - $start_task) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Run time: $mins_elapsed:$secs_elapsed"
    echo "***************************** Done *****************************"
  >>>
  runtime {
    preemptible: preemptibles
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File annotated_featuremap = "~{annotated_featuremap_vcf}"
    File annotated_featuremap_index = "~{annotated_featuremap_vcf}.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task FeatureMapSingleSubstitutions {
  input  {
    File featuremap
    File featuremap_index
    String output_basename
    String docker
    Int memory_gb
    Int cpus
    # Int disk_size
    Int preemptibles
    File monitoring_script
  }
  
  String single_sub_output = "~{output_basename}.single_substitution.vcf.gz"
  Int disk_size = ceil(size(featuremap, "GB") * 30) + 20

  command <<<
    set -xeo pipefail
    bash ~{monitoring_script} | tee monitoring.log >&2 &

    echo "$(date) ***************************** Importing vcf into SQLite database *****************************"

    vcflite import --vcf-in ~{featuremap} 

    echo "$(date) ***************************** Generating FeatureMap of single substitutions *****************************"

    vcflite query --group-by "chrom, pos" --having "COUNT(*) = 1" --vcf-out ~{single_sub_output}

    echo "$(date) ***************************** Indexing FeatureMap of single substitutions *****************************"
    tabix -p vcf ~{single_sub_output}

    echo "$(date) ***************************** Done *****************************"
  >>>
  runtime {
    preemptible: preemptibles
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
  }
  output{
    File featuremap_single_substitutions = "~{single_sub_output}"
    File featuremap_single_substitutions_index = "~{single_sub_output}.tbi"
    File monitoring_log = "monitoring.log"
  }
}

task FeatureMapMerge {
  input  {
    Array[File] featuremap_parts
    Array[File] featuremap_parts_indices
    String output_basename
    String tag
    Int memory_gb
    Int cpus
    Int disk_size
    String docker
    Int preemptibles
    File monitoring_script
  }

  command <<<
    set -eo pipefail
    set -x
    bash ~{monitoring_script} | tee monitoring.log >&2 &
    
    source ~/.bashrc
    start=$(date +%s)

    conda activate genomics.py3
    mkdir -p tmp
    # bypass for architectures where index files are not synced with vcf files (such as local miniwdl)
    for part in ~{sep=" " featuremap_parts}; do if ! test -f $part.tbi; then bcftools index -t $part; fi; done

    bcftools concat --threads ~{cpus} -a -Oz -o ~{output_basename}.~{tag}.vcf.gz ~{sep=" " featuremap_parts}
    bcftools index --threads ~{cpus} -f -t ~{output_basename}.~{tag}.vcf.gz

    end=$(date +%s)
    mins_elapsed=$(( ($end - $start) / 60))
    secs_elapsed=$(( ($end - $start) % 60 ))
    if [ $secs_elapsed -lt 10 ]; then
      secs_elapsed=0$secs_elapsed
    fi
    echo "Total run time: $mins_elapsed:$secs_elapsed"
  >>>
  runtime {
    preemptible: preemptibles
    cpu: "~{cpus}"
    memory: "~{memory_gb} GB"
    disks: "local-disk " + ceil(disk_size) + " HDD"
    docker: docker
    maxRetries: 0
  }
  output{
    File monitoring_log = "monitoring.log"
    File featuremap = "~{output_basename}.~{tag}.vcf.gz"
    File featuremap_index = "~{output_basename}.~{tag}.vcf.gz.tbi"
  }
}
