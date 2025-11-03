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
#   This workflow align reads to a reference using the Ultima Alignment (UA) algorithm.
# CHANGELOG

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/alignment_tasks.wdl" as AlignTasks

workflow UAAlignment {
    input {
        Array[File] input_files
        String base_file_name
        UaParameters ua_parameters
        References? references
        File cache_tarball
        Boolean no_address
        Int preemptible_tries

        #@wv not defined(ua_parameters['ua_index']) -> defined(references)
        #@wv suffix(ua_parameters['ref_alt']) == '.alt'
    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    if (! defined(ua_parameters.ua_index)){
        call AlignTasks.BuildUaIndex {
            input :
                references          = select_first([references]),
                preemptible_tries   = preemptible_tries,
                ua_docker           = global.ua_docker,
                monitoring_script   = global.monitoring_script, # !FileCoercion
                no_address          = no_address,
                dummy_input_for_call_caching = "",
        }
    }

    File ua_index = select_first([ua_parameters.ua_index, BuildUaIndex.ua_index])

    call AlignTasks.AlignWithUA {
        input :
            input_bams              = input_files,
            output_bam_basename     = base_file_name + ".ua.aln",
            ua_index                = ua_index,                   # pre-built in GCP
            ref_alt                 = ua_parameters.ref_alt,
            extra_args              = ua_parameters.ua_extra_args,
            use_v_aware_alignment   = ua_parameters.v_aware_alignment_flag,
            no_address              = no_address,
            cache_tarball           = cache_tarball,
            monitoring_script       = global.monitoring_script, # !FileCoercion
            preemptible_tries       = preemptible_tries,
            ua_docker               = global.ua_docker,
            cpu                     = ua_parameters.cpus,
            memory_gb               = ua_parameters.memory_gb,
    }
    output {
        Array[File] ua_output_json = AlignWithUA.ua_output_json
        File ua_output_bam         = AlignWithUA.ua_output_bam
        File? ua_index_build       = BuildUaIndex.ua_index
    }
}