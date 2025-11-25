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
#   This workflow align reads to a reference using the Ultima Alignment (UA) for methyl-seq algorithm.
# CHANGELOG

import "tasks/globals.wdl" as Globals
import "tasks/structs.wdl" as Structs
import "tasks/alignment_tasks.wdl" as AlignTasks

workflow UAMethAlignment {
    input {
        Array[File] input_files
        String base_file_name
        UaMethParameters ua_meth_parameters
        File cache_tarball
        References? references
        Boolean UaMethIntensiveMode = false  # If false, map against the three-letter C2T index, otherwise  map against both C2T and G2A indexes.

        Boolean no_address
        Int preemptible_tries

        #@wv not defined(ua_meth_parameters) -> defined(references)
    }
    call Globals.Globals as Globals
    GlobalVariables global = Globals.global_dockers

    Int cpu_default = 40 
    if (defined(ua_meth_parameters.index_c2t)) {
        File ua_index_c2t_ = select_first([ua_meth_parameters.index_c2t])
        File ua_index_g2a_ = select_first([ua_meth_parameters.index_g2a])
        
    } 
    if (!defined(ua_meth_parameters.index_c2t)) {
        call AlignTasks.BuildUaMethIndex {
            input :
                references          = select_first([references]),
                preemptible_tries   = preemptible_tries,
                ua_docker           = global.ua_docker,
                monitoring_script   = global.monitoring_script, # !FileCoercion
                no_address          = no_address,
        }
    }

    File ua_index_c2t = select_first([ua_index_c2t_, BuildUaMethIndex.index_c2t])
    File ua_index_g2a = select_first([ua_index_g2a_, BuildUaMethIndex.index_g2a])

    Int cpu = select_first([ua_meth_parameters.cpus, cpu_default])

    call AlignTasks.AlignWithUAMeth {
        input :
            input_bams              = input_files,
            cache_tarball           = cache_tarball,
            output_bam_basename     = base_file_name + ".ua.aln",
            index_c2t               = ua_index_c2t,
            index_g2a               = ua_index_g2a,
            ref_alt                 = ua_meth_parameters.ref_alt,
            extra_args              = ua_meth_parameters.ua_extra_args,
            use_v_aware_alignment   = ua_meth_parameters.v_aware_alignment_flag,
            UaMethIntensiveMode     = UaMethIntensiveMode,
            no_address              = no_address,
            monitoring_script       = global.monitoring_script,  # !FileCoercion
            preemptible_tries       = preemptible_tries,
            ua_docker               = global.ua_docker,
            memory_gb               = ua_meth_parameters.memory_gb,
            cpu                     = cpu
    }
    output {
        Array[File] ua_output_json = AlignWithUAMeth.ua_output_json
        File ua_output_bam         = AlignWithUAMeth.ua_output_bam
        File? ua_index_c2t_build   = BuildUaMethIndex.index_c2t
        File? ua_index_g2a_build   = BuildUaMethIndex.index_g2a
    }
}