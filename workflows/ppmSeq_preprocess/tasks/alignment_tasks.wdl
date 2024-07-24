version 1.0
import "structs.wdl"

# Note regarding dummy_input_for_call_caching
# When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your
# own workspace. This will prevent call caching from failing with "Cache Miss (10 failed copy attempts)". Outside of
# Terra this can be left as the default empty String. This dummy input is only needed for tasks that have no inputs
# specific to the sample being run (such as GetBwaVersion which does not take in any sample data).


# TASK DEFINITIONS
task SplitCram {
    input {
        File monitoring_script
        File input_cram_bam
        String base_file_name
        Int reads_per_file
        String docker
        Int disk_size = ceil(3 * size(input_cram_bam,"GB"))
        Int preemptible_tries
        Boolean no_address
    }

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        mkdir -p splitout
        /crammer/crammer --split --out splitout/~{base_file_name}-%d.cram \
        --progress --nreads-per-file ~{reads_per_file}  < ~{input_cram_bam}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " LOCAL"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        Array[File] split_outputs = glob("splitout/*.cram")
        File monitoring_log = "monitoring.log"
    }
}

task SamSplitter {
    input {
        File monitoring_script
        File input_bam
        Int n_reads
        Int preemptible_tries
        Int compression_level
        String docker
        Boolean no_address
    }
    Float unmapped_bam_size = size(input_bam, "GB")
    # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
    Float disk_multiplier = 2.5
    Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        mkdir output_dir

        total_reads=$(samtools view -c ~{input_bam})

        java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
        INPUT=~{input_bam} \
        OUTPUT=output_dir \
        SPLIT_TO_N_READS=~{n_reads} \
        TOTAL_READS_IN_INPUT=$total_reads
    >>>
    output {
        Array[File] split_outputs = glob("output_dir/*.bam")
        File monitoring_log = "monitoring.log"
    }
    runtime {
        docker: docker
        preemptible: preemptible_tries
        memory: "3.75 GB"
        disks: "local-disk " + disk_size + " HDD"
        noAddress: no_address
        maxRetries: 1
    }
}


task CreateReferenceCache {
    input {
        Array[File] references
        File cache_populate_script
        Int preemptible_tries
        Int disk_size = ceil(3 * size(references, "GB") + 20)

        String docker
        String dummy_input_for_call_caching # !UnusedDeclaration
    }

    command <<<
        set -xeo pipefail
        mkdir cache

        for file in ~{sep=' ' references}; do
            perl ~{cache_populate_script} -root cache ${file}
        done

        tar -zcf cache.tgz cache
    >>>

    output {
        File cache_tarball = "cache.tgz"
    }

    runtime {
        preemptible: preemptible_tries
        memory: "4 GB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        maxRetries: 1
    }
}


# Convert recalibrated cram to bam
task ConvertCramOrBamToUBam {
    input {
        File monitoring_script
        File input_file
        File cache_tarball
        String base_file_name
        String docker
        Int disk_size = ceil((8 * size(input_file,"GB")) + size(cache_tarball,"GB") + 20)
        Int preemptible_tries
        Boolean no_address
    }

    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &


        samtools view -H ~{input_file} | grep "^@RG" | \
            awk '{for (i=1;i<=NF;i++){if ($i ~/tp:/) {print substr($i,4)}}}' | head -1 \
            > tp_direction.txt

        a=`cat tp_direction.txt`

        if [ "$a" = "reference" ];
        then
            tp_flag="--ATTRIBUTE_TO_REVERSE tp --ATTRIBUTE_TO_REVERSE t0"
            sort_order_flag="--SO queryname"
        else
            tp_flag=""
            sort_order_flag="--SO unsorted"
        fi

        tar -zxf ~{cache_tarball}

        REF_CACHE=cache/%2s/%2s/ REF_PATH='.' samtools view -b -F 2048 -h ~{input_file} -o tmp.bam 
        java -Xmx11g -jar /usr/gitc/picard.jar RevertSam -I tmp.bam \
            -O ~{base_file_name}.u.bam \
            --MAX_DISCARD_FRACTION 0.005 \
            --ATTRIBUTE_TO_CLEAR XT \
            --ATTRIBUTE_TO_CLEAR XN \
            --ATTRIBUTE_TO_CLEAR AS \
            --ATTRIBUTE_TO_CLEAR OC \
            --ATTRIBUTE_TO_CLEAR OP \
            $tp_flag \
            --REMOVE_DUPLICATE_INFORMATION \
            --REMOVE_ALIGNMENT_INFORMATION \
            --VALIDATION_STRINGENCY LENIENT \
            $sort_order_flag
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "13 GB"
        cpu: "3"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File unmapped_bam = "~{base_file_name}.u.bam"
        File monitoring_log = "monitoring.log"
    }
}


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
    input {
        File input_bam
        String output_bam_basename
        AlignmentReferences bwa_references
        File monitoring_script
        Boolean no_address
        # The merged bam can be bigger than only the aligned bam,
        # so account for the output size by multiplying the input size by 3.5.
        Int preemptible_tries
        String docker
    }
    Int disk_size = ceil(3.5*size(input_bam,"GB") +
        size(bwa_references.references.ref_fasta,"GB") +
        size(bwa_references.ref_alt,"GB") +
        size(bwa_references.ref_amb,"GB") +
        size(bwa_references.ref_ann,"GB") +
        size(bwa_references.ref_alt,"GB") +
        size(bwa_references.ref_bwt,"GB") +
        size(bwa_references.ref_pac,"GB") +
        size(bwa_references.ref_sa,"GB") +
        20)
    command <<<
        /usr/gitc/bwa 2>bwa_help  # done before "set -o pipefail" because /bwa has a rc=1 and we dont want to allow rc=1
        set -eo pipefail
        export BWA_VERSION=$(grep -e '^Version' bwa_help | sed 's/Version: //')

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # set the bash variable needed for the command-line
        bash_ref_fasta=~{bwa_references.references.ref_fasta}
        bwa_commandline=$( echo "/usr/gitc/bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta")

        # if ref_alt has data in it,
        java -Xms5000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
        $bwa_commandline /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
        java -Xms3000m -jar /usr/gitc/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_RETAIN=tm \
            ATTRIBUTES_TO_RETAIN=tf \
            ATTRIBUTES_TO_RETAIN=RX \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ATTRIBUTES_TO_REVERSE=ti \
            ATTRIBUTES_TO_REVERSE=tp \
            ATTRIBUTES_TO_REVERSE=t0 \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=~{input_bam} \
            OUTPUT=~{output_bam_basename}.bam \
            REFERENCE_SEQUENCE=~{bwa_references.references.ref_fasta} \
            SORT_ORDER="queryname" \
            IS_BISULFITE_SEQUENCE=false \
            CLIP_ADAPTERS=false \
            ALIGNED_READS_ONLY=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwamem" \
            PROGRAM_GROUP_VERSION="$BWA_VERSION" \
            PROGRAM_GROUP_COMMAND_LINE="$bwa_commandline" \
            PROGRAM_GROUP_NAME="bwamem" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=false \
            ADD_PG_TAG_TO_READS=false

        echo "Piped SAM->FASTQ->BWA->MergeBamAlignment complete."

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "28 GB"
        cpu: "16"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File output_bam = "~{output_bam_basename}.bam"
        File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
        File monitoring_log = "monitoring.log"
    }
}


# Performs bwa-meth alignment for methyl-seq data
task SamToFastqAndBwaMeth {
        input {
            File input_bam
            String output_bam_basename
            File bwa_meth_ref # tar file
            File monitoring_script
            Boolean no_address
            Int preemptible_tries
            String docker
        }

        Int disk_size = ceil(5*size(input_bam,"GB") + (8 * size(bwa_meth_ref,"GB")) )


     command <<<
        /usr/gitc/bwa 2>bwa_help  # done before "set -o pipefail" because /bwa has a rc=1 and we dont want to allow rc=1
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        bwameth_version=$(bwameth.py --version)

        tar xvf ~{bwa_meth_ref}
        ref_basename_files=$(ls | grep -P "(.fa$)|(.fasta$)")
        bwa_meth_commandline=$( echo "bwameth.py -t 16 --reference $ref_basename_files")
        rm ~{bwa_meth_ref}

        picard \
         SamToFastq \
         INPUT=~{input_bam} \
         FASTQ=/dev/stdout \
         INTERLEAVE=true \
         NON_PF=true | \
        $bwa_meth_commandline /dev/stdin 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
        picard \
         SortSam -I /dev/stdin  -O /dev/stdout -SO queryname | \
        picard \
         MergeBamAlignment \
         VALIDATION_STRINGENCY=SILENT \
         EXPECTED_ORIENTATIONS=FR \
         ATTRIBUTES_TO_RETAIN=X0 \
         ATTRIBUTES_TO_RETAIN=tm \
         ATTRIBUTES_TO_RETAIN=tf \
         ATTRIBUTES_TO_RETAIN=RX \
         ATTRIBUTES_TO_REMOVE=NM \
         ATTRIBUTES_TO_REMOVE=MD \
         ATTRIBUTES_TO_REMOVE=RG \
         ATTRIBUTES_TO_REVERSE=ti \
         ATTRIBUTES_TO_REVERSE=tp \
         ATTRIBUTES_TO_REVERSE=t0 \
         ALIGNED_BAM=/dev/stdin \
         UNMAPPED_BAM=~{input_bam} \
         OUTPUT=~{output_bam_basename}.bam \
         REFERENCE_SEQUENCE="$ref_basename_files" \
         SORT_ORDER="queryname" \
         IS_BISULFITE_SEQUENCE=true \
         CLIP_ADAPTERS=false \
         ALIGNED_READS_ONLY=false \
         MAX_RECORDS_IN_RAM=2000000 \
         ADD_MATE_CIGAR=true \
         MAX_INSERTIONS_OR_DELETIONS=-1 \
         PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
         PROGRAM_RECORD_ID="bwameth" \
         PROGRAM_GROUP_VERSION="$bwameth_version" \
         PROGRAM_GROUP_COMMAND_LINE="$bwa_meth_commandline" \
         PROGRAM_GROUP_NAME="bwameth" \
         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
         ALIGNER_PROPER_PAIR_FLAGS=true \
         UNMAP_CONTAMINANT_READS=false \
         ADD_PG_TAG_TO_READS=false

        echo "Piped SAM->FASTQ->BWA_METH->MergeBamAlignment complete."

    >>>

    runtime {
        preemptible: preemptible_tries
        memory: "32 GB"
        cpu: "25"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File output_bam = "~{output_bam_basename}.bam"
        File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
        File monitoring_log = "monitoring.log"
    }
}

## Constructs the index for UA. We keep this task in the pipeline because UA index format keeps changing
task BuildUaIndex{
    input {
      References references
      File monitoring_script
      Boolean no_address
      Int disk_size = ceil((15 * size(references.ref_fasta,"GB")) + 20)
      Int preemptible_tries
      String ua_docker
      String dummy_input_for_call_caching # !UnusedDeclaration
    }

    String output_file = basename(references.ref_fasta) + ".uai"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # for compatibility with the old image where ua was in /ua/ua and not in PATH
        export PATH=$PATH:/ua

        ua \
            --build \
            --ref ~{references.ref_fasta} \
            --seed 20,200,5 \
            --index ~{output_file} \
            --progress

    >>>

    runtime {
        cpu : "1"
        cpuPlatform: "Intel Skylake"
        preemptible: preemptible_tries
        memory: "200 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: ua_docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File ua_index =  "~{output_file}"
    }
}


### Constructs the index for UA methylation
task BuildUaMethIndex {
  input {
    References references
    File monitoring_script
    Boolean no_address
    Int disk_size = ceil((60 * size(references.ref_fasta,"GB")) + 60)
    Int preemptible_tries
    String ua_docker

    }

    String output_file = basename(references.ref_fasta) + ".uai"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        # for compatibility with the old image where ua was in /ua/ua and not in PATH
        export PATH=$PATH:/ua

        ua \
        --methylation \
        --build \
        --ref ~{references.ref_fasta} \
        --seed 20,200,5 \
        --index ~{output_file} \
        --progress

    >>>

    runtime {
        cpu : "1"
        preemptible: preemptible_tries
        memory: "200 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: ua_docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File index_c2t = "~{output_file}.c2t"
        File index_g2a = "~{output_file}.g2a"
        File monitoring_log = "monitoring.log"
    }
}


# Run UA alignment. Since UA is fast, we run on all chunks together
# use_v_aware_alignment flag toggles running variant aware ua
# For now we pass the v_aware_vcf separately as it might be needed for
# pileup calculation

task AlignWithUA {
  input {
        Array[File] input_bams
        File? cache_tarball
        String output_bam_basename
        File ua_index
        File ref_alt
        Boolean use_v_aware_alignment
        File? v_aware_vcf
        String? extra_args
        File monitoring_script
        Int disk_size = ceil(size(ua_index, "GB") + 3*size(input_bams, "GB") + 20 + size(v_aware_vcf,"GB") + 2*size(cache_tarball, "GB"))
        Int preemptible_tries
        Boolean no_address
        String ua_docker
        Int cpu = 40
    }

    Int preemptible_tries_final = if (size(input_bams, "GB") < 250) then preemptible_tries else 0
    command <<<
    set -eo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &
    echo "~{sep='\n'input_bams}" > bam_list.txt


    ~{"tar -zxf "+cache_tarball}
    
    # for compatibility with the old image where ua was in /ua/ua and not in PATH
    export PATH=$PATH:/ua

    REF_CACHE=cache/%2s/%2s/ REF_PATH='.' samtools cat -b bam_list.txt | \
    samtools view -h -@ ~{cpu} - | \
    ua \
        --index ~{ua_index} \
        --align true \
        --progress \
        --tp reference \
        --alt=~{ref_alt} \
        --json=~{output_bam_basename}.%s.json \
        --nthread max \
        ~{"--vcf="+ v_aware_vcf}  \
        ~{true="--vcf-snps-only --vcf-af-threshold=0.01" false='' use_v_aware_alignment} \
        --sam-input - \
        --sam-output - \
        ~{extra_args} | \
    samtools view -@ ~{cpu} -o ~{output_bam_basename}.bam -
    >>>

    runtime {
        cpuPlatform: "Intel Skylake"
        preemptible: preemptible_tries_final
        memory: "100 GB"
        cpu: "~{cpu}"
        disks: "local-disk " + disk_size + " HDD"
        docker: ua_docker
        noAddress: no_address
        maxRetries: 1
    }

    output {
        Array[File] ua_output_json = glob("*.json")
        File ua_output_bam = "~{output_bam_basename}.bam"
        File monitoring_log = "monitoring.log"
    }
}


# Run UA alignmentfor methyl-seq
task AlignWithUAMeth {
  input {
    Array[File] input_bams
    String output_bam_basename
    File index_c2t
    File index_g2a
    File? cache_tarball
    Boolean UaMethIntensiveMode = false
    File monitoring_script
    Int disk_size = ceil( size(index_c2t, "GB") + size(index_g2a, "GB") + 2.5*size(input_bams, "GB") + 20 + 2*size(cache_tarball, "GB"))
    Int preemptible_tries
    Boolean no_address
    String ua_docker
    Int cpus

    }
    Int preemptible_tries_final = if (size(input_bams, "GB") < 250) then preemptible_tries else 0
    String ua_meth_mode = if (UaMethIntensiveMode) then "--methylation-intensive" else "--methylation"

    command <<<
    set -eo pipefail

    bash ~{monitoring_script} | tee monitoring.log >&2 &
    echo "~{sep='\n'input_bams}" > bam_list.txt

    ua_index_file_name=$(echo ~{index_c2t} | sed -r 's/'.c2t'//g')
    ls -lstr
    cat bam_list.txt

    ~{"tar -zxf "+cache_tarball}

    # for compatibility with the old image where ua was in /ua/ua and not in PATH
    export PATH=$PATH:/ua

    REF_CACHE=cache/%2s/%2s/ REF_PATH='.' samtools cat -b bam_list.txt | \
    samtools view -h -@ ~{cpus} - | \
    ua \
        ~{ua_meth_mode} \
        --align true \
        --index "${ua_index_file_name}" \
        --progress \
        --tp reference \
        --vector \
        --json ~{output_bam_basename}-%s.json \
        --nthread ~{cpus} \
        --sam-input - \
        --sam-output - | \
    samtools view -@ ~{cpus} -o "~{output_bam_basename}.bam" -

    >>>

    runtime {
        cpuPlatform: "Intel Skylake"
        preemptible: preemptible_tries_final
        memory: "200 GB"
        cpu: "~{cpus}"
        disks: "local-disk " + disk_size + " HDD"
        docker: ua_docker
        noAddress: no_address
        maxRetries: 1
    }

    output {
        Array[File] ua_output_json = glob("*.json")
        File ua_output_bam = "~{output_bam_basename}.bam"
        File monitoring_log = "monitoring.log"
    }
}

# Read unmapped BAM, convert to FASTQ, align with vg giraffe and cobvert to bam with vg surject, then stream to MergeBamAlignment
task SamToFastqAndGiraffeAndMba {
    input {
        File input_bam
        String output_bam_basename
        GiraffeReferences giraffe_references
        String in_prefix_to_strip = "GRCh38."
        File monitoring_script
        Boolean no_address
        # The merged bam can be bigger than only the aligned bam,
        # so account for the output size by multiplying the input size by 3.5.
        Int preemptible_tries
        String docker
    }
    Int disk_size = ceil(4.5*size(input_bam,"GB") +
        size(giraffe_references.references.ref_fasta,"GB") +
        size(giraffe_references.ref_gbz,"GB") +
        size(giraffe_references.ref_dist,"GB") +
        size(giraffe_references.ref_min,"GB") +
        20)
    Int threads = 16
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms5000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=~{output_bam_basename}.fq \
        INTERLEAVE=true \
        NON_PF=true

        echo "Writing FASTQ complete."

        vg giraffe \
        --progress \
        --sample ~{output_bam_basename} \
        --output-format gaf \
        -f ~{output_bam_basename}.fq \
        -Z ~{giraffe_references.ref_gbz} \
        -d ~{giraffe_references.ref_dist} \
        -m ~{giraffe_references.ref_min} \
        -t ~{threads} | gzip > ~{output_bam_basename}.gaf.gz

        echo "Giraffe alignment complete."

        vg surject \
          -F ~{giraffe_references.ref_path_list} \
          -x ~{giraffe_references.ref_gbz} \
          -t ~{threads} \
          --bam-output --gaf-input \
          --sample ~{output_bam_basename} \
          --prune-low-cplx \
          ~{output_bam_basename}.gaf.gz > ~{output_bam_basename}.vg.unsorted.bam
        #| samtools sort -@ ~{threads} -n -O BAM -o ~{output_bam_basename}.vg.bam
        #/dev/stdin

        echo "Conversion GAF to BAM complete."

        if [ ~{in_prefix_to_strip} != "" ]
        then
            # patch the SQ fields from the dict into a new header
            samtools view -H ~{output_bam_basename}.vg.unsorted.bam | grep ^@HD > new_header.sam
            grep ^@SQ ~{giraffe_references.references.ref_dict} | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
            samtools view -H ~{output_bam_basename}.vg.unsorted.bam  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam

            cat <(cat new_header.sam) <(samtools view ~{output_bam_basename}.vg.unsorted.bam) | \
                sed -e "s/~{in_prefix_to_strip}//g"  | \
                samtools sort --threads ~{threads} -n -O BAM -o ~{output_bam_basename}.vg.bam
        else
            samtools sort --threads ~{threads} ~{output_bam_basename}.vg.unsorted.bam -O BAM -o ~{output_bam_basename}.vg.bam

        fi

        echo "Fixing header + qsort complete"

        java -Xms3000m -jar /usr/gitc/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_RETAIN=tm \
            ATTRIBUTES_TO_RETAIN=tf \
            ATTRIBUTES_TO_RETAIN=RX \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ATTRIBUTES_TO_REVERSE=ti \
            ATTRIBUTES_TO_REVERSE=tp \
            ATTRIBUTES_TO_REVERSE=t0 \
            ALIGNED_BAM=~{output_bam_basename}.vg.bam \
            UNMAPPED_BAM=~{input_bam} \
            OUTPUT=~{output_bam_basename}.bam \
            REFERENCE_SEQUENCE=~{giraffe_references.references.ref_fasta} \
            SORT_ORDER="queryname" \
            IS_BISULFITE_SEQUENCE=false \
            CLIP_ADAPTERS=false \
            ALIGNED_READS_ONLY=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="vg_giraffe" \
            PROGRAM_GROUP_VERSION="v1.45.0-105-gd213c299c" \
            PROGRAM_GROUP_COMMAND_LINE="vg_commandline" \
            PROGRAM_GROUP_NAME="vg_giraffe" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=false \
            ADD_PG_TAG_TO_READS=false

        echo "Merge VG aligned and Ubam complete."

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "80 GB"
        cpu: "16"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        noAddress: no_address
        maxRetries: 1
    }
    output {
        File output_bam = "~{output_bam_basename}.bam"
        File monitoring_log = "monitoring.log"
    }
}
task MarkDuplicatesSpark {
    input {
        Array[File] input_bams
        String output_bam_basename
        File monitoring_script
        String docker
        Boolean no_address
        String gitc_path = "/usr/gitc/"
        String? args
    }
    parameter_meta {
        input_bams: {
            localization_optional: true
        }
    }

    Int md_disk_size = 6 * ceil(size(input_bams, "GB"))
    Int local_ssd_size = 375
    # add one local ssd disk (375 GB) when spare disk is < 50 GB
    # the actual size that will be required to google will be a multiple of 375
    Int rounded_disk_size = ceil(md_disk_size/local_ssd_size) * local_ssd_size
    Int total_md_disk_size = if (rounded_disk_size - md_disk_size) < 50 then md_disk_size + local_ssd_size else md_disk_size

    Int local_ssd_threshold = 3000
    Int mapped_bam_size_local_ssd = if total_md_disk_size < local_ssd_threshold then total_md_disk_size else 9000

    # increase memory size for large inputs
    Int increased_ram_size = 300
    Int mark_duplicates_memory_threshold = 600
    Int ram_size = if md_disk_size / 4 > mark_duplicates_memory_threshold then increased_ram_size else 208

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        bams_dirname=$(echo "~{sep='\n'input_bams}" | tail -1 | xargs dirname)

        java -Xmx190g -jar ~{gitc_path}GATK_ultima.jar MarkDuplicatesSpark \
        --spark-master local[24] \
        --conf spark.driver.memory=6g \
        --conf spark.executor.memory=5g \
        --spark-verbosity WARN \
        --input ~{sep=" --input " input_bams} \
        --output ~{output_bam_basename}.bam \
        --create-output-bam-index true \
        --verbosity WARNING \
        ~{args}
    >>>

    runtime {
        disks: "local-disk " + mapped_bam_size_local_ssd + " LOCAL"
        cpu: 32
        memory: ram_size + " GB"
        preemptible: 0
        docker: docker
        noAddress: no_address
    }
    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_index = "~{output_bam_basename}.bam.bai"
        File monitoring_log = "monitoring.log"
    }
}


# Unlike MarkDuplicatesSpark, MarkDuplicates task works only on a single file. It doesn't use parallelization, and requires a smaller machine.
task MarkDuplicates {
    input {
        File input_bam
        String base_file_name
        References references
        File monitoring_script
        String docker
        Boolean no_address
        Int preemptibles
        Int disk_size = ceil(3.5*size(input_bam,"GB") +
            size(references.ref_fasta,"GB") + 10)
    }
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xmx6g -jar /usr/gitc/picard.jar \
          MarkDuplicates \
          -I ~{input_bam} \
          -O ~{base_file_name}.markdup.bam \
          -M ~{base_file_name}.duplicate_metrics \
          -R ~{references.ref_fasta} \
          --VERBOSITY WARNING \
          --FLOW_MODE True
    >>>
    runtime {
        disks: "local-disk " + disk_size + " HDD"
        cpu: 2
        memory: "8 GB"
        preemptible: 1
        docker: docker
        noAddress: no_address
    }
    output {
        File output_bam = "~{base_file_name}.markdup.bam"
        File output_bam_index = "~{base_file_name}.markdup.bam.bai"
        File duplicate_metrics = "~{base_file_name}.duplicate_metrics"
        File monitoring_log = "monitoring.log"
    }
}


# Convert BAM file to CRAM format
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
    input {
        File input_bam
        References references
        String output_basename
        Int preemptible_tries
        String docker
        File monitoring_script
    }

    Float dynamic_convert_to_cram_disk_size = 2 * size(input_bam, "GB") + size(references.ref_fasta,"GB") + 20
    Float secure_disk_size_threshold = 510.0

    Float convert_to_cram_disk_size = if dynamic_convert_to_cram_disk_size > secure_disk_size_threshold then dynamic_convert_to_cram_disk_size else secure_disk_size_threshold
    Int disk_size = ceil(convert_to_cram_disk_size)

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        samtools view -C -T ~{references.ref_fasta} ~{input_bam} | \
        tee ~{output_basename}.cram | \
        md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

        # Create REF_CACHE. Used when indexing a CRAM
        seq_cache_populate.pl -root ./ref/cache ~{references.ref_fasta}
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s

        samtools index ~{output_basename}.cram
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        maxRetries: 1
    }
    output {
        File output_cram = "~{output_basename}.cram"
        File output_cram_index = "~{output_basename}.cram.crai"
        File output_cram_md5 = "~{output_basename}.cram.md5"
        File monitoring_log = "monitoring.log"
    }
}


task ValidateSamFile {
    input {
        File input_bam
        File? input_bam_index
        String report_filename
        References references
        Int? max_output
        Array[String]? ignore
        Int preemptible_tries
        String docker
        String gitc_path = "/usr/gitc/"
        Boolean is_methyl_seq
        File monitoring_script
    }
    Float secure_disk_size_threshold = 510.0

    Float cram_size = size(input_bam, "GB")
    Float dynamic_validate_cram_disk_size = cram_size + size(references.ref_fasta, "GB") + 20
    Float validate_cram_disk_size = if dynamic_validate_cram_disk_size > secure_disk_size_threshold then dynamic_validate_cram_disk_size else secure_disk_size_threshold
    Int disk_size = ceil(validate_cram_disk_size)

    command {
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        java -Xms6000m -jar ~{gitc_path}picard.jar \
        ValidateSamFile \
        INPUT=~{input_bam} \
        OUTPUT=~{report_filename} \
        REFERENCE_SEQUENCE=~{references.ref_fasta} \
        ~{"MAX_OUTPUT=" + max_output} \
        IGNORE=~{default="null" sep=" IGNORE=" ignore} \
        MODE=VERBOSE \
        SKIP_MATE_VALIDATION=true \
        IS_BISULFITE_SEQUENCED=is_methyl_seq
    }

    runtime {
        preemptible: preemptible_tries
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File report = "~{report_filename}"
        File monitoring_log = "monitoring.log"
    }
}

# Run STAR alignment. Since it is fast, we can run on all chunks together
task StarAlign {
  input {
    Array[File] input_bams
    String base_file_name
    File genome

    File? gtf_override
    String? extra_args

    Int cpu
    File monitoring_script
    Int disk_size = ceil(3*size(input_bams, "GB") + 3*size(genome, "GB") + 20 )
    Int preemptible_tries
    Boolean no_address
    String docker
    }
    String genome_dir = "genome_dir"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir ~{genome_dir}
        unzip ~{genome} -d ~{genome_dir}

        STAR \
            --readFilesIn ~{sep=',' input_bams} \
            --readFilesType SAM SE \
            --readFilesCommand samtools view -h -F0x200\
            --genomeDir ~{genome_dir} \
            --runThreadN ~{cpu} \
            --quantMode GeneCounts \
             ~{"--sjdbGTFfile " + gtf_override} \
            ~{default="" extra_args} \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix ~{base_file_name}.
    >>>

    runtime {
        preemptible: "~{preemptible_tries}"
        cpu: "~{cpu}"
        memory: "64 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File star_log_file = "~{base_file_name}.Log.final.out"
        File star_log_params_file = "~{base_file_name}.Log.out"
        File output_bam = "~{base_file_name}.Aligned.out.bam"
        File reads_per_gene_file = "~{base_file_name}.ReadsPerGene.out.tab"
        File sj_file = "~{base_file_name}.SJ.out.tab"
        File monitoring_log = "monitoring.log"
    }
}

# Run STAR in genome generate mode to generate the genome dir structure needed for alignment
task StarGenomeGenerate {
    input {
        Array[File] fasta_files
        File gtf_file
        String output_basename
        String? extra_args
        
        Int cpu
        Int disk_size = ceil(3*(size(fasta_files, "GB")) + size(gtf_file, "GB") + 100 )
        Int preemptible_tries
        Boolean no_address
        String docker
        File monitoring_script
    }

    String genome_dir = "genome_dir"
    String output_genome_zip_filename = "~{output_basename}.zip"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        mkdir ~{genome_dir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN ~{cpu} \
            --genomeDir ~{genome_dir} \
            --genomeFastaFiles ~{sep=" " fasta_files} \
            --sjdbGTFfile ~{gtf_file} \
            ~{default="" extra_args}

        # zip the generated genome to be used in future tasks (zip without comprassion for faster zip/unzip opretions)
        zip -0 -r -j ~{output_genome_zip_filename} ~{genome_dir}/*
    >>>

    runtime {
        preemptible: "~{preemptible_tries}"
        cpu: "~{cpu}"
        memory: "40 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File star_log_file = "~{genome_dir}/Log.out"
        File genome_zip = "~{output_genome_zip_filename}"
        File monitoring_log = "monitoring.log"
    }
}


# Parse STAR log.finle.out file to get stats in csv format
task StarAlignStats {
        input {
        File star_log_file
        String base_file_name

        Int disk_size = ceil(3*size(star_log_file, "GB") + 20 )
        Int preemptible_tries
        Boolean no_address
        String docker
        File monitoring_script
    }
    String stats_filename = "~{base_file_name}.star_stats.csv"

    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        conda activate genomics.py3

        python <<CODE
        import pandas as pd
        
        df = pd.read_csv("~{star_log_file}", header = None, sep = '\t')
        df.columns = ['metric', 'value']

        # parse metric description
        df['metric'] = df['metric'].str.replace('|', '').str.strip().str.replace(' ', '_')
        df.loc[:, "metric"] = df["metric"].str.replace(",", "").str.replace(":", "")

        # parse metric type: add metric type (general, unique_reads, multi_mapping_reads, unmapped_reads, chimeric_reads)
        df.loc[:, "metric_type"] = df["metric"].where(df["value"].isnull()).fillna(method="ffill").fillna("general").str.lower().str.replace(":", "")
        df = df.dropna(subset=["value"])

        # parse value: add value type (number, percentage, datetime)
        df['value_type'] = 'number'
        df.loc[df['value'].str.endswith('%'),'value_type'] = 'percentage'
        df.loc[df['value'].str.find(':') != -1,'value_type'] = 'datetime'
        df['value'] = df['value'].str.replace('%', '')

        df.to_csv("~{stats_filename}", index=False)
        
        CODE
    >>>
    runtime {
        preemptible: "~{preemptible_tries}"
        cpu: "1"
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
    }

    output {
        File star_stats = "~{stats_filename}"
        File monitoring_log = "monitoring.log"
    }
}

task SortBam {
    input {
    File input_bam 
    Int disk_size = ceil(3*size(input_bam, "GB"))
    String gitc_path
    String docker
    File monitoring_script
    }
    String base_file_name = basename(input_bam, ".bam")
    command <<<
        set -eo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        java -Xmx64G -jar ~{gitc_path}picard.jar SortSam \
            I=~{input_bam} \
            O=~{base_file_name}.sorted.bam \
            SORT_ORDER=queryname
    >>>
    runtime {
        cpu: "1"
        memory: "8 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
    }

    output {
        File sorted_bam = "~{base_file_name}.sorted.bam"
        File monitoring_log = "monitoring.log"
    }
}