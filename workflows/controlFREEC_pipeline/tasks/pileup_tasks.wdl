version 1.0 
import "structs.wdl"

task BcftoolsMpileupTask {
    input {
        References references
        File input_cram_bam
        File input_cram_bam_index
        File regions_bed
        Int min_mapq
        Boolean large_file=false
        String base_file_name
        File monitoring_script
        String docker
        Float disk_size = ceil(size(input_cram_bam, "GB")*2 + size(references.ref_fasta, "GB") + size(regions_bed, "GB") + 1)
        Float memory_gb = 4
        Int cpus = 2
        Int preemptibles
    }
    command <<<
        set -eo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &
        source ~/.bashrc
        conda activate genomics.py3

        bcftools mpileup ~{input_cram_bam} \
            -f ~{references.ref_fasta} \
            ~{true="-T" false="-R" large_file} ~{regions_bed} \
            -a FORMAT/AD -Oz \
            -o ~{base_file_name}.mpileup.vcf.gz \
            --skip-indels -q ~{min_mapq} \
            --threads ~{cpus} && \
            bcftools index -t ~{base_file_name}.mpileup.vcf.gz

        echo "*********** extract allele coverage ***********"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' ~{base_file_name}.mpileup.vcf.gz | \
            awk -F"\t" 'BEGIN{OFS="\t"}{split($5,ad,","); print $1,$2,$3,$4,ad[1],ad[2],$6}' > ~{base_file_name}.allele_coverage.csv
    >>>
    output {
        File mpileup_vcf = "~{base_file_name}.mpileup.vcf.gz"
        File mpileup_vcf_index = "~{base_file_name}.mpileup.vcf.gz.tbi"
        File allele_coverage_file = "~{base_file_name}.allele_coverage.csv"
    }
    runtime {
      preemptible: "~{preemptibles}"
      cpu: "~{cpus}"
      memory: "~{memory_gb} GB"
      disks: "local-disk " + ceil(disk_size) + " HDD"
      docker: docker
    }
}


task SamtoolsMpileup {
    input {
        File input_cram_bam
        File input_cram_bam_index
        References references
        String base_file_name
        File monitoring_script
        String docker
        Float disk_size = ceil(size(input_cram_bam, "GB")*2 + size(references.ref_fasta, "GB") + 1)
        Int preemptible
    }
    command <<<
        set -eo pipefail
        set -x 
        bash ~{monitoring_script} | tee monitoring.log >&2 &
   
        samtools mpileup -f ~{references.ref_fasta} \
            -Q 0 \
            -o ~{base_file_name}.mpileup \
            ~{input_cram_bam}
    >>>
    output {
        File mpileup = "~{base_file_name}.mpileup"
    }
    runtime {
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        preemptible: "~{preemptible}"
    }
}

task CreateMpileup {
    input{

        Array[File] input_bam_files
        Array[File] input_bam_files_index
        File reference_fasta
        File reference_fai
        File reference_dict
        Int max_depth = 8000
        Int min_MapQ
        Int min_BaseQ = 0
        File snp_file
        File? snp_file_index
        File interval
        String docker
        Boolean no_address
        Int preemptible_tries
        File monitoring_script
        String? cloud_provider
    }
    Int disk_size = ceil(size(reference_fasta,"GB") + size(snp_file,"GB") + 50)
    String base_input_name = basename(input_bam_files[0])
    Boolean is_aws = defined(cloud_provider)
    String snp_file_basename = basename(snp_file)


    parameter_meta {
      input_bam_files: {
          localization_optional: true
      }
    }

command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -xeou pipefail

    # if ~{is_aws}
    # then
    #     inputs_string="~{sep=' ' input_bam_files}"
    # else
    echo "DEBUG start PrintReads $(date)"
    # download only the region of interval
    gatk --java-options "-Xms1G" PrintReads \
        -I ~{sep=' -I ' input_bam_files} \
        -O input.bam \
        -L ~{interval} \
        -R ~{reference_fasta}
    inputs_string="input.bam"
    echo "DEBUG end PrintReads $(date)"
    # fi

    echo "inputs_string:"
    echo $inputs_string
    echo "DEBUG start intersect snp_file with current interval $(date)"
    # intersect snp_file with current interval
    gatk IntervalListToBed -I ~{interval} -O interval.bed
    
    cp ~{snp_file} ./~{snp_file_basename}
    file=./~{snp_file_basename}
    if [[ "$file" == *.vcf.gz ]]; then
        bcftools index -t $file
        pos_file="out.vcf.gz"
        bcftools view -R interval.bed -O z -o $pos_file $file && bcftools index -t $pos_file
    elif [[ "$file" == *.bed.gz ]]; then
        gunzip $file
        bed_file="${file%.gz}"
        pos_file="out.bed"
        bedtools intersect -a $bed_file -b interval.bed -header > $pos_file
    elif [[ "$file" == *.bed ]]; then
        pos_file="out.bed"
        bedtools intersect -a $file -b interval.bed -header > $pos_file
    else
        echo "positions file suffix is none of : [.vcf.gz , .bed.gz , .bed]"
    fi

    echo "DEBUG end intersect snp_file with current interval $(date)"

    echo "DEBUG start mpileup $(date)"
    #caclulate mpileup for current interval
    samtools mpileup -f ~{reference_fasta} \
    -d ~{max_depth} \
    -Q ~{min_BaseQ} \
    -q ~{min_MapQ} \
    -l $pos_file \
    $inputs_string \
    >  ~{base_input_name}_minipileup.pileup
    echo "DEBUG end mpileup $(date)"

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "2 GB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu:1
    }
    output {
        File out_pileup = "~{base_input_name}_minipileup.pileup"
        File monitoring_log = "monitoring.log"
    }
}