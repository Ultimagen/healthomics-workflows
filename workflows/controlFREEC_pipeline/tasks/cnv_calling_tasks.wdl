version 1.0
import "structs.wdl"

# TASK DEFINITIONS
task CnmopsGetReadCountsFromBam{
    input {
        File input_bam_file
        File input_bai_file
        File reference_genome
        File reference_genome_index
        Int mapq
        Array[String] ref_seq_names
        Int window_length
        String base_file_name
        Boolean save_hdf
        String docker
        Boolean no_address
        File monitoring_script
        Int preemptible_tries
    }

    Float input_bam_file_size = size(input_bam_file, "GB")
    Float additional_disk = 100
    Int disk_size = ceil(input_bam_file_size + input_bam_file_size + additional_disk)
    Int cpu = if size(input_bam_file, "GB")<8 then 2 else ceil(size(input_bam_file, "GB"))/8


    String out_bam_filtered='~{base_file_name}.MAPQ~{mapq}.bam'
    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        samtools view ~{input_bam_file} -O BAM -o ~{out_bam_filtered} -bq ~{mapq} -T ~{reference_genome}
        samtools index ~{out_bam_filtered}

        Rscript --vanilla  /home/ugbio/src/cnv/cnmops/get_reads_count_from_bam.R \
            -i ~{out_bam_filtered} \
            -refseq ~{sep="," ref_seq_names} \
            -wl ~{window_length} \
            -p ~{cpu} \
            -o ~{base_file_name} \
            ~{true="--save_hdf" false='' save_hdf}

        touch ~{base_file_name}.ReadCounts.hdf5
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "16 GiB"
        disks: "local-disk " + ceil(disk_size) + " LOCAL"
        docker: docker
        noAddress: no_address
        cpu: cpu
    }
    output {
        File out_reads_count="~{base_file_name}.ReadCounts.rds"
        File out_reads_count_hdf5="~{base_file_name}.ReadCounts.hdf5"
        String out_sample_name = "~{out_bam_filtered}"
        File monitoring_log = "monitoring.log"

    }
}

task ConvertBedGraphToGranges {
    input{
        String sample_name
        Array[File] input_bed_graph
        File genome_windows
        File genome_file
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }
    Int disk_size = ceil(5*size(input_bed_graph,"GB") + size(genome_windows,"GB") + 40)
    String input_bedGraph_basename = basename(input_bed_graph[0], ".bedGraph")
    Boolean multiple_bedgraph_files = length(input_bed_graph) > 1

    command <<<
        bash ~{monitoring_script} > monitoring.log &
        set -eo pipefail

        for bedgraph in ~{sep=" " input_bed_graph}; 
        do 
            file_basename=$(basename $bedgraph)
            if [[ $bedgraph =~ \.gz$ ]]; 
            then 
                gzip -d -c $bedgraph > "$file_basename.bedgraph"; 
            else
                cp $bedgraph "$file_basename.bedgraph";
            fi
            bedtools map -g ~{genome_file} -a ~{genome_windows} -b "$file_basename.bedgraph" -c 4 -o mean | \
            awk '{if($4=="."){print $1"\t"$2"\t"$3"\t"0}else{print $1"\t"$2"\t"$3"\t"$4}}' \
            > "$file_basename.bedgraph.mean";
        done
        
        if ~{multiple_bedgraph_files}
        then    
            bedtools unionbedg -i ./*.bedgraph.mean | \
            awk -v OFS="\t" -F "\t" 'NR>0{sum=0; for(i=4; i<=NF; i++) sum += $i; NF++; $NF=sum } 1' | \
            awk '{print $1"\t"$2"\t"$3"\t"$NF}' > ~{input_bedGraph_basename}.win.bedGraph
        else
            cp "$file_basename.bedgraph.mean" ~{input_bedGraph_basename}.win.bedGraph
        fi            

        Rscript --vanilla /home/ugbio/src/cnv/cnmops/convert_bedGraph_to_Granges.R \
        -i ~{input_bedGraph_basename}.win.bedGraph \
        -sample_name ~{sample_name}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "4 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File merged_bedGraph="~{input_bedGraph_basename}.win.bedGraph"
        File out_RC_Granges="~{sample_name}.ReadCounts.rds"
        File monitoring_log = "monitoring.log"
    }
}
task AddCountsToCohortMatrix {
    input {
        File sample_reads_count
        File cohort_reads_count_matrix
        Boolean save_hdf
        String docker

        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float cohort_reads_count_matrix_size = 3*size(cohort_reads_count_matrix, "GB")
    Float sample_reads_count_file_size = size(sample_reads_count, "GB")
    Float additional_disk = 50
    Int disk_size = ceil(cohort_reads_count_matrix_size + sample_reads_count_file_size + additional_disk)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        Rscript --vanilla /home/ugbio/src/cnv/cnmops/merge_reads_count_sample_to_cohort.R \
            -cohort_rc ~{cohort_reads_count_matrix} \
            -sample_rc ~{sample_reads_count} \
            ~{true="--save_hdf" false='' save_hdf}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "32 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        File merged_cohort_reads_count_matrix="merged_cohort_reads_count.rds"
        File monitoring_log = "monitoring.log"
    }
}
task CreateCohortReadsCountMatrix {
    input {
        Array[File] sample_reads_count_files
        Boolean save_csv
        Boolean save_hdf
        String docker

        Int parallel
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float sample_reads_count_file_size = 2*size(sample_reads_count_files, "GB")
    Float additional_disk = 10
    Int disk_size = ceil(sample_reads_count_file_size + additional_disk)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail
        echo "~{sep="\n" sample_reads_count_files}" >> samples_list.txt
        Rscript --vanilla /home/ugbio/src/cnv/cnmops/create_reads_count_cohort_matrix.R \
            -samples_read_count_files_list samples_list.txt \
            ~{true="--save_csv" false='' save_csv} \
            ~{true="--save_hdf" false='' save_hdf}
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "64 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: parallel
    }
    output {
        File merged_cohort_reads_count_matrix="merged_cohort_reads_count.rds"
        File monitoring_log = "monitoring.log"
    }
}

task RunCnmops {
    input {
        File merged_cohort_reads_count_matrix
        File? ploidy
        String? chrX_name
        String? chrY_name
        Boolean cap_coverage
        Int min_width_value
        Boolean save_hdf
        Boolean save_csv
        Boolean mod_cnv
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
        Int parallel
    }

    Float merged_cohort_reads_count_matrix_file_size = size(merged_cohort_reads_count_matrix, "GB")
    Float additional_disk = 25
    Int disk_size = ceil(merged_cohort_reads_count_matrix_file_size + additional_disk)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -eo pipefail

        # use touch for optional output files
        touch .estimate_gender
        touch chrX_mean_coverage_distribution.png

        Rscript --vanilla /home/ugbio/src/cnv/cnmops/normalize_reads_count.R \
             --cohort_reads_count_file ~{merged_cohort_reads_count_matrix} \
             ~{"--ploidy " + ploidy} \
             ~{"--chrX_name " + chrX_name} \
             ~{"--chrY_name " + chrY_name} \
             ~{true="--cap_coverage" false="" cap_coverage}

        Rscript --vanilla /home/ugbio/src/cnv/cnmops/cnv_calling_using_cnmops.R \
            -cohort_rc cohort_reads_count.norm.rds \
            -minWidth ~{min_width_value} \
            -p ~{parallel} \
            ~{true="--save_hdf" false='' save_hdf} \
            ~{true="--save_csv" false='' save_csv} \
            ~{true="--mod_cnv" false='' mod_cnv}

        touch cohort.cnmops_outputs.hdf5
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "32 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: parallel
    }
    output {
        File cohort_reads_count_norm = "cohort_reads_count.norm.rds"
        File cohort_cnvs_int_rds = "cohort.cnmops.resCNMOPS_Int.rds"
        File cohort_cnvs_csv = "cohort.cnmops.cnvs.csv"
        File cohort_cnmops_output_hdf5 = "cohort.cnmops_outputs.hdf5"
        File cohort_estimated_gender = ".estimate_gender"
        File cohort_mean_chrX_dist = "chrX_mean_coverage_distribution.png"
        File monitoring_log = "monitoring.log"
    }
}

task ExtractNormalizedReadCount {
    input {
        File cohort_reads_count_norm
        Array[String] sample_names
        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Float cohort_reads_count_norm_file_size = size(cohort_reads_count_norm, "GB")
    Float additional_disk = 10
    Int disk_size = ceil(cohort_reads_count_norm_file_size + additional_disk)

    command <<<
        bash ~{monitoring_script} | tee monitoring.log >&2 &
        set -exo pipefail
        for sample_name in ~{sep=" " sample_names}; do
            Rscript --vanilla /home/ugbio/src/cnv/cnmops/export_cohort_matrix_to_bed.R  \
                ~{cohort_reads_count_norm} --sample "$sample_name"
        done
        Rscript --vanilla  /home/ugbio/src/cnv/cnmops/export_cohort_matrix_to_bed.R  \
            ~{cohort_reads_count_norm} --mean
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 2
    }
    output {
        Array[File] sample_reads_count_bed = glob("*.cov.bed")
        File cohort_reads_count_bed = "coverage.cohort.bed"
        File monitoring_log = "monitoring.log"
    }
}

task ProcessCnmopsCnvs {
    input {
        File cohort_cnvs_csv
        Array[String] sample_names
        Int min_cnv_length
        Array[File] sample_norm_coverage_file
        File cohort_norm_avg_coverage_file
        Float intersection_cutoff
        File? cnv_lcr_file
        File ref_genome_file
        File germline_coverage_rds
        String docker
        Boolean skip_figure_generation
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    Boolean filter_lcr = defined(cnv_lcr_file)
    Float cohort_cnvs_csv_file_size = size(cohort_cnvs_csv, "GB")
    Float ref_genome_file_size = size(ref_genome_file, "GB")
    Float additional_disk = 25
    Int disk_size = ceil(cohort_cnvs_csv_file_size + ref_genome_file_size + additional_disk)
    
    command <<<
        set -xeo pipefail

        bash ~{monitoring_script} | tee monitoring.log >&2 &

        for f in ~{sep=' ' sample_norm_coverage_file}; do
            echo $f >> filelist.txt
        done

        #write all samples coverage to bed files
        Rscript --vanilla  /home/ugbio/src/cnv/cnmops/export_cohort_matrix_to_bed.R  ~{germline_coverage_rds}
        #get samples CNVs

        for sample_name in ~{sep=" " sample_names}
        do 
            if grep -q "$sample_name" ~{cohort_cnvs_csv};
                then grep "$sample_name" ~{cohort_cnvs_csv} > $sample_name.cnvs.csv ;
                awk -F "," '{print $1"\t"$2-1"\t"$3"\t"$NF}' $sample_name.cnvs.csv > $sample_name.cnvs.bed;

                ## looking for the depth file that matches the sample
                while read -r line; do
                    bn=$(basename "$line" .cov.bed)
                    if [[ "$bn" == "$sample_name" ]]; then
                        desired_file="$line"
                        break
                    fi
                done < filelist.txt
                
                echo "$sample_name coverage file: $desired_file" >&2
                if ~{filter_lcr}; then 
                    process_cnvs \
                        --sample_name "$sample_name" \
                        --input_bed_file $sample_name.cnvs.bed \
                        --intersection_cutoff ~{intersection_cutoff} \
                        --cnv_lcr_file ~{cnv_lcr_file} \
                        --min_cnv_length ~{min_cnv_length} \
                        --out_directory . \
                        --sample_norm_coverage_file "$desired_file" \
                        --cohort_avg_coverage_file ~{cohort_norm_avg_coverage_file} \
                        --fasta_index_file ~{ref_genome_file}

                                        
                else
                    #convert bed file to vcf
                    process_cnvs \
                        --sample_name "$sample_name" \
                        --input_bed_file $sample_name.cnvs.bed \
                        --intersection_cutoff ~{intersection_cutoff} \
                        --min_cnv_length ~{min_cnv_length} \
                        --out_directory . \
                        --sample_norm_coverage_file "$desired_file" \
                        --cohort_avg_coverage_file ~{cohort_norm_avg_coverage_file} \
                        --fasta_index_file ~{ref_genome_file}
                fi

                bcftools view -i "INFO/SVTYPE='DEL'" $sample_name.cnvs.annotate.vcf.gz | \
                bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO;FILTER=%FILTER\n' - > "$sample_name".DEL.bed
                
                bcftools view -i "INFO/SVTYPE='DUP'" $sample_name.cnvs.annotate.vcf.gz | \
                bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO;FILTER=%FILTER\n' - > "$sample_name".DUP.bed
                
                #generate figure for each sample
                if ~{skip_figure_generation}; then
                    echo "skip figure generation"; 
                else
                    cut -f1-3,5 $sample_name.cov.bed > $sample_name.cov.for_plot.bed
                    plot_cnv_results \
                        --germline_coverage $sample_name.cov.for_plot.bed \
                        --duplication_cnv_calls $sample_name.DUP.bed \
                        --deletion_cnv_calls $sample_name.DEL.bed \
                        --sample_name $sample_name \
                        --out_directory CNV_figures \
                        --vcf-like
                fi

            else 
                echo "$sample_name not found in ~{cohort_cnvs_csv}"
            fi
        done

    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "8 GiB"
        disks: "local-disk " + ceil(disk_size) + " HDD"
        docker: docker
        noAddress: no_address
        cpu: 4
    }
    output {
        Array[File] csvs = glob("*.csv")
        Array[File] sample_cnvs_vcf = glob("*.cnvs.annotate.vcf.gz")
        Array[File] sample_cnvs_vcf_index = glob("*.cnvs.annotate.vcf.gz.tbi")
        Array[File] coverage_plot = glob("*/*.CNV.coverage.jpeg")
        Array[File] dup_del_plot = glob("*/*.dup_del.calls.jpeg")
        Array[File] copy_number_plot = glob("*/*.CNV.calls.jpeg")
        File monitoring_log = "monitoring.log"
    }
}


task UniqReads {
    input {
        String sampleId
        String tempPrefix = "~{sampleId}_"
        Array[String] tempSeqsPaths
        File finalBam
        File finalBam_index
        Int memoryGb
        String docker
    }
    Int diskSize = ceil( size(finalBam, "GB") ) + 100
    command {
        #adding first command to support cram as input
        samtools view -b -h ~{finalBam} | \
        /samtools-0.1.7a_getUnique-0.1.3/samtools \
        view \
        -U "BWA,~{tempPrefix},N,N" \
        -
    }

    output {
        Array[File] tempSeqs = tempSeqsPaths
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GiB"
        docker : docker
        runtime_minutes: "500"
    }
}

task Bicseq2Norm {
    input {
        Int memoryGb
        Int readLength
        Int medianInsertSize
        String GCvsRDPath = "~{sampleId}.GCvsRD.pdf"
        String paramsPath = "~{sampleId}.params.out"
        String sampleId

        Array[File] tempSeqs
        Array[String] tempNormPaths

        File bicseq2ConfigFile
        String configFilePath = "~{sampleId}.bicseq2.config"

        Array[File] uniqCoordsFiles
        Array[File] chromFastasFiles
        Int diskSize = 70
        String docker
    }

    command {
        set -e -o pipefail

        python3 \
        /bicseq2_config_writer.py \
        --fa-files ~{sep=" " chromFastasFiles} \
        --mappability-files ~{sep=" " uniqCoordsFiles} \
        --temp-seqs ~{sep=" " tempSeqs} \
        --temp-norm-paths ~{sep=" " tempNormPaths} \
        --norm-bicseq2-config ~{bicseq2ConfigFile} \
        --sample-id ~{sampleId} \
        --out-file ~{configFilePath}

        mkdir -p ~{sampleId}

        /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
        -l=~{readLength} \
        -s=~{medianInsertSize} \
        -fig=~{GCvsRDPath} \
        -tmp=~{sampleId} \
        ~{configFilePath} \
        ~{paramsPath}
    }

    output {
        Array[File] tempNorm = tempNormPaths
        File params = paramsPath
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GiB"
        docker : docker
        runtime_minutes: "500"
    }
}


task Bicseq2Wgs {
    input {
        Int memoryGb
        String pairName
        String bicseq2PngPath = "~{pairName}.bicseq2.png"
        String bicseq2Path = "~{pairName}.bicseq2.txt"
        Array[File] tempTumorNorms
        Array[File] tempNormalNorms
        File bicseq2SegConfigFile
        String segConfigFilePath = "~{pairName}.bicseq2.seg.config"
        Int lambda = 4
        Int diskSize = 10
        String docker
    }

    command {
        set -e -o pipefail

        mkdir -p ~{pairName}

        python3 \
        /bicseq2_seg_config_writer.py \
        --tumor-norms ~{sep=" " tempTumorNorms} \
        --normal-norms ~{sep=" " tempNormalNorms} \
        --seg-bicseq2-config ~{bicseq2SegConfigFile} \
        --out-file ~{segConfigFilePath} \
        --pair-id ~{pairName}

        perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
        --control \
        --tmp ~{pairName} \
        --fig=~{bicseq2PngPath} \
        --title=~{pairName} \
        --lambda=4 \
        ~{segConfigFilePath} \
        ~{bicseq2Path}

    }

    output {
        File bicseq2Png = "~{pairName}.bicseq2.png"
        File bicseq2 = "~{pairName}.bicseq2.txt"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GiB"
        docker : docker
        runtime_minutes: "60"
    }
}

task Bicseq2PostProcessing {
    input {
        File input_bicseq2_txt_file
        Float? ratio_DEL_cutoff
        Float? ratio_DUP_cutoff
        String docker
    }
    Int diskSize = ceil( size(input_bicseq2_txt_file, "GB") ) + 5
    String out_basename = basename(input_bicseq2_txt_file, ".txt")
    Int memoryGb = 8

    command {
        set -eo pipefail

        bicseq2_post_processing \
        --input_bicseq2_txt_file ~{input_bicseq2_txt_file} \
        ~{"--ratio_DUP_cutoff " + ratio_DUP_cutoff} \
        ~{"--ratio_DEL_cutoff " + ratio_DEL_cutoff}
    }

    output {
        File out_bed_file = "~{out_basename}.bed"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GiB"
        docker : docker
        runtime_minutes: "60"
    }
}


task CnvVcfToBed {
    input {
        File input_cnv_vcf
        String base_file_name

        String docker
        File monitoring_script
        Boolean no_address
        Int preemptible_tries
    }

    command <<<
        set -xeo pipefail
        bash ~{monitoring_script} | tee monitoring.log >&2 &

        bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%INFO;FILTER=%FILTER\n' ~{input_cnv_vcf} > ~{base_file_name}.cnv.bed
    >>>
    runtime {
        preemptible: preemptible_tries
        memory: "4 GiB"
        disks: "local-disk 10 LOCAL"
        docker: docker
        noAddress: no_address
        cpu: 2
    }

    output {
        File output_cnv_bed = "~{base_file_name}.cnv.bed"
    }
}

