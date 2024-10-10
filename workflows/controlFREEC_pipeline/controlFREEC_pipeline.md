# SomaticCNVCallingControlFREEC
Runs single sample somatic CNV calling workflow based on \<a href="https://boevalab.inf.ethz.ch/FREEC/"\>ControlFREEC</a>.

CNVs are called based on both coverage and allele frequencies in the tumor and the matched germline sample.

The pipeline will gather coverage and allele frequencies, run controlFREEC and filter called CNVs by length and low confidence regions.

coverage will be calculted based on the input cram/bam. Alternativley, it can recieve coverage input as one of:bedGraph, cpn formats.

Allele frequencies will be calculated based on the input cram/bam and a given vcf file to specify locations. Alternativley, it can recieve precalculated frequencies as mpileup format.

The pipeline outputs: 

&nbsp;&nbsp;- calculated coverage for tumor and normal samples

&nbsp;&nbsp;- calculated mpileup for tumor and normal samples

&nbsp;&nbsp;- called CNVs + filtered called CNVs

&nbsp;&nbsp;- controlFREEC run-summary

&nbsp;&nbsp;-coverage plot that shows normalized (log scale) coverage along the genome for the germline and tumor samples.

&nbsp;&nbsp;-duplications and deletions figure - showing gains and losses along the genome.

&nbsp;&nbsp;-copy-number figure  shows the copy number along the genome.



## Inputs

### Required inputs
<p name="SomaticCNVCallingControlFREEC.base_file_name">
        <b>SomaticCNVCallingControlFREEC.base_file_name</b><br />
        <i>String </i> &mdash; 
         Base file name used for some output files <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.input_tumor_cram_bam_file">
        <b>SomaticCNVCallingControlFREEC.input_tumor_cram_bam_file</b><br />
        <i>Array[File]+ </i> &mdash; 
         Input tumor BAM/CRAM files <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.input_tumor_cram_bam_file_index">
        <b>SomaticCNVCallingControlFREEC.input_tumor_cram_bam_file_index</b><br />
        <i>Array[File]+ </i> &mdash; 
         Input tumor BAI/CRAI index files <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.input_normal_cram_bam_file">
        <b>SomaticCNVCallingControlFREEC.input_normal_cram_bam_file</b><br />
        <i>Array[File]+ </i> &mdash; 
         Input normal BAM/CRAM files <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.input_normal_cram_bam_file_index">
        <b>SomaticCNVCallingControlFREEC.input_normal_cram_bam_file_index</b><br />
        <i>Array[File]+ </i> &mdash; 
         Input normal BAI/CRAI index files <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.SNPfile">
        <b>SomaticCNVCallingControlFREEC.SNPfile</b><br />
        <i>File </i> &mdash; 
         Vcf file holding locations of the common variants to calculate pileup statistics on <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.SNPfile_index">
        <b>SomaticCNVCallingControlFREEC.SNPfile_index</b><br />
        <i>File </i> &mdash; 
         Vcf.tbi index file for SNPfile <br /> 
</p>

### Required parameters
<p name="SomaticCNVCallingControlFREEC.interval_list">
        <b>SomaticCNVCallingControlFREEC.interval_list</b><br />
        <i>File </i> &mdash; 
         Interval list defining the regions to gather allele frequencies on, recommended value set in the template <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.genome_windows">
        <b>SomaticCNVCallingControlFREEC.genome_windows</b><br />
        <i>File </i> &mdash; 
         Bed file of the genome binned to equal sized windows <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.window">
        <b>SomaticCNVCallingControlFREEC.window</b><br />
        <i>Int </i> &mdash; 
         The size of the window over which the coverage is aggregated <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.min_cnv_length">
        <b>SomaticCNVCallingControlFREEC.min_cnv_length</b><br />
        <i>Int </i> &mdash; 
         Minimum length for reported CNVs. Default is 10,000 <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.intersection_cutoff">
        <b>SomaticCNVCallingControlFREEC.intersection_cutoff</b><br />
        <i>Float </i> &mdash; 
         Intersection cutoff with UG-CNV-LCR regions to filter out CNV calls. Default is  0.5 <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.cnv_lcr_file">
        <b>SomaticCNVCallingControlFREEC.cnv_lcr_file</b><br />
        <i>File </i> &mdash; 
         UG-CNV-LCR bed file <br /> 
</p>

### Required references
<p name="SomaticCNVCallingControlFREEC.references">
        <b>SomaticCNVCallingControlFREEC.references</b><br />
        <i>References </i> &mdash; 
         Struct of reference objects holding reference fasta with corresponding fai and dict files <br /> 
</p>

### Optional inputs
<p name="SomaticCNVCallingControlFREEC.normal_mpileup_override">
        <b>SomaticCNVCallingControlFREEC.normal_mpileup_override</b><br />
        <i>File? </i> &mdash; 
         Pre-calculated mpileup for normal sample <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_mpileup_override">
        <b>SomaticCNVCallingControlFREEC.tumor_mpileup_override</b><br />
        <i>File? </i> &mdash; 
         Pre-calculated mpileup for tumor sample <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.normal_coverage_cpn">
        <b>SomaticCNVCallingControlFREEC.normal_coverage_cpn</b><br />
        <i>File? </i> &mdash; 
         Pre-calculated binned coverage for the normal sample in the format needed for controlFREEC (cpn) <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_coverage_cpn">
        <b>SomaticCNVCallingControlFREEC.tumor_coverage_cpn</b><br />
        <i>File? </i> &mdash; 
         Pre-calculated binned coverage for the tumor sample in the format needed for controlFREEC (cpn) <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.normal_sorter_zipped_bed_graph">
        <b>SomaticCNVCallingControlFREEC.normal_sorter_zipped_bed_graph</b><br />
        <i>Array[File]? </i> &mdash; 
         Pre-calculated bedGraph files containing per-base coverage for the normal sample <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_sorter_zipped_bed_graph">
        <b>SomaticCNVCallingControlFREEC.tumor_sorter_zipped_bed_graph</b><br />
        <i>Array[File]? </i> &mdash; 
         Pre-calculated bedGraph files containing per-base coverage for the tumor sample <br /> 
</p>

### Optional parameters
<p name="SomaticCNVCallingControlFREEC.mapq_override">
        <b>SomaticCNVCallingControlFREEC.mapq_override</b><br />
        <i>Int? </i> &mdash; 
         Reads mapping-quality cutoff for coverage calculation. Default is 1 <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.collect_coverage_region">
        <b>SomaticCNVCallingControlFREEC.collect_coverage_region</b><br />
        <i>Array[String]? </i> &mdash; 
         Genomic region to limit the CNV calling to <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.contaminationFraction">
        <b>SomaticCNVCallingControlFREEC.contaminationFraction</b><br />
        <i>Float? </i> &mdash; 
         a priori known value of tumor sample contamination by normal cells. Default: contamination=0 <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.chrLenFile_override">
        <b>SomaticCNVCallingControlFREEC.chrLenFile_override</b><br />
        <i>File? </i> &mdash; 
         Chromosome lengths file focusing controlFREEC regions. file is expected to be tab-delimited where the first column indicates the chromosome name and the second column indicates the chromosome length. By default, the reference.fai file will be used. <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.ploidy">
        <b>SomaticCNVCallingControlFREEC.ploidy</b><br />
        <i>Int? </i> &mdash; 
         Average sample ploidy (if known) <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.maxThreads_override">
        <b>SomaticCNVCallingControlFREEC.maxThreads_override</b><br />
        <i>Int? </i> &mdash; 
         maximal threads for controlFREEC. Default is 8 <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.sex">
        <b>SomaticCNVCallingControlFREEC.sex</b><br />
        <i>String? </i> &mdash; 
         Sample's sex value, should be 'XX' or 'XY' <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.gemMappabilityFile">
        <b>SomaticCNVCallingControlFREEC.gemMappabilityFile</b><br />
        <i>File? </i> &mdash; 
         Gem file holding mappablity biased regions.  <br /> 
</p>
<p name="SomaticCNVCallingControlFREEC.preemptible_tries_override">
        <b>SomaticCNVCallingControlFREEC.preemptible_tries_override</b><br />
        <i>Int? </i> &mdash; 
         Number of preemptible tries <br /> 
</p>
</details>


## Outputs
<p name="SomaticCNVCallingControlFREEC.tumor_mpileup">
        <b>SomaticCNVCallingControlFREEC.tumor_mpileup</b><br />
        <i>File</i><br />
        mpileup file for tumor sample
</p>
<p name="SomaticCNVCallingControlFREEC.normal_mpileup">
        <b>SomaticCNVCallingControlFREEC.normal_mpileup</b><br />
        <i>File</i><br />
        mpileup file for normal sample
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_coverage">
        <b>SomaticCNVCallingControlFREEC.tumor_coverage</b><br />
        <i>File</i><br />
        Coverage file for tumor sample
</p>
<p name="SomaticCNVCallingControlFREEC.normal_coverage">
        <b>SomaticCNVCallingControlFREEC.normal_coverage</b><br />
        <i>File</i><br />
        Coverage file for normal sample
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_CNVs">
        <b>SomaticCNVCallingControlFREEC.tumor_CNVs</b><br />
        <i>File</i><br />
        controlFREEC predicted copy number alterations for tumor sample
</p>
<p name="SomaticCNVCallingControlFREEC.controlFREEC_info">
        <b>SomaticCNVCallingControlFREEC.controlFREEC_info</b><br />
        <i>File</i><br />
        controlFREEC run summary
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_ratio_bedgraph">
        <b>SomaticCNVCallingControlFREEC.tumor_ratio_bedgraph</b><br />
        <i>File</i><br />
        controlFREEC ratios in BedGraph format
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_CNVs_annotated_bed_file">
        <b>SomaticCNVCallingControlFREEC.tumor_CNVs_annotated_bed_file</b><br />
        <i>File</i><br />
        Called CNVs for tumor sample with LCR and LENGTH annotations
</p>
<p name="SomaticCNVCallingControlFREEC.tumor_CNVs_filtered_bed_file">
        <b>SomaticCNVCallingControlFREEC.tumor_CNVs_filtered_bed_file</b><br />
        <i>File</i><br />
        Filtered called CNVs for tumor sample
</p>
<p name="SomaticCNVCallingControlFREEC.coverage_plot">
        <b>SomaticCNVCallingControlFREEC.coverage_plot</b><br />
        <i>File</i><br />
        Coverage plot that shows normalized (log scale) coverage along the genome for the germline and tumor samples
</p>
<p name="SomaticCNVCallingControlFREEC.dup_del_plot">
        <b>SomaticCNVCallingControlFREEC.dup_del_plot</b><br />
        <i>File</i><br />
        Duplications and deletions figure - showing gains and losses along the genome
</p>
<p name="SomaticCNVCallingControlFREEC.copy_number_plot">
        <b>SomaticCNVCallingControlFREEC.copy_number_plot</b><br />
        <i>File</i><br />
        Copy-number figure  shows the copy number along the genome
</p>

<hr />

> Generated using WDL AID (1.0.0)
