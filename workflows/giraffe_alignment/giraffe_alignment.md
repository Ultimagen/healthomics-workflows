# GiraffeAlignment
Single-sample alignment + sort + mark-duplicates using vg giraffe only. This is a simplified variant of AlignCramOrBam with all non-giraffe aligner options removed.

## Inputs

### Required inputs
<p name="GiraffeAlignment.input_cram_list">
        <b>GiraffeAlignment.input_cram_list</b><br />
        <i>Array[File] </i> &mdash;
         List of input CRAM/BAM files to align using vg giraffe. All files in the list are processed together. <br />
</p>
<p name="GiraffeAlignment.base_file_name">
        <b>GiraffeAlignment.base_file_name</b><br />
        <i>String </i> &mdash;
         Base name for all output files. Should not contain spaces, #, or comma characters. <br />
</p>
<p name="GiraffeAlignment.ref_files_for_tarball">
        <b>GiraffeAlignment.ref_files_for_tarball</b><br />
        <i>Array[File] </i> &mdash;
         List of reference FASTA files used to build the reference cache. <br />
</p>
<p name="GiraffeAlignment.sorter_params">
        <b>GiraffeAlignment.sorter_params</b><br />
        <i>SorterParams </i> &mdash;
         Parameters for sorting and marking duplicates. <br />
</p>

### Required references
<p name="GiraffeAlignment.references">
        <b>GiraffeAlignment.references</b><br />
        <i>References </i> &mdash;
         Linear reference genome (FASTA, FAI, DICT) used for read processing and downstream sorting/markdup. <br />
</p>
<p name="GiraffeAlignment.giraffe_parameters">
        <b>GiraffeAlignment.giraffe_parameters</b><br />
        <i>GiraffeReferences </i> &mdash;
         Giraffe graph reference bundle (GBZ/DIST/MIN/ZIPCODES/PATH_LIST) plus linear reference (FASTA/FAI/DICT). <br />
</p>

### Optional parameters
<p name="GiraffeAlignment.preemptible_tries">
        <b>GiraffeAlignment.preemptible_tries</b><br />
        <i>Int? </i> &mdash;
         Number of preemptible retries for scatter tasks. <br />
</p>
<p name="GiraffeAlignment.output_haplotypes_cram">
        <b>GiraffeAlignment.output_haplotypes_cram</b><br />
        <i>Boolean </i> &mdash;
         When true, run haplotype sampling and output a haplotype CRAM. <br />
</p>
<p name="GiraffeAlignment.reads_per_split">
        <b>GiraffeAlignment.reads_per_split</b><br />
        <i>Int </i> &mdash;
         Approximate number of reads per split chunk for parallel processing. <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.kmer_length">
        <b>GiraffeAlignment.HaplotypeSampling.kmer_length</b><br />
        <i>Int </i> &mdash;
         K-mer length for KMC counting (default: 29) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.num_haplotypes">
        <b>GiraffeAlignment.HaplotypeSampling.num_haplotypes</b><br />
        <i>Int </i> &mdash;
         Number of haplotypes to sample (default: 32) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.include_reference">
        <b>GiraffeAlignment.HaplotypeSampling.include_reference</b><br />
        <i>Boolean </i> &mdash;
         Include reference in sampled haplotypes (default: true) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.diploid_sampling">
        <b>GiraffeAlignment.HaplotypeSampling.diploid_sampling</b><br />
        <i>Boolean </i> &mdash;
         Use diploid sampling strategy (default: true) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.window_size">
        <b>GiraffeAlignment.HaplotypeSampling.window_size</b><br />
        <i>Int </i> &mdash;
         Sliding window size for seqkit (default: 50000) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.step_size">
        <b>GiraffeAlignment.HaplotypeSampling.step_size</b><br />
        <i>Int </i> &mdash;
         Sliding window step size for seqkit (default: 50000) <br />
</p>
<p name="GiraffeAlignment.HaplotypeSampling.minimap2_preset">
        <b>GiraffeAlignment.HaplotypeSampling.minimap2_preset</b><br />
        <i>String </i> &mdash;
         Minimap2 preset for alignment (default: asm5) <br />
</p>
</details>


## Outputs
<p name="GiraffeAlignment.output_cram">
        <b>GiraffeAlignment.output_cram</b><br />
        <i>File?</i><br />
        Giraffe-aligned cram file.
</p>
<p name="GiraffeAlignment.output_cram_index">
        <b>GiraffeAlignment.output_cram_index</b><br />
        <i>File?</i><br />
        Index file for the Giraffe-aligned cram file.
</p>
<p name="GiraffeAlignment.sorter_stats_csv">
        <b>GiraffeAlignment.sorter_stats_csv</b><br />
        <i>File?</i><br />
        Sorter stats in csv format.
</p>
<p name="GiraffeAlignment.sorter_stats_json">
        <b>GiraffeAlignment.sorter_stats_json</b><br />
        <i>File?</i><br />
        Sorter stats in json format.
</p>
<p name="GiraffeAlignment.output_cram_bam_list">
        <b>GiraffeAlignment.output_cram_bam_list</b><br />
        <i>Array[File?]?</i><br />
        Output file list after the pipeline is executed in multiple outputs mode
</p>
<p name="GiraffeAlignment.output_cram_bam_index_list">
        <b>GiraffeAlignment.output_cram_bam_index_list</b><br />
        <i>Array[File?]?</i><br />
        Index files for the output cram files (when running in multiple outputs mode).
</p>
<p name="GiraffeAlignment.sorter_stats_csv_list">
        <b>GiraffeAlignment.sorter_stats_csv_list</b><br />
        <i>Array[File?]?</i><br />
        Sorter stats in csv format (when running in multiple outputs mode).
</p>
<p name="GiraffeAlignment.sorter_stats_json_list">
        <b>GiraffeAlignment.sorter_stats_json_list</b><br />
        <i>Array[File?]?</i><br />
        Sorter stats in json format (when running in multiple outputs mode).
</p>
<p name="GiraffeAlignment.haplotypes_cram">
        <b>GiraffeAlignment.haplotypes_cram</b><br />
        <i>File?</i><br />
        Haplotype-sampled CRAM output (when enabled).
</p>
<p name="GiraffeAlignment.haplotypes_cram_index">
        <b>GiraffeAlignment.haplotypes_cram_index</b><br />
        <i>File?</i><br />
        Index for haplotype-sampled CRAM output (when enabled).
</p>

<hr />

> Generated using WDL AID (1.0.1)