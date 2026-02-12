# STRGenotyper
Alignment-based STR genotype caller using Smith-Waterman alignment

## Inputs

### Required inputs
<p name="STRGenotyper.base_file_name">
        <b>STRGenotyper.base_file_name</b><br />
        <i>String </i> &mdash;
         Prefix for name of all output files <br />
</p>
<p name="STRGenotyper.cram_file">
        <b>STRGenotyper.cram_file</b><br />
        <i>File </i> &mdash;
         Input CRAM file for STR genotyping <br />
</p>
<p name="STRGenotyper.cram_index">
        <b>STRGenotyper.cram_index</b><br />
        <i>File </i> &mdash;
         CRAM index file (.crai) <br />
</p>
<p name="STRGenotyper.variant_catalog">
        <b>STRGenotyper.variant_catalog</b><br />
        <i>File </i> &mdash;
         JSON file containing STR variant catalog with locus definitions <br />
</p>
<p name="STRGenotyper.references">
        <b>STRGenotyper.references</b><br />
        <i>References </i> &mdash;
         Reference genome files (fasta, fasta.fai, dict) as References struct <br />
</p>

### Optional parameters
<p name="STRGenotyper.ref_padding">
        <b>STRGenotyper.ref_padding</b><br />
        <i>Int </i> &mdash;
         Number of bases to extend around the STR repeat region when building auxiliary references for alignment. Larger values provide more flanking sequence context for accurate alignment. <br />
</p>
<p name="STRGenotyper.min_repeat">
        <b>STRGenotyper.min_repeat</b><br />
        <i>Int </i> &mdash;
         Minimum number of repeat units to include in auxiliary reference sequences <br />
</p>
<p name="STRGenotyper.max_repeat">
        <b>STRGenotyper.max_repeat</b><br />
        <i>Int </i> &mdash;
         Maximum number of repeat units to include in auxiliary reference sequences <br />
</p>
<p name="STRGenotyper.min_score_ratio">
        <b>STRGenotyper.min_score_ratio</b><br />
        <i>Float </i> &mdash;
         Minimum ratio of alignment score to the theoretical maximum score (read_length * match_score). Alignments below this threshold are filtered out. Range: 0.0-1.0, where 1.0 requires perfect alignment. <br />
</p>
<p name="STRGenotyper.spanning_flank_bases">
        <b>STRGenotyper.spanning_flank_bases</b><br />
        <i>Int </i> &mdash;
         Minimum number of bases that must align on each side of the STR repeat region for a read to be considered 'spanning' the locus <br />
</p>
<p name="STRGenotyper.min_mapping_quality">
        <b>STRGenotyper.min_mapping_quality</b><br />
        <i>Int </i> &mdash;
         Minimum mapping quality for reads to be included in analysis <br />
</p>
<p name="STRGenotyper.threads">
        <b>STRGenotyper.threads</b><br />
        <i>Int </i> &mdash;
         Number of threads for parallel processing <br />
</p>
<p name="STRGenotyper.output_detailed_csv">
        <b>STRGenotyper.output_detailed_csv</b><br />
        <i>Boolean </i> &mdash;
         Whether to output detailed per-read CSV file. Set to false for large catalogs to reduce I/O. <br />
</p>
<p name="STRGenotyper.output_summary_csv">
        <b>STRGenotyper.output_summary_csv</b><br />
        <i>Boolean </i> &mdash;
         Whether to output summary per-locus CSV file. Set to false for large catalogs to reduce I/O. <br />
</p>
<p name="STRGenotyper.haploid">
        <b>STRGenotyper.haploid</b><br />
        <i>Boolean </i> &mdash;
         Enable haploid mode: report single allele instead of diploid pairs. Use for X/Y chromosomes in males or haploid organisms. <br />
</p>
</details>


## Outputs
<p name="STRGenotyper.detailed_csv_files">
        <b>STRGenotyper.detailed_csv_files</b><br />
        <i>Array[File]</i><br />
        Detailed per-read alignment results in CSV format, containing alignment scores, repeat counts, and read metadata for each alignment. Empty array if output_detailed_csv=false.
</p>
<p name="STRGenotyper.summary_csv_files">
        <b>STRGenotyper.summary_csv_files</b><br />
        <i>Array[File]</i><br />
        Per-locus summary statistics in CSV format, aggregating alignment results across all reads for each STR locus. Empty array if output_summary_csv=false.
</p>
<p name="STRGenotyper.genotypes_bed">
        <b>STRGenotyper.genotypes_bed</b><br />
        <i>File</i><br />
        Final genotype calls in BED format for visualization in genome browsers (IGV, UCSC). Contains chromosome, start, end, and genotype information
</p>
<p name="STRGenotyper.genotypes_vcf">
        <b>STRGenotyper.genotypes_vcf</b><br />
        <i>File</i><br />
        Final genotype calls in compressed VCF format with per-allele support counts (ADSP, ADFL). Compatible with standard VCF tools.
</p>
<p name="STRGenotyper.genotypes_vcf_index">
        <b>STRGenotyper.genotypes_vcf_index</b><br />
        <i>File</i><br />
        Tabix index for the genotypes VCF file
</p>

<hr />

> Generated using WDL AID (1.0.1)