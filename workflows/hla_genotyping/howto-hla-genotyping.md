# Human Leukocyte Antigen (HLA) and Killer Immunoglobulin Receptor (KIR) Genotyping

The HLA genotyping pipeline uses T1K tool for identifying HLA and KIR alleles
The pipeline takes a CRAM/BAM file as input and produces genotyping results in text/TSV format.

### Overview
T1K (The ONE genotyper for KIR and HLA) is a computational tool that infers alleles for polymorphic genes including KIR and HLA. It calculates allele abundances based on RNA-seq/WES/WGS read alignments on provided allele reference sequences.


### Requirements
1. Input CRAM/BAM files and respective indexes
2. T1K coordinate files (pre-generated, versioned by database release)
   - HLA coordinate file: `hlaidx_dna_coord.fa`
   - KIR coordinate file: `kiridx_dna_coord.fa` (optional, only if running KIR genotyping)
  Download from `s3://ultimagen-workflow-resources-us-east-1/hla/t1k_index_hladb_v3.63.0.tar.gz`
3. Docker for t1k ultimagenomics/ugbio_t1k:1.28.0

### T1K Output Files

T1K produces the following output files:

#### HLA Genotyping Outputs
- `{base_file_name}_hla_genotype.tsv`: Main genotyping results
  - One row per HLA gene
  - Columns: gene name, number of different alleles, allele calls with abundance and quality scores
  - **Recommendation**: Ignore alleles with quality ≤ 0
- `{base_file_name}_hla_allele.tsv`: Detailed per-allele information
  - Per-allele abundance and coverage metrics

#### KIR Genotyping Outputs (if enabled)
- `{base_file_name}_kir_genotype.tsv`: KIR genotyping results (same format as HLA)
- `{base_file_name}_kir_allele.tsv`: Detailed KIR allele information

### Output Format Details

The genotype TSV file format:
```
gene_name  num_diff_alleles  allele_1  abundance_1  quality_1  allele_2  abundance_2  quality_2  secondary_alleles
```

- Missing or homozygous alleles: shown as `. 0 -1` as placeholders
- Secondary alleles: pipe-separated fields (allele;abundance;quality) that met abundance criteria but were filtered by tie-breaking

### Docker Image

The T1K Docker image includes:
- T1K v1.0.9 with CRAM support via htslib 1.15.1
- All required compression libraries (libdeflate, lzma, bz2)

### Running T1K 

If running T1K directly (not via the workflow):

```bash
# Set REF_PATH for CRAM decoding
export REF_PATH=/path/to/hg38/fasta/dir/

# HLA genotyping
perl run-t1k \
    -b sample.cram \
    -f hlaidx/hlaidx_dna_seq.fa \
    -c hlaidx_dna_coord.fa \
    --preset hla-wgs \
    --skipPostAnalysis \
    -t 16 \
    -o sample_hla \
    --od output/

# KIR genotyping
perl run-t1k \
    -b sample.cram \
    -f kiridx/kiridx_dna_seq.fa \
    -c kiridx_dna_coord.fa \
    --preset kir-wgs \
    --skipPostAnalysis \
    -t 16 \
    -o sample_kir \
    --od output/
```

---
