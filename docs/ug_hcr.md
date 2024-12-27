# UG-High confidence region
## Version Log
* v1.0 Initial version (BIOIN-286)
* v1.1 Fixed a mistake in annotation of the exclusion regions: for details see analysis/220419_fix_high_confidence_annotations.ipynb. 
* v1.2 ug_cnv_lcr added. see description below.
* v1.2.1 annotation tag updated per 
* v1.3 HMER11 was split off GAP annotation (that annotates Ns in the reference +- 4 bases) 
* v1.4 ug_cnv_lcr updated: annotation tags were merged and additional regions (telomere-centromere) were added. 
* v2.1 UG-HCR updated for the data from RM6 ePCR formulation and basecalling 5.0. 
* v2.1.1  ug_cnv_lcr updated: telomere annotation tags were updated to include chromosome edges.
* v2.1.2 sorted bed files based on chromosome order in reference genome
* v3.0 ug_hcr is developed based on gVCFs of multiple runs.

## Definition of `ug_lcr.bed` and `ug_hcr.bed`

Genomic areas where UG performance is consistently poor are annotated in ug_lcr.bed. The regions were defined as follows: 

gVCFs were generated from 56 GIAB samples sequenced to a coverage of ~40X (HG001-7). 

Run HCR was defined as any block (reference or variant) that has GQ>=20

Run HCRs were intersected and any block that occured in less than 35 samples was assigned to ug_lcr.bed

In addition, reference homopolymers of length 15 or longer were extended by three bases to the left and two bases to the right and also assigned to `ug_lcr.bed`. 

 Total excluded BED size is 272 Mb. Majority of the excluded bases lie outside of GIAB high confidence regions where less than 1% is excluded from each GIAB sample high confidence region (version 4.2.1).  

`ug_hcr.bed` are the regions complementary to `ug_lcr.bed`. 



## Intersection of UG-HCR with various genomic regions

| tag | HG001 | HG002 | HG003 | HG004 | HG005 | HG006 | HG007 | exome | WG |
| --- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | -- |
excluded bases | 21,913,516 | 23,639,890 | 23,175,722 | 22,500,237 | 19,207,155 | 21,368,823 | 21,324,483 | 1,050,817 | 272,782,058 |
% excluded bases | 0.77 | 0.83 | 0.81 | 0.79 | 0.68 | 0.75 | 0.75 | 0.03 | 8.48 |

## Description ug_cnv_lcr 

We defined regions where coverage biases and mappability issues are likely to generate noisy copy number calls (UG-CNV-LCR). These regions were defined by the following criteria based on the coverage profile of 90 unrelated samples sequenced to approximate coverage of 50x on the UG sequencer.  

* Coverage-Mappability: 
   * Low mappability - defined as 50 bp windows that, on average, are covered by at least 20 reads, but less than 10% of these reads are aligned with mapping quality > 20. 
   * High coverage variability - 50 bp windows with coverage that is highly variable between samples (std/mean > 0.5) 
   * Low coverage - defined as 50 bp windows that, on average, are covered by 10 reads or less. 
* Telomere_Centromere:  250kb at the end of each chromosome, 250kb flanking sequence around the centromeres 
* Clusters: regions that were empirically found to generate noisy CNV calls and thus were marked in the callset. 

SV calls that had larger than 50% overlap with excluded areas were marked as “UG-CNV-LCR” in the callset.  

The following table details each UG-CNV-LCR criteria by size

| tag                  | GIAB        | exome       | %GIAB        | %exome | %hg38 |
| -------------------- | ----------- | ------------ |-------------|--------|-------|
| Coverage-Mappability | 157,351,620 | 42,620,663   | 977,296 | 1.7% | 2.9% | 5.1% |
| Telomere-Centromere  | 66,802,541  | 5,151,741    | 148,638 | 0.2% | 0.4% | 2.2% |
| Clusters             | 7,087,267   | 1,595,921    | 73,661 | 0.06% | 0.2% | 0.2% |
