# UG-High confidence region
## Version Log
* v1.0 Initial version (BIOIN-286)
* v1.1 Fixed a mistake in annotation of the exclusion regions: for details see analysis/220419_fix_high_confidence_annotations.ipynb. 
* v1.2 ug_cnv_lcr added. see description below.
* v1.2.1 annotation tag updated per 
* v1.3 HMER11 was split off GAP annotation (that annotates Ns in the reference +- 4 bases) 
* v1.4 ug_cnv_lcr updated: annotation tags were merged and additional regions (telomere-centromere) were added. 
* v2.1 UG-HCR updated for the data from RM6 ePCR formulation and basecalling 5.0. 


## Definition of `ug_lcr.bed` and `ug_hcr.bed`

Genomic areas where UG performance is poor are annotated in `ug_lcr.bed`. The regions and their annotations are: 

* HMER11 - Homopolymer runs of length 11 bp and above, padded with three bases to the left and two bases to the right around the genomic coordinates (total 17.7 Mb).  
* Tandem- Tandem repeats that are not spanned (end-to-end) by at least 10-reads on average between five samples (total 4.6 Mb). Tandem repeats in hg38 were annotated using trf requiring perfect match only.
* Coverage-Mappability- (total 111 Mb) 100bp windows resulting from one of the following analyses:
   - windows that, on average, are covered by at least 10 reads, but less than 5% of these reads are aligned with mapping quality >= 20.
   - windows with reads (mapq>=20) coverage that is highly variable between samples (std/mean > 0.55).
   - windows with reads (mapq>=20) coverage smaller than 10. 

Coverage-Mappability track was generated by collecting coverage statistics from 98 samples of mean coverage 50x.  

Total excluded BED size is 131.3 Mb. Majority of the excluded bases lie outside of GIAB high confidence regions where less than 1% is excluded from each GIAB sample high confidence region (version 4.2.1).  

`ug_hcr.bed` is the genomic complement of `ug_lcr.bed`. 

## Intersection of UG-HCR with various genomic regions

Detailed intersection statistics. Table 1: sizes of `ug_lcr` intersections with genome, exome and NIST HCRs in bp Table 2: per cent of the region excluded. GIAB samples high confidence region is NIST v4.2.1 high confidence region .  

tag | genome | exome | 
----|--------| -----|
| HMER | 17,719,256 | 508 |
| Tandem | 4,573,628 | 12,959 |
| Coverage-Mappability | 110,955,208 | 760,591 | 
| all_lcr | 131,334,542 | 768,900 | 

| tag                  | HG001      | HG002      | HG003      | HG004      | HG005      | HG006      | HG007      |
| -------------------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| HMER                 | 7,312,081  | 8,023,149  | 7,879,037  | 8,083,603  | 4,020,636  | 6,329,343  | 6,502,421  |
| Tandem               | 1,464,371  | 1,755,125  | 1,543,470  | 1,543,511  | 1,420,919  | 1,394,484  | 1,399,610  |
| Coverage-Mappability | 14,906,802 | 15,619,905 | 15,567,759 | 14,835,967 | 14,553,842 | 15,012,885 | 14,882,908 |
| all_lcr              | 23,419,815 | 25,101,766 | 24,756,209 | 24,225,272 | 19,784,405 | 22,518,391 | 22,559,871 |

 Table 1. Number of bases excluded from genome, exome and GIAB high confidence regions

| tag                  | exome | hg38 | HG001 | HG002 | HG003 | HG004 | HG005 | HG006 | HG007 |
| -------------------- | ----- | ---- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| HMER                 | 0.001 | 0.57 | 0.29  | 0.32  | 0.31  | 0.32  | 0.16  | 0.25  | 0.26  |
| Tandem               | 0.04  | 0.15 | 0.06  | 0.07  | 0.06  | 0.06  | 0.06  | 0.06  | 0.06  |
| Coverage-Mappability | 2.29  | 3.59 | 0.59  | 0.61  | 0.62  | 0.59  | 0.58  | 0.60  | 0.59  |
| all_lcr              | 2.32  | 4.25 | 0.93  | 0.99  | 0.98  | 0.96  | 0.79  | 0.89  | 0.90  |

Table 2. % of the region excluded from genome, exome and GIAB high confidence regions

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
