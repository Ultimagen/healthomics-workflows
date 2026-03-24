# FFPEv1.9 Model Features Summary

## Model Overview
The **FFPEv1.9** is an XGBoost model trained on FFPE tumor-normal paired samples to classify somatic single-nucleotide variants (SNVs). The model distinguishes true somatic mutations from germline variants and sequencing artifacts by learning FFPE-specific coverage patterns, SRSNV aggregated qualities and  Mapping properties.

**Total Features:** 58
**Prefix Convention:**
- `t_` = Tumor/Case sample features
- `n_` = Normal/Control sample features

---

## Feature Categories

### 1. Coverage & Variant Allele Frequency (VAF)
Critical metrics for variant quality and somatic status.

| Feature | Description | Notes |
|---------|-------------|-------|
| `{t,n}_dp` | **Depth** - Total read coverage at variant position | Higher depth = more confident calls |
| `{t,n}_vaf` | **Variant Allele Frequency** - Fraction of reads supporting ALT allele | Tumor VAF > Normal VAF suggests somatic variant |
| `{t,n}_raw_vaf` | **Raw VAF** - Unfiltered VAF before quality filters | Comparison with filtered VAF shows filtering impact |
| `{t,n}_alt_reads` | **ALT read count** - Number of reads supporting variant | Absolute count supporting the variant |
| `{t,n}_pass_alt_reads` | **Passing ALT reads** - ALT reads passing quality filters | High-quality evidence for variant |

**Key Insight:** True somatic variants typically show high tumor VAF with low/zero normal VAF.

---

### 2. SRSNV Model Quality (MQUAL)
Per-read variant confidence scores from the Single-Read SNV (SRSNV) model. MQUAL represents the model's confidence that each read truly supports the variant, aggregated across all supporting reads.

| Feature | Description | Why It Matters | Valid Range |
|---------|-------------|----------------|-------------|
| `{t,n}_mqual_mean` | **Mean SRSNV quality** | Average per-read variant confidence | 0-10 (typically 2-5) |
| `{t,n}_mqual_max` | **Max SRSNV quality** | Highest confidence supporting read | 0-10 |
| `{t,n}_mqual_min` | **Min SRSNV quality** | Lowest confidence supporting read | 0-10 (typically 0.5-2) |

**Key Insight:** True variants have higher, more consistent MQUAL scores. Low MQUAL suggests reads with ambiguous or low-confidence variant support. Note: MQUAL is NOT mapping quality (MAPQ)—it's the SRSNV model's per-read confidence score. 

---

### 3. Pileup Context (5-base window)
Read pileup counts of reference and non-reference alleles in the variant's vicinity (±2 bases). These features capture the local sequence context and allele distribution around the variant position.

**Position Mapping:** The feature indices 0-4 correspond to positions **-2, -1, 0, +1, +2** relative to the variant:
- `ref0` / `nonref0` = position -2 (2 bases upstream)
- `ref1` / `nonref1` = position -1 (1 base upstream)
- `ref2` / `nonref2` = **position 0 (the variant position itself)**
- `ref3` / `nonref3` = position +1 (1 base downstream)
- `ref4` / `nonref4` = position +2 (2 bases downstream)

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_ref0-4` | **REF allele pileup counts** | Number of reads supporting REF at each position in the 5-base window |
| `{t,n}_nonref0-4` | **Non-REF allele pileup counts** | Number of reads supporting non-REF alleles at each position |

**Key Insight:** The central position (index 2 = ref2/nonref2) represents the variant position itself and is particularly informative for variant quality. High nonref counts at flanking positions may indicate alignment errors or homopolymer-related issues.

---

### 4. Strand Bias
Directionality of supporting reads.

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_reverse_count` | **Reverse strand read count** | Reads aligned to reverse strand |
| `{t,n}_forward_count` | **Forward strand read count** | Reads aligned to forward strand |

**Key Insight:** True variants appear on both strands roughly equally. Strong bias suggests strand-specific artifacts (e.g., FFPE damage, oxidation).

---

### 5. Soft-Clipped Reads at Read Boundaries
Count of reads with soft-clipping at read start/end positions. Soft-clipping indicates poor alignment or structural variation at read boundaries.

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_scst_num_reads` | **Soft-clipped at start** - Reads with soft-clipping at read start | High counts suggest alignment artifacts or structural variants |
| `{t,n}_sced_num_reads` | **Soft-clipped at end** - Reads with soft-clipping at read end | High counts indicate poor alignment quality at read terminus |

**Key Insight:** High scst/sced counts suggest alignment issues, structural variation, or low-quality read boundaries → likely artifacts or complex regions.

---

### 6. Mapping Ambiguity

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_map0_count` | **MAPQ=0 reads** - Count of reads with mapping quality = 0 | High count = repetitive region or multi-mapping reads |
| `t_tr_distance` | **Tandem repeat distance** | Distance to nearest tandem repeat region |

**Key Insight:** True variants in unique regions have low map0 counts. High map0 = alignment uncertainty. Variants near tandem repeats (low tr_distance) are more prone to alignment errors.

---

### 7. Edit Distance
Measure of read-to-reference similarity.

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_edist_mean` | **Mean edit distance** - Average mismatches/indels per read | Lower = better alignment quality |
| `{t,n}_edist_max` | **Max edit distance** | Worst-aligned supporting read |
| `{t,n}_edist_min` | **Min edit distance** | Best-aligned supporting read |

**Key Insight:** True variants have low edit distance (few mismatches). High edist suggests noisy/problematic region.

---

### 8. Read Length (RL)
Supporting read size distribution.

| Feature | Description | Why It Matters |
|---------|-------------|----------------|
| `{t,n}_rl_mean` | **Mean read length** | Average size of supporting reads |
| `{t,n}_rl_max` | **Max read length** | Longest supporting read |
| `{t,n}_rl_min` | **Min read length** | Shortest supporting read |

**Key Insight:** Short reads can indicate fragmentation artifacts (common in FFPE). Very long reads = high-quality DNA.

---

## Feature Interactions for Somatic Calling

### True Somatic Variant Profile:
✅ **High tumor VAF**, low/zero normal VAF

✅ **High MQUAL** (confident per-read variant support)

✅ **Low edit distance** (good alignment quality)

✅ **Balanced strand counts** (no bias)

✅ **Moderate pileup context** (not in extreme homopolymer runs)

✅ **Low scst/sced** (minimal soft-clipping)


### False Positive (Artifact) Profile:
❌ **Elevated normal VAF** (germline or systematic error)

❌ **Low MQUAL** or high map0 (alignment/quality issue)

❌ **Strand bias** (forward or reverse only)

❌ **High edit distance** (noisy region)

❌ **Abnormal pileup patterns** (unexpected local reference/non-reference read-count distributions can indicate sequencing or alignment artifacts)

❌ **Short reads** (FFPE fragmentation)

❌ **High scst/sced** (soft-clipping artifacts)


---

## Top 10 Most Influential Features (by XGBoost Gain)

Based on XGBoost feature importance analysis (higher gain = more influential in model decisions):

| Rank | Feature | Importance | Category | Description |
|------|---------|------------|----------|-------------|
| 1 | **n_vaf** | 1569.6 | Coverage & VAF | Normal VAF - most critical for germline filtering |
| 2 | **t_vaf** | 995.6 | Coverage & VAF | Tumor VAF - primary somatic signal |
| 3 | **n_raw_vaf** | 866.2 | Coverage & VAF | Normal raw VAF - unfiltered germline detection |
| 4 | **t_pass_alt_reads** | 681.9 | Coverage & VAF | High-quality tumor ALT evidence |
| 5 | **t_ref2** | 612.8 | Pileup Context | Tumor REF pileup at variant position (index 2 = position 0) |
| 6 | **t_map0_count** | 535.4 | Mapping Ambiguity | Tumor reads with MAPQ=0 |
| 7 | **n_nonref2** | 482.7 | Pileup Context | Normal non-REF pileup at variant position (index 2 = position 0) |
| 8 | **t_mqual_min** | 481.7 | SRSNV Quality | Worst tumor SRSNV quality - confidence floor |
| 9 | **n_ref2** | 385.0 | Pileup Context | Normal REF pileup at variant position (index 2 = position 0) |
| 10 | **t_edist_mean** | 367.2 | Edit Distance | Average tumor edit distance - alignment quality |

**Key Insight:** The model heavily relies on **VAF comparison** (n_vaf, t_vaf, n_raw_vaf) to distinguish somatic from germline variants, then uses **pileup context** (ref2, nonref2) and **quality metrics** (SRSNV quality, edit distance) to filter sequencing artifacts.

---

## Model Training Strategy

The **FFPEv1.9** model was trained on FFPE samples with known somatic variants to learn:
1. **FFPE-specific artifacts** (C>T/G>A damage patterns)
2. **Flow-based error signatures** (homopolymer indels)
3. **Somatic vs. germline discrimination**
4. **Context-dependent error rates**

**Scoring:** XGBoost outputs probability (0-1) that a variant is a **true somatic variant**.
- Score ≥ 0.6 → High confidence somatic (PASS)
- Score < 0.6 → Likely artifact or germline (FILTERED)

---

*Generated: 2026-03-02*
*Model: gs://concordanz/somatic_snvfind/sfm.xgb_model.FFPEV1.9.json*
