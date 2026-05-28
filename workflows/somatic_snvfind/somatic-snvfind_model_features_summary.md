# Somatic SNV Find Model Features Summary (FFPEv1.9 & FreshFrozenv1.21)

## Model Overview
The **FFPEv1.9** and **FreshFrozenv1.21** are XGBoost models trained on tumor-normal paired samples to classify somatic single-nucleotide variants (SNVs). Each model distinguishes true somatic mutations from germline variants and sequencing artifacts by learning sample-type-specific coverage patterns, SRSNV aggregated qualities and Mapping properties.

- **FFPEv1.9** — trained on FFPE tumor samples vs. BuffyCoat normal.
- **FreshFrozenv1.21** — trained on Fresh Frozen tumor samples vs. BuffyCoat normal.

**Feature counts:**
- Common to both models: **40**
- Private to **FFPEv1.9**: **18**
- Private to **FreshFrozenv1.21**: **13**

**Prefix Convention:**
- `t_` = Tumor/Case sample features
- `n_` = Normal/Control sample features

**Model availability legend (used in tables below):**
- **Common** — used by both FFPE and FreshFrozen models
- **FFPE only** — used by FFPEv1.9 only
- **FF only** — used by FreshFrozenv1.21 only

---

## Feature Categories

### 1. Coverage & Variant Allele Frequency (VAF)
Critical metrics for variant quality and somatic status.

| Feature | Description | Model | Notes |
|---------|-------------|-------|-------|
| `{t,n}_dp` | **Depth** - Total read coverage at variant position | Common | Higher depth = more confident calls |
| `{t,n}_vaf` | **Variant Allele Frequency** - Fraction of reads supporting ALT allele | Common | Tumor VAF > Normal VAF suggests somatic variant |
| `n_raw_vaf` | **Raw VAF (normal)** - Unfiltered VAF before quality filters | Common | Comparison with filtered VAF shows filtering impact |
| `t_raw_vaf` | **Raw VAF (tumor)** - Unfiltered VAF before quality filters | **FF only** | Comparison with filtered VAF shows filtering impact |
| `{t,n}_alt_reads` | **ALT read count** - Number of reads supporting variant | Common | Absolute count supporting the variant |
| `{t,n}_pass_alt_reads` | **Passing ALT reads** - ALT reads passing quality filters | Common | High-quality evidence for variant |
| `{t,n}_dp_filt` | **Filtered depth** - Total read coverage after applying quality filters | **FF only** | Effective coverage available to the caller after filters |

**Key Insight:** True somatic variants typically show high tumor VAF with low/zero normal VAF.

---

### 2. SRSNV Model Quality (MQUAL)
Per-read variant confidence scores from the Single-Read SNV (SRSNV) model. MQUAL represents the model's confidence that each read truly supports the variant, aggregated across all supporting reads.

| Feature | Description | Model | Why It Matters | Valid Range |
|---------|-------------|-------|----------------|-------------|
| `{t,n}_mqual_mean` | **Mean SRSNV quality** | Common | Average per-read variant confidence | 0-10 (typically 2-5) |
| `{t,n}_mqual_max` | **Max SRSNV quality** | Common | Highest confidence supporting read | 0-10 |
| `{t,n}_mqual_min` | **Min SRSNV quality** | Common | Lowest confidence supporting read | 0-10 (typically 0.5-2) |

**Key Insight:** True variants have higher, more consistent MQUAL scores. Low MQUAL suggests reads with ambiguous or low-confidence variant support. Note: MQUAL is NOT mapping quality (MAPQ)—it's the SRSNV model's per-read confidence score.

---

### 3. Mapping Quality (MAPQ)
Aligner-reported mapping quality statistics across supporting reads. Distinct from SRSNV MQUAL.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_mapq_mean` | **Mean MAPQ** - Average aligner mapping quality across supporting reads | **FF only** | Lower mean MAPQ → more ambiguous alignment |
| `{t,n}_mapq_max` | **Max MAPQ** | **FF only** | Best-mapped supporting read |
| `{t,n}_mapq_min` | **Min MAPQ** | **FF only** | Worst-mapped supporting read |
| `{t,n}_map0_count` | **MAPQ=0 reads** - Count of reads with mapping quality = 0 | **FFPE only** | High count = repetitive region or multi-mapping reads |

**Key Insight:** True variants in unique regions show consistently high MAPQ and low map0 counts. Low/variable MAPQ or high map0 across supporting reads indicates alignment uncertainty.

---

### 4. Pileup Context (5-base window)
Read pileup counts of reference and non-reference alleles in the variant's vicinity (±2 bases). These features capture the local sequence context and allele distribution around the variant position.

**Position Mapping:** The feature indices 0-4 correspond to positions **-2, -1, 0, +1, +2** relative to the variant:
- `ref0` / `nonref0` = position -2 (2 bases upstream)
- `ref1` / `nonref1` = position -1 (1 base upstream)
- `ref2` / `nonref2` = **position 0 (the variant position itself)**
- `ref3` / `nonref3` = position +1 (1 base downstream)
- `ref4` / `nonref4` = position +2 (2 bases downstream)

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_ref0-4` | **REF allele pileup counts** | Common | Number of reads supporting REF at each position in the 5-base window |
| `{t,n}_nonref0-4` | **Non-REF allele pileup counts** | Common | Number of reads supporting non-REF alleles at each position |

**Key Insight:** The central position (index 2 = ref2/nonref2) represents the variant position itself and is particularly informative for variant quality. High nonref counts at flanking positions may indicate alignment errors or homopolymer-related issues.

---

### 5. Strand Bias
Directionality of supporting reads.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_reverse_count` | **Reverse strand read count** | Common | Reads aligned to reverse strand |
| `{t,n}_forward_count` | **Forward strand read count** | Common | Reads aligned to forward strand |

**Key Insight:** True variants appear on both strands roughly equally. Strong bias suggests strand-specific artifacts (e.g., FFPE damage, oxidation).

---

### 6. Soft-Clipped Reads at Read Boundaries
Count of reads with soft-clipping at read start/end positions. Soft-clipping indicates poor alignment or structural variation at read boundaries.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_scst_num_reads` | **Soft-clipped at start** - Reads with soft-clipping at read start | **FFPE only** | High counts suggest alignment artifacts or structural variants |
| `{t,n}_sced_num_reads` | **Soft-clipped at end** - Reads with soft-clipping at read end | **FFPE only** | High counts indicate poor alignment quality at read terminus |

**Key Insight:** High scst/sced counts suggest alignment issues, structural variation, or low-quality read boundaries → likely artifacts or complex regions.

---

### 7. Tandem Repeats

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `t_tr_distance` | **Tandem repeat distance** - Distance to nearest tandem repeat region | Common | Variants near tandem repeats are more prone to alignment errors |
| `t_tr_length` | **Tandem repeat length** - Length of the nearest tandem repeat | **FF only** | Longer repeats increase alignment ambiguity |
| `t_tr_seq_unit_length` | **Tandem repeat unit length** - Length of the repeat motif | **FF only** | Shorter motifs (mono/di-nucleotide) are most error-prone |

**Key Insight:** Variants near tandem repeats (low tr_distance, high tr_length, short tr_seq_unit_length) are more prone to alignment errors.

---

### 8. Edit Distance
Measure of read-to-reference similarity.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_edist_mean` | **Mean edit distance** - Average mismatches/indels per read | **FFPE only** | Lower = better alignment quality |
| `{t,n}_edist_max` | **Max edit distance** | **FFPE only** | Worst-aligned supporting read |
| `{t,n}_edist_min` | **Min edit distance** | **FFPE only** | Best-aligned supporting read |

**Key Insight:** True variants have low edit distance (few mismatches). High edist suggests noisy/problematic region.

---

### 9. Read Length (RL)
Supporting read size distribution.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `{t,n}_rl_mean` | **Mean read length** | **FFPE only** | Average size of supporting reads |
| `{t,n}_rl_max` | **Max read length** | **FFPE only** | Longest supporting read |
| `{t,n}_rl_min` | **Min read length** | **FFPE only** | Shortest supporting read |

**Key Insight:** Short reads can indicate fragmentation artifacts (common in FFPE). Very long reads = high-quality DNA.

---

### 10. Homopolymer Context
Local homopolymer length around the variant — flow-based sequencing is most error-prone in long homopolymer runs.

| Feature | Description | Model | Why It Matters |
|---------|-------------|-------|----------------|
| `x_hmer_ref` | **REF homopolymer length** at the variant position | **FF only** | Long REF homopolymers are error-prone for flow-based indels |
| `x_hmer_alt` | **ALT homopolymer length** induced by the variant | **FF only** | Variants creating/extending homopolymers are more likely to be artifacts |

**Key Insight:** Variants in/adjacent to long homopolymers carry higher artifact risk; the model uses these features to down-weight such calls.

---

## Per-Model Feature Inventory

### Common features (40)
`n_alt_reads`, `n_dp`, `n_forward_count`, `n_mqual_max`, `n_mqual_mean`, `n_mqual_min`,
`n_nonref0`, `n_nonref1`, `n_nonref2`, `n_nonref3`, `n_nonref4`,
`n_pass_alt_reads`, `n_raw_vaf`, `n_ref0`, `n_ref1`, `n_ref2`, `n_ref3`, `n_ref4`,
`n_reverse_count`, `n_vaf`,
`t_alt_reads`, `t_dp`, `t_forward_count`, `t_mqual_max`, `t_mqual_mean`, `t_mqual_min`,
`t_nonref0`, `t_nonref1`, `t_nonref2`, `t_nonref3`, `t_nonref4`,
`t_pass_alt_reads`, `t_ref0`, `t_ref1`, `t_ref2`, `t_ref3`, `t_ref4`,
`t_reverse_count`, `t_tr_distance`, `t_vaf`

### Private to FFPEv1.9 (18)
`n_edist_max`, `n_edist_mean`, `n_edist_min`,
`n_map0_count`,
`n_rl_max`, `n_rl_mean`, `n_rl_min`,
`n_sced_num_reads`, `n_scst_num_reads`,
`t_edist_max`, `t_edist_mean`, `t_edist_min`,
`t_map0_count`,
`t_rl_max`, `t_rl_mean`, `t_rl_min`,
`t_sced_num_reads`, `t_scst_num_reads`

### Private to FreshFrozenv1.21 (13)
`n_dp_filt`, `n_mapq_max`, `n_mapq_mean`, `n_mapq_min`,
`t_dp_filt`, `t_mapq_max`, `t_mapq_mean`, `t_mapq_min`,
`t_raw_vaf`, `t_tr_length`, `t_tr_seq_unit_length`,
`x_hmer_alt`, `x_hmer_ref`

---

## Feature Interactions for Somatic Calling

### True Somatic Variant Profile:
- **High tumor VAF**, low/zero normal VAF
- **High MQUAL** (confident per-read variant support)
- **High MAPQ** (FreshFrozen model — unambiguous alignment)
- **Low edit distance** (FFPE model — good alignment quality)
- **Balanced strand counts** (no bias)
- **Moderate pileup context** (not in extreme homopolymer runs)
- **Low scst/sced** (FFPE model — minimal soft-clipping)
- **Short hmer context** (FreshFrozen model — variant not inside a long homopolymer)

### False Positive (Artifact) Profile:
- **Elevated normal VAF** (germline or systematic error)
- **Low MQUAL** or high map0 (FFPE) / low MAPQ (FreshFrozen) — alignment/quality issue
- **Strand bias** (forward or reverse only)
- **High edit distance** (FFPE — noisy region)
- **Abnormal pileup patterns** (unexpected local reference/non-reference read-count distributions can indicate sequencing or alignment artifacts)
- **Short reads** (FFPE — fragmentation)
- **High scst/sced** (FFPE — soft-clipping artifacts)
- **Long hmer context** (FreshFrozen — flow-based homopolymer artifact)
- **Long tandem repeat / short repeat unit** (FreshFrozen — alignment ambiguity)

---

## Top 10 Most Influential Features — FFPEv1.9 (by XGBoost Gain)

Based on XGBoost feature importance analysis (higher gain = more influential in model decisions):

| Rank | Feature | Importance | Category | Description |
|------|---------|------------|----------|-------------|
| 1 | **n_vaf** | 1569.6 | Coverage & VAF | Normal VAF - most critical for germline filtering |
| 2 | **t_vaf** | 995.6 | Coverage & VAF | Tumor VAF - primary somatic signal |
| 3 | **n_raw_vaf** | 866.2 | Coverage & VAF | Normal raw VAF - unfiltered germline detection |
| 4 | **t_pass_alt_reads** | 681.9 | Coverage & VAF | High-quality tumor ALT evidence |
| 5 | **t_ref2** | 612.8 | Pileup Context | Tumor REF pileup at variant position (index 2 = position 0) |
| 6 | **t_map0_count** *(FFPE only)* | 535.4 | Mapping Ambiguity | Tumor reads with MAPQ=0 |
| 7 | **n_nonref2** | 482.7 | Pileup Context | Normal non-REF pileup at variant position (index 2 = position 0) |
| 8 | **t_mqual_min** | 481.7 | SRSNV Quality | Worst tumor SRSNV quality - confidence floor |
| 9 | **n_ref2** | 385.0 | Pileup Context | Normal REF pileup at variant position (index 2 = position 0) |
| 10 | **t_edist_mean** *(FFPE only)* | 367.2 | Edit Distance | Average tumor edit distance - alignment quality |

**Key Insight:** The model heavily relies on **VAF comparison** (n_vaf, t_vaf, n_raw_vaf) to distinguish somatic from germline variants, then uses **pileup context** (ref2, nonref2) and **quality metrics** (SRSNV quality, edit distance) to filter sequencing artifacts.

---

## Top 10 Most Influential Features — FreshFrozenv1.21 (by XGBoost Gain)

Based on XGBoost feature importance analysis (higher gain = more influential in model decisions):

| Rank | Feature | Importance | Category | Description |
|------|---------|------------|----------|-------------|
| 1 | **n_alt_reads** | 80393.4 | Coverage & VAF | Normal ALT read count - dominant germline / contamination signal |
| 2 | **n_mapq_min** *(FF only)* | 60232.3 | Mapping Quality | Worst normal MAPQ - flags low-confidence normal alignments |
| 3 | **t_pass_alt_reads** | 7360.8 | Coverage & VAF | High-quality tumor ALT evidence |
| 4 | **n_vaf** | 7073.5 | Coverage & VAF | Normal VAF - core germline filter |
| 5 | **n_nonref2** | 1772.2 | Pileup Context | Normal non-REF pileup at variant position (index 2 = position 0) |
| 6 | **t_mapq_min** *(FF only)* | 1439.3 | Mapping Quality | Worst tumor MAPQ - alignment confidence floor |
| 7 | **n_pass_alt_reads** | 1265.8 | Coverage & VAF | High-quality normal ALT evidence |
| 8 | **t_mapq_max** *(FF only)* | 735.0 | Mapping Quality | Best tumor MAPQ |
| 9 | **t_mapq_mean** *(FF only)* | 721.9 | Mapping Quality | Average tumor MAPQ |
| 10 | **t_mqual_min** | 453.6 | SRSNV Quality | Worst tumor SRSNV quality - confidence floor |

**Key Insight:** Compared to FFPEv1.9, the FreshFrozen model places much heavier weight on **MAPQ statistics** (4 of the top 10 features are MAPQ-based) and on **normal-sample evidence** (`n_alt_reads`, `n_vaf`, `n_pass_alt_reads`, `n_nonref2`) for germline/artifact rejection. SRSNV `t_mqual_min` remains a shared top-tier signal.

---

## Scoring:
Each XGBoost model outputs probability (0-1) that a variant is a true somatic variant. The PASS threshold is model-specific:

| Model | PASS threshold (`xgb_proba_threshold`) |
|-------|----------------------------------------|
| **FFPEv1.9** | **≥ 0.8** |
| **FreshFrozenv1.21** | **≥ 0.95** |

Variants with score below the threshold are FILTERED (likely artifact or germline). Thresholds are configured per use-case in the SomaticSNVfind input templates.

---

*Generated: 2026-03-02*
*Updated: 2026-05-24 — added FreshFrozenv1.21 features and per-model markings*
*Models:*
- *FFPE: gs://concordanz/somatic_snvfind/sfm.xgb_model.FFPEV1.9.json*
- *FreshFrozen: gs://concordanz/somatic_snvfind/somatic_snvfind.xgb_model.UGSB4v1.21_FreshFrozen.json*
