# CLAUDE.md — KD Transcriptomics Pipeline
## Sanjay Senthilnathan | Kawasaki Disease Subgroup Meta-Analysis

---

## Project Overview

This is a bioinformatics research pipeline for a **peer-review-targeted publication** on transcriptomic subgroup identification in Kawasaki Disease (KD). The pipeline performs:

1. Multi-cohort microarray data loading and preprocessing
2. Batch effect correction via ComBat (sva package)
3. Unsupervised consensus clustering (ConsensusClusterPlus)
4. Differential gene expression analysis (limma)
5. Pathway enrichment analysis (GSEA via fgsea, GSVA)
6. Clinical phenotype association testing

**Target journals:** Frontiers in Pediatrics, Scientific Reports, PLOS ONE, BMJ Open

---

## Repository Structure

```
kd_transcriptomics/
├── .git/
├── .gitignore
├── kd_transcriptomics.Rproj
├── CLAUDE.md                       # This file
├── README.md
│
├── pipeline.r                      # MAIN PIPELINE — entry point (at ROOT, NOT in r_code/)
│
├── r_code/
│   ├── illumina_loader.R           # Loads raw Illumina BeadStudio text files
│   ├── meta_filter.R               # Filters expr/pval matrices by GEO metadata
│   ├── process_study.R             # Collapses probes to Entrez IDs per study
│   ├── clean_dx.R                  # Maps diagnosis strings to KD/Control/Uncertain
│   └── gene_collapse.R             # DEAD CODE — do not use, scheduled for deletion
│
├── transcriptome_data/
│   ├── uncompressed/               # Extracted GEO series matrix files + raw files
│   ├── compressed/                 # GSE63881 non-normalized gzipped file
│   └── platforms/
│       ├── v4/                     # NEEDS to contain GPL10558.txt (extracted from root)
│       └── v3/                     # NEEDS to contain GPL6947.txt (extracted from root)
│
├── GPL10558.soft.gz                # Illumina HT-12 V4 platform — COMPRESSED at ROOT
├── GPL6480.soft.gz                 # Agilent platform — COMPRESSED at ROOT
├── GSE68004_series_matrix.txt.gz
├── GSE18606_series_matrix.txt.gz
│
├── final_results/
├── kd_clusters/
├── old_results/
├── Poster/
├── Biogeneous/
├── ISAS_forms/
├── ISEF_submissions/
├── Senthilnathan_2026_final.pdf
└── KD sci fair speech script.pdf
```

---

## Critical Setup Issue — Platform Files Path Mismatch

The code in pipeline.r expects platform annotation files at:
- `transcriptome_data/platforms/v4/GPL10558.txt`
- `transcriptome_data/platforms/v3/GPL6947.txt`

But the actual files are at the repo ROOT as compressed files:
- `GPL10558.soft.gz`
- `GPL6480.soft.gz`

**Fix this BEFORE running pipeline.r (Option A recommended):**

### Option A: Extract and organize files
```bash
mkdir -p transcriptome_data/platforms/v4
mkdir -p transcriptome_data/platforms/v3
gunzip -c GPL10558.soft.gz > transcriptome_data/platforms/v4/GPL10558.txt
gunzip -c GPL6480.soft.gz  > transcriptome_data/platforms/v3/GPL6947.txt
```

### Option B: Update pipeline.r to read from root
Change lines ~181-188 in pipeline.r:
```r
# FROM:
gplv4 <- read.delim("transcriptome_data/platforms/v4/GPL10558.txt", ...)
gplv3 <- read.delim("transcriptome_data/platforms/v3/GPL6947.txt", ...)

# TO:
gplv4 <- read.delim("GPL10558.soft.gz", comment.char = "#", stringsAsFactors = FALSE)
gplv3 <- read.delim("GPL6480.soft.gz",  comment.char = "#", stringsAsFactors = FALSE)
```

---

## Datasets

### Currently Integrated (4 datasets, ~559 samples)

| GSE ID   | Platform  | Center           | Scale on GEO                    | IVIG labels | CAA labels |
|----------|-----------|------------------|---------------------------------|-------------|------------|
| GSE73461 | HT-12 V4  | UK (Kaforou lab) | LINEAR (range: -58 to 68529)    | No          | No         |
| GSE73462 | HT-12 V3  | UK (Kaforou lab) | LOG2 (range: 4.19 to 15.3)      | No          | No         |
| GSE73463 | HT-12 V4  | UK (Kaforou lab) | LOG2 (range: -24.88 to 16)      | No          | No         |
| GSE63881 | HT-12 V4  | Burns lab UCSD   | LINEAR (raw Genome Studio)      | YES         | YES        |

**Critical:** GSE73461/62/63 are discovery + validation splits of ONE study (Wright/Levin/Kaforou),
not three independent studies. GSE73462 and 73463 had batch effects removed before upload to GEO.
This must be stated in the Methods section.

### Planned Addition (Priority 1)

| GSE ID   | Platform  | Center                               | Scale on GEO | KD n           | IVIG labels        | CAA labels |
|----------|-----------|--------------------------------------|--------------|----------------|--------------------|------------|
| GSE68004 | HT-12 V4  | Nationwide Children's / Baylor       | LINEAR       | 89 cKD + 13 iKD| In paper supp only | Same       |

IVIG/CAA labels for GSE68004 are NOT in GEO characteristics fields. They are in the supplementary
table of PMID 29813106 (PLOS One 2018). Email corresponding authors Ramilo/Mejias at Nationwide
Children's Hospital to request the per-sample annotation file.

---

## Known Bugs — Fix In Priority Order

### BUG 0: Platform file paths don't match actual locations [CRASHES PIPELINE]
**File:** pipeline.r lines ~181-188
**Problem:** Code reads from transcriptome_data/platforms/ but files are compressed at repo root.
**Fix:** See Critical Setup Issue section above. Do this FIRST.

### BUG 1: `res61` used before creation [CRASHES PIPELINE]
**File:** pipeline.r, CLEANING section (around line 100)
**Problem:** `res61$metadata <- res61_meta_clean[match(colnames(res61$expr), ...)]`
appears before `res61 <- filter_by_metadata(...)`. res61 does not exist yet.
**Fix:** Delete the premature res61$metadata assignment line. The filter_by_metadata
block 30 lines later handles this correctly and should not be changed.

### BUG 2: `master_expr_filtered` never defined [CRASHES PIPELINE]
**File:** pipeline.r, line ~802
**Problem:** `colnames(master_expr_filtered)` is referenced in the cluster-vs-controls
GSEA block but this variable is never created anywhere in the pipeline.
The pipeline has master_expr_final, master_norm, and master_expr_combat but not
master_expr_filtered.
**Fix:** Replace all instances of master_expr_filtered with master_expr_combat.

### BUG 3: No set.seed() before consensus clustering [RESULTS NOT REPRODUCIBLE]
**File:** pipeline.r, before ConsensusClusterPlus() call (around line 393)
**Problem:** ConsensusClusterPlus uses stochastic resampling (pItem=0.8, reps=1000).
Without a seed every run produces different cluster assignments. All figures, gene
lists, and p-values in the paper are run-specific and cannot be reproduced by anyone.
**Fix:** Add `set.seed(42)` on the line immediately before ConsensusClusterPlus().

### BUG 4: filter_by_metadata output immediately overwritten [PROBE FILTERING LOST]
**File:** pipeline.r, GSE73461/62/63 processing blocks (lines ~207-251)
**Problem:** res61 <- filter_by_metadata(...) correctly filters probes detected in
>50% of samples. But the very next line res61$expr <- expr73461_log2[, common61]
replaces the filtered matrix with the full unfiltered matrix. GSE63881 correctly
keeps its filtered matrix. Result: 3 of 4 datasets silently skip probe filtering.
**Fix:** Remove the three lines that reassign res6X$expr after filter_by_metadata.
Keep only the metadata realignment line:
```r
res61 <- filter_by_metadata("GSE73461", expr73461_log2, pval73461, threshold = 0.5)
res61$metadata <- res61_meta_clean[intersect(colnames(res61$expr), rownames(res61_meta_clean)), ]
# Do NOT add res61$expr <- ... after this
```

### BUG 5: rownames() called on a list object [RETURNS NULL SILENTLY]
**File:** pipeline.r, probe_ids section (around line 170)
**Problem:** probe_ids61 <- rownames(res61) — res61 is a list, not a matrix.
rownames() on a list returns NULL. These four variables are also never used anywhere
downstream.
**Fix:** Delete all four probe_ids lines entirely.

### BUG 6: Volcano plots use raw P.Value instead of adj.P.Val [SCIENTIFIC ERROR]
**File:** pipeline.r, make_poster_volcano() function (around line 685)
**Problem:** y = -log10(P.Value) and coloring threshold uses df$P.Value < p_thresh.
With n=559 samples, raw p-values are extremely small due to sample size alone, not
effect size. The paper states FDR=0.05 is used but the figures show uncorrected values.
**Fix:** Change to y = -log10(adj.P.Val) and color threshold to adj.P.Val < 0.05.

### BUG 7: GSEA duplicate-symbol deduplication is biased [SCIENTIFIC ERROR]
**File:** pipeline.r, GSEA loop (around line 742)
**Problem:** ranks <- ranks[!duplicated(names(ranks))] keeps the FIRST duplicate,
which is whichever Entrez ID sorts first. This biases GSEA toward lower/older
Entrez IDs (better-characterized genes), not the most informative probe.
**Fix:** For each duplicated symbol, keep the entry with highest absolute t-statistic:
```r
ranks_df <- data.frame(symbol = names(ranks), t = ranks)
ranks_df <- ranks_df[order(-abs(ranks_df$t)), ]
ranks_df <- ranks_df[!duplicated(ranks_df$symbol), ]
ranks <- setNames(ranks_df$t, ranks_df$symbol)
ranks <- sort(ranks, decreasing = TRUE)
```

---

## Normalization — The Critical Scientific Issue

### Verified scale of each dataset
```
GSE73461: range -58 to 68529  → LINEAR  (floor at 1, log2, within-study quantile normalize)
GSE73462: range  4.19 to 15.3 → LOG2    (already RSN normalized — no transform)
GSE73463: range -24.88 to 16  → LOG2    (already RSN normalized — clip negatives to 0 only)
GSE63881: raw linear Genome Studio export (floor at 1, log2, within-study quantile normalize)
GSE68004: raw linear GenomeStudio avg norm (floor at 1, log2, within-study quantile normalize)
```

### Why log2(x) not log2(x+1)
log2(x+1) is a convention for RNA-seq COUNT data where zero is a true observed value.
Microarray intensity is continuous — zeros are instrument noise, not true absence.
Floor at 1 so log2(1)=0, then log2() directly. The +1 pseudocount introduces a small
upward bias across the low-intensity range and should not be used for microarray data.

### Why within-study quantile normalization before ComBat
ComBat assumes its input has already been normalized within each batch.
Within-study QN before ComBat is the standard approach.
Do NOT run a second quantile normalization on the pooled merged matrix after ComBat.
Remove the pre-ComBat cross-study normalize.quantiles() call currently in pipeline.r
(around line 345) — this is double-normalizing and removes biological signal.

### Correct preprocessing code per dataset
```r
# GSE73461 — raw linear
expr73461_lin[expr73461_lin < 1] <- 1
expr73461_log2 <- log2(expr73461_lin)
expr73461_log2 <- normalize.quantiles(as.matrix(expr73461_log2))
rownames(expr73461_log2) <- rownames(expr73461_lin)
colnames(expr73461_log2) <- colnames(expr73461_lin)
rm(expr73461_lin); gc()

# GSE73462 — already log2 RSN
expr73462_log2 <- gse73462$expr_lin

# GSE73463 — already log2 RSN, has -24 artifacts
expr73463_log2 <- gse73463$expr_lin
expr73463_log2[expr73463_log2 < 0] <- 0

# GSE63881 — raw linear
expr81_lin[expr81_lin < 1] <- 1
expr81_log2 <- log2(expr81_lin)
expr81_log2 <- normalize.quantiles(as.matrix(expr81_log2))
rownames(expr81_log2) <- rownames(expr81_lin)
colnames(expr81_log2) <- colnames(expr81_lin)
rm(expr81_lin); gc()

# GSE68004 (when added) — same treatment as 73461 and 63881
expr68004_lin[expr68004_lin < 1] <- 1
expr68004_log2 <- log2(expr68004_lin)
expr68004_log2 <- normalize.quantiles(as.matrix(expr68004_log2))
rownames(expr68004_log2) <- rownames(expr68004_lin)
colnames(expr68004_log2) <- colnames(expr68004_lin)
rm(expr68004_lin); gc()
```

---

## Annotation Mapping Notes

### Key asymmetry between datasets (NOT a bug — must be preserved)
GSE73461/62/63 and GSE68004 expression matrices use Array_Address_Id as row keys.
Map these using annov4$PROBEID (which holds Array_Address_Id values).

GSE63881 expression matrix uses ILMN_ probe IDs as row keys.
Map these using annov4$ID (which holds ILMN_ values).

### GPL file column names expected by pipeline.r
```r
# annov4 columns after renaming (from GPL10558):
# PROBEID = Array_Address_Id (numeric)
# ID      = ILMN_ probe ID
# ENTREZID = Entrez_Gene_ID
# SYMBOL  = Symbol

# annov3 columns after renaming (from GPL6947):
# PROBEID = Array_Address_Id
# ENTREZID = Entrez_Gene_ID
# SYMBOL  = Symbol
```

---

## Planned Improvements (in priority order)

### Must fix before any journal submission
- [x] BUG 0: Extract platform files to expected paths
- [x] BUG 1: Delete premature res61$metadata assignment
- [x] BUG 2: Replace master_expr_filtered with master_expr_combat (was already absent)
- [x] BUG 3: Add set.seed(42) before ConsensusClusterPlus
- [x] BUG 4: Remove filter_by_metadata output overwrite
- [x] BUG 5: Delete dead probe_ids lines
- [x] Fix normalization chain (floor+log2+QN for 61 and 81; remove pooled QN at line ~345)
- [x] BUG 6: Switch volcano plots to adj.P.Val
- [x] BUG 7: Fix GSEA deduplication by max |t|

### Strengthens paper significantly
- [ ] Add GSE68004 as fifth dataset
- [ ] Contact Jaggi/Mejias for per-sample IVIG/CAL annotation (PMID 29813106)
- [ ] Add power calculation for phenotype tests using pwr package
- [ ] Add per-cluster n labels to phenotype bar charts (Figures 12-14)
- [ ] Add cluster stability sensitivity analysis (re-run at 3k/5k/8k MAD genes, report adjusted Rand index)
- [ ] Delete gene_collapse.R (dead code with wrong collapse key — uses SYMBOL not ENTREZID)

### Adds analytical impact
- [ ] CIBERSORTx or xCell immune cell deconvolution
- [ ] WGCNA co-expression modules correlated with cluster membership
- [ ] Cross-validate cluster marker genes against GSE64486 (coronary artery tissue)

---

## What NOT to Change

These are correct and must not be modified:
- The ComBat model matrix: `mod <- model.matrix(~as.factor(Diagnosis), data = master_metadata_final)`
- The k=4 cluster selection (justified by silhouette, PAC, delta area convergence)
- The one-vs-rest limma contrast design for multi-class DE
- pItem=0.8, reps=1000 for ConsensusClusterPlus (standard values)
- MaxMean method in collapseRows (correct for microarray multi-probe collapse)
- FDR threshold of 0.05 with Benjamini-Hochberg

---

## Key Package Roles

- `preprocessCore` — within-study normalize.quantiles()
- `sva` — ComBat batch correction
- `ConsensusClusterPlus` — consensus clustering (REQUIRES set.seed before call)
- `limma` — lmFit > contrasts.fit > eBayes > topTable
- `fgsea` — GSEA on ranked gene lists from limma t-statistics
- `GSVA` — sample-level pathway scoring on full expression matrix
- `WGCNA::collapseRows` — probe-to-gene collapse using MaxMean
- `org.Hs.eg.db` — Entrez ID to gene symbol for GSEA input

---

## R Environment

- R version 4.x
- Machine: MacBook Pro M1 Pro 16GB RAM
- RStudio Version 2025.09.0+387
- Data files are local (pipeline assumes files already downloaded)
- GitHub: https://github.com/splayfern3/kd_transcriptomics

---

## Claude Code Permissions

Claude Code may freely edit:
- `pipeline.r` (main pipeline at repo ROOT)
- `r_code/` — all .R source files
- `CLAUDE.md` (update checklists as bugs are fixed)

Do NOT edit:
- `.git/` directory
- `transcriptome_data/` contents (data files, not code)

---

## Testing Protocol After Each Fix

After every change run:
```r
source("r_code/illumina_loader.R")
source("r_code/meta_filter.R")
source("r_code/process_study.R")
source("r_code/clean_dx.R")
message("All source files loaded successfully")
```

If any fail, stop and report the error. Do NOT commit broken code.

---

## Git Workflow

After each fix:
```bash
git add -A
git commit -m "Fix: [describe what was fixed]"
git log --oneline -3
```

If something breaks:
```bash
git reset --hard HEAD~1
git status
```

---

## Questions Claude Code Should Answer Before Starting

1. `git branch` — confirm on fix/preprocessing-and-bugs, NOT main
2. `ls transcriptome_data/platforms/v4/` — check if GPL10558.txt exists
3. `ls transcriptome_data/platforms/v3/` — check if GPL6947.txt exists
4. If either platform file is missing, fix BUG 0 first before touching anything else
