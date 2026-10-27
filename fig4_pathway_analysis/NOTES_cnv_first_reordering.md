# CNV-first three-way convergence workflow

This workflow uses the **uncorrected 22q11.2 classification**, exactly as requested.

## Goal
Reverse the original Kosmicki-first analysis order while leaving the underlying statistical tests unchanged:

1. Discover pathways enriched in **asymmetric CNVs**.
2. Among those pathways, test enrichment in the **207 Kosmicki height genes**.
3. Among those, test **negative-effect bias** in Kosmicki genes.

With the uncorrected 22q11.2 classification, **GO:0008284 Positive Regulation Of Cell Population Proliferation** has asymmetric-CNV P=0.031, so a standard nominal discovery threshold of **P<0.05 and OR>1** retains it.

## Fastest reproducible route
Use `kosmicki_pathway_cnv_master_results.csv`. It already contains the exact pathway-level statistics from the uncorrected analysis across 5,712 tested pathways. Reordering the analysis does **not** require recomputing the enrichment tests; it only changes the sequence of filtering and multiple-testing families.

Run:

```bash
Rscript run_CNV_first_three_way_analysis.R \
  kosmicki_pathway_cnv_master_results.csv \
  CNV_first_results
```

Required R packages:

```r
install.packages(c("data.table", "ggplot2", "patchwork", "ggrepel"))
```

## Sequential tests

### Step 1 — Asymmetric-CNV discovery
Retain pathways with:
- `asymmetric_cnv_odds_ratio > 1`
- `asymmetric_cnv_p < 0.05`

### Step 2 — Kosmicki enrichment
Within Step-1 pathways, retain:
- `kosmicki_odds_ratio > 1`
- `kosmicki_p < 0.05`

Also calculate BH FDR across only the Step-1 pathways (`kosmicki_q_among_step1`).

### Step 3 — Negative-effect bias
Within Step-2 pathways, retain:
- `direction_odds_ratio > 1`
- `direction_fisher_p_one_sided < 0.05`

Also calculate BH FDR across only the Step-2 pathways (`direction_q_among_step2`).

## Expected GO:0008284 values from the uncorrected analysis
- Asymmetric-CNV enrichment: OR 7.04, P 0.031
- Kosmicki enrichment: 17/176 mapped genes, OR 3.11, P 1.05e-4
- Negative-effect bias: 16 negative / 1 positive, OR 8.70, P 0.009
- Locus-block permutation in the original analysis: P 0.033

## Main figure
The script creates a three-panel sequential discovery figure:
- **A:** all pathways tested for asymmetric-CNV enrichment, highlighting GO:0008284
- **B:** Step-1 pathways tested for Kosmicki enrichment
- **C:** Step-2 pathways tested for negative-effect bias

The title line reports the funnel counts at each stage. This makes the logic of the reversed analysis visible in one figure.

## Raw inputs for rebuilding from scratch
The following are included for audit/reconstruction:
- `kosmicki_207_genes_annotated.csv`
- `cnv_genes_annotated.csv` — this is the **uncorrected** CNV classification used in the original analysis
- `kosmicki_pathway_cnv_pipeline.R` — original Kosmicki-first pipeline

The original pipeline additionally expects:
- `GO_Biological_Process_2023.gmt`
- `KEGG_2021_Human.gmt`

Those exact GMT files were part of the prior reproducibility bundle, but they are not separately exposed in the current file surface. They are **not needed** for the CNV-first reordering if you start from `kosmicki_pathway_cnv_master_results.csv`.

## Recommended Claude Code prompt

> Run the CNV-first three-way convergence analysis using `kosmicki_pathway_cnv_master_results.csv`. Do not alter the 22q11.2 classifications. First identify pathways enriched in asymmetric CNVs using one-sided P<0.05 and OR>1. Then test only those pathways for Kosmicki enrichment using P<0.05 and OR>1, with BH FDR across that selected family. Then test only the survivors for negative-effect bias in Kosmicki genes using the existing one-sided Fisher test, again reporting BH FDR across that family. Produce the three-panel main figure specified in `run_CNV_first_three_way_analysis.R`, a compact results table, and an audit line for GO:0008284. Do not recompute or relabel the underlying 22q11.2 CNV classes.
