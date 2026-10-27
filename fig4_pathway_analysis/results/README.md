# Three-way convergence analysis — results

Generated 2026-09-20, R 4.6.1 (data.table, ggplot2, patchwork, ggrepel). No pathway statistic recomputed.

**The full write-up is `../FINAL_REPORT.pdf`.** This folder holds the figures, tables and logs it cites.

## Headline

Across an analysis set of **102 pathways**, **GO:0008284 — Positive Regulation Of Cell Population
Proliferation** is the only pathway significant in all three tests of the asymmetric analysis. It is
**Bonferroni significant in test C** (P = 1.05e-4 vs 0.05/102 = 4.9e-4; BH q = 0.0045). The same analysis run
with mirror CNVs yields **no convergent pathway**.

| Test | GO:0008284 | P | q (BH) | P<0.05 | Bonferroni |
|---|---|---|---|---|---|
| A — enrichment in mirror CNVs | 3 genes, OR 0.44 | 0.94 | 1.00 | 3 | 0 |
| B — enrichment in asymmetric CNVs | 9/10 universe genes, OR 7.04 | 0.031 | 1.00 | 2 | 0 |
| C — enrichment in Kosmicki height genes | 17/176 genes, OR 3.11 | **1.05e-4** | **0.0045** | 22 | 3 |
| D — negative-effect bias vs all 207 | 16 neg / 1 pos, OR 8.70 | 0.0090 | 0.70 | 3 | 0 |
| **Convergent (asymmetric: B, C, D)** | — | — | — | **1** | — |
| **Convergent (mirror: A, C, D)** | — | — | — | **0** | — |

## Analysis set

5,727 GMT terms → 5,712 with statistics (retained at 5–500 genes) → **102 analysed**.

The set is the **union** of the pathways attainable for the asymmetric CNV test (≥6 genes in the CNV universe,
48 pathways) and those attainable for the mirror test (≥4, 102 pathways). The asymmetric set is a strict subset
of the mirror set, so the union is the mirror set: 102 pathways (68 GO:BP + 34 KEGG).

**This one set is used for all four tests** and fixes the multiple-testing burden throughout:
Bonferroni = 0.05/102 = 4.9e-4, BH FDR within the set. Tests C and D do not involve CNVs, so they are computed
once and shared by the two convergence analyses; only the CNV test differs between them.

## Rank of GO:0008284 within the set

| Test | Rank by P | Rank by OR |
|---|---|---|
| A — mirror-CNV | 88 of 102 | 83 of 102 |
| B — asymmetric-CNV | **1 of 102** | 3 of 102 |
| C — Kosmicki | 2 of 102 | 9 of 102 |
| D — negative-effect bias | **1 of 102** | 2 of 102 |

## Intersections

| Analysis | Intersection | Pathway | Database |
|---|---|---|---|
| Asymmetric | **B\|C\|D** | **Pos Reg of Cell Population Proliferation (GO:0008284)** | GO:BP |
| Asymmetric | B\|C | *none* | — |
| Asymmetric | C\|D | PI3K-Akt signaling pathway | KEGG |
| Asymmetric | C\|D | Pos Reg of MAPK Cascade (GO:0043410) | GO:BP |
| Mirror | **A\|C\|D** | ***none*** | — |
| Mirror | A\|C | Protein Ubiquitination (GO:0016567) | GO:BP |
| Mirror | C\|D | GO:0008284, PI3K-Akt, GO:0043410 | GO:BP, KEGG, GO:BP |

No gene in any of these pathways is both a Kosmicki height gene and a Sacks CNV gene.

## Files

| File | Contents |
|---|---|
| `00_analysis_log.txt` | Full log: set derivation, the four tests, convergence, ranks |
| `Figure_1_three_way_convergence.png` / `.pdf` | **Figure 1**, 6 panels (A–D pathway level, E–F gene level), 400 dpi |
| `analysis_set_102_pathways.csv` | The 102 pathways: all four tests, q-values, Bonferroni and convergence flags |
| `analysis_summary.csv` | Headline counts and thresholds |
| `analysis_set_attainability.csv` | Best-case attainable P by genes-in-universe, both CNV classes |
| `convergent_pathways.csv` | Pathways significant in all three tests of either analysis |
| `convergent_pathway_ranks.csv` | Rank of GO:0008284 by P and by OR in each test |
| `intersection_pathways.csv` | Pathways at each intersection, with gene counts and shared genes |
| `testA_mirror_significant.csv` etc. | Pathways significant in each test |

## Re-run

From the parent `Final/` directory, in order:

```bash
Rscript run_three_way_convergence.R
Rscript make_main_figure.R
Rscript make_intersection_table.R
```
