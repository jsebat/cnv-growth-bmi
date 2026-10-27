# Figure 4: pathway convergence of asymmetric CNVs and single-gene height associations

Pathway analysis behind Figure 4 of Sacks et al. (2026), "Context-dependent effects of CNVs on growth".
Jonathan Sebat (UCSD), analyses carried out with Claude, July to October 2026.

## Question

Recurrent CNVs fall into two dose-response classes for height (Figure 3): "mirror" loci, where deletions and
duplications shift height in opposite directions, and "asymmetric" loci, where the duplication effect is buffered or
reversed. Do the genes inside asymmetric CNVs share biology with the single genes that reduce height in exome
association studies?

## Inputs (`inputs/`)

| File | Contents |
|---|---|
| `cnv_genes_annotated.csv` | 942 gene-by-interval rows (708 unique genes) for the 34 recurrent CNV intervals with reciprocal height effects, each labelled with its interval's dose-response class (13 asymmetric, 21 mirror intervals; the nested 22q11.2 intervals are entered separately, so 13 genes belong to both classes). GRCh38 coordinates. |
| `kosmicki_207_genes_annotated.csv` | The 207 exome-wide significant height genes of Kosmicki et al. (medRxiv 2026, doi 10.64898/2026.06.22.26355163) with effect size (s.d. per allele) and direction. |
| `GO_Biological_Process_2023.gmt` | GO Biological Process 2023 gene-set library (Enrichr release). |
| `KEGG_2021_Human.gmt` | KEGG 2021 Human gene-set library (Enrichr release). |

## Scripts, in run order

1. `kosmicki_pathway_cnv_pipeline.R` computes every pathway statistic once and writes the master table
   `kosmicki_pathway_cnv_master_results.csv` (5,712 gene sets of 5 to 500 genes). Per gene set: enrichment of the
   207 Kosmicki genes (one-sided hypergeometric against the library universe), the height-decreasing bias of the
   set's Kosmicki genes relative to all mapped Kosmicki genes (one-sided hypergeometric), enrichment in asymmetric
   and in mirror CNV genes (one-sided hypergeometric within the CNV gene universe), and a locus-block permutation
   for the asymmetric contrast (5,000 permutations that reassign the asymmetric label among the 34 intervals within
   gene-count tertiles, seed 20260719).

   ```bash
   Rscript kosmicki_pathway_cnv_pipeline.R inputs pipeline_results
   ```

   The master table shipped here is the one the figure and tables were made from; rerunning regenerates it.

2. `run_three_way_convergence.R` reads the master table, defines the 102-pathway analysis set (gene sets with at
   least four genes in the CNV universe, the minimum at which the mirror test can reach P < 0.05; the asymmetric set,
   at least six genes, is a subset), applies Bonferroni (0.05/102) and BH FDR within that set, and reports the
   pathways significant in all three tests. Writes `results/analysis_set_102_pathways.csv`,
   `results/00_analysis_log.txt` and the `test*_significant.csv` tables.

   ```bash
   Rscript run_three_way_convergence.R kosmicki_pathway_cnv_master_results.csv results
   ```

3. `make_main_figure.R` draws Figure 4 (`results/Figure_4_convergence.{png,pdf}`). Panel A is the asymmetric-versus-
   mirror contrast with a two-sided Fisher exact P computed from the stored 2 x 2 counts
   (`results/panelA_asym_vs_mirror_two_sided.csv`); that is the only statistic computed outside the pipeline.
   `make_main_figure_v1_6panel.R` is the earlier six-panel layout.

   ```bash
   Rscript make_main_figure.R kosmicki_pathway_cnv_master_results.csv inputs/kosmicki_207_genes_annotated.csv inputs/cnv_genes_annotated.csv results
   ```

4. `make_intersection_table.R` lists the pathways at each intersection of the three tests and checks whether any
   gene is both a Kosmicki height gene and a CNV gene of the relevant class (none is). It rebuilds gene memberships
   from the GMT files, so point it at `inputs/`:

   ```bash
   Rscript make_intersection_table.R kosmicki_pathway_cnv_master_results.csv inputs results
   ```

## Result

GO:0008284, positive regulation of cell population proliferation, is the only gene set nominal in all three tests:
asymmetric-CNV enrichment (9 of its 10 CNV-universe genes asymmetric, OR 7.0, two-sided P = 0.047, locus-block
permutation P = 0.033), Kosmicki height-gene enrichment (17 of 176 mapped genes, OR 3.1, P = 1.1 x 10^-4, Bonferroni
significant) and height-decreasing bias (16 of 17 decreasing, OR 8.7, P = 0.009). The same analysis with mirror CNVs
yields no convergent pathway. `FINAL_REPORT.pdf` is the full write-up; `results/README.md` and
`results/Figure_4_legend.txt` give the headline table and the legend; `NOTES_cnv_first_reordering.md` describes the
CNV-first ordering of the three tests, which gives the same answer.

## Requirements

R 4.6 with data.table, ggplot2, patchwork and ggrepel. Gene-level inputs only; no participant data.
