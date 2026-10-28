# cnv-growth-bmi

Code and summary results for Figures 3, 4 and 6 of Sacks et al. (2026), "Context-dependent effects of CNVs on growth"
(companion to the biobank meta-analysis code at https://github.com/mollysacks/Sacks2026).
Jonathan Sebat (UCSD), analyses carried out with Claude, July to October 2026.

## Part 1. Figures 1-2: CNV main effects on height and BMI (association tests)

The biobank association tests and the meta-analysis are in Molly Sacks's repository, https://github.com/mollysacks/Sacks2026
(`SharedAnalysisCode/` for height and BMI, `MetaAnalysisCode/`). The recurrent-CNV genotyping tutorial is also in this
repository under `recurrentCNVgenotyping/`.

## Part 2. Figure 3: mirror and asymmetric CNV effects (oligogenic model)

`fig3_mirror_asymmetric/`: the perpendicular offset of each reciprocal CNV from the mirror line and its bimodality
(measurement-error-aware mixture, parametric bootstrap), the size and gene-count comparisons of mirror and asymmetric loci,
and simulations comparing the mirror-only, continuous and single-hit models. The figure is
`perp_offset_bimodality/figures/Figure_height_main_v2_me`; see that folder's `README.md` for the run order.

## Part 3. Figure 4: pathway convergence of asymmetric CNVs and single-gene height associations

`fig4_pathway_analysis/`: gene-set enrichment of the genes inside asymmetric versus mirror recurrent CNVs, of the 207
exome height genes of Kosmicki et al., and of their height-decreasing bias, on one 102-pathway analysis set with a
locus-block permutation for the CNV layer. Inputs (gene tables and GMT libraries), scripts in run order, the master
results table, Figure 4 and the full report are all in that folder; see its `README.md`.

## Part 4. Figure 6 and supplement: body mass as a modifier of growth and autistic traits in CNV carriers

Analyses of body mass, growth (height, head circumference) and social-communication / cognitive phenotypes in carriers of
16p11.2 and 22q11.2 copy-number variants (G2MH, SVIP, Searchlight, Cardiff) versus unaffected SSC relatives and SSC probands.

**Start here:** `HANDOFF_2026-10-26.md` (pipeline, run order, caveats, Expanse resources, next steps) and
`manuscript_materials_2026-10-24/Analysis_report_2026-10-26.docx` (full methods, results, figure legends, table index, paths, PMID references).

### Layout
- `00_*` to `33_*.R`: analysis scripts, numbered in run order (longitudinal chain 01 -> 05 -> 07 -> 11 -> 29; cross-sectional chain 21/23 -> 25 -> 26 -> 27 -> 28 -> 30; 31 builds the report).
- `manuscript_materials_2026-10-24/`: frozen main figure, supplementary figures (incl. `fig_supp_asd_by_cnv`: height, head circumference and BMI by autism diagnosis within CNV group, with a PGS sensitivity analysis), Tables S1-S8 (`supp_tables_bmi_effects.xlsx` with a notes sheet), legends, QC rules, reports.
- Top-level `fig_*.png/.pdf`, `*.csv`, `*.md`: live outputs and write-ups (`unaffected_controls_report.md` is the running analysis log).
- `Molly/`, `aou/`, `ukbb_bodysize/`: side analyses (head-circumference figures, All of Us and UK Biobank anchors, UKBB body size at 10 vs adult).
- `resources/`: copy of the lab GWAS summary-statistics sheet.

### What is not here
Participant-level data stay in the lab Dropbox / Expanse (paths in the handoff and the report): the growth genomatrix, phenotype
releases, `analysis_table.rds`, `hc_visits.rds`, PGS score files, the flagged-record list (`implausible_height_weight*.csv`, which
carries sample IDs) and `tbx1_subjects_zscores*`. Figures show de-identified points; tables are aggregate statistics.

### Requirements
R 4.6 with data.table, mgcv, lme4, lmerTest, childsds, ggplot2, patchwork, openxlsx, officer; plink2 and Slurm on Expanse for the
genotype scoring (`30a_build_ssc_pgs.R` reads the scores back from the data folder). Paths at the top of each script point at the
lab Dropbox; adjust `HERE`, `DB` and `DATA` to run elsewhere.
