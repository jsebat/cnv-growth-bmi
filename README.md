# cnv-growth-bmi

Analyses of body mass, growth (height, head circumference) and social-communication / cognitive phenotypes in carriers of
16p11.2 and 22q11.2 copy-number variants (G2MH, SVIP, Searchlight, Cardiff) versus unaffected SSC relatives and SSC probands.
Jonathan Sebat (UCSD), analyses carried out with Claude, October 2026. Companion to the Sacks et al. 2026 biobank study.

**Start here:** `HANDOFF_2026-10-26.md` (pipeline, run order, caveats, Expanse resources, next steps) and
`manuscript_materials_2026-10-24/Analysis_report_2026-10-26.docx` (full methods, results, figure legends, table index, paths, PMID references).

## Layout
- `00_*` to `31_*.R`: analysis scripts, numbered in run order (longitudinal chain 01 -> 05 -> 07 -> 11 -> 29; cross-sectional chain 21/23 -> 25 -> 26 -> 27 -> 28 -> 30; 31 builds the report).
- `manuscript_materials_2026-10-24/`: frozen main figure, supplementary figures, Tables S1-S6 (`supp_tables_bmi_effects.xlsx` with a notes sheet), legends, QC rules, reports.
- Top-level `fig_*.png/.pdf`, `*.csv`, `*.md`: live outputs and write-ups (`unaffected_controls_report.md` is the running analysis log).
- `Molly/`, `aou/`, `ukbb_bodysize/`: side analyses (head-circumference figures, All of Us and UK Biobank anchors, UKBB body size at 10 vs adult).
- `resources/`: copy of the lab GWAS summary-statistics sheet.

## What is not here
Participant-level data stay in the lab Dropbox / Expanse (paths in the handoff and the report): the growth genomatrix, phenotype
releases, `analysis_table.rds`, `hc_visits.rds`, PGS score files, the flagged-record list (`implausible_height_weight*.csv`, which
carries sample IDs) and `tbx1_subjects_zscores*`. Figures show de-identified points; tables are aggregate statistics.

## Requirements
R 4.6 with data.table, mgcv, lme4, lmerTest, childsds, ggplot2, patchwork, openxlsx, officer; plink2 and Slurm on Expanse for the
genotype scoring (`30a_build_ssc_pgs.R` reads the scores back from the data folder). Paths at the top of each script point at the
lab Dropbox; adjust `HERE`, `DB` and `DATA` to run elsewhere.
