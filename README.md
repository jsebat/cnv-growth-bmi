# cnv-growth-bmi

Code and summary results for the Sacks et al. (2026) paper on body mass as a modifier of growth and autistic traits in carriers
of 16p11.2 and 22q11.2 copy-number variants (G2MH, SVIP, Searchlight, Cardiff, CHU Sainte-Justine) versus unaffected SSC
relatives and SSC probands. Jonathan Sebat (UCSD), analyses carried out with Claude, July to November 2026.

The companion paper, "Gene dosage architecture of anthropometric traits" (mirror, constraint and single-hit components of CNV
effects on height and BMI, CNV GWAS in UK Biobank and All of Us), has its own repository,
https://github.com/jsebat/cnv-dosage-architecture. The two were split on 2026-11-19; the folders `concordance/`,
`fig1_main_effects/`, `fig3_mirror_asymmetric/`, `fig4_pathway_analysis/` and `recurrentCNVgenotyping/` moved there, as did the
`cnv-gwas` branch. The biobank association tests are in Molly Sacks's repository, https://github.com/mollysacks/Sacks2026.

**Start here:** `HANDOFF_2026-11-20_v2.md` (v2 draft: paper restructured around CNV x BMI effects on height, head growth, autistic traits, polygenic score and sex; the control set was corrected on 2026-11-20, see that file; all v2 results are in `run_v2_2026-11-20/`, scripts `80`-`89`). Earlier state: `HANDOFF_2026-11-14.md` (current state: Box-Cox scaling of SRS/SCQ and adult BMI, combined cohorts, SVIP
replication, three-way and PGS analyses, data hub), then `HANDOFF_2026-10-26.md` (pipeline, run order, caveats, Expanse
resources) and `manuscript_materials_2026-10-24/Analysis_report_2026-10-26.docx` (methods, results, legends, table index).

### Layout
- `00_*` to `45_*.R`: analysis scripts, numbered in run order (longitudinal chain 01 -> 05 -> 07 -> 11 -> 29; cross-sectional chain 21/23 -> 25 -> 26 -> 27 -> 28 -> 30; 31 builds the report).
- `46_*` to `70_*.R`, `35b`-`35f`, `51`/`51b`, `00_boxcox_srs_scq.R`: the November 2026 round (Box-Cox SRS/SCQ, main figures v6-v8, interaction forests, SVIP replication, combined-cohort and three-way models, sex and PGS checks, data-hub manifest).
- `manuscript_materials_2026-10-24/`: frozen pre-split package (main figure, supplementary figures, Tables S1-S8 and the shared supplementary tables, legends, QC rules, reports). Its T29/T30, S23 and S25-S26 tables belong to the dosage architecture paper.
- Top-level `fig_*.png/.pdf`, `*.csv`, `*.md`: live outputs and write-ups (`unaffected_controls_report.md` is the running analysis log).
- `Molly/`, `aou/`, `ukbb_bodysize/`: side analyses (head-circumference figures, All of Us and UK Biobank anchors, UKBB body size at 10 vs adult).
- `resources/`: copy of the lab GWAS summary-statistics sheet.

### Pubertal timing and Mendelian randomization (Supplementary Fig. 5, Supplementary Tables 24, 27, 28)

Reviewer request: is body mass causal? `fig_supp5_combined.R` (runs on Expanse) builds Supplementary Fig. 5 from UK Biobank: adult
height by recalled body size at age 10 and by adult BMI group (a), age at menarche by genotype (b), adult height per year of
menarche (c), polygenic-score instruments (d) and the Mendelian randomization of age at menarche on adult height (e).
`ukbb_mr_childhood_obesity.R` (childhood obesity, childhood BMI and adult BMI scores as instruments) and
`ukbb_mr_puberty_pathway.R` (adiposity scores on age at menarche, two-step mediation, multivariable model) give
Supplementary Table 27. In the clinical cohorts, `mr_reduced_form.R`, `mr_mediation.R` and `mr_first_stage.R` (Supplementary
Table 28) and `45_menarche_pgs_mr.R` (menarche score) test the same scores in G2MH carriers and SSC controls.
`44_supp_tables_27_29.R` writes Supplementary Tables 24 and 27-29. Caveat: the menarche (Guan et al. 2025), adult BMI and height
GWAS include UK Biobank participants, so those instruments are not independent of the UK Biobank sample.

### What is not here
Participant-level data stay in the lab Dropbox / Expanse (paths in the handoff and the report): the growth genomatrix, phenotype
releases, `analysis_table.rds`, `hc_visits.rds`, PGS score files, the flagged-record list (`implausible_height_weight*.csv`, which
carries sample IDs) and `tbx1_subjects_zscores*`. Figures show de-identified points; tables are aggregate statistics.

### Requirements
R 4.6 with data.table, mgcv, lme4, lmerTest, childsds, ggplot2, patchwork, openxlsx, officer; plink2 and Slurm on Expanse for the
genotype scoring (`30a_build_ssc_pgs.R` reads the scores back from the data folder). Paths at the top of each script point at the
lab Dropbox; adjust `HERE`, `DB` and `DATA` to run elsewhere.
