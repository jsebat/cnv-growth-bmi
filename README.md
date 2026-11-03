# cnv-growth-bmi

Code and summary results for Figures 1F-G, 2, 3 and 5 and Supplementary Fig. 5 of Sacks et al. (2026), "Context-dependent effects of CNVs on growth"
(companion to the biobank meta-analysis code at https://github.com/mollysacks/Sacks2026).
Jonathan Sebat (UCSD), analyses carried out with Claude, July to November 2026. Figure numbers follow the November 2026 renumbering (old Fig. 0 is now Fig. 1A; old Figs 3, 4, 6 are now Figs 2, 3, 5).

## Part 1. Figure 1: CNV main effects on height and BMI (association tests)

The biobank association tests and the meta-analysis are in Molly Sacks's repository, https://github.com/mollysacks/Sacks2026
(`SharedAnalysisCode/` for height and BMI, `MetaAnalysisCode/`). The recurrent-CNV genotyping tutorial is also in this
repository under `recurrentCNVgenotyping/`.

`fig1_main_effects/`: Figure 1 with the federated-analysis schematic as panel A, and `46_fig1FG_dose_response_whiskers.R`, which
redraws the duplication-versus-deletion panels (F, G) from Supplementary Table 3 with ±1 s.e. whiskers on both effects (the original
plotting code was not recovered).

## Part 2. Figure 2: mirror and asymmetric CNV effects (oligogenic model)

`fig3_mirror_asymmetric/`: the perpendicular offset of each reciprocal CNV from the mirror line and its bimodality
(measurement-error-aware mixture, parametric bootstrap), the size and gene-count comparisons of mirror and asymmetric loci,
and simulations comparing the mirror-only, continuous and single-hit models. The current figure is
`make_height_figure_v4.py` -> `Figure_height_main_v4_ad` (panels A-B bimodality, C-D gene count, E-F TWAS concordance,
G-J size and effect differences, K-O models); `perp_offset_bimodality/` holds the analysis code, outputs and earlier figure
versions (see its `README.md` for the run order), `CNV_height_BMI_corrected_directional_models_v2/` the model fits and
simulations, `CNV_burden_behavior_figure/` the model-behaviour supplementary figure, and `PGC_asym_model_diff/` the comparison
with Molly Sacks's model implementation.

### Common-variant TWAS concordance (Figure 2E-F, Supplementary Note 4, Supplementary Table 29)

`concordance/`: tests whether the genes inside each CNV predict its effects. Gene-level height and BMI TWAS z scores (Mancuso et
al. 2017, mean over expression panels) are summed over the contiguous run of genes in each CNV and compared with the CNV's dosage
effect, (DUP - DEL)/2, which cancels a shift shared by both alleles and halves one on the duplication axis. Significance comes from
circular genomic permutation (10,000 rotations, locus clusters moved as units), with a null of random windows matched for gene
count; the main analysis uses non-overlapping CNVs (spanning intervals kept, 22q11.2 A-D). Also: rare-variant concordance with
exome height and BMI genes (Kosmicki, Akbari), per-CNV sign tests and FDR, leave-one-out, separate-lines and asymmetric-offset
tests, and the Supplementary Table 29 builder (`build_supp_table_s29.py`). `verify_numbers.py` checks every effect size against
Supplementary Table 3 and the Figure 2 inputs. See `concordance/README.md`. The raw TWAS input files (public, Mancuso et al.)
are not included because of their size.

## Part 3. Figure 3: pathway convergence of asymmetric CNVs and single-gene height associations

`fig4_pathway_analysis/`: gene-set enrichment of the genes inside asymmetric versus mirror recurrent CNVs, of the 207
exome height genes of Kosmicki et al., and of their height-decreasing bias, on one 102-pathway analysis set with a
locus-block permutation for the CNV layer. Inputs (gene tables and GMT libraries), scripts in run order, the master
results table, the figure and the full report are all in that folder; see its `README.md`.

## Part 4. Figure 5 and supplement: body mass as a modifier of growth and autistic traits in CNV carriers

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

### Extended two-component models (Supplementary Note 5, Supplementary Table 30)

`fig3_mirror_asymmetric/CNV_height_BMI_corrected_directional_models_v2/code/models_grid.py` extends the Fig. 2K-O models in a 2 x 2
design (along the mirror line / perpendicular; diffuse / sparse hits) with summed TWAS z on the mirror component: 21 architectures per
trait (`models_twas_burden.py`, `models_mirror_hits.py` hold the step-by-step fits; outputs in `py_outputs/`). Height: mirror + TWAS +
sparse asymmetric hits is best; no mirror hits. BMI: mirror + continuous upper-right shift with constrained genes, weakly supported.
Set TWAS_SOURCE=brain to use the PsychENCODE cortex TWAS (`concordance/twas_psychencode.py`, gene tables in
`concordance/psychencode_twas/`).

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
