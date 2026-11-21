# BMI modifier paper v2 (2026-11-20): shared brief for the analysis work packages

Jon Sebat asked for a major revision of the Sacks "BMI modifier" paper into a v2 draft with this outline:

1. CNV x BMI effects on height
2. CNV x BMI effects on trajectories of head circumference and autistic traits
3. Three-way interactions CNV x BMI x polygenic score -> autism
4. Two-way and three-way interactions CNV x BMI x sex -> autism (exploratory)

The coordinating session writes the Google Doc. Work packages produce analyses, result tables and figure blocks (PNG).
Do not edit Google Docs, do not commit or push to git, do not message Jon. Report back to the coordinator.

## Where things are

- **Code (working copy, HERE):** `C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories`
  Numbered R scripts (00_*.R helpers, 17_main_with_ukbb.R main growth figure, 32_supp_asd_by_cnv.R autism by CNV,
  51b/60/61/62/64/66 SRS-SCQ interaction scripts, 03/05/07/11 longitudinal). Read the header comments of a script before using it.
- **Participant-level data (controlled, never copy elsewhere, never into git):** `C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data`
  - `analysis_table_chop_2026-11-19.rds`: all clinical cohorts, all visits, including CHOP 22q11.2 A-D records (built by 01_build.R + 00_add_chop.R; the 279 G2MH 22q carriers of site 203 / penn_chop are already removed as presumed duplicates of CHOP).
  - `hc_visits.rds`: Simons VIP and Cardiff head-circumference visits including infancy.
  - `combined_srs_scq.rds`: SRS / SCQ sample (G2MH, Simons VIP, Searchlight, SSC) with polygenic scores.
- **Current run and its notes:** `HERE/run_chop_fillgaps_2026-11-19/README_selection_criteria.md` (visit rule `fill_gaps`, CHOP overlap rule), outputs in `expanse_v9e/`, `asd_v9b_random/`.
- **New outputs for this pass:** `HERE/run_v2_2026-11-20/<package>/` (your own subfolder; never write into older run folders).
- **Local R:** `"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"`. Write scripts with the Write tool (bash heredocs mangle `\n`, `\t` and non-ASCII).
- **Expanse (only if your package says so):** `wsl ssh -o BatchMode=yes expanse '<cmd>'` (ignore the Lmod gcc message). Do not let WSL expand `$(...)` locally: write a script file, `scp` it, run it. Jobs: `sbatch` in `ind-shared`, account `ddp195`, template `/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v9e_chop_genotype_bands/run_v9e.sbatch` (Rscript `/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript`, `R_LIBS_USER=$HOME/R/dnmt3a-r`). Small jobs (few cores, < 1 h) need no approval. Every job variant gets its own new output folder. Never cancel a running job. UK Biobank individual data never leave Expanse; copy back only figures and aggregate tables.

## Jon's analysis rules (binding)

- **One analysis per genotype** (16p11.2 del, 16p11.2 dup, 22q11.2 del, 22q11.2 dup) on all data, then a pooled "all four CNVs" estimate adjusted for CNV, then a **heterogeneity test** across genotypes. No discovery / replication split. No "failure" framing.
- Figure scope: 16p11.2 BP4-BP5 and 22q11.2 A-D only.
- One visit per person in cross-sectional analyses, chosen by `00_select_visit.R` rule `fill_gaps` (as in the current run). Longitudinal analyses use all cleaned visits with a child random intercept.
- Scoring through `00_scoring.R` (children: WHO height and BMI, UK90 head circumference; adults: Box-Cox relative to the non-carrier adults of the same cohort). If you use another scorer, report the comparison.
- Fix errors; do not swap in alternative methods for an existing analysis. New analyses are additions and are labelled as new.
- Every script: header comment saying what it does, inputs, outputs, date, and a line `# null model: ...` (or `# null model: none; ...` with the test used).
- Say "trait", not "risk". Label which cohort / device / data source every number comes from.
- Report sample size per genotype and per cohort for every analysis. A cell with zero or near-zero rows means a join error; check it.
- Sensitivity analyses that matter go in the same output folder as CSVs, each named for what it varies.
- Autism outcomes: **SRS and SCQ scores are primary; autism diagnosis is secondary** (diagnosis also brings in CHOP).
- Sex and polygenic-score analyses of autism are exploratory; say so.

## Figure blocks

Final figures are composed by the coordinator from PNG blocks, so keep one style:
ggplot2 `theme_classic(base_size = 15)`, bold lowercase panel tags (size 22), BMI tails red `#F0645A` (>= 70th centile) and teal `#00A5AD` (<= 30th), carriers `#4B6EA9`, non-carriers `#8C9A9A`, females `#C2528B`, males `#3C7FB1`; P values printed as "p = ...", bold when < 0.05; n shown in each panel. Save PNG at 200 dpi, white background, plus the PDF. Full figure width is 20 in; a half-width block is 10 in. Also save a CSV with every plotted number and test next to each PNG.

## What to report back (keep it under 60 lines)

1. Files written (scripts, CSVs, PNGs), with full paths.
2. The key numbers, each with n and its CSV file.
3. Sensitivity analyses run and whether each changed the conclusion.
4. Anything wrong, surprising or unverified, stated plainly. Do not smooth over a result that did not hold.
