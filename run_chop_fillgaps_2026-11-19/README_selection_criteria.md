# CHOP 22q11.2 records and age-gap-filling visit selection (run 2026-11-19)

Jon's request (2026-11-19): add the new samples to the height and BMI analysis. The 22q11.2 deletion group has almost no children
under 7, and most groups have a gap between 18 and 40. Fill these gaps to get a full trajectory, choosing the youngest visit
rather than the oldest for some people where the longitudinal data allow. Do the same for Simons VIP and the other cohorts with
repeated visits. Record the selection criteria.

## Samples added
- **CHOP 22q11.2 growth records** (T. Blaine Crowley, received 2026-11-19), loaded by `00_add_chop.R`.
  - Only standard **A–D** carriers: 696 deletions and 71 duplications, 24,501 encounters. The figure covers 16p11.2 BP4–BP5 and 22q11.2 A–D only.
  - Height and weight are paired when recorded at the same age (same encounter); BMI = weight / height².
  - Head circumference is joined at the same encounter.
  - ASD diagnosis comes from the CHOP index (Yes/No; the sheet does not say whether "No" means evaluated).
- **Overlap.** CHOP is also a G2MH recruitment site (`penn_chop`, MS ID prefix 203), and no CHOP–G2MH identifier map exists yet.
  - Rule `drop_g2mh_chop_site` (default): the 279 G2MH 22q11.2 carriers from that site are removed from the growth tables, because their growth is expected to be in the CHOP records.
  - Non-carrier controls are untouched.
  - `CHOP_OVERLAP=keep` gives the sensitivity analysis.

## Visit selection (`00_select_visit.R`, rule `fill_gaps`)
1. **Eligible visits:** outcome (height or head circumference) and BMI both measured, age 2–73.
2. **Age bands:** 2–7, 7–13, 13–19, 19–40, 40–73 years.
3. **Visits in one band only:** the person keeps the latest visit in that band. This is the previous rule.
4. **Visits in several bands:** the person goes to the band with the fewest people of the same genotype so far.
   - People with fewer candidate bands are placed first, then by ID, so the result is deterministic.
   - Ties go to the oldest band, so the previous choice is kept unless a younger band is emptier.
   - Within the chosen band, one visit is drawn at random (seed 20261119), so ages spread across the band. (First version used the latest visit of the band; that piled ages at the band edges 7, 13 and 19 y, and was replaced on 2026-11-19.)
5. **One record per person,** as before; only which visit represents a person changes.

`VISIT_RULE=latest` restores the previous rule.

People per band (height; `expanse_v9/visit_selection_HEIGHT.csv`). UK Biobank adults 40–73 are added on top of these.

| Genotype | Rule | 2–7 | 7–13 | 13–19 | 19–40 | 40–73 |
|---|---|---|---|---|---|---|
| 16p11.2 del | latest visit | 65 | 160 | 111 | 42 | 15 |
| 16p11.2 del | fill gaps | 106 | 123 | 107 | 42 | 15 |
| 16p11.2 dup | latest visit | 34 | 77 | 44 | 43 | 33 |
| 16p11.2 dup | fill gaps | 45 | 67 | 45 | 41 | 33 |
| 22q11.2 del | latest visit | 79 | 251 | 269 | 334 | 32 |
| 22q11.2 del | fill gaps | 233 | 233 | 233 | 234 | 32 |
| 22q11.2 dup | latest visit | 10 | 56 | 57 | 40 | 17 |
| 22q11.2 dup | fill gaps | 41 | 41 | 41 | 40 | 17 |

Both rows use the new build, which includes CHOP. The 19–40 gap for 16p11.2 cannot be filled from the data we have: none of its carriers with repeated visits were measured in that range.

## Runs
**Current: v9e** (Expanse `v9e_chop_genotype_bands`, jobs 55894696 and 55894742; copies in `expanse_v9e/`).
- Same data and statistics as v9c and v9d. Panel d shows 400 of 684 points.
- Panels a′ and f′ now print, per age band, whether the BMI effect differs by genotype (F test of genotype × BMI); saved in `main_with_ukbb_bmi_genotype_by_band.csv`.
- Supplementary Fig. 6: `fig_supp_main_effects_by_band.png`, from `43_supp_main_effects_by_band.R`.
- Supplementary Tables 17, 18, 21 and 22: `supplementary_tables/` here, and `manuscript_materials_2026-10-24/Sacks Supplementary Tables v3_chop.xlsx` (built by `47_supp_tables_chop_2026-11-19.R`).

Earlier figure, v9c: (Expanse `v9c_chop_smooth`, job 55894514; copies in `expanse_v9c/`). Same data and tests as v9b. Descriptive BMI-tail curves in panels with fewer than 250 people now use k = 5 and gamma = 1.5, which removes the left-edge curl Jon saw in the 16p11.2 del and 22q11.2 dup head-circumference panels.

**v9b: random visit within the chosen band (random visit within the chosen band, stiffer curves for small panels).**
- Expanse folder: `v9b_chop_fillgaps_random`, job 55894440.
- Copies: `expanse_v9b/`; autism-by-CNV outputs in `asd_v9b_random/`.
- In panels with fewer than 150 people, the descriptive BMI-tail curves use k = 5 instead of 8. Tests are not affected.
- The first version (`v9_chop_fillgaps`, `expanse_v9/`, autism outputs in this folder's root) is superseded; it is kept for comparison.

Earlier run (v9):
- **Main figure (height and head-circumference trajectories by BMI, with UK Biobank):** run on Expanse.
  - Folder: `/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v9_chop_fillgaps`
  - Job 55894135; script `35d_main_figure_v9_chop.R` via `run35d_v9.sbatch`
  - Settings: `ANALYSIS_TABLE=analysis_table_chop_2026-11-19.rds`, `VISIT_RULE=fill_gaps`
  - Aggregate outputs copied to `expanse_v9/`; the v8 tables, for comparison, are in `expanse_v8_reference/`.
- **Autism by CNV (supplementary figure and Tables S7–S8):** `32_supp_asd_by_cnv.R` run locally with `ADD_CHOP=TRUE VISIT_RULE=fill_gaps OUT_DIR=<this folder>`.
- **Participant-level table:** `Jon/Claude/CNV_growth_trajectories_data/analysis_table_chop_2026-11-19.rds`, built by `01_build.R`. The previous build, `analysis_table.rds`, is kept.
