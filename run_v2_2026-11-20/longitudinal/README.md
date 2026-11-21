# Longitudinal ("next visit") analyses, BMI modifier paper v2 (work package A, 2026-11-20)

> **Update 2026-11-21 (JS): cohort adjustment, final rule.** Cohort and genotype are partly confounded (CHOP has only 22q11.2 carriers,
> Simons VIP and Searchlight only 16p11.2). A model of one genotype carries a cohort term (unchanged). A model that pools the four
> CNVs carries one intercept per cohort within genotype, in place of separate cohort and genotype terms. The same rule is used for
> every pooled model of the paper (autistic traits, autism diagnosis, polygenic score, sex). The three possible adjustments agree:
>
> | Pooled estimate | cohort within genotype (primary) | cohort + genotype (until 2026-11-20) | genotype only |
> |---|---|---|---|
> | Height, 2-8 y | +0.089, P 1.4e-14 | +0.089, P 1.9e-14 | +0.088, P 3.7e-14 |
> | Height, 8-13 y | +0.030, P 0.0047 | +0.029, P 0.0055 | +0.027, P 0.012 |
> | Height, 13-18 y | -0.023, P 0.076 | -0.022, P 0.092 | -0.023, P 0.073 |
> | Height, BMI x window | 4.5e-11 | 8.2e-11 | 6.8e-11 |
> | Head circumference, all ages | +0.029, P 0.010 | +0.030, P 0.0089 | +0.028, P 0.013 |
> | Head circumference, heterogeneity across genotypes | 0.046 | 0.030 | 0.034 |
> | Head circumference, one interval per child | +0.035, P 0.18 | +0.037, P 0.17 | +0.031, P 0.24 |
>
> The CSV files in this folder are current (both alternatives are rows of `sensitivity.csv`; see also `cohort_adjustment_check.csv`,
> script 91). Earlier outputs: `previous_cohort_adjusted_2026-11-20/` and `previous_genotype_only_2026-11-21/`. In the tables below, pooled
> cells still show earlier values where they differ in the third decimal.


Question: in CNV-carrier children, does BMI Z at a visit predict the change in (1) head-circumference Z and (2) height Z by
the next visit? Rebuilt inside this code base so the paper no longer depends on the grant-figure code
(`Jon/Claude/COMET_ASD_ARPAH/analysis` A1, A4, A7, A8).

Scripts (code-base root): `80_longitudinal_v2.R` (analysis, log in `80_run.log`), `81_fig_longitudinal_v2.R` (figure blocks).
Participant-level panels for the figures are in `Jon/Claude/CNV_growth_trajectories_data/longitudinal_v2_panels.rds` (controlled;
nothing participant-level is in this folder).

## Data

| Outcome | Cohorts that contribute | Source |
|---|---|---|
| Height, 2-18 y | Simons VIP, Searchlight, Cardiff (16p11.2 and 22q11.2), CHOP 22q11.2 A-D | `analysis_table_chop_2026-11-19.rds` (raw height and BMI) |
| Head circumference, 0.1-18 y | Simons VIP (v10, every visit from infancy), Cardiff, CHOP 22q11.2 A-D | `hc_visits.rds` (Simons VIP, Cardiff); CHOP workbook through `chop_rows()` of `00_add_chop.R` |

- Every cohort was checked for repeated child visits (`carrier_repeat_visits_by_cohort.csv`). G2MH (290 carrier children), CHU
  Sainte-Justine (37) and SSC carriers (15) have one child visit each and cannot contribute; G2MH repeats exist only in adults.
  Searchlight has repeated height (27 children with 2+ visits) but no repeated head circumference.
- **No non-carrier control has repeated measurements** (`controls_repeat_check.csv`: SSC 8,478, G2MH 1,023, Simons VIP 226
  non-carriers, all with one visit). Every estimate here is therefore a within-carrier estimate; there is no control slope.
- Scoring: `00_scoring.R` for every cohort, recomputed here from raw values (height and BMI: WHO 2006 under 5 y, WHO 2007 from 5 y;
  head circumference: UK90).
- CHOP head circumference is read from the workbook, not from `analysis_table_chop`: rule 2 of `00_qc_growth.R` (height outside
  60-230 cm -> height and BMI set to missing) deletes every valid infant length under 60 cm, which removes the BMI of 378 CHOP
  head-circumference visits before about 4 months (3,133 visits with BMI in the workbook, 2,711 in the table). Simons VIP infancy
  visits never went through that rule. The table-based version is a sensitivity row.

## Cleaning rules (`cleaning_log.csv`)

- Height: `clean_visits()` of `00_helpers.R` (per child, against the last kept visit: under 0.25 y later; identical height before
  16 y; more than 2 cm shorter or more than 20 cm/y after 3 y; Z jump over 2.5 SD). Visits need height and BMI, |Z| <= 5.
  CHOP: 17,819 -> 6,787 visits (10,807 dropped as closer than 0.25 y: dense EHR). Pairs = consecutive kept visits >= 0.5 y apart.
- Head circumference: rule of `07_hc_trajectories.R` (under 0.04 y later before 2 y, 0.25 y after; identical value before 3 y;
  more than 1.5 cm smaller or more than 15 cm/y after 1 y; Z jump over 2.5 SD); 30-70 cm, BMI 10-70, |Z| <= 5. CHOP: 3,451 ->
  2,466 visits. Simons VIP and Cardiff arrive already cleaned in `hc_visits.rds` (no further visit dropped). Pairs = consecutive
  visits with BMI at the first visit, >= 0.1 y apart.

## Models

Per genotype, then all four CNVs adjusted for genotype, then heterogeneity:

`change in Z ~ Z at the visit + interval + sex (+ cohort | + genotype) + ns(age at the visit, 3) + BMI Z at the visit + (1 | child)`
(cohort in the models of one genotype; one intercept per cohort within genotype in the pooled model)

- lmerTest, REML, Satterthwaite P for the BMI coefficient. Heterogeneity = likelihood-ratio test of BMI x genotype (ML fits).
- Age windows (those already used in the paper): the same model with `window + window:BMI`; window = age at the first visit of
  the pair. BMI x window = likelihood-ratio test. Separate-model-per-window fits (the grant code's form) are in
  `results_by_window_separate_models.csv`.
- Sex (exploratory): the same model with `BMI + BMI:sex`.
- Estimates from fewer than 15 children are not reported (NA).
- null model: none; these are the models' own Wald and likelihood-ratio tests.

## Sample sizes (children / visit pairs)

Head circumference (`sample_head_circumference_by_genotype_cohort.csv`):

| Genotype | Simons VIP | Cardiff | CHOP | Total | Median age at first visit of pair (IQR) |
|---|---|---|---|---|---|
| 16p11.2 del | 72 / 274 | 22 / 22 | - | 94 / 296 | 1.1 y (0.5-2.2) |
| 16p11.2 dup | 41 / 158 | 6 / 6 | - | 47 / 164 | 1.2 y (0.5-2.1) |
| 22q11.2 del | - | 50 / 56 | 281 / 1,355 | 331 / 1,411 | 2.2 y (1.0-7.0) |
| 22q11.2 dup | - | 6 / 7 | 40 / 225 | 46 / 232 | 2.0 y (0.8-5.9) |
| All four CNVs | | | | 518 / 2,103 | 1.8 y (0.8-5.7) |

Height (`sample_height_by_genotype_cohort.csv`; by window in `sample_by_genotype_cohort_window.csv`):

| Genotype | Simons VIP | Searchlight | Cardiff | CHOP | Total | Median age |
|---|---|---|---|---|---|---|
| 16p11.2 del | 71 / 172 | 17 / 26 | 29 / 30 | - | 117 / 228 | 8.8 y |
| 16p11.2 dup | 31 / 66 | 5 / 6 | 12 / 12 | - | 48 / 84 | 10.0 y |
| 22q11.2 del | - | - | 73 / 102 | 468 / 2,380 | 541 / 2,482 | 8.8 y |
| 22q11.2 dup | - | - | 11 / 12 | 57 / 301 | 68 / 313 | 8.8 y |
| All four CNVs | | | | | 774 / 3,107 | 8.9 y |

Genotype and cohort are largely confounded: 16p11.2 comes from Simons VIP (+ Searchlight, Cardiff), 22q11.2 mostly from CHOP.

## Results

Head circumference: change in Z to the next visit per SD of BMI Z (`results_head_circumference.csv`)

| Genotype | All ages | P | <2 y | 2-5 y | 5-18 y | BMI x window P |
|---|---|---|---|---|---|---|
| 16p11.2 del | +0.087 (SE 0.037) | 0.023 | +0.080 (P 0.069) | +0.121 (0.15) | +0.073 (0.40) | 0.87 |
| 16p11.2 dup | +0.052 (0.053) | 0.33 | +0.049 (0.46) | +0.075 (0.43) | NA (10 children) | 0.51 |
| 22q11.2 del | +0.012 (0.013) | 0.33 | -0.003 (0.88) | +0.016 (0.60) | +0.028 (0.15) | 0.49 |
| 22q11.2 dup | +0.049 (0.042) | 0.24 | +0.025 (0.63) | -0.005 (0.96) | +0.116 (0.12) | 0.44 |
| All four CNVs | +0.029 (0.011) | 0.010 | +0.024 (0.12) | +0.027 (0.31) | +0.039 (0.039) | 0.79 |

Heterogeneity across genotypes (all ages): P = 0.046. No evidence that the BMI effect differs between age windows.

Height: change in Z to the next visit per SD of BMI Z (`results_height.csv`)

| Genotype | All ages 2-18 y | P | 2-8 y | 8-13 y | 13-18 y | BMI x window P |
|---|---|---|---|---|---|---|
| 16p11.2 del | +0.088 (0.042) | 0.039 | +0.233 (P 3.5e-4) | +0.008 (0.90) | -0.024 (0.77) | 0.0082 |
| 16p11.2 dup | +0.022 (0.064) | 0.73 | +0.153 (0.20) | +0.111 (0.22) | NA (12 children) | 0.015 |
| 22q11.2 del | +0.012 (0.007) | 0.071 | +0.052 (1.7e-6) | +0.006 (0.56) | -0.022 (0.074) | 1.6e-5 |
| 22q11.2 dup | +0.065 (0.023) | 0.0041 | +0.127 (4.9e-4) | +0.045 (0.14) | +0.011 (0.80) | 0.065 |
| All four CNVs | +0.033 (0.008) | 2.0e-5 | +0.089 (1.4e-14) | +0.030 (0.0047) | -0.023 (0.076) | 4.5e-11 |

Heterogeneity across genotypes: all ages P = 0.43; 2-8 y P = 0.083; 8-13 y P = 0.070; 13-18 y P = 0.17.
The earlier finding holds: BMI predicts faster height gain before about 8 y (16p11.2 del, 22q11.2 del, and also 22q11.2 dup) and
not after. The 16p11.2 dup BMI x window P (0.015) rests on a 13-18 y window with 12 children and should not be leaned on.

Sex, exploratory (`sex_interaction.csv`): BMI x sex is not significant for head circumference in any genotype or pooled
(pooled: girls +0.042, boys +0.023, P = 0.34; 22q11.2 dup not estimable, 13 girls). For height, all ages, the pooled slope is larger
in boys (boys +0.047, girls +0.012, BMI x sex P = 0.024; per genotype P = 0.11-0.64; heterogeneity of BMI x sex across genotypes
P = 0.71). In the 2-8 y window the slopes are similar (girls +0.069, boys +0.093, P = 0.36); the difference sits at 8-13 y
(P = 0.053). Uncorrected for the number of tests.

## Sensitivity (`sensitivity.csv`; every row has all genotypes and windows)

Head circumference, all ages, BMI coefficient (P):

| Analysis | 16p del | 16p dup | 22q del | 22q dup | All four | Heterogeneity P |
|---|---|---|---|---|---|---|
| Primary | 0.087 (0.023) | 0.052 (0.33) | 0.012 (0.33) | 0.049 (0.24) | 0.030 (0.0089) | 0.030 |
| One interval per child (first) | 0.047 (0.57) | 0.002 (0.98) | 0.017 (0.57) | 0.138 (0.20) | 0.037 (0.17) | 0.21 |
| Intervals 0.5-3 y only | 0.144 (0.024) | 0.093 (0.29) | 0.025 (0.24) | 0.037 (0.66) | 0.051 (0.0052) | 0.020 |
| CHOP excluded | 0.087 (0.023) | 0.052 (0.33) | 0.192 (0.098; Cardiff, 50 children) | NA | 0.102 (5.7e-4) | 0.20 |
| CHOP only (22q11.2) | - | - | 0.005 (0.66) | 0.043 (0.32) | 0.009 (0.43) | 0.34 |
| Candidate duplicates removed (height + weight rule, 2 children) | same | same | 0.012 | 0.049 | 0.029 (0.010) | 0.029 |
| Candidate duplicates removed (loose upper bound, 66 Cardiff children) | same | same | 0.007 (0.60) | 0.048 (0.26) | 0.026 (0.023) | 0.020 |
| childsds scoring for CHOP | same | same | 0.012 (0.34) | 0.047 (0.27) | 0.030 (0.0097) | 0.027 |
| childsds scoring for every cohort | 0.082 (0.039) | 0.050 (0.35) | 0.012 (0.34) | 0.047 (0.27) | 0.029 (0.013) | 0.059 |
| CHOP from the analysis table (infant lengths under 60 cm lost) | same | same | 0.011 (0.42) | 0.053 (0.22) | 0.030 (0.010) | 0.015 |

Height, 2-8 y window, BMI coefficient (P):

| Analysis | 16p del | 22q del | 22q dup | All four |
|---|---|---|---|---|
| Primary | 0.233 (3.5e-4) | 0.052 (1.7e-6) | 0.127 (4.9e-4) | 0.089 (1.9e-14) |
| One interval per child (first) | 0.191 (0.024) | 0.045 (0.056) | 0.095 (0.18) | 0.070 (0.0017) |
| Intervals 0.5-3 y only | 0.087 (0.24) | 0.052 (4.4e-7) | 0.157 (9.1e-6) | 0.083 (4.4e-14) |
| CHOP excluded | 0.233 (3.5e-4) | NA | NA | 0.127 (0.0064) |
| CHOP only | - | 0.048 (2.3e-6) | 0.127 (3.1e-4) | 0.054 (2.3e-8) |
| Duplicates removed (any rule) | same | 0.048-0.056 | 0.121-0.127 | 0.089-0.097 |
| childsds scoring (CHOP or all) | 0.233-0.234 | 0.054 | 0.114 | 0.090-0.091 |
| Separate model per window (grant-code form) | 0.168 (0.018) | 0.044 (1.6e-4) | 0.170 (2.7e-4) | see CSV |

## The three known weaknesses

(a) **childsds vs `00_scoring.R`** (`scoring_agreement_childsds_vs_00_scoring.csv`): where both give a score they are the same
number (correlation 1.000; mean difference under 1e-5 SD, largest difference 0.003 SD, in every cohort and measure). The only
difference is coverage: childsds returns no WHO 2007 score at exactly 5.0 y, so it drops 16 Simons VIP BMI values at
head-circumference visits (and 11 CHOP, 17 Simons VIP height visits). `00_scoring.R` scores them. That is the whole reason the
16p11.2 numbers moved from the grant version (91 -> 94 children; 0.082, P 0.039 -> 0.087, P 0.023).

(b) **CHOP vs Cardiff double counting** (`duplicate_check_chop_cardiff.csv`; 108 Cardiff and 767 CHOP 22q11.2 children, no ID map):
- Rule as specified (same sex and genotype, a visit within 0.02 y with head circumference or height within 0.2 cm): 62 child
  pairs (41 Cardiff children), none matching at 2 or more visits. The same rule between children of opposite sex, who cannot be
  the same child, gives 72 pairs, so these matches are at chance level.
- Cardiff ages are whole years, so a looser rule was also run (age within -0.5 to +1 y): 14 same-sex pairs match at 3 or more
  Cardiff visits (5 in the opposite-sex control), but their weights disagree by 3-50 kg, so they are different children.
- Requiring height and weight to agree at the same visit: 2 pairs agree at 2 visits (0 in the control).
- Conclusion: no sign of systematic overlap; at most 2 candidate children. Removing them, or the upper-bound set of 66, changes
  nothing (tables above).

(c) **Ages covered** (`age_coverage.csv`): the head-circumference pairs are young. Half of all pairs start before 1.8 y, 52%
before 2 y, 72% before 5 y; 16p11.2 pairs are 69% before 2 y (median 1.1-1.2 y, 90% before 6.4 y and 4.0 y). Only 14% (16p11.2
del) and 7% (16p11.2 dup) of pairs start at 5 y or later, and those come mostly from Cardiff with 3-year gaps. The
head-circumference result is a statement about infancy and early childhood; it says little about school age in 16p11.2.

## Comparison with the grant-figure version (`compare_with_grant_version.csv`)

| | Grant version (A7/A8) | v2 primary |
|---|---|---|
| Children / pairs, head circumference | 513 / 2,070 | 518 / 2,103 |
| 16p11.2 del | 0.082, P 0.039 | 0.087, P 0.023 |
| 16p11.2 dup | 0.050, P 0.35 | 0.052, P 0.33 |
| 22q11.2 del | 0.012, P 0.34 | 0.012, P 0.33 |
| 22q11.2 dup | 0.048, P 0.26 | 0.049, P 0.24 |
| All carriers | 0.029, P 0.012 | 0.030, P 0.0089 |
| Heterogeneity P | 0.060 | 0.030 |

With childsds scoring for every cohort the rebuild returns the grant numbers (513 children, 2,069 pairs, 0.029, P 0.013,
heterogeneity P 0.059), so the code is reproduced and the differences are the age-5.0 scoring gap alone.

## Caveats

- Within-carrier estimates only; no control children with repeated visits.
- Head circumference: the pooled association is small (0.03 SD per SD of BMI per interval) and is carried by 16p11.2 deletion
  (Simons VIP infants) and Cardiff. In CHOP alone (321 children, 1,580 pairs) it is 0.009, P = 0.43. It is not significant with
  one interval per child (pooled P = 0.17; 16p11.2 del 0.047, P = 0.57). The heterogeneity P ranges 0.015-0.06 across reasonable
  data choices (0.20-0.34 inside non-CHOP or CHOP alone), and genotype is confounded with cohort, so the heterogeneity cannot be
  assigned to genotype rather than cohort or age.
- Height: the 2-8 y effect is robust in 22q11.2 (CHOP). In 16p11.2 del it is 0.23 in the primary model, 0.17 with a separate
  model per window, 0.19 with one interval per child, and 0.087 (P = 0.24) when intervals are limited to 0.5-3 y.
- Cardiff ages are whole years and its intervals are long (median 3 y); CHOP is clinical EHR data with visit timing driven by care.
- `00_qc_growth.R` rule 2 (height under 60 cm removed) also affects any other use of `analysis_table_chop` below about 4 months
  of age (for example a 0-2 y cross-sectional view): infant BMI and length are missing there for CHOP.
- P values are not corrected for the number of genotypes, windows or sex tests.

## Files

- Figure blocks: `block_hc_longitudinal.{png,pdf,csv}` (tags g, h, i), `block_height_longitudinal.{png,pdf,csv}` (tags g, h).
- Results: `results_head_circumference.csv`, `results_height.csv`, `results_by_window_separate_models.csv`, `sex_interaction.csv`.
- Samples: `sample_head_circumference_by_genotype_cohort.csv`, `sample_height_by_genotype_cohort.csv`,
  `sample_by_genotype_cohort_window.csv`, `age_coverage.csv`, `carrier_repeat_visits_by_cohort.csv`, `controls_repeat_check.csv`.
- Checks: `sensitivity.csv`, `scoring_agreement_childsds_vs_00_scoring.csv`, `duplicate_check_chop_cardiff.csv`,
  `cleaning_log.csv`, `compare_with_grant_version.csv`, `80_run.log`.
