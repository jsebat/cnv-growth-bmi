# Work package C: autistic traits and autism diagnosis (CNV x BMI, x polygenic score, x sex)

> **Update 2026-11-21 (JS): cohort adjustment and Fig. 4 panels.** (1) Pooled autism-diagnosis models (`85_v2_asd_by_cnv.R`, Fig. 3)
> now carry one intercept per cohort within CNV instead of additive cohort + CNV terms, because cohort and CNV are partly confounded
> (CHOP has only 22q11.2 carriers). Per-CNV models are unchanged. Pooled results: head circumference by diagnosis +0.21 (P 0.043) ->
> +0.20 (P 0.051; family-clustered 0.061); odds ratio per SD of head circumference 1.13 (P 0.026) -> 1.13 (P 0.036; clustered 0.053);
> BMI by diagnosis +0.17 (P 0.056 -> 0.057); odds ratio per SD of BMI 1.08 (P 0.094 -> 0.097, n 1,482 -> 1,474 because cohort x CNV
> cells without both outcomes drop out). The earlier model is a sensitivity row in both CSVs and the earlier outputs are in
> `previous_additive_cohort_2026-11-20/`. Numbers quoted further down in this README for those pooled models are the earlier ones.
> (2) Fig. 4 (`86b`): panel a prints the BMI x score and three-way tests and is narrower; panel c is a bar plot of the three-way
> term by CNV and pooled; the full forest is `fig_supp_threeway_forest.png` (Supplementary Fig. 6).


BMI modifier paper v2, run of 2026-11-20. Everything here is aggregate; participant-level tables and caches stay in
`Jon/Claude/CNV_growth_trajectories_data` (cache folder `v2_autism_cache/`).

**Bottom line**

- **CNV x BMI on autistic traits (primary): unchanged.** Re-running the old code on the current data reproduces every number.
  SRS rises 0.07 SD more per SD of BMI in carriers than in relatives (P 0.013); SCQ 0.06 (P 0.063).
- **Autism diagnosis (secondary) does not clearly rise with BMI.** Pooled odds ratio 1.08 per SD of BMI (P 0.094), and the signal
  comes from CHOP alone. Diagnosis does rise with head circumference (odds ratio 1.13, P 0.026).
- **CNV x BMI x autism polygenic score holds for SRS** (-0.075, P 0.012; permutation P 0.018) **but not for SCQ**, rests on G2MH
  (81% of the sample), and differs between CNVs (heterogeneity P 0.009; 16p11.2 duplication goes the other way).
- **Mendelian randomisation is inconclusive in carriers.**
- **The sex difference does not survive a robustness check.** The pooled CNV x BMI x sex term for SRS (P 0.009) falls to P 0.12
  when each sex gets its own cohort x CNV intercepts, and to P 0.28 when the four per-CNV estimates are averaged.

Sections 3 and 4 are exploratory; no P value is corrected for multiple testing.

## Files

| Deliverable | Figure | Tables | Script |
|---|---|---|---|
| 1. SRS / SCQ block (tags j-s) | `block_srs_scq.png` | `block_srs_scq.csv`, `block_srs_scq_binned_means.csv`, `block_srs_scq_n_by_cohort.csv` | `84_v2_block_srs_scq.R` |
| 2. Autism diagnosis by CNV (tags a-e) | `fig_asd_by_cnv_v2.png` | `asd_by_cnv_tests_v2.csv`, `asd_by_bmi_stratum.csv`, `asd_by_bmi_logistic.csv`, `asd_by_cnv_n_by_cohort_sex.csv`, `fig_asd_by_cnv_v2_plotted.csv` | `85_v2_asd_by_cnv.R` |
| 3. Polygenic score, main (tags a-d) | `fig_cnv_bmi_pgs_main.png` | `threeway_pgs_robustness.csv`, `pgs_bmi_slope_by_score_fifth.csv`, `mr_bmi_autistic_traits_v2.csv`, `asd_dx_bmi_x_pgs_logistic.csv`, `fig_cnv_bmi_pgs_main_plotted.csv` | `86_v2_cnv_bmi_pgs_stats.R`, `86b_v2_fig_cnv_bmi_pgs_main.R` |
| 3. Polygenic score, supplement (tags a-g) | `fig_supp_cnv_bmi_pgs_full.png` | `fig_supp_cnv_bmi_pgs_full_plotted.csv` | `86a_v2_rerun_49_pgs_adjustment.R`, `86c_v2_fig_supp_cnv_bmi_pgs_full.R` |
| 4. Sex block (tags e-g) | `block_sex_autism.png` | `sex_cnv_bmi_autism.csv`, `block_sex_autism_plotted.csv` | `87_v2_sex_cnv_bmi_autism.R` |
| Samples | | `sample_n_by_cnv_cohort_sex.csv` | `87_v2_sex_cnv_bmi_autism.R` |

Every figure also has a PDF. `83_v2_autism_common.R` holds the shared helpers; `83b_v2_verify_combined_build.R` checks the SRS / SCQ
table against a fresh build. `aux_rerun/` holds the outputs of the old scripts as re-run here (nothing in the code folder or in
older run folders was overwritten) and the run logs.

Run order: 84 (builds the cache, about 4 min), 85, 86a, 86, 86b, 86c, 87. `V2_REBUILD=TRUE` forces a rebuild of the cache.

## Samples

Full counts per outcome x CNV x cohort x sex are in `sample_n_by_cnv_cohort_sex.csv`. Totals:

| Analysis | 16p11.2 del | 16p11.2 dup | 22q11.2 del | 22q11.2 dup | Cohorts | Reference |
|---|---|---|---|---|---|---|
| SRS | 279 | 140 | 317 | 78 | G2MH 522, Simons VIP 159, Searchlight 116, SSC 17 | 5,609 SSC unaffected relatives |
| SCQ | 259 | 98 | 275 | 57 | G2MH 423, Simons VIP 125, Searchlight 131, SSC 10 | 1,174 SSC siblings |
| SRS with autism score | 145 | 80 | 289 | 69 | G2MH 471, Simons VIP 99, SSC 13 | 5,120 |
| SCQ with autism score | 122 | 45 | 250 | 51 | G2MH 382, Simons VIP 77, SSC 9 | 1,068 |
| Diagnosis with BMI and height | 252 | 174 | 891 | 165 | CHOP 728, G2MH 444, Simons VIP 194, Cardiff 74, CHUSJ 24, SSC 18 | 3,727 SSC probands and siblings |
| Diagnosis with head circumference | 213 | 137 | 550 | 122 | CHOP 454, G2MH 310, Simons VIP 162, Cardiff 73, SSC 18, CHUSJ 5 | 3,707 |

- 22q11.2 carriers have SRS / SCQ in G2MH only (plus one SSC carrier); Simons VIP and Searchlight carry 16p11.2 only.
- The SRS / SCQ sample includes the G2MH carriers recruited at the CHOP site. The growth table used for diagnosis drops them
  (CHOP overlap rule of the 2026-11-19 run) and adds the CHOP records. CHOP has no SRS, SCQ or polygenic score.
- Searchlight has no autism diagnosis in the 32 pipeline and no head circumference; a sensitivity row adds Searchlight v14
  diagnoses (clinical or reported).
- The SRS / SCQ table (`combined_srs_scq.rds`, written 2026-11-14) was rebuilt in memory from the current code and source files and
  is identical: 8,537 rows, zero difference in trait score, BMI Z and age (`aux_rerun/verify_combined_build.txt`).
- One visit per person: SRS / SCQ with the BMI nearest the questionnaire (within 3 y); growth table by rule `fill_gaps`.

## Models

- **Scale.** SRS and SCQ raw totals are Box-Cox transformed and z-scored against SSC unaffected relatives (as 51b / 60). The ordinal
  check is a cumulative-logit model of the raw total.
- **Carrier BMI slope:** `gam y ~ s(age) + BMI (+ cohort x CNV)`.
- **CNV x BMI:** `bam y ~ carrier + s(age, by = carrier) + s(age, by = BMI) + carrier:BMI (+ cohort x CNV)`, reference SSC
  unaffected relatives (61 / 51b).
- **CNV x BMI x autism score:** `lm y ~ carrier * (age + age^2 + sex) + cohort x CNV + carrier * BMI * score`, HC1 robust SE (62 /
  64). Score = April 2026 autism GWAS, full weights for G2MH and Simons VIP, leave-SSC-out weights for SSC; residualised on 10 PCs
  and z-scored within cohort.
- **Diagnosis by CNV (panels a-c):** `gam y ~ diagnosis + s(age) + sex (+ cohort) (+ CNV)` (32).
- **Diagnosis by BMI or head circumference (new):** logistic `diagnosis ~ Z + ns(age, 3) + sex + cohort (+ CNV)`. A cohort without
  both outcomes drops out once it has its own intercept.
- **Sex (new, one set):** BMI slope per sex; CNV x BMI per sex against same-sex relatives; one model with carrier x BMI x sex and all
  lower-order terms (`bam`, sex-specific age smooths); logistic `diagnosis ~ BMI x sex` within carriers.
- **Heterogeneity across the four CNVs:** Cochran Q on the four per-CNV estimates.
- **Permutation null (three-way term):** the autism score is permuted within cohort x CNV among carriers and, separately, among
  relatives; 10,000 permutations; statistic = robust t.
- **Mendelian randomisation:** instruments = polygenic scores for childhood obesity, childhood BMI and adult BMI. Two-stage least
  squares with HC1 SE, Wald ratio per score, leave-one-instrument-out, Anderson-Rubin 95% sets, Sargan over-identification test.

## Results

### 1. SRS / SCQ against BMI (`block_srs_scq.csv`)

No number moved from the 2026-11-14 values (column `moved` is FALSE in every row).

| Group | n | SRS slope | SRS CNV x BMI (P Box-Cox / ordinal) | n | SCQ slope | SCQ CNV x BMI (P Box-Cox / ordinal) |
|---|---|---|---|---|---|---|
| SSC unaffected relatives | 5,609 | +0.02 | reference | 1,174 | +0.00 | reference |
| SSC autistic probands, no CNV | 2,065 | +0.01 | +0.00 (0.88) | 1,503 | +0.00 | +0.00 (0.99) |
| All four CNVs | 814 | +0.10 | +0.072 (0.013 / 0.088) | 689 | +0.06 | +0.061 (0.063 / 0.086) |
| 16p11.2 del | 279 | +0.05 | +0.033 (0.44 / 0.97) | 259 | +0.00 | -0.001 (0.98 / 0.55) |
| 16p11.2 dup | 140 | +0.09 | +0.075 (0.25 / 0.20) | 98 | +0.11 | +0.102 (0.15 / 0.006) |
| 22q11.2 del | 317 | +0.15 | +0.110 (0.009 / 0.011) | 275 | +0.09 | +0.088 (0.082 / 0.12) |
| 22q11.2 dup | 78 | -0.02 | -0.055 (0.49 / 0.66) | 57 | +0.19 | +0.179 (0.053 / 0.010) |

- The pooled SRS result is nominal on the Box-Cox scale and weaker on the ordinal check (P 0.088).
- 22q11.2 deletion is the only single CNV where both scales agree for SRS.

### 2. Autism diagnosis within CNV group (`asd_by_cnv_tests_v2.csv`, `asd_by_bmi_logistic.csv`)

**Panels a-c, carriers with vs without a diagnosis (adjusted difference in SD):**

| Outcome | All four CNVs (n no dx / dx) | P | Heterogeneity P | SSC probands vs siblings |
|---|---|---|---|---|
| Height | +0.13 (1,137 / 345) | 0.11 | 0.44 | -0.08, P 0.072 |
| Head circumference | +0.21 (736 / 286) | 0.043 | 0.76 | +0.14, P 0.0016 |
| BMI | +0.17 (1,137 / 345) | 0.056 | 0.73 | +0.12, P 0.018 |

- Per-CNV rows match the 2026-11-19 run exactly. Overgrowth / undergrowth groupings are in the table only.
- BMI by cohort: CHOP +0.34 (P 0.015); without CHOP +0.08 (P 0.50).
- Head circumference without each cohort in turn: +0.18 to +0.26 (P 0.04 to 0.11).

**Panel d, diagnosis by BMI band (pooled):** 21.6% (85 / 393) at or below the 30th centile, 18.7% (77 / 411) in the middle, 27.0%
(183 / 678) at or above the 70th.

| Analysis | n (diagnosed) | Odds ratio per SD of BMI (95% CI) | P |
|---|---|---|---|
| All four CNVs, cohort terms | 1,482 (345) | 1.08 (0.99-1.18) | 0.094 |
| 16p11.2 del / dup | 248 / 174 | 1.04 / 1.05 | 0.70 / 0.75 |
| 22q11.2 del / dup | 891 / 161 | 1.13 / 1.02 | 0.059 / 0.88 |
| Heterogeneity across CNVs | | | 0.83 |
| CHOP only | 728 (112) | 1.19 (1.02-1.39) | 0.027 |
| G2MH only | 444 (139) | 1.01 (0.88-1.16) | 0.90 |
| Simons VIP only | 194 (35) | 0.91 (0.70-1.19) | 0.48 |
| Cardiff only | 74 (41) | 1.18 (0.83-1.66) | 0.36 |
| Without CHOP | 754 (233) | 1.02 (0.91-1.14) | 0.78 |
| Without cohort terms | 1,482 | 1.04 (0.95-1.13) | 0.40 |
| Age under 18 | 1,040 (256) | 1.11 (1.00-1.23) | 0.059 |
| With Searchlight diagnoses | 1,639 (426) | 1.06 (0.98-1.15) | 0.14 |

- **Depends on one cohort.** Any BMI signal is in CHOP; the other cohorts together show none. The CHOP index sheet does not say
  whether "No" means evaluated, so a recording difference by BMI cannot be excluded.

**Panel e, diagnosis by head-circumference band (pooled):** 24.1% (133 / 551), 30.8% (77 / 250), 34.4% (76 / 221). Odds ratio per SD
1.13 (1.02-1.27), P 0.026; heterogeneity across CNVs P 0.88; 1.11 to 1.16 with each cohort left out (P 0.03 to 0.09).

### 3. CNV x BMI x autism polygenic score (`threeway_pgs_robustness.csv`)

The re-run of 66 and of 49 is identical to the saved tables. No error was found in the old statistics.

| SRS three-way term | n carriers | Estimate (SE) | P |
|---|---|---|---|
| All carriers vs SSC relatives | 583 | -0.075 (0.030) | 0.012; permutation 0.018 |
| Within carriers, BMI x score | 583 | -0.085 (0.027) | 0.0019; permutation 0.0038 |
| 16p11.2 del | 145 | -0.072 (0.062) | 0.24 |
| 16p11.2 dup | 80 | +0.119 (0.058) | 0.039 |
| 22q11.2 del | 289 | -0.100 (0.038) | 0.009 |
| 22q11.2 dup | 69 | -0.135 (0.077) | 0.079 |
| Heterogeneity across CNVs | | | 0.009 |
| G2MH only / Simons VIP only | 471 / 99 | -0.077 / -0.063 | 0.015 / 0.44 |
| Without G2MH | 112 | -0.064 (0.079) | 0.42 |
| Without Simons VIP / SSC carriers | 484 / 570 | -0.077 / -0.075 | 0.014 / 0.013 |
| Without 16p11.2 del / dup | 438 / 503 | -0.071 / -0.084 | 0.042 / 0.007 |
| Without 22q11.2 del / dup | 294 / 514 | -0.065 / -0.077 | 0.11 / 0.013 |
| SCQ, all carriers | 468 | -0.026 (0.037) | 0.49; permutation 0.49 |

- BMI slope on SRS: +0.20 in carriers below the median score, +0.00 above. By fifth of the score: +0.30, +0.17, +0.05, -0.04, +0.01
  (`pgs_bmi_slope_by_score_fifth.csv`). Relatives: flat (BMI x score -0.010, P 0.45).
- **Flags.** Not seen for SCQ. The estimate keeps its sign and size without G2MH but only 112 carriers remain. 16p11.2 duplication
  has the opposite sign. Without 22q11.2 deletion the term is no longer significant.

**Autism diagnosis (`asd_dx_bmi_x_pgs_logistic.csv`).** Carriers of the SRS / SCQ sample with a score (n = 586, 162 diagnosed): BMI
x score odds ratio 0.92 (0.82-1.04), P 0.19. Same direction as SRS, not significant. Growth-table carriers (n = 571): 0.91, P 0.14.
16p11.2 duplication again goes the other way (odds ratio 2.8, P 0.008, 27 diagnosed of 81); heterogeneity P 0.033.

**Mendelian randomisation (`mr_bmi_autistic_traits_v2.csv`), effect of BMI in SD of trait per SD of BMI:**

| Stratum | n | Observational | MR, three scores (SE) | First-stage F | Instruments disagree, P | Robust 95% set | Detectable at 80% power |
|---|---|---|---|---|---|---|---|
| Carriers, all genotyped, SRS | 583 | +0.10 | +0.21 (0.12), P 0.072 | 14.7 | 0.094 | 0.00 to 0.47 | 0.33 |
| Carriers, G2MH, same model, SRS | 471 | +0.12 | +0.18 (0.13) | 9.6 | 0.019 | not reliable | 0.37 |
| Carriers, G2MH, model of 49, SRS | 471 | +0.12 | +0.02 (0.13) | 9.5 | 0.007 | empty | 0.36 |
| Carriers, all genotyped, SCQ | 468 | +0.07 | +0.04 (0.09) | 13.1 | 0.98 | -0.27 to 0.33 | 0.25 |
| SSC unaffected relatives, SRS | 5,120 | +0.03 | +0.05 (0.04) | 212 | 0.55 | -0.05 to 0.14 | 0.11 |
| SSC autistic probands, no CNV, SRS | 1,888 | +0.01 | +0.02 (0.04) | 59 | 0.059 | -0.04 to 0.08 | 0.11 |

- **Inconclusive in carriers.** The smallest effect detectable at 80% power (0.33) is three times the observed slope (0.10); power to
  detect an effect of the observed size is 14%.
- The carrier estimate depends on the covariate set: +0.21 with height and autism scores as covariates, +0.02 with the seven score
  covariates of 49 (G2MH only, where the instruments also disagree).
- The childhood BMI score alone gives a much larger estimate (+0.50) than the other two (+0.18 each), which is what pleiotropy of
  that score would look like.
- Adding Simons VIP and SSC carriers raises the first-stage F from 9.6 to 14.7. Simons VIP scores come from Anjali's Searchlight
  set (same weights as G2MH).
- In relatives and in probands without a CNV the instruments are strong and the estimate is near zero, with an upper bound of 0.14
  and 0.08.

### 4. Sex (`sex_cnv_bmi_autism.csv`)

| SRS, all four CNVs | Female | Male | Male minus female | P |
|---|---|---|---|---|
| BMI slope, carriers (n 397 / 417) | +0.16 (P 0.0001) | +0.06 (P 0.085) | -0.10 (separate fits) | 0.049 |
| BMI slope, relatives (n 2,875 / 2,734) | +0.03 | +0.01 | -0.02 | 0.38 |
| CNV x BMI vs same-sex relatives | +0.126 (P 0.004) | +0.041 (P 0.30) | | |
| CNV x BMI x sex, one model | | | -0.137 (0.052) | 0.009 |
| Same, linear model with robust SE | | | -0.129 (0.054) | 0.018 |
| Same, sex-specific cohort x CNV intercepts | | | -0.086 (0.055) | 0.12 |
| Mean of the four per-CNV estimates | | | -0.056 (0.051) | 0.28 |

- Per CNV the three-way term is -0.10, -0.04, -0.03 and +0.04 (none significant; heterogeneity P 0.89). The pooled estimate is
  larger than every per-CNV estimate, which points to between-CNV differences in the sex gap leaking into the pooled term.
- **Does not survive the robustness check.** Direction is consistent (steeper in females); the evidence is weak.
- Leaving one cohort out: -0.17 (without G2MH), -0.12 (Searchlight), -0.13 (SSC), -0.14 (Simons VIP); P 0.016 to 0.039.
- SCQ: female +0.096 (P 0.056), male +0.032 (P 0.46); three-way -0.070 (P 0.26).
- **Diagnosis** (n = 1,482): odds ratio per SD of BMI 1.16 in females (P 0.086) and 1.06 in males (P 0.32); BMI x sex P 0.10.
  Without CHOP P 0.035; without G2MH P 0.58. In the genotyped SRS / SCQ sample (n = 586) BMI x sex P 0.029. Nominal at best.
- **Polygenic score.** Adjusting for the score, BMI x score and their carrier interactions leaves the sex term unchanged (-0.167 to
  -0.162 in the genotyped sample). BMI x score is similar in female (-0.095) and male (-0.077) carriers; difference P 0.77.

## Sensitivity analyses and whether they changed the conclusion

| Analysis | Changed the conclusion? |
|---|---|
| Ordinal instead of Box-Cox scale (CNV x BMI) | Weaker for pooled SRS (P 0.088); 22q11.2 deletion holds |
| Diagnosis by BMI: per cohort, without CHOP, without cohort terms | Yes: no signal outside CHOP |
| Diagnosis by head circumference: leave one cohort out | No |
| Three-way score term: permutation | No (P 0.018) |
| Three-way score term: leave one cohort out | Sign and size hold; not significant without G2MH (n = 112) |
| Three-way score term: leave one CNV out | Not significant without 22q11.2 deletion (P 0.11) |
| Three-way score term: per CNV | Yes: heterogeneous, 16p11.2 duplication opposite |
| Mendelian randomisation: covariate set, leave one instrument out, robust sets | Yes: estimate moves from +0.21 to +0.02 |
| Sex term: sex-specific cohort x CNV intercepts, per-CNV mean | Yes: P 0.009 becomes 0.12 and 0.28 |
| Sex term: adjustment for the autism score | No |

## Family-clustered check (`family_clustered_check.csv`, `87b_v2_family_clustered_check.R`)

No conclusion changes: no headline P crosses 0.05 when standard errors are clustered by family. Clustered SEs are 0.98 to 1.06
times the reported ones. Figures are unchanged.

| Test | Estimate | Reported SE, P | Family-clustered SE, P | Carriers sharing a family / all |
|---|---|---|---|---|
| CNV x BMI, SRS, all four CNVs | +0.072 | 0.029, 0.013 | 0.031, 0.018 | 180 / 814 |
| CNV x BMI, SCQ, all four CNVs | +0.061 | 0.033, 0.063 | 0.032, 0.058 | 79 / 689 |
| CNV x BMI, SRS, 22q11.2 del | +0.110 | 0.042, 0.009 | 0.045, 0.014 | 12 / 317 |
| CNV x BMI x autism score, SRS, vs relatives | -0.075 | 0.030, 0.012 | 0.030, 0.012 | 117 / 583 |
| BMI x autism score, SRS, within carriers | -0.085 | 0.027, 0.0019 | 0.027, 0.0018 | 106 / 583 |
| CNV x BMI x autism score, SCQ | -0.026 | 0.037, 0.49 | 0.038, 0.49 | 40 / 468 |
| CNV x BMI x sex, SRS, one model | -0.137 | 0.052, 0.009 | 0.055, 0.013 | 180 / 814 |
| Same, sex-specific cohort x CNV intercepts | -0.086 | 0.055, 0.12 | 0.057, 0.13 | 180 / 814 |
| Head circumference by diagnosis, pooled | +0.21 | 0.104, 0.043 | 0.107, 0.0495 | 167 / 1,022 |
| Diagnosis on head-circumference Z (log odds ratio) | +0.126 | 0.057, 0.026 | 0.060, 0.035 | 167 / 1,022 |
| BMI by diagnosis, pooled | +0.17 | 0.091, 0.056 | 0.092, 0.058 | 220 / 1,482 |
| Diagnosis on BMI Z (log odds ratio) | +0.078 | 0.046, 0.094 | 0.046, 0.094 | 220 / 1,482 |

- "Sharing a family" counts carriers whose family has another person (carrier or relative) in the same analysis.
- Family = SSC family number; Simons VIP / Searchlight family number; G2MH and Cardiff `rarecnv_family_id` of the G2MH release
  (else site + `local_family_id`). Every G2MH, Cardiff, Simons and SSC carrier had a family identifier.
- CHOP (728) and CHU Sainte-Justine (24) records carry no family identifier, so each is its own cluster; relatedness within CHOP is
  not accounted for.
- The head-circumference contrast sits at P 0.0495 after clustering: it stays under 0.05 but only just.
- Method: `sandwich::vcovCL` for the linear and logistic models; for the penalised `bam` / `gam` models, the sandwich built from
  the model's penalised bread and per-family score sums.

## Caveats

- Reported standard errors are not clustered by family; the check above shows this does not change any headline conclusion.
- Diagnosis is recorded differently in every cohort: research assessment in SSC (probands by design), clinical summary in Simons
  VIP, subject table in G2MH and Cardiff, index sheet in CHOP. Cohort intercepts absorb level differences, not differences in how
  recording relates to BMI.
- Polygenic-score analyses cover G2MH, Simons VIP and SSC carriers only. Searchlight, Cardiff, CHUSJ and CHOP carriers have no score.
- 22q11.2 results for SRS / SCQ and for the score are G2MH only, so CNV and cohort cannot be separated for 22q11.2.
- The 49 re-run uses `ADD_CHOP=FALSE`, as when 49 was first run: its only growth input is the carriers' own measured height.
- The Wald-ratio standard error is the reduced-form SE divided by the first stage (as in 49), which understates uncertainty when
  the instrument is weak; the Anderson-Rubin sets are the ones to read.
