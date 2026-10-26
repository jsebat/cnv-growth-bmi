# Growth and cognitive trajectories by BMI: carriers vs unaffected (non-ASD) relatives

JS / Claude, 2026-10-24. Script `23_unaffected_controls.R` (sources the data build of `21_pheno_bmi_gam.R`).
Outputs: `unaffected_controls_tests.csv`, `fig_unaffected_controls.{png,pdf}`, `fig_carriers_only_trajectories.{png,pdf}`.
Run log: `run_unaffected.log` (the superseded first run, with wrong sibling ages, is `run_unaffected.snapshotage.log`).

## 1. Design

* **Carriers:** 16p11.2 proximal del/dup and 22q11.2 A-D del/dup from the clinical growth genomatrix (G2MH, Cardiff,
  SVIP, Searchlight; Cardiff-masked G2MH rows removed), ages 2-65, latest visit with the outcome and BMI. ASD status from
  `Jon/Claude/G2MH_CNV_growth/subjects.csv` (G2MH) and SVIP `diagnosis_summary.csv`.
* **Unaffected controls (no ASD):** SSC designated siblings and parents (non-carriers in the genomatrix; sibling
  height/weight/HC from the SSC sibling table) and SVIP non-carrier family members with clinical ASD = false.
  No SSC probands, no UK Biobank (UKBB can be added back to the height/HC controls later).
* **Outcomes:** height and head circumference (WHO / UK90 Z), SRS total (Z vs SSC unaffected), SCQ lifetime
  (Z vs SSC siblings); carriers only: RBS-R (29 items), FSIQ, Penn CNB accuracy and speed (no controls with BMI).
* **BMI:** child WHO Z; adult centile vs the pooled non-carrier adult reference, mapped to Z. Plots split at the
  70th / 30th centiles; tests use continuous BMI.
* **Tests (mgcv GAM/bam):** within group, BMI main effect (`y ~ s(age) + BMI`) and age x BMI (`ti(age, BMI)`).
  Carriers vs unaffected: CNV x BMI constant shift (`y ~ cnv + s(age, by = cnv) + s(age, by = BMI) + cnv:BMI`) and
  CNV x age x BMI (`ti(age, BMI, by = cnv)`). Each run for all carriers and for non-ASD carriers only, per locus and
  pooled (overgrowth = 16p del + 22q dup, undergrowth = 16p dup + 22q del, all), with locus as a covariate in pooled fits.

## 2. FLAG: age at measurement in the SSC phenotype release

The SSC v15 age columns (`Age_individual_v14_parent_age.txt` `age`; `SSC_Age_of_Subjects_&_ParentalAgeAtBirth.csv`
`subject_age`) are ages **as of a single snapshot date, September 2012**, not ages at the study visit (2008-2011).
They overstate visit age by a median of 3.2 years (range 1.4-6.6). Consequences and fixes:

| Where | Effect of the wrong age | Status |
|---|---|---|
| SSC siblings, height/HC/BMI Z (this analysis, first run) | Height Z of child controls ~ -2.5 SD; HC Z ~ -1.5 | Fixed: sibling visit age = proband age at evaluation (`ssc_case_HC.csv`) + (proband birth month - sibling birth month). Reproduces SSC's own sibling height Z (r = 0.96). |
| SSC parents in Molly's genomatrix (`age` column) | Identical to the snapshot age, i.e. ~3 y too old. Adult height Z and adult HC Z do not depend on age; adult BMI centile does weakly. Age axis of all parent points shifted +3 y. | Fixed in 21/23 (`gm` ages replaced by visit age for 4452 parent rows). **Not fixed in the genomatrix itself, in `analysis_table.rds`, or in the adult BMI reference `REFad` used by `bmiZ()`.** Small effect, but Molly should know. |
| SSC probands | Genomatrix age = `age_at_ados` = age at evaluation. Correct. | OK |
| Cognitive traits (SRS, SCQ) of SSC siblings and parents | Scores are not re-normed by age here (Z vs the unaffected pool), so the trait values are unchanged, but the age axis and the age smooth were off by ~3 y for all SSC controls in the `21_pheno_bmi_gam.R` outputs (`pheno_bmi_gam_tests.csv`, `fig_supp_pheno_bmi_gam.*`), and sibling BMI Z (the exposure) was scored at the wrong age. Any age-normed cognitive score computed by us from SSC raw data would be affected directly. | Fixed in the code; **21 has not been rerun** (~70 min for both passes). Its SSC-proband-control results are superseded by this analysis anyway. |
| `18_ssc_sibpairs.R` | Already used the visit-age derivation. | OK |
| Main UKBB figure (`17_main_with_ukbb.R`) | Clinical No CNV adults include SSC parents at snapshot age; adult height/HC Z are age-invariant, adult BMI centile weakly age-adjusted. | Expected to be negligible; not rerun. |

Documented for the data folder in `...\Autism\SSC_Version_15_Phenotype_Data_Set_9\README_ages_and_line_endings.md`.
Age-dependent measures from SSC (anything scored with age norms by us rather than by SSC) must use the visit age,
never these columns.

## 3. Sample (with BMI; unaffected = non-ASD)

| Outcome | 16p del | 16p dup | 22q del | 22q dup | Unaffected |
|---|---|---|---|---|---|
| Height | 359 | 200 | 361 | 119 | 5718 (1696 sibs, ~3800 parents, ~220 SVIP) |
| Head circumference | 213 | 126 | 243 | 92 | 5691 |
| SRS | 79 | 34 | 209 | 67 | 5493 |
| SCQ | 60 | 21 | 171 | 50 | 1134 (siblings only) |
| RBS-R / FSIQ | 70 / 56 | 31 / 25 | 216 / 111 | 63 / 36 | none |
| CNB accuracy / speed | 141 / 143 | 57 / 57 | 348 / 354 | 142 / 142 | none (CNB controls have no BMI) |

Non-ASD carriers are roughly half to two thirds of each group (height: 149, 102, 217, 82).

## 4. Results, carriers vs unaffected controls (`fig_unaffected_controls`)

CNV x BMI = age-adjusted difference in the BMI slope (SD of outcome per SD of BMI) between carriers and unaffected;
3-way = CNV x age x BMI.

**Height**
* No locus shows a CNV x BMI shift (16p del p = 0.71, 16p dup 0.17, 22q del 0.065, 22q dup 0.52). Undergrowth pooled
  -0.09, p = 0.019, driven by 22q del and absent in non-ASD carriers (p = 0.35).
* Age-dependence is where the signal is: CNV x age x BMI p = 4.4e-6 (16p del), < 1e-16 (16p dup), 0.013 (22q del),
  0.063 (22q dup). Within carriers the age x BMI term is significant for 16p del (p = 0.007), 22q del (0.046),
  overgrowth pooled (6.7e-4) and all pooled (2.8e-6). It is also very strong within the unaffected (p ~ 0) because the
  BMI-height relation changes sharply across childhood. The 16p dup three-way p of 0 is a bam fit of 200 carriers
  against 5700 controls and should be read as "very small" rather than quoted.
* Non-ASD carriers only: 16p del CNV x BMI +0.13 (p = 0.056), 16p dup three-way p = 0.038, 22q del three-way n.s.
  (0.90). The 22q del age x BMI pattern therefore comes mainly from the ASD subgroup.

**Head circumference**
* 16p del: carriers' HC rises more steeply with BMI than in unaffected relatives (CNV x BMI +0.16, p = 0.0053; three-way
  p = 0.015), and this holds in non-ASD carriers (+0.19, p = 0.010). Overgrowth pooled +0.11, p = 0.020.
* 16p dup, 22q del, 22q dup: no CNV x BMI shift (p = 0.22, 0.72, 0.65). Within every carrier group HC increases with BMI
  (slopes 0.26-0.50 SD/SD, all p < 0.005), as it does in the unaffected (0.37, p ~ 1e-139): the BMI-HC relation is
  general, and only 16p del exceeds it.

**SRS** (higher = more impaired)
* Carriers' SRS rises with BMI more than unaffected: all CNVs +0.20 SD/SD, p = 6e-6; 22q del +0.20 (5e-5);
  16p dup +0.50 (4e-5, n = 34); 16p del +0.17 (0.044); 22q dup none (p = 0.66). No age dependence (all three-way
  p > 0.1).
* Restricting to non-ASD carriers weakens it (all CNVs +0.08, p = 0.11; 22q del +0.17, p = 0.004) and reverses it in
  16p del (-0.20, p = 0.049, n = 43). The within-unaffected BMI effect on SRS is nil (0.01, p = 0.28).

**SCQ**
* Same picture: all CNVs +0.33, p = 5e-5; 22q del +0.27 (1e-3); 22q dup +0.39 (2e-3, and +0.50 in non-ASD, p = 4e-4);
  16p del +0.26 (0.034, nil in non-ASD). 16p dup too small (n = 21). No age dependence.
* Caveat: SCQ (and SRS) contributed to ASD case assignment in some cohorts, so "non-ASD carriers" is a selection on the
  outcome; both views are reported.

## 5. Results, carriers only (`fig_carriers_only_trajectories`; BMI and age x BMI within group)

* **SRS / SCQ:** within-carrier BMI effects are modest (SRS all CNVs +0.22, p = 0.006; SCQ +0.33, p = 0.009) and come
  from 22q del and 16p dup; none survive in the non-ASD subsets. No age x BMI anywhere.
* **RBS-R:** nothing (all p > 0.1).
* **FSIQ:** no BMI main effect; 22q del age x BMI p = 0.0034 (n = 111; higher-BMI carriers decline with age in the
  plot), not in the non-ASD subset (p = 0.14).
* **CNB accuracy:** 16p dup lower accuracy with higher BMI (-0.15 SD/SD, p = 0.015; non-ASD -0.18, p = 0.0036);
  22q dup the opposite in non-ASD carriers (+0.13, p = 0.038). 22q del age x BMI p = 0.059 (0.035 non-ASD).
* **CNB speed:** no BMI main effects; scattered age x BMI in non-ASD subsets (16p dup p = 0.0035, 22q del 0.033).
* The "vs unaffected" rows for CNB in `unaffected_controls_tests.csv` are not meaningful (the 216 CNB controls have
  BMI for 2 individuals) and should be ignored.

## 6. Bottom line

With a clean non-ASD relative control group and no UKBB: (i) 16p11.2 deletion is the only locus whose head
circumference-BMI coupling exceeds that of unaffected relatives, and it does so in carriers without ASD too;
(ii) height shows no constant BMI shift at any locus, only age-dependent differences, largest in 16p11.2 dup and
16p11.2 del; (iii) social-communication measures (SRS, SCQ) increase with BMI in carriers but not in unaffected
relatives, mostly at 22q11.2 del, and attenuate when ASD carriers are removed; (iv) no robust BMI effects on RBS-R,
FSIQ or CNB beyond a 16p11.2 dup accuracy decrease.

Multiple testing: ~190 tests are tabulated; only the 16p del HC, the pooled and 22q del SRS/SCQ, and the height
three-way terms survive a Bonferroni-type threshold (p < 2.6e-4).

## 7. Addendum (2026-10-24, later): BMI slopes on SRS / SCQ without the age dimension

Script `24_cognition_bmi_slopes.R`; outputs `cognition_bmi_slopes.csv`, `fig_cognition_bmi_slopes.{png,pdf}`.
Because the SRS/SCQ effects are not age-dependent and the 70/30 tails are sparse, the trait is plotted against
continuous BMI Z (age-adjusted within group) with the fitted slope and binned means. The non-carrier row shows
unaffected SSC relatives and, as the ASD comparison, SSC non-carrier probands.

* SSC non-carrier probands (ASD, n = 2066 SRS / 1504 SCQ) have no BMI slope (+0.02, p = 0.4 and 0.7), the same as
  siblings and parents (slopes 0.00-0.02). The carrier slopes therefore are CNV-specific, not an ASD phenomenon.
* Carriers: SRS +0.22 SD/SD (p = 0.006), SCQ +0.33 (p = 0.009), monotonic across BMI bins. Per locus: 22q del
  (SRS +0.22, SCQ +0.29), 16p dup SRS (+0.53, n = 34), 22q dup SCQ (+0.44, in non-ASD carriers +0.52); 16p del is
  positive but not significant within carriers (SRS +0.18, SCQ +0.44) while differing from the unaffected slope.
* Confounders checked: sex (no change), FSIQ (no change; BMI and FSIQ uncorrelated in carriers), age group (same slope
  in children and adults), cohort (not applicable: all carriers with SRS/SCQ are G2MH). ASD status matters: carriers
  with ASD have higher BMI at every locus; adjusting for ASD halves the SRS slope (+0.12, p = 0.09) and leaves SCQ at
  +0.26 (p = 0.02). Non-ASD carriers alone: SRS +0.09, SCQ +0.15 (n.s.). Medication (antipsychotics) is unmeasured.

## 8. Medication sensitivity analysis (2026-10-24; `27_medication_sensitivity.R`)

Source: G2MH medications release (one row per drug; `meds_current` 1 = current). ATC codes are empty, so classes were assigned
from drug names (antipsychotic, stimulant/ADHD incl. atomoxetine, ADHD alpha-agonist, antidepressant, antiepileptic/mood
stabiliser, anxiolytic; weight-gain = antipsychotic + valproate + lithium + mirtazapine + gabapentinoids; weight-loss =
stimulants + atomoxetine + topiramate + metformin + GLP-1 agonists). Medication status is known for 87-92% of carriers with
SRS, SCQ or CNB. Outputs: `medication_prevalence.csv`, `medication_sensitivity.csv`.

Prevalence (current, all carriers): any psychotropic 29-33%, stimulants 13-16%, antidepressants 13-14%, antipsychotics
5-7% (22q11.2 del 9%, 16p11.2 dup SRS subset 17%), weight-gain medication 6-10%, weight-loss medication 15-17%.

Medication classes are strongly associated with the traits (antipsychotic: SRS +2.3 SD, SCQ +3.5 SD, CNB accuracy -0.6 SD;
any psychotropic: SRS +1.2, SCQ +1.3, CNB accuracy -0.2), but not with BMI in this sample (antipsychotic -0.2 SD, p = 0.5;
weight-gain medication -0.1, p = 0.7; weight-loss medication +0.3, p = 0.15; any psychotropic +0.3, p = 0.09). Without a
medication-BMI association there is little scope for medication to confound the BMI-trait slopes, and the slopes bear this out:

| Trait, all carriers | unadjusted (med-known subset) | + medication classes | no current psychotropic | + classes + ASD |
|---|---|---|---|---|
| SRS   | +0.20 (p = 0.019), n = 345 | +0.19 (0.024) | +0.23 (0.041), n = 231 | +0.08 (0.26) |
| SCQ   | +0.34 (0.011), n = 263 | +0.40 (0.003) | +0.38 (0.024), n = 177 | +0.33 (0.005) |
| CNB accuracy | -0.01 (0.81), n = 635 | 0.00 (0.86) | +0.01 (0.75) | -0.01 (0.59) |
| CNB accuracy, 16p11.2 dup | -0.10 (0.087), n = 52 | -0.10 (0.10) | -0.11 (0.086), n = 44 | -0.12 (0.055) |

Per locus the picture is the same: 22q11.2 del SRS +0.23 -> +0.24 with classes, +0.29 in unmedicated carriers; SCQ +0.36 ->
+0.40 / +0.42. The 16p11.2 dup SRS slope shrinks when the five carriers without medication data are dropped (+0.53 -> +0.40)
and further with adjustment (+0.28), but this is a 30-person group. The 16p11.2 dup CNB accuracy slope is unchanged by any
adjustment. ASD status, not medication, remains the only covariate that materially reduces the SRS slope.

Conclusion: psychotropic medication (including antipsychotics and other weight-affecting drugs) does not explain the BMI
effects on SRS and SCQ, nor the 16p11.2 duplication BMI effect on CNB accuracy.

## 9. Plausibility QC of growth measures (2026-10-24; `00_qc_growth.R`) and its effect on the results

A check of the 16p11.2 duplication CNB result exposed implausible raw values in the inputs, so a plausibility screen is now
applied in the shared data build (`21_pheno_bmi_gam.R`, hence 23-27 and the main figure v3) before any scoring: heights
recorded in metres (27 G2MH site-213 rows) are multiplied by 100; heights outside 60-230 cm (4; e.g. 602 cm) are set to
missing together with their BMI; BMI outside 10-70 (9; e.g. 9.2, 165, 181), |height Z| > 5 (22), child |BMI Z| > 5 (7)
and head circumference outside 30-70 cm or |Z| > 5 (17) are set to missing. The adult BMI reference and the G2MH
medical-history height/weight used for CNB BMI get the same bounds. Flagged records: `implausible_height_weight.{csv,md}`
(66 visit rows, 55 individuals; mostly SSC parents, Searchlight and the site-213 block; 1-6 carriers per CNV group).
Pre-QC outputs are kept in `preQC_20261024/`. Scripts not yet rerun with the QC: 17 (UKBB main figure, Expanse), 18, 19,
20, 22 and 10_supplement.

**What changed.** SRS and SCQ results are unchanged (all slopes within 0.03 SD/SD; e.g. all carriers vs unaffected SRS
+0.21, p = 1.6e-6; SCQ +0.34, p = 4.8e-5). The 16p11.2 deletion head-circumference CNV x BMI effect persists (+0.13,
p = 0.021; three-way p = 0.020) with its 16p-specificity (16p dup, 22q del, 22q dup CNV x BMI p = 0.39, 0.73, 0.66). The
16p11.2 duplication CNB accuracy slope is -0.14 (p = 0.074) and the across-CNV heterogeneity for accuracy p = 0.11
(was 0.027); no CNB measure now shows significant heterogeneity. The pooled undergrowth height CNV x BMI shift
(-0.09, p = 0.019) is gone (-0.05, p = 0.15).

**A caution on the CNV x age x BMI terms.** Removing a handful of extreme height values changed several height three-way
p-values by orders of magnitude in both directions (16p del 4.4e-6 -> 0.10; 16p dup < 1e-16 -> 0.012; 22q del 0.013 ->
0.99; 22q dup 0.063 -> 0.0055; HC 22q del 0.80 -> 0.045; HC 22q dup 0.998 -> 0.0014). These tensor-product
interactions in `bam` with a few hundred carriers against ~5700 controls were already known to be sensitive to row
order; they are evidently also sensitive to single extreme observations and should not be headlined. The CNV x BMI
(constant shift) tests and the SRS/SCQ slopes were stable.

**Longitudinal (prospective) analysis after QC** (`01_build.R` -> `05_age_varying.R` -> `07_hc_trajectories.R` -> `11_supp_prospective.R`
-> `29_supp_table_S3_prospective.R`; `01_build.R` now also replaces SSC parent snapshot ages by visit age, 4300 rows). The height
results are identical before and after QC: none of the flagged records survived 05's existing longitudinal cleaning (age-varying
p 16p11.2 del 0.006, pooled 0.022; < 8 y window 16p11.2 del +0.23 SD per SD, t = 3.9, p = 1.6e-4). Head circumference lost 5 SVIP
intervals and the all-ages coefficient rose: 16p11.2 del +0.04 -> +0.08 (t = 2.1, p = 0.038), pooled +0.07 -> +0.10 (t = 3.4,
p = 0.0007; age-varying p 0.0039 -> 0.0007). The SVIP v10 folder had moved (`G2MH/Simons_Searchlight/SVIP/SVIP/...`); paths in 07
and 09 are updated, and `lmerTest` was installed for R 4.6.1. Effects, SE, t and p are in Table S3 (`supp_table_S3_prospective.csv`,
sheet "S3 prospective" of `supp_tables_bmi_effects.xlsx`).
