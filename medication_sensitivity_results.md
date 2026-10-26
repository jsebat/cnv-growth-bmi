# Medication sensitivity analysis: BMI effects on SRS, SCQ and Penn CNB in CNV carriers

(Updated 2026-10-24 after plausibility QC of height, weight and head circumference; see the QC note at the end. Pre-QC outputs are kept in `preQC_20261024/`.)

## Methods

Height, weight and head circumference were screened for implausible values before scoring (heights recorded in metres corrected; heights outside 60-230 cm, BMI outside 10-70 kg/m2, height or child-BMI Z beyond +/-5, and head circumference outside 30-70 cm or Z beyond +/-5 set to missing; `00_qc_growth.R`, flagged records in `implausible_height_weight.csv`). Current medication use was taken from the G2MH medications release (February 2026). ATC codes were not populated, so drugs were classed by name as antipsychotics, ADHD medication (stimulants, atomoxetine), ADHD alpha-agonists, antidepressants, antiepileptics/mood stabilisers, and anxiolytics; "any psychotropic" is their union. Status was known for 87-92% of carriers with SRS, SCQ or CNB; carriers absent from the medications file were excluded from adjusted models. For each trait and CNV group the age-adjusted BMI slope (GAM, trait ~ s(age) + BMI; locus as a covariate when pooled) was refit in this subset unadjusted, with medication-class indicators, in carriers on no psychotropic, and with medication classes plus ASD diagnosis. Associations of each class with BMI and with the trait were tested in the same models. Slopes are SD of trait per SD of BMI.

## Results

Psychotropic use was common (29-33% of carriers; antipsychotics 5-7%, 9% in 22q11.2 deletion) and strongly associated with the traits (antipsychotic users: SRS +2.3 SD, SCQ +3.5 SD, CNB accuracy -0.6 SD; all p < 1e-4), but no class was associated with BMI (antipsychotics -0.2 SD, p = 0.5; weight-gain-associated drugs -0.1, p = 0.7; any psychotropic +0.3, p = 0.09), leaving little scope for confounding. The BMI slopes were unchanged by adjustment or by excluding medicated carriers (Table): across all carriers SRS +0.22, +0.19 and +0.23; SCQ +0.35, +0.37 and +0.38; CNB accuracy 0.00 throughout. The 22q11.2 deletion SRS and SCQ slopes persisted in unmedicated carriers. The negative BMI slope on CNB accuracy in 16p11.2 duplication carriers, -0.15 (p = 0.015) before QC, fell to -0.14 (p = 0.07) once one carrier with a recorded BMI of 9.2 kg/m2 was removed, and to -0.05 to -0.10 under medication adjustment or in unmedicated carriers; it should be regarded as suggestive. ASD diagnosis, not medication, remained the only covariate that attenuated the SRS slope (+0.09), while SCQ was robust to all adjustments (+0.32, p = 0.008). Medication therefore does not account for the BMI effects on social-communication measures.

## Table. Age-adjusted BMI slope (SD of trait per SD of BMI; estimate (SE), p) by model, after QC

n = carriers with medication information; % psychotropic / % antipsychotic = current use. "No psychotropic" is the slope in carriers on no current psychotropic (n in brackets; omitted where n < 20). "+ classes + ASD" adds ASD diagnosis to the medication-adjusted model. Source: `27_medication_sensitivity.R`; full output in `medication_sensitivity.csv` and `medication_prevalence.csv`.

| Trait | Group | n | % psychotropic | % antipsychotic | Unadjusted | + medication classes | No psychotropic (n) | + classes + ASD |
|---|---|---|---|---|---|---|---|---|
| SRS | 16p11.2 del | 69 | 22 | 3 | +0.09 (0.22), p = 0.68 | +0.09 (0.23), p = 0.69 | +0.05 (0.28), p = 0.85 (54) | -0.01 (0.17), p = 0.96 |
| SRS | 16p11.2 dup | 30 | 40 | 17 | +0.40 (0.32), p = 0.22 | +0.28 (0.34), p = 0.41 | n/a | +0.22 (0.27), p = 0.41 |
| SRS | 22q11.2 del | 184 | 37 | 9 | +0.26 (0.11), p = 0.01 | +0.23 (0.10), p = 0.03 | +0.29 (0.13), p = 0.02 (116) | +0.12 (0.10), p = 0.23 |
| SRS | 22q11.2 dup | 61 | 30 | 2 | +0.08 (0.24), p = 0.73 | +0.10 (0.27), p = 0.71 | +0.05 (0.34), p = 0.88 (43) | -0.01 (0.21), p = 0.97 |
| SRS | All CNVs | 344 | 33 | 7 | +0.22 (0.09), p = 0.01 | +0.19 (0.09), p = 0.03 | +0.23 (0.11), p = 0.04 (231) | +0.09 (0.08), p = 0.23 |
| SCQ | 16p11.2 del | 52 | 15 | 2 | +0.45 (0.35), p = 0.21 | +0.26 (0.35), p = 0.45 | +0.13 (0.36), p = 0.73 (44) | +0.32 (0.26), p = 0.23 |
| SCQ | 22q11.2 del | 147 | 38 | 9 | +0.38 (0.17), p = 0.02 | +0.37 (0.16), p = 0.02 | +0.42 (0.19), p = 0.03 (91) | +0.25 (0.15), p = 0.09 |
| SCQ | 22q11.2 dup | 45 | 36 | 0 | +0.37 (0.37), p = 0.33 | +0.54 (0.43), p = 0.22 | +0.82 (0.55), p = 0.15 (29) | +0.66 (0.38), p = 0.09 |
| SCQ | All CNVs | 262 | 32 | 7 | +0.35 (0.14), p = 0.01 | +0.37 (0.14), p = 0.008 | +0.38 (0.17), p = 0.02 (177) | +0.32 (0.12), p = 0.008 |
| CNB accuracy | 16p11.2 del | 133 | 20 | 2 | +0.03 (0.04), p = 0.43 | +0.04 (0.04), p = 0.36 | +0.04 (0.05), p = 0.41 (107) | +0.01 (0.04), p = 0.74 |
| CNB accuracy | 16p11.2 dup | 51 | 16 | 0 | -0.05 (0.08), p = 0.51 | -0.05 (0.08), p = 0.57 | -0.06 (0.09), p = 0.50 (43) | -0.10 (0.08), p = 0.24 |
| CNB accuracy | 22q11.2 del | 318 | 36 | 9 | -0.02 (0.03), p = 0.53 | -0.02 (0.03), p = 0.51 | +0.03 (0.04), p = 0.39 (205) | -0.01 (0.03), p = 0.74 |
| CNB accuracy | 22q11.2 dup | 130 | 28 | 2 | +0.07 (0.06), p = 0.20 | +0.09 (0.06), p = 0.11 | +0.04 (0.07), p = 0.53 (94) | +0.08 (0.06), p = 0.17 |
| CNB accuracy | All CNVs | 632 | 29 | 5 | +0.00 (0.02), p = 0.85 | +0.01 (0.02), p = 0.78 | +0.03 (0.03), p = 0.33 (449) | -0.00 (0.02), p = 0.93 |
| CNB speed | 16p11.2 del | 135 | 20 | 2 | -0.02 (0.05), p = 0.73 | -0.01 (0.05), p = 0.84 | +0.01 (0.05), p = 0.83 (108) | +0.01 (0.05), p = 0.83 |
| CNB speed | 16p11.2 dup | 51 | 16 | 0 | +0.01 (0.11), p = 0.96 | -0.01 (0.11), p = 0.96 | -0.03 (0.12), p = 0.83 (43) | -0.02 (0.12), p = 0.84 |
| CNB speed | 22q11.2 del | 324 | 35 | 9 | -0.00 (0.03), p = 0.98 | -0.01 (0.03), p = 0.78 | -0.01 (0.04), p = 0.87 (210) | -0.00 (0.04), p = 0.96 |
| CNB speed | 22q11.2 dup | 130 | 28 | 2 | +0.05 (0.09), p = 0.57 | +0.06 (0.09), p = 0.51 | +0.10 (0.08), p = 0.19 (94) | +0.05 (0.09), p = 0.60 |
| CNB speed | All CNVs | 640 | 29 | 5 | +0.00 (0.03), p = 0.97 | -0.00 (0.03), p = 0.98 | +0.00 (0.03), p = 0.87 (455) | -0.00 (0.03), p = 0.99 |

## Addendum: the 16p11.2 duplication CNB accuracy slope and the plausibility QC

Before QC the 16p11.2 duplication slope on CNB accuracy was -0.15 (p = 0.015, n = 57) and the heterogeneity of the BMI slope across the four CNV groups was p = 0.027. The slope was not explained by age (present in children, -0.18, p = 0.009; absent in adults; no age x BMI interaction), sex, ASD, site, parental or own education, number of valid subtests, timing of the BMI measurement, or medication (-0.10 to -0.12 under all adjustments). It was, however, carried by one 10-year-old with a recorded BMI of 9.2 kg/m2 (BMI Z -9.6; from the G2MH medical-history height and weight) and above-average accuracy. That record is removed by the plausibility QC now applied to all analyses (BMI outside 10-70 set to missing). After QC: slope -0.14 (SE 0.08, p = 0.074, n = 56); with medication classes -0.05 (p = 0.57); in unmedicated carriers -0.06 (p = 0.50); with classes + ASD -0.10 (p = 0.24); heterogeneity across CNVs F = 2.02, p = 0.11. All ten accuracy subtests still slope negative in this group, so the direction is consistent, but the result is suggestive only.

The QC changed little else in the medication analysis (SRS and SCQ slopes within 0.03 of their pre-QC values). Its effects on the other analyses are summarised in section 9 of `unaffected_controls_report.md`. The screening rules and the 66 flagged records (55 individuals) are in `00_qc_growth.R` and `implausible_height_weight.{csv,md}`.
