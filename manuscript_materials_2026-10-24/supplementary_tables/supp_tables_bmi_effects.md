# Supplementary tables: BMI effects on growth and social-communication / cognitive measures in CNV carriers

- Table S1: height and head circumference (WHO / UK90 Z), clinical cohorts, ages 2-65, latest visit with the outcome and BMI (plausibility QC applied, 00_qc_growth.R). Groups: 16p11.2 proximal deletion / duplication, 22q11.2 A-D deletion / duplication; pooled overgrowth = 16p11.2 del + 22q11.2 dup, undergrowth = 16p11.2 dup + 22q11.2 del, all CNVs (locus as a covariate). Reference groups: unaffected controls = SSC designated siblings and parents (ages at the study visit) plus SVIP non-carrier family members without ASD; SSC ASD probands = SSC probands with no CNV.
- Table S2: SRS total (Z vs SSC unaffected relatives), SCQ lifetime (Z vs SSC siblings), Penn CNB accuracy and speed (Z vs G2MH non-carriers; speed = -RT z). Reference groups as in S1 (SCQ: siblings only; CNB: the G2MH non-carriers have no BMI, so no reference contrast).
- BMI: WHO child BMI Z; adults, centile vs pooled non-carrier adults mapped to Z. All tests use continuous BMI.
- BMI slope: within-group age-adjusted BMI effect, SD of trait per SD of BMI (GAM y ~ BMI + s(age) [+ locus]). age x BMI p: tensor-product age x BMI interaction within the group (GAM y ~ s(age) + s(BMI) + ti(age, BMI) [+ locus]).
- CNV main effect: age-adjusted mean difference carriers - reference (SD) from the CNV x BMI model. CNV x BMI: difference between carriers and the reference group in the BMI slope, adjusted for age (bam y ~ [locus +] CNV + s(age, by CNV) + s(age, by BMI) + CNV:BMI). CNV x age x BMI: whether that difference varies with age (bam y ~ [locus +] CNV + s(age) + s(BMI) + s(age, by CNV) + ti(age, BMI) + ti(age, BMI, by CNV); p of the carrier-specific tensor product). The three-way terms are sensitive to extreme observations (see unaffected_controls_report.md section 9).
- Reference-group rows give the within-group BMI slope and age x BMI for that reference group itself.
- Generated 2026-10-24 by 28_supp_tables.R.

## Table S1. Height and head circumference

| Trait | Group | Members | Reference group | n | n reference | BMI slope (SE) | BMI slope p | age x BMI p | CNV main effect (SD) | CNV main p | CNV x BMI (SE) | CNV x BMI p | CNV x age x BMI p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Height | Unaffected controls | Unaffected controls |  | 5708 |  | +0.038 (0.013) | 0.0043 | 2.7e-07 |  |  |  |  |  |
| Height | SSC ASD probands | SSC ASD probands |  | 2074 |  | +0.127 (0.018) | 6.6e-12 | 5.2e-04 |  |  |  |  |  |
| Height | 16p11.2 del | 16p11.2 del | SSC ASD probands |  358 | 2074 | +0.105 (0.048) | 0.029 | 0.0027 | -0.37 | 2.8e-05 | +0.009 (0.044) | 0.84 | 0.87 |
| Height | 16p11.2 del | 16p11.2 del | Unaffected controls |  358 | 5708 | +0.105 (0.048) | 0.029 | 0.0027 | -0.70 | 1.3e-05 | +0.000 (0.042) | 0.99 | 0.1 |
| Height | 16p11.2 dup | 16p11.2 dup | SSC ASD probands |  199 | 2074 | +0.018 (0.070) | 0.8 | 0.54 | -0.49 | 2.8e-06 | -0.115 (0.068) | 0.091 | 0.068 |
| Height | 16p11.2 dup | 16p11.2 dup | Unaffected controls |  199 | 5708 | +0.018 (0.070) | 0.8 | 0.54 | +0.11 | 0.34 | -0.068 (0.058) | 0.24 | 0.012 |
| Height | 22q11.2 del | 22q11.2 del | SSC ASD probands |  357 | 2074 | +0.080 (0.044) | 0.069 | 0.16 | -0.94 | < 1e-16 | -0.022 (0.050) | 0.66 | 0.0025 |
| Height | 22q11.2 del | 22q11.2 del | Unaffected controls |  357 | 5708 | +0.080 (0.044) | 0.069 | 0.16 | -0.61 | 6.8e-12 | -0.031 (0.044) | 0.48 | 0.99 |
| Height | 22q11.2 dup | 22q11.2 dup | SSC ASD probands |  118 | 2074 | +0.106 (0.074) | 0.15 | 0.051 | +0.21 | 0.18 | +0.014 (0.090) | 0.88 | 0.0015 |
| Height | 22q11.2 dup | 22q11.2 dup | Unaffected controls |  118 | 5708 | +0.106 (0.074) | 0.15 | 0.051 | -0.32 | 0.012 | -0.001 (0.076) | 0.99 | 0.0055 |
| Height | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | SSC ASD probands |  476 | 2074 | +0.100 (0.040) | 0.014 | 9.0e-05 | -0.37 | 1.3e-05 | +0.008 (0.041) | 0.85 | 0.79 |
| Height | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | Unaffected controls |  476 | 5708 | +0.100 (0.040) | 0.014 | 9.0e-05 | -0.80 | 2.8e-11 | -0.000 (0.039) | 0.99 | 0.035 |
| Height | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | SSC ASD probands |  556 | 2074 | +0.055 (0.037) | 0.14 | 0.1 | -0.40 | 1.2e-04 | -0.052 (0.043) | 0.22 | 0.26 |
| Height | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | Unaffected controls |  556 | 5708 | +0.055 (0.037) | 0.14 | 0.1 | -0.01 | 0.88 | -0.054 (0.037) | 0.15 | 1.4e-04 |
| Height | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | SSC ASD probands | 1032 | 2074 | +0.069 (0.028) | 0.014 | < 1e-16 | -0.38 | 2.3e-06 | -0.007 (0.033) | 0.83 | 0.5 |
| Height | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | Unaffected controls | 1032 | 5708 | +0.069 (0.028) | 0.014 | < 1e-16 | -0.26 | 0.003 | -0.045 (0.030) | 0.13 | 0.066 |
| Head circumference | Unaffected controls | Unaffected controls |  | 5680 |  | +0.364 (0.014) | < 1e-16 | 0.0098 |  |  |  |  |  |
| Head circumference | SSC ASD probands | SSC ASD probands |  | 2060 |  | +0.295 (0.019) | < 1e-16 | 8.7e-05 |  |  |  |  |  |
| Head circumference | 16p11.2 del | 16p11.2 del | SSC ASD probands |  211 | 2060 | +0.474 (0.056) | 4.5e-15 | 0.082 | +0.02 | 0.88 | +0.172 (0.056) | 0.002 | 0.84 |
| Head circumference | 16p11.2 del | 16p11.2 del | Unaffected controls |  211 | 5680 | +0.474 (0.056) | 4.5e-15 | 0.082 | -0.05 | 0.79 | +0.131 (0.057) | 0.021 | 0.02 |
| Head circumference | 16p11.2 dup | 16p11.2 dup | SSC ASD probands |  125 | 2060 | +0.290 (0.078) | 2.9e-04 | 0.77 | -1.32 | < 1e-16 | -0.056 (0.089) | 0.53 | 0.11 |
| Head circumference | 16p11.2 dup | 16p11.2 dup | Unaffected controls |  125 | 5680 | +0.290 (0.078) | 2.9e-04 | 0.77 | -0.88 | 1.4e-09 | -0.069 (0.081) | 0.39 | 0.68 |
| Head circumference | 22q11.2 del | 22q11.2 del | SSC ASD probands |  242 | 2060 | +0.345 (0.055) | 2.2e-09 | 0.73 | -0.74 | 5.3e-11 | -0.011 (0.057) | 0.84 | 0.1 |
| Head circumference | 22q11.2 del | 22q11.2 del | Unaffected controls |  242 | 5680 | +0.345 (0.055) | 2.2e-09 | 0.73 | -0.80 | 1.8e-06 | -0.019 (0.055) | 0.73 | 0.045 |
| Head circumference | 22q11.2 dup | 22q11.2 dup | SSC ASD probands |   92 | 2060 | +0.417 (0.095) | 3.0e-05 | 0.13 | -0.62 | 0.0024 | +0.113 (0.097) | 0.25 | 0.08 |
| Head circumference | 22q11.2 dup | 22q11.2 dup | Unaffected controls |   92 | 5680 | +0.417 (0.095) | 3.0e-05 | 0.13 | -0.39 | 0.032 | +0.040 (0.089) | 0.66 | 0.0014 |
| Head circumference | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | SSC ASD probands |  303 | 2060 | +0.450 (0.048) | < 1e-16 | 0.13 | +0.01 | 0.93 | +0.155 (0.050) | 0.0019 | 0.066 |
| Head circumference | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | Unaffected controls |  303 | 5680 | +0.450 (0.048) | < 1e-16 | 0.13 | +0.04 | 0.78 | +0.094 (0.049) | 0.056 | 0.017 |
| Head circumference | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | SSC ASD probands |  367 | 2060 | +0.324 (0.045) | 3.6e-12 | 0.79 | -1.24 | < 1e-16 | -0.012 (0.050) | 0.81 | 0.11 |
| Head circumference | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | Unaffected controls |  367 | 5680 | +0.324 (0.045) | 3.6e-12 | 0.79 | -1.01 | 1.1e-15 | -0.040 (0.047) | 0.4 | 0.04 |
| Head circumference | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | SSC ASD probands |  670 | 2060 | +0.382 (0.033) | < 1e-16 | 0.075 | +0.07 | 0.46 | +0.077 (0.038) | 0.044 | 6.5e-05 |
| Head circumference | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | Unaffected controls |  670 | 5680 | +0.382 (0.033) | < 1e-16 | 0.075 | +0.24 | 0.037 | +0.018 (0.037) | 0.63 | 0.018 |

## Table S2. SRS, SCQ and Penn CNB

| Trait | Group | Members | Reference group | n | n reference | BMI slope (SE) | BMI slope p | age x BMI p | CNV main effect (SD) | CNV main p | CNV x BMI (SE) | CNV x BMI p | CNV x age x BMI p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SRS | Unaffected controls | Unaffected controls |  | 5486 |  | +0.013 (0.013) | 0.3 | 0.11 |  |  |  |  |  |
| SRS | SSC ASD probands | SSC ASD probands |  | 2065 |  | +0.025 (0.031) | 0.41 | 0.083 |  |  |  |  |  |
| SRS | 16p11.2 del | 16p11.2 del | SSC ASD probands |   79 | 2065 | +0.180 (0.202) | 0.38 | 0.77 | -1.68 | 1.4e-06 | +0.137 (0.151) | 0.36 | 0.68 |
| SRS | 16p11.2 del | 16p11.2 del | Unaffected controls |   79 | 5486 | +0.180 (0.202) | 0.38 | 0.77 | +1.73 | 4.6e-08 | +0.170 (0.085) | 0.044 | 0.24 |
| SRS | 16p11.2 dup | 16p11.2 dup | SSC ASD probands |   34 | 2065 | +0.528 (0.253) | 0.046 | 0.71 | -0.75 | 0.071 | +0.383 (0.236) | 0.1 | 1 |
| SRS | 16p11.2 dup | 16p11.2 dup | Unaffected controls |   34 | 5486 | +0.528 (0.253) | 0.046 | 0.71 | +2.45 | < 1e-16 | +0.502 (0.122) | 3.8e-05 | 0.28 |
| SRS | 22q11.2 del | 22q11.2 del | SSC ASD probands |  208 | 2065 | +0.254 (0.098) | 0.011 | 0.47 | -2.48 | < 1e-16 | +0.193 (0.098) | 0.05 | 0.5 |
| SRS | 22q11.2 del | 22q11.2 del | Unaffected controls |  208 | 5486 | +0.254 (0.098) | 0.011 | 0.47 | +2.14 | < 1e-16 | +0.232 (0.052) | 8.4e-06 | 0.7 |
| SRS | 22q11.2 dup | 22q11.2 dup | SSC ASD probands |   67 | 2065 | +0.073 (0.234) | 0.76 | 0.9 | -2.28 | 9.1e-10 | -0.036 (0.193) | 0.85 | 1 |
| SRS | 22q11.2 dup | 22q11.2 dup | Unaffected controls |   67 | 5486 | +0.073 (0.234) | 0.76 | 0.9 | +1.34 | 2.1e-06 | -0.047 (0.106) | 0.66 | 0.083 |
| SRS | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | SSC ASD probands |  146 | 2065 | +0.144 (0.151) | 0.34 | 1 | -1.72 | 7.3e-08 | +0.063 (0.123) | 0.61 | 1 |
| SRS | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | Unaffected controls |  146 | 5486 | +0.144 (0.151) | 0.34 | 1 | +1.92 | < 1e-16 | +0.075 (0.068) | 0.27 | 0.052 |
| SRS | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | SSC ASD probands |  242 | 2065 | +0.288 (0.093) | 0.0021 | 0.38 | -1.43 | 6.5e-05 | +0.200 (0.093) | 0.032 | 0.5 |
| SRS | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | Unaffected controls |  242 | 5486 | +0.288 (0.093) | 0.0021 | 0.38 | +2.94 | < 1e-16 | +0.268 (0.049) | 5.8e-08 | 0.74 |
| SRS | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | SSC ASD probands |  388 | 2065 | +0.235 (0.080) | 0.0036 | 0.6 | -2.12 | 1.2e-14 | +0.167 (0.080) | 0.036 | 1 |
| SRS | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | Unaffected controls |  388 | 5486 | +0.235 (0.080) | 0.0036 | 0.6 | +2.14 | < 1e-16 | +0.212 (0.044) | 1.6e-06 | 0.21 |
| SCQ | Unaffected controls | Unaffected controls |  | 1134 |  | +0.009 (0.023) | 0.71 | 0.72 |  |  |  |  |  |
| SCQ | SSC ASD probands | SSC ASD probands |  | 1503 |  | +0.026 (0.067) | 0.7 | 0.65 |  |  |  |  |  |
| SCQ | 16p11.2 del | 16p11.2 del | SSC ASD probands |   60 | 1503 | +0.440 (0.314) | 0.17 | 0.5 | -2.95 | 3.5e-05 | +0.400 (0.312) | 0.2 | 1 |
| SCQ | 16p11.2 del | 16p11.2 del | Unaffected controls |   60 | 1134 | +0.440 (0.314) | 0.17 | 0.5 | +7.07 | < 1e-16 | +0.261 (0.123) | 0.034 | 1 |
| SCQ | 16p11.2 dup | 16p11.2 dup | SSC ASD probands |   21 | 1503 | +0.351 (0.493) | 0.49 | 0.5 | -4.17 | 8.2e-06 | +0.283 (0.490) | 0.56 | 0.5 |
| SCQ | 16p11.2 dup | 16p11.2 dup | Unaffected controls |   21 | 1134 | +0.351 (0.493) | 0.49 | 0.5 | +9.14 | < 1e-16 | +0.478 (0.165) | 0.0039 | 1 |
| SCQ | 22q11.2 del | 22q11.2 del | SSC ASD probands |  170 | 1503 | +0.293 (0.152) | 0.056 | 0.47 | -5.76 | < 1e-16 | +0.241 (0.203) | 0.23 | 0.33 |
| SCQ | 22q11.2 del | 22q11.2 del | Unaffected controls |  170 | 1134 | +0.293 (0.152) | 0.056 | 0.47 | +3.86 | < 1e-16 | +0.274 (0.083) | 0.0011 | 0.052 |
| SCQ | 22q11.2 dup | 22q11.2 dup | SSC ASD probands |   50 | 1503 | +0.436 (0.370) | 0.24 | 0.58 | -5.17 | 5.8e-14 | +0.362 (0.367) | 0.32 | 0.38 |
| SCQ | 22q11.2 dup | 22q11.2 dup | Unaffected controls |   50 | 1134 | +0.436 (0.370) | 0.24 | 0.58 | +4.29 | < 1e-16 | +0.392 (0.128) | 0.0023 | 1 |
| SCQ | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | SSC ASD probands |  110 | 1503 | +0.440 (0.238) | 0.067 | 0.95 | -3.18 | 2.6e-07 | +0.388 (0.246) | 0.11 | 1 |
| SCQ | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | Unaffected controls |  110 | 1134 | +0.440 (0.238) | 0.067 | 0.95 | +6.32 | < 1e-16 | +0.418 (0.101) | 3.7e-05 | 0.47 |
| SCQ | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | SSC ASD probands |  191 | 1503 | +0.291 (0.144) | 0.045 | 0.5 | -4.86 | 3.8e-09 | +0.242 (0.192) | 0.21 | 0.59 |
| SCQ | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | Unaffected controls |  191 | 1134 | +0.291 (0.144) | 0.045 | 0.5 | +5.87 | < 1e-16 | +0.285 (0.081) | 4.5e-04 | 0.15 |
| SCQ | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | SSC ASD probands |  301 | 1503 | +0.333 (0.127) | 0.0095 | 0.49 | -3.63 | 1.6e-09 | +0.316 (0.161) | 0.049 | 1 |
| SCQ | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | Unaffected controls |  301 | 1134 | +0.333 (0.127) | 0.0095 | 0.49 | +6.06 | < 1e-16 | +0.338 (0.083) | 4.8e-05 | 0.29 |
| CNB accuracy | 16p11.2 del | 16p11.2 del | none (no reference group with BMI) |  140 |  | +0.024 (0.041) | 0.57 | 0.39 |  |  |  |  |  |
| CNB accuracy | 16p11.2 dup | 16p11.2 dup | none (no reference group with BMI) |   56 |  | -0.138 (0.076) | 0.074 | 0.48 |  |  |  |  |  |
| CNB accuracy | 22q11.2 del | 22q11.2 del | none (no reference group with BMI) |  347 |  | -0.018 (0.028) | 0.53 | 0.066 |  |  |  |  |  |
| CNB accuracy | 22q11.2 dup | 22q11.2 dup | none (no reference group with BMI) |  142 |  | +0.061 (0.051) | 0.24 | 1 |  |  |  |  |  |
| CNB accuracy | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | none (no reference group with BMI) |  282 |  | +0.042 (0.033) | 0.21 | 0.62 |  |  |  |  |  |
| CNB accuracy | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | none (no reference group with BMI) |  403 |  | -0.030 (0.027) | 0.26 | 0.1 |  |  |  |  |  |
| CNB accuracy | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | none (no reference group with BMI) |  685 |  | -0.004 (0.021) | 0.85 | 0.11 |  |  |  |  |  |
| CNB speed | 16p11.2 del | 16p11.2 del | none (no reference group with BMI) |  142 |  | -0.018 (0.045) | 0.69 | 0.54 |  |  |  |  |  |
| CNB speed | 16p11.2 dup | 16p11.2 dup | none (no reference group with BMI) |   56 |  | -0.035 (0.100) | 0.73 | 0.13 |  |  |  |  |  |
| CNB speed | 22q11.2 del | 22q11.2 del | none (no reference group with BMI) |  353 |  | -0.003 (0.033) | 0.92 | 0.1 |  |  |  |  |  |
| CNB speed | 22q11.2 dup | 22q11.2 dup | none (no reference group with BMI) |  142 |  | +0.026 (0.074) | 0.72 | 0.95 |  |  |  |  |  |
| CNB speed | Overgrowth CNVs | 16p11.2 del + 22q11.2 dup | none (no reference group with BMI) |  284 |  | +0.010 (0.044) | 0.82 | 0.92 |  |  |  |  |  |
| CNB speed | Undergrowth CNVs | 16p11.2 dup + 22q11.2 del | none (no reference group with BMI) |  409 |  | -0.004 (0.031) | 0.9 | 0.53 |  |  |  |  |  |
| CNB speed | All CNVs | 16p11.2 del + 16p11.2 dup + 22q11.2 del + 22q11.2 dup | none (no reference group with BMI) |  693 |  | -0.003 (0.026) | 0.9 | 0.2 |  |  |  |  |  |

- Table S3: prospective effect of BMI on subsequent growth in CNV carriers < 18 y (fig_supp_prospective). Outcome: change in height Z (WHO) or head-circumference Z (UK90) from one visit to the next, per 1 SD of BMI Z at the first visit, adjusted for Z at the visit, interval length, sex, locus (pooled) and cohort, with a child random intercept. Height: cleaned longitudinal visits of SVIP, Searchlight and Cardiff (intervals >= 0.5 y; 05_age_varying.R). Head circumference: full SVIP v10 series + Cardiff, G2MH, SSC (07_hc_trajectories.R). Plausibility QC (00_qc_growth.R) applied to the genomatrix and the SVIP v10 series. 22q11.2 dup has too few repeat measures and is not shown in the figure.
- Age-varying p: mixed GAM with s(age, by = BMI); p that the BMI effect differs from 0 at any age. All ages: the same model with a single constant BMI coefficient (effect, SE, t = effect / SE, p). Windows: linear mixed model (lmer) with BMI x age-window terms for height (< 8, 8-13, 13-18 y) and GAM-based windows for head circumference (2-5, 5-18 y); t = estimate / SE.


| Outcome | Group | Term | Children | Intervals | Effect per SD of BMI Z (SE) | t | p | Age-varying p |
|---|---|---|---|---|---|---|---|---|
| Height | 16p11.2 del | all ages (constant BMI coefficient) | 116 | 224 | +0.086 (0.041) | 2.12 | 0.04 | 0.006 |
| Height | 16p11.2 del | window 13-18 y |  28 |  41 | -0.007 (0.083) | -0.08 | 0.93 |  |
| Height | 16p11.2 del | window 8-13 y |  60 |  90 | +0.003 (0.063) | 0.05 | 0.96 |  |
| Height | 16p11.2 del | window <8 y |  61 |  93 | +0.234 (0.061) | 3.86 | 1.6e-04 |  |
| Height | 16p11.2 dup | all ages (constant BMI coefficient) |  48 |  84 | +0.003 (0.064) | 0.05 | 0.96 | 0.11 |
| Height | 16p11.2 dup | window 13-18 y |  12 |  24 | -0.202 (0.106) | -1.91 | 0.06 |  |
| Height | 16p11.2 dup | window 8-13 y |  25 |  40 | +0.107 (0.089) | 1.20 | 0.23 |  |
| Height | 16p11.2 dup | window <8 y |  20 |  20 | +0.051 (0.120) | 0.42 | 0.68 |  |
| Height | 22q11.2 del | all ages (constant BMI coefficient) |  73 | 102 | -0.038 (0.046) | -0.82 | 0.41 | 0.60 |
| Height | 22q11.2 del | window 13-18 y |  22 |  23 | -0.026 (0.119) | -0.22 | 0.83 |  |
| Height | 22q11.2 del | window 8-13 y |  52 |  65 | -0.057 (0.060) | -0.94 | 0.35 |  |
| Height | 22q11.2 del | window <8 y |  14 |  14 | +0.026 (0.097) | 0.26 | 0.79 |  |
| Height | Pooled | all ages (constant BMI coefficient) | 248 | 422 | +0.026 (0.026) | 1.00 | 0.32 | 0.02 |
| Height | Pooled | window 13-18 y |  63 |  89 | -0.066 (0.050) | -1.31 | 0.19 |  |
| Height | Pooled | window 8-13 y | 143 | 201 | +0.007 (0.036) | 0.20 | 0.84 |  |
| Height | Pooled | window <8 y | 100 | 132 | +0.127 (0.044) | 2.86 | 0.005 |  |
| Head circumference | 16p11.2 del | all ages (constant BMI coefficient) |  91 | 290 | +0.037 (0.034) | 1.08 | 0.28 | 0.49 |
| Head circumference | 16p11.2 del | window 2-5 y |  31 |  44 | +0.075 (0.088) | 0.85 | 0.40 |  |
| Head circumference | 16p11.2 del | window 5-18 y |  35 |  42 | +0.053 (0.086) | 0.61 | 0.54 |  |
| Head circumference | 16p11.2 del | window <2 y |  55 | 204 | +0.032 (0.040) | 0.79 | 0.43 |  |
| Head circumference | 16p11.2 dup | all ages (constant BMI coefficient) |  46 | 161 | +0.058 (0.049) | 1.17 | 0.24 | 0.29 |
| Head circumference | 16p11.2 dup | window 2-5 y |  25 |  37 | +0.098 (0.096) | 1.02 | 0.31 |  |
| Head circumference | 16p11.2 dup | window 5-18 y |  10 |  11 | -0.085 (0.206) | -0.41 | 0.68 |  |
| Head circumference | 16p11.2 dup | window <2 y |  30 | 113 | +0.061 (0.068) | 0.90 | 0.37 |  |
| Head circumference | 22q11.2 del | all ages (constant BMI coefficient) |  50 |  56 | +0.193 (0.111) | 1.75 | 0.09 | 0.23 |
| Head circumference | Pooled | all ages (constant BMI coefficient) | 193 | 514 | +0.069 (0.026) | 2.64 | 0.009 | 0.004 |
| Head circumference | Pooled | window 2-5 y |  56 |  81 | +0.091 (0.067) | 1.35 | 0.18 |  |
| Head circumference | Pooled | window 5-18 y | 101 | 116 | +0.143 (0.052) | 2.77 | 0.006 |  |
| Head circumference | Pooled | window <2 y |  85 | 317 | +0.032 (0.035) | 0.94 | 0.35 |  |

- Table S3: prospective effect of BMI on subsequent growth in CNV carriers < 18 y (fig_supp_prospective). Outcome: change in height Z (WHO) or head-circumference Z (UK90) from one visit to the next, per 1 SD of BMI Z at the first visit, adjusted for Z at the visit, interval length, sex, locus (pooled) and cohort, with a child random intercept. Height: cleaned longitudinal visits of SVIP, Searchlight and Cardiff (intervals >= 0.5 y; 05_age_varying.R). Head circumference: full SVIP v10 series + Cardiff, G2MH, SSC (07_hc_trajectories.R). Plausibility QC (00_qc_growth.R) applied to the genomatrix and the SVIP v10 series. 22q11.2 dup has too few repeat measures and is not shown in the figure.
- Age-varying p: mixed GAM with s(age, by = BMI); p that the BMI effect differs from 0 at any age. All ages: the same model with a single constant BMI coefficient (effect, SE, t = effect / SE, p). Windows: linear mixed model (lmer) with BMI x age-window terms for height (< 8, 8-13, 13-18 y) and GAM-based windows for head circumference (2-5, 5-18 y); t = estimate / SE.


| Outcome | Group | Term | Children | Intervals | Effect per SD of BMI Z (SE) | t | p | Age-varying p |
|---|---|---|---|---|---|---|---|---|
| Height | 16p11.2 del | all ages (constant BMI coefficient) | 116 | 224 | +0.086 (0.041) | 2.12 | 0.04 | 0.006 |
| Height | 16p11.2 del | window 13-18 y |  28 |  41 | -0.007 (0.083) | -0.08 | 0.93 |  |
| Height | 16p11.2 del | window 8-13 y |  60 |  90 | +0.003 (0.063) | 0.05 | 0.96 |  |
| Height | 16p11.2 del | window <8 y |  61 |  93 | +0.234 (0.061) | 3.86 | 1.6e-04 |  |
| Height | 16p11.2 dup | all ages (constant BMI coefficient) |  48 |  84 | +0.003 (0.064) | 0.05 | 0.96 | 0.11 |
| Height | 16p11.2 dup | window 13-18 y |  12 |  24 | -0.202 (0.106) | -1.91 | 0.06 |  |
| Height | 16p11.2 dup | window 8-13 y |  25 |  40 | +0.107 (0.089) | 1.20 | 0.23 |  |
| Height | 16p11.2 dup | window <8 y |  20 |  20 | +0.051 (0.120) | 0.42 | 0.68 |  |
| Height | 22q11.2 del | all ages (constant BMI coefficient) |  73 | 102 | -0.038 (0.046) | -0.82 | 0.41 | 0.60 |
| Height | 22q11.2 del | window 13-18 y |  22 |  23 | -0.026 (0.119) | -0.22 | 0.83 |  |
| Height | 22q11.2 del | window 8-13 y |  52 |  65 | -0.057 (0.060) | -0.94 | 0.35 |  |
| Height | 22q11.2 del | window <8 y |  14 |  14 | +0.026 (0.097) | 0.26 | 0.79 |  |
| Height | Pooled | all ages (constant BMI coefficient) | 248 | 422 | +0.026 (0.026) | 1.00 | 0.32 | 0.02 |
| Height | Pooled | window 13-18 y |  63 |  89 | -0.066 (0.050) | -1.31 | 0.19 |  |
| Height | Pooled | window 8-13 y | 143 | 201 | +0.007 (0.036) | 0.20 | 0.84 |  |
| Height | Pooled | window <8 y | 100 | 132 | +0.127 (0.044) | 2.86 | 0.005 |  |
| Head circumference | 16p11.2 del | all ages (constant BMI coefficient) |  91 | 290 | +0.037 (0.034) | 1.08 | 0.28 | 0.49 |
| Head circumference | 16p11.2 del | window 2-5 y |  31 |  44 | +0.075 (0.088) | 0.85 | 0.40 |  |
| Head circumference | 16p11.2 del | window 5-18 y |  35 |  42 | +0.053 (0.086) | 0.61 | 0.54 |  |
| Head circumference | 16p11.2 del | window <2 y |  55 | 204 | +0.032 (0.040) | 0.79 | 0.43 |  |
| Head circumference | 16p11.2 dup | all ages (constant BMI coefficient) |  46 | 161 | +0.058 (0.049) | 1.17 | 0.24 | 0.29 |
| Head circumference | 16p11.2 dup | window 2-5 y |  25 |  37 | +0.098 (0.096) | 1.02 | 0.31 |  |
| Head circumference | 16p11.2 dup | window 5-18 y |  10 |  11 | -0.085 (0.206) | -0.41 | 0.68 |  |
| Head circumference | 16p11.2 dup | window <2 y |  30 | 113 | +0.061 (0.068) | 0.90 | 0.37 |  |
| Head circumference | 22q11.2 del | all ages (constant BMI coefficient) |  50 |  56 | +0.193 (0.111) | 1.75 | 0.09 | 0.23 |
| Head circumference | Pooled | all ages (constant BMI coefficient) | 193 | 514 | +0.069 (0.026) | 2.64 | 0.009 | 0.004 |
| Head circumference | Pooled | window 2-5 y |  56 |  81 | +0.091 (0.067) | 1.35 | 0.18 |  |
| Head circumference | Pooled | window 5-18 y | 101 | 116 | +0.143 (0.052) | 2.77 | 0.006 |  |
| Head circumference | Pooled | window <2 y |  85 | 317 | +0.032 (0.035) | 0.94 | 0.35 |  |

- Table S3: prospective effect of BMI on subsequent growth in CNV carriers < 18 y (fig_supp_prospective). Outcome: change in height Z (WHO) or head-circumference Z (UK90) from one visit to the next, per 1 SD of BMI Z at the first visit, adjusted for Z at the visit, interval length, sex, locus (pooled) and cohort, with a child random intercept. Height: cleaned longitudinal visits of SVIP, Searchlight and Cardiff (intervals >= 0.5 y; 05_age_varying.R). Head circumference: full SVIP v10 series + Cardiff, G2MH, SSC (07_hc_trajectories.R). Plausibility QC (00_qc_growth.R) applied to the genomatrix and the SVIP v10 series. 22q11.2 dup has too few repeat measures and is not shown in the figure.
- Age-varying p: mixed GAM with s(age, by = BMI); p that the BMI effect differs from 0 at any age. All ages: the same model with a single constant BMI coefficient (effect, SE, t = effect / SE, p). Windows: linear mixed model (lmer) with BMI x age-window terms for height (< 8, 8-13, 13-18 y) and GAM-based windows for head circumference (2-5, 5-18 y); t = estimate / SE.

## Table S3. Prospective effect of BMI on subsequent growth (carriers < 18 y)

| Outcome | Group | Term | Children | Intervals | Effect per SD of BMI Z (SE) | t | p | Age-varying p |
|---|---|---|---|---|---|---|---|---|
| Height | 16p11.2 del | all ages (constant BMI coefficient) | 116 | 224 | +0.086 (0.041) | 2.12 | 0.04 | 0.006 |
| Height | 16p11.2 del | window 13-18 y |  28 |  41 | -0.007 (0.083) | -0.08 | 0.93 |  |
| Height | 16p11.2 del | window 8-13 y |  60 |  90 | +0.003 (0.063) | 0.05 | 0.96 |  |
| Height | 16p11.2 del | window <8 y |  61 |  93 | +0.234 (0.061) | 3.86 | 1.6e-04 |  |
| Height | 16p11.2 dup | all ages (constant BMI coefficient) |  48 |  84 | +0.003 (0.064) | 0.05 | 0.96 | 0.11 |
| Height | 16p11.2 dup | window 13-18 y |  12 |  24 | -0.202 (0.106) | -1.91 | 0.06 |  |
| Height | 16p11.2 dup | window 8-13 y |  25 |  40 | +0.107 (0.089) | 1.20 | 0.23 |  |
| Height | 16p11.2 dup | window <8 y |  20 |  20 | +0.051 (0.120) | 0.42 | 0.68 |  |
| Height | 22q11.2 del | all ages (constant BMI coefficient) |  73 | 102 | -0.038 (0.046) | -0.82 | 0.41 | 0.60 |
| Height | 22q11.2 del | window 13-18 y |  22 |  23 | -0.026 (0.119) | -0.22 | 0.83 |  |
| Height | 22q11.2 del | window 8-13 y |  52 |  65 | -0.057 (0.060) | -0.94 | 0.35 |  |
| Height | 22q11.2 del | window <8 y |  14 |  14 | +0.026 (0.097) | 0.26 | 0.79 |  |
| Height | Pooled | all ages (constant BMI coefficient) | 248 | 422 | +0.026 (0.026) | 1.00 | 0.32 | 0.02 |
| Height | Pooled | window 13-18 y |  63 |  89 | -0.066 (0.050) | -1.31 | 0.19 |  |
| Height | Pooled | window 8-13 y | 143 | 201 | +0.007 (0.036) | 0.20 | 0.84 |  |
| Height | Pooled | window <8 y | 100 | 132 | +0.127 (0.044) | 2.86 | 0.005 |  |
| Head circumference | 16p11.2 del | all ages (constant BMI coefficient) |  91 | 287 | +0.081 (0.039) | 2.09 | 0.04 | 0.11 |
| Head circumference | 16p11.2 del | window 2-5 y |  31 |  44 | +0.082 (0.088) | 0.93 | 0.35 |  |
| Head circumference | 16p11.2 del | window 5-18 y |  34 |  40 | +0.052 (0.095) | 0.55 | 0.59 |  |
| Head circumference | 16p11.2 del | window <2 y |  55 | 203 | +0.091 (0.046) | 1.97 | 0.05 |  |
| Head circumference | 16p11.2 dup | all ages (constant BMI coefficient) |  46 | 159 | +0.042 (0.049) | 0.87 | 0.38 | 0.42 |
| Head circumference | 16p11.2 dup | window 2-5 y |  25 |  36 | +0.081 (0.096) | 0.85 | 0.40 |  |
| Head circumference | 16p11.2 dup | window 5-18 y |  10 |  11 | -0.114 (0.205) | -0.56 | 0.58 |  |
| Head circumference | 16p11.2 dup | window <2 y |  30 | 112 | +0.055 (0.067) | 0.82 | 0.41 |  |
| Head circumference | 22q11.2 del | all ages (constant BMI coefficient) |  50 |  56 | +0.193 (0.111) | 1.75 | 0.09 | 0.23 |
| Head circumference | Pooled | all ages (constant BMI coefficient) | 193 | 509 | +0.096 (0.028) | 3.42 | 6.7e-04 | 7.5e-04 |
| Head circumference | Pooled | window 2-5 y |  56 |  80 | +0.087 (0.067) | 1.31 | 0.19 |  |
| Head circumference | Pooled | window 5-18 y | 100 | 114 | +0.153 (0.054) | 2.85 | 0.005 |  |
| Head circumference | Pooled | window <2 y |  85 | 315 | +0.074 (0.038) | 1.97 | 0.05 |  |

- Table S4: sensitivity of the main-figure BMI effects to sex and polygenic background. The main-figure tests are repeated within strata: sex (male / female; carriers and unaffected controls both stratified); childhood obesity PGS (Bradfield et al. 2019 trans-ancestral GWAS, SBayesRC weights), adult BMI PGS (Yengo et al. 2018), height PGS (Yengo et al. 2022) and autism PGS (April 2026 ASD GWAS; full weights for G2MH carriers, leave-SSC-out weights for SSC controls), each below / above the median. PGS were computed with the same weights on the G2MH WGS cohort (carriers) and on the SSC genotype set (controls), residualised on 10 ancestry principal components and z-scored within cohort; the median split is within cohort. Carriers from SVIP, Searchlight and Cardiff have no PGS in this pipeline and appear only in the sex strata.
- Columns: within-stratum BMI slope (SD of trait per SD of BMI, GAM y ~ BMI + s(age) [+ locus]) and age x BMI p; CNV x BMI and CNV x age x BMI versus the unaffected controls of the same stratum (models as in Tables S1-S2); BMI x stratum p: whether the within-carrier BMI slope differs between the two strata (age-adjusted, locus covariate when pooled); PGS main effect and BMI x PGS (continuous) p: within carriers, GAM y ~ s(age) + BMI * PGS [+ locus]. CNB has no control group with BMI.

## Table S4. Main-figure BMI effects within strata of sex, childhood obesity PGS, adult BMI PGS, height PGS and autism PGS

| Stratification | Stratum | Trait | Group | n | n controls | BMI slope (SE) | BMI slope p | age x BMI p | CNV x BMI vs controls (SE) | CNV x BMI p | CNV x age x BMI p | BMI x stratum (carriers) p | PGS main effect p | BMI x PGS (continuous) p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Adult BMI PGS | above median | Height | 16p11.2 del |  61 | 2504 | +0.503 (0.130) | 2.9e-04 | 0.005 | +0.374 (0.100) | 1.8e-04 | 0.16 | 0.056 | 0.99 | 0.048 |
| Adult BMI PGS | below median | Height | 16p11.2 del |  57 | 2506 | +0.185 (0.090) | 0.045 | 0.16 | +0.098 (0.097) | 0.31 | 0.2 | 0.056 | 0.99 | 0.048 |
| Adult BMI PGS | above median | Height | 16p11.2 dup |  31 | 2504 | +0.163 (0.168) | 0.34 | 0.056 | +0.043 (0.139) | 0.76 | 0.52 | 0.39 | 0.84 | 0.064 |
| Adult BMI PGS | below median | Height | 16p11.2 dup |  36 | 2506 | -0.045 (0.153) | 0.77 | 0.04 | -0.091 (0.144) | 0.53 | 0.0012 | 0.39 | 0.84 | 0.064 |
| Adult BMI PGS | above median | Height | 22q11.2 del | 128 | 2504 | +0.135 (0.076) | 0.078 | 0.0094 | +0.004 (0.074) | 0.96 | 1 | 0.27 | 0.14 | 0.67 |
| Adult BMI PGS | below median | Height | 22q11.2 del | 140 | 2506 | +0.015 (0.078) | 0.84 | 0.22 | -0.053 (0.072) | 0.47 | 0.81 | 0.27 | 0.14 | 0.67 |
| Adult BMI PGS | above median | Height | 22q11.2 dup |  51 | 2504 | +0.270 (0.134) | 0.05 | 0.64 | +0.150 (0.135) | 0.27 | 0.12 | 0.1 | 0.34 | 0.19 |
| Adult BMI PGS | below median | Height | 22q11.2 dup |  46 | 2506 | -0.043 (0.119) | 0.72 | 0.19 | -0.106 (0.122) | 0.38 | 0.27 | 0.1 | 0.34 | 0.19 |
| Adult BMI PGS | above median | Height | Overgrowth CNVs | 112 | 2504 | +0.406 (0.095) | 4.2e-05 | 0.021 | +0.259 (0.083) | 0.0019 | 0.18 | 0.018 | 0.46 | 0.019 |
| Adult BMI PGS | below median | Height | Overgrowth CNVs | 103 | 2506 | +0.107 (0.072) | 0.14 | 0.026 | +0.019 (0.078) | 0.81 | 0.99 | 0.018 | 0.46 | 0.019 |
| Adult BMI PGS | above median | Height | Undergrowth CNVs | 159 | 2504 | +0.144 (0.070) | 0.04 | 6.7e-04 | +0.011 (0.068) | 0.87 | 0.14 | 0.21 | 0.14 | 0.33 |
| Adult BMI PGS | below median | Height | Undergrowth CNVs | 176 | 2506 | +0.005 (0.069) | 0.94 | 0.37 | -0.059 (0.066) | 0.38 | 0.87 | 0.21 | 0.14 | 0.33 |
| Adult BMI PGS | above median | Height | All CNVs | 271 | 2504 | +0.231 (0.059) | 1.2e-04 | 1.4e-05 | +0.088 (0.058) | 0.13 | 0.18 | 0.0042 | 0.72 | 0.017 |
| Adult BMI PGS | below median | Height | All CNVs | 279 | 2506 | +0.022 (0.053) | 0.67 | 0.017 | -0.052 (0.054) | 0.33 | 0.17 | 0.0042 | 0.72 | 0.017 |
| Adult BMI PGS | above median | Head circumference | 16p11.2 del |  50 | 2492 | +0.650 (0.118) | 1.4e-06 | 0.067 | +0.285 (0.111) | 0.01 | 0.019 | 0.056 | 0.59 | 0.075 |
| Adult BMI PGS | below median | Head circumference | 16p11.2 del |  48 | 2493 | +0.337 (0.111) | 0.004 | 0.33 | +0.013 (0.106) | 0.9 | 0.0036 | 0.056 | 0.59 | 0.075 |
| Adult BMI PGS | above median | Head circumference | 16p11.2 dup |  23 | 2492 | +0.018 (0.163) | 0.91 | 0.62 | -0.339 (0.174) | 0.051 | 0.018 | 0.77 | 0.69 | 0.93 |
| Adult BMI PGS | below median | Head circumference | 16p11.2 dup |  31 | 2493 | +0.128 (0.164) | 0.44 | 0.014 | -0.218 (0.161) | 0.17 | 0.057 | 0.77 | 0.69 | 0.93 |
| Adult BMI PGS | above median | Head circumference | 22q11.2 del |  91 | 2492 | +0.397 (0.090) | 2.8e-05 | 0.27 | +0.035 (0.088) | 0.7 | 0.97 | 0.16 | 0.95 | 0.35 |
| Adult BMI PGS | below median | Head circumference | 22q11.2 del |  80 | 2493 | +0.197 (0.114) | 0.089 | 0.67 | -0.150 (0.105) | 0.15 | 8.3e-04 | 0.16 | 0.95 | 0.35 |
| Adult BMI PGS | above median | Head circumference | 22q11.2 dup |  38 | 2492 | +0.680 (0.165) | 2.2e-04 | 0.38 | +0.307 (0.163) | 0.06 | 1 | 0.17 | 0.82 | 0.15 |
| Adult BMI PGS | below median | Head circumference | 22q11.2 dup |  34 | 2493 | +0.376 (0.159) | 0.024 | 0.42 | -0.003 (0.145) | 0.99 | 0.0038 | 0.17 | 0.82 | 0.15 |
| Adult BMI PGS | above median | Head circumference | Overgrowth CNVs |  88 | 2492 | +0.654 (0.096) | 1.2e-09 | 0.44 | +0.286 (0.094) | 0.0024 | 0.019 | 0.041 | 0.92 | 0.027 |
| Adult BMI PGS | below median | Head circumference | Overgrowth CNVs |  82 | 2493 | +0.350 (0.090) | 2.1e-04 | 0.53 | +0.007 (0.088) | 0.94 | 8.7e-04 | 0.041 | 0.92 | 0.027 |
| Adult BMI PGS | above median | Head circumference | Undergrowth CNVs | 114 | 2492 | +0.338 (0.082) | 7.5e-05 | 0.65 | -0.023 (0.081) | 0.78 | 1 | 0.14 | 0.77 | 0.29 |
| Adult BMI PGS | below median | Head circumference | Undergrowth CNVs | 111 | 2493 | +0.160 (0.094) | 0.092 | 0.37 | -0.192 (0.089) | 0.032 | 0.002 | 0.14 | 0.77 | 0.29 |
| Adult BMI PGS | above median | Head circumference | All CNVs | 202 | 2492 | +0.468 (0.063) | 3.2e-12 | 0.43 | +0.103 (0.066) | 0.12 | 0.99 | 0.017 | 0.92 | 0.031 |
| Adult BMI PGS | below median | Head circumference | All CNVs | 193 | 2493 | +0.250 (0.065) | 1.7e-04 | 0.095 | -0.100 (0.067) | 0.13 | 0.0049 | 0.017 | 0.92 | 0.031 |
| Adult BMI PGS | above median | SRS | 16p11.2 del |  44 | 2494 | -0.028 (0.292) | 0.92 | 1 | -0.075 (0.109) | 0.49 | 0.21 | 0.3 | 0.38 | 0.18 |
| Adult BMI PGS | below median | SRS | 16p11.2 del |  30 | 2494 | +0.592 (0.310) | 0.068 | 0.4 | +0.597 (0.147) | 4.9e-05 | 1 | 0.3 | 0.38 | 0.18 |
| Adult BMI PGS | below median | SRS | 16p11.2 dup |  20 | 2494 | -0.022 (0.379) | 0.95 | 0.26 | -0.095 (0.207) | 0.64 | 1 | 0.44 | 0.73 | 0.15 |
| Adult BMI PGS | above median | SRS | 22q11.2 del |  88 | 2494 | +0.305 (0.153) | 0.05 | 0.78 | +0.290 (0.079) | 2.6e-04 | 0.2 | 0.92 | 0.37 | 0.6 |
| Adult BMI PGS | below median | SRS | 22q11.2 del |  98 | 2494 | +0.318 (0.152) | 0.04 | 0.59 | +0.307 (0.081) | 1.6e-04 | 0.99 | 0.92 | 0.37 | 0.6 |
| Adult BMI PGS | above median | SRS | 22q11.2 dup |  31 | 2494 | -0.309 (0.355) | 0.39 | 0.25 | -0.291 (0.167) | 0.081 | 0.18 | 0.75 | 0.33 | 0.83 |
| Adult BMI PGS | below median | SRS | 22q11.2 dup |  29 | 2494 | -0.029 (0.380) | 0.94 | 0.26 | -0.066 (0.195) | 0.74 | 1 | 0.75 | 0.33 | 0.83 |
| Adult BMI PGS | above median | SRS | Overgrowth CNVs |  75 | 2494 | -0.126 (0.225) | 0.58 | 0.44 | -0.148 (0.096) | 0.12 | 1 | 0.21 | 0.15 | 0.13 |
| Adult BMI PGS | below median | SRS | Overgrowth CNVs |  59 | 2494 | +0.385 (0.231) | 0.1 | 0.72 | +0.350 (0.118) | 0.0031 | 0.18 | 0.21 | 0.15 | 0.13 |
| Adult BMI PGS | above median | SRS | Undergrowth CNVs |  97 | 2494 | +0.335 (0.153) | 0.031 | 0.74 | +0.317 (0.077) | 4.3e-05 | 0.97 | 0.84 | 0.24 | 0.93 |
| Adult BMI PGS | below median | SRS | Undergrowth CNVs | 118 | 2494 | +0.273 (0.142) | 0.058 | 0.88 | +0.252 (0.077) | 0.0011 | 1 | 0.84 | 0.24 | 0.93 |
| Adult BMI PGS | above median | SRS | All CNVs | 172 | 2494 | +0.154 (0.131) | 0.24 | 0.66 | +0.122 (0.067) | 0.071 | 0.97 | 0.69 | 0.95 | 0.53 |
| Adult BMI PGS | below median | SRS | All CNVs | 177 | 2494 | +0.293 (0.120) | 0.016 | 0.88 | +0.282 (0.068) | 3.9e-05 | 1 | 0.69 | 0.95 | 0.53 |
| Adult BMI PGS | above median | SCQ | 16p11.2 del |  35 | 532 | +0.438 (0.473) | 0.36 | 0.68 | +0.361 (0.169) | 0.034 | 1 | 0.98 | 0.33 | 0.94 |
| Adult BMI PGS | below median | SCQ | 16p11.2 del |  21 | 498 | +0.422 (0.411) | 0.32 | 0.61 | +0.448 (0.151) | 0.0032 | 0.27 | 0.98 | 0.33 | 0.94 |
| Adult BMI PGS | above median | SCQ | 22q11.2 del |  72 | 532 | +0.332 (0.235) | 0.16 | 0.61 | +0.299 (0.120) | 0.013 | 0.16 | 0.52 | 0.83 | 0.36 |
| Adult BMI PGS | below median | SCQ | 22q11.2 del |  79 | 498 | +0.518 (0.247) | 0.039 | 0.5 | +0.645 (0.145) | 1.1e-05 | 0.85 | 0.52 | 0.83 | 0.36 |
| Adult BMI PGS | above median | SCQ | 22q11.2 dup |  26 | 532 | +0.044 (0.563) | 0.94 | 0.39 | +0.085 (0.194) | 0.66 | 0.97 | 0.42 | 0.18 | 0.47 |
| Adult BMI PGS | below median | SCQ | 22q11.2 dup |  18 | 498 | +0.831 (0.689) | 0.25 | 0.18 | +0.856 (0.258) | 9.8e-04 | 1 | 0.42 | 0.18 | 0.47 |
| Adult BMI PGS | above median | SCQ | Overgrowth CNVs |  61 | 532 | +0.278 (0.357) | 0.44 | 0.47 | +0.208 (0.154) | 0.18 | 0.52 | 0.39 | 0.89 | 0.91 |
| Adult BMI PGS | below median | SCQ | Overgrowth CNVs |  39 | 498 | +0.501 (0.364) | 0.18 | 0.72 | +0.594 (0.173) | 6.6e-04 | 1 | 0.39 | 0.89 | 0.91 |
| Adult BMI PGS | above median | SCQ | Undergrowth CNVs |  78 | 532 | +0.308 (0.242) | 0.21 | 0.69 | +0.260 (0.124) | 0.037 | 1 | 0.52 | 0.44 | 0.47 |
| Adult BMI PGS | below median | SCQ | Undergrowth CNVs |  90 | 498 | +0.469 (0.240) | 0.054 | 0.43 | +0.574 (0.145) | 8.9e-05 | 0.5 | 0.52 | 0.44 | 0.47 |
| Adult BMI PGS | above median | SCQ | All CNVs | 139 | 532 | +0.227 (0.208) | 0.28 | 0.63 | +0.224 (0.127) | 0.077 | 0.6 | 0.4 | 0.86 | 0.6 |
| Adult BMI PGS | below median | SCQ | All CNVs | 129 | 498 | +0.482 (0.199) | 0.017 | 0.64 | +0.560 (0.132) | 2.6e-05 | 0.62 | 0.4 | 0.86 | 0.6 |
| Adult BMI PGS | above median | CNB accuracy | 16p11.2 del |  57 |  | +0.025 (0.063) | 0.69 | 0.32 |  |  |  | 0.74 | 0.42 | 0.43 |
| Adult BMI PGS | below median | CNB accuracy | 16p11.2 del |  44 |  | +0.015 (0.065) | 0.82 | 0.19 |  |  |  | 0.74 | 0.42 | 0.43 |
| Adult BMI PGS | below median | CNB accuracy | 16p11.2 dup |  26 |  | -0.105 (0.097) | 0.29 | 0.21 |  |  |  | 0.31 | 0.39 | 0.57 |
| Adult BMI PGS | above median | CNB accuracy | 22q11.2 del | 129 |  | +0.015 (0.047) | 0.75 | 0.57 |  |  |  | 0.11 | 0.93 | 0.54 |
| Adult BMI PGS | below median | CNB accuracy | 22q11.2 del | 153 |  | -0.088 (0.044) | 0.05 | 0.21 |  |  |  | 0.11 | 0.93 | 0.54 |
| Adult BMI PGS | above median | CNB accuracy | 22q11.2 dup |  63 |  | +0.139 (0.085) | 0.11 | 0.67 |  |  |  | 0.34 | 0.71 | 0.77 |
| Adult BMI PGS | below median | CNB accuracy | 22q11.2 dup |  47 |  | +0.025 (0.077) | 0.75 | 0.5 |  |  |  | 0.34 | 0.71 | 0.77 |
| Adult BMI PGS | above median | CNB accuracy | Overgrowth CNVs | 120 |  | +0.083 (0.053) | 0.12 | 0.38 |  |  |  | 0.55 | 0.61 | 0.67 |
| Adult BMI PGS | below median | CNB accuracy | Overgrowth CNVs |  91 |  | +0.013 (0.050) | 0.8 | 0.094 |  |  |  | 0.55 | 0.61 | 0.67 |
| Adult BMI PGS | above median | CNB accuracy | Undergrowth CNVs | 142 |  | -0.002 (0.044) | 0.97 | 0.81 |  |  |  | 0.11 | 0.78 | 0.48 |
| Adult BMI PGS | below median | CNB accuracy | Undergrowth CNVs | 179 |  | -0.086 (0.041) | 0.037 | 0.15 |  |  |  | 0.11 | 0.78 | 0.48 |
| Adult BMI PGS | above median | CNB accuracy | All CNVs | 262 |  | +0.033 (0.034) | 0.34 | 0.37 |  |  |  | 0.12 | 0.45 | 0.35 |
| Adult BMI PGS | below median | CNB accuracy | All CNVs | 270 |  | -0.051 (0.032) | 0.11 | 0.025 |  |  |  | 0.12 | 0.45 | 0.35 |
| Adult BMI PGS | above median | CNB speed | 16p11.2 del |  57 |  | +0.038 (0.062) | 0.54 | 0.59 |  |  |  | 0.16 | 0.76 | 0.13 |
| Adult BMI PGS | below median | CNB speed | 16p11.2 del |  46 |  | -0.114 (0.101) | 0.27 | 0.6 |  |  |  | 0.16 | 0.76 | 0.13 |
| Adult BMI PGS | below median | CNB speed | 16p11.2 dup |  26 |  | +0.116 (0.124) | 0.36 | 0.99 |  |  |  | 0.47 | 0.34 | 0.49 |
| Adult BMI PGS | above median | CNB speed | 22q11.2 del | 133 |  | -0.012 (0.043) | 0.78 | 0.27 |  |  |  | 0.69 | 0.14 | 0.55 |
| Adult BMI PGS | below median | CNB speed | 22q11.2 del | 154 |  | -0.035 (0.056) | 0.53 | 0.24 |  |  |  | 0.69 | 0.14 | 0.55 |
| Adult BMI PGS | above median | CNB speed | 22q11.2 dup |  63 |  | +0.118 (0.164) | 0.47 | 0.85 |  |  |  | 0.48 | 0.67 | 0.42 |
| Adult BMI PGS | below median | CNB speed | 22q11.2 dup |  47 |  | -0.034 (0.079) | 0.67 | 0.0092 |  |  |  | 0.48 | 0.67 | 0.42 |
| Adult BMI PGS | above median | CNB speed | Overgrowth CNVs | 120 |  | +0.097 (0.088) | 0.27 | 0.88 |  |  |  | 0.15 | 0.79 | 0.14 |
| Adult BMI PGS | below median | CNB speed | Overgrowth CNVs |  93 |  | -0.081 (0.061) | 0.19 | 0.49 |  |  |  | 0.15 | 0.79 | 0.14 |
| Adult BMI PGS | above median | CNB speed | Undergrowth CNVs | 146 |  | +0.011 (0.044) | 0.8 | 0.3 |  |  |  | 0.47 | 0.33 | 0.85 |
| Adult BMI PGS | below median | CNB speed | Undergrowth CNVs | 180 |  | -0.018 (0.050) | 0.71 | 0.54 |  |  |  | 0.47 | 0.33 | 0.85 |
| Adult BMI PGS | above median | CNB speed | All CNVs | 266 |  | +0.037 (0.045) | 0.42 | 0.53 |  |  |  | 0.2 | 0.54 | 0.45 |
| Adult BMI PGS | below median | CNB speed | All CNVs | 273 |  | -0.040 (0.039) | 0.31 | 0.23 |  |  |  | 0.2 | 0.54 | 0.45 |
| Autism PGS (2026) | above median | Height | 16p11.2 del |  59 | 2504 | +0.482 (0.130) | 4.9e-04 | 0.017 | +0.363 (0.094) | 1.1e-04 | 0.0017 | 0.079 | 0.016 | 0.04 |
| Autism PGS (2026) | below median | Height | 16p11.2 del |  59 | 2506 | +0.200 (0.086) | 0.023 | 0.046 | +0.072 (0.100) | 0.47 | 0.99 | 0.079 | 0.016 | 0.04 |
| Autism PGS (2026) | above median | Height | 16p11.2 dup |  41 | 2504 | +0.047 (0.125) | 0.71 | 0.51 | -0.043 (0.111) | 0.7 | 0.81 | 0.45 | 0.1 | 0.75 |
| Autism PGS (2026) | below median | Height | 16p11.2 dup |  26 | 2506 | +0.221 (0.237) | 0.36 | 0.78 | +0.119 (0.209) | 0.57 | 4.6e-05 | 0.45 | 0.1 | 0.75 |
| Autism PGS (2026) | above median | Height | 22q11.2 del | 131 | 2504 | +0.007 (0.077) | 0.92 | 0.13 | -0.097 (0.070) | 0.17 | 1 | 0.3 | 0.86 | 0.52 |
| Autism PGS (2026) | below median | Height | 22q11.2 del | 137 | 2506 | +0.111 (0.073) | 0.13 | 0.28 | +0.008 (0.072) | 0.91 | 1.7e-04 | 0.3 | 0.86 | 0.52 |
| Autism PGS (2026) | above median | Height | 22q11.2 dup |  46 | 2504 | -0.013 (0.152) | 0.93 | 0.047 | -0.089 (0.137) | 0.52 | 1 | 0.26 | 0.086 | 0.011 |
| Autism PGS (2026) | below median | Height | 22q11.2 dup |  51 | 2506 | +0.192 (0.113) | 0.095 | 0.3 | +0.072 (0.123) | 0.56 | 0.017 | 0.26 | 0.086 | 0.011 |
| Autism PGS (2026) | above median | Height | Overgrowth CNVs | 105 | 2504 | +0.323 (0.100) | 0.0016 | 2.4e-04 | +0.221 (0.079) | 0.0055 | 0.0017 | 0.28 | 0.57 | 0.99 |
| Autism PGS (2026) | below median | Height | Overgrowth CNVs | 110 | 2506 | +0.203 (0.067) | 0.0033 | 0.044 | +0.074 (0.080) | 0.36 | 0.22 | 0.28 | 0.57 | 0.99 |
| Autism PGS (2026) | above median | Height | Undergrowth CNVs | 172 | 2504 | +0.014 (0.065) | 0.83 | 0.25 | -0.081 (0.062) | 0.19 | 0.46 | 0.18 | 0.45 | 0.55 |
| Autism PGS (2026) | below median | Height | Undergrowth CNVs | 163 | 2506 | +0.123 (0.070) | 0.079 | 0.17 | +0.019 (0.069) | 0.78 | 1.1e-04 | 0.18 | 0.45 | 0.55 |
| Autism PGS (2026) | above median | Height | All CNVs | 277 | 2504 | +0.132 (0.058) | 0.023 | 3.5e-04 | +0.015 (0.053) | 0.77 | 0.5 | 0.63 | 0.66 | 0.58 |
| Autism PGS (2026) | below median | Height | All CNVs | 273 | 2506 | +0.127 (0.052) | 0.015 | 1.2e-04 | +0.012 (0.056) | 0.83 | 7.6e-05 | 0.63 | 0.66 | 0.58 |
| Autism PGS (2026) | above median | Head circumference | 16p11.2 del |  46 | 2494 | +0.529 (0.124) | 1.1e-04 | 0.34 | +0.188 (0.107) | 0.079 | 0.03 | 0.56 | 0.29 | 0.45 |
| Autism PGS (2026) | below median | Head circumference | 16p11.2 del |  52 | 2491 | +0.437 (0.105) | 1.2e-04 | 0.45 | +0.063 (0.107) | 0.55 | 0.0087 | 0.56 | 0.29 | 0.45 |
| Autism PGS (2026) | above median | Head circumference | 16p11.2 dup |  32 | 2494 | +0.143 (0.119) | 0.24 | 0.84 | -0.206 (0.137) | 0.13 | 0.06 | 0.76 | 0.8 | 0.5 |
| Autism PGS (2026) | below median | Head circumference | 16p11.2 dup |  22 | 2491 | +0.036 (0.219) | 0.87 | 0.44 | -0.314 (0.222) | 0.16 | 0.0091 | 0.76 | 0.8 | 0.5 |
| Autism PGS (2026) | above median | Head circumference | 22q11.2 del |  83 | 2494 | +0.276 (0.115) | 0.018 | 0.22 | -0.071 (0.099) | 0.47 | 0.061 | 0.61 | 0.42 | 0.79 |
| Autism PGS (2026) | below median | Head circumference | 22q11.2 del |  88 | 2491 | +0.362 (0.081) | 2.6e-05 | 0.82 | -0.019 (0.087) | 0.83 | 0.99 | 0.61 | 0.42 | 0.79 |
| Autism PGS (2026) | above median | Head circumference | 22q11.2 dup |  36 | 2494 | +0.430 (0.189) | 0.03 | 0.058 | +0.060 (0.167) | 0.72 | 0.11 | 0.6 | 0.99 | 0.9 |
| Autism PGS (2026) | below median | Head circumference | 22q11.2 dup |  36 | 2491 | +0.549 (0.144) | 5.5e-04 | 0.12 | +0.159 (0.146) | 0.28 | 0.99 | 0.6 | 0.99 | 0.9 |
| Autism PGS (2026) | above median | Head circumference | Overgrowth CNVs |  82 | 2494 | +0.492 (0.103) | 8.4e-06 | 0.34 | +0.135 (0.092) | 0.14 | 0.12 | 0.93 | 0.49 | 0.72 |
| Autism PGS (2026) | below median | Head circumference | Overgrowth CNVs |  88 | 2491 | +0.476 (0.084) | 2.1e-07 | 0.73 | +0.099 (0.089) | 0.26 | 0.56 | 0.93 | 0.49 | 0.72 |
| Autism PGS (2026) | above median | Head circumference | Undergrowth CNVs | 115 | 2494 | +0.213 (0.089) | 0.019 | 0.13 | -0.140 (0.082) | 0.089 | 0.077 | 0.4 | 0.33 | 0.8 |
| Autism PGS (2026) | below median | Head circumference | Undergrowth CNVs | 110 | 2491 | +0.315 (0.078) | 9.6e-05 | 0.71 | -0.064 (0.082) | 0.44 | 0.98 | 0.4 | 0.33 | 0.8 |
| Autism PGS (2026) | above median | Head circumference | All CNVs | 197 | 2494 | +0.359 (0.068) | 3.9e-07 | 0.16 | +0.004 (0.065) | 0.95 | 0.1 | 0.39 | 0.64 | 0.96 |
| Autism PGS (2026) | below median | Head circumference | All CNVs | 198 | 2491 | +0.390 (0.057) | 1.2e-10 | 0.77 | +0.015 (0.064) | 0.81 | 0.61 | 0.39 | 0.64 | 0.96 |
| Autism PGS (2026) | above median | SRS | 16p11.2 del |  40 | 2493 | -0.111 (0.254) | 0.67 | 0.47 | -0.041 (0.116) | 0.72 | 0.98 | 0.21 | 0.17 | 0.058 |
| Autism PGS (2026) | below median | SRS | 16p11.2 del |  34 | 2495 | +0.478 (0.350) | 0.18 | 0.95 | +0.413 (0.125) | 1.0e-03 | 0.71 | 0.21 | 0.17 | 0.058 |
| Autism PGS (2026) | above median | SRS | 16p11.2 dup |  17 | 2493 | +0.661 (0.446) | 0.16 | 0.48 | +0.483 (0.228) | 0.034 | 0.047 | 0.031 | 0.076 | 0.18 |
| Autism PGS (2026) | above median | SRS | 22q11.2 del |  89 | 2493 | -0.033 (0.173) | 0.85 | 0.76 | -0.023 (0.083) | 0.78 | 0.072 | 0.0031 | 0.3 | 0.074 |
| Autism PGS (2026) | below median | SRS | 22q11.2 del |  97 | 2495 | +0.607 (0.119) | 1.8e-06 | 0.85 | +0.561 (0.074) | 3.4e-14 | 1 | 0.0031 | 0.3 | 0.074 |
| Autism PGS (2026) | above median | SRS | 22q11.2 dup |  31 | 2493 | -0.346 (0.399) | 0.39 | 0.62 | -0.532 (0.184) | 0.0039 | 0.97 | 0.6 | 0.13 | 0.65 |
| Autism PGS (2026) | below median | SRS | 22q11.2 dup |  29 | 2495 | -0.060 (0.368) | 0.87 | 0.81 | -0.108 (0.175) | 0.54 | 1 | 0.6 | 0.13 | 0.65 |
| Autism PGS (2026) | above median | SRS | Overgrowth CNVs |  71 | 2493 | -0.164 (0.211) | 0.44 | 0.51 | -0.142 (0.102) | 0.16 | 0.1 | 0.076 | 0.036 | 0.026 |
| Autism PGS (2026) | below median | SRS | Overgrowth CNVs |  63 | 2495 | +0.291 (0.256) | 0.26 | 0.96 | +0.215 (0.106) | 0.043 | 1 | 0.076 | 0.036 | 0.026 |
| Autism PGS (2026) | above median | SRS | Undergrowth CNVs | 106 | 2493 | +0.065 (0.162) | 0.69 | 0.22 | +0.048 (0.079) | 0.55 | 0.36 | 0.018 | 0.21 | 0.18 |
| Autism PGS (2026) | below median | SRS | Undergrowth CNVs | 109 | 2495 | +0.559 (0.123) | 1.6e-05 | 0.94 | +0.495 (0.073) | 1.8e-11 | 1 | 0.018 | 0.21 | 0.18 |
| Autism PGS (2026) | above median | SRS | All CNVs | 177 | 2493 | +0.000 (0.129) | 1 | 0.39 | -0.007 (0.068) | 0.91 | 0.1 | 0.003 | 0.036 | 0.019 |
| Autism PGS (2026) | below median | SRS | All CNVs | 172 | 2495 | +0.449 (0.118) | 2.1e-04 | 0.78 | +0.396 (0.066) | 2.3e-09 | 1 | 0.003 | 0.036 | 0.019 |
| Autism PGS (2026) | above median | SCQ | 16p11.2 del |  30 | 492 | +0.161 (0.387) | 0.68 | 0.64 | +0.182 (0.156) | 0.24 | 0.5 | 0.44 | 0.44 | 0.46 |
| Autism PGS (2026) | below median | SCQ | 16p11.2 del |  26 | 538 |  | n/a | n/a | +0.656 (0.213) | 0.0022 | 1 | 0.44 | 0.44 | 0.46 |
| Autism PGS (2026) | above median | SCQ | 22q11.2 del |  75 | 492 | +0.200 (0.289) | 0.49 | 0.26 | +0.081 (0.156) | 0.6 | 0.54 | 0.36 | 0.79 | 0.68 |
| Autism PGS (2026) | below median | SCQ | 22q11.2 del |  76 | 538 | +0.556 (0.191) | 0.0049 | 0.48 | +0.470 (0.105) | 9.2e-06 | 0.19 | 0.36 | 0.79 | 0.68 |
| Autism PGS (2026) | above median | SCQ | 22q11.2 dup |  22 | 492 | +0.742 (0.675) | 0.29 | 0.55 | +0.712 (0.237) | 0.0028 | 0.5 | 0.24 | 0.45 | 0.12 |
| Autism PGS (2026) | below median | SCQ | 22q11.2 dup |  22 | 538 | -0.325 (0.565) | 0.57 | 0.82 | -0.464 (0.220) | 0.035 | 1 | 0.24 | 0.45 | 0.12 |
| Autism PGS (2026) | above median | SCQ | Overgrowth CNVs |  52 | 492 | +0.366 (0.342) | 0.29 | 0.71 | +0.347 (0.153) | 0.024 | 1 | 0.81 | 0.27 | 0.88 |
| Autism PGS (2026) | below median | SCQ | Overgrowth CNVs |  48 | 538 | +0.126 (0.380) | 0.74 | 0.85 | -0.119 (0.163) | 0.47 | 1 | 0.81 | 0.27 | 0.88 |
| Autism PGS (2026) | above median | SCQ | Undergrowth CNVs |  84 | 492 | +0.196 (0.281) | 0.49 | 0.23 | +0.102 (0.156) | 0.51 | 0.9 | 0.47 | 0.4 | 0.57 |
| Autism PGS (2026) | below median | SCQ | Undergrowth CNVs |  84 | 538 | +0.500 (0.194) | 0.012 | 0.45 | +0.398 (0.106) | 1.9e-04 | 1 | 0.47 | 0.4 | 0.57 |
| Autism PGS (2026) | above median | SCQ | All CNVs | 136 | 492 | +0.262 (0.222) | 0.24 | 0.96 | +0.230 (0.144) | 0.11 | 1 | 0.39 | 0.12 | 0.81 |
| Autism PGS (2026) | below median | SCQ | All CNVs | 132 | 538 | +0.341 (0.182) | 0.063 | 0.27 | +0.358 (0.111) | 0.0014 | 0.93 | 0.39 | 0.12 | 0.81 |
| Autism PGS (2026) | above median | CNB accuracy | 16p11.2 del |  49 |  | -0.026 (0.073) | 0.72 | 0.63 |  |  |  | 0.28 | 0.3 | 0.44 |
| Autism PGS (2026) | below median | CNB accuracy | 16p11.2 del |  52 |  | +0.061 (0.056) | 0.28 | 0.29 |  |  |  | 0.28 | 0.3 | 0.44 |
| Autism PGS (2026) | above median | CNB accuracy | 16p11.2 dup |  24 |  | -0.240 (0.075) | 0.0043 | 0.84 |  |  |  | 0.081 | 0.61 | 0.24 |
| Autism PGS (2026) | below median | CNB accuracy | 16p11.2 dup |  15 |  | +0.043 (0.170) | 0.8 | 0.58 |  |  |  | 0.081 | 0.61 | 0.24 |
| Autism PGS (2026) | above median | CNB accuracy | 22q11.2 del | 132 |  | -0.032 (0.047) | 0.5 | 0.096 |  |  |  | 0.98 | 0.98 | 0.92 |
| Autism PGS (2026) | below median | CNB accuracy | 22q11.2 del | 150 |  | -0.034 (0.041) | 0.4 | 0.79 |  |  |  | 0.98 | 0.98 | 0.92 |
| Autism PGS (2026) | above median | CNB accuracy | 22q11.2 dup |  50 |  | +0.035 (0.097) | 0.72 | 0.89 |  |  |  | 0.84 | 0.05 | 0.71 |
| Autism PGS (2026) | below median | CNB accuracy | 22q11.2 dup |  60 |  | +0.022 (0.072) | 0.77 | 0.18 |  |  |  | 0.84 | 0.05 | 0.71 |
| Autism PGS (2026) | above median | CNB accuracy | Overgrowth CNVs |  99 |  | +0.010 (0.058) | 0.86 | 0.45 |  |  |  | 0.25 | 0.011 | 0.38 |
| Autism PGS (2026) | below median | CNB accuracy | Overgrowth CNVs | 112 |  | +0.048 (0.047) | 0.31 | 0.058 |  |  |  | 0.25 | 0.011 | 0.38 |
| Autism PGS (2026) | above median | CNB accuracy | Undergrowth CNVs | 156 |  | -0.051 (0.042) | 0.22 | 0.23 |  |  |  | 0.64 | 0.92 | 0.74 |
| Autism PGS (2026) | below median | CNB accuracy | Undergrowth CNVs | 165 |  | -0.025 (0.040) | 0.52 | 0.83 |  |  |  | 0.64 | 0.92 | 0.74 |
| Autism PGS (2026) | above median | CNB accuracy | All CNVs | 255 |  | -0.035 (0.034) | 0.3 | 0.21 |  |  |  | 0.56 | 0.1 | 0.84 |
| Autism PGS (2026) | below median | CNB accuracy | All CNVs | 277 |  | +0.006 (0.030) | 0.85 | 0.12 |  |  |  | 0.56 | 0.1 | 0.84 |
| Autism PGS (2026) | above median | CNB speed | 16p11.2 del |  49 |  | -0.046 (0.065) | 0.49 | 0.16 |  |  |  | 0.96 | 0.37 | 0.68 |
| Autism PGS (2026) | below median | CNB speed | 16p11.2 del |  54 |  | -0.020 (0.089) | 0.82 | 0.66 |  |  |  | 0.96 | 0.37 | 0.68 |
| Autism PGS (2026) | above median | CNB speed | 16p11.2 dup |  24 |  | +0.274 (0.169) | 0.12 | 0.37 |  |  |  | 0.2 | 0.81 | 0.018 |
| Autism PGS (2026) | below median | CNB speed | 16p11.2 dup |  15 |  | -0.028 (0.173) | 0.87 | 0.82 |  |  |  | 0.2 | 0.81 | 0.018 |
| Autism PGS (2026) | above median | CNB speed | 22q11.2 del | 136 |  | +0.015 (0.058) | 0.8 | 0.36 |  |  |  | 0.32 | 0.31 | 0.91 |
| Autism PGS (2026) | below median | CNB speed | 22q11.2 del | 151 |  | -0.041 (0.039) | 0.3 | 0.44 |  |  |  | 0.32 | 0.31 | 0.91 |
| Autism PGS (2026) | above median | CNB speed | 22q11.2 dup |  50 |  | +0.084 (0.068) | 0.22 | 0.61 |  |  |  | 0.18 | 0.41 | 0.95 |
| Autism PGS (2026) | below median | CNB speed | 22q11.2 dup |  60 |  | -0.074 (0.150) | 0.63 | 0.6 |  |  |  | 0.18 | 0.41 | 0.95 |
| Autism PGS (2026) | above median | CNB speed | Overgrowth CNVs |  99 |  | +0.012 (0.045) | 0.8 | 0.5 |  |  |  | 0.6 | 0.23 | 0.72 |
| Autism PGS (2026) | below median | CNB speed | Overgrowth CNVs | 114 |  | -0.026 (0.093) | 0.78 | 0.44 |  |  |  | 0.6 | 0.23 | 0.72 |
| Autism PGS (2026) | above median | CNB speed | Undergrowth CNVs | 160 |  | +0.044 (0.053) | 0.41 | 0.43 |  |  |  | 0.16 | 0.17 | 0.24 |
| Autism PGS (2026) | below median | CNB speed | Undergrowth CNVs | 166 |  | -0.038 (0.038) | 0.32 | 0.6 |  |  |  | 0.16 | 0.17 | 0.24 |
| Autism PGS (2026) | above median | CNB speed | All CNVs | 259 |  | +0.031 (0.038) | 0.42 | 0.55 |  |  |  | 0.061 | 0.69 | 0.33 |
| Autism PGS (2026) | below median | CNB speed | All CNVs | 280 |  | -0.035 (0.044) | 0.42 | 0.32 |  |  |  | 0.061 | 0.69 | 0.33 |
| Childhood obesity PGS | above median | Height | 16p11.2 del |  55 | 2503 | +0.403 (0.144) | 0.0074 | 0.14 | +0.256 (0.109) | 0.019 | 0.003 | 0.65 | 0.59 | 0.75 |
| Childhood obesity PGS | below median | Height | 16p11.2 del |  63 | 2507 | +0.325 (0.092) | 7.8e-04 | 0.13 | +0.229 (0.091) | 0.012 | 1 | 0.65 | 0.59 | 0.75 |
| Childhood obesity PGS | above median | Height | 16p11.2 dup |  33 | 2503 | +0.270 (0.126) | 0.04 | 0.7 | +0.143 (0.121) | 0.24 | 9.9e-05 | 0.0095 | 0.8 | 0.052 |
| Childhood obesity PGS | below median | Height | 16p11.2 dup |  34 | 2507 | -0.372 (0.198) | 0.07 | 0.4 | -0.436 (0.178) | 0.014 | 0.014 | 0.0095 | 0.8 | 0.052 |
| Childhood obesity PGS | above median | Height | 22q11.2 del | 134 | 2503 | +0.063 (0.070) | 0.37 | 0.86 | -0.062 (0.068) | 0.36 | 0.27 | 0.9 | 0.28 | 0.81 |
| Childhood obesity PGS | below median | Height | 22q11.2 del | 134 | 2507 | +0.050 (0.082) | 0.54 | 0.019 | -0.025 (0.077) | 0.74 | 0.5 | 0.9 | 0.28 | 0.81 |
| Childhood obesity PGS | above median | Height | 22q11.2 dup |  54 | 2503 | +0.177 (0.109) | 0.11 | 0.21 | +0.030 (0.119) | 0.8 | 2.9e-04 | 0.37 | 0.028 | 0.34 |
| Childhood obesity PGS | below median | Height | 22q11.2 dup |  43 | 2507 | -0.014 (0.143) | 0.92 | 0.2 | -0.074 (0.139) | 0.59 | 0.5 | 0.37 | 0.028 | 0.34 |
| Childhood obesity PGS | above median | Height | Overgrowth CNVs | 109 | 2503 | +0.295 (0.091) | 0.0016 | 0.073 | +0.136 (0.083) | 0.1 | 2.0e-04 | 0.98 | 0.42 | 0.39 |
| Childhood obesity PGS | below median | Height | Overgrowth CNVs | 106 | 2507 | +0.241 (0.078) | 0.0027 | 0.095 | +0.156 (0.078) | 0.045 | 0.005 | 0.98 | 0.42 | 0.39 |
| Childhood obesity PGS | above median | Height | Undergrowth CNVs | 167 | 2503 | +0.112 (0.062) | 0.072 | 0.62 | -0.007 (0.062) | 0.91 | 0.0019 | 0.21 | 0.55 | 0.37 |
| Childhood obesity PGS | below median | Height | Undergrowth CNVs | 168 | 2507 | -0.014 (0.076) | 0.86 | 0.032 | -0.087 (0.072) | 0.23 | 1 | 0.21 | 0.55 | 0.37 |
| Childhood obesity PGS | above median | Height | All CNVs | 276 | 2503 | +0.141 (0.054) | 0.0098 | 0.011 | +0.006 (0.053) | 0.91 | 3.2e-04 | 0.92 | 0.41 | 0.49 |
| Childhood obesity PGS | below median | Height | All CNVs | 274 | 2507 | +0.102 (0.057) | 0.073 | 6.4e-07 | +0.013 (0.057) | 0.82 | 1 | 0.92 | 0.41 | 0.49 |
| Childhood obesity PGS | above median | Head circumference | 16p11.2 del |  46 | 2489 | +0.615 (0.131) | 2.8e-05 | 0.18 | +0.255 (0.121) | 0.035 | 0.015 | 0.2 | 0.72 | 0.19 |
| Childhood obesity PGS | below median | Head circumference | 16p11.2 del |  52 | 2496 | +0.404 (0.103) | 2.6e-04 | 0.19 | +0.063 (0.098) | 0.52 | 0.048 | 0.2 | 0.72 | 0.19 |
| Childhood obesity PGS | above median | Head circumference | 16p11.2 dup |  26 | 2489 | +0.109 (0.151) | 0.48 | 0.067 | -0.254 (0.142) | 0.074 | 0.41 | 0.9 | 0.19 | 0.85 |
| Childhood obesity PGS | below median | Head circumference | 16p11.2 dup |  28 | 2496 | +0.137 (0.213) | 0.53 | 0.3 | -0.187 (0.222) | 0.4 | 0.69 | 0.9 | 0.19 | 0.85 |
| Childhood obesity PGS | above median | Head circumference | 22q11.2 del |  84 | 2489 | +0.163 (0.097) | 0.097 | 0.44 | -0.194 (0.089) | 0.03 | 0.067 | 0.019 | 0.09 | 0.024 |
| Childhood obesity PGS | below median | Head circumference | 22q11.2 del |  87 | 2496 | +0.477 (0.091) | 1.1e-06 | 0.75 | +0.115 (0.095) | 0.23 | 0.13 | 0.019 | 0.09 | 0.024 |
| Childhood obesity PGS | above median | Head circumference | 22q11.2 dup |  43 | 2489 | +0.558 (0.137) | 2.2e-04 | 0.61 | +0.186 (0.134) | 0.17 | 0.35 | 0.49 | 0.43 | 0.87 |
| Childhood obesity PGS | below median | Head circumference | 22q11.2 dup |  29 | 2496 | +0.405 (0.212) | 0.068 | 0.78 | +0.012 (0.185) | 0.95 | 0.052 | 0.49 | 0.43 | 0.87 |
| Childhood obesity PGS | above median | Head circumference | Overgrowth CNVs |  89 | 2489 | +0.587 (0.094) | 1.5e-08 | 0.29 | +0.219 (0.093) | 0.018 | 0.0043 | 0.28 | 0.57 | 0.33 |
| Childhood obesity PGS | below median | Head circumference | Overgrowth CNVs |  81 | 2496 | +0.409 (0.092) | 3.2e-05 | 0.56 | +0.058 (0.088) | 0.51 | 0.12 | 0.28 | 0.57 | 0.33 |
| Childhood obesity PGS | above median | Head circumference | Undergrowth CNVs | 110 | 2489 | +0.148 (0.084) | 0.079 | 0.25 | -0.216 (0.078) | 0.0057 | 0.069 | 0.032 | 0.044 | 0.038 |
| Childhood obesity PGS | below median | Head circumference | Undergrowth CNVs | 115 | 2496 | +0.416 (0.084) | 2.9e-06 | 0.82 | +0.059 (0.088) | 0.5 | 0.55 | 0.032 | 0.044 | 0.038 |
| Childhood obesity PGS | above median | Head circumference | All CNVs | 199 | 2489 | +0.323 (0.064) | 9.5e-07 | 0.099 | -0.048 (0.064) | 0.45 | 0.44 | 0.19 | 0.035 | 0.3 |
| Childhood obesity PGS | below median | Head circumference | All CNVs | 196 | 2496 | +0.412 (0.062) | 2.5e-10 | 0.54 | +0.060 (0.066) | 0.36 | 0.0081 | 0.19 | 0.035 | 0.3 |
| Childhood obesity PGS | above median | SRS | 16p11.2 del |  34 | 2494 | +0.346 (0.293) | 0.25 | 0.29 | +0.274 (0.135) | 0.042 | 0.015 | 0.4 | 0.33 | 0.35 |
| Childhood obesity PGS | below median | SRS | 16p11.2 del |  40 | 2494 | +0.020 (0.301) | 0.95 | 0.96 | -0.012 (0.118) | 0.92 | 0.74 | 0.4 | 0.33 | 0.35 |
| Childhood obesity PGS | below median | SRS | 16p11.2 dup |  15 | 2494 | -0.398 (0.521) | 0.46 | 0.22 | -0.442 (0.267) | 0.098 | 1 | 0.24 | 0.86 | 0.14 |
| Childhood obesity PGS | above median | SRS | 22q11.2 del |  93 | 2494 | +0.216 (0.142) | 0.13 | 0.93 | +0.164 (0.077) | 0.034 | 0.98 | 0.35 | 0.53 | 0.78 |
| Childhood obesity PGS | below median | SRS | 22q11.2 del |  93 | 2494 | +0.415 (0.159) | 0.011 | 0.53 | +0.404 (0.082) | 9.2e-07 | 1 | 0.35 | 0.53 | 0.78 |
| Childhood obesity PGS | above median | SRS | 22q11.2 dup |  39 | 2494 | +0.052 (0.366) | 0.89 | 0.76 | -0.153 (0.163) | 0.35 | 0.011 | 0.42 | 0.49 | 0.18 |
| Childhood obesity PGS | below median | SRS | 22q11.2 dup |  21 | 2494 | -0.440 (0.409) | 0.3 | 0.67 | -0.437 (0.205) | 0.033 | 0.61 | 0.42 | 0.49 | 0.18 |
| Childhood obesity PGS | above median | SRS | Overgrowth CNVs |  73 | 2494 | +0.236 (0.225) | 0.3 | 0.36 | +0.030 (0.108) | 0.78 | 0.033 | 0.12 | 0.9 | 0.014 |
| Childhood obesity PGS | below median | SRS | Overgrowth CNVs |  61 | 2494 | -0.105 (0.244) | 0.67 | 0.83 | -0.098 (0.105) | 0.35 | 0.86 | 0.12 | 0.9 | 0.014 |
| Childhood obesity PGS | above median | SRS | Undergrowth CNVs | 107 | 2494 | +0.247 (0.139) | 0.079 | 0.19 | +0.208 (0.074) | 0.0052 | 0.99 | 0.68 | 0.9 | 0.69 |
| Childhood obesity PGS | below median | SRS | Undergrowth CNVs | 108 | 2494 | +0.342 (0.154) | 0.028 | 0.9 | +0.326 (0.080) | 4.7e-05 | 1 | 0.68 | 0.9 | 0.69 |
| Childhood obesity PGS | above median | SRS | All CNVs | 180 | 2494 | +0.234 (0.120) | 0.053 | 0.7 | +0.138 (0.067) | 0.041 | 1 | 0.55 | 0.62 | 0.052 |
| Childhood obesity PGS | below median | SRS | All CNVs | 169 | 2494 | +0.206 (0.134) | 0.13 | 0.87 | +0.212 (0.070) | 0.0024 | 0.88 | 0.55 | 0.62 | 0.052 |
| Childhood obesity PGS | above median | SCQ | 16p11.2 del |  27 | 504 |  | n/a | n/a | +0.147 (0.178) | 0.41 | 0.5 | 0.16 | 0.79 | 0.34 |
| Childhood obesity PGS | below median | SCQ | 16p11.2 del |  29 | 526 | -0.094 (0.451) | 0.84 | 0.69 | -0.087 (0.161) | 0.59 | 1 | 0.16 | 0.79 | 0.34 |
| Childhood obesity PGS | above median | SCQ | 22q11.2 del |  75 | 504 | +0.312 (0.227) | 0.17 | 0.41 | +0.232 (0.122) | 0.058 | 0.32 | 0.71 | 0.14 | 0.99 |
| Childhood obesity PGS | below median | SCQ | 22q11.2 del |  76 | 526 | +0.459 (0.256) | 0.077 | 0.91 | +0.441 (0.136) | 0.0012 | 0.75 | 0.71 | 0.14 | 0.99 |
| Childhood obesity PGS | above median | SCQ | 22q11.2 dup |  29 | 504 | -0.260 (0.619) | 0.68 | 0.49 | -0.588 (0.222) | 0.0082 | 1 | 0.15 | 0.65 | 0.22 |
| Childhood obesity PGS | below median | SCQ | 22q11.2 dup |  15 | 526 | +1.113 (0.696) | 0.14 | 0.87 | +1.072 (0.231) | 4.4e-06 | 1 | 0.15 | 0.65 | 0.22 |
| Childhood obesity PGS | above median | SCQ | Overgrowth CNVs |  56 | 504 | +0.385 (0.371) | 0.3 | 0.42 | +0.109 (0.165) | 0.51 | 1 | 0.48 | 0.89 | 0.46 |
| Childhood obesity PGS | below median | SCQ | Overgrowth CNVs |  44 | 526 | +0.259 (0.379) | 0.5 | 0.87 | +0.247 (0.154) | 0.11 | 1 | 0.48 | 0.89 | 0.46 |
| Childhood obesity PGS | above median | SCQ | Undergrowth CNVs |  81 | 504 | +0.294 (0.227) | 0.2 | 0.058 | +0.252 (0.121) | 0.037 | 0.35 | 0.64 | 0.24 | 0.82 |
| Childhood obesity PGS | below median | SCQ | Undergrowth CNVs |  87 | 526 | +0.431 (0.250) | 0.088 | 0.86 | +0.398 (0.135) | 0.0034 | 0.91 | 0.64 | 0.24 | 0.82 |
| Childhood obesity PGS | above median | SCQ | All CNVs | 137 | 504 | +0.296 (0.198) | 0.14 | 0.95 | +0.225 (0.127) | 0.077 | 1 | 0.91 | 0.65 | 0.64 |
| Childhood obesity PGS | below median | SCQ | All CNVs | 131 | 526 | +0.367 (0.214) | 0.089 | 0.6 | +0.362 (0.131) | 0.0061 | 1 | 0.91 | 0.65 | 0.64 |
| Childhood obesity PGS | above median | CNB accuracy | 16p11.2 del |  51 |  | -0.030 (0.068) | 0.66 | 0.42 |  |  |  | 0.56 | 0.73 | 0.43 |
| Childhood obesity PGS | below median | CNB accuracy | 16p11.2 del |  50 |  | +0.068 (0.063) | 0.29 | 0.14 |  |  |  | 0.56 | 0.73 | 0.43 |
| Childhood obesity PGS | above median | CNB accuracy | 16p11.2 dup |  20 |  | -0.154 (0.099) | 0.14 | 0.32 |  |  |  | 0.62 | 0.41 | 0.2 |
| Childhood obesity PGS | below median | CNB accuracy | 16p11.2 dup |  19 |  | -0.218 (0.147) | 0.16 | 0.54 |  |  |  | 0.62 | 0.41 | 0.2 |
| Childhood obesity PGS | above median | CNB accuracy | 22q11.2 del | 133 |  | +0.003 (0.046) | 0.96 | 0.45 |  |  |  | 0.37 | 0.002 | 0.032 |
| Childhood obesity PGS | below median | CNB accuracy | 22q11.2 del | 149 |  | -0.055 (0.042) | 0.19 | 0.19 |  |  |  | 0.37 | 0.002 | 0.032 |
| Childhood obesity PGS | above median | CNB accuracy | 22q11.2 dup |  65 |  | +0.125 (0.078) | 0.11 | 0.66 |  |  |  | 0.94 | 0.17 | 0.74 |
| Childhood obesity PGS | below median | CNB accuracy | 22q11.2 dup |  45 |  | +0.104 (0.080) | 0.2 | 0.26 |  |  |  | 0.94 | 0.17 | 0.74 |
| Childhood obesity PGS | above median | CNB accuracy | Overgrowth CNVs | 116 |  | +0.052 (0.055) | 0.34 | 0.25 |  |  |  | 0.86 | 0.29 | 0.66 |
| Childhood obesity PGS | below median | CNB accuracy | Overgrowth CNVs |  95 |  | +0.073 (0.050) | 0.15 | 0.16 |  |  |  | 0.86 | 0.29 | 0.66 |
| Childhood obesity PGS | above median | CNB accuracy | Undergrowth CNVs | 153 |  | -0.002 (0.042) | 0.96 | 0.12 |  |  |  | 0.34 | 0.012 | 0.049 |
| Childhood obesity PGS | below median | CNB accuracy | Undergrowth CNVs | 168 |  | -0.068 (0.040) | 0.09 | 0.21 |  |  |  | 0.34 | 0.012 | 0.049 |
| Childhood obesity PGS | above median | CNB accuracy | All CNVs | 269 |  | +0.012 (0.033) | 0.72 | 0.11 |  |  |  | 0.58 | 0.0029 | 0.2 |
| Childhood obesity PGS | below median | CNB accuracy | All CNVs | 263 |  | -0.021 (0.032) | 0.5 | 0.048 |  |  |  | 0.58 | 0.0029 | 0.2 |
| Childhood obesity PGS | above median | CNB speed | 16p11.2 del |  51 |  | +0.004 (0.076) | 0.95 | 0.74 |  |  |  | 0.91 | 0.63 | 0.61 |
| Childhood obesity PGS | below median | CNB speed | 16p11.2 del |  52 |  | -0.033 (0.084) | 0.69 | 0.54 |  |  |  | 0.91 | 0.63 | 0.61 |
| Childhood obesity PGS | above median | CNB speed | 16p11.2 dup |  20 |  | +0.178 (0.193) | 0.37 | 0.28 |  |  |  | 0.12 | 0.64 | 0.007 |
| Childhood obesity PGS | below median | CNB speed | 16p11.2 dup |  19 |  | -0.140 (0.133) | 0.31 | 0.084 |  |  |  | 0.12 | 0.64 | 0.007 |
| Childhood obesity PGS | above median | CNB speed | 22q11.2 del | 135 |  | -0.071 (0.053) | 0.18 | 0.5 |  |  |  | 0.091 | 0.57 | 0.43 |
| Childhood obesity PGS | below median | CNB speed | 22q11.2 del | 152 |  | +0.037 (0.046) | 0.42 | 0.034 |  |  |  | 0.091 | 0.57 | 0.43 |
| Childhood obesity PGS | above median | CNB speed | 22q11.2 dup |  65 |  | +0.076 (0.150) | 0.61 | 0.93 |  |  |  | 0.93 | 0.69 | 0.79 |
| Childhood obesity PGS | below median | CNB speed | 22q11.2 dup |  45 |  | +0.055 (0.098) | 0.58 | 0.072 |  |  |  | 0.93 | 0.69 | 0.79 |
| Childhood obesity PGS | above median | CNB speed | Overgrowth CNVs | 116 |  | +0.073 (0.091) | 0.42 | 0.78 |  |  |  | 0.91 | 0.92 | 0.68 |
| Childhood obesity PGS | below median | CNB speed | Overgrowth CNVs |  97 |  | +0.006 (0.062) | 0.93 | 0.13 |  |  |  | 0.91 | 0.92 | 0.68 |
| Childhood obesity PGS | above median | CNB speed | Undergrowth CNVs | 155 |  | -0.031 (0.050) | 0.53 | 0.97 |  |  |  | 0.39 | 0.87 | 0.75 |
| Childhood obesity PGS | below median | CNB speed | Undergrowth CNVs | 171 |  | +0.029 (0.043) | 0.49 | 0.095 |  |  |  | 0.39 | 0.87 | 0.75 |
| Childhood obesity PGS | above median | CNB speed | All CNVs | 271 |  | -0.018 (0.047) | 0.7 | 0.99 |  |  |  | 0.5 | 0.86 | 0.76 |
| Childhood obesity PGS | below median | CNB speed | All CNVs | 268 |  | +0.022 (0.036) | 0.53 | 0.015 |  |  |  | 0.5 | 0.86 | 0.76 |
| Height PGS | above median | Height | 16p11.2 del |  54 | 2504 | +0.535 (0.114) | 2.3e-05 | 0.043 | +0.397 (0.092) | 1.7e-05 | 0.0052 | 0.02 | 0.8 | 0.098 |
| Height PGS | below median | Height | 16p11.2 del |  64 | 2506 | +0.165 (0.101) | 0.11 | 0.33 | +0.036 (0.084) | 0.67 | 0.65 | 0.02 | 0.8 | 0.098 |
| Height PGS | above median | Height | 16p11.2 dup |  36 | 2504 | +0.237 (0.166) | 0.16 | 0.6 | +0.161 (0.146) | 0.27 | 0.36 | 0.08 | 4.1e-05 | 0.14 |
| Height PGS | below median | Height | 16p11.2 dup |  31 | 2506 | -0.139 (0.126) | 0.28 | 0.24 | -0.267 (0.117) | 0.023 | 0.99 | 0.08 | 4.1e-05 | 0.14 |
| Height PGS | above median | Height | 22q11.2 del | 136 | 2504 | +0.101 (0.077) | 0.19 | 0.068 | -0.006 (0.068) | 0.93 | 0.99 | 0.52 | 7.0e-13 | 0.24 |
| Height PGS | below median | Height | 22q11.2 del | 132 | 2506 | +0.031 (0.062) | 0.62 | 0.12 | -0.094 (0.062) | 0.13 | 1 | 0.52 | 7.0e-13 | 0.24 |
| Height PGS | above median | Height | 22q11.2 dup |  39 | 2504 | +0.217 (0.109) | 0.055 | 0.19 | +0.091 (0.114) | 0.43 | 0.24 | 0.42 | 3.4e-06 | 0.88 |
| Height PGS | below median | Height | 22q11.2 dup |  58 | 2506 | +0.101 (0.130) | 0.44 | 0.037 | +0.024 (0.127) | 0.85 | 0.99 | 0.42 | 3.4e-06 | 0.88 |
| Height PGS | above median | Height | Overgrowth CNVs |  93 | 2504 | +0.402 (0.081) | 3.7e-06 | 0.012 | +0.267 (0.074) | 3.2e-04 | 0.0038 | 0.075 | 5.5e-04 | 0.42 |
| Height PGS | below median | Height | Overgrowth CNVs | 122 | 2506 | +0.166 (0.078) | 0.035 | 0.088 | +0.058 (0.072) | 0.42 | 0.65 | 0.075 | 5.5e-04 | 0.42 |
| Height PGS | above median | Height | Undergrowth CNVs | 172 | 2504 | +0.121 (0.069) | 0.083 | 0.15 | +0.018 (0.063) | 0.77 | 0.99 | 0.21 | < 1e-16 | 0.11 |
| Height PGS | below median | Height | Undergrowth CNVs | 163 | 2506 | -0.005 (0.056) | 0.92 | 0.59 | -0.116 (0.056) | 0.037 | 1 | 0.21 | < 1e-16 | 0.11 |
| Height PGS | above median | Height | All CNVs | 265 | 2504 | +0.223 (0.056) | 9.0e-05 | 1.3e-04 | +0.094 (0.052) | 0.072 | 0.0048 | 0.04 | < 1e-16 | 0.31 |
| Height PGS | below median | Height | All CNVs | 285 | 2506 | +0.030 (0.048) | 0.54 | 8.8e-05 | -0.079 (0.047) | 0.094 | 0.62 | 0.04 | < 1e-16 | 0.31 |
| Height PGS | above median | Head circumference | 16p11.2 del |  42 | 2490 | +0.516 (0.118) | 9.1e-05 | 0.12 | +0.147 (0.108) | 0.17 | 0.25 | 0.72 | 0.82 | 0.32 |
| Height PGS | below median | Head circumference | 16p11.2 del |  56 | 2495 | +0.455 (0.112) | 1.6e-04 | 0.18 | +0.111 (0.104) | 0.28 | 0.0016 | 0.72 | 0.82 | 0.32 |
| Height PGS | above median | Head circumference | 16p11.2 dup |  28 | 2490 | +0.091 (0.175) | 0.61 | 0.22 | -0.252 (0.202) | 0.21 | 0.24 | 0.81 | 0.24 | 0.28 |
| Height PGS | below median | Head circumference | 16p11.2 dup |  26 | 2495 | +0.154 (0.152) | 0.32 | 0.61 | -0.207 (0.144) | 0.15 | 0.0017 | 0.81 | 0.24 | 0.28 |
| Height PGS | above median | Head circumference | 22q11.2 del |  90 | 2490 | +0.216 (0.091) | 0.019 | 0.42 | -0.162 (0.086) | 0.061 | 0.025 | 0.081 | 0.22 | 0.028 |
| Height PGS | below median | Head circumference | 22q11.2 del |  81 | 2495 | +0.448 (0.100) | 2.7e-05 | 0.39 | +0.096 (0.095) | 0.31 | 0.78 | 0.081 | 0.22 | 0.028 |
| Height PGS | above median | Head circumference | 22q11.2 dup |  29 | 2490 | +0.528 (0.107) | 3.9e-05 | 0.094 | +0.145 (0.145) | 0.32 | 0.024 | 0.85 | 0.0016 | 0.86 |
| Height PGS | below median | Head circumference | 22q11.2 dup |  43 | 2495 | +0.581 (0.212) | 0.0093 | 0.44 | +0.196 (0.177) | 0.27 | 0.26 | 0.85 | 0.0016 | 0.86 |
| Height PGS | above median | Head circumference | Overgrowth CNVs |  71 | 2490 | +0.522 (0.083) | 2.6e-08 | 0.84 | +0.150 (0.089) | 0.092 | 0.022 | 0.94 | 0.018 | 0.64 |
| Height PGS | below median | Head circumference | Overgrowth CNVs |  99 | 2495 | +0.486 (0.101) | 5.7e-06 | 0.74 | +0.132 (0.091) | 0.15 | 0.0033 | 0.94 | 0.018 | 0.64 |
| Height PGS | above median | Head circumference | Undergrowth CNVs | 118 | 2490 | +0.193 (0.080) | 0.017 | 0.48 | -0.181 (0.080) | 0.024 | 0.024 | 0.26 | 0.51 | 0.048 |
| Height PGS | below median | Head circumference | Undergrowth CNVs | 107 | 2495 | +0.352 (0.088) | 1.2e-04 | 0.11 | -0.006 (0.081) | 0.94 | 0.0017 | 0.26 | 0.51 | 0.048 |
| Height PGS | above median | Head circumference | All CNVs | 189 | 2490 | +0.335 (0.059) | 5.9e-08 | 0.24 | -0.039 (0.063) | 0.54 | 0.23 | 0.53 | 0.029 | 0.28 |
| Height PGS | below median | Head circumference | All CNVs | 206 | 2495 | +0.419 (0.066) | 1.6e-09 | 0.69 | +0.061 (0.064) | 0.35 | 0.0039 | 0.53 | 0.029 | 0.28 |
| Height PGS | above median | SRS | 16p11.2 del |  35 | 2490 | +0.274 (0.242) | 0.27 | 0.4 | +0.228 (0.106) | 0.032 | 0.99 | 0.39 | 0.73 | 0.26 |
| Height PGS | below median | SRS | 16p11.2 del |  39 | 2498 | -0.127 (0.385) | 0.74 | 0.52 | -0.060 (0.147) | 0.69 | 0.56 | 0.39 | 0.73 | 0.26 |
| Height PGS | above median | SRS | 16p11.2 dup |  18 | 2490 | -0.201 (0.615) | 0.75 | 0.61 | -0.058 (0.258) | 0.82 | 0.19 | 0.36 | 0.26 | 0.39 |
| Height PGS | above median | SRS | 22q11.2 del |  92 | 2490 | +0.364 (0.161) | 0.026 | 0.91 | +0.371 (0.081) | 5.3e-06 | 1 | 0.72 | 0.27 | 0.81 |
| Height PGS | below median | SRS | 22q11.2 del |  94 | 2498 | +0.291 (0.139) | 0.039 | 0.65 | +0.242 (0.077) | 0.0018 | 0.72 | 0.72 | 0.27 | 0.81 |
| Height PGS | above median | SRS | 22q11.2 dup |  23 | 2490 | -0.367 (0.456) | 0.43 | 0.51 | -0.356 (0.196) | 0.069 | 0.23 | 0.84 | 0.21 | 0.79 |
| Height PGS | below median | SRS | 22q11.2 dup |  37 | 2498 | -0.123 (0.335) | 0.72 | 0.76 | -0.158 (0.171) | 0.36 | 0.43 | 0.84 | 0.21 | 0.79 |
| Height PGS | above median | SRS | Overgrowth CNVs |  58 | 2490 | +0.105 (0.201) | 0.6 | 0.77 | +0.091 (0.094) | 0.33 | 0.61 | 0.5 | 0.66 | 0.25 |
| Height PGS | below median | SRS | Overgrowth CNVs |  76 | 2498 | -0.106 (0.251) | 0.67 | 0.62 | -0.101 (0.115) | 0.38 | 0.5 | 0.5 | 0.66 | 0.25 |
| Height PGS | above median | SRS | Undergrowth CNVs | 110 | 2490 | +0.309 (0.157) | 0.052 | 1 | +0.318 (0.079) | 6.4e-05 | 0.31 | 0.86 | 0.46 | 0.95 |
| Height PGS | below median | SRS | Undergrowth CNVs | 105 | 2498 | +0.321 (0.134) | 0.018 | 0.38 | +0.258 (0.075) | 6.0e-04 | 0.61 | 0.86 | 0.46 | 0.95 |
| Height PGS | above median | SRS | All CNVs | 168 | 2490 | +0.256 (0.125) | 0.042 | 0.73 | +0.247 (0.067) | 2.1e-04 | 1 | 0.48 | 0.25 | 0.44 |
| Height PGS | below median | SRS | All CNVs | 181 | 2498 | +0.201 (0.125) | 0.11 | 0.54 | +0.130 (0.069) | 0.059 | 0.81 | 0.48 | 0.25 | 0.44 |
| Height PGS | above median | SCQ | 16p11.2 del |  27 | 512 | +0.316 (0.404) | 0.44 | 0.96 | -0.125 (0.196) | 0.53 | 1 | 0.64 | 0.11 | 0.29 |
| Height PGS | below median | SCQ | 16p11.2 del |  29 | 518 | +0.272 (0.525) | 0.61 | 0.54 | +0.051 (0.184) | 0.78 | 0.87 | 0.64 | 0.11 | 0.29 |
| Height PGS | above median | SCQ | 22q11.2 del |  74 | 512 | +0.462 (0.256) | 0.075 | 0.47 | +0.319 (0.135) | 0.018 | 0.98 | 0.79 | 0.73 | 0.87 |
| Height PGS | below median | SCQ | 22q11.2 del |  77 | 518 | +0.323 (0.222) | 0.15 | 0.91 | +0.359 (0.122) | 0.0033 | 0.5 | 0.79 | 0.73 | 0.87 |
| Height PGS | above median | SCQ | 22q11.2 dup |  17 | 512 | -0.619 (0.574) | 0.3 | 0.56 | -0.738 (0.211) | 5.1e-04 | 1 | 0.17 | 0.14 | 0.51 |
| Height PGS | below median | SCQ | 22q11.2 dup |  27 | 518 | +0.802 (0.652) | 0.23 | 0.75 | +1.081 (0.242) | 9.7e-06 | 0.75 | 0.17 | 0.14 | 0.51 |
| Height PGS | above median | SCQ | Overgrowth CNVs |  44 | 512 | +0.193 (0.325) | 0.56 | 0.61 | -0.127 (0.156) | 0.42 | 0.5 | 0.44 | 0.89 | 0.84 |
| Height PGS | below median | SCQ | Overgrowth CNVs |  56 | 518 | +0.446 (0.405) | 0.28 | 0.9 | +0.564 (0.172) | 0.0011 | 0.72 | 0.44 | 0.89 | 0.84 |
| Height PGS | above median | SCQ | Undergrowth CNVs |  85 | 512 | +0.380 (0.258) | 0.15 | 0.93 | +0.301 (0.139) | 0.031 | 1 | 0.85 | 0.71 | 0.86 |
| Height PGS | below median | SCQ | Undergrowth CNVs |  83 | 518 | +0.303 (0.218) | 0.17 | 0.8 | +0.328 (0.121) | 0.007 | 1 | 0.85 | 0.71 | 0.86 |
| Height PGS | above median | SCQ | All CNVs | 129 | 512 | +0.293 (0.204) | 0.15 | 0.55 | +0.232 (0.127) | 0.069 | 0.5 | 0.9 | 0.92 | 0.99 |
| Height PGS | below median | SCQ | All CNVs | 139 | 518 | +0.322 (0.201) | 0.11 | 0.87 | +0.355 (0.130) | 0.0065 | 0.5 | 0.9 | 0.92 | 0.99 |
| Height PGS | above median | CNB accuracy | 16p11.2 del |  50 |  | -0.052 (0.056) | 0.36 | 0.32 |  |  |  | 0.069 | 0.58 | 0.17 |
| Height PGS | below median | CNB accuracy | 16p11.2 del |  51 |  | +0.128 (0.074) | 0.092 | 0.23 |  |  |  | 0.069 | 0.58 | 0.17 |
| Height PGS | above median | CNB accuracy | 16p11.2 dup |  21 |  | -0.282 (0.119) | 0.029 | 0.5 |  |  |  | 0.13 | 0.19 | 0.42 |
| Height PGS | below median | CNB accuracy | 16p11.2 dup |  18 |  | +0.041 (0.099) | 0.69 | 0.24 |  |  |  | 0.13 | 0.19 | 0.42 |
| Height PGS | above median | CNB accuracy | 22q11.2 del | 156 |  | -0.011 (0.045) | 0.81 | 0.14 |  |  |  | 0.49 | 0.47 | 0.22 |
| Height PGS | below median | CNB accuracy | 22q11.2 del | 126 |  | -0.053 (0.042) | 0.21 | 0.38 |  |  |  | 0.49 | 0.47 | 0.22 |
| Height PGS | above median | CNB accuracy | 22q11.2 dup |  51 |  | +0.069 (0.083) | 0.41 | 0.84 |  |  |  | 0.94 | 0.61 | 0.81 |
| Height PGS | below median | CNB accuracy | 22q11.2 dup |  59 |  | +0.084 (0.084) | 0.32 | 0.22 |  |  |  | 0.94 | 0.61 | 0.81 |
| Height PGS | above median | CNB accuracy | Overgrowth CNVs | 101 |  | +0.003 (0.049) | 0.96 | 0.45 |  |  |  | 0.18 | 0.43 | 0.28 |
| Height PGS | below median | CNB accuracy | Overgrowth CNVs | 110 |  | +0.093 (0.056) | 0.097 | 0.098 |  |  |  | 0.18 | 0.43 | 0.28 |
| Height PGS | above median | CNB accuracy | Undergrowth CNVs | 177 |  | -0.026 (0.043) | 0.55 | 0.23 |  |  |  | 0.6 | 0.71 | 0.26 |
| Height PGS | below median | CNB accuracy | Undergrowth CNVs | 144 |  | -0.047 (0.038) | 0.22 | 0.67 |  |  |  | 0.6 | 0.71 | 0.26 |
| Height PGS | above median | CNB accuracy | All CNVs | 278 |  | -0.018 (0.032) | 0.57 | 0.064 |  |  |  | 0.6 | 0.89 | 0.81 |
| Height PGS | below median | CNB accuracy | All CNVs | 254 |  | +0.004 (0.032) | 0.9 | 0.099 |  |  |  | 0.6 | 0.89 | 0.81 |
| Height PGS | above median | CNB speed | 16p11.2 del |  50 |  | +0.004 (0.086) | 0.96 | 0.85 |  |  |  | 0.77 | 0.26 | 0.83 |
| Height PGS | below median | CNB speed | 16p11.2 del |  53 |  | -0.046 (0.069) | 0.51 | 0.66 |  |  |  | 0.77 | 0.26 | 0.83 |
| Height PGS | above median | CNB speed | 16p11.2 dup |  21 |  | -0.220 (0.127) | 0.1 | 0.89 |  |  |  | 0.066 | 0.37 | 0.06 |
| Height PGS | below median | CNB speed | 16p11.2 dup |  18 |  | +0.195 (0.149) | 0.21 | 0.4 |  |  |  | 0.066 | 0.37 | 0.06 |
| Height PGS | above median | CNB speed | 22q11.2 del | 157 |  | +0.038 (0.043) | 0.37 | 0.054 |  |  |  | 0.06 | 0.92 | 0.06 |
| Height PGS | below median | CNB speed | 22q11.2 del | 130 |  | -0.074 (0.054) | 0.18 | 0.89 |  |  |  | 0.06 | 0.92 | 0.06 |
| Height PGS | above median | CNB speed | 22q11.2 dup |  51 |  | +0.050 (0.170) | 0.77 | 0.97 |  |  |  | 0.99 | 0.86 | 0.96 |
| Height PGS | below median | CNB speed | 22q11.2 dup |  59 |  | +0.048 (0.098) | 0.62 | 0.14 |  |  |  | 0.99 | 0.86 | 0.96 |
| Height PGS | above median | CNB speed | Overgrowth CNVs | 101 |  | +0.029 (0.095) | 0.76 | 0.59 |  |  |  | 0.83 | 0.78 | 1 |
| Height PGS | below median | CNB speed | Overgrowth CNVs | 112 |  | +0.002 (0.058) | 0.97 | 0.7 |  |  |  | 0.83 | 0.78 | 1 |
| Height PGS | above median | CNB speed | Undergrowth CNVs | 178 |  | +0.032 (0.041) | 0.43 | 0.29 |  |  |  | 0.3 | 0.65 | 0.19 |
| Height PGS | below median | CNB speed | Undergrowth CNVs | 148 |  | -0.044 (0.051) | 0.4 | 0.72 |  |  |  | 0.3 | 0.65 | 0.19 |
| Height PGS | above median | CNB speed | All CNVs | 279 |  | +0.025 (0.043) | 0.57 | 0.58 |  |  |  | 0.55 | 0.92 | 0.25 |
| Height PGS | below median | CNB speed | All CNVs | 260 |  | -0.029 (0.038) | 0.44 | 0.92 |  |  |  | 0.55 | 0.92 | 0.25 |
| Sex | female | Height | 16p11.2 del | 161 | 2917 | +0.166 (0.075) | 0.029 | 0.32 | +0.045 (0.068) | 0.51 | 0.99 | 0.86 |  |  |
| Sex | male | Height | 16p11.2 del | 197 | 2772 | +0.089 (0.062) | 0.15 | 0.14 | -0.006 (0.055) | 0.91 | 0.098 | 0.86 |  |  |
| Sex | female | Height | 16p11.2 dup |  94 | 2917 | +0.073 (0.110) | 0.51 | 0.072 | +0.013 (0.085) | 0.88 | 0.38 | 0.6 |  |  |
| Sex | male | Height | 16p11.2 dup | 103 | 2772 | -0.018 (0.087) | 0.83 | 0.78 | -0.144 (0.080) | 0.07 | 1.5e-05 | 0.6 |  |  |
| Sex | female | Height | 22q11.2 del | 180 | 2917 | +0.086 (0.075) | 0.25 | 0.02 | -0.010 (0.067) | 0.88 | 0.5 | 0.8 |  |  |
| Sex | male | Height | 22q11.2 del | 175 | 2772 | +0.062 (0.051) | 0.23 | 0.34 | -0.075 (0.058) | 0.2 | 0.12 | 0.8 |  |  |
| Sex | female | Height | 22q11.2 dup |  54 | 2917 | +0.033 (0.138) | 0.81 | 0.23 | -0.085 (0.152) | 0.58 | 0.78 | 0.054 |  |  |
| Sex | male | Height | 22q11.2 dup |  62 | 2772 | +0.233 (0.085) | 0.0084 | 0.72 | +0.070 (0.094) | 0.46 | 1.5e-04 | 0.054 |  |  |
| Sex | female | Height | Overgrowth CNVs | 215 | 2917 | +0.127 (0.066) | 0.055 | 0.11 | +0.013 (0.064) | 0.84 | 0.99 | 0.71 |  |  |
| Sex | male | Height | Overgrowth CNVs | 259 | 2772 | +0.121 (0.052) | 0.021 | 0.1 | +0.010 (0.050) | 0.84 | 0.052 | 0.71 |  |  |
| Sex | female | Height | Undergrowth CNVs | 274 | 2917 | +0.068 (0.062) | 0.27 | 0.084 | -0.015 (0.054) | 0.78 | 0.94 | 0.58 |  |  |
| Sex | male | Height | Undergrowth CNVs | 278 | 2772 | +0.036 (0.045) | 0.43 | 0.083 | -0.096 (0.051) | 0.059 | 0.039 | 0.58 |  |  |
| Sex | female | Height | All CNVs | 489 | 2917 | +0.071 (0.046) | 0.12 | 4.9e-05 | -0.026 (0.046) | 0.57 | 0.99 | 0.65 |  |  |
| Sex | male | Height | All CNVs | 537 | 2772 | +0.081 (0.035) | 0.021 | 6.7e-04 | -0.054 (0.041) | 0.18 | 0.068 | 0.65 |  |  |
| Sex | female | Head circumference | 16p11.2 del |  99 | 2910 | +0.547 (0.104) | 8.1e-07 | 0.53 | +0.132 (0.101) | 0.19 | 0.032 | 0.2 |  |  |
| Sex | male | Head circumference | 16p11.2 del | 112 | 2752 | +0.421 (0.063) | 1.4e-09 | 0.18 | +0.129 (0.065) | 0.048 | 0.99 | 0.2 |  |  |
| Sex | female | Head circumference | 16p11.2 dup |  57 | 2910 | +0.297 (0.137) | 0.035 | 0.21 | -0.099 (0.132) | 0.45 | 1 | 0.73 |  |  |
| Sex | male | Head circumference | 16p11.2 dup |  66 | 2752 | +0.222 (0.094) | 0.021 | 0.19 | -0.087 (0.100) | 0.39 | 0.3 | 0.73 |  |  |
| Sex | female | Head circumference | 22q11.2 del | 119 | 2910 | +0.340 (0.093) | 4.0e-04 | 0.3 | -0.069 (0.090) | 0.45 | 0.99 | 0.89 |  |  |
| Sex | male | Head circumference | 22q11.2 del | 123 | 2752 | +0.352 (0.066) | 5.7e-07 | 0.75 | +0.031 (0.066) | 0.64 | 0.2 | 0.89 |  |  |
| Sex | female | Head circumference | 22q11.2 dup |  39 | 2910 | +0.608 (0.215) | 0.0077 | 0.66 | +0.178 (0.199) | 0.37 | 0.042 | 0.43 |  |  |
| Sex | male | Head circumference | 22q11.2 dup |  52 | 2752 | +0.407 (0.103) | 2.6e-04 | 0.061 | +0.073 (0.096) | 0.45 | 1 | 0.43 |  |  |
| Sex | female | Head circumference | Overgrowth CNVs | 138 | 2910 | +0.543 (0.092) | 3.1e-08 | 0.4 | +0.126 (0.091) | 0.17 | 0.033 | 0.072 |  |  |
| Sex | male | Head circumference | Overgrowth CNVs | 164 | 2752 | +0.411 (0.054) | 1.9e-12 | 0.86 | +0.106 (0.056) | 0.059 | 0.13 | 0.072 |  |  |
| Sex | female | Head circumference | Undergrowth CNVs | 176 | 2910 | +0.323 (0.078) | 5.5e-05 | 0.42 | -0.086 (0.077) | 0.26 | 1 | 0.99 |  |  |
| Sex | male | Head circumference | Undergrowth CNVs | 189 | 2752 | +0.323 (0.054) | 1.2e-08 | 0.54 | +0.006 (0.057) | 0.92 | 0.83 | 0.99 |  |  |
| Sex | female | Head circumference | All CNVs | 314 | 2910 | +0.413 (0.060) | 2.3e-11 | 0.072 | -0.007 (0.062) | 0.91 | 0.99 | 0.58 |  |  |
| Sex | male | Head circumference | All CNVs | 353 | 2752 | +0.367 (0.038) | < 1e-16 | 0.82 | +0.057 (0.044) | 0.19 | 0.41 | 0.58 |  |  |
| Sex | female | SRS | 16p11.2 del |  38 | 2796 | +0.223 (0.301) | 0.46 | 0.53 | +0.235 (0.115) | 0.042 | 0.17 | 0.8 |  |  |
| Sex | male | SRS | 16p11.2 del |  41 | 2673 | +0.182 (0.261) | 0.49 | 0.85 | +0.205 (0.118) | 0.083 | 0.28 | 0.8 |  |  |
| Sex | female | SRS | 16p11.2 dup |  18 | 2796 | +0.585 (0.505) | 0.27 | 0.73 | +0.601 (0.223) | 0.007 | 0.14 | 0.85 |  |  |
| Sex | male | SRS | 16p11.2 dup |  16 | 2673 | +0.610 (0.255) | 0.033 | 0.48 | +0.631 (0.182) | 5.2e-04 | 0.0065 | 0.85 |  |  |
| Sex | female | SRS | 22q11.2 del | 109 | 2796 | +0.307 (0.157) | 0.054 | 0.25 | +0.249 (0.079) | 0.0017 | 0.99 | 0.61 |  |  |
| Sex | male | SRS | 22q11.2 del |  98 | 2673 | +0.204 (0.131) | 0.12 | 0.62 | +0.219 (0.072) | 0.0023 | 0.3 | 0.61 |  |  |
| Sex | female | SRS | 22q11.2 dup |  34 | 2796 | +0.207 (0.444) | 0.64 | 0.35 | +0.142 (0.218) | 0.52 | 0.018 | 0.44 |  |  |
| Sex | male | SRS | 22q11.2 dup |  32 | 2673 | -0.179 (0.287) | 0.54 | 0.72 | -0.131 (0.146) | 0.37 | 0.28 | 0.44 |  |  |
| Sex | female | SRS | Overgrowth CNVs |  72 | 2796 | +0.229 (0.238) | 0.34 | 0.38 | +0.159 (0.104) | 0.13 | 0.26 | 0.089 |  |  |
| Sex | male | SRS | Overgrowth CNVs |  73 | 2673 | +0.004 (0.192) | 0.98 | 0.87 | +0.030 (0.096) | 0.76 | 0.3 | 0.089 |  |  |
| Sex | female | SRS | Undergrowth CNVs | 127 | 2796 | +0.380 (0.154) | 0.015 | 0.13 | +0.318 (0.075) | 2.5e-05 | 0.99 | 0.45 |  |  |
| Sex | male | SRS | Undergrowth CNVs | 114 | 2673 | +0.239 (0.119) | 0.047 | 0.46 | +0.254 (0.068) | 2.1e-04 | 0.3 | 0.45 |  |  |
| Sex | female | SRS | All CNVs | 199 | 2796 | +0.331 (0.131) | 0.012 | 0.34 | +0.241 (0.067) | 3.4e-04 | 0.99 | 0.048 |  |  |
| Sex | male | SRS | All CNVs | 187 | 2673 | +0.160 (0.103) | 0.12 | 0.56 | +0.172 (0.061) | 0.0046 | 0.25 | 0.048 |  |  |
| Sex | female | SCQ | 16p11.2 del |  28 | 608 | +0.720 (0.430) | 0.11 | 0.55 | +0.695 (0.148) | 3.0e-06 | 0.84 | 0.39 |  |  |
| Sex | male | SCQ | 16p11.2 del |  32 | 526 | +0.163 (0.474) | 0.73 | 0.65 | +0.192 (0.181) | 0.29 | 1 | 0.39 |  |  |
| Sex | female | SCQ | 22q11.2 del |  91 | 608 | +0.275 (0.203) | 0.18 | 0.23 | +0.134 (0.112) | 0.23 | 0.73 | 0.74 |  |  |
| Sex | male | SCQ | 22q11.2 del |  78 | 526 | +0.178 (0.235) | 0.45 | 0.63 | +0.180 (0.126) | 0.15 | 0.83 | 0.74 |  |  |
| Sex | female | SCQ | 22q11.2 dup |  27 | 608 | +0.173 (0.609) | 0.78 | 0.57 | -0.273 (0.219) | 0.21 | 0.55 | 0.92 |  |  |
| Sex | male | SCQ | 22q11.2 dup |  23 | 526 | +0.146 (0.538) | 0.79 | 0.6 | +0.220 (0.198) | 0.27 | 0.5 | 0.92 |  |  |
| Sex | female | SCQ | Overgrowth CNVs |  55 | 608 | +0.533 (0.360) | 0.14 | 0.29 | +0.178 (0.159) | 0.26 | 0.99 | 0.35 |  |  |
| Sex | male | SCQ | Overgrowth CNVs |  55 | 526 | +0.244 (0.343) | 0.48 | 0.49 | +0.208 (0.154) | 0.18 | 0.018 | 0.35 |  |  |
| Sex | female | SCQ | Undergrowth CNVs | 104 | 608 | +0.304 (0.202) | 0.14 | 0.12 | +0.142 (0.112) | 0.21 | 1 | 0.83 |  |  |
| Sex | male | SCQ | Undergrowth CNVs |  86 | 526 | +0.287 (0.222) | 0.2 | 0.88 | +0.324 (0.122) | 0.0083 | 0.72 | 0.83 |  |  |
| Sex | female | SCQ | All CNVs | 159 | 608 | +0.376 (0.182) | 0.04 | 0.52 | +0.224 (0.115) | 0.051 | 1 | 0.33 |  |  |
| Sex | male | SCQ | All CNVs | 141 | 526 | +0.253 (0.189) | 0.18 | 0.55 | +0.287 (0.123) | 0.02 | 0.91 | 0.33 |  |  |
| Sex | female | CNB accuracy | 16p11.2 del |  68 |  | +0.079 (0.068) | 0.25 | 0.92 |  |  |  | 0.43 |  |  |
| Sex | male | CNB accuracy | 16p11.2 del |  72 |  | +0.001 (0.052) | 0.99 | 0.27 |  |  |  | 0.43 |  |  |
| Sex | female | CNB accuracy | 16p11.2 dup |  30 |  | -0.124 (0.110) | 0.27 | 0.25 |  |  |  | 0.71 |  |  |
| Sex | male | CNB accuracy | 16p11.2 dup |  25 |  | -0.073 (0.119) | 0.55 | 0.015 |  |  |  | 0.71 |  |  |
| Sex | female | CNB accuracy | 22q11.2 del | 190 |  | +0.012 (0.042) | 0.78 | 0.12 |  |  |  | 0.49 |  |  |
| Sex | male | CNB accuracy | 22q11.2 del | 154 |  | -0.030 (0.041) | 0.46 | 0.45 |  |  |  | 0.49 |  |  |
| Sex | female | CNB accuracy | 22q11.2 dup |  64 |  | +0.291 (0.096) | 0.0037 | 0.25 |  |  |  | 0.021 |  |  |
| Sex | male | CNB accuracy | 22q11.2 dup |  72 |  | -0.012 (0.062) | 0.84 | 0.75 |  |  |  | 0.021 |  |  |
| Sex | female | CNB accuracy | Overgrowth CNVs | 132 |  | +0.158 (0.057) | 0.0064 | 0.36 |  |  |  | 0.064 |  |  |
| Sex | male | CNB accuracy | Overgrowth CNVs | 144 |  | -0.002 (0.042) | 0.95 | 0.59 |  |  |  | 0.064 |  |  |
| Sex | female | CNB accuracy | Undergrowth CNVs | 220 |  | -0.011 (0.040) | 0.77 | 0.26 |  |  |  | 0.61 |  |  |
| Sex | male | CNB accuracy | Undergrowth CNVs | 179 |  | -0.033 (0.039) | 0.4 | 0.34 |  |  |  | 0.61 |  |  |
| Sex | female | CNB accuracy | All CNVs | 352 |  | +0.044 (0.033) | 0.18 | 0.19 |  |  |  | 0.26 |  |  |
| Sex | male | CNB accuracy | All CNVs | 323 |  | -0.019 (0.028) | 0.5 | 0.35 |  |  |  | 0.26 |  |  |
| Sex | female | CNB speed | 16p11.2 del |  69 |  | +0.025 (0.082) | 0.76 | 0.26 |  |  |  | 0.38 |  |  |
| Sex | male | CNB speed | 16p11.2 del |  73 |  | -0.045 (0.047) | 0.34 | 0.14 |  |  |  | 0.38 |  |  |
| Sex | female | CNB speed | 16p11.2 dup |  30 |  | -0.063 (0.174) | 0.72 | 0.11 |  |  |  | 0.6 |  |  |
| Sex | male | CNB speed | 16p11.2 dup |  25 |  | +0.047 (0.127) | 0.71 | 0.15 |  |  |  | 0.6 |  |  |
| Sex | female | CNB speed | 22q11.2 del | 191 |  | -0.077 (0.057) | 0.18 | 0.082 |  |  |  | 0.089 |  |  |
| Sex | male | CNB speed | 22q11.2 del | 159 |  | +0.017 (0.038) | 0.65 | 0.21 |  |  |  | 0.089 |  |  |
| Sex | female | CNB speed | 22q11.2 dup |  64 |  | -0.067 (0.196) | 0.73 | 0.78 |  |  |  | 0.46 |  |  |
| Sex | male | CNB speed | 22q11.2 dup |  72 |  | +0.038 (0.056) | 0.5 | 0.28 |  |  |  | 0.46 |  |  |
| Sex | female | CNB speed | Overgrowth CNVs | 133 |  | +0.010 (0.099) | 0.92 | 0.84 |  |  |  | 0.71 |  |  |
| Sex | male | CNB speed | Overgrowth CNVs | 145 |  | +0.005 (0.038) | 0.89 | 0.094 |  |  |  | 0.71 |  |  |
| Sex | female | CNB speed | Undergrowth CNVs | 221 |  | -0.080 (0.055) | 0.14 | 0.19 |  |  |  | 0.072 |  |  |
| Sex | male | CNB speed | Undergrowth CNVs | 184 |  | +0.021 (0.036) | 0.57 | 0.7 |  |  |  | 0.072 |  |  |
| Sex | female | CNB speed | All CNVs | 354 |  | -0.056 (0.050) | 0.26 | 0.66 |  |  |  | 0.1 |  |  |
| Sex | male | CNB speed | All CNVs | 329 |  | +0.014 (0.026) | 0.61 | 0.31 |  |  |  | 0.1 |  |  |

- Table S5: the effect of BMI with sex and polygenic scores as covariates (the confounding check). For each trait and group, the carrier BMI slope (GAM y ~ BMI + s(age) [+ locus]) and the CNV x BMI shift versus unaffected controls (model as in Tables S1-S2) are refit (i) unadjusted, (ii) adding sex and its interaction with BMI, (iii) adding one PGS and its interaction with BMI (childhood obesity, adult BMI, height, or autism; the latter from the April 2026 ASD GWAS, full weights for G2MH carriers and leave-SSC-out weights for SSC controls), (iv) adding all four PGS with their BMI interactions, and (v) all four PGS plus sex. All rows use the same individuals (G2MH carriers and SSC controls with all four scores). PGS x BMI p: whether the BMI slope depends on that PGS within carriers.

## Table S5. BMI effects with sex and polygenic scores as covariates

| Trait | Group | Model | n carriers | n controls | Carrier BMI slope (SE) | BMI slope p | PGS x BMI p | CNV x BMI vs controls (SE) | CNV x BMI p |
|---|---|---|---|---|---|---|---|---|---|
| Height | 16p11.2 del | unadjusted (same individuals) | 118 | 5010 | +0.345 (0.079) | 2.7e-05 |  | +0.213 (0.068) | 0.0019 |
| Height | 16p11.2 del | + sex | 118 | 5010 | +0.387 (0.124) | 0.0024 | 0.67 | +0.215 (0.068) | 0.0017 |
| Height | 16p11.2 del | + childhood obesity PGS | 118 | 5010 | +0.353 (0.081) | 3.0e-05 | 0.75 | +0.215 (0.068) | 0.0017 |
| Height | 16p11.2 del | + adult BMI PGS | 118 | 5010 | +0.343 (0.083) | 7.4e-05 | 0.048 | +0.216 (0.068) | 0.0016 |
| Height | 16p11.2 del | + height PGS | 118 | 5010 | +0.368 (0.079) | 8.3e-06 | 0.098 | +0.233 (0.059) | 7.3e-05 |
| Height | 16p11.2 del | + autism PGS (2026) | 118 | 5010 | +0.348 (0.078) | 1.9e-05 | 0.04 | +0.213 (0.069) | 0.0019 |
| Height | 16p11.2 del | + all four PGS | 118 | 5010 | +0.350 (0.082) | 4.9e-05 |  | +0.237 (0.059) | 5.1e-05 |
| Height | 16p11.2 del | + all four PGS + sex | 118 | 5010 | +0.364 (0.126) | 0.0048 |  | +0.243 (0.058) | 3.2e-05 |
| Height | 16p11.2 dup | unadjusted (same individuals) |  67 | 5010 | +0.072 (0.111) | 0.52 |  | -0.025 (0.099) | 0.8 |
| Height | 16p11.2 dup | + sex |  67 | 5010 | +0.164 (0.181) | 0.37 | 0.49 | -0.017 (0.099) | 0.86 |
| Height | 16p11.2 dup | + childhood obesity PGS |  67 | 5010 | -0.006 (0.117) | 0.96 | 0.052 | -0.026 (0.099) | 0.79 |
| Height | 16p11.2 dup | + adult BMI PGS |  67 | 5010 | +0.052 (0.116) | 0.65 | 0.064 | -0.025 (0.099) | 0.8 |
| Height | 16p11.2 dup | + height PGS |  67 | 5010 | +0.064 (0.097) | 0.51 | 0.14 | -0.053 (0.084) | 0.53 |
| Height | 16p11.2 dup | + autism PGS (2026) |  67 | 5010 | +0.108 (0.124) | 0.39 | 0.75 | -0.021 (0.100) | 0.83 |
| Height | 16p11.2 dup | + all four PGS |  67 | 5010 | +0.012 (0.113) | 0.91 |  | -0.041 (0.085) | 0.62 |
| Height | 16p11.2 dup | + all four PGS + sex |  67 | 5010 | +0.118 (0.158) | 0.46 |  | -0.031 (0.084) | 0.72 |
| Height | 22q11.2 del | unadjusted (same individuals) | 268 | 5010 | +0.063 (0.053) | 0.23 |  | -0.044 (0.050) | 0.38 |
| Height | 22q11.2 del | + sex | 268 | 5010 | +0.093 (0.081) | 0.25 | 0.57 | -0.046 (0.050) | 0.36 |
| Height | 22q11.2 del | + childhood obesity PGS | 268 | 5010 | +0.054 (0.053) | 0.32 | 0.81 | -0.044 (0.050) | 0.39 |
| Height | 22q11.2 del | + adult BMI PGS | 268 | 5010 | +0.086 (0.055) | 0.12 | 0.67 | -0.043 (0.050) | 0.39 |
| Height | 22q11.2 del | + height PGS | 268 | 5010 | +0.054 (0.048) | 0.26 | 0.24 | -0.055 (0.043) | 0.2 |
| Height | 22q11.2 del | + autism PGS (2026) | 268 | 5010 | +0.061 (0.053) | 0.25 | 0.52 | -0.044 (0.050) | 0.39 |
| Height | 22q11.2 del | + all four PGS | 268 | 5010 | +0.047 (0.051) | 0.35 |  | -0.048 (0.043) | 0.26 |
| Height | 22q11.2 del | + all four PGS + sex | 268 | 5010 | +0.065 (0.076) | 0.39 |  | -0.048 (0.043) | 0.26 |
| Height | 22q11.2 dup | unadjusted (same individuals) |  97 | 5010 | +0.116 (0.090) | 0.2 |  | +0.005 (0.090) | 0.96 |
| Height | 22q11.2 dup | + sex |  97 | 5010 | -0.183 (0.148) | 0.22 | 0.016 | +0.003 (0.090) | 0.97 |
| Height | 22q11.2 dup | + childhood obesity PGS |  97 | 5010 | +0.070 (0.090) | 0.43 | 0.34 | +0.004 (0.090) | 0.97 |
| Height | 22q11.2 dup | + adult BMI PGS |  97 | 5010 | +0.096 (0.091) | 0.29 | 0.19 | +0.006 (0.090) | 0.94 |
| Height | 22q11.2 dup | + height PGS |  97 | 5010 | +0.137 (0.078) | 0.082 | 0.88 | +0.011 (0.077) | 0.89 |
| Height | 22q11.2 dup | + autism PGS (2026) |  97 | 5010 | +0.060 (0.091) | 0.51 | 0.011 | +0.002 (0.090) | 0.98 |
| Height | 22q11.2 dup | + all four PGS |  97 | 5010 | +0.038 (0.076) | 0.62 |  | +0.018 (0.077) | 0.82 |
| Height | 22q11.2 dup | + all four PGS + sex |  97 | 5010 | -0.174 (0.127) | 0.18 |  | +0.024 (0.077) | 0.75 |
| Height | Overgrowth CNVs | unadjusted (same individuals) | 215 | 5010 | +0.263 (0.059) | 1.5e-05 |  | +0.139 (0.056) | 0.014 |
| Height | Overgrowth CNVs | + sex | 215 | 5010 | +0.225 (0.093) | 0.017 | 0.6 | +0.140 (0.056) | 0.013 |
| Height | Overgrowth CNVs | + childhood obesity PGS | 215 | 5010 | +0.256 (0.060) | 3.3e-05 | 0.39 | +0.140 (0.056) | 0.013 |
| Height | Overgrowth CNVs | + adult BMI PGS | 215 | 5010 | +0.245 (0.061) | 7.9e-05 | 0.019 | +0.141 (0.056) | 0.012 |
| Height | Overgrowth CNVs | + height PGS | 215 | 5010 | +0.290 (0.057) | 7.4e-07 | 0.42 | +0.159 (0.048) | 1.0e-03 |
| Height | Overgrowth CNVs | + autism PGS (2026) | 215 | 5010 | +0.268 (0.060) | 1.3e-05 | 0.99 | +0.138 (0.056) | 0.014 |
| Height | Overgrowth CNVs | + all four PGS | 215 | 5010 | +0.272 (0.059) | 8.3e-06 |  | +0.164 (0.048) | 6.5e-04 |
| Height | Overgrowth CNVs | + all four PGS + sex | 215 | 5010 | +0.248 (0.091) | 0.0072 |  | +0.170 (0.048) | 4.1e-04 |
| Height | Undergrowth CNVs | unadjusted (same individuals) | 335 | 5010 | +0.065 (0.048) | 0.18 |  | -0.040 (0.046) | 0.39 |
| Height | Undergrowth CNVs | + sex | 335 | 5010 | +0.092 (0.073) | 0.21 | 0.6 | -0.040 (0.046) | 0.39 |
| Height | Undergrowth CNVs | + childhood obesity PGS | 335 | 5010 | +0.055 (0.048) | 0.26 | 0.37 | -0.040 (0.046) | 0.39 |
| Height | Undergrowth CNVs | + adult BMI PGS | 335 | 5010 | +0.084 (0.050) | 0.092 | 0.33 | -0.039 (0.046) | 0.4 |
| Height | Undergrowth CNVs | + height PGS | 335 | 5010 | +0.053 (0.043) | 0.21 | 0.11 | -0.053 (0.039) | 0.18 |
| Height | Undergrowth CNVs | + autism PGS (2026) | 335 | 5010 | +0.067 (0.048) | 0.16 | 0.55 | -0.039 (0.046) | 0.4 |
| Height | Undergrowth CNVs | + all four PGS | 335 | 5010 | +0.046 (0.045) | 0.32 |  | -0.049 (0.039) | 0.22 |
| Height | Undergrowth CNVs | + all four PGS + sex | 335 | 5010 | +0.064 (0.068) | 0.34 |  | -0.045 (0.039) | 0.26 |
| Height | All CNVs | unadjusted (same individuals) | 550 | 5010 | +0.123 (0.039) | 0.0016 |  | +0.009 (0.039) | 0.81 |
| Height | All CNVs | + sex | 550 | 5010 | +0.085 (0.059) | 0.15 | 0.41 | +0.010 (0.039) | 0.8 |
| Height | All CNVs | + childhood obesity PGS | 550 | 5010 | +0.121 (0.039) | 0.0023 | 0.49 | +0.010 (0.039) | 0.8 |
| Height | All CNVs | + adult BMI PGS | 550 | 5010 | +0.124 (0.040) | 0.0024 | 0.017 | +0.011 (0.038) | 0.78 |
| Height | All CNVs | + height PGS | 550 | 5010 | +0.131 (0.036) | 3.2e-04 | 0.31 | +0.007 (0.033) | 0.83 |
| Height | All CNVs | + autism PGS (2026) | 550 | 5010 | +0.124 (0.039) | 0.0016 | 0.58 | +0.009 (0.039) | 0.81 |
| Height | All CNVs | + all four PGS | 550 | 5010 | +0.119 (0.038) | 0.0017 |  | +0.013 (0.033) | 0.7 |
| Height | All CNVs | + all four PGS + sex | 550 | 5010 | +0.076 (0.056) | 0.17 |  | +0.014 (0.033) | 0.66 |
| Head circumference | 16p11.2 del | unadjusted (same individuals) |  98 | 4985 | +0.486 (0.080) | 2.5e-08 |  | +0.130 (0.076) | 0.086 |
| Head circumference | 16p11.2 del | + sex |  98 | 4985 | +0.598 (0.124) | 5.9e-06 | 0.24 | +0.139 (0.075) | 0.065 |
| Head circumference | 16p11.2 del | + childhood obesity PGS |  98 | 4985 | +0.496 (0.082) | 3.0e-08 | 0.19 | +0.133 (0.076) | 0.079 |
| Head circumference | 16p11.2 del | + adult BMI PGS |  98 | 4985 | +0.501 (0.085) | 6.0e-08 | 0.075 | +0.132 (0.076) | 0.082 |
| Head circumference | 16p11.2 del | + height PGS |  98 | 4985 | +0.497 (0.081) | 1.8e-08 | 0.32 | +0.136 (0.075) | 0.068 |
| Head circumference | 16p11.2 del | + autism PGS (2026) |  98 | 4985 | +0.486 (0.080) | 3.0e-08 | 0.45 | +0.131 (0.076) | 0.083 |
| Head circumference | 16p11.2 del | + all four PGS |  98 | 4985 | +0.515 (0.087) | 6.8e-08 |  | +0.142 (0.075) | 0.057 |
| Head circumference | 16p11.2 del | + all four PGS + sex |  98 | 4985 | +0.550 (0.133) | 8.4e-05 |  | +0.151 (0.074) | 0.041 |
| Head circumference | 16p11.2 dup | unadjusted (same individuals) |  54 | 4985 | +0.109 (0.115) | 0.35 |  | -0.251 (0.116) | 0.03 |
| Head circumference | 16p11.2 dup | + sex |  54 | 4985 | +0.128 (0.203) | 0.53 | 0.92 | -0.216 (0.115) | 0.061 |
| Head circumference | 16p11.2 dup | + childhood obesity PGS |  54 | 4985 | +0.112 (0.123) | 0.37 | 0.85 | -0.251 (0.116) | 0.031 |
| Head circumference | 16p11.2 dup | + adult BMI PGS |  54 | 4985 | +0.093 (0.124) | 0.46 | 0.93 | -0.250 (0.116) | 0.031 |
| Head circumference | 16p11.2 dup | + height PGS |  54 | 4985 | +0.070 (0.121) | 0.57 | 0.28 | -0.256 (0.114) | 0.026 |
| Head circumference | 16p11.2 dup | + autism PGS (2026) |  54 | 4985 | +0.083 (0.131) | 0.53 | 0.5 | -0.242 (0.116) | 0.037 |
| Head circumference | 16p11.2 dup | + all four PGS |  54 | 4985 | +0.076 (0.151) | 0.62 |  | -0.246 (0.115) | 0.032 |
| Head circumference | 16p11.2 dup | + all four PGS + sex |  54 | 4985 | +0.075 (0.241) | 0.76 |  | -0.212 (0.114) | 0.064 |
| Head circumference | 22q11.2 del | unadjusted (same individuals) | 171 | 4985 | +0.322 (0.067) | 3.7e-06 |  | -0.042 (0.064) | 0.52 |
| Head circumference | 22q11.2 del | + sex | 171 | 4985 | +0.375 (0.098) | 1.8e-04 | 0.46 | -0.042 (0.064) | 0.51 |
| Head circumference | 22q11.2 del | + childhood obesity PGS | 171 | 4985 | +0.307 (0.068) | 1.0e-05 | 0.024 | -0.042 (0.064) | 0.52 |
| Head circumference | 22q11.2 del | + adult BMI PGS | 171 | 4985 | +0.308 (0.071) | 2.8e-05 | 0.35 | -0.043 (0.064) | 0.51 |
| Head circumference | 22q11.2 del | + height PGS | 171 | 4985 | +0.357 (0.069) | 6.4e-07 | 0.028 | -0.051 (0.064) | 0.42 |
| Head circumference | 22q11.2 del | + autism PGS (2026) | 171 | 4985 | +0.324 (0.068) | 3.7e-06 | 0.79 | -0.042 (0.064) | 0.51 |
| Head circumference | 22q11.2 del | + all four PGS | 171 | 4985 | +0.321 (0.073) | 2.2e-05 |  | -0.053 (0.064) | 0.41 |
| Head circumference | 22q11.2 del | + all four PGS + sex | 171 | 4985 | +0.346 (0.101) | 7.9e-04 |  | -0.053 (0.063) | 0.41 |
| Head circumference | 22q11.2 dup | unadjusted (same individuals) |  72 | 4985 | +0.518 (0.114) | 2.5e-05 |  | +0.137 (0.108) | 0.2 |
| Head circumference | 22q11.2 dup | + sex |  72 | 4985 | +0.537 (0.199) | 0.0088 | 0.83 | +0.153 (0.107) | 0.15 |
| Head circumference | 22q11.2 dup | + childhood obesity PGS |  72 | 4985 | +0.498 (0.121) | 1.1e-04 | 0.87 | +0.133 (0.108) | 0.22 |
| Head circumference | 22q11.2 dup | + adult BMI PGS |  72 | 4985 | +0.496 (0.115) | 5.4e-05 | 0.15 | +0.140 (0.108) | 0.19 |
| Head circumference | 22q11.2 dup | + height PGS |  72 | 4985 | +0.564 (0.107) | 1.6e-06 | 0.86 | +0.153 (0.106) | 0.15 |
| Head circumference | 22q11.2 dup | + autism PGS (2026) |  72 | 4985 | +0.518 (0.119) | 4.7e-05 | 0.9 | +0.138 (0.108) | 0.2 |
| Head circumference | 22q11.2 dup | + all four PGS |  72 | 4985 | +0.518 (0.118) | 4.7e-05 |  | +0.159 (0.106) | 0.13 |
| Head circumference | 22q11.2 dup | + all four PGS + sex |  72 | 4985 | +0.660 (0.204) | 0.002 |  | +0.175 (0.105) | 0.097 |
| Head circumference | Overgrowth CNVs | unadjusted (same individuals) | 170 | 4985 | +0.492 (0.066) | 3.4e-12 |  | +0.130 (0.064) | 0.042 |
| Head circumference | Overgrowth CNVs | + sex | 170 | 4985 | +0.611 (0.104) | 2.4e-08 | 0.14 | +0.141 (0.063) | 0.027 |
| Head circumference | Overgrowth CNVs | + childhood obesity PGS | 170 | 4985 | +0.485 (0.066) | 9.3e-12 | 0.33 | +0.132 (0.064) | 0.038 |
| Head circumference | Overgrowth CNVs | + adult BMI PGS | 170 | 4985 | +0.492 (0.067) | 1.2e-11 | 0.027 | +0.133 (0.064) | 0.037 |
| Head circumference | Overgrowth CNVs | + height PGS | 170 | 4985 | +0.514 (0.065) | 2.6e-13 | 0.64 | +0.141 (0.063) | 0.025 |
| Head circumference | Overgrowth CNVs | + autism PGS (2026) | 170 | 4985 | +0.497 (0.066) | 3.7e-12 | 0.72 | +0.131 (0.064) | 0.039 |
| Head circumference | Overgrowth CNVs | + all four PGS | 170 | 4985 | +0.512 (0.067) | 1.8e-12 |  | +0.147 (0.063) | 0.019 |
| Head circumference | Overgrowth CNVs | + all four PGS + sex | 170 | 4985 | +0.614 (0.105) | 2.6e-08 |  | +0.158 (0.062) | 0.011 |
| Head circumference | Undergrowth CNVs | unadjusted (same individuals) | 225 | 4985 | +0.269 (0.059) | 9.7e-06 |  | -0.097 (0.058) | 0.094 |
| Head circumference | Undergrowth CNVs | + sex | 225 | 4985 | +0.323 (0.089) | 3.6e-04 | 0.41 | -0.090 (0.057) | 0.12 |
| Head circumference | Undergrowth CNVs | + childhood obesity PGS | 225 | 4985 | +0.266 (0.059) | 1.3e-05 | 0.038 | -0.096 (0.058) | 0.097 |
| Head circumference | Undergrowth CNVs | + adult BMI PGS | 225 | 4985 | +0.253 (0.063) | 7.8e-05 | 0.29 | -0.097 (0.058) | 0.093 |
| Head circumference | Undergrowth CNVs | + height PGS | 225 | 4985 | +0.287 (0.060) | 3.1e-06 | 0.048 | -0.105 (0.057) | 0.066 |
| Head circumference | Undergrowth CNVs | + autism PGS (2026) | 225 | 4985 | +0.273 (0.060) | 8.3e-06 | 0.8 | -0.095 (0.058) | 0.1 |
| Head circumference | Undergrowth CNVs | + all four PGS | 225 | 4985 | +0.274 (0.064) | 2.6e-05 |  | -0.103 (0.057) | 0.07 |
| Head circumference | Undergrowth CNVs | + all four PGS + sex | 225 | 4985 | +0.295 (0.092) | 0.0015 |  | -0.096 (0.057) | 0.09 |
| Head circumference | All CNVs | unadjusted (same individuals) | 395 | 4985 | +0.367 (0.044) | 1.6e-15 |  | +0.001 (0.046) | 0.98 |
| Head circumference | All CNVs | + sex | 395 | 4985 | +0.422 (0.067) | 1.0e-09 | 0.28 | +0.010 (0.045) | 0.83 |
| Head circumference | All CNVs | + childhood obesity PGS | 395 | 4985 | +0.358 (0.044) | 1.1e-14 | 0.3 | +0.002 (0.046) | 0.96 |
| Head circumference | All CNVs | + adult BMI PGS | 395 | 4985 | +0.357 (0.046) | 8.1e-14 | 0.031 | +0.002 (0.046) | 0.97 |
| Head circumference | All CNVs | + height PGS | 395 | 4985 | +0.371 (0.044) | 7.2e-16 | 0.28 | +0.000 (0.045) | 0.99 |
| Head circumference | All CNVs | + autism PGS (2026) | 395 | 4985 | +0.367 (0.044) | 1.9e-15 | 0.96 | +0.003 (0.046) | 0.96 |
| Head circumference | All CNVs | + all four PGS | 395 | 4985 | +0.358 (0.046) | 7.2e-14 |  | +0.004 (0.045) | 0.94 |
| Head circumference | All CNVs | + all four PGS + sex | 395 | 4985 | +0.397 (0.068) | 1.3e-08 |  | +0.012 (0.045) | 0.78 |
| SRS | 16p11.2 del | unadjusted (same individuals) |  74 | 4988 | +0.169 (0.211) | 0.43 |  | +0.161 (0.086) | 0.063 |
| SRS | 16p11.2 del | + sex |  74 | 4988 | +0.213 (0.293) | 0.47 | 0.86 | +0.161 (0.086) | 0.062 |
| SRS | 16p11.2 del | + childhood obesity PGS |  74 | 4988 | +0.157 (0.206) | 0.45 | 0.35 | +0.163 (0.086) | 0.059 |
| SRS | 16p11.2 del | + adult BMI PGS |  74 | 4988 | +0.250 (0.222) | 0.26 | 0.18 | +0.167 (0.086) | 0.053 |
| SRS | 16p11.2 del | + height PGS |  74 | 4988 | +0.128 (0.216) | 0.56 | 0.26 | +0.162 (0.086) | 0.061 |
| SRS | 16p11.2 del | + autism PGS (2026) |  74 | 4988 | +0.185 (0.209) | 0.38 | 0.058 | +0.164 (0.086) | 0.057 |
| SRS | 16p11.2 del | + all four PGS |  74 | 4988 | +0.203 (0.223) | 0.37 |  | +0.179 (0.086) | 0.038 |
| SRS | 16p11.2 del | + all four PGS + sex |  74 | 4988 | +0.157 (0.329) | 0.63 |  | +0.179 (0.086) | 0.038 |
| SRS | 16p11.2 dup | unadjusted (same individuals) |  29 | 4988 | +0.208 (0.346) | 0.55 |  | +0.159 (0.161) | 0.32 |
| SRS | 16p11.2 dup | + sex |  29 | 4988 | +0.366 (0.517) | 0.49 | 0.95 | +0.169 (0.161) | 0.29 |
| SRS | 16p11.2 dup | + childhood obesity PGS |  29 | 4988 | -0.072 (0.395) | 0.86 | 0.14 | +0.151 (0.161) | 0.35 |
| SRS | 16p11.2 dup | + adult BMI PGS |  29 | 4988 | +0.378 (0.363) | 0.31 | 0.15 | +0.158 (0.161) | 0.33 |
| SRS | 16p11.2 dup | + height PGS |  29 | 4988 | +0.252 (0.351) | 0.48 | 0.39 | +0.160 (0.161) | 0.32 |
| SRS | 16p11.2 dup | + autism PGS (2026) |  29 | 4988 | +0.277 (0.377) | 0.47 | 0.18 | +0.184 (0.161) | 0.25 |
| SRS | 16p11.2 dup | + all four PGS |  29 | 4988 | +0.745 (0.501) | 0.15 |  | +0.171 (0.161) | 0.29 |
| SRS | 16p11.2 dup | + all four PGS + sex |  29 | 4988 | +1.092 (0.924) | 0.25 |  | +0.182 (0.161) | 0.26 |
| SRS | 22q11.2 del | unadjusted (same individuals) | 186 | 4988 | +0.313 (0.104) | 0.0031 |  | +0.293 (0.055) | 1.3e-07 |
| SRS | 22q11.2 del | + sex | 186 | 4988 | +0.430 (0.169) | 0.012 | 0.3 | +0.292 (0.055) | 1.4e-07 |
| SRS | 22q11.2 del | + childhood obesity PGS | 186 | 4988 | +0.305 (0.107) | 0.0047 | 0.78 | +0.292 (0.055) | 1.4e-07 |
| SRS | 22q11.2 del | + adult BMI PGS | 186 | 4988 | +0.338 (0.108) | 0.002 | 0.6 | +0.293 (0.055) | 1.4e-07 |
| SRS | 22q11.2 del | + height PGS | 186 | 4988 | +0.318 (0.105) | 0.0028 | 0.81 | +0.293 (0.055) | 1.3e-07 |
| SRS | 22q11.2 del | + autism PGS (2026) | 186 | 4988 | +0.294 (0.105) | 0.0055 | 0.074 | +0.293 (0.055) | 1.2e-07 |
| SRS | 22q11.2 del | + all four PGS | 186 | 4988 | +0.304 (0.110) | 0.0063 |  | +0.292 (0.055) | 1.4e-07 |
| SRS | 22q11.2 del | + all four PGS + sex | 186 | 4988 | +0.374 (0.177) | 0.036 |  | +0.292 (0.055) | 1.5e-07 |
| SRS | 22q11.2 dup | unadjusted (same individuals) |  60 | 4988 | -0.122 (0.272) | 0.66 |  | -0.254 (0.125) | 0.043 |
| SRS | 22q11.2 dup | + sex |  60 | 4988 | +0.019 (0.431) | 0.97 | 0.2 | -0.264 (0.125) | 0.035 |
| SRS | 22q11.2 dup | + childhood obesity PGS |  60 | 4988 | -0.151 (0.280) | 0.59 | 0.18 | -0.255 (0.125) | 0.042 |
| SRS | 22q11.2 dup | + adult BMI PGS |  60 | 4988 | -0.129 (0.289) | 0.66 | 0.83 | -0.249 (0.125) | 0.047 |
| SRS | 22q11.2 dup | + height PGS |  60 | 4988 | -0.150 (0.274) | 0.59 | 0.79 | -0.254 (0.125) | 0.043 |
| SRS | 22q11.2 dup | + autism PGS (2026) |  60 | 4988 | -0.167 (0.275) | 0.55 | 0.65 | -0.259 (0.125) | 0.039 |
| SRS | 22q11.2 dup | + all four PGS |  60 | 4988 | -0.211 (0.307) | 0.49 |  | -0.252 (0.125) | 0.045 |
| SRS | 22q11.2 dup | + all four PGS + sex |  60 | 4988 | +0.113 (0.448) | 0.8 |  | -0.263 (0.125) | 0.036 |
| SRS | Overgrowth CNVs | unadjusted (same individuals) | 134 | 4988 | +0.071 (0.164) | 0.66 |  | +0.018 (0.074) | 0.81 |
| SRS | Overgrowth CNVs | + sex | 134 | 4988 | +0.389 (0.233) | 0.097 | 0.063 | +0.017 (0.074) | 0.81 |
| SRS | Overgrowth CNVs | + childhood obesity PGS | 134 | 4988 | +0.039 (0.161) | 0.81 | 0.014 | +0.019 (0.074) | 0.8 |
| SRS | Overgrowth CNVs | + adult BMI PGS | 134 | 4988 | +0.131 (0.172) | 0.45 | 0.13 | +0.024 (0.074) | 0.74 |
| SRS | Overgrowth CNVs | + height PGS | 134 | 4988 | +0.059 (0.163) | 0.72 | 0.25 | +0.018 (0.074) | 0.8 |
| SRS | Overgrowth CNVs | + autism PGS (2026) | 134 | 4988 | +0.073 (0.162) | 0.65 | 0.026 | +0.019 (0.074) | 0.79 |
| SRS | Overgrowth CNVs | + all four PGS | 134 | 4988 | +0.108 (0.167) | 0.52 |  | +0.032 (0.074) | 0.67 |
| SRS | Overgrowth CNVs | + all four PGS + sex | 134 | 4988 | +0.334 (0.242) | 0.17 |  | +0.031 (0.074) | 0.67 |
| SRS | Undergrowth CNVs | unadjusted (same individuals) | 215 | 4988 | +0.301 (0.101) | 0.0033 |  | +0.281 (0.054) | 2.0e-07 |
| SRS | Undergrowth CNVs | + sex | 215 | 4988 | +0.446 (0.162) | 0.0063 | 0.22 | +0.282 (0.054) | 1.8e-07 |
| SRS | Undergrowth CNVs | + childhood obesity PGS | 215 | 4988 | +0.292 (0.104) | 0.0053 | 0.69 | +0.279 (0.054) | 2.4e-07 |
| SRS | Undergrowth CNVs | + adult BMI PGS | 215 | 4988 | +0.331 (0.105) | 0.0018 | 0.93 | +0.280 (0.054) | 2.1e-07 |
| SRS | Undergrowth CNVs | + height PGS | 215 | 4988 | +0.303 (0.102) | 0.0032 | 0.95 | +0.281 (0.054) | 1.9e-07 |
| SRS | Undergrowth CNVs | + autism PGS (2026) | 215 | 4988 | +0.303 (0.101) | 0.0032 | 0.18 | +0.284 (0.054) | 1.4e-07 |
| SRS | Undergrowth CNVs | + all four PGS | 215 | 4988 | +0.320 (0.107) | 0.0032 |  | +0.282 (0.054) | 1.9e-07 |
| SRS | Undergrowth CNVs | + all four PGS + sex | 215 | 4988 | +0.418 (0.170) | 0.015 |  | +0.283 (0.054) | 1.7e-07 |
| SRS | All CNVs | unadjusted (same individuals) | 349 | 4988 | +0.219 (0.088) | 0.013 |  | +0.192 (0.048) | 5.8e-05 |
| SRS | All CNVs | + sex | 349 | 4988 | +0.456 (0.133) | 6.8e-04 | 0.013 | +0.194 (0.048) | 5.0e-05 |
| SRS | All CNVs | + childhood obesity PGS | 349 | 4988 | +0.190 (0.089) | 0.033 | 0.052 | +0.190 (0.048) | 6.9e-05 |
| SRS | All CNVs | + adult BMI PGS | 349 | 4988 | +0.226 (0.091) | 0.013 | 0.53 | +0.194 (0.048) | 5.2e-05 |
| SRS | All CNVs | + height PGS | 349 | 4988 | +0.219 (0.088) | 0.013 | 0.44 | +0.193 (0.048) | 5.8e-05 |
| SRS | All CNVs | + autism PGS (2026) | 349 | 4988 | +0.224 (0.088) | 0.011 | 0.019 | +0.196 (0.048) | 4.1e-05 |
| SRS | All CNVs | + all four PGS | 349 | 4988 | +0.205 (0.090) | 0.024 |  | +0.196 (0.048) | 4.1e-05 |
| SRS | All CNVs | + all four PGS + sex | 349 | 4988 | +0.396 (0.135) | 0.0035 |  | +0.198 (0.048) | 3.6e-05 |
| SCQ | 16p11.2 del | unadjusted (same individuals) |  56 | 1030 | +0.411 (0.318) | 0.2 |  | +0.362 (0.116) | 0.0018 |
| SCQ | 16p11.2 del | + sex |  56 | 1030 | +0.646 (0.457) | 0.16 | 0.47 | +0.365 (0.116) | 0.0016 |
| SCQ | 16p11.2 del | + childhood obesity PGS |  56 | 1030 | +0.283 (0.320) | 0.38 | 0.34 | +0.360 (0.116) | 0.0019 |
| SCQ | 16p11.2 del | + adult BMI PGS |  56 | 1030 | +0.484 (0.330) | 0.15 | 0.94 | +0.362 (0.116) | 0.0019 |
| SCQ | 16p11.2 del | + height PGS |  56 | 1030 | +0.352 (0.320) | 0.28 | 0.29 | +0.357 (0.116) | 0.0021 |
| SCQ | 16p11.2 del | + autism PGS (2026) |  56 | 1030 | +0.397 (0.323) | 0.22 | 0.46 | +0.364 (0.116) | 0.0017 |
| SCQ | 16p11.2 del | + all four PGS |  56 | 1030 | +0.228 (0.347) | 0.52 |  | +0.352 (0.116) | 0.0024 |
| SCQ | 16p11.2 del | + all four PGS + sex |  56 | 1030 | +0.368 (0.501) | 0.47 |  | +0.355 (0.116) | 0.0022 |
| SCQ | 22q11.2 del | unadjusted (same individuals) | 151 | 1030 | +0.369 (0.166) | 0.028 |  | +0.350 (0.091) | 1.2e-04 |
| SCQ | 22q11.2 del | + sex | 151 | 1030 | +0.335 (0.255) | 0.19 | 0.95 | +0.339 (0.091) | 1.9e-04 |
| SCQ | 22q11.2 del | + childhood obesity PGS | 151 | 1030 | +0.404 (0.168) | 0.017 | 0.99 | +0.344 (0.091) | 1.6e-04 |
| SCQ | 22q11.2 del | + adult BMI PGS | 151 | 1030 | +0.387 (0.171) | 0.026 | 0.36 | +0.350 (0.091) | 1.2e-04 |
| SCQ | 22q11.2 del | + height PGS | 151 | 1030 | +0.368 (0.167) | 0.03 | 0.87 | +0.350 (0.091) | 1.2e-04 |
| SCQ | 22q11.2 del | + autism PGS (2026) | 151 | 1030 | +0.354 (0.173) | 0.042 | 0.68 | +0.351 (0.091) | 1.2e-04 |
| SCQ | 22q11.2 del | + all four PGS | 151 | 1030 | +0.401 (0.181) | 0.028 |  | +0.347 (0.091) | 1.4e-04 |
| SCQ | 22q11.2 del | + all four PGS + sex | 151 | 1030 | +0.410 (0.274) | 0.14 |  | +0.337 (0.091) | 2.2e-04 |
| SCQ | 22q11.2 dup | unadjusted (same individuals) |  44 | 1030 | +0.218 (0.452) | 0.63 |  | +0.114 (0.153) | 0.46 |
| SCQ | 22q11.2 dup | + sex |  44 | 1030 | +0.273 (0.652) | 0.68 | 0.57 | +0.101 (0.153) | 0.51 |
| SCQ | 22q11.2 dup | + childhood obesity PGS |  44 | 1030 | +0.286 (0.460) | 0.54 | 0.22 | +0.114 (0.153) | 0.46 |
| SCQ | 22q11.2 dup | + adult BMI PGS |  44 | 1030 | +0.425 (0.514) | 0.41 | 0.47 | +0.126 (0.154) | 0.41 |
| SCQ | 22q11.2 dup | + height PGS |  44 | 1030 | +0.150 (0.454) | 0.74 | 0.51 | +0.108 (0.153) | 0.48 |
| SCQ | 22q11.2 dup | + autism PGS (2026) |  44 | 1030 | +0.036 (0.445) | 0.94 | 0.12 | +0.095 (0.153) | 0.54 |
| SCQ | 22q11.2 dup | + all four PGS |  44 | 1030 | +0.046 (0.580) | 0.94 |  | +0.105 (0.155) | 0.5 |
| SCQ | 22q11.2 dup | + all four PGS + sex |  44 | 1030 | +0.604 (0.778) | 0.44 |  | +0.093 (0.155) | 0.55 |
| SCQ | Overgrowth CNVs | unadjusted (same individuals) | 100 | 1030 | +0.349 (0.258) | 0.18 |  | +0.316 (0.109) | 0.0038 |
| SCQ | Overgrowth CNVs | + sex | 100 | 1030 | +0.580 (0.357) | 0.11 | 0.27 | +0.314 (0.109) | 0.0039 |
| SCQ | Overgrowth CNVs | + childhood obesity PGS | 100 | 1030 | +0.306 (0.265) | 0.25 | 0.46 | +0.316 (0.109) | 0.0037 |
| SCQ | Overgrowth CNVs | + adult BMI PGS | 100 | 1030 | +0.348 (0.263) | 0.19 | 0.91 | +0.318 (0.109) | 0.0036 |
| SCQ | Overgrowth CNVs | + height PGS | 100 | 1030 | +0.354 (0.262) | 0.18 | 0.84 | +0.312 (0.109) | 0.0043 |
| SCQ | Overgrowth CNVs | + autism PGS (2026) | 100 | 1030 | +0.316 (0.260) | 0.23 | 0.88 | +0.308 (0.109) | 0.0048 |
| SCQ | Overgrowth CNVs | + all four PGS | 100 | 1030 | +0.277 (0.275) | 0.32 |  | +0.307 (0.109) | 0.0051 |
| SCQ | Overgrowth CNVs | + all four PGS + sex | 100 | 1030 | +0.576 (0.389) | 0.14 |  | +0.306 (0.109) | 0.0051 |
| SCQ | Undergrowth CNVs | unadjusted (same individuals) | 168 | 1030 | +0.330 (0.165) | 0.047 |  | +0.305 (0.090) | 7.9e-04 |
| SCQ | Undergrowth CNVs | + sex | 168 | 1030 | +0.386 (0.252) | 0.13 | 0.59 | +0.294 (0.090) | 0.0012 |
| SCQ | Undergrowth CNVs | + childhood obesity PGS | 168 | 1030 | +0.354 (0.166) | 0.035 | 0.82 | +0.298 (0.091) | 0.001 |
| SCQ | Undergrowth CNVs | + adult BMI PGS | 168 | 1030 | +0.366 (0.169) | 0.032 | 0.47 | +0.302 (0.091) | 9.1e-04 |
| SCQ | Undergrowth CNVs | + height PGS | 168 | 1030 | +0.330 (0.165) | 0.048 | 0.86 | +0.306 (0.091) | 7.6e-04 |
| SCQ | Undergrowth CNVs | + autism PGS (2026) | 168 | 1030 | +0.330 (0.168) | 0.052 | 0.57 | +0.309 (0.091) | 6.6e-04 |
| SCQ | Undergrowth CNVs | + all four PGS | 168 | 1030 | +0.394 (0.176) | 0.027 |  | +0.303 (0.091) | 9.0e-04 |
| SCQ | Undergrowth CNVs | + all four PGS + sex | 168 | 1030 | +0.466 (0.268) | 0.084 |  | +0.293 (0.091) | 0.0013 |
| SCQ | All CNVs | unadjusted (same individuals) | 268 | 1030 | +0.311 (0.143) | 0.031 |  | +0.318 (0.091) | 5.0e-04 |
| SCQ | All CNVs | + sex | 268 | 1030 | +0.468 (0.207) | 0.025 | 0.17 | +0.311 (0.091) | 6.6e-04 |
| SCQ | All CNVs | + childhood obesity PGS | 268 | 1030 | +0.312 (0.145) | 0.033 | 0.64 | +0.314 (0.092) | 6.2e-04 |
| SCQ | All CNVs | + adult BMI PGS | 268 | 1030 | +0.322 (0.146) | 0.028 | 0.6 | +0.319 (0.092) | 5.2e-04 |
| SCQ | All CNVs | + height PGS | 268 | 1030 | +0.311 (0.144) | 0.031 | 0.99 | +0.317 (0.091) | 5.3e-04 |
| SCQ | All CNVs | + autism PGS (2026) | 268 | 1030 | +0.316 (0.144) | 0.029 | 0.81 | +0.322 (0.091) | 4.3e-04 |
| SCQ | All CNVs | + all four PGS | 268 | 1030 | +0.329 (0.149) | 0.028 |  | +0.318 (0.092) | 5.6e-04 |
| SCQ | All CNVs | + all four PGS + sex | 268 | 1030 | +0.473 (0.213) | 0.027 |  | +0.312 (0.092) | 7.1e-04 |
| CNB accuracy | 16p11.2 del | unadjusted (same individuals) | 101 | 0 | +0.014 (0.045) | 0.75 |  |  |  |
| CNB accuracy | 16p11.2 del | + sex | 101 | 0 | +0.067 (0.064) | 0.29 | 0.23 |  |  |
| CNB accuracy | 16p11.2 del | + childhood obesity PGS | 101 | 0 | +0.014 (0.045) | 0.76 | 0.43 |  |  |
| CNB accuracy | 16p11.2 del | + adult BMI PGS | 101 | 0 | +0.007 (0.047) | 0.88 | 0.43 |  |  |
| CNB accuracy | 16p11.2 del | + height PGS | 101 | 0 | +0.013 (0.046) | 0.78 | 0.17 |  |  |
| CNB accuracy | 16p11.2 del | + autism PGS (2026) | 101 | 0 | +0.009 (0.045) | 0.85 | 0.44 |  |  |
| CNB accuracy | 16p11.2 del | + all four PGS | 101 | 0 | -0.007 (0.049) | 0.88 |  |  |  |
| CNB accuracy | 16p11.2 del | + all four PGS + sex | 101 | 0 | +0.079 (0.070) | 0.27 |  |  |  |
| CNB accuracy | 16p11.2 dup | unadjusted (same individuals) |  39 | 0 | -0.168 (0.078) | 0.038 |  |  |  |
| CNB accuracy | 16p11.2 dup | + sex |  39 | 0 | -0.327 (0.110) | 0.0053 | 0.037 |  |  |
| CNB accuracy | 16p11.2 dup | + childhood obesity PGS |  39 | 0 | -0.193 (0.086) | 0.031 | 0.2 |  |  |
| CNB accuracy | 16p11.2 dup | + adult BMI PGS |  39 | 0 | -0.155 (0.083) | 0.07 | 0.57 |  |  |
| CNB accuracy | 16p11.2 dup | + height PGS |  39 | 0 | -0.153 (0.079) | 0.063 | 0.42 |  |  |
| CNB accuracy | 16p11.2 dup | + autism PGS (2026) |  39 | 0 | -0.116 (0.084) | 0.18 | 0.24 |  |  |
| CNB accuracy | 16p11.2 dup | + all four PGS |  39 | 0 | -0.080 (0.099) | 0.42 |  |  |  |
| CNB accuracy | 16p11.2 dup | + all four PGS + sex |  39 | 0 | -0.297 (0.123) | 0.022 |  |  |  |
| CNB accuracy | 22q11.2 del | unadjusted (same individuals) | 282 | 0 | -0.032 (0.031) | 0.3 |  |  |  |
| CNB accuracy | 22q11.2 del | + sex | 282 | 0 | -0.021 (0.049) | 0.67 | 0.8 |  |  |
| CNB accuracy | 22q11.2 del | + childhood obesity PGS | 282 | 0 | -0.008 (0.031) | 0.79 | 0.032 |  |  |
| CNB accuracy | 22q11.2 del | + adult BMI PGS | 282 | 0 | -0.031 (0.032) | 0.34 | 0.54 |  |  |
| CNB accuracy | 22q11.2 del | + height PGS | 282 | 0 | -0.036 (0.031) | 0.24 | 0.22 |  |  |
| CNB accuracy | 22q11.2 del | + autism PGS (2026) | 282 | 0 | -0.032 (0.031) | 0.31 | 0.92 |  |  |
| CNB accuracy | 22q11.2 del | + all four PGS | 282 | 0 | -0.016 (0.032) | 0.63 |  |  |  |
| CNB accuracy | 22q11.2 del | + all four PGS + sex | 282 | 0 | -0.005 (0.050) | 0.93 |  |  |  |
| CNB accuracy | 22q11.2 dup | unadjusted (same individuals) | 110 | 0 | +0.077 (0.057) | 0.18 |  |  |  |
| CNB accuracy | 22q11.2 dup | + sex | 110 | 0 | +0.272 (0.110) | 0.015 | 0.037 |  |  |
| CNB accuracy | 22q11.2 dup | + childhood obesity PGS | 110 | 0 | +0.096 (0.058) | 0.1 | 0.74 |  |  |
| CNB accuracy | 22q11.2 dup | + adult BMI PGS | 110 | 0 | +0.081 (0.058) | 0.17 | 0.77 |  |  |
| CNB accuracy | 22q11.2 dup | + height PGS | 110 | 0 | +0.083 (0.058) | 0.16 | 0.81 |  |  |
| CNB accuracy | 22q11.2 dup | + autism PGS (2026) | 110 | 0 | +0.062 (0.060) | 0.3 | 0.71 |  |  |
| CNB accuracy | 22q11.2 dup | + all four PGS | 110 | 0 | +0.095 (0.063) | 0.13 |  |  |  |
| CNB accuracy | 22q11.2 dup | + all four PGS + sex | 110 | 0 | +0.276 (0.121) | 0.025 |  |  |  |
| CNB accuracy | Overgrowth CNVs | unadjusted (same individuals) | 211 | 0 | +0.048 (0.037) | 0.19 |  |  |  |
| CNB accuracy | Overgrowth CNVs | + sex | 211 | 0 | +0.141 (0.057) | 0.015 | 0.032 |  |  |
| CNB accuracy | Overgrowth CNVs | + childhood obesity PGS | 211 | 0 | +0.053 (0.037) | 0.15 | 0.66 |  |  |
| CNB accuracy | Overgrowth CNVs | + adult BMI PGS | 211 | 0 | +0.049 (0.037) | 0.19 | 0.67 |  |  |
| CNB accuracy | Overgrowth CNVs | + height PGS | 211 | 0 | +0.046 (0.037) | 0.22 | 0.28 |  |  |
| CNB accuracy | Overgrowth CNVs | + autism PGS (2026) | 211 | 0 | +0.028 (0.037) | 0.45 | 0.38 |  |  |
| CNB accuracy | Overgrowth CNVs | + all four PGS | 211 | 0 | +0.027 (0.038) | 0.49 |  |  |  |
| CNB accuracy | Overgrowth CNVs | + all four PGS + sex | 211 | 0 | +0.119 (0.060) | 0.05 |  |  |  |
| CNB accuracy | Undergrowth CNVs | unadjusted (same individuals) | 321 | 0 | -0.040 (0.029) | 0.16 |  |  |  |
| CNB accuracy | Undergrowth CNVs | + sex | 321 | 0 | -0.044 (0.045) | 0.33 | 0.91 |  |  |
| CNB accuracy | Undergrowth CNVs | + childhood obesity PGS | 321 | 0 | -0.027 (0.029) | 0.34 | 0.049 |  |  |
| CNB accuracy | Undergrowth CNVs | + adult BMI PGS | 321 | 0 | -0.037 (0.030) | 0.22 | 0.48 |  |  |
| CNB accuracy | Undergrowth CNVs | + height PGS | 321 | 0 | -0.044 (0.029) | 0.13 | 0.26 |  |  |
| CNB accuracy | Undergrowth CNVs | + autism PGS (2026) | 321 | 0 | -0.041 (0.029) | 0.16 | 0.74 |  |  |
| CNB accuracy | Undergrowth CNVs | + all four PGS | 321 | 0 | -0.033 (0.030) | 0.28 |  |  |  |
| CNB accuracy | Undergrowth CNVs | + all four PGS + sex | 321 | 0 | -0.036 (0.046) | 0.44 |  |  |  |
| CNB accuracy | All CNVs | unadjusted (same individuals) | 532 | 0 | -0.009 (0.023) | 0.68 |  |  |  |
| CNB accuracy | All CNVs | + sex | 532 | 0 | +0.006 (0.035) | 0.85 | 0.54 |  |  |
| CNB accuracy | All CNVs | + childhood obesity PGS | 532 | 0 | -0.000 (0.023) | 0.99 | 0.2 |  |  |
| CNB accuracy | All CNVs | + adult BMI PGS | 532 | 0 | -0.006 (0.023) | 0.8 | 0.35 |  |  |
| CNB accuracy | All CNVs | + height PGS | 532 | 0 | -0.009 (0.023) | 0.69 | 0.81 |  |  |
| CNB accuracy | All CNVs | + autism PGS (2026) | 532 | 0 | -0.010 (0.023) | 0.65 | 0.84 |  |  |
| CNB accuracy | All CNVs | + all four PGS | 532 | 0 | -0.001 (0.023) | 0.95 |  |  |  |
| CNB accuracy | All CNVs | + all four PGS + sex | 532 | 0 | +0.015 (0.035) | 0.67 |  |  |  |
| CNB speed | 16p11.2 del | unadjusted (same individuals) | 103 | 0 | -0.022 (0.054) | 0.68 |  |  |  |
| CNB speed | 16p11.2 del | + sex | 103 | 0 | -0.010 (0.077) | 0.9 | 0.76 |  |  |
| CNB speed | 16p11.2 del | + childhood obesity PGS | 103 | 0 | -0.022 (0.055) | 0.68 | 0.61 |  |  |
| CNB speed | 16p11.2 del | + adult BMI PGS | 103 | 0 | -0.053 (0.056) | 0.35 | 0.13 |  |  |
| CNB speed | 16p11.2 del | + height PGS | 103 | 0 | -0.004 (0.056) | 0.95 | 0.83 |  |  |
| CNB speed | 16p11.2 del | + autism PGS (2026) | 103 | 0 | -0.031 (0.055) | 0.57 | 0.68 |  |  |
| CNB speed | 16p11.2 del | + all four PGS | 103 | 0 | -0.045 (0.059) | 0.45 |  |  |  |
| CNB speed | 16p11.2 del | + all four PGS + sex | 103 | 0 | +0.009 (0.085) | 0.92 |  |  |  |
| CNB speed | 16p11.2 dup | unadjusted (same individuals) |  39 | 0 | +0.143 (0.112) | 0.21 |  |  |  |
| CNB speed | 16p11.2 dup | + sex |  39 | 0 | +0.121 (0.180) | 0.5 | 0.87 |  |  |
| CNB speed | 16p11.2 dup | + childhood obesity PGS |  39 | 0 | -0.008 (0.105) | 0.94 | 0.007 |  |  |
| CNB speed | 16p11.2 dup | + adult BMI PGS |  39 | 0 | +0.195 (0.118) | 0.11 | 0.49 |  |  |
| CNB speed | 16p11.2 dup | + height PGS |  39 | 0 | +0.114 (0.106) | 0.29 | 0.06 |  |  |
| CNB speed | 16p11.2 dup | + autism PGS (2026) |  39 | 0 | +0.031 (0.112) | 0.78 | 0.018 |  |  |
| CNB speed | 16p11.2 dup | + all four PGS |  39 | 0 | -0.159 (0.126) | 0.22 |  |  |  |
| CNB speed | 16p11.2 dup | + all four PGS + sex |  39 | 0 | -0.120 (0.178) | 0.51 |  |  |  |
| CNB speed | 22q11.2 del | unadjusted (same individuals) | 287 | 0 | -0.013 (0.034) | 0.7 |  |  |  |
| CNB speed | 22q11.2 del | + sex | 287 | 0 | -0.101 (0.055) | 0.068 | 0.044 |  |  |
| CNB speed | 22q11.2 del | + childhood obesity PGS | 287 | 0 | -0.020 (0.035) | 0.58 | 0.43 |  |  |
| CNB speed | 22q11.2 del | + adult BMI PGS | 287 | 0 | -0.030 (0.036) | 0.41 | 0.55 |  |  |
| CNB speed | 22q11.2 del | + height PGS | 287 | 0 | -0.020 (0.034) | 0.57 | 0.06 |  |  |
| CNB speed | 22q11.2 del | + autism PGS (2026) | 287 | 0 | -0.014 (0.034) | 0.68 | 0.91 |  |  |
| CNB speed | 22q11.2 del | + all four PGS | 287 | 0 | -0.040 (0.037) | 0.28 |  |  |  |
| CNB speed | 22q11.2 del | + all four PGS + sex | 287 | 0 | -0.124 (0.057) | 0.031 |  |  |  |
| CNB speed | 22q11.2 dup | unadjusted (same individuals) | 110 | 0 | +0.046 (0.094) | 0.62 |  |  |  |
| CNB speed | 22q11.2 dup | + sex | 110 | 0 | -0.153 (0.183) | 0.41 | 0.22 |  |  |
| CNB speed | 22q11.2 dup | + childhood obesity PGS | 110 | 0 | +0.057 (0.096) | 0.56 | 0.79 |  |  |
| CNB speed | 22q11.2 dup | + adult BMI PGS | 110 | 0 | +0.055 (0.095) | 0.56 | 0.42 |  |  |
| CNB speed | 22q11.2 dup | + height PGS | 110 | 0 | +0.047 (0.096) | 0.63 | 0.96 |  |  |
| CNB speed | 22q11.2 dup | + autism PGS (2026) | 110 | 0 | +0.033 (0.100) | 0.74 | 0.95 |  |  |
| CNB speed | 22q11.2 dup | + all four PGS | 110 | 0 | +0.050 (0.106) | 0.64 |  |  |  |
| CNB speed | 22q11.2 dup | + all four PGS + sex | 110 | 0 | -0.276 (0.205) | 0.18 |  |  |  |
| CNB speed | Overgrowth CNVs | unadjusted (same individuals) | 213 | 0 | +0.013 (0.055) | 0.82 |  |  |  |
| CNB speed | Overgrowth CNVs | + sex | 213 | 0 | -0.049 (0.088) | 0.58 | 0.4 |  |  |
| CNB speed | Overgrowth CNVs | + childhood obesity PGS | 213 | 0 | +0.013 (0.056) | 0.82 | 0.68 |  |  |
| CNB speed | Overgrowth CNVs | + adult BMI PGS | 213 | 0 | +0.005 (0.056) | 0.93 | 0.14 |  |  |
| CNB speed | Overgrowth CNVs | + height PGS | 213 | 0 | +0.016 (0.057) | 0.78 | 1 |  |  |
| CNB speed | Overgrowth CNVs | + autism PGS (2026) | 213 | 0 | -0.002 (0.057) | 0.97 | 0.72 |  |  |
| CNB speed | Overgrowth CNVs | + all four PGS | 213 | 0 | -0.012 (0.059) | 0.84 |  |  |  |
| CNB speed | Overgrowth CNVs | + all four PGS + sex | 213 | 0 | -0.092 (0.093) | 0.32 |  |  |  |
| CNB speed | Undergrowth CNVs | unadjusted (same individuals) | 326 | 0 | +0.003 (0.032) | 0.93 |  |  |  |
| CNB speed | Undergrowth CNVs | + sex | 326 | 0 | -0.075 (0.051) | 0.15 | 0.052 |  |  |
| CNB speed | Undergrowth CNVs | + childhood obesity PGS | 326 | 0 | +0.004 (0.033) | 0.9 | 0.75 |  |  |
| CNB speed | Undergrowth CNVs | + adult BMI PGS | 326 | 0 | -0.007 (0.034) | 0.84 | 0.85 |  |  |
| CNB speed | Undergrowth CNVs | + height PGS | 326 | 0 | -0.001 (0.033) | 0.97 | 0.19 |  |  |
| CNB speed | Undergrowth CNVs | + autism PGS (2026) | 326 | 0 | +0.001 (0.032) | 0.98 | 0.24 |  |  |
| CNB speed | Undergrowth CNVs | + all four PGS | 326 | 0 | -0.012 (0.035) | 0.72 |  |  |  |
| CNB speed | Undergrowth CNVs | + all four PGS + sex | 326 | 0 | -0.089 (0.053) | 0.095 |  |  |  |
| CNB speed | All CNVs | unadjusted (same individuals) | 539 | 0 | +0.002 (0.029) | 0.95 |  |  |  |
| CNB speed | All CNVs | + sex | 539 | 0 | -0.074 (0.045) | 0.097 | 0.029 |  |  |
| CNB speed | All CNVs | + childhood obesity PGS | 539 | 0 | +0.003 (0.029) | 0.93 | 0.76 |  |  |
| CNB speed | All CNVs | + adult BMI PGS | 539 | 0 | -0.003 (0.030) | 0.93 | 0.45 |  |  |
| CNB speed | All CNVs | + height PGS | 539 | 0 | +0.002 (0.029) | 0.94 | 0.25 |  |  |
| CNB speed | All CNVs | + autism PGS (2026) | 539 | 0 | +0.003 (0.029) | 0.92 | 0.33 |  |  |
| CNB speed | All CNVs | + all four PGS | 539 | 0 | -0.000 (0.030) | 0.99 |  |  |  |
| CNB speed | All CNVs | + all four PGS + sex | 539 | 0 | -0.071 (0.046) | 0.12 |  |  |  |

- Table S6: power to detect modifier x BMI interactions within carriers. For each modifier (sex; each PGS as a continuous z-score) the interaction term BMI x modifier was fitted in GAM y ~ s(age) + BMI * modifier [+ locus]. Power (normal approximation from the term's SE) is given for an interaction as large as the carriers' own BMI effect on that trait, i.e. the BMI effect disappearing in one sex or changing by its full size per SD of PGS; 'detectable interaction' is the smallest effect detectable with 80% power (2.8 x SE). Observed interaction estimates and p-values are listed for completeness.

## Table S6. Power to detect modifier x BMI interactions

| Trait | Group | Modifier | n carriers | BMI slope (reference effect) | Interaction SE | Detectable interaction (80% power) | Power for interaction = BMI slope | Observed interaction (SE) | Interaction p |
|---|---|---|---|---|---|---|---|---|---|
| Height | 16p11.2 del | adult BMI PGS |  118 | +0.343 | 0.063 | 0.18 | 100% | +0.126 (0.063) | 0.048 |
| Height | 16p11.2 del | autism PGS (2026) |  118 | +0.348 | 0.088 | 0.25 | 98% | +0.184 (0.088) | 0.04 |
| Height | 16p11.2 del | childhood obesity PGS |  118 | +0.353 | 0.081 | 0.23 | 99% | -0.026 (0.081) | 0.75 |
| Height | 16p11.2 del | height PGS |  118 | +0.368 | 0.076 | 0.21 | 100% | +0.127 (0.076) | 0.098 |
| Height | 16p11.2 del | sex (male vs female) |  358 | +0.111 | 0.096 | 0.27 | 21% | -0.017 (0.096) | 0.86 |
| Height | 16p11.2 dup | adult BMI PGS |   67 | +0.052 | 0.135 | 0.38 | 7% | +0.254 (0.135) | 0.064 |
| Height | 16p11.2 dup | autism PGS (2026) |   67 | +0.108 | 0.109 | 0.30 | 17% | +0.034 (0.109) | 0.75 |
| Height | 16p11.2 dup | childhood obesity PGS |   67 | -0.006 | 0.109 | 0.30 | 5% | +0.216 (0.109) | 0.052 |
| Height | 16p11.2 dup | height PGS |   67 | +0.064 | 0.112 | 0.31 | 9% | +0.166 (0.112) | 0.14 |
| Height | 16p11.2 dup | sex (male vs female) |  197 | +0.052 | 0.137 | 0.38 | 7% | -0.072 (0.137) | 0.6 |
| Height | 22q11.2 del | adult BMI PGS |  268 | +0.086 | 0.055 | 0.15 | 35% | +0.024 (0.055) | 0.67 |
| Height | 22q11.2 del | autism PGS (2026) |  268 | +0.061 | 0.055 | 0.15 | 20% | -0.036 (0.055) | 0.52 |
| Height | 22q11.2 del | childhood obesity PGS |  268 | +0.054 | 0.054 | 0.15 | 17% | -0.013 (0.054) | 0.81 |
| Height | 22q11.2 del | height PGS |  268 | +0.054 | 0.045 | 0.12 | 22% | +0.052 (0.045) | 0.24 |
| Height | 22q11.2 del | sex (male vs female) |  355 | +0.087 | 0.089 | 0.25 | 16% | -0.023 (0.089) | 0.8 |
| Height | 22q11.2 dup | adult BMI PGS |   97 | +0.096 | 0.080 | 0.22 | 22% | +0.105 (0.080) | 0.19 |
| Height | 22q11.2 dup | autism PGS (2026) |   97 | +0.060 | 0.089 | 0.25 | 10% | -0.230 (0.089) | 0.011 |
| Height | 22q11.2 dup | childhood obesity PGS |   97 | +0.070 | 0.081 | 0.23 | 14% | +0.078 (0.081) | 0.34 |
| Height | 22q11.2 dup | height PGS |   97 | +0.137 | 0.095 | 0.27 | 30% | +0.015 (0.095) | 0.88 |
| Height | 22q11.2 dup | sex (male vs female) |  116 | -0.105 | 0.167 | 0.47 | 10% | +0.325 (0.167) | 0.054 |
| Height | Overgrowth CNVs | adult BMI PGS |  215 | +0.245 | 0.049 | 0.14 | 100% | +0.115 (0.049) | 0.019 |
| Height | Overgrowth CNVs | autism PGS (2026) |  215 | +0.268 | 0.059 | 0.16 | 100% | +0.001 (0.059) | 0.99 |
| Height | Overgrowth CNVs | childhood obesity PGS |  215 | +0.256 | 0.054 | 0.15 | 100% | -0.046 (0.054) | 0.39 |
| Height | Overgrowth CNVs | height PGS |  215 | +0.290 | 0.058 | 0.16 | 100% | +0.047 (0.058) | 0.42 |
| Height | Overgrowth CNVs | sex (male vs female) |  474 | +0.080 | 0.081 | 0.23 | 17% | +0.031 (0.081) | 0.71 |
| Height | Undergrowth CNVs | adult BMI PGS |  335 | +0.084 | 0.050 | 0.14 | 39% | +0.049 (0.050) | 0.33 |
| Height | Undergrowth CNVs | autism PGS (2026) |  335 | +0.067 | 0.048 | 0.13 | 29% | -0.028 (0.048) | 0.55 |
| Height | Undergrowth CNVs | childhood obesity PGS |  335 | +0.055 | 0.047 | 0.13 | 21% | +0.042 (0.047) | 0.37 |
| Height | Undergrowth CNVs | height PGS |  335 | +0.053 | 0.041 | 0.12 | 25% | +0.066 (0.041) | 0.11 |
| Height | Undergrowth CNVs | sex (male vs female) |  552 | +0.071 | 0.075 | 0.21 | 16% | -0.041 (0.075) | 0.58 |
| Height | All CNVs | adult BMI PGS |  550 | +0.124 | 0.037 | 0.10 | 92% | +0.088 (0.037) | 0.017 |
| Height | All CNVs | autism PGS (2026) |  550 | +0.124 | 0.037 | 0.10 | 92% | -0.020 (0.037) | 0.58 |
| Height | All CNVs | childhood obesity PGS |  550 | +0.121 | 0.036 | 0.10 | 92% | -0.025 (0.036) | 0.49 |
| Height | All CNVs | height PGS |  550 | +0.131 | 0.035 | 0.10 | 97% | +0.035 (0.035) | 0.31 |
| Height | All CNVs | sex (male vs female) | 1026 | +0.053 | 0.054 | 0.15 | 17% | +0.025 (0.054) | 0.65 |
| Head circumference | 16p11.2 del | adult BMI PGS |   98 | +0.501 | 0.065 | 0.18 | 100% | +0.116 (0.065) | 0.075 |
| Head circumference | 16p11.2 del | autism PGS (2026) |   98 | +0.486 | 0.092 | 0.26 | 100% | +0.071 (0.092) | 0.45 |
| Head circumference | 16p11.2 del | childhood obesity PGS |   98 | +0.496 | 0.083 | 0.23 | 100% | +0.109 (0.083) | 0.19 |
| Head circumference | 16p11.2 del | height PGS |   98 | +0.497 | 0.081 | 0.23 | 100% | +0.081 (0.081) | 0.32 |
| Head circumference | 16p11.2 del | sex (male vs female) |  211 | +0.564 | 0.114 | 0.32 | 100% | -0.146 (0.114) | 0.2 |
| Head circumference | 16p11.2 dup | adult BMI PGS |   54 | +0.093 | 0.139 | 0.39 | 10% | -0.013 (0.139) | 0.93 |
| Head circumference | 16p11.2 dup | autism PGS (2026) |   54 | +0.083 | 0.109 | 0.31 | 12% | +0.074 (0.109) | 0.5 |
| Head circumference | 16p11.2 dup | childhood obesity PGS |   54 | +0.112 | 0.129 | 0.36 | 14% | +0.025 (0.129) | 0.85 |
| Head circumference | 16p11.2 dup | height PGS |   54 | +0.070 | 0.153 | 0.43 | 7% | -0.168 (0.153) | 0.28 |
| Head circumference | 16p11.2 dup | sex (male vs female) |  123 | +0.309 | 0.162 | 0.45 | 48% | -0.056 (0.162) | 0.73 |
| Head circumference | 22q11.2 del | adult BMI PGS |  171 | +0.308 | 0.067 | 0.19 | 100% | +0.062 (0.067) | 0.35 |
| Head circumference | 22q11.2 del | autism PGS (2026) |  171 | +0.324 | 0.072 | 0.20 | 99% | +0.020 (0.072) | 0.79 |
| Head circumference | 22q11.2 del | childhood obesity PGS |  171 | +0.307 | 0.065 | 0.18 | 100% | -0.148 (0.065) | 0.024 |
| Head circumference | 22q11.2 del | height PGS |  171 | +0.357 | 0.062 | 0.17 | 100% | -0.137 (0.062) | 0.028 |
| Head circumference | 22q11.2 del | sex (male vs female) |  242 | +0.336 | 0.112 | 0.31 | 85% | +0.015 (0.112) | 0.89 |
| Head circumference | 22q11.2 dup | adult BMI PGS |   72 | +0.496 | 0.114 | 0.32 | 99% | +0.166 (0.114) | 0.15 |
| Head circumference | 22q11.2 dup | autism PGS (2026) |   72 | +0.518 | 0.122 | 0.34 | 99% | -0.015 (0.122) | 0.9 |
| Head circumference | 22q11.2 dup | childhood obesity PGS |   72 | +0.498 | 0.112 | 0.31 | 99% | +0.019 (0.112) | 0.87 |
| Head circumference | 22q11.2 dup | height PGS |   72 | +0.564 | 0.130 | 0.36 | 99% | -0.023 (0.130) | 0.86 |
| Head circumference | 22q11.2 dup | sex (male vs female) |   91 | +0.575 | 0.225 | 0.63 | 72% | -0.177 (0.225) | 0.43 |
| Head circumference | Overgrowth CNVs | adult BMI PGS |  170 | +0.492 | 0.056 | 0.16 | 100% | +0.125 (0.056) | 0.027 |
| Head circumference | Overgrowth CNVs | autism PGS (2026) |  170 | +0.497 | 0.068 | 0.19 | 100% | +0.024 (0.068) | 0.72 |
| Head circumference | Overgrowth CNVs | childhood obesity PGS |  170 | +0.485 | 0.062 | 0.17 | 100% | +0.060 (0.062) | 0.33 |
| Head circumference | Overgrowth CNVs | height PGS |  170 | +0.514 | 0.067 | 0.19 | 100% | +0.031 (0.067) | 0.64 |
| Head circumference | Overgrowth CNVs | sex (male vs female) |  302 | +0.571 | 0.098 | 0.27 | 100% | -0.176 (0.098) | 0.072 |
| Head circumference | Undergrowth CNVs | adult BMI PGS |  225 | +0.253 | 0.060 | 0.17 | 99% | +0.064 (0.060) | 0.29 |
| Head circumference | Undergrowth CNVs | autism PGS (2026) |  225 | +0.273 | 0.059 | 0.17 | 100% | +0.015 (0.059) | 0.8 |
| Head circumference | Undergrowth CNVs | childhood obesity PGS |  225 | +0.266 | 0.056 | 0.16 | 100% | -0.116 (0.056) | 0.038 |
| Head circumference | Undergrowth CNVs | height PGS |  225 | +0.287 | 0.057 | 0.16 | 100% | -0.113 (0.057) | 0.048 |
| Head circumference | Undergrowth CNVs | sex (male vs female) |  365 | +0.322 | 0.091 | 0.26 | 94% | +0.001 (0.091) | 0.99 |
| Head circumference | All CNVs | adult BMI PGS |  395 | +0.357 | 0.041 | 0.11 | 100% | +0.089 (0.041) | 0.031 |
| Head circumference | All CNVs | autism PGS (2026) |  395 | +0.367 | 0.041 | 0.12 | 100% | +0.002 (0.041) | 0.96 |
| Head circumference | All CNVs | childhood obesity PGS |  395 | +0.358 | 0.041 | 0.11 | 100% | -0.042 (0.041) | 0.3 |
| Head circumference | All CNVs | height PGS |  395 | +0.371 | 0.043 | 0.12 | 100% | -0.046 (0.043) | 0.28 |
| Head circumference | All CNVs | sex (male vs female) |  667 | +0.408 | 0.063 | 0.18 | 100% | -0.035 (0.063) | 0.58 |
| SRS | 16p11.2 del | adult BMI PGS |   74 | +0.250 | 0.248 | 0.69 | 17% | -0.338 (0.248) | 0.18 |
| SRS | 16p11.2 del | autism PGS (2026) |   74 | +0.185 | 0.221 | 0.62 | 13% | -0.427 (0.221) | 0.058 |
| SRS | 16p11.2 del | childhood obesity PGS |   74 | +0.157 | 0.230 | 0.64 | 11% | +0.215 (0.230) | 0.35 |
| SRS | 16p11.2 del | height PGS |   74 | +0.128 | 0.208 | 0.58 | 9% | +0.237 (0.208) | 0.26 |
| SRS | 16p11.2 del | sex (male vs female) |   79 | +0.225 | 0.392 | 1.10 | 9% | -0.102 (0.392) | 0.8 |
| SRS | 16p11.2 dup | adult BMI PGS |   29 | +0.378 | 0.437 | 1.22 | 14% | +0.645 (0.437) | 0.15 |
| SRS | 16p11.2 dup | autism PGS (2026) |   29 | +0.277 | 0.283 | 0.79 | 16% | +0.393 (0.283) | 0.18 |
| SRS | 16p11.2 dup | childhood obesity PGS |   29 | -0.072 | 0.364 | 1.02 | 5% | +0.552 (0.364) | 0.14 |
| SRS | 16p11.2 dup | height PGS |   29 | +0.252 | 0.437 | 1.22 | 9% | -0.387 (0.437) | 0.39 |
| SRS | 16p11.2 dup | sex (male vs female) |   34 | +0.726 | 0.559 | 1.56 | 25% | -0.109 (0.559) | 0.85 |
| SRS | 22q11.2 del | adult BMI PGS |  186 | +0.338 | 0.117 | 0.33 | 82% | -0.062 (0.117) | 0.6 |
| SRS | 22q11.2 del | autism PGS (2026) |  186 | +0.294 | 0.110 | 0.31 | 76% | -0.197 (0.110) | 0.074 |
| SRS | 22q11.2 del | childhood obesity PGS |  186 | +0.305 | 0.114 | 0.32 | 76% | -0.033 (0.114) | 0.78 |
| SRS | 22q11.2 del | height PGS |  186 | +0.318 | 0.104 | 0.29 | 86% | +0.025 (0.104) | 0.81 |
| SRS | 22q11.2 del | sex (male vs female) |  207 | +0.300 | 0.204 | 0.57 | 31% | -0.105 (0.204) | 0.61 |
| SRS | 22q11.2 dup | adult BMI PGS |   60 | -0.129 | 0.308 | 0.86 | 7% | -0.066 (0.308) | 0.83 |
| SRS | 22q11.2 dup | autism PGS (2026) |   60 | -0.167 | 0.294 | 0.82 | 9% | -0.136 (0.294) | 0.65 |
| SRS | 22q11.2 dup | childhood obesity PGS |   60 | -0.151 | 0.278 | 0.78 | 8% | +0.379 (0.278) | 0.18 |
| SRS | 22q11.2 dup | height PGS |   60 | -0.150 | 0.280 | 0.78 | 8% | +0.074 (0.280) | 0.79 |
| SRS | 22q11.2 dup | sex (male vs female) |   66 | +0.241 | 0.534 | 1.50 | 7% | -0.411 (0.534) | 0.44 |
| SRS | Overgrowth CNVs | adult BMI PGS |  134 | +0.131 | 0.184 | 0.52 | 11% | -0.279 (0.184) | 0.13 |
| SRS | Overgrowth CNVs | autism PGS (2026) |  134 | +0.073 | 0.150 | 0.42 | 8% | -0.338 (0.150) | 0.026 |
| SRS | Overgrowth CNVs | childhood obesity PGS |  134 | +0.039 | 0.151 | 0.42 | 6% | +0.377 (0.151) | 0.014 |
| SRS | Overgrowth CNVs | height PGS |  134 | +0.059 | 0.152 | 0.43 | 7% | +0.175 (0.152) | 0.25 |
| SRS | Overgrowth CNVs | sex (male vs female) |  145 | +0.435 | 0.295 | 0.83 | 31% | -0.506 (0.295) | 0.089 |
| SRS | Undergrowth CNVs | adult BMI PGS |  215 | +0.331 | 0.114 | 0.32 | 83% | +0.009 (0.114) | 0.93 |
| SRS | Undergrowth CNVs | autism PGS (2026) |  215 | +0.303 | 0.100 | 0.28 | 85% | -0.134 (0.100) | 0.18 |
| SRS | Undergrowth CNVs | childhood obesity PGS |  215 | +0.292 | 0.108 | 0.30 | 77% | +0.043 (0.108) | 0.69 |
| SRS | Undergrowth CNVs | height PGS |  215 | +0.303 | 0.102 | 0.29 | 84% | +0.007 (0.102) | 0.95 |
| SRS | Undergrowth CNVs | sex (male vs female) |  241 | +0.374 | 0.190 | 0.53 | 50% | -0.145 (0.190) | 0.45 |
| SRS | All CNVs | adult BMI PGS |  349 | +0.226 | 0.093 | 0.26 | 68% | -0.059 (0.093) | 0.53 |
| SRS | All CNVs | autism PGS (2026) |  349 | +0.224 | 0.078 | 0.22 | 81% | -0.184 (0.078) | 0.019 |
| SRS | All CNVs | childhood obesity PGS |  349 | +0.190 | 0.085 | 0.24 | 61% | +0.167 (0.085) | 0.052 |
| SRS | All CNVs | height PGS |  349 | +0.219 | 0.083 | 0.23 | 75% | +0.064 (0.083) | 0.44 |
| SRS | All CNVs | sex (male vs female) |  386 | +0.422 | 0.157 | 0.44 | 77% | -0.311 (0.157) | 0.048 |
| SCQ | 16p11.2 del | adult BMI PGS |   56 | +0.484 | 0.234 | 0.66 | 54% | +0.017 (0.234) | 0.94 |
| SCQ | 16p11.2 del | autism PGS (2026) |   56 | +0.397 | 0.336 | 0.94 | 22% | -0.251 (0.336) | 0.46 |
| SCQ | 16p11.2 del | childhood obesity PGS |   56 | +0.283 | 0.425 | 1.19 | 10% | +0.406 (0.425) | 0.34 |
| SCQ | 16p11.2 del | height PGS |   56 | +0.352 | 0.364 | 1.02 | 16% | +0.389 (0.364) | 0.29 |
| SCQ | 16p11.2 del | sex (male vs female) |   60 | +0.720 | 0.637 | 1.78 | 20% | -0.557 (0.637) | 0.39 |
| SCQ | 16p11.2 dup | sex (male vs female) |   21 | +0.324 | 1.073 | 3.00 | 6% | +0.793 (1.073) | 0.47 |
| SCQ | 22q11.2 del | adult BMI PGS |  151 | +0.387 | 0.186 | 0.52 | 55% | -0.169 (0.186) | 0.36 |
| SCQ | 22q11.2 del | autism PGS (2026) |  151 | +0.354 | 0.182 | 0.51 | 49% | -0.075 (0.182) | 0.68 |
| SCQ | 22q11.2 del | childhood obesity PGS |  151 | +0.404 | 0.169 | 0.47 | 67% | -0.002 (0.169) | 0.99 |
| SCQ | 22q11.2 del | height PGS |  151 | +0.368 | 0.159 | 0.45 | 64% | +0.027 (0.159) | 0.87 |
| SCQ | 22q11.2 del | sex (male vs female) |  169 | +0.298 | 0.312 | 0.87 | 16% | -0.103 (0.312) | 0.74 |
| SCQ | 22q11.2 dup | adult BMI PGS |   44 | +0.425 | 0.545 | 1.53 | 12% | -0.400 (0.545) | 0.47 |
| SCQ | 22q11.2 dup | autism PGS (2026) |   44 | +0.036 | 0.418 | 1.17 | 5% | +0.670 (0.418) | 0.12 |
| SCQ | 22q11.2 dup | childhood obesity PGS |   44 | +0.286 | 0.453 | 1.27 | 10% | -0.565 (0.453) | 0.22 |
| SCQ | 22q11.2 dup | height PGS |   44 | +0.150 | 0.531 | 1.49 | 6% | -0.357 (0.531) | 0.51 |
| SCQ | 22q11.2 dup | sex (male vs female) |   50 | +0.420 | 0.801 | 2.24 | 8% | -0.080 (0.801) | 0.92 |
| SCQ | Overgrowth CNVs | adult BMI PGS |  100 | +0.348 | 0.208 | 0.58 | 39% | -0.024 (0.208) | 0.91 |
| SCQ | Overgrowth CNVs | autism PGS (2026) |  100 | +0.316 | 0.235 | 0.66 | 27% | -0.035 (0.235) | 0.88 |
| SCQ | Overgrowth CNVs | childhood obesity PGS |  100 | +0.306 | 0.263 | 0.74 | 21% | +0.196 (0.263) | 0.46 |
| SCQ | Overgrowth CNVs | height PGS |  100 | +0.354 | 0.279 | 0.78 | 25% | -0.058 (0.279) | 0.84 |
| SCQ | Overgrowth CNVs | sex (male vs female) |  110 | +0.665 | 0.466 | 1.31 | 30% | -0.439 (0.466) | 0.35 |
| SCQ | Undergrowth CNVs | adult BMI PGS |  168 | +0.366 | 0.185 | 0.52 | 50% | -0.133 (0.185) | 0.47 |
| SCQ | Undergrowth CNVs | autism PGS (2026) |  168 | +0.330 | 0.163 | 0.46 | 53% | -0.092 (0.163) | 0.57 |
| SCQ | Undergrowth CNVs | childhood obesity PGS |  168 | +0.354 | 0.165 | 0.46 | 57% | -0.037 (0.165) | 0.82 |
| SCQ | Undergrowth CNVs | height PGS |  168 | +0.330 | 0.158 | 0.44 | 55% | +0.027 (0.158) | 0.86 |
| SCQ | Undergrowth CNVs | sex (male vs female) |  190 | +0.295 | 0.296 | 0.83 | 17% | -0.062 (0.296) | 0.83 |
| SCQ | All CNVs | adult BMI PGS |  268 | +0.322 | 0.137 | 0.38 | 65% | -0.072 (0.137) | 0.6 |
| SCQ | All CNVs | autism PGS (2026) |  268 | +0.316 | 0.126 | 0.35 | 71% | -0.030 (0.126) | 0.81 |
| SCQ | All CNVs | childhood obesity PGS |  268 | +0.312 | 0.141 | 0.39 | 60% | +0.065 (0.141) | 0.64 |
| SCQ | All CNVs | height PGS |  268 | +0.311 | 0.138 | 0.39 | 62% | -0.002 (0.138) | 0.99 |
| SCQ | All CNVs | sex (male vs female) |  300 | +0.444 | 0.249 | 0.70 | 43% | -0.242 (0.249) | 0.33 |
| CNB accuracy | 16p11.2 del | adult BMI PGS |  101 | +0.007 | 0.051 | 0.14 | 5% | +0.040 (0.051) | 0.43 |
| CNB accuracy | 16p11.2 del | autism PGS (2026) |  101 | +0.009 | 0.048 | 0.13 | 5% | -0.037 (0.048) | 0.44 |
| CNB accuracy | 16p11.2 del | childhood obesity PGS |  101 | +0.014 | 0.049 | 0.14 | 6% | -0.038 (0.049) | 0.43 |
| CNB accuracy | 16p11.2 del | height PGS |  101 | +0.013 | 0.049 | 0.14 | 6% | -0.067 (0.049) | 0.17 |
| CNB accuracy | 16p11.2 del | sex (male vs female) |  140 | +0.066 | 0.081 | 0.23 | 13% | -0.065 (0.081) | 0.43 |
| CNB accuracy | 16p11.2 dup | adult BMI PGS |   39 | -0.155 | 0.109 | 0.30 | 30% | -0.063 (0.109) | 0.57 |
| CNB accuracy | 16p11.2 dup | autism PGS (2026) |   39 | -0.116 | 0.070 | 0.20 | 38% | -0.084 (0.070) | 0.24 |
| CNB accuracy | 16p11.2 dup | childhood obesity PGS |   39 | -0.193 | 0.089 | 0.25 | 58% | +0.117 (0.089) | 0.2 |
| CNB accuracy | 16p11.2 dup | height PGS |   39 | -0.153 | 0.105 | 0.29 | 31% | -0.086 (0.105) | 0.42 |
| CNB accuracy | 16p11.2 dup | sex (male vs female) |   55 | -0.132 | 0.160 | 0.45 | 13% | +0.060 (0.160) | 0.71 |
| CNB accuracy | 22q11.2 del | adult BMI PGS |  282 | -0.031 | 0.034 | 0.09 | 15% | +0.021 (0.034) | 0.54 |
| CNB accuracy | 22q11.2 del | autism PGS (2026) |  282 | -0.032 | 0.032 | 0.09 | 16% | +0.003 (0.032) | 0.92 |
| CNB accuracy | 22q11.2 del | childhood obesity PGS |  282 | -0.008 | 0.029 | 0.08 | 6% | +0.062 (0.029) | 0.032 |
| CNB accuracy | 22q11.2 del | height PGS |  282 | -0.036 | 0.029 | 0.08 | 24% | +0.036 (0.029) | 0.22 |
| CNB accuracy | 22q11.2 del | sex (male vs female) |  344 | +0.010 | 0.059 | 0.17 | 5% | -0.041 (0.059) | 0.49 |
| CNB accuracy | 22q11.2 dup | adult BMI PGS |  110 | +0.081 | 0.047 | 0.13 | 41% | +0.014 (0.047) | 0.77 |
| CNB accuracy | 22q11.2 dup | autism PGS (2026) |  110 | +0.062 | 0.070 | 0.20 | 14% | +0.026 (0.070) | 0.71 |
| CNB accuracy | 22q11.2 dup | childhood obesity PGS |  110 | +0.096 | 0.061 | 0.17 | 35% | -0.021 (0.061) | 0.74 |
| CNB accuracy | 22q11.2 dup | height PGS |  110 | +0.083 | 0.079 | 0.22 | 18% | +0.019 (0.079) | 0.81 |
| CNB accuracy | 22q11.2 dup | sex (male vs female) |  136 | +0.261 | 0.118 | 0.33 | 60% | -0.278 (0.118) | 0.021 |
| CNB accuracy | Overgrowth CNVs | adult BMI PGS |  211 | +0.049 | 0.032 | 0.09 | 33% | +0.014 (0.032) | 0.67 |
| CNB accuracy | Overgrowth CNVs | autism PGS (2026) |  211 | +0.028 | 0.037 | 0.10 | 12% | -0.033 (0.037) | 0.38 |
| CNB accuracy | Overgrowth CNVs | childhood obesity PGS |  211 | +0.053 | 0.035 | 0.10 | 33% | -0.016 (0.035) | 0.66 |
| CNB accuracy | Overgrowth CNVs | height PGS |  211 | +0.046 | 0.042 | 0.12 | 19% | -0.046 (0.042) | 0.28 |
| CNB accuracy | Overgrowth CNVs | sex (male vs female) |  276 | +0.123 | 0.065 | 0.18 | 47% | -0.122 (0.065) | 0.064 |
| CNB accuracy | Undergrowth CNVs | adult BMI PGS |  321 | -0.037 | 0.032 | 0.09 | 22% | +0.022 (0.032) | 0.48 |
| CNB accuracy | Undergrowth CNVs | autism PGS (2026) |  321 | -0.041 | 0.029 | 0.08 | 29% | -0.010 (0.029) | 0.74 |
| CNB accuracy | Undergrowth CNVs | childhood obesity PGS |  321 | -0.027 | 0.027 | 0.08 | 17% | +0.053 (0.027) | 0.049 |
| CNB accuracy | Undergrowth CNVs | height PGS |  321 | -0.044 | 0.028 | 0.08 | 35% | +0.032 (0.028) | 0.26 |
| CNB accuracy | Undergrowth CNVs | sex (male vs female) |  399 | -0.012 | 0.055 | 0.15 | 6% | -0.028 (0.055) | 0.61 |
| CNB accuracy | All CNVs | adult BMI PGS |  532 | -0.006 | 0.022 | 0.06 | 6% | +0.021 (0.022) | 0.35 |
| CNB accuracy | All CNVs | autism PGS (2026) |  532 | -0.010 | 0.022 | 0.06 | 8% | -0.004 (0.022) | 0.84 |
| CNB accuracy | All CNVs | childhood obesity PGS |  532 | -0.000 | 0.021 | 0.06 | 5% | +0.027 (0.021) | 0.2 |
| CNB accuracy | All CNVs | height PGS |  532 | -0.009 | 0.022 | 0.06 | 7% | +0.005 (0.022) | 0.81 |
| CNB accuracy | All CNVs | sex (male vs female) |  675 | +0.025 | 0.040 | 0.11 | 9% | -0.045 (0.040) | 0.26 |
| CNB speed | 16p11.2 del | adult BMI PGS |  103 | -0.053 | 0.060 | 0.17 | 14% | +0.092 (0.060) | 0.13 |
| CNB speed | 16p11.2 del | autism PGS (2026) |  103 | -0.031 | 0.058 | 0.16 | 8% | -0.024 (0.058) | 0.68 |
| CNB speed | 16p11.2 del | childhood obesity PGS |  103 | -0.022 | 0.058 | 0.16 | 7% | -0.030 (0.058) | 0.61 |
| CNB speed | 16p11.2 del | height PGS |  103 | -0.004 | 0.059 | 0.16 | 5% | -0.013 (0.059) | 0.83 |
| CNB speed | 16p11.2 del | sex (male vs female) |  142 | +0.032 | 0.089 | 0.25 | 6% | -0.078 (0.089) | 0.38 |
| CNB speed | 16p11.2 dup | adult BMI PGS |   39 | +0.195 | 0.153 | 0.43 | 25% | +0.108 (0.153) | 0.49 |
| CNB speed | 16p11.2 dup | autism PGS (2026) |   39 | +0.031 | 0.095 | 0.27 | 6% | +0.236 (0.095) | 0.018 |
| CNB speed | 16p11.2 dup | childhood obesity PGS |   39 | -0.008 | 0.111 | 0.31 | 5% | +0.322 (0.111) | 0.007 |
| CNB speed | 16p11.2 dup | height PGS |   39 | +0.114 | 0.142 | 0.40 | 13% | -0.278 (0.142) | 0.06 |
| CNB speed | 16p11.2 dup | sex (male vs female) |   55 | -0.066 | 0.214 | 0.60 | 6% | +0.113 (0.214) | 0.6 |
| CNB speed | 22q11.2 del | adult BMI PGS |  287 | -0.030 | 0.035 | 0.10 | 14% | -0.021 (0.035) | 0.55 |
| CNB speed | 22q11.2 del | autism PGS (2026) |  287 | -0.014 | 0.037 | 0.10 | 7% | +0.004 (0.037) | 0.91 |
| CNB speed | 22q11.2 del | childhood obesity PGS |  287 | -0.020 | 0.033 | 0.09 | 9% | -0.026 (0.033) | 0.43 |
| CNB speed | 22q11.2 del | height PGS |  287 | -0.020 | 0.033 | 0.09 | 9% | +0.062 (0.033) | 0.06 |
| CNB speed | 22q11.2 del | sex (male vs female) |  350 | -0.080 | 0.068 | 0.19 | 22% | +0.116 (0.068) | 0.089 |
| CNB speed | 22q11.2 dup | adult BMI PGS |  110 | +0.055 | 0.077 | 0.22 | 11% | +0.062 (0.077) | 0.42 |
| CNB speed | 22q11.2 dup | autism PGS (2026) |  110 | +0.033 | 0.116 | 0.33 | 6% | +0.008 (0.116) | 0.95 |
| CNB speed | 22q11.2 dup | childhood obesity PGS |  110 | +0.057 | 0.101 | 0.28 | 9% | -0.027 (0.101) | 0.79 |
| CNB speed | 22q11.2 dup | height PGS |  110 | +0.047 | 0.130 | 0.36 | 7% | -0.006 (0.130) | 0.96 |
| CNB speed | 22q11.2 dup | sex (male vs female) |  136 | -0.075 | 0.177 | 0.50 | 7% | +0.131 (0.177) | 0.46 |
| CNB speed | Overgrowth CNVs | adult BMI PGS |  213 | +0.005 | 0.048 | 0.14 | 5% | +0.072 (0.048) | 0.14 |
| CNB speed | Overgrowth CNVs | autism PGS (2026) |  213 | -0.002 | 0.057 | 0.16 | 5% | -0.020 (0.057) | 0.72 |
| CNB speed | Overgrowth CNVs | childhood obesity PGS |  213 | +0.013 | 0.053 | 0.15 | 6% | -0.022 (0.053) | 0.68 |
| CNB speed | Overgrowth CNVs | height PGS |  213 | +0.016 | 0.064 | 0.18 | 6% | -0.000 (0.064) | 1 |
| CNB speed | Overgrowth CNVs | sex (male vs female) |  278 | -0.012 | 0.088 | 0.25 | 5% | +0.033 (0.088) | 0.71 |
| CNB speed | Undergrowth CNVs | adult BMI PGS |  326 | -0.007 | 0.033 | 0.09 | 6% | -0.006 (0.033) | 0.85 |
| CNB speed | Undergrowth CNVs | autism PGS (2026) |  326 | +0.001 | 0.033 | 0.09 | 5% | +0.039 (0.033) | 0.24 |
| CNB speed | Undergrowth CNVs | childhood obesity PGS |  326 | +0.004 | 0.031 | 0.09 | 5% | +0.010 (0.031) | 0.75 |
| CNB speed | Undergrowth CNVs | height PGS |  326 | -0.001 | 0.032 | 0.09 | 5% | +0.042 (0.032) | 0.19 |
| CNB speed | Undergrowth CNVs | sex (male vs female) |  405 | -0.079 | 0.064 | 0.18 | 23% | +0.116 (0.064) | 0.072 |
| CNB speed | All CNVs | adult BMI PGS |  539 | -0.003 | 0.027 | 0.08 | 5% | +0.021 (0.027) | 0.45 |
| CNB speed | All CNVs | autism PGS (2026) |  539 | +0.003 | 0.028 | 0.08 | 5% | +0.027 (0.028) | 0.33 |
| CNB speed | All CNVs | childhood obesity PGS |  539 | +0.003 | 0.027 | 0.08 | 5% | -0.008 (0.027) | 0.76 |
| CNB speed | All CNVs | height PGS |  539 | +0.002 | 0.029 | 0.08 | 5% | +0.033 (0.029) | 0.25 |
| CNB speed | All CNVs | sex (male vs female) |  683 | -0.057 | 0.050 | 0.14 | 21% | +0.082 (0.050) | 0.1 |

