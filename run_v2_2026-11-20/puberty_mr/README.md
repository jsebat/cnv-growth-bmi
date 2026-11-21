# Childhood body size, puberty and adult height: Mendelian-randomisation block (work package E, 2026-11-20)

Lower part of main Figure 1 (panels i-l; tags set by `TAGS` in `88_fig_puberty_mr_v2.R`). All numbers are from UK Biobank
(European-ancestry cluster, QC'd, as in the earlier scripts). Individual data stayed on Expanse; only aggregates are here.

## Files

- `block_puberty_mr.png` / `.pdf` (20 x 7 in), `block_puberty_mr.csv` (every plotted number), `panel_i.png` ... `panel_l.png` (no tags).
- `block_puberty_mr_ukb_trained_instrument.png` and `panel_k_ukb_trained_instrument.png`: same block with the path diagram drawn
  from the menarche score used in the previous draft (weights trained in UK Biobank). Numbers differ in the third decimal only.
- `mr_checks.csv`: 19 checks, each with the numbers and whether the conclusion holds.
- `expanse_results/`: aggregate tables from Expanse (`panel_j_means.csv`, `panel_j_trends.csv`, `panel_k_menarche.csv`,
  `mr_onesample.csv`, `mr_twosample.csv`, `mr_twosample_snps.csv`, `mr_controls.csv`, `mr_sample_counts.csv`,
  `mr_score_correlations.csv`, `reproduction_check.csv`, logs). File names keep the first lettering (j, k).
- Scripts (in the code folder): `88_fig_puberty_mr_v2.R` (figure, local), `88a_puberty_mr_ukbb_expanse.R` (analysis, Expanse),
  `88b_puberty_mr_select_instruments_expanse.R` (SNP instruments, Expanse); runners in `expanse_scripts/`.
- Expanse run folder: `/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120/`
  (`01_verify` job 55917353, `02_instruments` + `03_mr` job 55917415, `03_mr_r2` job 55917464 = final, two variants added;
  `00_external_sumstats/perry2014` = the one download, 23 MB, straight onto Expanse).

## Data and n

| Set | n |
|---|---|
| Body-size panel (i): height, BMI and a body-size-at-10 answer | 16p11.2 del 106, 16p11.2 dup 140, 22q11.2 dup 262, non-carriers 444,852 |
| Women with age at menarche (j) | 16p11.2 del 51, 16p11.2 dup 75, 22q11.2 dup 121, non-carriers 238,681 |
| MR sample (k, l): non-carrier women with menarche, height, all scores | 238,681 (235,757 with a body-size answer) |
| Non-carrier men (reduced forms only) | 207,136 |

22q11.2 deletion: 5 carriers, not analysed.

## What kind of Mendelian randomisation this is

The analysis in the draft is a **one-sample MR inside UK Biobank with genome-wide polygenic scores as single instruments**
(SBayesRC weights over millions of SNPs). It is not a SNP-level MR, and the first step (score -> menarche) is a score-outcome
association, not an effect per unit of childhood body size.

| Score | Source GWAS | UK Biobank in the GWAS? |
|---|---|---|
| Childhood obesity | Bradfield et al. 2019 (PMID 31504550), EGG trans-ancestral | No |
| Childhood BMI | Vogelezang et al. 2020 (PMID 33045005), EGG | No |
| Age at menarche (draft) | Guan et al. 2025 (PMID 40065166), direct-effect column, effective N about 12,500 | **Yes, UK Biobank only: complete overlap** |
| Adult BMI (covariate / comparison) | Yengo et al. 2018 (PMID 30124842) | Yes |
| Height (covariate) | Yengo et al. 2022 (PMID 36224396) | Yes |
| Age at menarche (new) | Perry et al. 2014 (PMID 25231870), ReproGen, up to 182,416 women | No |

New in this pass: (1) a menarche allele score from 92 independent Perry 2014 SNPs with external weights; (2) a two-sample
SNP-level MR in which SNP-exposure effects come from the GWAS without UK Biobank and SNP-outcome effects are estimated in the
238,681 women (inverse-variance weighted, weighted median, MR-Egger, multivariable IVW). The path diagram uses (1).

## Models

- Height: residual of height on sex, age, age squared and 10 PCs, in SD (as before). Body-size panel: Box-Cox within sex (as before).
- Scores residualised on 10 PCs and scaled to SD. Covariates in every MR regression: height score, age, 10 PCs.
- Score -> menarche: linear regression. Menarche -> height and the direct effect: two-stage least squares with menarche
  instrumented by the menarche score and the childhood-obesity score as covariate. Through menarche = product of the two, delta-method SE.
- SNP instruments: P < 5e-8, one per 1 Mb, then r2 < 0.01 in UK Biobank; palindromic SNPs with MAF > 0.42 dropped.

## Results (non-carrier women, n = 238,681; `expanse_results/mr_onesample.csv`, `mr_twosample.csv`)

| Step | Estimate (95% CI) | P | Instrument strength |
|---|---|---|---|
| Childhood-obesity score -> age at menarche | -0.12 years per SD of score (-0.13 to -0.12) | < 1e-300 | F 1,417 |
| Age at menarche -> adult height (Perry allele score) | +0.096 SD per year later (0.081 to 0.110) | 3e-39 | F 6,246 |
| Direct effect of the score on height | +0.014 SD per SD (0.010 to 0.018) | 8e-13 | |
| Through menarche | -0.012 SD per SD (-0.014 to -0.010) | 3e-35 | |
| Total | +0.002 SD per SD (-0.001 to 0.005) | 0.21 | |
| Two-sample: childhood BMI -> menarche, IVW (16 SNPs) | -0.59 years per SD (-0.76 to -0.42) | 2e-11 | mean F 48 |
| Two-sample: menarche -> height, IVW (92 SNPs) | +0.18 SD per year (0.13 to 0.23) | 3e-12 | mean F 53 |
| Multivariable: childhood BMI -> height, direct (95 SNPs) | +0.10 SD per SD (0.01 to 0.19) | 0.028 | |

Carriers (observational; `panel_j_trends.csv`, `panel_k_menarche.csv`): 16p11.2 del -0.28 SD of height per category plumper
(p = 0.036, n = 106), 16p11.2 dup -0.17 (p = 0.20, n = 140), 22q11.2 dup -0.18 (p = 0.027, n = 262), non-carriers -0.07.
Menarche: 16p11.2 del 1.3 years earlier (p = 1e-8), 16p11.2 dup 0.9 years later (p = 3e-7), 22q11.2 dup no difference.
All previously reported numbers reproduce exactly (`reproduction_check.csv`).

## What would not survive review as written

1. **The menarche score in the draft was trained in UK Biobank and used as an instrument in UK Biobank.** The result does not
   depend on it (Perry allele score: 0.096 vs 0.094 SD per year), so the fix is to report the Perry-based numbers. The supplementary
   table note says the Guan GWAS "includes" UK Biobank; it is UK Biobank only.
2. **"The direct effect is offset by the effect through puberty" holds only with the height-score covariate.** Without it the
   total effect of the childhood-obesity score on height is +0.008 SD (p = 5e-5); with the adult-BMI score added it is +0.014
   (p = 1e-14). The direct effect is positive and the path through menarche negative in every specification; "partly offset" is safe, "cancel" is not.
3. **The size of menarche -> height is not pinned down.** 0.09-0.10 SD per year with the height score adjusted, 0.18-0.21
   without (one-sample and two-sample agree), 0.46 by MR-Egger with a significant intercept (p = 0.002) and very strong
   heterogeneity. Adjusting for a height score is unusual in MR and can remove part of the effect being estimated; a reviewer may ask for the unadjusted number as primary.
4. **Genome-wide scores are not clean instruments.** The childhood-obesity and childhood-BMI scores are associated with the
   mother smoking around birth (p = 2e-4, 0.002), a variable fixed before the participant's birth; the menarche instruments
   also lower recalled childhood body size (p = 8e-127). The SNP-level and multivariable analyses address this only in part (16 and 9 SNPs for childhood adiposity).
5. Step 1 is per SD of a score, not per unit of childhood body size; scaling by the recalled 3-level answer gives an implausible
   -1.3 years per category (observational -0.42) and should not be reported.
6. The MR is in non-carrier women only. It says nothing directly about carriers (51-121 women per genotype), and the carrier
   trends in panel i are both sexes combined; within women alone none reaches p < 0.05 (`panel_j_trends.csv`, rows sexes = women).
7. Not done: exclusion of the UK Biobank participants behind the Guan weights (not identifiable; replaced by the Perry score),
   Steiger / conditional-F statistics for the multivariable model, and a height GWAS without UK Biobank as outcome.
