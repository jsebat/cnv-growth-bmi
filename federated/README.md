# federated: run the BMI-modifier analyses at your site and send back tables only

Version 1.0.0, written 2026-11-20. Companion to Sacks et al. 2026 (the "BMI modifier" paper).

## Purpose

The paper asks whether BMI changes the effect of four CNVs (16p11.2 BP4-BP5 deletion and duplication, 22q11.2 A-D
deletion and duplication) on height, head circumference and autistic traits. This folder lets another site run the
paper's non-genetic analyses on its own carriers and non-carrier controls. The site keeps its data. It sends back
**aggregate tables only**. The tables from all sites are then combined.

No polygenic scores are needed. A paediatric cohort with few adults is fine.

## What you need

- R (4.2 or later) with packages `data.table`, `mgcv`, `MASS`, `jsonlite`; `lmerTest` for the longitudinal part;
  `ordinal` for one optional check; `ggplot2` and `patchwork` only for the figures.
- This folder, complete. The scripts must stay together.
- One table of visits (below) and a config file (copy `config_example.json`).

## Commands

```
# 1. at each site (about 2 to 5 minutes)
Rscript site_export.R --config my_config.json

# 2. at the coordinating site, on the folders the sites sent back
Rscript meta_combine.R --in siteA_out siteB_out lab_out --out meta_out

# 3. figures
Rscript plot_federated.R --meta meta_out --out figures

# test on synthetic data (about 3 minutes; writes to a temporary folder)
Rscript tests/test_federated.R
```

## Input: one long table of visits

One row per visit. Tab- or comma-separated text, or `.rds`. Column names are free; you name them in the config.

| column | required | notes |
|---|---|---|
| id | yes | person id; never written to any output |
| family id | no | if given, family-clustered standard errors are added |
| sex | yes | any coding; list the male and female values in the config |
| age at visit | yes | years (or months / days: set `age_unit`) |
| height (cm) | per analysis | |
| weight (kg) or BMI | yes | BMI is computed from weight and height when no BMI column is given |
| head circumference (cm) | per analysis | |
| genotype group | yes | one of the four CNVs or control; `genotype_map` in the config translates your labels; other rows are left out |
| autism diagnosis | no | yes / no per person |
| SRS total, SCQ total | no | raw totals, on any one row of the person; controls need the same questionnaire |
| age at questionnaire | no | if missing, the age of the row that holds the score is used |
| cohort / clinic | no | added as a covariate where the paper adds cohort (longitudinal models: within one genotype only; pooled models are adjusted for genotype, not cohort) |

Repeated visits are welcome: they feed the longitudinal models. Controls should not carry any of the four CNVs.

## What to send back

The CSV files in `out_dir`, and nothing else. Please open them before sending.

| file | content |
|---|---|
| `sample_n.csv` | persons by genotype x age band x sex, for height and for head circumference |
| `estimates.csv` | every model estimate with its standard error (schema below) |
| `tests.csv` | F tests and likelihood-ratio tests (heterogeneity, BMI x window, genotype by band) |
| `autism_by_bmi_band.csv` | share with an autism diagnosis by BMI band and by head-circumference band |
| `trajectory_bins.csv` | mean Z by fixed age bin, genotype and BMI tail |
| `trajectory_spline.csv`, `trajectory_spline_matrices.csv` | the curve per genotype and BMI tail: coefficients, covariance, n, supported ages |
| `scoring_info.csv` | how adults were scored; how many visits the cleaning rules removed |
| `run_info.csv` | settings, version, counts of excluded rows |

## Data rules

- **Nothing individual-level leaves the site.** The outputs hold no ids, no dates and no single measurements.
- **Minimum cell size** (`min_cell`, default 5). A count from 1 to 4 is written `<5` and its estimate is left empty.
- **No recovery by subtraction.** When a hidden count could be worked out from a total, a second count in the same
  line is written `>=5`. Its estimate is kept. The pooled "All four CNVs" count is masked when a genotype is hidden.
- **Models need more than the minimum.** A band model needs 15 carriers and 30 controls; a longitudinal model needs
  15 children; a curve needs 20 persons in the BMI tail. Rows that fall short say so in `status`.
- **The synthetic test data are the only data in this repository.** Never add real input or output files to it.

## Output schema: `estimates.csv`

One row per estimate. `site`, `script_version`, `run_date`, `scoring` and `estimator` are on every row of every table.

| column | meaning |
|---|---|
| analysis | cross-sectional band; cross-sectional all ages; sex; longitudinal; longitudinal sex (exploratory); autism diagnosis; autistic traits |
| outcome | Height, Head circumference, BMI, BMI Z, Head-circumference Z, SRS, SCQ |
| group | one genotype, `All four CNVs` (pooled), or `No CNV` (controls) |
| stratum | age band (2-9, 10-18, 19-39, 40-73), `all ages`, or longitudinal window |
| sex | all, male, female, or `male - female` for the three-way term |
| term | what is estimated, for example `CNV x BMI`, `BMI slope, carriers`, `BMI Z at the visit` |
| model | the model in words |
| estimate, se, df, p | estimate in SD units (log odds for diagnosis and for the ordinal check), standard error, degrees of freedom, P |
| se_family_clustered, p_family_clustered, n_families | sandwich standard error clustered on family (band models, when a family id is given) |
| n_carriers, n_controls, n_pairs, n_asd | sample sizes; `<5` suppressed, `>=5` masked |
| status | ok, or why there is no estimate |
| same_model_as_paper | yes, or what is new |

## The analyses

All follow the paper's scripts in the root of this repository. One analysis per genotype, then all four CNVs pooled
with a genotype term, then a heterogeneity test.

- **Cross-sectional, by age band** (`17_main_with_ukbb.R`, `82_growth_blocks_v2.R`). One visit per person. Within each
  band: `Z ~ CNV * BMI Z + sex + age`, carriers of one genotype plus controls. Reported: BMI slope in controls, BMI
  slope in carriers, and the difference (`CNV x BMI`). Pooled: the same with a genotype term. Tests: BMI slope differs
  across genotypes (carriers only), and BMI effect differs by genotype (controls plus genotypes).
- **Cross-sectional, all ages.** `Z ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z` (mgcv `bam`, REML).
  The BMI slope may change with age; the `CNV x BMI` term is the carriers' constant shift from the controls' slope.
- **Sex.** `CNV x BMI` within each sex, and the `CNV x BMI x sex` term, per band and for all ages.
- **Longitudinal** (`80_longitudinal_v2.R`). Children under 18. All cleaned visits.
  `change in Z to the next visit ~ Z at the visit + interval + sex (+ cohort | + genotype) + ns(age, 3) + BMI Z + (1 | child)`:
  cohort in the models of one genotype, genotype (and not cohort) in the pooled model (changed 2026-11-21, because cohort and genotype are confounded).
  Height pairs are at least 0.5 y apart, windows 2-8, 8-13, 13-18 y. Head-circumference pairs are at least 0.1 y apart,
  windows under 2, 2-5, 5-18 y. Controls with repeated visits get the same model (new; the paper has none).
- **Autism diagnosis** (`85_v2_asd_by_cnv.R`), if given. Height, head circumference and BMI Z in carriers with and without
  a diagnosis; logistic regression of diagnosis on BMI Z and on head-circumference Z; share with a diagnosis by band.
- **Autistic traits** (`51b_cog_panels_boxcox_combined.R`), if SRS or SCQ are given for carriers and controls. BMI slope
  per group, `CNV x BMI` against controls, and an ordinal check that does not depend on the scale.

## Known choices

- **One visit per person** in the cross-sectional models, chosen by `00_select_visit.R`, rule `fill_gaps` (it fills thin
  age bands; fixed seed). Set `"visit_rule": "latest"` for the latest visit.
- **Children** are scored on growth references: WHO 2006 / 2007 for height and BMI, UK 1990 for head circumference
  (`00_scoring.R`, `growth_refs_lms.csv`).
- **Adults** (height and BMI from 19 y, head circumference from 18 y) have no growth reference. They are scored against
  the site's own non-carrier adults of the same sex: Box-Cox transform, then the residual of `y ~ age + age^2`, divided
  by the residual SD. The Box-Cox exponents are **fixed** at the UK Biobank values (height 0.65 men, 0.60 women; BMI
  -0.75 men, -1.05 women; head circumference not transformed). A site needs 100 non-carrier adults per sex
  (`min_adult_ref`); otherwise its adults are left out and the adult bands are empty. That is expected at a paediatric site.
- **SRS and SCQ** raw totals (+ 1) are Box-Cox transformed with fixed exponents (SRS child form 0.32, SRS adult form 0.29,
  SCQ -0.16) and expressed as Z against the site's controls with the same form. The BMI used is the one measured nearest
  the questionnaire, within 3 years. The SRS form is taken from age (adult form from 18 y).
- **BMI tails** are fixed cut-offs on the BMI Z scale: at or above the 70th centile (Z >= 0.524) and at or below the
  30th (Z <= -0.524).
- **Meta-analysis.** Fixed effects, inverse-variance weights, Cochran Q and I2. Sites are strata. A suppressed or
  not-estimated cell is never pooled. Heterogeneity across the four genotypes is Cochran Q on the pooled estimates.
- **Traits, not risks.** Estimates describe how a trait changes with BMI.

## How the federated trajectory works

The paper's Figure 1a-e and 2a-e show height or head circumference against age, one curve for people with high BMI
and one for low BMI. The published figures draw every person as a point. A federated figure cannot do that.

So each site sends two summaries per genotype and BMI tail:

1. **Binned means.** Ages are cut into the same 12 bins at every site (2-4, 4-6, 6-8, 8-10, 10-13, 13-16, 16-19, 19-30,
   30-40, 40-50, 50-60, 60-73 years). The site sends the number of people, the mean Z and its standard error per bin.
2. **A curve on a shared ruler.** Every site fits the same curve shape: a natural cubic spline of log(age) with knots at
   6, 11, 16 and 30 years and ends at 2 and 73 years. The knots are fixed in `00_federated_common.R` and written into the
   output. The fit is plain least squares, with no smoothing choice. The site sends the six coefficients, their
   covariance matrix, the number of people and the age range that its data cover.

Because every site uses the same six building blocks, the coefficient sets can be averaged. `meta_combine.R` weights each
site by how precisely it knows each coefficient (generalised least squares with each site's covariance). A paediatric
site knows the childhood part well and the adult part not at all; it then counts for childhood and adds nothing for
adults. The pooled curve gets a 95% band and is drawn only at ages that at least one site covers. Binned means are
pooled by weighting each site's mean by its number of people.

This is exact, not approximate: with one common residual variance (`--trajectory_weights common_sigma`) the pooled curve
equals the curve fitted to all sites' data stacked together, to 13 decimal places in the test. The default lets each
site have its own residual variance, which equals weighted least squares on the stacked data.

Two limits. The curves are not smoothed, so a tail with few people gives a wavy curve with a wide band. And a site's
own curve is drawn only between its 5th-youngest and 5th-oldest person, rounded outward to half a year.

## Files

| file | purpose |
|---|---|
| `site_export.R` | all analyses at one site; writes the aggregate tables |
| `meta_combine.R` | combines the sites |
| `plot_federated.R` | trajectory panels, band forest, longitudinal forest |
| `00_federated_common.R` | shared constants: genotype labels, age bands, age bins, spline basis. Do not edit at a site |
| `00_scoring.R`, `growth_refs_lms.csv`, `00_select_visit.R` | copied unchanged from the repository root |
| `config_example.json` | config template |
| `tests/make_synthetic.R`, `tests/test_federated.R` | two simulated sites with known effects; end-to-end test |

## Test

`tests/test_federated.R` simulates two sites (ids start with `FAKE_`): one with children and adults, one paediatric
clinic with age-matched controls and few adults. It runs export, combine and plot and checks that the simulated effects
come back within sampling error, that the pooled curve equals the stacked-data curve, that small cells are suppressed
and cannot be recovered by subtraction, and that no output contains an id. Result on 2026-11-20: 36 of 36 checks pass.

## Validation against the paper run (lab clinical cohorts, 2026-11-20)

`site_export.R` was run on the lab's clinical table. UK Biobank adults were not available for this check.

- Age bands 2-9 and 10-18, height and head circumference, every genotype and pooled: identical to the paper run
  (estimates, standard errors, family-clustered standard errors, sample sizes, F tests).
- Longitudinal height, and longitudinal head circumference on the same visit sources: identical.
- Head circumference at 19 y and over, and all ages: differences up to 0.003 (bands) and 0.006 (all-ages sex term).
  The paper run estimated the adult BMI exponent in the clinical adults for this table; the package fixes it at the
  UK Biobank value. With `"adult_boxcox": "estimate"` the numbers are identical.
- Height at 19-39 y: differences up to 0.0004. The paper run has one UK Biobank control under 40.
- Height at 40-73 y and all ages: not comparable. The paper run adds UK Biobank adults.
