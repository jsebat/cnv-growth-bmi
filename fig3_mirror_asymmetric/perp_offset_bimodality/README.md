# Is the perpendicular offset from the mirror line bimodal?

2026-10-10. Height and BMI, exact Sacks reciprocal-CNV effects
(`CNV_height_BMI_corrected_directional_models_v2/{height,bmi}/inputs/`).

## Question

The mirror model places every locus on the line `DUP = -rho*DEL`, with position along the
line set by gene count. Each locus's signed **perpendicular distance** from that line is
the part the mirror model does not explain. If asymmetry is a discrete class, those
distances should be a two-component mixture; if it is a continuous process, a single
Gaussian.

## Method

1. Fit the mirror model (`fit_mirror`, ported exactly from the published package) and take
   its `rho`.
2. For each locus, `offset = (rho*DEL + DUP) / sqrt(1 + rho^2)`, with
   `SE = sqrt(rho^2*seDEL^2 + seDUP^2) / sqrt(1 + rho^2)`.
3. Fit a one-component and a two-component Gaussian to the offsets; statistic is
   `LRT = 2*(logL2 - logL1)`.
4. Parametric bootstrap: simulate 5,000 datasets of the same size from the fitted single
   Gaussian, refit both models to each, and count how often the null reaches the observed
   LRT. Mixture LRTs do not follow chi-square, so the bootstrap **is** the test.

Two nulls:

| | |
|---|---|
| **A** | plain single Gaussian fitted to the observed offsets |
| **B** | the same latent Gaussian **plus each locus's own measurement error**, so the null carries the heteroscedastic SEs the loci actually have |

Two alternatives:

| | extra par. | |
|---|---|---|
| **equal variance** | 2 | the mechanistically faithful one: under a sparse single hit, a hit *displaces* a locus but does not change its dispersion, so the components differ in mean only |
| **free variance** | 3 | more flexible, but at n~33 it can pin a component on one outlier |

The observed data and every bootstrap replicate are fitted with the **identical** EM
procedure and the same number of random starts (25). This matters: an earlier R attempt
fitted the observed data effectively harder than the nulls and returned p = 0.0002, which
was an artefact of the null replicates being under-fitted, not a real signal.

## Mirror calibration

Reported both ways, because the choice changes the answer.

- **mirror-class only** (primary) — the line is defined by the 21 mirror-class height loci,
  and the 12 asymmetric loci are held out of that fit, so their offsets are measured
  against a line they played no part in shaping. `rho = 0.351`.
- **full data** — all 33 loci. The asymmetric loci tilt the very line their displacement is
  measured against: `rho` flattens to `0.218`.

The BMI locus table carries no `Dose_response_class`, so all 32 loci are used and the two
calibrations coincide.

## Result

5,000 bootstrap replicates, 25 EM starts (`outputs/bimodality_test.csv`).

| trait | mirror fit | alternative | LRT | p (null A) | p (null B) |
|---|---|---|---|---|---|
| **Height** | **mirror-class** | **equal var** | **13.20** | **0.0034** | **0.0278** |
| **Height** | **mirror-class** | **free var** | **13.20** | **0.0238** | **0.0904** |
| Height | full data | equal var | 5.22 | 0.1030 | 0.1580 |
| Height | full data | free var | 7.59 | 0.1958 | 0.2865 |
| BMI | (either) | equal var | 3.88 | 0.1908 | 0.4511 |
| BMI | (either) | free var | 9.40 | 0.1034 | 0.3495 |

**Height is bimodal, BMI is not** — with the mirror calibrated on mirror-class loci.
The height mixture puts ~37% of loci in a component displaced ~0.28 SD off the line, which
is the same perpendicular displacement the sparse single-hit model estimates
(`PGC_asym_model_diff` reports perp ~0.27-0.31 SD at P(hit|G=23) ~0.35).

### The two modes are the published classes

Under the mirror-class calibration the separation is **complete**:

| class | n | mean offset | min | max |
|---|---|---|---|---|
| asymmetric | 12 | -0.306 | -0.425 | -0.224 |
| mirror | 21 | -0.023 | -0.145 | +0.135 |

All 12 asymmetric loci fall below -0.18; all 21 mirror loci fall above it. Not one
misclassification.

That is the strongest and the weakest thing about this result. It confirms the published
`Dose_response_class` labels correspond to a real geometric displacement. But because the
line was fitted on the mirror-class loci and the labels themselves encode dose-response
behaviour, the bimodality is **partly a restatement of the labels**, not an independent
discovery of structure in unlabelled data. The label-free version of the test is the
full-data calibration, and that one is not significant (p = 0.10-0.20).

Both belong in any write-up: the mirror-class number says the classes are geometrically
real; the full-data number says 33 loci are not enough to recover the classes without
being told them.

Two further caveats:

- **It depends on the mirror calibration.** Fitting the mirror to all 33 loci drops height
  to p = 0.10-0.20. The result is conditional on holding the asymmetric loci out of the
  line fit, which is defensible but is an assumption, not a free lunch.
- **It weakens under null B.** Carrying the per-locus SEs into the null moves height from
  0.0034 to 0.0278 (equal var) and from 0.0238 to 0.0904 (free var). Part of the apparent
  separation is measurement error that the plain-Gaussian null does not model.

BMI's free-variance LRT (9.40) is *larger* than height's under the full-data calibration,
but it is not evidence of bimodality: that fit pins a component on a single locus
(weight 0.031, component SD exactly at the 5%-of-sd floor). The equal-variance constraint
removes the artefact and BMI drops to 3.88. This is why equal variance is the primary
alternative.

## Files

```
code/perp_offset_bimodality.py   mirror fit, offsets, mixture fits, parametric bootstrap
code/make_figure.py              the bimodality figure on its own
outputs/bimodality_test.csv      one row per trait x mirror-fit x alternative
outputs/perp_offsets_per_locus.csv   per-locus offsets and their SEs
outputs/null_LRT_draws.npz       the 5,000 null LRT draws per configuration
figures/Figure_perp_offset_bimodality.{pdf,png}

code/score_models_on_criteria.py scores three models against the five criteria
code/make_height_figure.py       the height manuscript figure (panels A-M)
outputs/model_criteria_scores.csv    simulated medians, intervals and percentiles
outputs/model_criteria_observed.csv  the observed value of each criterion
outputs/model_criteria_draws.npz     the per-replicate draws
figures/Figure_height_main.{pdf,png}
Figure_height_main_LEGEND.md     figure legend + how the scoring criteria work
```

Reproduce:

```bash
python code/perp_offset_bimodality.py --nboot 5000 --starts 25
python code/make_figure.py
python code/score_models_on_criteria.py --nsim 2000
python code/make_height_figure.py
```

`score_models_on_criteria.py` reads the model fits written by
`../CNV_height_BMI_corrected_directional_models_v2/code/run_axis_mixture_models.py`, so
run that first if the fits have changed.

## Relation to the other folders

- `PGC_asym_model_diff/` decomposes the *fitted model's* shift vector into perpendicular
  and along-line parts and shows the clustering in simulated clouds. That is a property of
  the fitted models.
- This folder tests the *observed* offsets directly, with no asymmetry model fitted.
