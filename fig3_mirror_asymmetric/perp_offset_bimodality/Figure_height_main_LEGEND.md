# Figure — legend and scoring methods

Figure file: `perp_offset_bimodality/figures/Figure_height_main.pdf` (and `.png`)
Built by:    `perp_offset_bimodality/code/make_height_figure.py`
Scores by:   `perp_offset_bimodality/code/score_models_on_criteria.py`
Model fits:  `CNV_height_BMI_corrected_directional_models_v2/code/run_axis_mixture_models.py`

---

## Figure legend

**Height: asymmetry is a rare single hit acting on the duplication axis.**

All panels use the 33 reciprocal CNV loci with published height effects. Blue, the 21
*mirror* loci; red, the 12 *asymmetric* loci. Effects are in trait SD units; genome-wide
burden is in cm per Mb.

**(A) The perpendicular offset from the mirror line is bimodal.** For each locus the
offset is its signed distance from the fitted mirror line DUP = −ρ·DEL, computed as
(ρ·β_DEL + β_DUP)/√(1+ρ²) with ρ = 0.351. Histogram, observed offsets; dashed grey, the
best single Gaussian; solid blue, the best two-component Gaussian mixture with a common
variance. Rug marks below the axis give each locus's published class. The mixture places
12 loci in the displaced component (mean −0.30) and 21 in the undisplaced component
(mean −0.02); that partition agrees with the published class labels for 33 of 33 loci.

**(B) The bimodality is stronger than a unimodal process produces.** Grey, the null
distribution of the 1-vs-2 component likelihood ratio over 5,000 parametric bootstrap
replicates drawn from a single Gaussian with per-locus measurement error added (null B);
blue line, the observed ratio of 13.20. *p* = 0.028. *Caveat:* the mirror line is
calibrated on the mirror-class loci only, so this test is partly circular — fitted on
all 33 loci the same statistic gives *p* = 0.10–0.29. It is reported as supporting, not
primary, evidence.

**(C–F) Mirror and asymmetric CNVs differ by size and by duplication effect, not by
deletion effect.** Each locus is one point; the horizontal bar is the group median.
*p* values are Welch two-sample *t* tests, on the log scale for length and gene count
and on the raw scale for effect sizes. Asymmetric loci are longer (median 1.19 vs
0.46 Mb, *p* = 0.006) and gene-richer (29.5 vs 19.0 genes, *p* = 0.019). Their deletion
effects are indistinguishable from mirror loci (*p* = 0.651), whereas their duplication
effects are strongly shifted downward (median −0.30 vs +0.02 SD, *p* = 1 × 10⁻⁶). The
asymmetry is therefore a duplication-side phenomenon, not a general increase in
severity.

**(G, H) Effect size scales with gene count only after stratifying by class.** Points are
loci; solid lines are Theil–Sen fits within class (heavy where *p* < 0.05, dashed where
not); the black dotted line is the pooled fit across all 33 loci. Spearman ρ and *p* are
inset. Pooled across classes the duplication effect shows no relationship with gene count
(ρ = +0.01, *p* = 0.968); within mirror loci it rises clearly (ρ = +0.51, *p* = 0.018).
The deletion effect behaves the same way (pooled ρ = −0.33, *p* = 0.064; mirror
ρ = −0.61, *p* = 0.003). Mixing the two classes masks a real dose response — a Simpson's
paradox. Within asymmetric loci neither slope is detectable (*p* = 0.713, 0.966), which is
what a single-hit mechanism predicts, though with *n* = 12 this is weak evidence: |ρ|
would have to exceed 0.59 to reach significance.

**(I–M) Only mirror + single hit satisfies all five criteria.** The criterion numbers
1–5 enumerate the five tests and are independent of the panel order: criterion 3,
bimodality, is the test shown in panels A and B; criteria 1 and 2 are the comparisons in
panels C–F; criterion 5 is the within-class slope in panel H. Rows are models; all three share
the same fitted mirror layer, and the two directional models share the same angular
freedom, so occupancy — continuous versus a rare single hit — is the only thing that
differs between them. Black vertical line, the observed value. Blue bar, the central 95%
of 2,000 simulated replicates; filled circle, the simulated median. The number above each
bar is the percentile of the observed value within that model's simulations; blue when it
falls inside the central 90%, red (open circle) when it falls outside. **(L)** is
different in kind: it shows a deterministic propagation to genome-wide burden rather than
a simulation. Grey band, the published 95% CI for the duplication burden (−0.134 to
−0.023 cm/Mb); black line, the published point estimate (−0.078 cm/Mb; Mac&eacute; et al.
2017). The deletion burden is the calibration anchor (−0.413 cm/Mb), so the duplication
burden is an out-of-sample prediction. Mirror alone predicts +0.144 cm/Mb — duplications
would make people *taller* — and cannot do otherwise: its per-gene duplication effect is
−ρ times its deletion effect, so the two burdens always carry opposite signs and a
U-shaped burden is structurally impossible. Both directional models reproduce the
U shape, and both land inside the published CI (continuous −0.087, single hit −0.074).

---

## The three models

Every model is conditional on the same mirror layer, fitted first and held fixed
(p = 0.528, m = 0.0527, ρ = 0.351, τ_DEL = 0, τ_DUP = 0.0507). Under the mirror layer a
locus with G genes has expected effects (G·δ, −ρ·G·δ): deletions and duplications push
the trait in opposite directions and both grow with gene count.

| | added parameters | what it adds | ΔAIC | ΔLOO |
|---|---|---|---|---|
| Mirror alone | 0 | nothing | 77.2 | 38.8 |
| Mirror + continuous | 3 | *every* gene carries a small directional effect, so the displacement is G·q·**b** and grows with gene count | 2.6 | 0.8 |
| Mirror + single hit | 3 | each gene independently carries the effect with rare probability q; one hit displaces the whole locus by a fixed vector **b**, and most loci are not hit at all | 0 | 0 |

Fitted values:

* **continuous** q = 0.076, **b** = (−0.042, +0.073), |**b**| = 0.085, θ = 119.8°
* **single hit** q = 0.019, **b** = (−0.195, +0.348), |**b**| = 0.399, θ = 119.3°

The two independently converge on essentially the same direction — θ ≈ 119°, i.e. the
duplication component is large and negative while the deletion component is *positive*,
pulling the deletion effect back toward null rather than making it worse. (This is the
geometry referred to elsewhere as "buffering": triplosensitive, deletion-tolerant. The
angle is a free parameter in both models, bounded only to (5°, 175°).) What separates the
two models is occupancy, not direction: the single-hit model affects 36% of loci on
average (9% at G = 5, 36% at G = 23, 67% at G = 58) against an observed asymmetric
fraction of 12/33 = 36%, while the continuous model touches 76% of loci with a much
smaller displacement each.

> **A note on naming.** Earlier drafts called one geometry "sparse single hit" and the
> other "buffering", which read as two different mechanisms. They are not. Both are
> single-hit mixtures with identical occupancy machinery; they differed only in whether
> the shift angle was capped at 85°. The capped version is a constrained special case
> that is not trying to reproduce the data, so it has been dropped from the main
> comparison, and the models are now named for the two things that actually vary:
> occupancy (*single hit* vs *continuous shift*) and direction (*constrained angle* vs
> *free angle*).

---

## How the scoring criteria work

### The simulation (criteria 1, 2, 3, 5)

For each model we generate 2,000 replicate datasets. Each replicate has the **real** gene
counts and the **real** per-locus standard errors of the 33 loci, so sampling noise and
the gene-count distribution are held fixed and only the model's own structure varies:

1. Compute the mirror mean vector and covariance for every locus from the fixed mirror
   parameters, adding each locus's own measurement variance.
2. Apply the model's directional component. For the continuous model this shifts the mean
   by G·q·**b** and inflates the covariance by G(2q − q²)·**bb**ᵀ. For the single-hit
   model, each locus is independently "hit" with probability 1 − (1 − q)^G and, if hit,
   its mean is displaced by **b**; the covariance is unchanged.
3. Draw one bivariate normal (β_DEL, β_DUP) per locus.

Each replicate is then analysed **exactly as the real data are**, with no access to the
simulation's hidden state:

* Offsets are computed from the simulated effects using the same ρ.
* The two classes are **re-inferred** from those offsets by the same equal-variance
  2-component mixture used in panel A, with a hard assignment at the midpoint between the
  component means.
* The five statistics below are computed from the re-inferred classes.

This blind re-inference is the point of the design. A continuous model has no classes by
construction, but noise plus a gene-count gradient will still make a mixture fit find two
groups; scoring it this way means it is never penalised by fiat for "not having classes".
Replicates where the inferred split is degenerate (fewer than 2 loci on either side) are
dropped.

The observed values are computed the same way — from mixture-inferred classes, not from
the published labels. For the real height data the two coincide exactly (33/33), so the
observed numbers are directly comparable to the simulated ones.

### The five criteria

| # | panel | statistic | observed |
|---|---|---|---|
| 1 | I | size difference: mean log(1+genes) in the displaced class minus the majority class | **+0.549** |
| 2 | J | duplication-effect separation: mean β_DUP, displaced minus majority | **−0.321** |
| 3 | K | bimodality: likelihood ratio, 2-component vs 1-component fit to the offsets | **13.20** |
| 4 | L | U-shaped burden: predicted genome-wide duplication burden per Mb, with the deletion burden calibrated | **−0.078** (published) |
| 5 | M | dose response: inverse-variance weighted slope of β_DUP on gene count, **within the majority class only** | **+0.0041** |

Criterion 2 is also scored on the deletion axis, where the observed separation is near
zero (+0.062) — the "not DEL" half of section 1. It is not given its own panel; the values
are in the results table below.

### Reading the percentile

The number above each bar is the fraction of that model's simulated replicates falling at
or below the observed value. A value near 0.5 means the model routinely produces what we
see. Values near 0 or 1 mean the observed data sit in a tail of what the model produces.
Blue marks 0.05–0.95; red marks outside that range.

The percentile is a posterior-predictive check, not a hypothesis test, and no multiplicity
correction is applied across the five criteria.

### Why criterion 4 is handled differently

The burden criterion is not a property of the 33 loci and cannot be simulated from them.
It is an algebraic propagation of the fitted per-gene effects to genome-wide CNV burden:

1. From the fitted model, compute the expected per-affected-gene deletion effect and
   duplication effect, averaged over the observed gene-count distribution.
2. Choose the effective dosage-sensitive gene density per Mb that makes the predicted
   **deletion** burden equal the published value, −0.413 cm/Mb. This is the single
   calibration constant and it fixes the scale.
3. The predicted **duplication** burden then follows with no further freedom. It is
   compared against the published −0.078 cm/Mb (95% CI −0.134 to −0.023).

A U-shaped burden — both deletion and duplication burden reducing height — requires both
predictions to be negative. Because the mirror layer sets the per-gene duplication effect
to exactly −ρ times the per-gene deletion effect, mirror alone always predicts burdens of
opposite sign. No parameter value can rescue it. Both directional models predict negative
duplication burdens inside the published CI.

Reference values are in `CNV_height_BMI_corrected_directional_models_v2/code/mace_burden_reference.csv`,
transcribed from Mac&eacute; et al. 2017 Supplementary Table 1 (confidence intervals derived
from the reported β and two-sided *p* by Wald inversion).

---

## Results table

Simulated median [2.5%, 97.5%] with the percentile of the observed value in parentheses.

| criterion | observed | Mirror alone | Mirror + continuous | Mirror + single hit |
|---|---|---|---|---|
| 1. size difference (log genes) | +0.549 | +0.009 [−0.83, +0.81] (0.90) | +0.612 [−0.33, +1.03] (0.40) | +0.477 [−0.38, +0.90] (0.62) |
| 2. DUP separation (SD) | −0.321 | −0.190 [−0.44, −0.03] (0.12) | −0.298 [−0.61, −0.16] (0.39) | −0.312 [−0.52, −0.18] (0.44) |
| 2b. DEL separation (SD) | +0.062 | −0.215 [−0.79, +0.21] (0.91) | −0.127 [−0.60, +0.24] (0.87) | −0.058 [−0.59, +0.25] (0.79) |
| 3. bimodality (LR) | 13.20 | 4.48 [0.14, 20.1] (0.91) | 2.92 [0.08, 14.6] (0.96) | 4.12 [0.20, 16.3] (0.94) |
| 4. DUP burden (cm/Mb) | −0.078 [−0.134, −0.023] | **+0.144** (wrong sign) | −0.087 ✓ | −0.074 ✓ |
| 5. within-class DUP slope | +0.0041 | +0.0043 [−0.0004, +0.0109] (0.47) | +0.0015 [−0.0048, +0.0091] (0.78) | +0.0033 [−0.0030, +0.0089] (0.61) |

### Pass / fail

A criterion is **passed** when the observed value falls inside the central 90% of that
model's simulations (criteria 1, 2, 2b, 3, 5) or inside the published 95% CI
(criterion 4).

| | 1. size | 2. DUP sep | 2b. DEL sep | 3. bimodality | 4. U-shaped burden | 5. dose response | total |
|---|---|---|---|---|---|---|---|
| Mirror alone | ✓ (0.90) | ✓ (0.12) | ✓ (0.91) | ✓ (0.91) | **✗ wrong sign** | ✓ (0.47) | 5/6 |
| Mirror + continuous | ✓ (0.40) | ✓ (0.39) | ✓ (0.87) | **✗ (0.96)** | ✓ | ✓ (0.78) | 5/6 |
| **Mirror + single hit** | **✓ (0.62)** | **✓ (0.44)** | **✓ (0.79)** | **✓ (0.94)** | **✓** | **✓ (0.61)** | **6/6** |

**Mirror + single hit is the only model that satisfies every criterion.** Each of the
other two fails exactly one, and they fail different ones:

* **Mirror alone fails criterion 4,** and fails it structurally rather than marginally.
  Its per-gene duplication effect is fixed at −ρ times its deletion effect, so the
  deletion and duplication burdens always carry opposite signs. No parameter value
  produces a U shape. This is the result that makes the published U-shaped burden
  (−0.132 cm/Mb, *p* = 4.5 × 10⁻⁷) incompatible with a purely mirror architecture.
* **Mirror + continuous fails criterion 3.** Spreading the directional effect across every
  gene produces a smooth gradient rather than two groups, so its offsets are the *least*
  bimodal of the three (median LR 2.9 against an observed 13.2).

### How strong is each of those conclusions

* **Criterion 4 against mirror alone is decisive**, because it is algebraic. It is also
  the only criterion that is.
* **Criterion 4 does not separate continuous from single hit.** Both reproduce the U shape
  and both land inside the published CI. The U shape argues for *some* directional
  component, not specifically for a sparse one.
* **Criteria 1, 2, 2b and 5 favour the single-hit model directionally, not decisively.**
  With 33 loci the simulated intervals are wide, and on these four criteria the observed
  value lies inside the 95% interval for all three models. The evidence is in where within
  those intervals it sits: mirror alone is near the edge on 1 (0.90) and 2b (0.91), the
  continuous model is off-centre on 5 (0.78) because spreading the effect over every gene
  dilutes the within-class dose response, and the single-hit model is the most centred on
  all four.
* **Criterion 3 discriminates weakly and should not be leaned on.** The observed
  likelihood ratio of 13.20 sits in the upper tail for all three models (0.91, 0.96,
  0.94); the continuous model's failure is a matter of 0.96 versus 0.94. The offsets are
  more sharply separated than even the single-hit model comfortably generates, which may
  reflect the partial circularity noted for panel B, heterogeneity none of the models
  capture, or chance at *n* = 33.
* **Model comparison on the likelihood itself** (ΔAIC 77.2 / 2.6 / 0, ΔLOO 38.8 / 0.8 / 0)
  gives the same ordering, with the single-hit model preferred over the continuous one by
  2.6 AIC at equal parameter count.

**Summary.** The genome-wide burden criterion rules out mirror alone on structural
grounds, and the single-hit model is the only one of the three that satisfies all five
criteria. The margin over the continuous model rests on the convergent but individually
weak evidence of criteria 1, 2, 3 and 5 together with the likelihood comparison, and
should be presented at that strength.
