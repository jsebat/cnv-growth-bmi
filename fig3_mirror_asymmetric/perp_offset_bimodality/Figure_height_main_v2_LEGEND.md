# Figure (v2) — legend

Figure file: `perp_offset_bimodality/figures/Figure_height_main_v2_me.pdf` (and `.png`)
Built by:    `perp_offset_bimodality/code/make_height_figure_v2.py --bimodality-dir me/all_loci --criteria-dir me/all_loci --name Figure_height_main_v2_me --note "..."`
Bimodality:  `perp_offset_bimodality/code/me_bimodality.py` → `outputs/me/all_loci/`
Simulations: `perp_offset_bimodality/code/score_models_on_criteria.py --me --out-dir outputs/me/all_loci`
Panel J–K numbers: `outputs/me/all_loci/model_expected_bimodality_v2.csv`

Changes from the original figure (`Figure_height_main_LEGEND.md`): Arial throughout;
bimodality in A, B, J and K uses the measurement-error-aware mixture; the five-criteria
panels I–M are replaced by I–K, which show only the two places the models differ.
The descriptions of the three models, their fitted parameters, and the simulation design
are unchanged and remain in `Figure_height_main_LEGEND.md`.

---

## Figure legend

**Height: asymmetry is a rare single hit acting on the duplication axis.**

All panels use the 33 reciprocal CNV loci with published height effects. Blue, the 21
*mirror* loci; red, the 12 *asymmetric* loci. Effects are in trait SD units; genome-wide
burden is in cm per Mb.

**(A) The perpendicular offset from the mirror line is bimodal.** For each locus the
offset is its signed distance from the fitted mirror line DUP = −ρ·DEL, computed as
(ρ·β_DEL + β_DUP)/√(1+ρ²) with ρ = 0.351. Histogram, observed offsets. Both fitted
curves carry each locus's own measurement error: a locus's offset is modelled as its
group mean plus true scatter τ plus its own sampling SE, so noisy loci count for less.
Dashed grey, the best single group; solid blue, the best two-group mixture with a common
τ (each curve is the average of the 33 loci's predictive densities). Rug marks below the
axis give each locus's published class. The mixture places 34% of loci in a displaced
group (mean −0.27 SD) and the rest in an undisplaced group (mean −0.03 SD), with true
scatter τ = 0.037 SD. Assigning each locus to its more probable group agrees with the
published class labels for 32 of 33 loci. The exception is an asymmetric locus whose
offset (−0.29) is typical of its class but whose SE (0.24) is so large that the fit
cannot place it.

**(B) The bimodality is stronger than a unimodal process produces.** Grey, the null
distribution of the two-group vs one-group likelihood ratio over 5,000 parametric
bootstrap replicates drawn from the fitted single group, each locus with its own SE;
blue line, the observed ratio of 9.9. *p* = 0.016. *Caveat:* the mirror line is
calibrated on the mirror-class loci only, so this test is partly circular. With the line
fitted on all 33 loci, the plain (non-measurement-error) version of the test gives
*p* = 0.10–0.16; the measurement-error-aware version has not been run on that
calibration.

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

**(I–K) Where the models differ.** Three models are compared. All share the same fitted
mirror layer, and the two directional models share the same angular freedom, so
occupancy — every gene contributing a little (*continuous*) versus a rare per-gene hit
that displaces the whole locus (*single hit*) — is the only thing that differs between
them.

**(I) Only the directional models predict that duplication burden reduces height.**
Predicted change in height against genome-wide duplication burden. Each model is
calibrated so that its predicted deletion burden equals the published −0.413 cm/Mb, so
its duplication slope is an out-of-sample prediction. Black line and grey band, the
published duplication burden and its 95% CI (−0.078 cm/Mb, −0.134 to −0.023, *p* = 0.006;
Macé et al. 2017). Mirror alone predicts +0.144 cm/Mb — duplications would make people
*taller* — and cannot do otherwise: its per-gene duplication effect is −ρ times its
deletion effect, so the deletion and duplication burdens always carry opposite signs.
Both directional models predict a negative slope inside the published CI (continuous
−0.087, single hit −0.074 cm/Mb). This panel separates mirror alone from the other two
but not continuous from single hit.

**(J) Only the single-hit model expects two groups of loci.** The offset distribution
each model expects, pooled over its 2,000 simulated datasets (each with the real gene
counts and per-locus SEs of the 33 loci), drawn over the observed offsets. Mirror alone
expects one narrow peak at zero; the continuous model one broad, skewed peak; the
single-hit model two peaks, at about −0.27 and 0 SD, matching the observed groups.

**(K) Only the single-hit model typically produces bimodality as strong as observed.**
Grey, the unimodal null from B. Curves, the distribution of the same likelihood ratio
across each model's 2,000 simulated datasets, each analysed exactly as the real data.
Blue line, the observed ratio. Legend values are the median *p*, against the unimodal
null, of each model's simulated datasets — the result a typical dataset from that model
would give. Mirror alone (median *p* = 0.42) and continuous (0.31) produce likelihood
ratios close to the null; the single-hit model (0.014) produces bimodality of the
observed strength, and the observed ratio sits at the centre of its distribution.

---

## Numbers behind I–K

| model | DUP burden (cm/Mb) | median LR | median *p* vs unimodal null | datasets with *p* < 0.05 |
|---|---|---|---|---|
| Mirror alone | +0.144 | 2.3 | 0.42 | 4.8% |
| Mirror + continuous | −0.087 | 3.0 | 0.31 | 15% |
| Mirror + single hit | −0.074 | 10.2 | 0.014 | 73% |
| observed | −0.078 (Macé 2017) | 9.9 | 0.016 | |

Each simulated dataset's *p* is its likelihood ratio referred to the observed null in B,
rather than to a bootstrap null of its own. This is an approximation: each dataset's own
null would be fitted to that dataset.

The simulated datasets in J and K are regenerated in `make_height_figure_v2.py` from
`score_models_on_criteria.simulate` with its seed and model order, so they are the same
replicates whose likelihood ratios are stored in `model_criteria_draws.npz`. Checked by
refitting a subset of replicates per model: plain-test run r ≥ 0.997 (60 replicates),
measurement-error run r ≥ 0.999 (40 replicates).
