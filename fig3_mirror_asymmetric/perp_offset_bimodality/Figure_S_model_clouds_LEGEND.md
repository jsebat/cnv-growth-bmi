# Supplementary figure — continuous vs single-hit model clouds

Figure file: `perp_offset_bimodality/figures/Figure_S_model_clouds.pdf` (and `.png`)
Built by:    `perp_offset_bimodality/code/make_supp_model_clouds.py`
Typical replicates: `outputs/me/all_loci/supp_model_clouds_typical.csv`

---

## Figure legend

**Continuous versus single-hit asymmetry: simulated and observed height effects.**

Each panel plots duplication effect against deletion effect (trait SD units) for CNV
loci. Orange, the 33 observed loci: circles, the 21 mirror loci; triangles, the 12
asymmetric loci. Dashed purple, the fitted mirror line DUP = −ρ·DEL (ρ = 0.35). Both
models share this fitted mirror layer and differ only in how the additional, directional
effect is distributed across genes. Simulations use the real gene counts and the real
per-locus standard errors of the 33 loci, so the simulated points carry the same
sampling noise as the observed estimates.

**(A) Mirror + continuous.** Every gene carries a small directional effect (q = 0.076 of
genes), so a locus's displacement from the mirror line grows smoothly with its gene count
— about (+0.07, −0.13) SD for a 23-gene locus. Points, 150 simulated datasets pooled
(4,950 loci). The model produces a single diffuse cloud that broadens below the mirror
line; it contains no separate group of loci. ΔAIC 2.6, ΔLOO 0.8 relative to the
single-hit model.

**(B) Mirror + single hit.** Each gene independently carries the effect with low
probability (q = 0.019); a single hit displaces the whole locus by a fixed vector,
+0.20 SD on deletion and −0.35 SD on duplication, and most loci carry no hit (probability
of a hit 0.36 for a 23-gene locus). Points, 150 simulated datasets pooled; blue, loci
without a hit; red, loci with a hit. The model produces two parallel bands, the lower
displaced mainly on the duplication axis, and the observed asymmetric loci fall in the
lower band. ΔAIC 0.0, ΔLOO 0.0 (best-fitting model).

**(C, D) One typical dataset from each model.** A single simulated dataset of 33 loci
(open diamonds), shown over the observed loci. Each is the replicate whose bimodality
likelihood ratio (measurement-error-aware test, as in the main figure) is closest to that
model's median over 2,000 replicates, so it shows what one study of these loci would
typically see if the model were true. **(C)** Continuous model, likelihood ratio 3.0: the
loci spread around and below the mirror line with no distinct group. **(D)** Single-hit
model, likelihood ratio 10.2; 10 of 33 loci carry a hit (red). The hit loci form a
separate group below the line, as the observed asymmetric loci do. The observed
likelihood ratio is 9.9.

---

## Two-panel version (A, B only)

Figure file: `perp_offset_bimodality/figures/Figure_S_model_clouds_AB.pdf` (and `.png`),
written by the same script from the same simulations.

**Continuous versus single-hit asymmetry: simulated and observed height effects.**

Each panel plots duplication effect against deletion effect (trait SD units) for CNV
loci. Orange, the 33 observed loci: circles, the 21 mirror loci; triangles, the 12
asymmetric loci. Dashed purple, the fitted mirror line DUP = −ρ·DEL (ρ = 0.35). Both
models share this fitted mirror layer and differ only in how the additional, directional
effect is distributed across genes. Points are 150 simulated datasets per model, pooled
(4,950 loci), generated with the real gene counts and per-locus standard errors of the
33 loci, so the simulated points carry the same sampling noise as the observed estimates.

**(A) Mirror + continuous.** Every gene carries a small directional effect (q = 0.076 of
genes), so a locus's displacement from the mirror line grows smoothly with its gene count
— about (+0.07, −0.13) SD for a 23-gene locus. The model produces a single diffuse cloud
that broadens below the mirror line; it contains no separate group of loci. ΔAIC 2.6,
ΔLOO 0.8 relative to the single-hit model.

**(B) Mirror + single hit.** Each gene independently carries the effect with low
probability (q = 0.019); a single hit displaces the whole locus by a fixed vector,
+0.20 SD on deletion and −0.35 SD on duplication, and most loci carry no hit (probability
of a hit 0.36 for a 23-gene locus). Blue, loci without a hit; red, loci with a hit. The
model produces two parallel bands, the lower one displaced mainly on the duplication
axis, and the observed asymmetric loci fall in the lower band. ΔAIC 0.0, ΔLOO 0.0
(best-fitting model).

---

## Notes

* Model parameters are the fits in
  `CNV_height_BMI_corrected_directional_models_v2/height/outputs/height_axis_mixture_model_comparison.csv`
  ("Mirror + continuous shift, free angle" and "Mirror + single hit, free angle"), the same
  models compared in panels I–K of `Figure_height_main_v2_me`. ΔAIC and ΔLOO are the
  conditional values from that table (mirror layer fitted on the mirror-class loci and
  held fixed).
* The pooled clouds in A and B use a separate random seed from the 2,000 scoring
  replicates. The typical datasets in C and D are regenerated from the scoring run's seed
  and order; refitting a subset reproduces the stored likelihood ratios (r ≥ 0.999).
* Earlier versions of this comparison (`PGC_asym_model_diff/code/loo_and_clouds.R`) fitted
  every parameter jointly on a 34-locus table and showed latent clouds without sampling
  noise. Their AIC and LOO values are on a different scale and are not comparable with
  the values here.
