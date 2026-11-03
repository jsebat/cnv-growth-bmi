# Why Molly's v4 asymmetric model does not reproduce Jon's clustered single-hit shift (height)

2026-10-09. Inputs: Molly v4 `inputs/height_loci.csv` (34 loci, incl. 22q11.21_A-D).

## Method
Molly's CNV-hit likelihood (`cnv_hit_deleterious_model.R`) generalised with one toggle per
difference between her model and Jon's v2 single-hit model:

| toggle | Molly v4 | Jon v2 |
|---|---|---|
| shift direction | fixed along (1, rho) | free angle (v2 bounded 5-85 deg; hit the bound) |
| hit size | N(beta, tau^2), tau free | fixed vector (tau = 0) |
| mirror slope rho | fixed log2(3/2) = 0.585 | fitted (0.35 on mirror loci; ~0.31 jointly) |
| off-line residual tauU | none | fitted (0.05) |
| mirror fit | joint, all loci | on the 21 mirror-class loci, then fixed (conditional AIC) |

`code/factor_test.R`: all 32 combinations, hit and matched continuous model each ->
`outputs/height_factor_grid.csv`. Validation: row (joint, fixed, free, fixed, zero) reproduces
Molly (logL 13.00, pi 0.087, beta -0.191, tau 0.111); row (cond, free, zero, free, free)
reproduces Jon's mirror (rho 0.351, tauU 0.051) and q = 0.019.

`code/loo_and_clouds.R`: common-scale comparison (same 34 loci, everything fit jointly, all
parameters counted, leave-one-locus-out with full refit) -> `outputs/height_common_scale_comparison.csv`,
`outputs/height_loo_pointwise.csv`, `figures/Figure_height_factor_clouds.{png,pdf}`.

## Result
A minority-hit, parallel, perpendicular-offset cluster (P(hit | G = 23) ~0.35, shift ~0.27 SD
off the line) appears only when **all three** hold: rho fitted (~0.3), direction free, tau = 0.
Reverting any one gives P(hit) 0.77-0.88: the whole cloud shifts and no cluster forms.
- rho = 0.585 puts the large mirror loci (1q21.1, 3q29, 16p11.2 BP4-5) ~0.2 SD below the line,
  so they also need the shift -> nearly every CNV is "hit".
- direction (1, rho) moves DEL 1.7x as much as DUP; the asymmetric loci are DUP-only shifts.
- free tau lets one wide hit distribution absorb everything (joint fit then degenerates to rho = 0).
Conditional vs joint mirror fitting is not what drives the clustering (joint and conditional
Jon-type fits agree).

Common-scale fit (height): Molly continuous AIC -23.5 / LOO 10.9; Jon-type single hit (joint)
AIC -22.6 (dAIC 0.9) / LOO 4.9 (dLOO 6.0), of which 5.0 comes from the single locus 1q21.1.

## Follow-up 2026-10-10: two clusters in Molly's model WITHOUT tauU
Scripts: `code/her_model_two_clusters.R` (8 variants, pi profiled), `code/her_model_rho_profile.R`,
`code/her_model_theta_profile.R` (hit direction profiled); outputs `outputs/her_model_*`, figure
`figures/Figure_her_model_no_tauU_clouds.*`. Base = her continuous mirror + continuous deleterious (sigma_d),
plus a fixed-size CNV-level single hit; no tauU anywhere.
- Minimal change (her model exactly, rho 0.585, deleterious along (1,rho), + fixed-size hit in a free direction):
  logL 18.7, k 7, AIC -23.5 = her continuous (-23.5); P(hit|G=23) 0.44; shift (-0.08, -0.10).
- Deleterious direction freed: logL 20.9, AIC -23.8 (k 9) or -25.8 if rho is not counted. Rho is then not
  identified (two free G-scaled components give an arbitrary covariance); bands are parallel, slope ~ -0.75.
- Hit direction is not identified: logL 19.6-20.9 for every direction from DEL-left (-180 deg) to DUP-down
  (-90 deg). A pure DUP-down hit: logL 19.6, AIC -21.3, P(hit) 0.44-0.61, shift -0.16 DUP. Jon's
  geometry (minority hit, ~0.3 SD shift right-and-down, -60 deg): logL 17.1.
So tauU is not what makes a two-cluster model competitive: sigma_d does the same job. But no variant beats her
continuous model by more than ~2 AIC units, and the direction and frequency of the hit are not pinned down.

## Power simulation 2026-10-10 (`code/power_simulation.R`, outputs `outputs/power_simulation_*`)
1000 simulated 34-locus datasets per scenario (real gene counts and SEs); each fit with CONT (Molly continuous,
k 4) and HIT (CONT + fixed-size single hit, free direction, k 7), no tauU. dAIC = AIC_HIT - AIC_CONT; real data -0.02.
| truth | AIC picks HIT | by >2 | by >4 | median dAIC | P(|dAIC|<1) |
|---|---|---|---|---|---|
| A two clusters as fitted (P 0.44, shift 0.13 SD) | 0.75 | 0.59 | 0.44 | -3.1 | 0.16 |
| B Jon geometry (P 0.35, shift 0.3 SD, tauU) | 0.87 | 0.77 | 0.67 | -8.1 | 0.09 |
| C no clusters (Molly continuous) | 0.37 | 0.18 | 0.09 | +0.9 | 0.26 |
The observed tie is ~1.6x more likely with no clusters than with clusters as fitted, ~3x more likely than under
Jon's geometry: weak evidence against clustering, moderate against the large-shift version. AIC is liberal here
(18% false positives at dAIC < -2). Caveat: B's data include tauU, which neither analysis model has.
