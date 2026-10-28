Corrected directional models package (revised BMI simulations)
===============================================================

Height
------
Single-locus and burden conclusions are unchanged from the corrected directional analysis. Mirror alone cannot explain the lower-left displaced height loci or the adverse direction of DUP burden. Both continuous and sparse asymmetric models correct the burden direction, but the sparse single-hit model remains the most parsimonious because it also matches the occupancy, size, and class-structure observations among recurrent CNVs.

BMI
---
Single-locus and burden model-comparison results are unchanged: the BMI data are best described by an obesity-directed mirror relationship plus an upper-right directional shift, with the continuous model preferred over sparse single hits by conditional AIC and leave-one-out prediction, while BIC favors the simpler mirror-only model because the extra shift is small.

What changed
------------
The original BMI *simulation clouds* used the literal mirror point estimate p≈0.02 and m≈0.28 SD as a gene-occupancy model, which produced an artificial vertical streak of points at DEL≈0. The revised simulations now generate the mirror component from a moment-matched aggregate Gaussian approximation that preserves the fitted CNV-level means, variances, and DEL–DUP covariance as a function of gene count without creating the zero-inflated artifact. This affects the visual simulations and the figure-generation code; it does not change the inferential model-comparison tables or burden predictions.

Updated files
-------------
- behavior_figure_revised/: full revised figure package
- BMI_SIMULATION_UPDATE.txt: explanation of the revision
