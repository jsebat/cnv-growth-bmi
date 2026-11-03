CNV burden behavior figure
==========================
This package contains a main figure summarizing how the fitted mirror, continuous, and sparse single-hit models propagate from recurrent CNV effects to genome-wide CNV burden predictions.

Files
-----
- figures/Figure_burden_model_behavior.pdf|png : rendered figure
- code/make_burden_behavior_figure.R : R script to reproduce the figure
- inputs/*.csv : burden-prediction and literature-target tables
- outputs/simulated_burden_clouds.csv : Python reference simulation cloud used for rendering here

Interpretation
--------------
The observed panels summarize the literature burden targets: exact deletion-burden coefficient from Macé et al. and the overall adverse burden direction from Auwerx et al. The model panels show how varying the effective dosage-sensitive gene density around the Macé calibration point generates predicted DEL and DUP burden coefficients. Mirror alone predicts the wrong DUP burden direction for both height and BMI; both asymmetric models correct that direction.
