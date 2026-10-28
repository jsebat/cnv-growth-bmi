CNV continuous-vs-single-hit behavior figure (revised BMI simulations)
=====================================================================
This package regenerates the main figure comparing continuous and sparse single-hit models for reciprocal CNV effects.

Revision note
-------------
The original BMI simulation panels showed an artificial vertical streak at DEL≈0 because the fitted mirror point estimate (very low contributor probability with a large per-gene effect) was used too literally as a generative gene-occupancy model. The revised version uses a moment-matched aggregate Gaussian approximation for the mirror component, preserving the fitted CNV-level means, variances, and DEL–DUP covariance while avoiding the zero-inflated artifact.

Files
-----
- code/make_model_behavior_figure.R : R script to recreate the figure
- code/python_reference_generator.py : matched Python reference implementation used to render the included figure here
- figures/Figure_main_continuous_vs_single_hit.pdf|png : main figure
- figures/Figure_S1_gene_count_mechanisms.pdf|png : companion mechanistic figure
- outputs/*.csv : simulated clouds, directional residuals, and companion-figure source data
- inputs/*.csv : exact locus-level inputs and fitted parameter tables

Run
---
1. Rscript code/install_packages.R
2. Rscript code/make_model_behavior_figure.R
