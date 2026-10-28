# Corrected directional CNV models for height and BMI

This package reruns the full analysis after correcting the direction of asymmetric effects:

- **Height:** asymmetric/generalized-dysfunction effects move both copy states toward the **lower-left**.
- **BMI:** asymmetric/generalized-dysfunction effects move both copy states toward the **upper-right**.

For each trait, the package compares:

1. Mirror alone
2. Mirror + continuous directional asymmetry
3. Mirror + sparse single hit

Results fitted to the **Sacks et al. single-locus CNVs** are separated from predictions calibrated to the **Macé et al. deletion-burden coefficient** and checked against the global burden direction in **Auwerx et al.**

Primary report: `report/CNV_height_BMI_corrected_models_report.pdf`

Python reference implementation: `code/run_corrected_analysis.py`
R implementation: `code/run_corrected_analysis.R` and `code/corrected_model_functions.R`

Key corrected point estimates:

- Height continuous: q = 8.17%, shift = 0.056 SD.
- Height single hit: q = 1.95%, shift = 0.283 SD.
- BMI continuous: diffuse small upper-right load; q reaches the broad bound and should not be interpreted literally.
- BMI single hit: q = 4.26%, shift = 0.078 SD.

The numerical outputs were generated with the Python reference implementation. R was unavailable in the artifact environment, so the R workflow was structurally checked but not executed.
