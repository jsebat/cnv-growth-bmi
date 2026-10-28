# Continuous vs single-hit models of asymmetric reciprocal CNV effects

Working directory for the height/BMI directional-model analysis
(mirror alone vs continuous generalized dysfunction vs sparse single-hit asymmetry).

Originally delivered as two zips in Google Drive ("Simulations of Assymetric Effects").
Both have been unpacked here, the pipeline has been repaired so that it runs outside
the authoring sandbox, and the results have been verified against the shipped outputs.

## Layout

```
CNV_height_BMI_corrected_directional_models_v2/   main analysis package
  code/          drivers, model functions, QA scripts, environment files
  height/        inputs, outputs, figures  (33 reciprocal loci: 21 mirror, 12 asymmetric)
  bmi/           inputs, outputs, figures  (32 reciprocal loci)
  report/        report .pdf / .docx, summary figure
  behavior_figure_revised/   Figure_main_continuous_vs_single_hit + Figure_S1
  py_outputs/    full rerun of the Python pipeline (gitignore-able, regenerable)
  R_outputs/     R cross-check outputs (regenerable)
  qa/            R-vs-Python and rerun-vs-shipped comparison tables
CNV_burden_behavior_figure/                        Figure_burden_model_behavior (separate figure)
```

`py_outputs/`, `R_outputs/` and `qa/` are all generated. Only `code/`, the trait
`inputs/`, and the shipped `outputs/` + `report/` are source material.

## Canonical pipeline: Python

`code/run_corrected_analysis.py` is the pipeline that produced every number, figure
and table in the report. It fits the mirror model, the two asymmetry models, the
leave-one-locus-out comparison, the q-profiles, the posterior predictive checks, the
height occupancy simulation, and the burden propagation, then renders the figures
and the .docx report.

```bash
python code/run_corrected_analysis.py --project-dir . 
```

Runtime ~40 s. Writes to `py_outputs/` by default; `--out-dir DIR` to redirect, or
`--in-place` to overwrite the shipped outputs (destructive — it also rewrites
`README.md` and `SOURCES.md` inside the package).

Environment: `code/requirements.txt`. Verified with Python 3.14.7, numpy 2.5.1,
pandas 3.0.5, scipy 1.18.0, matplotlib 3.11.1, python-docx 1.2.0.

## Independent cross-check: R

`code/run_corrected_analysis.R` is a separate implementation of the core inference —
mirror fit, both asymmetry models, conditional LOO, locus predictions, burden
propagation, and the height size/gradient tests. It is a genuine second implementation,
not a wrapper, so agreement between the two is meaningful evidence.

```bash
Rscript code/install_packages.R
Rscript code/run_corrected_analysis.R --project-dir . --formats html
```

Runtime ~6 s. Writes to `R_outputs/`. Verified with R 4.6.1, ggplot2 4.0.3,
dplyr 1.2.1, tidyr 1.3.2, rmarkdown 2.31, knitr 1.51 (`code/R_sessionInfo.txt`).

Rendering the .Rmd report additionally needs pandoc, which is not installed here;
without it the CSV and figure outputs are still produced and the render step is
skipped with a warning.

The R pipeline does **not** implement the q-profiles, the posterior predictive checks,
or the height occupancy simulation. Those exist only in Python.

## Extended models (axis-aligned geometry and aggregate drift)

`code/run_axis_mixture_models.py` adds three models to the published three, all fitted
conditional on the same calibrated mirror component so the information criteria are
directly comparable: axis-aligned sparse and continuous mixtures (each affected locus
displaced along one axis only, deletion-side with probability `w`), and a reparameterised
aggregate drift that drops the unidentified event frequency. It also runs conditional LOO,
burden propagation, and posterior predictive checks for all six models.

```bash
python code/run_axis_mixture_models.py --project-dir .
```

Runtime ~15 s with LOO. `--max-q-continuous` raises the frequency bound (used to show that
BMI's `q` is pinned, not estimated); `--tag` suffixes the output filenames. The script
self-validates by refitting the two published models and comparing to the shipped fits —
it recovers them to 1.8e-15.

Outputs per trait: `{trait}_axis_mixture_model_comparison.csv`, `_conditional_LOO.csv`,
`_burden_predictions.csv`, `_posterior_predictive_checks.csv`, and `{trait}_ppc_clouds.npz`.

## Final report

`code/build_final_report.py` composes methods, results, five figures and conclusions into
`report/CNV_asymmetry_final_report.pdf` (9 pages, vector). It reads the tables and
simulation clouds produced by the two scripts above, so run both first.

```bash
python code/build_final_report.py --project-dir .
```

## Verifying a run

```bash
python code/compare_outputs.py --a . --b py_outputs --out qa/python_rerun_vs_shipped.csv
Rscript code/compare_r_vs_python.R --project-dir .
```

Both write a per-file table with the largest absolute and relative difference.
Files with no shared numeric columns are reported as `NO_NUMERIC_OVERLAP`, not as passes.

## Replication status

Verified on the repair run:

- **Python rerun vs shipped outputs**: 24 of 25 CSVs clean at 1e-4. The one flag is
  `bmi_continuous_q_profile.csv` at 4.0e-4 absolute / 5.9e-6 relative — optimizer
  tolerance on a profile curve, not a discrepancy.
- **R vs Python**: all fitted parameters, locus predictions, burden predictions and
  size/gradient tests agree. Largest difference anywhere is 3e-3 absolute, in the
  leave-one-out sums, which accumulate small optimizer differences across 33 refits.
- **Conclusions identical in both implementations**: height ΔAIC 71.96 / 2.44 / 0 and
  ΔLOO 36.90 / 0.87 / 0 (sparse single hit preferred); BMI ΔAIC 1.96 / 0 / 1.73 and
  ΔLOO 2.13 / 0 / 1.04 (continuous preferred).

See `CNV_height_BMI_corrected_directional_models_v2/REPAIR_LOG.md` for what was
broken and what changed.

## Known gaps

- `{height,BMI}_frequency_effect_profile.csv` were read from upstream `_rework_*`
  packages that were never shipped. The ridge is now **reconstructed in-package** from
  the same mirror likelihood (`frequency_effect_profile()` in the Python driver). It is
  consistent with the reported fit, but it is not byte-identical to whatever curve
  produced the original Figure 1 panel B. If the upstream packages surface, drop the
  CSVs into `{height,bmi}/outputs/` and the driver will prefer them automatically.
- The two `python_reference_generator.py` files in the figure packages are stubs
  (75 and 114 bytes) pointing at a "package creation notebook" that was not shipped.
  The behavior figures can only be regenerated from their R scripts.
- The R pipeline covers the core inference only (see above).
