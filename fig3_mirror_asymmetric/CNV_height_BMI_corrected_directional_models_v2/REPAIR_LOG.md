# Repair log

Changes made to get the delivered package running and verified outside the
authoring sandbox. Every change is listed; nothing about the models themselves
was altered, and no reported result changed.

## Python driver — `code/run_corrected_analysis.py`

The shipped script could not run anywhere except the machine it was written on.

1. **Hardcoded sandbox paths.** `ROOT` was fixed to
   `/mnt/data/CNV_height_BMI_corrected_directional_models`, and inputs were copied
   from `/mnt/data/_rework_height/...` and `/mnt/data/_rework_bmi/...`, neither of
   which ships with the package. Those copies were guarded by `if src.exists()`, so
   they failed silently and the script then died at the first unguarded
   `pd.read_csv`. Replaced with `--project-dir` / `--out-dir` / `--in-place`
   arguments; inputs are now read from the package itself, with `--height-src` /
   `--bmi-src` available to point at external sources. Missing inputs now raise a
   named error instead of failing silently.

2. **`shutil.rmtree(ROOT)` on startup.** Harmless when ROOT was a scratch directory;
   catastrophic once ROOT can point at the package. Now only ever clears a dedicated
   output directory, never `--project-dir`.

3. **`shutil.copy2(__file__, ROOT/'code/run_corrected_analysis.py')`** would raise
   `SameFileError` when running in place. Now skipped when source and destination
   resolve to the same file.

4. **Missing upstream input: frequency-effect ridge.** Figure 1 panel B for both
   traits read `{height,BMI}_frequency_effect_profile.csv` from the unshipped
   `_rework_*` packages. Added `frequency_effect_profile()`, which profiles the
   mirror likelihood over contributor frequency `p`, optimising `(m, rho, tauD, tauU)`
   at each fixed `p`. `load_frequency_effect_profile()` prefers a shipped file if one
   is present and otherwise reconstructs and saves the curve. **This is the one place
   where output is not guaranteed identical to the original** — the reconstruction is
   consistent with the reported mirror fit, but the original curve is unavailable for
   comparison.

5. **Overflow in the new profile.** Nelder-Mead is unconstrained and pushed the
   log-parameters far enough to overflow `exp`, raising
   `OverflowError: (34, 'Result too large')`. Clamped the log-parameters to
   [-50, 50] and made non-finite objective values return the sentinel penalty.

6. **Cloud-sync file locks.** `rmtree` failed with `WinError 32` because Dropbox
   holds handles on files it is syncing. Added `_rmtree_resilient()`: retries, then
   degrades to overwrite-in-place with a warning rather than aborting the run.

## R model functions — `code/corrected_model_functions.R`

7. **LOO folds were labelled by row index.** `held_out_row = i` with column
   `log_predictive_density`, where Python writes `held_out_locus` (e.g.
   `10q11.21_11.23`) and `held_out_log_predictive_density`. The integer index was
   meaningless outside the run and made the two outputs impossible to join. Now uses
   `CNV_ID` and the Python column names.

8. **`single_hit_probability` was `NA` for the mirror model**, where Python writes
   `0`. NA would silently propagate through any sum or mean of that column.
   Continuous stays `NA`, matching Python.

9. **Missing columns.** Added `mirror_expected_DEL` / `mirror_expected_DUP` to the
   locus-prediction tables and `max_q` to the model-comparison table (with
   `max_q = NA_real_` on the mirror branch, which had no such field).

## R driver — `code/run_corrected_analysis.R`

10. **Output filename mismatch.** Wrote `*_locus_predictions_single.csv` where Python
    writes `*_locus_predictions_single_hit.csv`, so the two trees did not line up.

## QA scripts (new)

11. `code/compare_r_vs_python.R` and `code/compare_outputs.py` — per-file numeric
    diffs between output trees. The Python driver created a `qa/` directory but never
    wrote anything to it, so no cross-check had ever been shipped.

12. **False pass in the first version of the R comparison.** Because R's LOO columns
    were named differently (defect 7), there were zero shared numeric columns and the
    file was reported `OK` without ever being compared. Zero overlap is now a hard
    `NO_NUMERIC_OVERLAP` failure.

## Environment (new)

13. `code/requirements.txt` and `code/R_sessionInfo.txt` — nothing pinned any
    versions previously. `code/install_packages.R` still installs latest-from-CRAN;
    the session file records what the verified run actually used.

## Not changed

- No model, likelihood, prior, seed, or parameter was modified.
- The seeds were already correct and explicit throughout, in both implementations.
- The line in the generated README stating that "R was unavailable in the artifact
  environment, so the R workflow was structurally checked but not executed" is now
  out of date — the R workflow has been executed and cross-checked. It is left as-is
  because it is generated text describing the original run.
