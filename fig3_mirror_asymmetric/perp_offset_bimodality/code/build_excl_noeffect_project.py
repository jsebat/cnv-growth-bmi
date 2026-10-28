#!/usr/bin/env python
"""Build a copy of the v2 model project without the no-effect CNVs.

A locus is dropped when neither its deletion nor its duplication effect reaches p < 0.05
(two-sided, from beta/SE) -- the same rule as perp_offset_bimodality_excl_noeffect.py.

Writes outputs/excl_noeffect_p05/v2_project/<trait>/{inputs,outputs}/ with
  - the filtered locus table, under the original file name
  - <trait>_mirror_parameters.csv refitted on the filtered loci with the same fit_mirror
    and calibration set as the shipped one (height: mirror-class loci; BMI: all loci).
    On the full data this reproduces the shipped parameters.
  - an empty <trait>_single_locus_model_comparison.csv, so run_axis_mixture_models.py
    skips its check against the shipped (unfiltered) fits instead of reporting a mismatch

Then, from this folder:
  python ../CNV_height_BMI_corrected_directional_models_v2/code/run_axis_mixture_models.py \
      --project-dir outputs/excl_noeffect_p05/v2_project
  python code/score_models_on_criteria.py --project-dir outputs/excl_noeffect_p05/v2_project \
      --out-dir outputs/excl_noeffect_p05
  python code/make_height_figure.py --variant excl_noeffect_p05
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import norm

from perp_offset_bimodality import HERE, fit_mirror

V2 = HERE.parent / 'CNV_height_BMI_corrected_directional_models_v2'
OUT = HERE / 'outputs' / 'excl_noeffect_p05' / 'v2_project'
P_CUTOFF = 0.05
TRAITS = [('height', 'sacks_height_all_loci_exact.csv', 0.30, True),
          ('bmi', 'sacks_bmi_complete_reciprocal_loci_exact.csv', 11 / 14, False)]

for trait, loci_file, pi_pos, mirror_class_only in TRAITS:
    dat = pd.read_csv(V2 / trait / 'inputs' / loci_file)
    pD = 2 * norm.sf(np.abs(dat.Deletion_beta / dat.Deletion_SE))
    pU = 2 * norm.sf(np.abs(dat.Duplication_beta / dat.Duplication_SE))
    kept = dat[(pD < P_CUTOFF) | (pU < P_CUTOFF)].reset_index(drop=True)
    fit_dat = (kept[kept.Dose_response_class.str.lower() == 'mirror']
               if mirror_class_only else kept)
    base = fit_mirror(fit_dat, pi_pos)
    base.update(pi_positive=pi_pos, n_loci_calibration=len(fit_dat))

    (OUT / trait / 'inputs').mkdir(parents=True, exist_ok=True)
    (OUT / trait / 'outputs').mkdir(parents=True, exist_ok=True)
    kept.to_csv(OUT / trait / 'inputs' / loci_file, index=False)
    pd.DataFrame([base]).to_csv(OUT / trait / 'outputs' / f'{trait}_mirror_parameters.csv',
                                index=False)
    pd.DataFrame(columns=['model', 'negative_log_likelihood']).to_csv(
        OUT / trait / 'outputs' / f'{trait}_single_locus_model_comparison.csv', index=False)
    print(f'{trait}: kept {len(kept)} of {len(dat)} loci; mirror layer on {len(fit_dat)} loci: '
          + ', '.join(f'{k}={v:.4g}' for k, v in base.items() if isinstance(v, float)))
