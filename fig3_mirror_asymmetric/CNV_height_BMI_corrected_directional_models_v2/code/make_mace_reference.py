#!/usr/bin/env python
"""Build the observed CNV-burden reference table from Mace et al. 2017.

Source: Mace A et al. CNV-association meta-analysis in 191,161 European adults
reveals new loci associated with anthropometric traits. Nat Commun 8:744 (2017),
doi:10.1038/s41467-017-00556-x -- Supplementary Table 1, "CNV burden analysis
results per trait and CNV categories". Betas are per Mb.

Supplementary Table 1 reports beta and P only; no standard errors are given. SEs
below are DERIVED from the two-sided Wald relation SE = |beta| / z, with
z = Phi^-1(1 - P/2). This is exact if the published P came from a Wald test on
these coefficients, which is the standard for this kind of burden regression, but
it is a reconstruction and is labelled as such wherever it is plotted.

The deletion coefficients here reproduce the values quoted in the abstract
(-0.41 cm/Mb height, +0.14 kg/m2/Mb BMI), confirming the calibration anchors used
throughout this analysis.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm

ap = argparse.ArgumentParser()
ap.add_argument('--project-dir', default=str(Path(__file__).resolve().parent.parent))
args = ap.parse_args()
P = Path(args.project_dir).resolve()

# Supplementary Table 1, verbatim. Columns: mirror, deletion, duplication, U-shaped.
ROWS = [
    # trait, unit, model, beta, pvalue
    ('height', 'cm/Mb', 'mirror', -8.37e-02, 8.56e-05),
    ('height', 'cm/Mb', 'deletion', -4.13e-01, 2.51e-10),
    ('height', 'cm/Mb', 'duplication', -7.84e-02, 6.00e-03),
    ('height', 'cm/Mb', 'U-shaped', -1.32e-01, 4.46e-07),
    ('bmi', 'kg/m2/Mb', 'mirror', 3.33e-02, 3.74e-02),
    ('bmi', 'kg/m2/Mb', 'deletion', 1.44e-01, 2.85e-03),
    ('bmi', 'kg/m2/Mb', 'duplication', 3.15e-02, 1.40e-01),
    ('bmi', 'kg/m2/Mb', 'U-shaped', 4.98e-02, 1.06e-02),
    ('weight', 'kg/Mb', 'mirror', 2.35e-02, 6.22e-01),
    ('weight', 'kg/Mb', 'deletion', 1.12e-01, 4.46e-01),
    ('weight', 'kg/Mb', 'duplication', 1.12e-02, 8.61e-01),
    ('weight', 'kg/Mb', 'U-shaped', 2.73e-02, 6.42e-01),
    ('whr', 'WHR/Mb', 'mirror', 6.93e-04, 3.37e-03),
    ('whr', 'WHR/Mb', 'deletion', 2.94e-03, 6.03e-05),
    ('whr', 'WHR/Mb', 'duplication', 8.21e-04, 9.18e-03),
    ('whr', 'WHR/Mb', 'U-shaped', 1.15e-03, 6.90e-05),
]

df = pd.DataFrame(ROWS, columns=['trait', 'units', 'burden_model', 'beta', 'p_value'])
z = norm.isf(df.p_value / 2)
df['z_derived'] = z
df['se_derived'] = np.abs(df.beta) / z
df['ci_low_derived'] = df.beta - 1.959963985 * df.se_derived
df['ci_high_derived'] = df.beta + 1.959963985 * df.se_derived
df['source'] = 'Mace et al. 2017 Nat Commun 8:744, Supplementary Table 1'
df['se_note'] = 'SE derived from beta and P via two-sided Wald; not published'

out = P / 'code' / 'mace_burden_reference.csv'
df.to_csv(out, index=False)
pd.set_option('display.width', 200)
print(df[df.trait.isin(['height', 'bmi'])]
      [['trait', 'burden_model', 'beta', 'p_value', 'se_derived',
        'ci_low_derived', 'ci_high_derived']]
      .to_string(index=False, float_format=lambda v: f'{v:.5g}'))
print(f'\nwrote {out}')
