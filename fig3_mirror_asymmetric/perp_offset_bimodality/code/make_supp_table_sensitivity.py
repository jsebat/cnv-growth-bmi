#!/usr/bin/env python
"""Supplementary table: the excl_noeffect_p05 sensitivity analysis.

Every headline number in the height figure, recomputed on the reduced locus set, beside
the value from the full set. The reduced set drops the 9 height loci where NEITHER the
deletion nor the duplication effect reaches p < 0.05 (33 -> 24 loci) and refits the
mirror layer and all three models on the remainder.

This is a sensitivity analysis and belongs only in the supplement. The main figure uses
all 33 loci.

Sources:
  full set   ../CNV_height_BMI_corrected_directional_models_v2   (model fits)
             outputs/              and outputs/me/all_loci       (bimodality)
  reduced    outputs/excl_noeffect_p05/v2_project                (model fits)
             outputs/excl_noeffect_p05/ and outputs/me/excl_noeffect_p05  (bimodality)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
MAIN = ROOT / 'CNV_height_BMI_corrected_directional_models_v2'
VAR = HERE / 'outputs' / 'excl_noeffect_p05'

MODELS = {'Mirror alone': 'Mirror alone',
          'Mirror + continuous': 'Mirror + continuous shift, free angle',
          'Mirror + single hit': 'Mirror + single hit, free angle (buffering)'}

SETS = [('all 33 loci', MAIN, HERE / 'outputs', HERE / 'outputs' / 'me' / 'all_loci'),
        ('24-locus refit', VAR / 'v2_project', VAR,
         HERE / 'outputs' / 'me' / 'excl_noeffect_p05')]

rows = []


def add(section, quantity, values):
    rows.append(dict(section=section, quantity=quantity,
                     **{k: v for k, v in values.items()}))


def fmt(v, nd=3):
    return '' if v is None or (isinstance(v, float) and not np.isfinite(v)) else f'{v:.{nd}f}'


vals = {}
for name, proj, out, me in SETS:
    h = pd.read_csv(proj / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
    h['cls'] = h.Dose_response_class.str.lower()
    cmp_ = pd.read_csv(proj / 'height' / 'outputs' /
                       'height_axis_mixture_model_comparison.csv').set_index('model')
    bur = pd.read_csv(proj / 'height' / 'outputs' /
                      'height_axis_mixture_burden_predictions.csv').set_index('model')
    mir = pd.read_csv(proj / 'height' / 'outputs' /
                      'height_mirror_parameters.csv').iloc[0]
    bim = pd.read_csv(out / 'bimodality_test.csv')
    bim = bim[(bim.trait == 'Height') & (bim.alternative == 'equal_var')
              & (bim.mirror_fit == 'mirror-class only')].iloc[0]
    mbim = pd.read_csv(me / 'bimodality_test.csv')
    mbim = mbim[(mbim.trait == 'Height') & (mbim.alternative == 'equal_var')
                & (mbim.mirror_fit == 'mirror-class only')].iloc[0]
    exp = out / 'model_expected_bimodality_v2.csv'
    exp = pd.read_csv(exp).set_index('model') if exp.exists() else None
    vals[name] = dict(h=h, cmp=cmp_, bur=bur, mir=mir, bim=bim, mbim=mbim, exp=exp)

NAMES = [s[0] for s in SETS]
V = lambda f: {n: f(vals[n]) for n in NAMES}

# -- A. cohort -------------------------------------------------------------------
add('A. Loci', 'Loci analysed', V(lambda d: len(d['h'])))
add('A. Loci', 'Mirror class', V(lambda d: int((d['h'].cls == 'mirror').sum())))
add('A. Loci', 'Asymmetric class', V(lambda d: int((d['h'].cls == 'asymmetric').sum())))

# -- B. mirror layer --------------------------------------------------------------
for k, lab in [('rho', 'Mirror slope rho'), ('p', 'Fraction of genes with an effect, p'),
               ('m', 'Per-gene effect scale, m'), ('tauU', 'Duplication dispersion tau_DUP')]:
    add('B. Mirror layer', lab, V(lambda d, k=k: float(d['mir'][k])))

# -- C. model comparison ------------------------------------------------------------
for short, full in MODELS.items():
    for col, lab in [('delta_conditional_AIC', 'dAIC'), ('delta_conditional_BIC', 'dBIC'),
                     ('delta_conditional_LOO', 'dLOO')]:
        add('C. Model comparison', f'{short}: {lab}',
            V(lambda d, f=full, c=col: float(d['cmp'].loc[f, c])))
for short, full in MODELS.items():
    if short == 'Mirror alone':
        continue
    for col, lab in [('q', 'q'), ('shift_DEL', 'shift DEL'), ('shift_DUP', 'shift DUP')]:
        add('D. Fitted directional parameters', f'{short}: {lab}',
            V(lambda d, f=full, c=col: float(d['cmp'].loc[f, c])))
    # computed, not read: the reduced-set project was written before shift_angle_deg
    # was added to the comparison table
    add('D. Fitted directional parameters', f'{short}: shift angle (deg)',
        V(lambda d, f=full: float(np.degrees(np.arctan2(
            d['cmp'].loc[f, 'shift_DUP'], d['cmp'].loc[f, 'shift_DEL'])))))

# -- E. genome-wide burden ----------------------------------------------------------
for short, full in MODELS.items():
    add('E. Predicted DUP burden (cm/Mb)', short,
        V(lambda d, f=full: float(
            d['bur'].loc[f, 'predicted_DUP_effect_per_Mb_at_matched_density'])))
add('E. Predicted DUP burden (cm/Mb)', 'Published (Mace et al. 2017)',
    {n: -0.0784 for n in NAMES})

# -- F. bimodality ------------------------------------------------------------------
add('F. Bimodality, measurement-error-aware', 'Observed likelihood ratio',
    V(lambda d: float(d['mbim'].LRT)))
add('F. Bimodality, measurement-error-aware', 'p vs unimodal null',
    V(lambda d: float(d['mbim'].p_bootstrap_nullB)))
add('F. Bimodality, plain mixture', 'Observed likelihood ratio',
    V(lambda d: float(d['bim'].LRT)))
add('F. Bimodality, plain mixture', 'p vs unimodal null',
    V(lambda d: float(d['bim'].p_bootstrap_nullB)))

# -- G. class differences and dose response ------------------------------------------
PANELS = [('CNV_length_Mb', 'CNV length', np.log),
          ('Gene_count', 'Gene count', np.log1p),
          ('Deletion_beta', 'Deletion effect', None),
          ('Duplication_beta', 'Duplication effect', None)]
for v, lab, tr in PANELS:
    def pv(d, v=v, tr=tr):
        a = d['h'][d['h'].cls == 'asymmetric'][v].to_numpy(float)
        m = d['h'][d['h'].cls == 'mirror'][v].to_numpy(float)
        a, m = (tr(a), tr(m)) if tr is not None else (a, m)
        return float(stats.ttest_ind(a, m, equal_var=False).pvalue)
    add('G. Mirror vs asymmetric (Welch p)', lab, V(pv))

for eff, lab in [('Deletion_beta', 'Deletion'), ('Duplication_beta', 'Duplication')]:
    for cls_ in ['mirror', 'asymmetric']:
        def rho(d, eff=eff, cls_=cls_):
            z = d['h'][d['h'].cls == cls_]
            return float(stats.spearmanr(z[eff], z.Gene_count)[0])
        add('H. Effect vs gene count (Spearman rho)', f'{lab}, {cls_}', V(rho))

tbl = pd.DataFrame(rows)[['section', 'quantity'] + NAMES]
OUTDIR = HERE / 'tables'
OUTDIR.mkdir(exist_ok=True)
csv = OUTDIR / 'TableS_sensitivity_excl_noeffect_p05.csv'
tbl.to_csv(csv, index=False)

# readable markdown
lines = ['# Table S — sensitivity analysis: excluding loci with no detected effect', '',
         'The reduced set drops the 9 height loci where neither the deletion nor the',
         'duplication effect reaches *p* < 0.05, leaving 24, and refits the mirror layer',
         'and all three models on the remainder. **The main figure uses all 33 loci; this',
         'is a sensitivity analysis only.**', '',
         '| | all 33 loci | 24-locus refit |', '|---|---|---|']
cur = None
for _, r in tbl.iterrows():
    if r.section != cur:
        cur = r.section
        lines.append(f'| **{cur}** | | |')
    nd = 0 if r.section.startswith('A.') else (4 if 'burden' in r.section or 'q' in str(r.quantity) else 3)
    lines.append(f'| {r.quantity} | {fmt(r[NAMES[0]], nd)} | {fmt(r[NAMES[1]], nd)} |')
md = OUTDIR / 'TableS_sensitivity_excl_noeffect_p05.md'
md.write_text('\n'.join(lines) + '\n', encoding='utf-8')
print(tbl.to_string(index=False))
print('\nwrote', csv)
print('wrote', md)
