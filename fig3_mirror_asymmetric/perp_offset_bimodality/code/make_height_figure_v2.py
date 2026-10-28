#!/usr/bin/env python
"""Height main figure, v2 (Arial; section 4 reduced to where the models differ).

1  A-B  the perpendicular offset from the mirror line is bimodal
2  C-F  mirror and asymmetric CNVs differ by size and by DUP effect, not DEL
3  G-H  effect size vs gene count, pooled and stratified by class
4  I-K  where the three models differ:
        I  predicted height change against genome-wide duplication burden, each model
           calibrated so its deletion burden matches Mace et al. 2017
        J  the offset distribution each model expects (pooled over its 2,000 simulated
           datasets) over the observed offsets, with the p value each model's datasets
           get from the same bimodality test against the unimodal null

The three models compared in section 4 share the same fitted mirror layer, and the two
directional models share the same angular freedom, so occupancy -- continuous versus a
rare single hit -- is the only thing that differs between them.

The simulated datasets in J are regenerated with score_models_on_criteria.simulate and
its seed, in its model order, so they are the very replicates whose LRTs are stored in
model_criteria_draws.npz (verified 2026-10-14: refitting the first 60 replicates of each
model reproduces the stored LRTs, r >= 0.997; the residue is EM random-start variation). Each replicate's p is its LRT referred to the
observed null-B bootstrap distribution (panel B).

Section headings are placed from the realised axes positions, so they cannot collide
with the panel titles above them.

--variant NAME draws from outputs/NAME/ (and its v2_project/ copy of the model project)
and writes Figure_height_main_NAME; e.g. excl_noeffect_p05, built by
build_excl_noeffect_project.py. With no variant the original figure is drawn.
--criteria-variant NAME takes panels I-M from outputs/NAME/ instead (e.g.
excl_noeffect_p05_matched, from score_models_on_criteria.py --match-filter 0.05).
--bimodality-dir, --loci-dir, --criteria-dir (paths under outputs/, or the v2 project for
loci), --name and --note set every source explicitly and override the variant options;
e.g. the measurement-error-aware figures read outputs/me/... (see me_bimodality.py).
"""
from __future__ import annotations
import argparse, pickle, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))   # python -I drops the script directory from sys.path
HERE = CODE.parent
ROOT = HERE.parent
IN = ROOT / 'CNV_height_BMI_corrected_directional_models_v2'
ALT = 'equal_var'

ap = argparse.ArgumentParser()
ap.add_argument('--variant', default='')
ap.add_argument('--criteria-variant', default=None)
_a, _ = ap.parse_known_args()      # the source overrides below are parsed fully later
VARIANT = _a.variant
CRIT = _a.criteria_variant if _a.criteria_variant is not None else VARIANT
OUT = HERE / 'outputs' / VARIANT if VARIANT else HERE / 'outputs'
OUT_CRIT = HERE / 'outputs' / CRIT if CRIT else HERE / 'outputs'
LOCI = OUT / 'v2_project' if VARIANT else IN
for opt in ('--bimodality-dir', '--loci-dir', '--criteria-dir', '--name'):
    ap.add_argument(opt, default=None)
ap.add_argument('--note', action='append', default=None)
ap.add_argument('--with-clouds', action='store_true',
                help='Insert a row of simulated (DEL, DUP) clouds for the two '
                     'directional models directly above section 4, and write '
                     'Figure_height_main_v2_clouds. The clouds show what each model '
                     'PRODUCES; section 4 then shows where those models differ.')
_a = ap.parse_args()
if _a.bimodality_dir: OUT = HERE / 'outputs' / _a.bimodality_dir
if _a.loci_dir: LOCI = HERE / 'outputs' / _a.loci_dir
if _a.criteria_dir: OUT_CRIT = HERE / 'outputs' / _a.criteria_dir
VARIANT_NOTE = {'excl_noeffect_p05': 'Excludes loci where neither DEL nor DUP reaches p < 0.05; '
                                     'the mirror layer and all three models are refitted on the rest.',
                'excl_noeffect_p05_matched': 'I-M: models fitted on all 33 loci; the same p < 0.05 '
                                             'exclusion is applied to the observed data and to every '
                                             'simulated replicate before scoring.'}

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 8,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'axes.titlesize': 8.4, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GOLD = '#c47b1a'; PURPLE = '#6a3d9a'

h = pd.read_csv(LOCI / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
po = pd.read_csv(OUT / 'perp_offsets_per_locus.csv')
po = po[po.trait == 'Height'].tail(len(h))     # the mirror-class calibration comes last
h = h.merge(po[['CNV_ID', 'perp_offset']], on='CNV_ID')
h['cls'] = h.Dose_response_class.str.lower()
MIR = h[h.cls == 'mirror']; ASY = h[h.cls == 'asymmetric']
GRP = [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]

res = pickle.loads((OUT / '_res.pkl').read_bytes())
R = [r for r in res if r['mirror_subset'] and r['name'] == 'Height'][0]
tab = pd.read_csv(OUT / 'bimodality_test.csv')
ROW = tab[(tab.trait == 'Height') & (tab.mirror_fit == 'mirror-class only')
          & (tab.alternative == ALT)].iloc[0]
sc = pd.read_csv(OUT_CRIT / 'model_criteria_scores.csv')
obs = pd.read_csv(OUT_CRIT / 'model_criteria_observed.csv', index_col=0).iloc[:, 0]
BIMODAL = ROW.p_bootstrap_nullB < .05
# a model passes when every simulated criterion's percentile is inside the central 90%
# and its predicted burden is U-shaped
PASSES = sc.model[sc[[c + '_obs_pct' for c in ('sizeG', 'sepU', 'LRT', 'slope')]]
                  .apply(lambda r: r.between(.05, .95).all(), axis=1) & sc.burden_U_shaped]
CRITERIA_HEADING = {
    ('Mirror + single hit',): '4.  Only mirror + single hit satisfies all five criteria',
    ('Mirror + continuous', 'Mirror + single hit'):
        '4.  Both directional models satisfy all five criteria; mirror alone does not',
}.get(tuple(PASSES), '4.  How the three models score on the five criteria')
mace = pd.read_csv(IN / 'code' / 'mace_burden_reference.csv')
mace = mace[mace.trait == 'height'].set_index('burden_model')


def dens(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


CLOUDS = _a.with_clouds
# One gridspec per row. A single grid cannot do this: the rows need 2, 4, 2 (, 2) and 3
# panels, and wspace is a fraction of the average column width, so a shared grid caps
# the gap between panels far too narrow for the y-axis labels in the wide rows.
# (n panels, left, right, wspace, height ratio)
SPECS = [(2, .095, .985, .22, 1.15),      # bimodality
         (4, .075, .985, .55, 1.00),      # class differences
         (2, .095, .985, .22, 1.15)]      # dose response
if CLOUDS:
    SPECS.append((2, .095, .985, .22, 1.25))   # simulated clouds
SPECS.append((3, .085, .985, .50, 1.15))       # where the models differ
S4_ROW = len(SPECS) - 1
# Gap above every row after the first. With the clouds row the figure is taller, so the
# gaps are set in INCHES (a fixed fraction would shrink them); without it the original
# fractions are kept exactly.
if CLOUDS:
    FIGH = 17.0
    GAP_IN = [1.30, 1.30, 1.30, 1.45]
else:
    FIGH = 12.6
    GAP_IN = [.080 * 12.6, .080 * 12.6, .085 * 12.6]
fig = plt.figure(figsize=(7.5, FIGH))
ROWS = {}
_top, _bot = .915, .095
_gaps = [g / FIGH for g in GAP_IN]
_h = [sp[4] for sp in SPECS]
_unit = (_top - _bot - sum(_gaps)) / sum(_h)
_y = _top
for _i, (_n, _l, _r, _w, _) in enumerate(SPECS):
    _b = _y - _unit * _h[_i]
    ROWS[_i] = fig.add_gridspec(1, _n, left=_l, right=_r, top=_y, bottom=_b, wspace=_w)
    _y = _b - (_gaps[_i] if _i < len(_gaps) else 0)
AX = {}
# panel letters: the clouds row takes I-J, so section 4 shifts to K-M
CLOUD_LETTERS = 'IJ' if CLOUDS else ''
S4 = 'KLM' if CLOUDS else 'IJK'

# ============================== 2. class differences (C-F) ==============
rng = np.random.default_rng(3)
PANELS = [('CNV_length_Mb', 'CNV length (Mb)', np.log, True),
          ('Gene_count', 'Gene count', np.log1p, True),
          ('Deletion_beta', 'Deletion effect (SD)', None, False),
          ('Duplication_beta', 'Duplication effect (SD)', None, False)]
for i, (v, lab, tr, from_zero) in enumerate(PANELS):
    ax = fig.add_subplot(ROWS[1][0, i]); AX['CDEF'[i]] = ax
    a = ASY[v].to_numpy(float); m = MIR[v].to_numpy(float)
    aa, mm = (tr(a), tr(m)) if tr is not None else (a, m)
    p = stats.ttest_ind(aa, mm, equal_var=False).pvalue
    for j, (vals, cc) in enumerate([(m, BLUE), (a, RED)]):
        ax.scatter(np.full(len(vals), j) + rng.normal(0, .085, len(vals)), vals, s=9,
                   color=cc, alpha=.8, lw=0)
        ax.plot([j - .30, j + .30], [np.median(vals)] * 2, color=INK, lw=1.5, zorder=4)
    span = max(a.max(), m.max()) - min(a.min(), m.min())
    if from_zero:
        ax.set_ylim(0, max(a.max(), m.max()) * 1.12)
    else:
        ax.set_ylim(min(a.min(), m.min()) - .07 * span, max(a.max(), m.max()) + .12 * span)
        ax.axhline(0, color=LIGHT, lw=.6, zorder=1)
    ax.set_xlim(-.62, 1.62); ax.set_xticks([0, 1]); ax.set_xticklabels(['mir', 'asy'])
    ax.set_ylabel(lab, size=7.6)
    ptxt = f'$p$ = {p:.3f}' if p >= 1e-4 else f'$p$ = {p:.0e}'
    ax.set_title(ptxt, size=8, pad=5, color=RED if p < .05 else INK,
                 weight='bold' if p < .01 else 'normal')

# ============================== 1. bimodality (A-B) =====================
ax = fig.add_subplot(ROWS[0][0, 0]); AX['A'] = ax
x = R['x']; f = R['fits'][ALT]
lo, hi = x.min() - .10, x.max() + .10
ax.hist(x, bins=np.linspace(lo, hi, 15), density=True, color=BLUE, alpha=.20,
        edgecolor='white', lw=.5)
g = np.linspace(lo, hi, 500)
if 'me' in R:
    # measurement-error-aware fit: each locus has variance tau^2 + s_i^2, so the curve
    # drawn is the average of the loci's own predictive densities
    S = R['me']['S'][:, None]; mu1, t1 = R['me']['one']; pi_, m1, m2, t2 = R['me']['two']
    one = dens(g[None, :], mu1, np.sqrt(t1 ** 2 + S ** 2)).mean(0)
    two = (pi_ * dens(g[None, :], m1, np.sqrt(t2 ** 2 + S ** 2))
           + (1 - pi_) * dens(g[None, :], m2, np.sqrt(t2 ** 2 + S ** 2))).mean(0)
    ax.plot(g, one, color=GREY, lw=1.5, ls='--', label='1 component')
    ax.plot(g, two, color=BLUE, lw=2.0, label='2 components')
else:
    ax.plot(g, dens(g, R['mu'], R['sd']), color=GREY, lw=1.5, ls='--', label='1 component')
    pi_, m1, s1, m2, s2 = f['par']
    ax.plot(g, pi_ * dens(g, m1, s1) + (1 - pi_) * dens(g, m2, s2), color=BLUE, lw=2.0,
            label='2 components')
cls = np.asarray([str(s).lower() for s in R['cls']])
ymax = ax.get_ylim()[1]
for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
    s = cls == lab
    ax.plot(x[s], np.full(int(s.sum()), -ymax * .06), '|', color=cc, ms=7, mew=1.3,
            clip_on=False, label=f'{lab} ({int(s.sum())})')
ax.set_ylim(-ymax * .11, ymax * 1.05)
ax.set_xlabel('Perpendicular offset from mirror line (SD)')
ax.set_ylabel('Density')
ax.set_title('Offsets are bimodal' if BIMODAL else 'Offsets: bimodality not significant',
             size=8.4, pad=5, weight='bold')
ax.legend(frameon=False, loc='upper left', handletextpad=.5, labelspacing=.22, fontsize=6.8)

ax = fig.add_subplot(ROWS[0][0, 1]); AX['B'] = ax
nullB, o = f['nullB'], f['obs']
top = max(np.quantile(nullB, .995), o) * 1.10
ax.hist(nullB, bins=np.linspace(0, top, 40), density=True, color=GREY, alpha=.45,
        edgecolor='white', lw=.4)
ax.axvline(o, color=BLUE, lw=2.2)
ax.set_xlim(0, top)
ax.set_xlabel('Likelihood ratio, 2 vs 1 component')
ax.set_ylabel('Density')
ax.set_title(f'vs unimodal null:  $p$ = {ROW.p_bootstrap_nullB:.3f}', size=8.4, pad=5,
             weight='bold', color=RED if BIMODAL else INK)
# keep the label clear of the observed line: move it left of the line when the line
# would run through the label's usual upper-right spot
_r = o / top
_x, _y, _ha = (.97, .92, 'right') if not .66 < _r < .99 else (_r - .02, .62, 'right')
ax.text(_x, _y, f'observed {o:.1f}\n{len(nullB):,} null draws', transform=ax.transAxes,
        ha=_ha, va='top', size=7, color=INK)

# ============================== 3. dose response (G-H) ==================
for c, (eff, elab) in enumerate([('Deletion_beta', 'Deletion effect (SD)'),
                                 ('Duplication_beta', 'Duplication effect (SD)')]):
    ax = fig.add_subplot(ROWS[2][0, c]); AX['GH'[c]] = ax
    v = 'Gene_count'
    for cls_, col, mk in GRP:
        z = h[h.cls == cls_]
        ax.scatter(z[v], z[eff], s=20, marker=mk, color=col, alpha=.85,
                   edgecolor='white', lw=.4)
        sl, ic = stats.theilslopes(z[eff], z[v])[:2]
        xx = np.linspace(z[v].min(), z[v].max(), 40)
        rr, pp = stats.spearmanr(z[eff], z[v])
        ax.plot(xx, ic + sl * xx, color=col, lw=1.9 if pp < .05 else 1.0,
                ls='-' if pp < .05 else (0, (4, 2)), zorder=4)
    sl, ic = stats.theilslopes(h[eff], h[v])[:2]
    xx = np.linspace(h[v].min(), h[v].max(), 40)
    ra, pa = stats.spearmanr(h[eff], h[v])
    ax.plot(xx, ic + sl * xx, color=INK, lw=1.4, ls=(0, (1, 1.6)), zorder=2)
    ax.axhline(0, color=LIGHT, lw=.7)
    ax.set_xlabel('Gene count'); ax.set_ylabel(elab)
    ents = [('all ', INK, ra, pa)]
    for cls_, col, _ in GRP:
        z = h[h.cls == cls_]
        rr, pp = stats.spearmanr(z[eff], z[v])
        ents.append((cls_[:4], col, rr, pp))
    y0 = .30 if eff == 'Deletion_beta' else .99
    for k, (lbl, col, rr, pp) in enumerate(ents):
        ax.text(.03, y0 - k * .092, f'{lbl}  $\\rho$={rr:+.2f}  $p$={pp:.3f}',
                transform=ax.transAxes, ha='left', va='top', size=6.9, color=col,
                weight='bold' if pp < .05 else 'normal')
    ax.set_title('deletion' if c == 0 else 'duplication', size=8.4, pad=5)

# ============================== 4. where the models differ (I-K) =========
import score_models_on_criteria as SC          # noqa: E402  (simulate + its seed)
MODEL_STYLE = [('mirror', 'Mirror alone', '#4d4d4d', (0, (4, 2))),
               ('open_continuous', '+ continuous', GOLD, '-'),
               ('open_single', '+ single hit', PURPLE, '-')]

# ---- I: predicted height change against duplication burden -----------------------
ax = fig.add_subplot(ROWS[S4_ROW][0, 0]); AX[S4[0]] = ax
mb = float(mace.loc['duplication'].beta)
mlo = float(mace.loc['duplication'].ci_low_derived)
mhi = float(mace.loc['duplication'].ci_high_derived)
mp = float(mace.loc['duplication'].p_value)
XMAX = 4.0
xb = np.linspace(0, XMAX, 50)
ax.fill_between(xb, mlo * xb, mhi * xb, color=GREY, alpha=.22, lw=0, zorder=1)
ax.plot(xb, mb * xb, color=INK, lw=1.8, zorder=3)
ax.axhline(0, color=LIGHT, lw=.8, zorder=0)
_scm = sc.set_index('key')
hd = []
for key, lab, col, ls in MODEL_STYLE:
    s = float(_scm.loc[key, 'burden_DUP_per_Mb'])
    ax.plot(xb, s * xb, color=col, lw=2.0, ls=ls, zorder=4)
    hd.append(Line2D([], [], color=col, lw=2.0, ls=ls, label=f'{lab}  {s:+.3f}'))
hd.append(Line2D([], [], color=INK, lw=1.8, label=f'observed  {mb:+.3f}'))
hd.append(Patch(color=GREY, alpha=.22, lw=0, label='observed 95% CI'))
ax.legend(handles=hd, frameon=False, loc='upper left', fontsize=6.1, handlelength=2.0,
          labelspacing=.22, title='slope, cm/Mb', title_fontsize=6.1, alignment='left')
ax.set_xlim(0, XMAX)
ax.set_ylim(-.62, 1.15)
ax.set_xlabel('Duplication burden (Mb)')
ax.set_ylabel('Predicted change in height (cm)')
ax.set_title('DUP burden and short stature', size=8.4, pad=5, weight='bold')

# ---- simulated datasets, regenerated exactly as the scoring run drew them ----------
_cmp = pd.read_csv(LOCI / 'height' / 'outputs' /
                   'height_axis_mixture_model_comparison.csv').set_index('model')
_base = pd.read_csv(LOCI / 'height' / 'outputs' / 'height_mirror_parameters.csv').iloc[0].to_dict()
_rho = _base['rho']; _k = np.sqrt(1 + _rho * _rho)
assert abs(_rho - R['rho']) < 1e-6, 'offsets in A and the simulations use different lines'
_G = h.Gene_count.to_numpy(float)
_draws = np.load(OUT_CRIT / 'model_criteria_draws.npz')
_rng = np.random.default_rng(SC.SEED)
SIMX, EXPECTED = {}, []
for key, lab, col, ls in MODEL_STYLE:          # same order and rng as the scoring run
    pars = _cmp.loc[SC.LABEL_IN_CSV[key]].to_dict()
    Ds, Us = SC.simulate(key, pars, _G, h.Deletion_SE.to_numpy(float),
                         h.Duplication_SE.to_numpy(float), _base, -1,
                         len(_draws[key + '_lrt']), _rng)
    SIMX[key] = ((_rho * Ds + Us) / _k).ravel()
    lrt = _draws[key + '_lrt']
    pv = (nullB[None, :] >= lrt[:, None]).mean(1)
    EXPECTED.append(dict(model=lab, key=key, median_LRT=float(np.median(lrt)),
                         median_p_vs_unimodal=float(np.median(pv)),
                         frac_p_lt_05=float(np.mean(pv < .05)),
                         P_LRT_ge_observed=float(np.mean(lrt >= o)),
                         observed_percentile=float(np.mean(lrt <= o))))
E = {e['key']: e for e in EXPECTED}

# ---- J: the offset distribution each model expects --------------------------------
ax = fig.add_subplot(ROWS[S4_ROW][0, 1]); AX[S4[1]] = ax
lo, hi = -.62, .32
g = np.linspace(lo, hi, 400)
ax.hist(x, bins=np.linspace(x.min() - .10, x.max() + .10, 15), density=True,
        color=BLUE, alpha=.16, edgecolor='white', lw=.5, zorder=1)
for key, lab, col, ls in MODEL_STYLE:
    ax.plot(g, stats.gaussian_kde(SIMX[key], bw_method=.12)(g), color=col, lw=2.0, ls=ls,
            zorder=4)
ymax = ax.get_ylim()[1]
for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
    s = cls == lab
    ax.plot(x[s], np.full(int(s.sum()), -ymax * .06), '|', color=cc, ms=7, mew=1.3,
            clip_on=False)
ax.set_ylim(-ymax * .11, ymax * 1.05)
ax.set_xlim(lo, hi)
ax.set_xlabel('Perpendicular offset (SD)')
ax.set_ylabel('Density')
ax.set_title('Expected offset distribution', size=8.4, pad=5, weight='bold')
hd = [Line2D([], [], color=BLUE, lw=6, alpha=.25, label='observed')]
hd += [Line2D([], [], color=col, lw=2.0, ls=ls, label=lab) for key, lab, col, ls in MODEL_STYLE]
ax.legend(handles=hd, frameon=False, loc='upper left', fontsize=6.1, handlelength=2.0,
          labelspacing=.22)

# ---- K: likelihood ratio each model's datasets produce, against the unimodal null ---
ax = fig.add_subplot(ROWS[S4_ROW][0, 2]); AX[S4[2]] = ax
_all = np.concatenate([nullB] + [_draws[k + '_lrt'] for k in SIMX])
top = max(np.quantile(_all, .99), o) * 1.05
gl = np.linspace(0, top, 400)
ax.hist(nullB, bins=np.linspace(0, top, 36), density=True, color=GREY, alpha=.35,
        edgecolor='white', lw=.3, zorder=1)
for key, lab, col, ls in MODEL_STYLE:
    lr = _draws[key + '_lrt']
    kd = stats.gaussian_kde(np.concatenate([lr, -lr]), bw_method=.10)   # reflected at 0
    ax.plot(gl, 2 * kd(gl), color=col, lw=2.0, ls=ls, zorder=4)
ax.axvline(o, color=BLUE, lw=2.0, zorder=3.5)
ax.set_xlim(0, top)
ax.set_xlabel('Likelihood ratio, 2 vs 1 component')
ax.set_ylabel('Density')
ax.set_title('Expected bimodality', size=8.4, pad=5, weight='bold')
# legend carries each model's typical p against the unimodal null (median over its
# simulated datasets)
hd = [Line2D([], [], color=BLUE, lw=2.0,
             label='observed'),
      Patch(color=GREY, alpha=.35, lw=0, label='unimodal null')]
hd += [Line2D([], [], color=col, lw=2.0, ls=ls,
              label=f'{lab}: median $p$ = {E[key]["median_p_vs_unimodal"]:.3f}')
       for key, lab, col, ls in MODEL_STYLE]
ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
leg = ax.legend(handles=hd, frameon=True, loc='upper right', fontsize=6.0, handlelength=2.0,
                labelspacing=.22, facecolor='white', edgecolor='none', framealpha=1,
                borderpad=.3)
leg.set_zorder(10)
pd.DataFrame(EXPECTED).to_csv(OUT_CRIT / 'model_expected_bimodality_v2.csv', index=False)
print(pd.DataFrame(EXPECTED).to_string(index=False))

# ============================== clouds (I-J), optional ==================
# What each model PRODUCES, with the real gene counts and the real per-locus SEs,
# directly above the row showing where those same two models differ.
if CLOUDS:
    # make_supp_model_clouds reads the model fits from the MAIN project, so with
    # --variant the clouds would come from the full 33-locus fit while every other
    # panel came from the variant's refit. Refuse rather than draw mixed provenance.
    if VARIANT:
        raise SystemExit(
            f'--with-clouds cannot be combined with --variant {VARIANT}: the cloud '
            'panels would be simulated from the main 33-locus fit while the other '
            'panels use the variant refit. Teach cloud_data() to take a project '
            'directory before using these together.')
    import make_supp_model_clouds as MC        # no side effects on import
    CL = MC.cloud_data()
    for _c, (_key, _name) in enumerate(MC.SHOW):
        ax = fig.add_subplot(ROWS[3][0, _c]); AX[CLOUD_LETTERS[_c]] = ax
        MC.draw_cloud(ax, _key, *CL[_key])
        MC.observed_key(ax)       # no shared legend in this build, so name them here
        _row = MC.cmp_.loc[MC.SC.LABEL_IN_CSV[_key]]
        _sub = f'$\\Delta$AIC {_row.delta_conditional_AIC:.1f}'
        if _key == 'open_single':
            _sub += f';  P(hit | 23 genes) = {1 - (1 - _row.q) ** 23:.2f}'
        ax.set_title(f'{_name}\n{_sub}', size=8.6, pad=5, weight='bold', linespacing=1.5)

# ------------------------------- panel letters and section headings -----
fig.canvas.draw()
_rend = fig.canvas.get_renderer()
_F = fig.transFigure.inverted()
_SC = 12.6 / FIGH          # offsets stay fixed in INCHES as the figure grows
LETTER_TOP = {}
for L, ax in AX.items():
    p = ax.get_position()
    t = fig.text(p.x0 - .052, p.y1 + .013 * _SC, L, size=11, weight='bold', va='bottom')
    LETTER_TOP[L] = t.get_window_extent(_rend).transformed(_F).y1


def section_top(letters):
    """Highest rendered thing in a row: the axes, its title, or its panel letter.

    Measured rather than assumed, so a row that gains a two-line title (the clouds)
    cannot silently push its heading into the row above.
    """
    tops = []
    for L in letters:
        ax = AX[L]
        t = ax.get_position().y1
        if ax.title.get_text():
            t = max(t, ax.title.get_window_extent(_rend).transformed(_F).y1)
        tops.append(max(t, LETTER_TOP[L]))
    return max(tops)

SECTIONS = [('AB', 'The perpendicular offset from the mirror line is bimodal' if BIMODAL
             else 'Is the perpendicular offset from the mirror line bimodal?'),
            ('CDEF', 'Mirror and asymmetric CNVs differ by size and DUP effect, not DEL'),
            ('GH', 'Effect size scales with gene count only after stratifying by class')]
if CLOUDS:
    SECTIONS.append((CLOUD_LETTERS, 'What the two directional models produce'))
SECTIONS.append((S4, 'Where the models differ: duplication burden and bimodality'))
SECTIONS = [(ls, f'{n + 1}.  {t}') for n, (ls, t) in enumerate(SECTIONS)]
for letters, text in SECTIONS:
    fig.text(.045, section_top(letters) + .009 * _SC, text, size=8.8, weight='bold',
             color=INK, ha='left', va='bottom')

# The shared legend sits far from everything it describes (the classes in C-H), which
# makes it read as if it belonged to the panels just above it. In the clouds build every
# panel names its own marks in place, so it is dropped there.
if not CLOUDS:
    fig.legend(handles=[
        Line2D([], [], marker='o', color=BLUE, lw=0, ms=6, label=f'mirror ({len(MIR)})'),
        Line2D([], [], marker='s', color=RED, lw=0, ms=6, label=f'asymmetric ({len(ASY)})'),
        Line2D([], [], color=INK, lw=1.4, ls=(0, (1, 1.6)), label='all loci pooled')],
        frameon=False, ncol=3, loc='lower center', bbox_to_anchor=(.53, .044), fontsize=6.9)

NOTE = (f'{S4[0]}: each model is calibrated so its deletion burden equals the published '
        '-0.41 cm/Mb (Macé et al. 2017); its duplication slope is then an out-of-sample '
        'prediction. Grey band = published 95% CI.\n'
        f'{S4[1]}, {S4[2]}: each model’s 2,000 simulated datasets (real gene counts and SEs). '
        f'{S4[1]} pools their offsets; {S4[2]} shows their likelihood ratios. Median $p$ = the typical '
        'dataset’s $p$ against the unimodal null in B.')
if _a.note:
    NOTE += ''.join('\n' + t for t in _a.note)
elif VARIANT:
    NOTE += '\n' + VARIANT_NOTE.get(VARIANT, f'Variant: {VARIANT}.')
fig.text(.53, .038, NOTE, ha='center', va='top', size=5.9, color=GREY, linespacing=1.6)

fig.suptitle('Height: the asymmetry is a rare single hit acting on the duplication axis',
             size=12, weight='bold', y=.984)
# named after the panels I-M source when it differs (e.g. excl_noeffect_p05_matched)
out = HERE / 'figures' / (_a.name or ('Figure_height_main_v2'
                                      + ('_clouds' if CLOUDS else '')
                                      + (f'_{CRIT}' if CRIT else '')))
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=200)
print('wrote', out.with_suffix('.pdf'))
