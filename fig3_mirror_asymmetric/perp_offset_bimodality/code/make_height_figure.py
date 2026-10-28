#!/usr/bin/env python
"""Height main figure: four sections.

1  A-D  mirror and asymmetric CNVs differ by size and by DUP effect, not DEL
2  E-F  the perpendicular offset from the mirror line is bimodal
3  G-H  effect size vs gene count, pooled and stratified by class
4  I-M  how the three models score on the five criteria

The three models compared in section 4 share the same fitted mirror layer, and the two
directional models share the same angular freedom, so occupancy -- continuous versus a
rare single hit -- is the only thing that differs between them.

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
                help='Add a row of simulated (DEL, DUP) clouds for the two directional '
                     'models, between the dose-response row and the scorecard, and write '
                     'Figure_height_main_v2. The clouds show what each model PRODUCES; '
                     'the scorecard below scores those same models.')
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
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.titlesize': 8.4, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GOLD = '#c47b1a'

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
    ('Mirror + single hit',): 'Only mirror + single hit satisfies all five criteria',
    ('Mirror + continuous', 'Mirror + single hit'):
        'Both directional models satisfy all five criteria; mirror alone does not',
}.get(tuple(PASSES), 'How the three models score on the five criteria')
mace = pd.read_csv(IN / 'code' / 'mace_burden_reference.csv')
mace = mace[mace.trait == 'height'].set_index('burden_model')


def dens(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


CLOUDS = _a.with_clouds
# One gridspec per row. A single grid cannot do this: the rows need 2, 4, 2 (, 2) and 5
# panels, and wspace is a fraction of the average column width, so a shared 20-column
# grid caps the gap between panels at about a nineteenth of the figure -- too narrow
# for the y-axis labels in the wide rows.
# (n panels, left, right, wspace, height ratio)
SPECS = [(2, .095, .985, .22, 1.15),      # bimodality
         (4, .075, .985, .55, 1.00),      # class differences
         (2, .095, .985, .22, 1.15)]      # dose response
if CLOUDS:
    SPECS.append((2, .095, .985, .22, 1.25))   # simulated clouds
SPECS.append((5, .155, .985, .45, 1.05))       # scorecard
SCORE_ROW = len(SPECS) - 1
# Gap above every row after the first, in INCHES: it has to clear the x-axis label of
# the row above AND the section heading written into it. The scorecard's two-line panel
# titles need a wider one. (v1 keeps its original fractions exactly.)
if CLOUDS:
    FIGH = 17.0
    GAP_IN = [1.30, 1.30, 1.30, 1.85]
else:
    FIGH = 12.6
    GAP_IN = [.080 * 12.6, .080 * 12.6, .120 * 12.6]
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
# panel letters: the clouds row takes I-J, so the scorecard shifts to K-O
CLOUD_LETTERS = 'IJ' if CLOUDS else ''
SCORE = 'KLMNO' if CLOUDS else 'IJKLM'

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

# ============================== clouds (I-J), optional ==================
# What each model PRODUCES, drawn with the real gene counts and the real per-locus SEs,
# immediately above the row that scores those same two models.
if CLOUDS:
    import make_supp_model_clouds as MC        # no side effects on import
    CL = MC.cloud_data()
    for c, (key, name) in enumerate(MC.SHOW):
        ax = fig.add_subplot(ROWS[3][0, c]); AX[CLOUD_LETTERS[c]] = ax
        MC.draw_cloud(ax, key, *CL[key])
        row = MC.cmp_.loc[MC.SC.LABEL_IN_CSV[key]]
        sub = f'$\\Delta$AIC {row.delta_conditional_AIC:.1f}'
        if key == 'open_single':
            sub += f';  P(hit | 23 genes) = {1 - (1 - row.q) ** 23:.2f}'
        ax.set_title(f'{name}\n{sub}', size=8.4, pad=5, weight='bold', linespacing=1.5)

# ============================== 5. five criteria =========================
# Criteria 1, 2, 3, 5 are simulation percentiles; criterion 4 is an algebraic
# propagation to genome-wide burden and is drawn against the published CI instead.
SIM = [(0, 'sepG', '1. Size difference\nasym - mirror, log genes'),
       (1, 'sepU', '2. DUP effect difference\nasym - mirror (SD)'),
       (2, 'lrt', '3. Bimodality\nLR, 2 vs 1 component'),
       (4, 'slope', '5. DUP vs gene count\nwithin mirror class')]
BURDEN_POS = 3
FLD = {'sepG': 'sizeG', 'sepU': 'sepU', 'lrt': 'LRT', 'slope': 'slope'}
SHORT = {'Mirror alone': 'Mirror alone', 'Mirror + continuous': '+ continuous',
         'Mirror + single hit': '+ single hit'}
names = [SHORT[m] for m in sc.model]
yy = np.arange(len(names))[::-1]


def pct_colour(pc):
    """blue when the observed value sits inside the central 90% of the simulations"""
    return BLUE if .05 <= pc <= .95 else RED


for i, key, title in SIM:
    ax = fig.add_subplot(ROWS[SCORE_ROW][0, i]); AX[SCORE[i]] = ax
    p = FLD[key]
    ov = float(obs[key])
    ax.axvline(ov, color=INK, lw=1.6, zorder=1)
    for j, (_, r) in enumerate(sc.iterrows()):
        y = yy[j]
        pc = r[p + '_obs_pct']
        col = pct_colour(pc)
        ax.plot([r[p + '_lo'], r[p + '_hi']], [y, y], color=col, lw=2.4, alpha=.5,
                solid_capstyle='butt')
        ax.plot([r[p + '_med']], [y], 'o', ms=5,
                mfc=col if .05 <= pc <= .95 else 'white', mec=col, mew=1.3)
        ax.text(.5, y + .27, f'{pc:.2f}', transform=ax.get_yaxis_transform(),
                ha='center', va='bottom', size=5.8, color=col, zorder=6,
                bbox=dict(boxstyle='square,pad=.08', fc='white', ec='none'))
    ax.set_yticks(yy)
    ax.set_yticklabels(names if i == 0 else [], size=6.6)
    ax.set_title(title, size=6.9, pad=5, linespacing=1.3)
    ax.tick_params(labelsize=6.2)
    ax.set_ylim(-.8, len(names) - .15)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    if key == 'slope':
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        ax.xaxis.get_offset_text().set_size(5.8)

# criterion 4: genome-wide burden. DEL is the calibration anchor, DUP is predicted.
ax = fig.add_subplot(ROWS[SCORE_ROW][0, BURDEN_POS]); AX[SCORE[BURDEN_POS]] = ax
mlo = float(mace.loc['duplication'].ci_low_derived)
mhi = float(mace.loc['duplication'].ci_high_derived)
mb = float(mace.loc['duplication'].beta)
ax.axvspan(mlo, mhi, color=GREY, alpha=.22, lw=0, zorder=1)
ax.axvline(mb, color=INK, lw=1.6, zorder=2)
ax.axvline(0, color=LIGHT, lw=1.0, zorder=1)
for j, (_, r) in enumerate(sc.iterrows()):
    y = yy[j]
    val = r['burden_DUP_per_Mb']
    inside = mlo <= val <= mhi
    col = BLUE if inside else RED
    ax.plot([val], [y], 'o', ms=6, mfc=col if inside else 'white', mec=col, mew=1.4,
            zorder=4)
    ax.text(.5, y + .27, f'{val:+.3f}', transform=ax.get_yaxis_transform(),
            ha='center', va='bottom', size=5.8, color=col, zorder=6,
            bbox=dict(boxstyle='square,pad=.08', fc='white', ec='none'))
ax.set_yticks(yy); ax.set_yticklabels([])
ax.set_title('4. U-shaped burden\npredicted DUP (cm/Mb)', size=6.9, pad=5, linespacing=1.3)
ax.tick_params(labelsize=6.2)
ax.set_ylim(-.8, len(names) - .15)
ax.set_xlim(-.17, .21)
ax.xaxis.set_major_locator(plt.MaxNLocator(3))

# ------------------------------- panel letters and section headings -----
fig.canvas.draw()
_rend = fig.canvas.get_renderer()
_F = fig.transFigure.inverted()
_SC = 12.6 / FIGH          # offsets are fixed in INCHES, so they survive a taller figure
LETTER_TOP = {}
for L, ax in AX.items():
    p = ax.get_position()
    if L in SCORE:             # narrow panels with two-line titles
        t = fig.text(p.x0 - .026, p.y1 + .032 * _SC, L, size=10, weight='bold', va='bottom')
    else:
        t = fig.text(p.x0 - .052, p.y1 + .013 * _SC, L, size=11, weight='bold', va='bottom')
    LETTER_TOP[L] = t.get_window_extent(_rend).transformed(_F).y1


def section_top(letters):
    """Highest rendered thing in a row: the axes, its title, or its panel letter.

    Measured rather than assumed. A fixed offset silently collides as soon as a row
    gains a two-line title or the figure changes height.
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
SECTIONS.append((SCORE, CRITERIA_HEADING))
SECTIONS = [(ls, f'{n + 1}.  {t}') for n, (ls, t) in enumerate(SECTIONS)]
for letters, text in SECTIONS:
    fig.text(.045, section_top(letters) + .009 * _SC, text, size=8.8, weight='bold',
             color=INK, ha='left', va='bottom')

fig.legend(handles=[
    Line2D([], [], marker='o', color=BLUE, lw=0, ms=6, label=f'mirror ({len(MIR)})'),
    Line2D([], [], marker='s', color=RED, lw=0, ms=6, label=f'asymmetric ({len(ASY)})'),
    Line2D([], [], color=INK, lw=1.4, ls=(0, (1, 1.6)), label='all loci pooled'),
    Line2D([], [], color=INK, lw=1.6, label='observed'),
    Line2D([], [], color=BLUE, lw=2.4, alpha=.5, label='95% of 2,000 simulations')],
    frameon=False, ncol=5, loc='lower center', bbox_to_anchor=(.53, .052), fontsize=6.9)

NOTE = (f'{SCORE[0]}-{SCORE[2]}, {SCORE[4]}: the number above each bar is the percentile '
        'of the observed value in that model’s 2,000 simulations; blue = inside the '
        'central 90%.\n'
        f'{SCORE[BURDEN_POS]}: grey band is the published 95% CI (Macé et al. 2017); the '
        'deletion burden is the calibration anchor, so the duplication burden is an '
        'out-of-sample prediction.')
if CLOUDS:
    NOTE = (f'{CLOUD_LETTERS[0]}, {CLOUD_LETTERS[1]}: 150 simulated datasets per model, '
            'each locus drawn with its own standard error; observed loci in black.\n') + NOTE
if _a.note:
    NOTE += ''.join('\n' + t for t in _a.note)
elif VARIANT and CRIT != VARIANT:
    # panels A-H and I-M come from different runs, so say which note covers which
    NOTE += ('\nA-H: ' + VARIANT_NOTE.get(VARIANT, f'variant {VARIANT}.').split(';')[0] + '; the mirror line '
             'is refitted on the rest.\n' + VARIANT_NOTE.get(CRIT, f'I-M: variant {CRIT}.'))
elif VARIANT:
    NOTE += '\n' + VARIANT_NOTE.get(VARIANT, f'Variant: {VARIANT}.')
fig.text(.53, .038, NOTE, ha='center', va='top', size=5.9, color=GREY, linespacing=1.6)

fig.suptitle('Height: the asymmetry is a rare single hit acting on the duplication axis',
             size=12, weight='bold', y=.984)
# named after the panels I-M source when it differs (e.g. excl_noeffect_p05_matched)
out = HERE / 'figures' / (_a.name or ('Figure_height_main' + ('_clouds' if CLOUDS else '')
                                      + (f'_{CRIT}' if CRIT else '')))
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=200)
print('wrote', out.with_suffix('.pdf'))
