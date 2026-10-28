#!/usr/bin/env python
"""Manuscript figure: six panels.

A, B  perpendicular offsets from the mirror line, with 1- and 2-component fits
C, D  observed likelihood ratio against the bootstrap null
E     asymmetric height CNVs are larger than mirror CNVs
F     the offset itself does not scale with CNV size, within either class

Panels E and F together are the single-hit signature: a sparse hit is more likely in a
bigger CNV (E), but once a CNV is hit the displacement does not grow with how many more
genes it contains (F).
"""
from __future__ import annotations
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
ALT = 'equal_var'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.titlesize': 9, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GREEN = '#3f7d52'
COL = {'Height': BLUE, 'BMI': GREEN}

res = pickle.loads((HERE / 'outputs' / '_res.pkl').read_bytes())
tab = pd.read_csv(HERE / 'outputs' / 'bimodality_test.csv')
pick = {r['name']: r for r in res if r['mirror_subset']}

h = pd.read_csv(ROOT / 'CNV_height_BMI_corrected_directional_models_v2' /
                'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
po = pd.read_csv(HERE / 'outputs' / 'perp_offsets_per_locus.csv')
po = po[po.trait == 'Height'].tail(33)
h = h.merge(po[['CNV_ID', 'perp_offset']], on='CNV_ID')
h['cls'] = h.Dose_response_class.str.lower()
MIR = h[h.cls == 'mirror']; ASY = h[h.cls == 'asymmetric']


def dens(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def panel_letter(ax, s, dx=-0.20, dy=1.14):
    ax.text(dx, dy, s, transform=ax.transAxes, size=12, weight='bold', va='top')


fig = plt.figure(figsize=(7.5, 9.2))
gs = GridSpec(3, 2, figure=fig, left=.095, right=.975, top=.955, bottom=.055,
              hspace=.46, wspace=.26, height_ratios=[1, 1, 1.05])

# ------------------------------------------------------------------ A, B
for c, name in enumerate(['Height', 'BMI']):
    r = pick[name]
    sel = tab[(tab.trait == name) & (tab.mirror_fit == 'mirror-class only')]
    row = sel[sel.alternative == ALT].iloc[0]
    f = r['fits'][ALT]; x = r['x']; col = COL[name]
    ax = fig.add_subplot(gs[0, c])
    lo, hi = x.min() - .10, x.max() + .10
    ax.hist(x, bins=np.linspace(lo, hi, 14), density=True, color=col, alpha=.20,
            edgecolor='white', lw=.5)
    g = np.linspace(lo, hi, 500)
    ax.plot(g, dens(g, r['mu'], r['sd']), color=GREY, lw=1.4, ls='--')
    pi_, m1, s1, m2, s2 = f['par']
    ax.plot(g, pi_ * dens(g, m1, s1) + (1 - pi_) * dens(g, m2, s2), color=col, lw=1.8)
    cls = np.asarray([str(s).lower() for s in r['cls']])
    ymax = ax.get_ylim()[1]
    if (cls == 'mirror').any():
        for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
            s = cls == lab
            ax.plot(x[s], np.full(int(s.sum()), -ymax * .06), '|', color=cc, ms=6,
                    mew=1.2, clip_on=False)
    else:
        ax.plot(x, np.full(len(x), -ymax * .06), '|', color=GREY, ms=6, mew=1.2,
                clip_on=False)
    ax.set_ylim(-ymax * .11, ymax * 1.05)
    ax.set_xlabel('Offset from mirror line (SD)')
    ax.set_ylabel('Density')
    ax.set_title(f'{name}  (n = {int(row.n_loci)})', weight='bold', color=col)
    panel_letter(ax, 'AB'[c])
axs0 = fig.axes[0]
axs0.plot([], [], color=GREY, lw=1.4, ls='--', label='1 component')
axs0.plot([], [], color=BLUE, lw=1.8, label='2 components')
axs0.plot([], [], '|', color=BLUE, ms=6, mew=1.2, label='mirror')
axs0.plot([], [], '|', color=RED, ms=6, mew=1.2, label='asymmetric')
axs0.legend(frameon=False, loc='upper left', handletextpad=.5, labelspacing=.25)

# ------------------------------------------------------------------ C, D
for c, name in enumerate(['Height', 'BMI']):
    r = pick[name]
    sel = tab[(tab.trait == name) & (tab.mirror_fit == 'mirror-class only')]
    row = sel[sel.alternative == ALT].iloc[0]
    f = r['fits'][ALT]; col = COL[name]
    ax = fig.add_subplot(gs[1, c])
    nullB, obs = f['nullB'], f['obs']
    top = max(np.quantile(nullB, .995), obs) * 1.10
    ax.hist(nullB, bins=np.linspace(0, top, 40), density=True, color=GREY, alpha=.45,
            edgecolor='white', lw=.4)
    ax.axvline(obs, color=col, lw=2.2)
    ax.set_xlim(0, top)
    ax.set_xlabel('Likelihood ratio, 2 vs 1 component')
    ax.set_ylabel('Density')
    ax.set_title(f'{name}  $p$ = {row.p_bootstrap_nullB:.3f}', weight='bold', color=col)
    ax.text(.97, .93, f'observed {obs:.1f}', transform=ax.transAxes, ha='right',
            va='top', size=7, color=col)
    panel_letter(ax, 'CD'[c])

# --------------------------------------------------------------------- E
gE = GridSpecFromSubplotSpec(1, 2, gs[2, 0], wspace=.52)
rng = np.random.default_rng(3)
for i, (v, lab) in enumerate([('CNV_length_Mb', 'CNV length (Mb)'),
                              ('Gene_count', 'Genes')]):
    ax = fig.add_subplot(gE[i])
    a = ASY[v].to_numpy(float); m = MIR[v].to_numpy(float)
    tr = np.log if v == 'CNV_length_Mb' else np.log1p
    p = stats.ttest_ind(tr(a), tr(m), equal_var=False).pvalue
    for j, (vals, cc) in enumerate([(m, BLUE), (a, RED)]):
        ax.scatter(np.full(len(vals), j) + rng.normal(0, .07, len(vals)), vals, s=11,
                   color=cc, alpha=.75, lw=0)
        ax.plot([j - .26, j + .26], [np.median(vals)] * 2, color=INK, lw=1.6, zorder=4)
    ax.set_ylim(0, max(a.max(), m.max()) * 1.17)   # headroom for the p-value label
    ax.set_xlim(-.55, 1.55); ax.set_xticks([0, 1])
    ax.set_xticklabels(['mirror', 'asym.'])
    ax.set_ylabel(lab)
    ax.text(.04, .97, f'$p$ = {p:.3f}', transform=ax.transAxes, ha='left', va='top',
            size=7.4, color=RED if p < .05 else INK)
    if i == 0:
        panel_letter(ax, 'E', dx=-0.42)
fig.text(.285, .333, 'Asymmetric CNVs are larger', ha='center', size=8.8,
         weight='bold', color=INK)

# --------------------------------------------------------------------- F
gF = GridSpecFromSubplotSpec(1, 2, gs[2, 1], wspace=.40)
for i, (v, lab) in enumerate([('CNV_length_Mb', 'CNV length (Mb)'),
                              ('Gene_count', 'Genes')]):
    ax = fig.add_subplot(gF[i])
    txt = []
    for z, cc, nm in [(MIR, BLUE, 'mir'), (ASY, RED, 'asym')]:
        ax.scatter(z[v], z.perp_offset, s=11, color=cc, alpha=.8, lw=0)
        rho, p = stats.spearmanr(z[v], z.perp_offset)
        txt.append((f'{nm} $\\rho$={rho:+.2f}', cc))
    ax.set_xlabel(lab)
    if i == 0:
        ax.set_ylabel('Offset from\nmirror line (SD)')
        panel_letter(ax, 'F', dx=-0.38)
    else:
        ax.set_yticklabels([])
    ax.axhline(0, color=LIGHT, lw=.7)
    for k, (lbl, cc) in enumerate(txt):
        ax.text(.04, .97 - k * .115, lbl, transform=ax.transAxes, ha='left', va='top',
                size=6.9, color=cc)
ylo = min(MIR.perp_offset.min(), ASY.perp_offset.min()) - .05
yhi = max(MIR.perp_offset.max(), ASY.perp_offset.max()) + .14
for ax in fig.axes[-2:]:
    ax.set_ylim(ylo, yhi)
fig.text(.755, .333, 'Offset does not scale with size', ha='center', size=8.8,
         weight='bold', color=INK)

out = HERE / 'figures' / 'Figure_manuscript_bimodality_size'
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=300)
print('wrote', out.with_suffix('.pdf'))
