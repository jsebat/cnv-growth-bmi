#!/usr/bin/env python
"""Single combined figure.

A, B  perpendicular offsets from the mirror line, 1- vs 2-component fits
C, D  observed likelihood ratio against the bootstrap null
E     asymmetric height CNVs are larger than mirror CNVs
F     the offset itself does not scale with CNV size, within either class
G, H  height effect size vs gene count, pooled and stratified by class
I, J  BMI effect size vs gene count; BMI carries no class label, so only the pooled
      estimate exists, and there is no dose response to recover
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
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
ALT = 'equal_var'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
    'axes.titlesize': 8.8, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GREEN = '#3f7d52'
COL = {'Height': BLUE, 'BMI': GREEN}

res = pickle.loads((HERE / 'outputs' / '_res.pkl').read_bytes())
tab = pd.read_csv(HERE / 'outputs' / 'bimodality_test.csv')
pick = {r['name']: r for r in res if r['mirror_subset']}

IN = ROOT / 'CNV_height_BMI_corrected_directional_models_v2'
h = pd.read_csv(IN / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
po = pd.read_csv(HERE / 'outputs' / 'perp_offsets_per_locus.csv')
po = po[po.trait == 'Height'].tail(33)
h = h.merge(po[['CNV_ID', 'perp_offset']], on='CNV_ID')
h['cls'] = h.Dose_response_class.str.lower()
MIR = h[h.cls == 'mirror']; ASY = h[h.cls == 'asymmetric']
b = pd.read_csv(IN / 'bmi' / 'inputs' / 'sacks_bmi_complete_reciprocal_loci_exact.csv')
GRP = [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]
EFFS = [('Deletion_beta', 'Deletion effect (SD)'), ('Duplication_beta', 'Duplication effect (SD)')]


def dens(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def tag(ax, s, dx=-0.21, dy=1.09):
    ax.text(dx, dy, s, transform=ax.transAxes, size=11.5, weight='bold', va='top')


fig = plt.figure(figsize=(7.5, 13.2))
gs = GridSpec(6, 2, figure=fig, left=.095, right=.975, top=.966, bottom=.052,
              hspace=.44, wspace=.26)

# ------------------------------------------------------------------- A, B
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
    tag(ax, 'AB'[c])
ax0 = fig.axes[0]
ax0.plot([], [], color=GREY, lw=1.4, ls='--', label='1 component')
ax0.plot([], [], color=BLUE, lw=1.8, label='2 components')
ax0.legend(frameon=False, loc='upper left', handletextpad=.5, labelspacing=.22)

# ------------------------------------------------------------------- C, D
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
    tag(ax, 'CD'[c])

# ---------------------------------------------------------------------- E
gE = GridSpecFromSubplotSpec(1, 4, gs[2, :], wspace=.42)
rng = np.random.default_rng(3)
PANELS = [('CNV_length_Mb', 'Length (Mb)', np.log, True),
          ('Gene_count', 'Genes', np.log1p, True),
          ('Deletion_beta', 'DEL effect', None, False),
          ('Duplication_beta', 'DUP effect', None, False)]
for i, (v, lab, tr, from_zero) in enumerate(PANELS):
    ax = fig.add_subplot(gE[i])
    a = ASY[v].to_numpy(float); m = MIR[v].to_numpy(float)
    aa, mm = (tr(a), tr(m)) if tr is not None else (a, m)
    p = stats.ttest_ind(aa, mm, equal_var=False).pvalue
    for j, (vals, cc) in enumerate([(m, BLUE), (a, RED)]):
        ax.scatter(np.full(len(vals), j) + rng.normal(0, .08, len(vals)), vals, s=8,
                   color=cc, alpha=.8, lw=0)
        ax.plot([j - .30, j + .30], [np.median(vals)] * 2, color=INK, lw=1.4, zorder=4)
    span = max(a.max(), m.max()) - min(a.min(), m.min())
    if from_zero:
        ax.set_ylim(0, max(a.max(), m.max()) * 1.26)
    else:
        ax.set_ylim(min(a.min(), m.min()) - .06 * span, max(a.max(), m.max()) + .28 * span)
        ax.axhline(0, color=LIGHT, lw=.6, zorder=1)
    ax.set_xlim(-.62, 1.62); ax.set_xticks([0, 1])
    ax.set_xticklabels(['mir', 'asy'], size=6.4)
    ax.set_ylabel(lab, size=7.2)
    ax.tick_params(labelsize=6.4)
    ptxt = f'$p$ = {p:.3f}' if p >= 1e-4 else f'$p$ = {p:.0e}'
    ax.text(.5, .99, ptxt, transform=ax.transAxes, ha='center', va='top',
            size=6.6, color=RED if p < .05 else INK,
            weight='bold' if p < .01 else 'normal')
    if i == 0:
        tag(ax, 'E', dx=-0.30)

# ---------------------------------------------------------------------- F
gF = GridSpecFromSubplotSpec(1, 2, gs[3, :], wspace=.14)
ylo = h.perp_offset.min() - .05; yhi = h.perp_offset.max() + .19
for i, (v, lab) in enumerate([('CNV_length_Mb', 'CNV length (Mb)'), ('Gene_count', 'Genes')]):
    ax = fig.add_subplot(gF[i])
    for z, cc in [(MIR, BLUE), (ASY, RED)]:
        ax.scatter(z[v], z.perp_offset, s=10, color=cc, alpha=.8, lw=0)
    for k, (z, cc, nm) in enumerate([(MIR, BLUE, 'mirr'), (ASY, RED, 'asym')]):
        rho, _ = stats.spearmanr(z[v], z.perp_offset)
        ax.text(.04, .98 - k * .105, f'{nm} $\\rho$={rho:+.2f}', transform=ax.transAxes,
                ha='left', va='top', size=6.6, color=cc)
    ax.axhline(0, color=LIGHT, lw=.7)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel(lab)
    if i == 0:
        ax.set_ylabel('Offset from mirror line (SD)'); tag(ax, 'F', dx=-0.14)
    else:
        ax.set_yticklabels([])

# ------------------------------------------------------------- G, H  height
for c, (eff, elab) in enumerate(EFFS):
    ax = fig.add_subplot(gs[4, c]); v = 'Gene_count'
    for cls, col, mk in GRP:
        z = h[h.cls == cls]
        ax.scatter(z[v], z[eff], s=20, marker=mk, color=col, alpha=.85, edgecolor='white', lw=.4)
        sl, ic = stats.theilslopes(z[eff], z[v])[:2]
        xx = np.linspace(z[v].min(), z[v].max(), 40)
        rho, p = stats.spearmanr(z[eff], z[v])
        ax.plot(xx, ic + sl * xx, color=col, lw=1.8 if p < .05 else 1.0,
                ls='-' if p < .05 else (0, (4, 2)), zorder=4)
    sl, ic = stats.theilslopes(h[eff], h[v])[:2]
    xx = np.linspace(h[v].min(), h[v].max(), 40)
    rho_all, p_all = stats.spearmanr(h[eff], h[v])
    ax.plot(xx, ic + sl * xx, color=INK, lw=1.4, ls=(0, (1, 1.6)), zorder=2)
    ax.axhline(0, color=LIGHT, lw=.7)
    ax.set_xlabel('Genes'); ax.set_ylabel(elab)
    ents = [('all ', INK, rho_all, p_all)]
    for cls, col, _ in GRP:
        z = h[h.cls == cls]
        rr, pp = stats.spearmanr(z[eff], z[v])
        ents.append((cls[:4], col, rr, pp))
    y0 = .27 if eff == 'Deletion_beta' else .99
    for k, (lbl, col, rr, pp) in enumerate(ents):
        ax.text(.03, y0 - k * .088, f'{lbl}  $\\rho$={rr:+.2f}  $p$={pp:.3f}',
                transform=ax.transAxes, ha='left', va='top', size=6.8, color=col,
                weight='bold' if pp < .05 else 'normal')
    tag(ax, 'GH'[c])

# ---------------------------------------------------------------- I, J  BMI
for c, (eff, elab) in enumerate(EFFS):
    ax = fig.add_subplot(gs[5, c]); v = 'Gene_count'
    ax.scatter(b[v], b[eff], s=20, marker='o', color=GREEN, alpha=.85, edgecolor='white', lw=.4)
    sl, ic = stats.theilslopes(b[eff], b[v])[:2]
    xx = np.linspace(b[v].min(), b[v].max(), 40)
    rr, pp = stats.spearmanr(b[eff], b[v])
    ax.plot(xx, ic + sl * xx, color=INK, lw=1.4, ls=(0, (1, 1.6)), zorder=2)
    ax.axhline(0, color=LIGHT, lw=.7)
    ax.set_xlabel('Genes'); ax.set_ylabel(elab.replace('effect', 'effect, BMI'))
    ax.text(.03, .99, f'all   $\\rho$={rr:+.2f}  $p$={pp:.3f}', transform=ax.transAxes,
            ha='left', va='top', size=6.8, color=INK)
    ax.text(.03, .90, 'no class label to stratify', transform=ax.transAxes,
            ha='left', va='top', size=6.4, color=GREY, style='italic')
    tag(ax, 'IJ'[c])

# -------------------------------------------------------------- headings
fig.text(.145, .660, 'Classes differ by size and by DUP effect, not DEL', ha='left',
         size=8.8, weight='bold', color=INK)
fig.text(.145, .499, 'Offset does not scale with CNV size, within either class', ha='left',
         size=8.8, weight='bold', color=INK)
fig.text(.145, .339, 'Height: pooling the classes masks the gene-count dose response', ha='left',
         size=8.8, weight='bold', color=INK)
fig.text(.145, .178, 'BMI: no dose response, and no classes to pool', ha='left',
         size=8.8, weight='bold', color=INK)
fig.legend(handles=[
    Line2D([], [], marker='o', color=BLUE, lw=0, ms=6, label='mirror (21)'),
    Line2D([], [], marker='s', color=RED, lw=0, ms=6, label='asymmetric (12)'),
    Line2D([], [], marker='o', color=GREEN, lw=0, ms=6, label='BMI loci (32)'),
    Line2D([], [], color=INK, lw=1.4, ls=(0, (1, 1.6)), label='all loci pooled')],
    frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(.53, .010))

out = HERE / 'figures' / 'Figure_combined'
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=200)
print('wrote', out.with_suffix('.pdf'))
