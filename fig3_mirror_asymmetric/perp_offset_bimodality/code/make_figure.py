#!/usr/bin/env python
"""Figure: observed perpendicular offsets from the mirror line, and the bootstrap null.

Top row    distribution of each locus's signed perpendicular distance from the fitted
           mirror line, with the fitted one-component Gaussian and the two-component
           mixture overlaid, and a rug coloured by published dose-response class.
Bottom row null distribution of the likelihood-ratio statistic under the single-Gaussian
           null, with the observed value marked.

The mirror line is calibrated on the mirror-class loci, so the asymmetric loci are
measured against a line they played no part in defining. For BMI the locus table carries
no dose-response class, so all loci are used and the two calibrations coincide.

The equal-variance mixture is the primary alternative: under the sparse single-hit
mechanism a hit displaces a locus but does not change its dispersion, so the two
components should differ in mean only. The free-variance result is annotated alongside.

--variant NAME reads outputs/NAME/ (e.g. excl_noeffect_p05, written by
perp_offset_bimodality_excl_noeffect.py) and writes Figure_perp_offset_bimodality_NAME.
"""
from __future__ import annotations
import argparse, pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent.parent
PRIMARY = 'equal_var'
MIRROR_FIT = 'mirror-class only'

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 9,
    'axes.titlesize': 10, 'axes.labelsize': 9, 'axes.linewidth': .7,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7.5,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'
BLUE = '#1f5c8b'; RED = '#b03a2e'; GREEN = '#3f7d52'
COL = {'Height': BLUE, 'BMI': GREEN}

ap = argparse.ArgumentParser()
ap.add_argument('--variant', default='')
VARIANT = ap.parse_args().variant
OUT = HERE / 'outputs' / VARIANT if VARIANT else HERE / 'outputs'
VARIANT_NOTE = {'excl_noeffect_p05': 'Loci where neither DEL nor DUP reaches p < 0.05 are excluded before fitting.'}

res = pickle.loads((OUT / '_res.pkl').read_bytes())
tab = pd.read_csv(OUT / 'bimodality_test.csv')
# one result object per trait, using the mirror-class calibration
pick = {r['name']: r for r in res if r['mirror_subset']}


def dens(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


shapes = {}
fig, axs = plt.subplots(2, 2, figsize=(8.6, 7.0))
fig.subplots_adjust(left=.085, right=.975, top=.885, bottom=.145, hspace=.52, wspace=.235)

for c, name in enumerate(['Height', 'BMI']):
    r = pick[name]
    sel = tab[(tab.trait == name) & (tab.mirror_fit == MIRROR_FIT)]
    row = sel[sel.alternative == PRIMARY].iloc[0]
    rowf = sel[sel.alternative == 'free_var'].iloc[0]
    f = r['fits'][PRIMARY]
    x = r['x']
    col = COL[name]

    # ---------------------------------------------- top: observed offsets
    ax = axs[0, c]
    lo, hi = x.min() - .12, x.max() + .12
    ax.hist(x, bins=np.linspace(lo, hi, 15), density=True, color=col, alpha=.22,
            edgecolor='white', lw=.6)
    g = np.linspace(lo, hi, 600)
    ax.plot(g, dens(g, r['mu'], r['sd']), color=GREY, lw=1.6, ls='--',
            label='1 component (null)')
    pi_, m1, s1, m2, s2 = f['par']
    ax.plot(g, pi_ * dens(g, m1, s1) + (1 - pi_) * dens(g, m2, s2), color=col, lw=1.9,
            label='2 components (equal var.)')
    ax.plot(g, pi_ * dens(g, m1, s1), color=col, lw=.9, ls=':')
    ax.plot(g, (1 - pi_) * dens(g, m2, s2), color=col, lw=.9, ls=':')

    cls = np.asarray([str(s).lower() for s in r['cls']])
    ymax = ax.get_ylim()[1]
    for lab, cc in [('mirror', BLUE), ('asymmetric', RED)]:
        s = cls == lab
        if s.any():
            ax.plot(x[s], np.full(s.sum(), -ymax * .055), '|', color=cc, ms=7, mew=1.3,
                    clip_on=False, label=f'{lab} ({int(s.sum())})')
    oth = ~np.isin(cls, ['mirror', 'asymmetric'])
    if oth.any():
        ax.plot(x[oth], np.full(int(oth.sum()), -ymax * .055), '|', color=GREY, ms=7,
                mew=1.3, clip_on=False, label=f'all loci ({int(oth.sum())})')
    ax.axvline(0, color=LIGHT, lw=.8)
    ax.set_ylim(-ymax * .10, ymax * 1.04)
    ax.set_xlabel('Perpendicular offset from mirror line (SD)')
    ax.set_ylabel('Density')
    shape = 'bimodal' if row.p_bootstrap_nullA < .05 else 'unimodal'
    shapes[name] = shape
    ax.set_title(f'{name}: offsets are {shape}\n'
                 f'n = {int(row.n_loci)}, mirror line rho = {row.rho_full_fit:.3f}',
                 weight='bold', color=col)
    ax.legend(frameon=False, loc='upper left', handletextpad=.6)
    ax.text(.985, .97, f'minor component {row.minor_weight:.2f}\n'
                       f'separation {row.separation:.3f} SD',
            transform=ax.transAxes, ha='right', va='top', size=7.4, color=INK)

    # --------------------------------------------- bottom: bootstrap null
    ax = axs[1, c]
    nullA, nullB, obs = f['nullA'], f['nullB'], f['obs']
    top = max(np.quantile(nullA, .995), np.quantile(nullB, .995), obs) * 1.10
    bins = np.linspace(0, top, 46)
    ax.hist(nullA, bins=bins, density=True, color=GREY, alpha=.40, edgecolor='white',
            lw=.4, label='null A: single Gaussian')
    ax.hist(nullB, bins=bins, density=True, histtype='step', color=INK, lw=1.0,
            label='null B: + per-locus SE')
    ax.axvline(obs, color=col, lw=2.2, label=f'observed LRT = {obs:.2f}')
    ax.set_xlim(0, top)
    ax.set_xlabel('Likelihood-ratio statistic, 2 vs 1 component')
    ax.set_ylabel('Density')
    ax.set_title(f'{name}: {int(row.n_boot)} draws from the unimodal null', weight='bold')
    ax.legend(frameon=False, loc='upper right')
    sig = row.p_bootstrap_nullA < .05
    ax.text(.035, .62,
            f'equal var.  $p$ = {row.p_bootstrap_nullA:.3f} (A),  '
            f'{row.p_bootstrap_nullB:.3f} (B)\n'
            f'free var.   $p$ = {rowf.p_bootstrap_nullA:.3f} (A),  '
            f'{rowf.p_bootstrap_nullB:.3f} (B)',
            transform=ax.transAxes, ha='left', va='top', size=7.8,
            color=RED if sig else INK,
            bbox=dict(boxstyle='round,pad=.35', fc='white', ec=RED if sig else LIGHT, lw=.8))

fig.suptitle(f'Perpendicular offset from the mirror line: {shapes["Height"]} for height, '
             f'{shapes["BMI"]} for BMI', size=11.5, weight='bold', y=.975)
hcls = np.asarray([str(s).lower() for s in pick['Height']['cls']])
caption = (f'Mirror line calibrated on the {int((hcls == "mirror").sum())} mirror-class height loci, '
           f'so the {int((hcls == "asymmetric").sum())} asymmetric loci are '
           'held out of the fit that defines the line.\nThe BMI locus table carries no dose-response '
           f'class, so all {len(pick["BMI"]["x"])} loci are used and the two calibrations coincide. '
           'Rug shows individual loci.\nNull A simulates from the fitted single Gaussian; null B adds each locus'
           "'s own measurement error. 5,000 replicates, identical EM fitting throughout.")
if VARIANT:
    caption += '\n' + VARIANT_NOTE.get(VARIANT, f'Variant: {VARIANT}.')
fig.text(.5, .075, caption, ha='center', va='top', size=7.1, color=GREY, linespacing=1.6)
out = HERE / 'figures' / ('Figure_perp_offset_bimodality' + (f'_{VARIANT}' if VARIANT else ''))
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=300)
print('wrote', out.with_suffix('.pdf'))
