#!/usr/bin/env python
"""Build the final publication-quality PDF report.

Composes methods, results, figures and conclusions into a single vector PDF.
Requires that run_corrected_analysis.py and run_axis_mixture_models.py have both
been run, since it reads their output tables and simulation clouds.

Usage:
  python code/build_final_report.py --project-dir .
"""
from __future__ import annotations
import argparse, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import main_figure as MF

ap = argparse.ArgumentParser()
ap.add_argument('--project-dir', default=str(Path(__file__).resolve().parent.parent))
ap.add_argument('--out', default=None)
args = ap.parse_args()
P = Path(args.project_dir).resolve()
OUT = Path(args.out) if args.out else P / 'report' / 'CNV_asymmetry_final_report.pdf'
OUT.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8.5,
    'axes.titlesize': 9.5, 'axes.labelsize': 8.5, 'axes.linewidth': .7,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
    'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 200,
    'pdf.fonttype': 42,
})
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#d8d8d8'
RED = '#b03a2e'; BLUE = '#1f5c8b'; GREEN = '#3f7d52'; PURPLE = '#6b4c93'; ORANGE = '#c47b1a'
PAGE = (8.5, 11)

TEAL = '#1b7f79'; PINK = '#a3417a'
SHORT = {
    'Mirror alone': 'Mirror alone',
    'Mirror + continuous directional shift': 'Continuous\n(diagonal)',
    'Mirror + sparse single hit': 'Sparse single hit\n(diagonal)',
    'Mirror + continuous, axis-aligned': 'Continuous\n(axis-aligned)',
    'Mirror + sparse single hit, axis-aligned': 'Sparse single hit\n(axis-aligned)',
    'Mirror + continuous buffering': 'Continuous\n(buffering)',
    'Mirror + sparse single hit with buffering': 'Sparse single hit\n(buffering)',
    'Mirror + buffering, saturating occupancy': 'Buffering,\nsaturating occ.',
    'Mirror + aggregate directional shift (q dropped)': 'Aggregate drift\n(2 par)',
}
ORDER = list(SHORT.keys())
COLS = {ORDER[0]: GREY, ORDER[1]: BLUE, ORDER[2]: PURPLE, ORDER[3]: GREEN,
        ORDER[4]: RED, ORDER[5]: TEAL, ORDER[6]: PINK, ORDER[7]: '#7a5230',
        ORDER[8]: ORANGE}

# ------------------------------------------------------------------- load data
def load(trait, loci):
    d = pd.read_csv(P / trait / 'inputs' / loci)
    base = pd.read_csv(P / trait / 'outputs' / f'{trait}_mirror_parameters.csv').iloc[0].to_dict()
    cmp_ = pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_model_comparison.csv')
    bur = pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_burden_predictions.csv')
    ppc = pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_posterior_predictive_checks.csv')
    clouds = np.load(P / trait / 'outputs' / f'{trait}_ppc_clouds.npz')
    return dict(dat=d, base=base, cmp=cmp_, bur=bur, ppc=ppc, clouds=clouds)

H = load('height', 'sacks_height_all_loci_exact.csv')
B = load('bmi', 'sacks_bmi_complete_reciprocal_loci_exact.csv')
H['sign'] = -1; B['sign'] = +1
H['name'] = 'Height'; B['name'] = 'BMI'
H['unit'] = 'cm/Mb'; B['unit'] = 'kg/m$^2$/Mb'

def mirror_line(base, gmin, gmax):
    g = np.linspace(gmin, gmax, 60)
    d0 = base['p'] * base['sign_mean'] * base['m']
    return g * d0, -base['rho'] * g * d0

def residual_z(T):
    d = T['dat']; b = T['base']; s = T['sign']
    G = d.Gene_count.to_numpy(float)
    d0 = b['p'] * b['sign_mean'] * b['m']
    md = G * d0; mu = -b['rho'] * G * d0
    zD = s * (d.Deletion_beta.to_numpy(float) - md) / d.Deletion_SE.to_numpy(float)
    zU = s * (d.Duplication_beta.to_numpy(float) - mu) / d.Duplication_SE.to_numpy(float)
    return zD, zU

# --------------------------------------------------------------- text pages
def text_page(pdf, blocks, footer=None):
    fig = plt.figure(figsize=PAGE)
    y = 0.945
    for kind, txt in blocks:
        if kind == 'title':
            fig.text(.5, y, txt, ha='center', va='top', size=17, weight='bold', color=INK)
            y -= .050
        elif kind == 'subtitle':
            fig.text(.5, y, txt, ha='center', va='top', size=10.5, style='italic', color='#444')
            y -= .038
        elif kind == 'h1':
            y -= .012
            fig.text(.09, y, txt, ha='left', va='top', size=12, weight='bold', color=INK)
            y -= .028
        elif kind == 'h2':
            y -= .006
            fig.text(.09, y, txt, ha='left', va='top', size=9.8, weight='bold', color=INK)
            y -= .022
        elif kind == 'body':
            for line in textwrap.wrap(txt, 104):
                fig.text(.09, y, line, ha='left', va='top', size=8.8, color=INK)
                y -= .0165
            y -= .010
        elif kind == 'bullet':
            wrapped = textwrap.wrap(txt, 99)
            for i, line in enumerate(wrapped):
                fig.text(.115 if i else .098, y, ('•  ' + line) if i == 0 else line,
                         ha='left', va='top', size=8.8, color=INK)
                y -= .0165
            y -= .005
        elif kind == 'gap':
            y -= .018
    if footer:
        fig.text(.5, .035, footer, ha='center', size=7.2, color=GREY)
    pdf.savefig(fig); plt.close(fig)

def fig_page(pdf, fig, caption, number):
    fig.text(.07, .045, f'Figure {number}. ' + caption, ha='left', va='top', size=8,
             color=INK, wrap=True)
    pdf.savefig(fig); plt.close(fig)

# ------------------------------------------------------------------ figures
def figure1():
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.4))
    fig.subplots_adjust(left=.09, right=.97, top=.88, bottom=.28, wspace=.28)
    for ax, T in zip(axs, [H, B]):
        d = T['dat']
        has_cls = 'Dose_response_class' in d
        if has_cls:
            for cls, col, mk in [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]:
                z = d[d.Dose_response_class.str.lower() == cls]
                ax.errorbar(z.Deletion_beta, z.Duplication_beta, xerr=z.Deletion_SE,
                            yerr=z.Duplication_SE, fmt=mk, ms=4.2, color=col,
                            ecolor=LIGHT, elinewidth=.7, capsize=0, label=f'{cls} ({len(z)})',
                            lw=0, mew=.5, mec='white')
        else:
            ax.errorbar(d.Deletion_beta, d.Duplication_beta, xerr=d.Deletion_SE,
                        yerr=d.Duplication_SE, fmt='o', ms=4.2, color=BLUE, ecolor=LIGHT,
                        elinewidth=.7, capsize=0, label=f'reciprocal loci ({len(d)})',
                        lw=0, mew=.5, mec='white')
        xl, yl = mirror_line(T['base'], d.Gene_count.min(), d.Gene_count.max())
        ax.plot(xl, yl, color=INK, lw=1.6, label='fitted mirror model')
        ax.axhline(0, color=LIGHT, lw=.7); ax.axvline(0, color=LIGHT, lw=.7)
        arrow_to = (-.18, -.18) if T['sign'] < 0 else (.18, .18)
        ax.add_patch(FancyArrowPatch((0, 0), arrow_to, arrowstyle='->', mutation_scale=13,
                                     color=RED, lw=1.6, alpha=.85))
        ax.text(arrow_to[0], arrow_to[1] * 1.25,
                'lower-left\nasymmetry' if T['sign'] < 0 else 'upper-right\nasymmetry',
                color=RED, size=7.4, ha='center',
                va='top' if T['sign'] < 0 else 'bottom')
        ax.set_xlabel(f"Deletion effect on {T['name'].lower()} (SD)")
        ax.set_ylabel(f"Duplication effect on {T['name'].lower()} (SD)")
        ax.set_title(f"{T['name']}: observed reciprocal CNV effects", weight='bold')
        ax.legend(frameon=False, loc='upper right', handletextpad=.5)
    return fig

def figure2():
    fig, axs = plt.subplots(2, 2, figsize=(8.5, 6.6))
    fig.subplots_adjust(left=.09, right=.97, top=.93, bottom=.20, hspace=.42, wspace=.28)
    for col, T in enumerate([H, B]):
        zD, zU = residual_z(T)
        d = T['dat']
        ax = axs[0, col]
        has_cls = 'Dose_response_class' in d
        if has_cls:
            for cls, c, mk in [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]:
                s = d.Dose_response_class.str.lower() == cls
                ax.scatter(zD[s.to_numpy()], zU[s.to_numpy()], s=26, color=c, marker=mk,
                           edgecolor='white', lw=.5, label=cls, zorder=3)
        else:
            ax.scatter(zD, zU, s=26, color=BLUE, edgecolor='white', lw=.5,
                       label='reciprocal loci', zorder=3)
        lim = max(np.abs(np.r_[zD, zU])) * 1.15
        ax.plot([-lim, lim], [0, 0], color=LIGHT, lw=.8)
        ax.plot([0, 0], [-lim, lim], color=LIGHT, lw=.8)
        ax.plot([0, lim], [0, lim], color=GREY, lw=.8, ls='--', label='45$\\degree$ diagonal')
        cm = T['cmp'].set_index('model')
        for lab, c in [(SPARSE, PURPLE), (BUF_S, PINK)]:
            if lab in cm.index:
                bD, bU = cm.loc[lab, 'shift_DEL'], cm.loc[lab, 'shift_DUP']
                n = np.hypot(bD, bU)
                if n > 0:
                    ax.add_patch(FancyArrowPatch((0, 0), (bD / n * lim * .75, bU / n * lim * .75),
                                                 arrowstyle='->', mutation_scale=11, color=c,
                                                 lw=1.5, alpha=.9, zorder=4))
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel('Deletion residual (z, toward asymmetry)')
        ax.set_ylabel('Duplication residual (z)')
        ax.set_title(f"{T['name']}: per-locus residual geometry", weight='bold')
        ax.legend(frameon=False, loc='lower right', fontsize=6.8)

        ax = axs[1, col]
        ang = np.degrees(np.arctan2(zU, zD))
        bins = np.arange(-180, 181, 30)
        if has_cls:
            m = (d.Dose_response_class.str.lower() == 'mirror').to_numpy()
            ax.hist([ang[m], ang[~m]], bins=bins, stacked=True, color=[BLUE, RED],
                    edgecolor='white', lw=.5, label=['mirror', 'asymmetric'])
        else:
            ax.hist(ang, bins=bins, color=BLUE, edgecolor='white', lw=.5)
        for x, lab in [(0, 'DEL axis'), (90, 'DUP axis')]:
            ax.axvline(x, color=INK, lw=.9, ls=':')
            ax.text(x, 1.02, lab, size=6.6, ha='center', va='bottom', color=INK,
                    transform=ax.get_xaxis_transform())
        ax.axvline(45, color=GREY, lw=.9, ls='--')
        ax.text(45, 1.02, '45$\\degree$', size=6.6, ha='center', va='bottom', color=GREY,
                transform=ax.get_xaxis_transform())
        ax.set_xlabel('Residual displacement angle (degrees)')
        ax.set_ylabel('Loci')
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_title(f"{T['name']}: displacement direction", weight='bold')
        if has_cls: ax.legend(frameon=False, fontsize=6.8)
    return fig

def figure3():
    fig, axs = plt.subplots(2, 3, figsize=(8.5, 6.0))
    fig.subplots_adjust(left=.30, right=.97, top=.92, bottom=.22, hspace=.55, wspace=.22)
    metrics = [('delta_conditional_AIC', '$\\Delta$ conditional AIC'),
               ('delta_conditional_BIC', '$\\Delta$ conditional BIC'),
               ('delta_conditional_LOO', '$\\Delta$ conditional LOO')]
    for r, T in enumerate([H, B]):
        cm = T['cmp'].set_index('model').reindex(ORDER)
        for c, (key, lab) in enumerate(metrics):
            ax = axs[r, c]
            if key not in cm:
                ax.axis('off'); continue
            v = cm[key].to_numpy(float)
            y = np.arange(len(ORDER))[::-1]
            ax.barh(y, v, color=[COLS[m] for m in ORDER], height=.68, edgecolor='white', lw=.6)
            best = np.nanargmin(v)
            vmax = float(np.nanmax(v))
            ax.scatter([-vmax * .045], [y[best]], marker='*', s=80, color=INK, zorder=5,
                       clip_on=False)
            ax.set_yticks(y)
            ax.set_yticklabels([SHORT[m] for m in ORDER] if c == 0 else [], size=6.8)
            ax.set_xlabel(lab, size=7.8)
            ax.axvline(2, color=GREY, lw=.7, ls=':')
            if c == 0:
                ax.text(-.62, 1.16, T['name'], transform=ax.transAxes, size=11,
                        weight='bold', color=INK)
            ax.set_xlim(-vmax * .09, max(vmax * 1.16, 6))
            for yy, vv in zip(y, v):
                if np.isfinite(vv):
                    ax.text(vv + vmax * .025, yy, f'{vv:.1f}', va='center', size=6.4,
                            color=INK)
    return fig

def figure4():
    """Simulations that do and do not reproduce the observed data."""
    keys = ['mirror', 'single', 'open_single', 'aggregate']
    fig, axs = plt.subplots(2, 4, figsize=(8.5, 5.8))
    fig.subplots_adjust(left=.115, right=.98, top=.82, bottom=.26, hspace=.52, wspace=.18)
    T = H
    d = T['dat']; cl = T['clouds']
    ppc = T['ppc'].set_index(['model', 'statistic'])
    obs_corr = float(np.corrcoef(d.Deletion_beta, d.Duplication_beta)[0, 1])
    obs_dup = float(d.Duplication_beta.mean())
    for c, k in enumerate(keys):
        lab = [m for m in ORDER if SHORT[m] and m == {
            'mirror': MIRROR, 'single': SPARSE, 'open_single': BUF_S,
            'aggregate': AGG}[k]][0]
        ds = cl[f'{k}_DEL']; us = cl[f'{k}_DUP']
        ax = axs[0, c]
        ax.scatter(ds[:60].ravel(), us[:60].ravel(), s=2.2, color=COLS[lab], alpha=.10,
                   lw=0, rasterized=True)
        ax.scatter(d.Deletion_beta, d.Duplication_beta, s=16, color=INK, zorder=3,
                   edgecolor='white', lw=.4)
        ax.axhline(0, color=LIGHT, lw=.6); ax.axvline(0, color=LIGHT, lw=.6)
        ax.set_title(SHORT[lab], size=8, weight='bold', color=COLS[lab], pad=6)
        ax.set_xlabel('DEL effect (SD)', size=7.2)
        if c == 0:
            ax.set_ylabel('DUP effect (SD)', size=7.2)
        else:
            ax.set_yticklabels([])
        ax.set_xlim(-1.15, .45); ax.set_ylim(-.75, .55)

        ax = axs[1, c]
        cdf_corr = ppc.loc[(lab, 'DEL_DUP_correlation'), 'posterior_predictive_CDF']
        sim_corr = ppc.loc[(lab, 'DEL_DUP_correlation'), 'simulated_mean']
        lo = ppc.loc[(lab, 'DEL_DUP_correlation'), 'simulated_2.5%']
        hi = ppc.loc[(lab, 'DEL_DUP_correlation'), 'simulated_97.5%']
        cdf_dup = ppc.loc[(lab, 'DUP_mean'), 'posterior_predictive_CDF']
        dlo = ppc.loc[(lab, 'DUP_mean'), 'simulated_2.5%']
        dhi = ppc.loc[(lab, 'DUP_mean'), 'simulated_97.5%']
        dmean = ppc.loc[(lab, 'DUP_mean'), 'simulated_mean']
        for i, (nm, lo_, hi_, mn, ob, cdf) in enumerate(
                [('DUP mean', dlo, dhi, dmean, obs_dup, cdf_dup),
                 ('DEL-DUP corr.', lo, hi, sim_corr, obs_corr, cdf_corr)]):
            yy = 1 - i
            ax.plot([lo_, hi_], [yy, yy], color=COLS[lab], lw=4, alpha=.45,
                    solid_capstyle='butt')
            ax.plot([mn], [yy], 'o', color=COLS[lab], ms=5)
            ax.plot([ob], [yy], marker='|', color=INK, ms=13, mew=1.8)
            fail = (cdf < .025) or (cdf > .975)
            ax.text(.97, yy + .34, ('FAIL' if fail else 'pass') + f' {cdf:.3f}',
                    transform=ax.get_yaxis_transform(), va='center', ha='right', size=6.4,
                    color=RED if fail else GREEN, weight='bold' if fail else 'normal')
        ax.set_yticks([0, 1])
        if c == 0:
            ax.set_yticklabels(['DEL-DUP\ncorrelation', 'DUP mean'], size=6.8)
        else:
            ax.set_yticklabels([])
        ax.set_ylim(-.7, 1.75); ax.set_xlim(-.95, .35)
        ax.set_xlabel('statistic value', size=7.2)
        ax.axvline(0, color=LIGHT, lw=.6)
    fig.legend(handles=[
        Line2D([], [], marker='o', lw=0, color=INK, ms=5, label='observed loci'),
        Line2D([], [], marker='o', lw=0, color=GREY, ms=5, alpha=.5, label='simulated replicates'),
        Line2D([], [], color=GREY, lw=4, alpha=.45, label='95% simulated interval'),
        Line2D([], [], marker='|', color=INK, lw=0, ms=11, mew=1.8, label='observed statistic')],
        frameon=False, fontsize=7, ncol=4, loc='upper center', bbox_to_anchor=(.55, .90))
    fig.suptitle('Height: which simulations reproduce the observed data', size=11.5,
                 weight='bold', y=.975)
    return fig

def figure5():
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.3))
    fig.subplots_adjust(left=.30, right=.96, top=.88, bottom=.30, wspace=.62)
    for ax, T in zip(axs, [H, B]):
        b = T['bur'].set_index('model').reindex(ORDER)
        v = b['predicted_DUP_effect_per_Mb_at_matched_density'].to_numpy(float)
        y = np.arange(len(ORDER))[::-1]
        adverse = (v < 0) if T['sign'] < 0 else (v > 0)
        ax.barh(y, v, color=[COLS[m] for m in ORDER], height=.66, edgecolor='white', lw=.6,
                alpha=.95)
        for yy, vv, ok in zip(y, v, adverse):
            # DejaVu Serif lacks the tick/cross glyphs; DejaVu Sans has both.
            ax.text(vv + (.004 if vv >= 0 else -.004), yy, ('✓' if ok else '✗'),
                    va='center', ha='left' if vv >= 0 else 'right', size=9,
                    color=GREEN if ok else RED, weight='bold', family='DejaVu Sans')
        ax.axvline(0, color=INK, lw=.9)
        ax.set_yticks(y); ax.set_yticklabels([SHORT[m] for m in ORDER], size=6.8)
        ax.set_xlabel(f"Predicted duplication burden ({T['unit']})", size=7.8)
        ax.set_title(f"{T['name']}: burden propagation", weight='bold')
        note = ('Auwerx: greater CNV burden → lower height'
                if T['sign'] < 0 else 'Auwerx: greater CNV burden → higher adiposity')
        ax.text(.5, -.30, note + f"\ncalibrated to Macé deletion coefficient "
                f"{T['bur'].Mace_DEL_target_per_Mb.iloc[0]} {T['unit']}",
                transform=ax.transAxes, ha='center', size=6.8, color=GREY)
        pad = max(abs(v)) * .35
        ax.set_xlim(min(v) - pad, max(v) + pad)
    return fig

# --------------------------------------------------------------------- build
hc = H['cmp'].set_index('model'); bc = B['cmp'].set_index('model')
hb = H['bur'].set_index('model'); bb = B['bur'].set_index('model')
hp = H['ppc'].set_index(['model', 'statistic'])

def g(t, m, c):
    try: return float(t.loc[m, c])
    except Exception: return float('nan')

BUF_S = 'Mirror + sparse single hit with buffering'
BUF_C = 'Mirror + continuous buffering'
AX_S = 'Mirror + sparse single hit, axis-aligned'
SPARSE = 'Mirror + sparse single hit'
CONT = 'Mirror + continuous directional shift'
AGG = 'Mirror + aggregate directional shift (q dropped)'
MIRROR = 'Mirror alone'

def ANG(t, m):
    """Fitted shift angle in degrees, from the stored components."""
    return float(np.degrees(np.arctan2(g(t, m, 'shift_DUP'), g(t, m, 'shift_DEL'))))

_mace = pd.read_csv(P / 'code' / 'mace_burden_reference.csv')
mace_h = float(_mace[(_mace.trait == 'height') & (_mace.burden_model == 'duplication')].beta.iloc[0])
mace_b = float(_mace[(_mace.trait == 'bmi') & (_mace.burden_model == 'duplication')].beta.iloc[0])

with PdfPages(OUT) as pdf:
    text_page(pdf, [
        ('gap', ''), ('gap', ''),
        ('title', 'Geometry of asymmetric reciprocal CNV effects'),
        ('title', 'on height and BMI'),
        ('subtitle', 'Diagonal versus axis-aligned displacement, and what each implies '
                     'for genome-wide burden'),
        ('gap', ''),
        ('h1', 'Summary'),
        ('body',
         'Reciprocal copy-number variants displace quantitative traits in a partly mirrored '
         'fashion: deletion and duplication effects are negatively correlated, with duplication '
         'effects the weaker of the two. A mirror model alone cannot account for the observed '
         'loci, and predicts the wrong direction for genome-wide duplication burden. Directional '
         'asymmetry is required. This report asks what geometry that asymmetry takes.'),
        ('body',
         'The previously published models displace every affected locus by the same vector, with '
         'the deletion and duplication components both non-zero, so each locus moves diagonally. '
         'We add an axis-aligned alternative in which each affected locus is displaced along one '
         'axis only, deletion-side with probability w and duplication-side with probability 1-w. '
         'The aggregate cloud still drifts diagonally, but no individual locus does. We also add '
         'a reparameterised model that estimates only the aggregate per-gene drift, dropping an '
         'event-frequency parameter that the data cannot identify.'),
        ('h2', 'Principal findings'),
        ('bullet',
         f"For height the decisive constraint turned out to be the published one. Restricting the "
         f"shift to the lower-left quadrant forces the fitted angle onto its 85 degree boundary in "
         f"every model. Releasing that bound moves the estimate to "
         f"{ANG(hc, BUF_S):.0f} degrees, past 90, where the deletion effect is attenuated toward "
         f"the null while the duplication effect is driven past it. That is buffering in the sense "
         f"of Sacks et al., and it is not representable at any parameter value inside the "
         f"published window."),
        ('bullet',
         f"The sparse model with buffering takes every likelihood-based criterion, with three "
         f"parameters rather than the four the axis-aligned model needs: best conditional AIC, BIC "
         f"and leave-one-locus-out, the widest posterior predictive margin of any asymmetry model, "
         f"and a predicted duplication burden of "
         f"{g(hb, BUF_S,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f} cm per Mb against "
         f"an observed {float(mace_h):+.4f}, closer than any other model."),
        ('bullet',
         'It does not win everything. The buffered model implies that duplication effect '
         'decreases with gene count, because a large duplication displacement combined with rising '
         'occupancy overwhelms the positive trend contributed by the mirror component. The '
         'observed trend is positive and nominally significant. On that one diagnostic the '
         'axis-aligned sparse model is the best fit and buffering falls outside the observed '
         'interval. No single height model is therefore consistent with every line of evidence, '
         'and the two should be reported together.'),
        ('bullet',
         'All three lower-left-constrained models reproduce the duplication mean but break the '
         'deletion-duplication correlation, a mechanical consequence of a shared diagonal shift '
         'inflating the covariance. Buffering relieves that, because its two components have '
         'opposite sign and the induced covariance is negative rather than positive.'),
        ('bullet',
         f"For BMI the event frequency is not identified: the profile likelihood is flat and q "
         f"runs to whatever upper bound it is given. Estimating only the identified aggregate "
         f"drift uses two parameters instead of three and wins both conditional AIC and BIC, "
         f"resolving a criterion conflict present in the earlier analysis."),
        ('bullet',
         'The two traits therefore differ in mechanism. Height requires stochastic, discrete '
         'occupancy: a mean-only model is worse by more than 31 AIC units. BMI does not; a '
         'deterministic drift proportional to gene count is sufficient and preferred.'),
        ('bullet',
         'Every asymmetry model reverses the incorrect duplication-burden direction predicted by '
         'the mirror model alone. For height the axis-aligned sparse model predicts the strongest '
         'adverse duplication burden, roughly threefold that of the diagonal model.'),
    ], footer='Analysis conditional on the calibrated mirror component. '
              'Exact Sacks et al. reciprocal-CNV effects; burden compared with Macé et al. and Auwerx et al.')

    text_page(pdf, [
        ('h1', 'Methods'),
        ('h2', 'Data'),
        ('body',
         f"Single-locus analyses used {len(H['dat'])} independent reciprocal height loci "
         f"({int((H['dat'].Dose_response_class.str.lower()=='mirror').sum())} mirror, "
         f"{int((H['dat'].Dose_response_class.str.lower()=='asymmetric').sum())} asymmetric) and "
         f"{len(B['dat'])} reciprocal BMI loci, taken from the exact Sacks et al. supplementary "
         f"statistics. The overlapping 22q11.2 A-D interval was excluded. Each locus contributes a "
         f"deletion and duplication beta, their standard errors, and an annotated gene count. "
         f"Genome-wide burden coefficients are taken from Macé et al.; the overall burden "
         f"direction is taken from Auwerx et al."),
        ('h2', 'Mirror component'),
        ('body',
         'A CNV effect is modelled as the additive sum of latent per-gene effects. Kosmicki et al. '
         'informs the height sign prior (70 per cent of deletion or loss-of-function effects '
         'decrease height); Akbari et al. informs the BMI sign prior (posterior mean 78.6 per cent '
         'of predicted loss-of-function effects increase BMI). Contributor frequency, mean absolute '
         'per-gene effect, duplication-to-deletion magnitude ratio and residual dispersion were '
         'estimated from the Sacks CNVs themselves rather than imported from ascertainment-selected '
         'exome effects. All asymmetry models below are fitted conditional on this calibrated '
         'mirror component, so their information criteria are directly comparable.'),
        ('h2', 'Asymmetry geometries'),
        ('body',
         'Diagonal models displace an affected locus by b(cos theta, sin theta) with theta '
         'constrained to the interval 5 to 85 degrees, so both components are necessarily non-zero. '
         'Trait sign confines the displacement to the lower-left quadrant for height and the '
         'upper-right for BMI. Two occupancy structures are used: a continuous process in which '
         'many small effects accumulate with gene count, and a sparse process in which the presence '
         'of at least one rare gene produces a fixed displacement with probability 1-(1-q)^G.'),
        ('body',
         'The 5 to 85 degree window is itself a modelling choice, and it is the one that mattered. '
         'It forces both components of the shift to share a sign, so an affected locus must be '
         'displaced further from the null in both dosage states. Widening the window to 5 to 175 '
         'degrees keeps the duplication component pointing toward the asymmetric direction but '
         'lets the deletion component take either sign. Beyond 90 degrees the fitted displacement '
         'attenuates the deletion effect toward the null while carrying the duplication effect '
         'past it, which is buffering as defined in Sacks et al.: the locus is asymmetric because '
         'the duplication side carries the effect and the deletion side is buffered, not because '
         'both dosage states are damaging. The buffered models take the same three parameters as '
         'the models they generalise, and the published models are the restricted case.'),
        ('body',
         'Axis-aligned models replace the shared vector with a mixture. An affected locus is '
         'displaced by (b_DEL, 0) with probability w, or by (0, b_DUP) with probability 1-w. Both '
         'occupancy structures are implemented, giving four parameters (q, w, b_DEL, b_DUP) against '
         'three for the diagonal models. The discriminating structure is second-moment: a shared '
         'diagonal shift adds a positive term to the deletion-duplication covariance, whereas '
         'axis-aligned events add none, because no single event moves both coordinates.'),
        ('body',
         'The aggregate model estimates only the quantity the data identify. Because the continuous '
         'likelihood is flat in q, with the product q*b nearly constant, this model drops q and '
         'estimates the per-gene drift directly as a deterministic displacement with no added '
         'dispersion, using two parameters. Note that this collapses the diagonal versus '
         'axis-aligned distinction, which survives only in the dispersion term.'),
        ('h2', 'Model comparison and checks'),
        ('body',
         'Models are compared by conditional AIC and BIC and by conditional leave-one-locus-out '
         'log predictive density, refitting the asymmetry parameters within each fold. Posterior '
         'predictive checks simulate 1200 replicate datasets per model and compare nine summary '
         'statistics with their observed values; a posterior predictive CDF below 0.025 or above '
         '0.975 is treated as a failure. Each fitted gene-level model is propagated to expected '
         'deletion and duplication burden, calibrated to the Macé deletion coefficient by '
         'allowing an effective density of dosage-responsive genes per Mb, and the resulting '
         'duplication direction is compared with Auwerx.'),
        ('h2', 'Reproducibility'),
        ('body',
         'All random seeds are fixed. Two independent implementations were run: a Python reference '
         'pipeline and an R implementation of the core inference. Fitted parameters, locus '
         'predictions and burden predictions agree between them to within optimiser tolerance, the '
         'largest discrepancy anywhere being three parts in a thousand, in the leave-one-out sums. '
         'The axis-aligned code reproduces the published diagonal fits to fifteen decimal places.'),
    ])

    # ---- Main summary figure and Table 1 -------------------------------------
    MH, MB, mace = MF.load_all(P)
    mfig = MF.build_main_figure(MH, MB, mace)
    mfig.savefig(P / 'report' / 'Figure_main_summary.pdf')
    mfig.savefig(P / 'report' / 'Figure_main_summary.png', dpi=400)
    pdf.savefig(mfig); plt.close(mfig)

    text_page(pdf, [
        ('h1', 'Figure legends'),
        ('h2', 'Main figure'),
        ('body',
         'Main figure. Geometry of asymmetric reciprocal CNV effects and its consequences for '
         'genome-wide burden.'),
        ('bullet',
         'A. The motivating observation. Residuals of each height locus from the fitted mirror '
         'expectation, standardised by that locus’s own standard error and signed so the '
         'asymmetric direction is positive. Individual loci concentrate near the deletion and '
         'duplication axes rather than on the 45 degree diagonal that the published models assume, '
         'and the published asymmetric class label tracks the duplication axis specifically.'),
        ('bullet',
         'B, E. Reciprocal-effect plane for each model. Open symbols are the observed Sacks loci, '
         'sized by gene count; circles are mirror-class and squares asymmetric-class. Shading gives '
         'the 50 and 95 per cent predictive regions from 1200 simulated replicate datasets. Black '
         'lines are the model mean as a function of gene count, drawn as explicit mixture '
         'components with line weight proportional to component probability: one component for the '
         'mirror, continuous and aggregate models, two for the diagonal sparse model, and three for '
         'the axis-aligned sparse model, whose deletion-shifted and duplication-shifted branches '
         'are dashed and dotted.'),
        ('bullet',
         'C, F. Burden propagation expressed as a multiple of the observed Macé estimate. The '
         'deletion row is the calibration anchor and is fixed at one by construction, so it is not '
         'a test; the duplication row is a free out-of-sample prediction. Grey band is the 95 per '
         'cent confidence interval of the observed duplication coefficient. Values left of zero '
         'have the wrong sign.'),
        ('bullet',
         'D, G. Model comparison and posterior predictive checks, aligned under the model columns. '
         'Deltas are relative to the best model within each trait; posterior predictive values are '
         'CDFs, with failures below 0.025 or above 0.975 shaded red. The winning model is tinted, '
         'separately for each trait, because the winner differs between them.'),
        ('gap', ''),
        ('bullet',
         'Note on B and E. Only the mirror and the sparse models keep the classic reciprocal '
         'mean trajectory, in which duplication effect moves opposite to deletion effect as gene '
         'count rises. The continuous-type models, including the aggregate drift preferred '
         'for BMI, add enough per-gene duplication drift to reverse that trend, so their mean '
         'curve tilts the same way on both axes. Their predicted deletion-duplication correlation '
         'nonetheless remains negative, because across loci that correlation is dominated by the '
         'residual covariance rather than by the gene-count trajectory. Figure 6 separates the two.'),
        ('body',
         'Observed burden coefficients are from Macé et al. 2017 Supplementary Table 1, which '
         'reports all four burden models the study tested. Standard errors are not published and '
         'have been reconstructed from the reported coefficient and P value assuming a two-sided '
         'Wald test; the intervals shown should be read with that in mind.'),
        ('h2', 'Table 1'),
        ('body',
         'Model comparison, posterior predictive behaviour and burden prediction, by trait and '
         'model. Fitted parameters, added parameter count, deltas against the best model within '
         'each trait, posterior predictive CDFs for the duplication mean and the '
         'deletion-duplication correlation, the effective dosage-gene density each model requires '
         'to match the observed deletion burden, and its resulting out-of-sample duplication '
         'prediction. The implied affected fraction applies to the sparse models only and is the '
         'mean over loci of 1-(1-q)^G.'),
        ('h2', 'Supporting figures'),
        ('bullet',
         'Figure 1. Published reciprocal-CNV effects for both traits with standard errors, the '
         'fitted mirror model, and the direction of the residual asymmetry it leaves unexplained.'),
        ('bullet',
         'Figure 2. Per-locus displacement geometry and the distribution of displacement angles, '
         'with the caveat that mirror-model residuals are structurally anti-correlated.'),
        ('bullet',
         'Figure 3. Model comparison across all six models, three criteria and both traits.'),
        ('bullet',
         'Figure 4. Simulated replicate datasets against observed height data for four '
         'representative models, with posterior predictive pass and fail marked.'),
        ('bullet',
         'Figure 5. Burden propagation in native units against the Auwerx direction.'),
        ('bullet',
         'Figure 6. Effect versus gene count with exact model-implied mean curves, and the '
         'implied against observed gene-count trend in duplication effect.'),
    ])

    tdf = MF.table1(MH, MB, mace)
    tdf.to_csv(P / 'report' / 'Table1_model_summary.csv', index=False)
    tfig = MF.table1_page(tdf, figsize=(11, 8.5))
    tfig.savefig(P / 'report' / 'Table1_model_summary.pdf')
    pdf.savefig(tfig); plt.close(tfig)

    fig_page(pdf, figure1(),
             'Published reciprocal-CNV effects. Each point is one locus, with standard errors. '
             'The fitted mirror model (black) captures the negative reciprocal trend in both '
             'traits. Height loci are labelled by their published dose-response class. The red '
             'arrow marks the direction of the residual asymmetry that the mirror model does not '
             'explain: lower-left for height, upper-right for BMI.', 1)

    fig_page(pdf, figure2(),
             'Per-locus displacement geometry. Top: residuals from the fitted mirror expectation, '
             'standardised by each locus’s own standard error and sign-flipped so the '
             'asymmetric direction is positive on both axes, so a displacement drawn up and to the '
             'right here is lower-left in raw effect space for height. Arrows show the fitted '
             'displacement vectors of the published diagonal sparse model (purple) and the buffered '
             'sparse model (pink); the buffered vector points beyond 90 degrees, into the region '
             'the published 5 to 85 degree window excludes. Bottom: the distribution of '
             'displacement angles. Individual loci concentrate near the deletion and duplication '
             'axes rather than at 45 degrees. Note that residuals from a fitted mirror model are '
             'structurally anti-correlated, so angle alone is suggestive rather than decisive; the '
             'geometry is tested inside the likelihood in Figure 3.', 2)

    fig_page(pdf, figure3(),
             'Model comparison, conditional on the calibrated mirror component. Bars show the '
             'penalty relative to the best model on each criterion; a star marks the winner and '
             'the dotted line marks the conventional two-unit threshold. Height prefers sparse, '
             'discrete occupancy under every criterion, with the axis-aligned variant winning AIC '
             'and leave-one-out and the diagonal variant winning BIC. BMI prefers the '
             'two-parameter aggregate drift on both AIC and BIC, and the mean-only model that '
             'fails badly for height is the best model here.', 3)

    fig_page(pdf, figure4(),
             'Simulations that do and do not reproduce the observed height data. Top: 60 simulated '
             'replicate datasets (pale) under each fitted model, with the observed loci overlaid '
             '(black). Bottom: two posterior predictive statistics, showing the 95 per cent '
             'simulated interval, the simulated mean, and the observed value, with the posterior '
             'predictive CDF. The mirror model reproduces the deletion-duplication correlation but '
             'fails the duplication mean outright. The published diagonal sparse model corrects the '
             'mean and breaks the correlation, because a same-signed shift adds a positive term to '
             'the modelled covariance. The aggregate drift model fails the mean. The buffered '
             'sparse model fails neither and has the widest margin on both statistics, since its '
             'two components carry opposite signs and the covariance it induces is negative; even '
             'so its simulated correlation still falls short of the observed value.', 4)

    gfig = MF.gene_count_figure(MH, MB)
    gfig.savefig(P / 'report' / 'Figure_gene_count_diagnostic.pdf')
    fig_page(pdf, gfig,
             'Effect versus gene count, the trajectory the reciprocal-plane panels and the '
             'posterior predictive statistics are blind to. Top and middle: observed deletion and '
             'duplication effects against annotated gene count, with each model’s exact '
             'E[effect | gene count] curve. These are evaluated rather than linearised, because '
             'the sparse models enter through the occupancy 1-(1-q)^G, which is non-linear. '
             'Bottom: the inverse-variance weighted slope of duplication effect on gene count '
             'implied by each model, against the observed slope and its 95 per cent interval; '
             'filled markers fall inside the interval, open markers outside. For height the '
             'axis-aligned sparse model lands almost exactly on the observed trend, while mirror '
             'alone is roughly fourfold too steep. For BMI, where every other diagnostic is '
             'ambiguous, this is the one comparison that rejects mirror alone: it is the only '
             'model whose implied slope falls outside the observed interval.', 6)

    fig_page(pdf, figure5(),
             'Propagation to genome-wide burden. Each fitted gene-level model is calibrated to the '
             'published Macé deletion coefficient by allowing an effective density of '
             'dosage-responsive genes per Mb; bars show the duplication burden predicted at that '
             'same density. A tick marks agreement with the adverse genome-wide direction reported '
             'by Auwerx, a cross marks disagreement. The mirror model alone predicts a beneficial '
             'duplication burden for height and a protective one for BMI, in both cases the wrong '
             'sign. Every asymmetry model corrects it; for height the axis-aligned sparse model '
             'predicts the largest adverse effect.', 5)

    text_page(pdf, [
        ('h1', 'Results'),
        ('h2', 'Height: single-locus fit'),
        ('body',
         f"The mirror model alone is decisively rejected, trailing the best model by "
         f"{g(hc, MIRROR,'delta_conditional_AIC'):.1f} AIC units. Among the models restricted to "
         f"the lower-left quadrant, the axis-aligned sparse mixture is best on conditional AIC and "
         f"leave-one-out, with a deletion-axis share of w = {g(hc, AX_S,'w_deletion_share'):.3f}, "
         f"close to an even split. But every one of those models sits on the 85 degree boundary of "
         f"the permitted angular window, which is the signature of a constraint doing the work "
         f"rather than the data."),
        ('body',
         f"Releasing the constraint resolves it. The sparse model with buffering settles at "
         f"{ANG(hc, BUF_S):.1f} degrees, with a deletion component of "
         f"{g(hc, BUF_S,'shift_DEL'):+.4f} SD and a duplication component of "
         f"{g(hc, BUF_S,'shift_DUP'):+.4f} SD. Applied at the median gene count this attenuates "
         f"the expected deletion effect from about -0.28 to -0.08 SD while carrying the expected "
         f"duplication effect from about +0.10 to -0.25 SD: the deletion side is buffered and the "
         f"duplication side becomes the damaging one. It takes conditional AIC, BIC and "
         f"leave-one-out simultaneously, ahead of the axis-aligned model by "
         f"{g(hc, AX_S,'delta_conditional_AIC'):.2f}, {g(hc, AX_S,'delta_conditional_BIC'):.2f} "
         f"and {g(hc, AX_S,'delta_conditional_LOO'):.2f} units, and ahead of the published diagonal "
         f"sparse model by {g(hc, SPARSE,'delta_conditional_AIC'):.2f}, "
         f"{g(hc, SPARSE,'delta_conditional_BIC'):.2f} and "
         f"{g(hc, SPARSE,'delta_conditional_LOO'):.2f}. It does so with three parameters, one "
         f"fewer than the axis-aligned model. The continuous buffered variant is second on every "
         f"criterion, so the result is driven by the geometry rather than by the occupancy "
         f"structure."),
        ('body',
         f"The aggregate model, which removes the dispersion contributed by stochastic occupancy, "
         f"is worse by {g(hc, AGG,'delta_conditional_AIC'):.1f} AIC units and "
         f"{g(hc, AGG,'delta_conditional_LOO'):.1f} units of predictive density, and its "
         f"deletion component collapses to zero. Height asymmetry is therefore not a smooth drift; "
         f"the discrete occupancy structure carries real information."),
        ('h2', 'Height: posterior predictive behaviour'),
        ('body',
         f"The mirror model reproduces the observed deletion-duplication correlation of "
         f"{float(np.corrcoef(H['dat'].Deletion_beta, H['dat'].Duplication_beta)[0,1]):.3f} "
         f"(CDF {g(hp.xs('DEL_DUP_correlation', level=1), ORDER[0], 'posterior_predictive_CDF'):.3f}) "
         f"but fails the duplication mean outright "
         f"(CDF {g(hp.xs('DUP_mean', level=1), ORDER[0], 'posterior_predictive_CDF'):.3f}). "
         f"Both published asymmetry models correct the mean and then fail the correlation "
         f"(CDF {g(hp.xs('DEL_DUP_correlation', level=1), ORDER[2], 'posterior_predictive_CDF'):.3f} "
         f"for the diagonal sparse model), which follows mechanically from a shared diagonal shift "
         f"adding a positive term to the modelled covariance."),
        ('body',
         f"Buffering relieves this. Its two components carry opposite signs, so the covariance it "
         f"induces is negative rather than positive, working with the observed reciprocal "
         f"relationship instead of against it. The buffered sparse model reaches a duplication "
         f"mean CDF of "
         f"{g(hp.xs('DUP_mean', level=1), BUF_S, 'posterior_predictive_CDF'):.3f} and a "
         f"correlation CDF of "
         f"{g(hp.xs('DEL_DUP_correlation', level=1), BUF_S, 'posterior_predictive_CDF'):.3f}, the "
         f"widest margin on both statistics of any asymmetry model. The axis-aligned model clears "
         f"the same thresholds but only narrowly, at "
         f"{g(hp.xs('DUP_mean', level=1), AX_S, 'posterior_predictive_CDF'):.3f} and "
         f"{g(hp.xs('DEL_DUP_correlation', level=1), AX_S, 'posterior_predictive_CDF'):.3f}. Even "
         f"so, the buffered simulated correlation of "
         f"{g(hp.xs('DEL_DUP_correlation', level=1), BUF_S, 'simulated_mean'):+.3f} still falls "
         f"short of the observed value, so the correct statement is that buffering reduces the "
         f"discrepancy that a same-signed shift creates rather than removing it."),
        ('h2', 'BMI: identifiability and model choice'),
        ('body',
         'The continuous model is not identified for BMI. Along the profile, twice the log '
         'likelihood varies by only 0.34 units between q of 0.22 and q of 0.80, against a '
         'conventional threshold of 3.84, while the product of frequency and magnitude stays '
         'essentially constant. Raising the upper bound from 0.80 to 0.99 moves the estimate to '
         'the new bound and improves AIC by 0.03 units, leaving the model ranking unchanged. Only '
         'the aggregate drift is estimable.'),
        ('body',
         f"Estimating that quantity directly, with two parameters rather than three, gives the best "
         f"conditional AIC and the best conditional BIC, ahead of the continuous model by "
         f"{g(bc, CONT,'delta_conditional_AIC'):.2f} and "
         f"{g(bc, CONT,'delta_conditional_BIC'):.2f} units and ahead of the mirror model by "
         f"{g(bc, MIRROR,'delta_conditional_AIC'):.2f} and "
         f"{g(bc, MIRROR,'delta_conditional_BIC'):.2f}. Its leave-one-out score is within "
         f"{g(bc, AGG,'delta_conditional_LOO'):.2f} units of the best. The earlier conflict, "
         f"in which AIC and leave-one-out favoured the continuous model while BIC favoured mirror "
         f"alone, does not survive removing the unidentified parameter. The fitted per-gene drift "
         f"is {g(bc, AGG,'shift_DEL'):.4f} SD on the deletion axis and "
         f"{g(bc, AGG,'shift_DUP'):.4f} SD on the duplication axis, an upper-right direction "
         f"consistent with the earlier estimate."),
        ('body',
         'BMI posterior predictive checks discriminate nothing: every model, including mirror '
         'alone, lies inside the interval on all nine statistics. Nor does BMI provide an '
         'independent asymmetric class label, a size enrichment, or an occupancy result analogous '
         'to height. The diagonal and axis-aligned geometries differ by two AIC units but by only '
         '0.12 units of predictive density, so the geometry of the BMI shift remains undetermined.'),
        ('h2', 'Gene-count trajectory'),
        ('body',
         'The across-locus deletion-duplication correlation is dominated by the residual '
         'covariance and is therefore nearly uninformative: all six BMI models simulate a '
         'correlation between -0.56 and -0.59 against an observed -0.73. The models disagree '
         'instead about the mean trajectory in gene count, which no criterion above examines. '
         'Mirror and the two sparse models retain a negative duplication trend with gene count; '
         'all three continuous-type models reverse it, because the per-gene duplication drift '
         'exceeds the mirror component it has to overcome.'),
        ('body',
         'Compared against the observed inverse-variance weighted trend, this separates the models '
         'where other criteria cannot, and it is the one place the buffered model does badly. For '
         'height the observed slope of duplication effect on gene count is +0.00092, with a 95 per '
         'cent interval of +0.00008 to +0.00176. The axis-aligned sparse model implies +0.00108, '
         'close to the centre of that interval. Mirror alone implies +0.00391, roughly fourfold too '
         'steep. The buffered sparse model implies -0.00049, outside the interval and of the wrong '
         'sign: its large duplication displacement, combined with occupancy rising in gene count, '
         'overwhelms the positive trend the mirror component contributes.'),
        ('body',
         'This is a genuine conflict between lines of evidence rather than a defect in one of them. '
         'The buffered model is preferred by every likelihood-based criterion, by both posterior '
         'predictive statistics and by the burden prediction, while the axis-aligned model is the '
         'only height model whose gene-count trajectory matches the data and which also passes both '
         'posterior predictive checks. Both should be reported.'),
        ('body',
         'The obvious reconciliation was fitted and does not work. A displacement that is buffered '
         'on the deletion side but smaller on the duplication side is not a separate model, because '
         'the buffered parameterisation already estimates the two components freely; the only '
         'remaining lever is the occupancy function. Replacing 1-(1-q)^G with a saturating form, in '
         'which the effective number of eligible genes approaches a ceiling, adds one parameter and '
         'moves the implied gene-count slope from -0.00049 to -0.00002. That is very nearly the '
         'observed interval but still outside it, and it is bought at a cost: conditional AIC is '
         '1.75 units worse than plain buffering, BIC 3.25 worse, and the predicted duplication '
         'burden overshoots to -0.1054 against an observed -0.0784, where plain buffering was the '
         'closest of any model. The estimated saturation point, about 78 genes, lies beyond the '
         'largest observed CNV at 58 genes, so the data cannot locate it. For BMI the saturation '
         'parameter is unbounded and the model collapses back to plain buffering with a wasted '
         'parameter.'),
        ('body',
         'The conflict is therefore real and unresolved by this route. It is not that the buffered '
         'geometry is nearly right and needs a small correction to the occupancy; the gene-count '
         'trajectory and the reciprocal-plane likelihood are pulling in genuinely different '
         'directions at this sample size.'),
        ('body',
         'For BMI the observed slope is +0.00022 with an interval of -0.00064 to +0.00107, and '
         'mirror alone, implying -0.00160, is the only model of the eight that falls outside it. '
         'This is the one diagnostic in the whole analysis that rejects mirror alone for BMI, and '
         'it does so on a quantity the posterior predictive checks do not capture. Opening the '
         'angular bound changes nothing for BMI: the buffered variants converge to exactly the '
         'constrained solutions, because the BMI optimum is interior at about 70 degrees and the '
         'bound was never binding.'),
    ])

    text_page(pdf, [
        ('h1', 'Results, continued'),
        ('h2', 'Genome-wide burden'),
        ('body',
         f"Calibrated to the Macé deletion coefficient of "
         f"{hb.Mace_DEL_target_per_Mb.iloc[0]} cm per deleted Mb, the mirror model predicts that "
         f"duplication burden increases height by "
         f"{g(hb, MIRROR,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f} cm per Mb, "
         f"which contradicts the adverse genome-wide direction reported by Auwerx. Every asymmetry "
         f"model reverses the sign. The diagonal sparse model predicts "
         f"{g(hb, SPARSE,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f} and the "
         f"axis-aligned sparse model "
         f"{g(hb, AX_S,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}, roughly "
         f"threefold larger. The aggregate model, already rejected on fit for height, also fails "
         f"here, predicting "
         f"{g(hb, AGG,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}."),
        ('body',
         f"For BMI, calibrated to {bb.Mace_DEL_target_per_Mb.iloc[0]} kg/m2 per deleted Mb, the "
         f"mirror model predicts a protective duplication burden of "
         f"{g(bb, MIRROR,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}. All asymmetry "
         f"models move it to zero or positive; the preferred aggregate model predicts "
         f"{g(bb, AGG,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}, in the adverse "
         f"direction reported by Auwerx. Burden therefore rejects the mirror model in both traits."),
        ('body',
         f"Because Macé et al. report a duplication-burden coefficient as well as a deletion one, "
         f"the duplication comparison is quantitative rather than merely directional, and it does "
         f"discriminate among the asymmetry models for height. Against an observed "
         f"{float(mace_h):+.4f} cm per Mb, the buffered sparse model predicts "
         f"{g(hb, BUF_S,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}, the axis-aligned "
         f"model {g(hb, AX_S,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f} and the "
         f"published diagonal sparse model "
         f"{g(hb, SPARSE,'predicted_DUP_effect_per_Mb_at_matched_density'):+.4f}, the last being "
         f"roughly fourfold too small. For BMI the observed duplication coefficient is not "
         f"significant, so it remains a weak constraint there, and all asymmetry models are "
         f"correctly signed and undershoot."),
        ('h2', 'Independent evidence for discrete occupancy in height'),
        ('body',
         'Asymmetric-class height CNVs are longer than mirror-class CNVs (1.06 against 0.65 Mb, '
         'Welch p = 0.006) and contain more genes (31 against 22, p = 0.019), which is what a '
         'sparse occupancy process predicts and a within-class accumulation gradient does not. '
         'Consistent with that, the displacement score does not increase with size or gene count '
         'within either class (all p above 0.39).'),
        ('body',
         'Under the diagonal sparse model an event frequency of 1.95 per cent predicts an '
         'asymmetric prevalence of 36.40 per cent against an observed 36.36 per cent. Under the '
         'axis-aligned model the fitted frequency is higher and implies that 72.7 per cent of loci '
         'carry an event, which at first sight breaks this agreement. It does not: the axis model '
         'partitions that total into 35.3 per cent displaced along the deletion axis and 37.4 per '
         'cent along the duplication axis, each close to the observed prevalence. The published '
         'class label tracks the duplication axis, with 11 of 12 asymmetric-class loci displaced '
         'predominantly along it against 10 of 21 mirror-class loci (Fisher exact p = 0.022). The '
         'three largest deletion-axis movers are all labelled mirror. Under the axis model the '
         'prevalence match is evidence for one of two populations rather than for a single class.'),
        ('h1', 'Conclusions'),
        ('bullet',
         'Mirror symmetry is a major component of both traits but cannot be the whole story. It is '
         'rejected by single-locus fit in height, by burden direction in both traits, and by the '
         'duplication mean in posterior predictive checking.'),
        ('bullet',
         'Height asymmetry is sparse, discrete and stochastic. A deterministic drift is worse by '
         'more than 30 AIC units, so the occupancy structure carries real information.'),
        ('bullet',
         'The lower-left constraint in the published analysis was the binding one. Every model '
         'restricted to that quadrant sits on the 85 degree boundary. Released, the estimate moves '
         'past 90 degrees to buffering: the deletion effect is attenuated toward the null while the '
         'duplication effect carries the trait. This is a correction to the published '
         'parameterisation, not merely an additional model.'),
        ('bullet',
         'No single height model satisfies every line of evidence. Buffering wins conditional AIC, '
         'BIC and leave-one-out, both posterior predictive checks and the burden prediction; the '
         'axis-aligned mixture is the only model whose gene-count trajectory matches the observed '
         'trend. At 33 loci this should be reported as two competing descriptions rather than a '
         'settled winner, with the buffering result stated as the stronger of the two on '
         'likelihood grounds.'),
        ('bullet',
         'Adding saturating occupancy, the natural reconciliation, does not resolve it. It very '
         'nearly closes the gene-count gap but loses on every information criterion, overshoots the '
         'observed duplication burden, and places its saturation point beyond the largest observed '
         'CNV. The tension between the reciprocal-plane likelihood and the gene-count trajectory '
         'should be reported as an open question, not patched.'),
        ('bullet',
         'BMI asymmetry is a small aggregate drift proportional to gene count, adequately described '
         'by two parameters. Its frequency is not identifiable and should not be reported as an '
         'estimate. With the unidentified parameter removed, AIC and BIC agree for the first time.'),
        ('bullet',
         'The mechanistic implication of axis alignment is that asymmetric displacement reflects '
         'genes with one-sided dosage sensitivity, haploinsufficient or triplosensitive, rather '
         'than a single gene class pushing the trait in the same direction under both deletion and '
         'duplication. This is a different prediction about which genes drive asymmetry, and it is '
         'testable against dosage-sensitivity annotation.'),
        ('h1', 'Limitations'),
        ('body',
         'Mirror contributor frequency and per-gene magnitude are partially confounded, and the '
         'asymmetry frequency and displacement magnitude are weakly identified, severely so for '
         'BMI. The burden translation uses an effective density of dosage-responsive genes rather '
         'than a directly observed genome-wide count, and the BMI conversion assumes one BMI '
         'standard deviation equals 4.7 kg/m2. Macé et al. report duplication as well as deletion '
         'burden coefficients, so the duplication comparison is quantitative; their standard errors '
         'are not published and have been reconstructed from the coefficient and P value assuming a '
         'two-sided Wald test, so the intervals drawn against the observed values are themselves '
         'reconstructions. Auwerx is used only as a directional cross-check.'),
        ('body',
         'Model selection at 33 height loci is not decisive. The buffered model leads the '
         'axis-aligned model by about five AIC units and the published diagonal model by rather '
         'more, but it fails the gene-count trajectory, so the two descriptions are reported '
         'together rather than resolved. The buffered deletion component is large, and although the '
         'mirror component is calibrated on mirror-class loci only, so the asymmetric-class '
         'residuals driving the result are out of sample, a sensitivity analysis against mirror '
         'mis-specification has not been carried out. The frequency-effect ridge shown in the '
         'calibration figures is reconstructed in-package because the original upstream profile was '
         'not distributed with the analysis. These limitations affect exact parameter values more '
         'than the qualitative ordering of models.'),
    ])

    d = pdf.infodict()
    d['Title'] = 'Geometry of asymmetric reciprocal CNV effects on height and BMI'
    d['Subject'] = 'Diagonal versus axis-aligned displacement models and genome-wide burden'

print(f'wrote {OUT}')
