#!/usr/bin/env python
"""Main summary figure and companion table.

Layout, per trait band:
  response row   6 model columns, reciprocal DEL x DUP plane. Observed loci
                 (size proportional to gene count), the model's 95% predictive
                 region from simulated replicates, and the model's mean
                 prediction drawn as explicit mixture components whose line
                 weight tracks component probability as a function of gene count.
  burden row     model-predicted burden relative to the observed Mace estimate.
                 Deletion is the calibration anchor and is fixed at 1 by
                 construction; duplication is a free out-of-sample prediction.
  criteria strip conditional AIC / BIC / LOO and the two posterior predictive
                 checks, aligned under the model columns.

The winning model is tinted per band, because the winner differs between traits.

Run standalone to write report/Figure_main_summary.pdf and Table1 CSV, or import
build_main_figure() / table1() from the report builder.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; PALE = '#f2f2f2'
RED = '#b03a2e'; BLUE = '#1f5c8b'; GREEN = '#3f7d52'; PURPLE = '#6b4c93'; ORANGE = '#c47b1a'
WIN = '#fdf3d8'

TEAL = '#1b7f79'; PINK = '#a3417a'; BROWN = '#7a5230'

KEYS = ['mirror', 'continuous', 'single', 'axis_continuous', 'axis_sparse',
        'open_continuous', 'open_single', 'sat_single', 'aggregate']
LABELS = {
    'mirror': 'Mirror alone',
    'continuous': 'Mirror + continuous directional shift',
    'single': 'Mirror + sparse single hit',
    'axis_continuous': 'Mirror + continuous, axis-aligned',
    'axis_sparse': 'Mirror + sparse single hit, axis-aligned',
    'open_continuous': 'Mirror + continuous buffering',
    'open_single': 'Mirror + sparse single hit with buffering',
    'sat_single': 'Mirror + buffering, saturating occupancy',
    'aggregate': 'Mirror + aggregate directional shift (q dropped)',
}
HEAD = {'mirror': 'Mirror\nalone', 'continuous': 'Continuous\ndiagonal',
        'single': 'Sparse hit\ndiagonal', 'axis_continuous': 'Continuous\naxis-aligned',
        'axis_sparse': 'Sparse hit\naxis-aligned', 'open_continuous': 'Continuous\nbuffering',
        'open_single': 'Sparse hit\nbuffering', 'sat_single': 'Buffering,\nsaturating occ.',
        'aggregate': 'Aggregate\ndrift (2 par)'}
COLS = {'mirror': GREY, 'continuous': BLUE, 'single': PURPLE,
        'axis_continuous': GREEN, 'axis_sparse': RED, 'open_continuous': TEAL,
        'open_single': PINK, 'sat_single': BROWN, 'aggregate': ORANGE}


def load(P, trait, loci):
    d = pd.read_csv(P / trait / 'inputs' / loci)
    return dict(
        key=trait,
        dat=d,
        base=pd.read_csv(P / trait / 'outputs' / f'{trait}_mirror_parameters.csv').iloc[0].to_dict(),
        cmp=pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_model_comparison.csv'),
        bur=pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_burden_predictions.csv'),
        ppc=pd.read_csv(P / trait / 'outputs' / f'{trait}_axis_mixture_posterior_predictive_checks.csv'),
        clouds=np.load(P / trait / 'outputs' / f'{trait}_ppc_clouds.npz'))


def load_all(P):
    H = load(P, 'height', 'sacks_height_all_loci_exact.csv')
    B = load(P, 'bmi', 'sacks_bmi_complete_reciprocal_loci_exact.csv')
    H.update(sign=-1, name='Height', unit='cm/Mb')
    B.update(sign=+1, name='BMI', unit='kg/m$^2$/Mb')
    mace = pd.read_csv(P / 'code' / 'mace_burden_reference.csv')
    return H, B, mace


def components(T, key):
    """Mixture components of a model's mean prediction as a function of gene count.

    Returns a list of (DEL array, DUP array, weight array, style) tuples. Weight is
    the component probability at each gene count, used to modulate line weight.
    """
    d, b, s = T['dat'], T['base'], T['sign']
    G = np.linspace(d.Gene_count.min(), d.Gene_count.max(), 120)
    d0 = b['p'] * b['sign_mean'] * b['m']
    md = G * d0
    mu = -b['rho'] * G * d0
    row = T['cmp'].set_index('model').loc[LABELS[key]]
    one = np.ones_like(G)
    if key == 'mirror':
        return [(md, mu, one, 'solid')]
    bD, bU = float(row.shift_DEL), float(row.shift_DUP)
    if key in ('continuous', 'open_continuous'):
        q = float(row.q)
        return [(md + s * G * q * bD, mu + s * G * q * bU, one, 'solid')]
    if key == 'axis_continuous':
        q, w = float(row.q), float(row.w_deletion_share)
        return [(md + s * G * q * w * bD, mu + s * G * q * (1 - w) * bU, one, 'solid')]
    if key == 'aggregate':
        return [(md + s * G * bD, mu + s * G * bU, one, 'solid')]
    q = float(row.q)
    if key == 'sat_single':
        gs = float(row.gene_saturation)
        occ = 1 - (1 - q) ** (gs * G / (gs + G))
        return [(md, mu, 1 - occ, 'solid'),
                (md + s * bD, mu + s * bU, occ, 'solid')]
    occ = 1 - (1 - q) ** G
    if key in ('single', 'open_single'):
        return [(md, mu, 1 - occ, 'solid'),
                (md + s * bD, mu + s * bU, occ, 'solid')]
    w = float(row.w_deletion_share)
    return [(md, mu, 1 - occ, 'solid'),
            (md + s * bD, mu, occ * w, 'dashed'),
            (md, mu + s * bU, occ * (1 - w), 'dotted')]


def expected_given_G(T, key, G):
    """Exact model-implied E[DEL | G] and E[DUP | G], averaged over the mixture.

    For the sparse models the shift enters through the occupancy 1-(1-q)^G, which is
    non-linear in gene count, so these curves must be evaluated rather than
    linearised. That distinction matters: a single slope understates the
    gene-count dependence the sparse models actually imply.
    """
    b, s = T['base'], T['sign']
    G = np.asarray(G, float)
    d0 = b['p'] * b['sign_mean'] * b['m']
    md = G * d0
    mu = -b['rho'] * G * d0
    row = T['cmp'].set_index('model').loc[LABELS[key]]
    if key == 'mirror':
        return md, mu
    bD, bU = float(row.shift_DEL), float(row.shift_DUP)
    if key in ('continuous', 'open_continuous'):
        q = float(row.q)
        return md + s * G * q * bD, mu + s * G * q * bU
    if key == 'axis_continuous':
        q, w = float(row.q), float(row.w_deletion_share)
        return md + s * G * q * w * bD, mu + s * G * q * (1 - w) * bU
    if key == 'aggregate':
        return md + s * G * bD, mu + s * G * bU
    q = float(row.q)
    if key == 'sat_single':
        gs = float(row.gene_saturation)
        occ = 1 - (1 - q) ** (gs * G / (gs + G))
        return md + s * occ * bD, mu + s * occ * bU
    occ = 1 - (1 - q) ** G
    if key in ('single', 'open_single'):
        return md + s * occ * bD, mu + s * occ * bU
    w = float(row.w_deletion_share)
    return md + s * occ * w * bD, mu + s * occ * (1 - w) * bU


def ivw_slope(G, y, se=None):
    """Inverse-variance weighted regression slope of y on G (unweighted if se is None)."""
    G = np.asarray(G, float); y = np.asarray(y, float)
    w = np.ones_like(G) if se is None else 1.0 / np.asarray(se, float) ** 2
    Gm = np.average(G, weights=w); ym = np.average(y, weights=w)
    sxx = np.sum(w * (G - Gm) ** 2)
    return float(np.sum(w * (G - Gm) * (y - ym)) / sxx), float(np.sqrt(1.0 / sxx))


def gene_count_figure(H, B):
    """Observed effect versus gene count, with exact model-implied mean curves.

    The reciprocal-plane panels and the posterior predictive statistics are both
    largely blind to this: across-locus correlation is dominated by the residual
    covariance, so every BMI model reproduces it. The gene-count trajectory is
    where the models genuinely disagree.
    """
    fig, axs = plt.subplots(3, 2, figsize=(7.4, 7.8),
                            gridspec_kw=dict(height_ratios=[1, 1, .80], hspace=.48, wspace=.30))
    fig.subplots_adjust(left=.175, right=.985, top=.930, bottom=.130)
    for c, T in enumerate([H, B]):
        d = T['dat']
        Gobs = d.Gene_count.to_numpy(float)
        grid = np.linspace(Gobs.min(), Gobs.max(), 200)
        for r, (col, secol, nm) in enumerate([('Deletion_beta', 'Deletion_SE', 'Deletion'),
                                              ('Duplication_beta', 'Duplication_SE', 'Duplication')]):
            ax = axs[r, c]
            ax.errorbar(Gobs, d[col], yerr=d[secol], fmt='o', ms=3.2, color=INK,
                        ecolor=LIGHT, elinewidth=.6, capsize=0, lw=0, mfc='none', mew=.7,
                        zorder=3, label='observed loci')
            for k in KEYS:
                ed, eu = expected_given_G(T, k, grid)
                ax.plot(grid, ed if r == 0 else eu, color=COLS[k], lw=1.35,
                        label=HEAD[k].replace('\n', ' '), zorder=4,
                        ls='--' if k == 'aggregate' else '-')
            ax.axhline(0, color=LIGHT, lw=.6)
            ax.set_xlabel('Gene count', size=7.2)
            ax.set_ylabel(f'{nm} effect (SD)', size=7.2)
            ax.set_title(f"{T['name']}: {nm.lower()} effect vs gene count", size=8, weight='bold')
            ax.tick_params(labelsize=6.4)

        # implied vs observed slope, duplication only
        ax = axs[2, c]
        obs_b, obs_se = ivw_slope(Gobs, d.Duplication_beta, d.Duplication_SE)
        lo, hi = obs_b - 1.96 * obs_se, obs_b + 1.96 * obs_se
        ax.axvspan(lo, hi, color=INK, alpha=.13, zorder=1)
        ax.axvline(obs_b, color=INK, lw=1.2, zorder=3)
        y = np.arange(len(KEYS))[::-1]
        slopes = []
        for yy, k in zip(y, KEYS):
            _, eu = expected_given_G(T, k, Gobs)
            mb, _ = ivw_slope(Gobs, eu, d.Duplication_SE)
            slopes.append(mb)
            inside = lo <= mb <= hi
            ax.plot([mb], [yy], marker='o', ms=6,
                    mfc=COLS[k] if inside else 'white', mec=COLS[k], mew=1.3, zorder=4)
            ax.text(mb, yy + .30, f'{mb:+.5f}', ha='center', size=5.2, color=INK, zorder=5)
        span = [min(slopes + [lo]), max(slopes + [hi])]
        pad = (span[1] - span[0]) * .22
        ax.set_xlim(span[0] - pad, span[1] + pad)
        ax.set_yticks(y)
        ax.set_yticklabels([HEAD[k].replace('\n', ' ') for k in KEYS] if c == 0 else [], size=6)
        ax.set_xlabel('implied slope of DUP effect on gene count', size=6.8)
        ax.set_title(f"{T['name']}: implied vs observed gene-count trend", size=8, weight='bold')
        ax.tick_params(labelsize=6.0)
        ax.locator_params(axis='x', nbins=5)
        ax.set_ylim(-1.75, len(KEYS) - .1)   # headroom below the lowest marker for the note
        ax.text(.02, .015, f'observed {obs_b:+.5f}   95% CI [{lo:+.5f}, {hi:+.5f}]\n'
                           'filled = inside interval, open = outside',
                transform=ax.transAxes, ha='left', va='bottom', size=5.5, color=INK)
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=6.4, ncol=4, loc='lower center',
               bbox_to_anchor=(.55, .012))
    return fig


def region_contour(ax, ds, us, color, levels=(0.50, 0.95)):
    """Contours containing the given fractions of the simulated predictive mass.

    Smoothed before contouring; the raw 2-D histogram of a finite replicate sample
    produces speckle islands that read as structure but are Monte Carlo noise.
    """
    from scipy.ndimage import gaussian_filter
    x = ds.ravel(); y = us.ravel()
    hb, xe, ye = np.histogram2d(x, y, bins=70)
    hb = gaussian_filter(hb, 2.0)
    hs = np.sort(hb.ravel())[::-1]
    csum = np.cumsum(hs) / hs.sum()
    cuts = sorted(hs[min(np.searchsorted(csum, lv), len(hs) - 1)] for lv in levels)
    xc = .5 * (xe[1:] + xe[:-1]); yc = .5 * (ye[1:] + ye[:-1])
    ax.contourf(xc, yc, hb.T, levels=[cuts[0], hb.max() + 1], colors=[color],
                alpha=.13, zorder=1)
    ax.contour(xc, yc, hb.T, levels=cuts, colors=[color], linewidths=[.5, .8],
               alpha=.85, zorder=2)


def response_panel(ax, T, key, winner, xlim, ylim):
    d = T['dat']
    if key == winner:
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=WIN,
                               edgecolor='none', zorder=0))
    cl = T['clouds']
    region_contour(ax, cl[f'{key}_DEL'], cl[f'{key}_DUP'], COLS[key])
    for dd, uu, wt, ls in components(T, key):
        pts = np.column_stack([dd, uu]).reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lw = .5 + 2.4 * np.clip((wt[:-1] + wt[1:]) / 2, 0, 1)
        ax.add_collection(LineCollection(seg, linewidths=lw, colors=INK,
                                         linestyles=ls, zorder=4, alpha=.9))
    sz = 3 + 26 * (d.Gene_count / d.Gene_count.max())
    if 'Dose_response_class' in d:
        for cls, mk in [('mirror', 'o'), ('asymmetric', 's')]:
            z = d[d.Dose_response_class.str.lower() == cls]
            ax.scatter(z.Deletion_beta, z.Duplication_beta,
                       s=sz[z.index], marker=mk, facecolor='none',
                       edgecolor=INK, lw=.65, zorder=5)
    else:
        ax.scatter(d.Deletion_beta, d.Duplication_beta, s=sz, marker='o',
                   facecolor='none', edgecolor=INK, lw=.65, zorder=5)
    ax.axhline(0, color=LIGHT, lw=.5, zorder=1); ax.axvline(0, color=LIGHT, lw=.5, zorder=1)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.tick_params(labelsize=5.6, length=2, pad=1.5)


def burden_panel(ax, T, key, mace, winner):
    """Model burden as a multiple of the observed Mace estimate."""
    if key == winner:
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=WIN,
                               edgecolor='none', zorder=0))
    m = mace[mace.trait == T['key']].set_index('burden_model')
    obs_del = float(m.loc['deletion', 'beta'])
    obs_dup = float(m.loc['duplication', 'beta'])
    dup_lo = float(m.loc['duplication', 'ci_low_derived']) / obs_dup
    dup_hi = float(m.loc['duplication', 'ci_high_derived']) / obs_dup
    b = T['bur'].set_index('model').loc[LABELS[key]]
    pred_dup = float(b['predicted_DUP_effect_per_Mb_at_matched_density']) / obs_dup
    lo, hi = min(dup_lo, dup_hi), max(dup_lo, dup_hi)
    ax.add_patch(Rectangle((lo, .58), hi - lo, .84, facecolor=INK, alpha=.22,
                           edgecolor='none', zorder=1))
    ax.plot([1, 1], [.55, 1.45], color=INK, lw=1.1, zorder=3)
    ax.plot([1, 1], [-.45, .45], color=INK, lw=1.1, zorder=3)
    ax.plot([1], [0], marker='|', color=INK, ms=0)
    ax.barh([1], [pred_dup], height=.5, color=COLS[key], edgecolor='white', lw=.5, zorder=2)
    ax.barh([0], [1.0], height=.5, color=COLS[key], edgecolor='white', lw=.5, alpha=.35, zorder=2)
    ax.axvline(0, color=GREY, lw=.6, zorder=1)
    ax.set_xlim(-2.2, 2.0); ax.set_ylim(-.6, 1.6)
    ax.set_yticks([0, 1])
    ax.tick_params(labelsize=5.6, length=2, pad=1.5)
    ok = pred_dup > 0
    ax.text(.97, 1.30, f'{pred_dup:+.2f}', transform=ax.get_yaxis_transform(),
            va='center', ha='right', size=5.6, color=INK if ok else RED,
            weight='normal' if ok else 'bold', zorder=6)


def criteria_strip(ax, T, winner):
    cm = T['cmp'].set_index('model')
    pp = T['ppc'].set_index(['model', 'statistic'])
    rows = [('$\\Delta$AIC', [float(cm.loc[LABELS[k], 'delta_conditional_AIC']) for k in KEYS], 'delta'),
            ('$\\Delta$BIC', [float(cm.loc[LABELS[k], 'delta_conditional_BIC']) for k in KEYS], 'delta'),
            ('$\\Delta$LOO', [float(cm.loc[LABELS[k], 'delta_conditional_LOO']) for k in KEYS], 'delta'),
            ('PPC DUP mean', [float(pp.loc[(LABELS[k], 'DUP_mean'), 'posterior_predictive_CDF']) for k in KEYS], 'ppc'),
            ('PPC DEL-DUP corr.', [float(pp.loc[(LABELS[k], 'DEL_DUP_correlation'), 'posterior_predictive_CDF']) for k in KEYS], 'ppc')]
    n = len(rows)
    for r, (name, vals, kind) in enumerate(rows):
        yy = n - 1 - r
        for c, v in enumerate(vals):
            if kind == 'delta':
                best = v <= 1e-9
                shade = '#e7efe9' if v < 2 else ('#f6f6f6' if v < 6 else '#ffffff')
                face = '#cfe3d5' if best else shade
                txt = f'{v:.1f}'; col = INK
            else:
                fail = (v < .025) or (v > .975)
                face = '#f6dcd8' if fail else '#e7efe9'
                txt = f'{v:.3f}'; col = RED if fail else GREEN
            ax.add_patch(Rectangle((c, yy), 1, 1, facecolor=face, edgecolor='white', lw=.8))
            ax.text(c + .5, yy + .5, txt, ha='center', va='center', size=4.9, color=col,
                    weight='bold' if (kind == 'delta' and v <= 1e-9) else 'normal')
        ax.text(-.10, yy + .5, name, ha='right', va='center', size=5.4, color=INK)
    ax.set_xlim(0, len(KEYS)); ax.set_ylim(0, n); ax.axis('off')


def build_main_figure(H, B, mace):
    fig = plt.figure(figsize=(8.1, 10.4))
    outer = GridSpec(3, 1, figure=fig, height_ratios=[.74, 2.40, 2.40],
                     left=.120, right=.988, top=.965, bottom=.082, hspace=.46)

    # --- Panel A: the motivating observation -------------------------------
    top = GridSpecFromSubplotSpec(1, 3, outer[0], width_ratios=[1, 1, 1.15], wspace=.42)
    d = H['dat']; b = H['base']
    G = d.Gene_count.to_numpy(float); d0 = b['p'] * b['sign_mean'] * b['m']
    zD = -(d.Deletion_beta.to_numpy(float) - G * d0) / d.Deletion_SE.to_numpy(float)
    zU = -(d.Duplication_beta.to_numpy(float) + b['rho'] * G * d0) / d.Duplication_SE.to_numpy(float)
    ax = fig.add_subplot(top[0])
    for cls, c, mk in [('mirror', BLUE, 'o'), ('asymmetric', RED, 's')]:
        s = (d.Dose_response_class.str.lower() == cls).to_numpy()
        ax.scatter(zD[s], zU[s], s=14, color=c, marker=mk, edgecolor='white', lw=.4, label=cls)
    lim = np.abs(np.r_[zD, zU]).max() * 1.12
    ax.plot([0, lim], [0, lim], color=GREY, lw=.7, ls='--')
    ax.axhline(0, color=LIGHT, lw=.6); ax.axvline(0, color=LIGHT, lw=.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel('DEL residual (z)', size=6.4); ax.set_ylabel('DUP residual (z)', size=6.4)
    ax.tick_params(labelsize=5.6, length=2, pad=1.5)
    ax.set_title('Height residual geometry', size=7, weight='bold')
    ax.legend(frameon=False, fontsize=5.4, loc='lower right', handletextpad=.3)
    ax.text(.02, .97, 'A', transform=ax.transAxes, size=10, weight='bold', va='top')

    ax = fig.add_subplot(top[1])
    ang = np.degrees(np.arctan2(zU, zD))
    m = (d.Dose_response_class.str.lower() == 'mirror').to_numpy()
    ax.hist([ang[m], ang[~m]], bins=np.arange(-180, 181, 30), stacked=True,
            color=[BLUE, RED], edgecolor='white', lw=.4)
    for x, lab in [(0, 'DEL'), (90, 'DUP')]:
        ax.axvline(x, color=INK, lw=.7, ls=':')
        ax.text(x, .93, lab, size=5.2, ha='center', va='top', color=INK,
                transform=ax.get_xaxis_transform(),
                bbox=dict(boxstyle='round,pad=.15', fc='white', ec='none', alpha=.8))
    ax.axvline(45, color=GREY, lw=.7, ls='--')
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel('displacement angle ($\\degree$)', size=6.4)
    ax.set_ylabel('loci', size=6.4)
    ax.tick_params(labelsize=5.6, length=2, pad=1.5)
    ax.set_title('Individual loci move along one axis', size=7, weight='bold')

    ax = fig.add_subplot(top[2]); ax.axis('off')
    ax.text(0, 1.02, 'Published class vs displacement axis', size=6.5,
            weight='bold', va='top')
    tbl = ['             DUP-dom.  other',
           'asymmetric      11       1',
           'mirror          10      11']
    for i, line in enumerate(tbl):
        ax.text(0, .80 - .135 * i, line, size=5.6, va='top', family='monospace')
    ax.text(0, .36, 'Fisher exact $P$ = 0.022', size=6.0, va='top')
    ax.text(0, .19, 'The published asymmetric class tracks\nthe duplication axis; the three largest\n'
                    'deletion-axis movers are labelled mirror.', size=5.5, va='top', color='#444')

    xlims = {'height': (-1.15, .45), 'bmi': (-.55, .75)}
    ylims = {'height': (-.72, .52), 'bmi': (-.55, .55)}
    winners = {'height': 'open_single', 'bmi': 'aggregate'}
    panel_letters = {'height': ('B', 'C', 'D'), 'bmi': ('E', 'F', 'G')}

    for bi, T in enumerate([H, B]):
        band = GridSpecFromSubplotSpec(3, len(KEYS), outer[1 + bi],
                                       height_ratios=[1.55, .95, 1.05], hspace=.62, wspace=.16)
        w = winners[T['key']]
        L = panel_letters[T['key']]
        first = {}
        for c, k in enumerate(KEYS):
            ax = fig.add_subplot(band[0, c])
            response_panel(ax, T, k, w, xlims[T['key']], ylims[T['key']])
            ax.set_title(HEAD[k], size=5.8, weight='bold', color=COLS[k], pad=3)
            if c: ax.set_yticklabels([])
            else:
                ax.set_ylabel('DUP effect (SD)', size=6.2); first['resp'] = ax
            ax.set_xlabel('DEL effect (SD)', size=6.2, labelpad=1)

            ax = fig.add_subplot(band[1, c])
            burden_panel(ax, T, k, mace, w)
            if c:
                ax.set_yticklabels([])
            else:
                ax.set_yticklabels(['DEL\n(calibrated)', 'DUP\n(predicted)'], size=5.2)
                first['bur'] = ax
            ax.set_xlabel('pred / obs', size=5.0, labelpad=1)
        ax = fig.add_subplot(band[2, :])
        criteria_strip(ax, T, w)
        first['crit'] = ax

        # Panel letters and the band label live in the left margin, so they cannot
        # collide with axis labels or panel titles.
        pr = first['resp'].get_position()
        pc = first['crit'].get_position()
        fig.text(.014, pr.y1 + .004, L[0], size=10, weight='bold', va='top')
        fig.text(.014, first['bur'].get_position().y1 + .004, L[1], size=10,
                 weight='bold', va='top')
        fig.text(.014, pc.y1 + .004, L[2], size=10, weight='bold', va='top')
        fig.text(.048, (pr.y1 + pc.y0) / 2, T['name'].upper(), size=11, weight='bold',
                 rotation=90, va='center', ha='center')

    fig.legend(handles=[
        Line2D([], [], marker='o', mfc='none', mec=INK, lw=0, ms=5, label='observed locus (size $\\propto$ gene count)'),
        Line2D([], [], color=INK, lw=2, label='model mean, line weight $\\propto$ component probability'),
        Line2D([], [], color=INK, lw=1, ls='--', label='DEL-shifted component'),
        Line2D([], [], color=INK, lw=1, ls=':', label='DUP-shifted component'),
        Rectangle((0, 0), 1, 1, facecolor=BLUE, alpha=.13, edgecolor=BLUE, lw=.8,
                  label='50% / 95% predictive region'),
        Rectangle((0, 0), 1, 1, facecolor=INK, alpha=.28, edgecolor='none',
                  label='95% CI of observed Mace estimate'),
    ], frameon=False, fontsize=5.9, ncol=3, loc='lower center', bbox_to_anchor=(.55, .008))
    return fig


def table1(H, B, mace):
    rows = []
    for T in [H, B]:
        cm = T['cmp'].set_index('model')
        bu = T['bur'].set_index('model')
        pp = T['ppc'].set_index(['model', 'statistic'])
        m = mace[mace.trait == T['key']].set_index('burden_model')
        obs_dup = float(m.loc['duplication', 'beta'])
        G = T['dat'].Gene_count.to_numpy(float)
        for k in KEYS:
            L = LABELS[k]
            q = cm.loc[L, 'q']
            occ = float(np.mean(1 - (1 - q) ** G)) if (k in ('single', 'axis_sparse')
                                                       and np.isfinite(q)) else np.nan
            rows.append({
                'Trait': T['name'], 'Model': HEAD[k].replace('\n', ' '),
                'k': int(cm.loc[L, 'n_parameters_added']),
                'q': q, 'w_DEL': cm.loc[L, 'w_deletion_share'],
                'shift_DEL': cm.loc[L, 'shift_DEL'], 'shift_DUP': cm.loc[L, 'shift_DUP'],
                'dAIC': cm.loc[L, 'delta_conditional_AIC'],
                'dBIC': cm.loc[L, 'delta_conditional_BIC'],
                'dLOO': cm.loc[L, 'delta_conditional_LOO'],
                'PPC_DUP_mean': pp.loc[(L, 'DUP_mean'), 'posterior_predictive_CDF'],
                'PPC_corr': pp.loc[(L, 'DEL_DUP_correlation'), 'posterior_predictive_CDF'],
                'genes_per_Mb': bu.loc[L, 'effective_dosage_genes_per_Mb_to_match_DEL'],
                'pred_DUP_burden': bu.loc[L, 'predicted_DUP_effect_per_Mb_at_matched_density'],
                'obs_DUP_burden': obs_dup,
                'pred_over_obs': bu.loc[L, 'predicted_DUP_effect_per_Mb_at_matched_density'] / obs_dup,
                'implied_affected_fraction': occ,
            })
    return pd.DataFrame(rows)


def table1_page(df, figsize=(11, 8.5)):
    """Render Table 1 as a landscape page."""
    fig = plt.figure(figsize=figsize)
    fig.text(.5, .955, 'Table 1. Model comparison, posterior predictive behaviour and burden '
                       'prediction', ha='center', size=12, weight='bold')
    fig.text(.5, .925, 'All fits conditional on the calibrated mirror component. '
                       'Burden columns are calibrated on the observed deletion coefficient, so the '
                       'duplication prediction is out of sample.',
             ha='center', size=7.6, style='italic', color='#444')
    show = df.copy()
    fmt = {'q': '{:.4f}', 'w_DEL': '{:.3f}', 'shift_DEL': '{:.4f}', 'shift_DUP': '{:.4f}',
           'dAIC': '{:.2f}', 'dBIC': '{:.2f}', 'dLOO': '{:.2f}', 'PPC_DUP_mean': '{:.3f}',
           'PPC_corr': '{:.3f}', 'genes_per_Mb': '{:.2f}', 'pred_DUP_burden': '{:+.4f}',
           'obs_DUP_burden': '{:+.4f}', 'pred_over_obs': '{:.2f}',
           'implied_affected_fraction': '{:.3f}'}
    for c, f in fmt.items():
        show[c] = show[c].map(lambda v: '' if pd.isna(v) else f.format(v))
    hdr = ['Trait', 'Model', 'k', 'q', 'w$_{DEL}$', 'shift$_{DEL}$', 'shift$_{DUP}$',
           '$\\Delta$AIC', '$\\Delta$BIC', '$\\Delta$LOO', 'PPC\nDUP mean', 'PPC\ncorr.',
           'genes\nper Mb', 'pred. DUP\nburden', 'obs. DUP\nburden', 'pred/\nobs',
           'implied\naffected']
    ax = fig.add_axes([.02, .12, .96, .74]); ax.axis('off')
    # Explicit widths: the model names need roughly twice the share of a numeric column.
    widths = [.050, .130, .026] + [.0567] * (len(hdr) - 3)
    tab = ax.table(cellText=show.values, colLabels=hdr, loc='center', cellLoc='center',
                   colWidths=widths)
    tab.auto_set_font_size(False); tab.set_fontsize(6.4); tab.scale(1, 1.5)
    ncol = len(hdr)
    for (r, c), cell in tab.get_celld().items():
        cell.set_linewidth(.4); cell.set_edgecolor('#cccccc')
        if r == 0:
            cell.set_facecolor('#ececec'); cell.set_text_props(weight='bold')
        else:
            rec = df.iloc[r - 1]
            best = (rec.dAIC == 0) or (rec.dLOO == 0)
            cell.set_facecolor(WIN if (rec.dAIC == 0) else
                               ('#f7f7f7' if rec.Trait == 'BMI' else 'white'))
    fig.text(.02, .065,
             'k, added parameters. q, per-gene event frequency (not identified for the continuous '
             'BMI models). w$_{DEL}$, fraction of events on the deletion axis. Shifts in trait SD. '
             'Deltas are relative to the best model within each trait. PPC values are posterior '
             'predictive CDFs; below 0.025 or above 0.975 is a failure. Burden in cm/Mb (height) '
             'and kg/m$^2$/Mb (BMI); observed values from Mace et al. 2017 Supplementary Table 1. '
             'Implied affected fraction applies to the sparse models only.',
             size=6.4, color='#444', wrap=True)
    return fig


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-dir', default=str(Path(__file__).resolve().parent.parent))
    a = ap.parse_args()
    P = Path(a.project_dir).resolve()
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                         'axes.linewidth': .6, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    H, B, mace = load_all(P)
    fig = build_main_figure(H, B, mace)
    out = P / 'report' / 'Figure_main_summary.pdf'
    fig.savefig(out); fig.savefig(out.with_suffix('.png'), dpi=400)
    print('wrote', out)
    df = table1(H, B, mace)
    df.to_csv(P / 'report' / 'Table1_model_summary.csv', index=False)
    f2 = table1_page(df)
    f2.savefig(P / 'report' / 'Table1_model_summary.pdf')
    print('wrote Table1_model_summary.csv / .pdf')


