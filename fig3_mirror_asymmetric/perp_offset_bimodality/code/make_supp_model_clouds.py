#!/usr/bin/env python
"""Supplementary figure: what the continuous and single-hit models produce, against the
observed height loci.

Top row    simulated clouds: each model's (DEL, DUP) effects with each locus's real SE, pooled
           over N_CLOUD simulated datasets of the 33 real gene counts. The single-hit
           cloud is coloured by whether the locus carries a hit.
Bottom row one typical dataset per model: the simulated replicate whose likelihood ratio
           (measurement-error-aware bimodality test) is closest to that model's median,
           drawn with each locus's real SEs -- i.e. what a single study of these 33 loci
           would see if the model were true.

Both models are the ones compared in Figure_height_main_v2 (mirror layer fitted on the
mirror-class loci and held fixed; directional component with a free angle). The typical
replicates are regenerated with score_models_on_criteria.simulate, its seed and its model
order, so they are the replicates whose LRs are stored in outputs/me/all_loci.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score_models_on_criteria as SC  # noqa: E402

ROOT = HERE.parent
IN = SC.IN
CRIT = ROOT / 'outputs' / 'me' / 'all_loci'
N_CLOUD = 150

def apply_style():
    """Only called when this module builds its own figure.

    make_height_figure imports cloud_data/frame/observed/cloud_key to draw the same
    clouds in its own serif style, so the rcParams change must not happen at import.
    """
    plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 8,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'axes.titlesize': 8.6, 'axes.labelsize': 8, 'axes.linewidth': .7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
    })


INK = '#1a1a1a'; LIGHT = '#dcdcdc'
# The simulated clouds are context and the observed loci are the foreground, so the
# observed points take ink rather than a hue. Orange on the hit cloud was the problem:
# the clouds are drawn at low alpha, so the hit cloud reads as #f4c3c3 over white and
# orange sits only dE = 19 from it (17.7 under protan/deutan simulation) -- above the
# categorical floor, but that floor is for series side by side, not for small marks
# embedded in a dense field of the other colour. Ink is dE = 64 from both clouds.
CLOUD_ALPHA = .30
NOHIT = '#2b6cb0'; HIT = '#c0392b'; OBS = INK; LINE = '#7b3294'

h = pd.read_csv(IN / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
cmp_ = pd.read_csv(IN / 'height' / 'outputs' /
                   'height_axis_mixture_model_comparison.csv').set_index('model')
base = pd.read_csv(IN / 'height' / 'outputs' / 'height_mirror_parameters.csv').iloc[0].to_dict()
rho = base['rho']
G = h.Gene_count.to_numpy(float)
seD = h.Deletion_SE.to_numpy(float); seU = h.Duplication_SE.to_numpy(float)
D_obs = h.Deletion_beta.to_numpy(float); U_obs = h.Duplication_beta.to_numpy(float)
cls = h.Dose_response_class.str.lower().to_numpy()
draws = np.load(CRIT / 'model_criteria_draws.npz')
n = len(G)


def sim_with_hits(key, pars, sD, sU, nsim, rng):
    """SC.simulate, plus the hit indicator it draws internally (single hit only).

    simulate's first draw for the single-hit model is the hit matrix, so replaying the
    generator state gives exactly the same hits as the simulate call that follows.
    """
    hit = None
    if key == 'open_single':
        st = rng.bit_generator.state
        hit = rng.random((nsim, n)) < (1 - (1 - pars['q']) ** G)
        rng.bit_generator.state = st
    D, U = SC.simulate(key, pars, G, sD, sU, base, -1, nsim, rng)
    return D, U, hit


def typical_replicates():
    """The replicate whose likelihood ratio is closest to each model's median."""
    rng = np.random.default_rng(SC.SEED)
    out = {}
    for key in SC.MODELS:
        pars = cmp_.loc[SC.LABEL_IN_CSV[key]].to_dict()
        lrt = draws[key + '_lrt']
        D, U, hit = sim_with_hits(key, pars, seD, seU, len(lrt), rng)
        i = int(np.nanargmin(np.abs(lrt - np.nanmedian(lrt))))
        out[key] = dict(D=D[i], U=U[i], hit=None if hit is None else hit[i],
                        lr=float(lrt[i]), rep=i, med=float(np.nanmedian(lrt)))
    return out


def cloud_data(n_cloud=N_CLOUD):
    """Pooled (DEL, DUP) draws for each model, every locus with its real SEs."""
    rng = np.random.default_rng(SC.SEED + 101)
    out = {}
    for key in ('open_continuous', 'open_single'):
        pars = cmp_.loc[SC.LABEL_IN_CSV[key]].to_dict()
        D, U, hit = sim_with_hits(key, pars, seD, seU, n_cloud, rng)
        out[key] = (D.ravel(), U.ravel(), None if hit is None else hit.ravel())
    return out


def draw_cloud(ax, key, D, U, hit):
    """One cloud panel: simulated draws, the observed loci, and the in-panel key."""
    frame(ax)
    if hit is None:
        ax.scatter(D, U, s=3, color=NOHIT, alpha=CLOUD_ALPHA, lw=0, zorder=3)
    else:
        ax.scatter(D[~hit], U[~hit], s=3, color=NOHIT, alpha=CLOUD_ALPHA, lw=0, zorder=3)
        ax.scatter(D[hit], U[hit], s=3, color=HIT, alpha=CLOUD_ALPHA, lw=0, zorder=3)
    observed(ax)
    cloud_key(ax, key)


SHOW = [('open_continuous', 'Mirror + continuous'), ('open_single', 'Mirror + single hit')]
XL, YL = (-1.5, 0.7), (-0.85, 0.62)


def frame(ax):
    ax.axhline(0, color=INK, lw=.5, zorder=1); ax.axvline(0, color=INK, lw=.5, zorder=1)
    xx = np.array(XL)
    ax.plot(xx, -rho * xx, color=LINE, lw=1.3, ls=(0, (5, 3)), zorder=2)
    ax.set_xlim(*XL); ax.set_ylim(*YL)
    ax.set_xlabel('Deletion effect (SD)'); ax.set_ylabel('Duplication effect (SD)')
    ax.set_xticks([-1.5, -1.0, -0.5, 0, 0.5]); ax.set_yticks([-0.8, -0.4, 0, 0.4])
    ax.grid(color='#eeeeee', lw=.5, zorder=0); ax.set_axisbelow(True)


def observed(ax, size=26):
    for c, mk in [('mirror', 'o'), ('asymmetric', '^')]:
        s = cls == c
        ax.scatter(D_obs[s], U_obs[s], s=size, marker=mk, color=OBS, edgecolor='white',
                   lw=.8, zorder=6)


def observed_key(ax):
    """Name the observed marker shapes inside the panel, upper right.

    For a figure where these panels are not accompanied by a shared legend.
    """
    for i, (lab, mk) in enumerate([('observed, mirror', 'o'),
                                   ('observed, asymmetric', '^')]):
        y = .955 - i * .082
        ax.plot([.585], [y], marker=mk, ms=4.5, color=OBS, mec='white', mew=.6, lw=0,
                transform=ax.transAxes, clip_on=False, zorder=7)
        ax.text(.625, y, lab, transform=ax.transAxes, size=6.4, color=INK,
                ha='left', va='center', zorder=7,
                bbox=dict(boxstyle='square,pad=.18', fc='white', ec='none', alpha=.78))


def cloud_key(ax, key, locus=False):
    """Name the clouds inside the panel they appear in.

    A shared legend cannot do this honestly: in the continuous panel every simulated
    locus is blue because that model has no hit/no-hit split, so a single legend entry
    reading "simulated, not hit" would mislabel panel A.
    """
    rows = ([('simulated loci', NOHIT)] if key == 'open_continuous'
            else [('simulated, not hit', NOHIT), ('simulated, hit', HIT)])
    # the key's mark matches the mark in that panel: filled dots for the clouds,
    # open diamonds for the single-replicate panels
    mk = dict(marker='D', mfc='none', mew=1.1, ms=4.5) if locus else         dict(marker='o', ms=4, alpha=.85)
    for i, (lab, col) in enumerate(rows):
        y = .035 + (len(rows) - 1 - i) * .075
        ax.plot([.585], [y + .016], color=col, mec=col, lw=0,
                transform=ax.transAxes, clip_on=False, zorder=7, **mk)
        ax.text(.625, y, lab, transform=ax.transAxes, size=6.4, color=INK,
                ha='left', va='bottom', zorder=7,
                bbox=dict(boxstyle='square,pad=.18', fc='white', ec='none', alpha=.78))


if __name__ == '__main__':
  apply_style()
  TYP = typical_replicates()
  CLOUD = cloud_data()
  LR_OBS = float(pd.read_csv(CRIT / 'model_criteria_observed.csv', index_col=0).iloc[:, 0]['lrt'])
  for AB_ONLY in (False, True):     # full 2x2 figure, then panels A-B alone
    if AB_ONLY:
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9), squeeze=False)
        fig.subplots_adjust(left=.09, right=.98, top=.76, bottom=.25, wspace=.28)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0))
        fig.subplots_adjust(left=.09, right=.98, top=.865, bottom=.13, wspace=.28, hspace=.55)
    row_dAIC = {k: (float(cmp_.loc[SC.LABEL_IN_CSV[k], 'delta_conditional_AIC']),
                    float(cmp_.loc[SC.LABEL_IN_CSV[k], 'delta_conditional_LOO'])) for k, _ in SHOW}

    for c, (key, name) in enumerate(SHOW):
        pars = cmp_.loc[SC.LABEL_IN_CSV[key]].to_dict()
        dA, dL = row_dAIC[key]
        # top: simulated cloud
        ax = axes[0, c]
        draw_cloud(ax, key, *CLOUD[key])
        sub = f'ΔAIC {dA:.1f}, ΔLOO {dL:.1f}'
        if key == 'open_single':
            sub += f';  P(hit | 23 genes) = {1 - (1 - pars["q"]) ** 23:.2f}'
        ax.set_title(name, weight='bold', loc='left', pad=24)
        ax.text(0, 1.015, f'{N_CLOUD} simulated datasets, with locus SE noise\n{sub}',
                transform=ax.transAxes, size=6.6, color=INK, va='bottom', linespacing=1.3)
        if AB_ONLY:
            continue
        # bottom: one typical dataset
        ax = axes[1, c]; frame(ax)
        t = TYP[key]
        if t['hit'] is None:
            ax.scatter(t['D'], t['U'], s=26, marker='D', facecolor='none', edgecolor=NOHIT,
                       lw=1.1, zorder=5)
        else:
            for sel, col in [(~t['hit'], NOHIT), (t['hit'], HIT)]:
                ax.scatter(t['D'][sel], t['U'][sel], s=26, marker='D', facecolor='none',
                           edgecolor=col, lw=1.1, zorder=5)
        observed(ax, size=20); cloud_key(ax, key, locus=True)
        ax.set_title(f'{name}: one typical dataset', weight='bold', loc='left', pad=24)
        ax.text(0, 1.015, f'33 loci, with locus SE noise\nlikelihood ratio {t["lr"]:.1f} '
                f'(model median {t["med"]:.1f}; observed {LR_OBS:.1f})',
                transform=ax.transAxes, size=6.6, color=INK, va='bottom', linespacing=1.3)

    for ax, L in zip(axes.ravel(), 'ABCD'):
        p = ax.get_position()
        fig.text(p.x0 - .06, p.y1 + (.10 if AB_ONLY else .055), L, size=12, weight='bold',
                 va='bottom')

    # The clouds are named in-panel by cloud_key, so the legend carries only what a
    # reader cannot get from the panel itself: the two observed marker shapes and the
    # mirror line.
    HANDLES = [
        Line2D([], [], marker='o', color=OBS, mec='white', mew=.8, lw=0, ms=6,
               label='observed, mirror (21)'),
        Line2D([], [], marker='^', color=OBS, mec='white', mew=.8, lw=0, ms=6,
               label='observed, asymmetric (12)'),
        Line2D([], [], color=LINE, lw=1.3, ls=(0, (5, 3)), label=f'mirror line, ρ = {rho:.2f}')]
    fig.legend(handles=HANDLES, frameon=False, ncol=3, loc='lower center',
               bbox_to_anchor=(.5, .015), fontsize=6.9)

    fig.suptitle('Continuous versus single-hit asymmetry: simulated and observed height effects',
                 size=10.5, weight='bold', y=.975 if not AB_ONLY else .965)
    out = ROOT / 'figures' / ('Figure_S_model_clouds' + ('_AB' if AB_ONLY else ''))
    fig.savefig(out.with_suffix('.pdf'))
    fig.savefig(out.with_suffix('.png'), dpi=200)
    print('wrote', out.with_suffix('.pdf'))
    plt.close(fig)
  pd.DataFrame([dict(model=k, replicate=v['rep'], LR=v['lr'], model_median_LR=v['med'],
                     n_hit=None if v['hit'] is None else int(v['hit'].sum()))
                for k, v in TYP.items()]).to_csv(CRIT / 'supp_model_clouds_typical.csv',
                                                 index=False)
  print({k: (v['rep'], round(v['lr'], 2), round(v['med'], 2)) for k, v in TYP.items()})
