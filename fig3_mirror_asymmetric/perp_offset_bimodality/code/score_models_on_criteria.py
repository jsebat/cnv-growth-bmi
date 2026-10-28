#!/usr/bin/env python
"""Score three models against the five empirical criteria for height.

The three models are all conditional on the same fitted mirror layer, and the two
directional models use the SAME angular freedom, so the only thing that differs
between them is occupancy -- how many genes in a CNV carry the directional effect:

    Mirror alone      no directional component at all
    + continuous      every gene contributes, displacement grows with gene count
    + single hit      a rare per-gene frequency q; one hit displaces the whole locus
                      by a fixed vector, and most loci are not hit at all

Criteria 1, 2, 3 and 5 are scored by simulation: for every model we generate replicate
datasets with the real gene counts and the real standard errors, then compute on each
replicate exactly the quantity we compute on the observed data. Nothing uses the
simulation's hidden state -- classes are re-inferred from the simulated offsets by the
same 2-component mixture used on the real data, so a continuous model is never
penalised for "not having classes" by fiat.

  1  size            asymmetric minus mirror difference in log gene count
  2  DUP effect      class separation on the duplication axis (and on the deletion
                     axis, which the observed data say should be near zero)
  3  bimodality      1-vs-2 component likelihood ratio for the perpendicular offsets
  5  dose response   inverse-variance weighted slope of duplication effect on gene
                     count WITHIN the majority (mirror-like) class

Criterion 4 is different in kind and is handled separately, by algebraic propagation
rather than simulation: the fitted per-gene effects are scaled so that the predicted
genome-wide DELETION burden matches Mace et al. 2017, and the predicted DUPLICATION
burden is then read off out of sample. A U-shaped burden requires both to be negative.
For the mirror model that is impossible by construction, not merely unlikely: its
per-gene duplication effect is -rho times its deletion effect, so the two burdens
always carry opposite signs, whatever the data.

--project-dir and --out-dir point the script at another copy of the v2 project (e.g. the
filtered one written by build_excl_noeffect_project.py); the defaults are the originals.

--me uses the measurement-error-aware mixture (me_ll1/me_ll2: each locus's offset SE
enters both fits) for the bimodality LR and for inferring the two classes, matching
me_bimodality.py. Without it the original fits are used.

--match-filter P applies the no-effect exclusion to the observed data AND to every
simulated replicate before scoring: a locus is dropped when neither its DEL nor its DUP
effect reaches p < P (two-sided, effect / its own SE). The models are still fitted on,
and simulate, all loci, so selection acts identically on both sides.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from perp_offset_bimodality import ll1, ll2, me_ll1, me_ll2, me_posterior_low  # noqa: E402

ROOT = HERE.parent.parent
IN = ROOT / 'CNV_height_BMI_corrected_directional_models_v2'
SEED = 20261011

# Both directional models use the free angle, so occupancy is the only difference.
MODELS = ['mirror', 'open_continuous', 'open_single']
NICE = {'mirror': 'Mirror alone',
        'open_continuous': 'Mirror + continuous',
        'open_single': 'Mirror + single hit'}
LABEL_IN_CSV = {
    'mirror': 'Mirror alone',
    'open_continuous': 'Mirror + continuous shift, free angle',
    'open_single': 'Mirror + single hit, free angle (buffering)',
}


def mirror_components(G, seD, seU, z):
    p, m, rho, tD, tU, mus = z['p'], z['m'], z['rho'], z['tauD'], z['tauU'], z['sign_mean']
    mean_gene = p * mus * m
    var_gene = max(p * 2 * m * m - (p * mus * m) ** 2, 1e-12)
    v = G * var_gene
    return (G * mean_gene, -rho * G * mean_gene,
            v + seD ** 2 + tD ** 2, rho * rho * v + seU ** 2 + tU ** 2, -rho * v)


def simulate(key, pars, G, seD, seU, base, sign, nsim, rng):
    """Replicate datasets under one fitted model. Returns (nsim, n) DEL and DUP."""
    n = len(G)
    md, mu, vd, vu, cv = mirror_components(G, seD, seU, base)
    mD = np.tile(md, (nsim, 1)); mU = np.tile(mu, (nsim, 1))
    VD = np.tile(vd, (nsim, 1)); VU = np.tile(vu, (nsim, 1)); CV = np.tile(cv, (nsim, 1))
    bD, bU = pars['shift_DEL'], pars['shift_DUP']
    if key in ('continuous', 'open_continuous'):
        q = pars['q']; vv = G * (2 * q - q * q)
        mD += sign * G * q * bD; mU += sign * G * q * bU
        VD += vv * bD * bD; VU += vv * bU * bU; CV += vv * bD * bU
    elif key in ('single', 'open_single', 'dup_only'):
        hit = rng.random((nsim, n)) < (1 - (1 - pars['q']) ** G)
        mD += sign * bD * hit; mU += sign * bU * hit
    a = np.sqrt(VD); b = CV / a; c = np.sqrt(np.maximum(VU - b * b, 1e-300))
    z1 = rng.standard_normal((nsim, n)); z2 = rng.standard_normal((nsim, n))
    return mD + a * z1, mU + b * z1 + c * z2


def score(D, U, G, seU, rho, n_starts=20, seD=None, me=False):
    """The four simulation-scored criteria, from observable quantities only.

    D, U are (nrep, n); G and seU are (n,) or, when replicates keep different loci,
    (nrep, n). Classes are inferred from the offsets by a 2-component mixture
    with a common variance, the same alternative used in the bimodality test. With me=True
    (needs seD) both fits carry each locus's offset SE and loci are assigned to the
    lower class when their posterior probability of it exceeds 0.5.
    """
    k = np.sqrt(1 + rho * rho)
    X = (rho * D + U) / k
    if me:
        S = np.sqrt(rho * rho * np.asarray(seD, float) ** 2 + np.asarray(seU, float) ** 2) / k
        S = np.broadcast_to(S, X.shape)
        a, _, _ = me_ll1(X, S)
        b, par = me_ll2(X, S, n_starts=n_starts, seed=SEED)
        lrt = 2 * (b - a)
        low = me_posterior_low(X, S, par) > .5
        pib = np.where(par[:, 1] <= par[:, 2], par[:, 0], 1 - par[:, 0])
    else:
        a, _, _ = ll1(X)
        b, par = ll2(X, n_starts=n_starts, seed=SEED, equal_var=True)
        lrt = 2 * (b - a)
        pi_, m1, s1, m2, s2 = (par[:, i] for i in range(5))
        # orient so component 1 is the displaced (more negative) one
        swap = m1 > m2
        m1b = np.where(swap, m2, m1); m2b = np.where(swap, m1, m2)
        pib = np.where(swap, 1 - pi_, pi_)
        # hard assignment at the midpoint between component means
        cut = 0.5 * (m1b + m2b)
        low = X < cut[:, None]
    nlow = low.sum(1)
    ok = (nlow >= 2) & (nlow <= X.shape[1] - 2)

    def sep(M):
        """displaced-class mean minus majority-class mean, per replicate"""
        M = np.broadcast_to(M, low.shape)
        lo = np.where(low, M, np.nan)
        hi = np.where(~low, M, np.nan)
        with np.errstate(invalid='ignore'):
            return np.nanmean(lo, 1) - np.nanmean(hi, 1)

    G = np.broadcast_to(G, X.shape)
    sepG = sep(np.log1p(G))                   # criterion 1, size
    sepD, sepU = sep(D), sep(U)               # criterion 2, effect-size separation
    # criterion 5: dose response within the majority (mirror-like) class
    slope = np.full(len(X), np.nan)
    w = 1.0 / np.broadcast_to(seU, X.shape) ** 2
    for i in range(len(X)):
        sel = ~low[i]
        if sel.sum() >= 5:
            gg, yy, ww = G[i][sel], U[i][sel], w[i][sel]
            Gm = np.average(gg, weights=ww)
            ym = np.average(yy, weights=ww)
            sxx = np.sum(ww * (gg - Gm) ** 2)
            if sxx > 0:
                slope[i] = np.sum(ww * (gg - Gm) * (yy - ym)) / sxx

    def msk(v):
        return np.where(ok, v, np.nan)

    return dict(lrt=lrt, sepG=msk(sepG), sepD=msk(sepD), sepU=msk(sepU), slope=slope,
                minor_w=np.minimum(pib, 1 - pib))


def has_effect(D, U, seD, seU, p_cut):
    """loci where DEL or DUP reaches p < p_cut, the complement of the no-effect rule"""
    return ((2 * norm.sf(np.abs(D / seD)) < p_cut) | (2 * norm.sf(np.abs(U / seU)) < p_cut))


def score_filtered(Ds, Us, G, seD, seU, rho, p_cut, n_starts, me=False):
    """score() on each replicate after dropping its own no-effect loci.

    Which loci survive differs between replicates, so replicates are batched by how many
    loci they keep and each batch is scored on its own kept loci.
    Returns the same fields as score(), plus n_kept per replicate.
    """
    keep = has_effect(Ds, Us, seD[None, :], seU[None, :], p_cut)
    n_kept = keep.sum(1)
    fields = ('lrt', 'sepG', 'sepD', 'sepU', 'slope', 'minor_w')
    out = {f: np.full(len(Ds), np.nan) for f in fields}
    for n in np.unique(n_kept):
        if n < 6:                            # too few loci to infer two classes
            continue
        rows = np.flatnonzero(n_kept == n)
        idx = np.vstack([np.flatnonzero(keep[i]) for i in rows])
        pick = lambda M: np.take_along_axis(M[rows], idx, 1)
        Gb = np.broadcast_to(G, Ds.shape); Sb = np.broadcast_to(seU, Ds.shape)
        SDb = np.broadcast_to(seD, Ds.shape)
        sc = score(pick(Ds), pick(Us), pick(Gb), pick(Sb), rho, n_starts=n_starts,
                   seD=pick(SDb), me=me)
        for f in fields:
            out[f][rows] = sc[f]
    out['n_kept'] = n_kept.astype(float)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nsim', type=int, default=2000)
    ap.add_argument('--starts', type=int, default=20)
    ap.add_argument('--project-dir', default=str(IN))
    ap.add_argument('--out-dir', default=str(HERE.parent / 'outputs'))
    ap.add_argument('--match-filter', type=float, default=None, metavar='P')
    ap.add_argument('--me', action='store_true')
    a = ap.parse_args()
    IN = Path(a.project_dir).resolve()

    h = pd.read_csv(IN / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
    cmp_ = pd.read_csv(IN / 'height' / 'outputs' /
                       'height_axis_mixture_model_comparison.csv').set_index('model')
    bur = pd.read_csv(IN / 'height' / 'outputs' /
                      'height_axis_mixture_burden_predictions.csv').set_index('model')
    base = pd.read_csv(IN / 'height' / 'outputs' /
                       'height_mirror_parameters.csv').iloc[0].to_dict()
    rho = base['rho']
    G = h.Gene_count.to_numpy(float)
    seD = h.Deletion_SE.to_numpy(float)
    seU = h.Duplication_SE.to_numpy(float)
    D = h.Deletion_beta.to_numpy(float)
    U = h.Duplication_beta.to_numpy(float)

    if a.match_filter is None:
        obs = score(D[None, :], U[None, :], G, seU, rho, n_starts=40, seD=seD, me=a.me)
        obs = {k: float(np.asarray(v)[0]) for k, v in obs.items()}
    else:
        k = has_effect(D, U, seD, seU, a.match_filter)
        obs = score(D[k][None, :], U[k][None, :], G[k], seU[k], rho, n_starts=40,
                    seD=seD[k], me=a.me)
        obs = {f: float(np.asarray(v)[0]) for f, v in obs.items()}
        obs['n_kept'] = int(k.sum())
        print(f'match filter p < {a.match_filter}: observed keeps {k.sum()} of {len(k)} loci')
    print('OBSERVED  sizeDiff %+.3f | sep DEL %+.3f | sep DUP %+.3f | LRT %.2f | '
          'mirror-class DUP slope %+.5f | minor wt %.2f'
          % (obs['sepG'], obs['sepD'], obs['sepU'], obs['lrt'], obs['slope'], obs['minor_w']))

    rng = np.random.default_rng(SEED)
    rows, store = [], {}
    for key in MODELS:
        pars = cmp_.loc[LABEL_IN_CSV[key]].to_dict()
        Ds, Us = simulate(key, pars, G, seD, seU, base, -1, a.nsim, rng)
        sc = (score(Ds, Us, G, seU, rho, n_starts=a.starts, seD=seD, me=a.me)
              if a.match_filter is None else
              score_filtered(Ds, Us, G, seD, seU, rho, a.match_filter, a.starts, me=a.me))
        store[key] = sc

        def pct(v, o):
            v = v[np.isfinite(v)]
            return float(np.mean(v <= o)) if len(v) else np.nan

        r = dict(model=NICE[key], key=key, n_sim=a.nsim)
        if 'n_kept' in sc:
            r.update(n_kept_med=float(np.median(sc['n_kept'])),
                     n_kept_lo=float(np.quantile(sc['n_kept'], .025)),
                     n_kept_hi=float(np.quantile(sc['n_kept'], .975)),
                     n_kept_obs_pct=float(np.mean(sc['n_kept'] <= obs['n_kept'])))
        for nm, fld in [('sizeG', 'sepG'), ('sepD', 'sepD'), ('sepU', 'sepU'),
                        ('LRT', 'lrt'), ('slope', 'slope')]:
            v = sc[fld]
            r[nm + '_med'] = float(np.nanmedian(v))
            r[nm + '_lo'] = float(np.nanquantile(v, .025))
            r[nm + '_hi'] = float(np.nanquantile(v, .975))
            r[nm + '_obs_pct'] = pct(v, obs[fld])
        # criterion 4 is algebraic, not simulated
        b = bur.loc[LABEL_IN_CSV[key]]
        r['burden_DEL_per_Mb'] = float(b.Mace_DEL_target_per_Mb)
        r['burden_DUP_per_Mb'] = float(b.predicted_DUP_effect_per_Mb_at_matched_density)
        r['burden_U_shaped'] = bool(r['burden_DEL_per_Mb'] < 0 and r['burden_DUP_per_Mb'] < 0)
        rows.append(r)
        print('%-22s size %+.3f | sepDEL %+.3f | sepDUP %+.3f [%+.3f,%+.3f] | LRT %5.1f '
              '| slope %+.5f | DUP burden %+.4f | U-shaped %s'
              % (NICE[key], r['sizeG_med'], r['sepD_med'], r['sepU_med'], r['sepU_lo'],
                 r['sepU_hi'], r['LRT_med'], r['slope_med'], r['burden_DUP_per_Mb'],
                 'YES' if r['burden_U_shaped'] else 'NO'))

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / 'model_criteria_scores.csv', index=False)
    pd.Series(obs).to_csv(out / 'model_criteria_observed.csv')
    np.savez_compressed(out / 'model_criteria_draws.npz',
                        **{k + '_' + f: store[k][f] for k in store
                           for f in ('lrt', 'sepG', 'sepD', 'sepU', 'slope')})
    print('\nwrote', out / 'model_criteria_scores.csv')
