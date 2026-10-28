#!/usr/bin/env python
"""Is the perpendicular offset from the mirror line bimodal?

The mirror model places every locus on the line DUP = -rho*DEL, with position along the
line set by gene count. Each locus's signed PERPENDICULAR distance from that line is the
part the mirror model does not explain. If asymmetry is a discrete class, those distances
should be a two-component mixture; if it is a continuous process, a single Gaussian.

The mirror is fitted to the FULL locus set, not the mirror-class subset. Fitting it on
mirror-class loci only pre-selects against exactly the loci that would form the second
component, which can manufacture the result.

Test: likelihood ratio of a 2-component against a 1-component Gaussian, with the null
distribution obtained by parametric bootstrap (simulate from the fitted single Gaussian,
refit both models, repeat). Mixture LRTs do not follow chi-square, so the bootstrap null
is the test, not an approximation to it.

Two nulls are reported:
  A  plain single Gaussian fitted to the observed distances (as specified)
  B  the same latent Gaussian plus each locus's own measurement error, so the null carries
     the heteroscedastic SEs the loci actually have

Both the observed data and every bootstrap replicate are fitted with the identical EM
procedure and the same number of random starts. That matters: fitting the observed data
harder than the nulls biases p downward.
"""
from __future__ import annotations
import argparse, math, pickle
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2

SEED = 20261010
HERE = Path(__file__).resolve().parent.parent


# ------------------------------------------------------------------ mirror fit
def bvlogpdf(x, y, m1, m2, v1, v2, cov):
    det = v1 * v2 - cov * cov
    if np.any(det <= 1e-14) or np.any(~np.isfinite(det)):
        return None
    dx, dy = x - m1, y - m2
    quad = (v2 * dx * dx + v1 * dy * dy - 2 * cov * dx * dy) / det
    return -np.log(2 * np.pi) - .5 * np.log(det) - .5 * quad


def fit_mirror(dat, pi_positive, nstarts=14):
    """Exact port of fit_mirror() from the published package."""
    D = dat.Deletion_beta.to_numpy(float)
    U = dat.Duplication_beta.to_numpy(float)
    G = dat.Gene_count.to_numpy(float)
    seD = dat.Deletion_SE.to_numpy(float)
    seU = dat.Duplication_SE.to_numpy(float)
    mus = 2 * pi_positive - 1

    def unpack(x):
        return dict(p=float(expit(x[0])), m=float(np.exp(x[1])), rho=float(np.exp(x[2])),
                    tauD=float(np.exp(x[3])), tauU=float(np.exp(x[4])), sign_mean=mus)

    def nll(x):
        z = unpack(x)
        p, m, rho, tD, tU = z['p'], z['m'], z['rho'], z['tauD'], z['tauU']
        mg = p * mus * m
        vg = p * 2 * m * m - (p * mus * m) ** 2
        if vg <= 1e-12 or not np.isfinite(vg):
            return 1e100
        v = G * vg
        ll = bvlogpdf(D, U, G * mg, -rho * G * mg, v + seD ** 2 + tD ** 2,
                      rho * rho * v + seU ** 2 + tU ** 2, -rho * v)
        return 1e100 if ll is None else float(-ll.sum())

    p0 = .4 if pi_positive < .5 else .03
    m0 = max(.02, abs(D.mean()) / (max(abs(mus), .1) * G.mean() * p0))
    rho0 = max(.15, min(1.0, abs(np.polyfit(D, U, 1)[0])))
    x0 = np.array([math.log(p0 / (1 - p0)), math.log(m0), math.log(rho0),
                   math.log(.02), math.log(.05)])
    best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 107 * i + int(pi_positive * 1000))
        r = minimize(nll, x0 + rg.normal(0, [1.4, 1, .6, 1, 1]), method='Nelder-Mead',
                     options={'maxiter': 8000, 'xatol': 1e-8, 'fatol': 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    return unpack(best.x)


# --------------------------------------------------- vectorised mixture fitting
def ll1(X):
    """Log-likelihood of the MLE single Gaussian, per row."""
    mu = X.mean(1, keepdims=True)
    sd = np.maximum(np.sqrt(((X - mu) ** 2).mean(1, keepdims=True)), 1e-12)
    ll = (-0.5 * np.log(2 * np.pi) - np.log(sd) - 0.5 * ((X - mu) / sd) ** 2).sum(1)
    return ll, mu[:, 0], sd[:, 0]


def _norm(X, m, s):
    return np.exp(-0.5 * ((X - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))


def ll2(X, n_starts=25, iters=200, seed=0, equal_var=False, tol=1e-8):
    """Best 2-component Gaussian mixture per row, by EM from many random starts.

    Identical treatment for observed and bootstrap rows. Each component sd is floored at
    5% of that row's own sd, so a component cannot collapse onto a single point and send
    the likelihood to infinity.

    equal_var=True constrains the two components to a common sd. That is the alternative
    implied by the sparse single-hit mechanism: a hit displaces a locus along the
    perpendicular axis but does not change its residual dispersion, so the components
    differ in mean only. It costs 2 extra parameters rather than 3, and is better behaved
    at n~33 because it cannot chase a tight cluster of outliers.
    """
    B, n = X.shape
    rg = np.random.default_rng(seed)
    floor = 0.05 * X.std(1, keepdims=True)
    sd0 = X.std(1, keepdims=True)
    best = np.full(B, -np.inf)
    bpar = np.zeros((B, 5))
    for _ in range(n_starts):
        i1 = rg.integers(0, n, B)
        i2 = rg.integers(0, n, B)
        m1 = np.take_along_axis(X, i1[:, None], 1).astype(float).copy()
        m2 = np.take_along_axis(X, i2[:, None], 1).astype(float).copy()
        same = (m1 == m2)[:, 0]
        if same.any():
            m1[same, 0] = X[same].min(1)
            m2[same, 0] = X[same].max(1)
        s1 = sd0.copy()
        s2 = sd0.copy()
        p = np.full((B, 1), 0.5)
        prev = None
        for it in range(iters):
            d1 = p * _norm(X, m1, s1)
            d2 = (1 - p) * _norm(X, m2, s2)
            w = d1 / np.maximum(d1 + d2, 1e-300)
            sw = np.clip(w.sum(1, keepdims=True), 1e-4 * n, (1 - 1e-4) * n)
            sw2 = n - sw
            p = sw / n
            m1 = (w * X).sum(1, keepdims=True) / sw
            m2 = ((1 - w) * X).sum(1, keepdims=True) / sw2
            if equal_var:
                sp = np.sqrt(((w * (X - m1) ** 2).sum(1, keepdims=True)
                              + ((1 - w) * (X - m2) ** 2).sum(1, keepdims=True)) / n)
                s1 = s2 = np.maximum(sp, floor)
            else:
                s1 = np.maximum(np.sqrt((w * (X - m1) ** 2).sum(1, keepdims=True) / sw), floor)
                s2 = np.maximum(np.sqrt(((1 - w) * (X - m2) ** 2).sum(1, keepdims=True) / sw2), floor)
            if it % 10 == 9:                      # periodic convergence check
                cur = np.log(np.maximum(p * _norm(X, m1, s1)
                                        + (1 - p) * _norm(X, m2, s2), 1e-300)).sum(1)
                if prev is not None and np.max(np.abs(cur - prev)) < tol:
                    break
                prev = cur
        ll = np.log(np.maximum(p * _norm(X, m1, s1) + (1 - p) * _norm(X, m2, s2), 1e-300)).sum(1)
        better = ll > best
        best = np.where(better, ll, best)
        cur = np.column_stack([p[:, 0], m1[:, 0], s1[:, 0], m2[:, 0], s2[:, 0]])
        bpar[better] = cur[better]
    return best, bpar


def lrt(X, n_starts, seed, equal_var=False):
    a, mu, sd = ll1(X)
    b, par = ll2(X, n_starts=n_starts, seed=seed, equal_var=equal_var)
    return 2 * (b - a), mu, sd, par


# ------------------------------------- measurement-error-aware mixture fitting
# Each offset x_i carries its own known sampling SE s_i, so the variance of locus i about
# its group mean is tau^2 + s_i^2: tau is the true scatter of loci about their line,
# shared by all loci. Noisy loci therefore count for less, and their scatter is not
# mistaken for a second group. Same vectorised layout as ll1/ll2: one row per dataset.
#
# tau^2 is updated by the ML fixed point for heteroscedastic random effects,
#   tau^2 = sum((R_i - s_i^2) / v_i^2) / sum(1 / v_i^2),  v_i = tau^2 + s_i^2,
# where R_i is the (expected) squared residual; the final log-likelihood is exact.

def _me_logpdf(X, m, v):
    return -0.5 * (np.log(2 * np.pi * v) + (X - m) ** 2 / v)


def _tau2_step(R, S2, t2):
    v = t2 + S2
    return np.maximum(((R - S2) / v ** 2).sum(1, keepdims=True)
                      / (1 / v ** 2).sum(1, keepdims=True), 0.0)


def me_ll1(X, S, iters=300):
    """Best single group with measurement error, per row: returns ll, mean, tau."""
    S2 = np.broadcast_to(np.asarray(S, float) ** 2, X.shape)
    t2 = np.maximum(X.var(1, keepdims=True) - S2.mean(1, keepdims=True), 1e-8)
    for _ in range(iters):
        w = 1 / (t2 + S2)
        a = (w * X).sum(1, keepdims=True) / w.sum(1, keepdims=True)
        t2 = _tau2_step((X - a) ** 2, S2, t2)
    w = 1 / (t2 + S2)
    a = (w * X).sum(1, keepdims=True) / w.sum(1, keepdims=True)
    return _me_logpdf(X, a, t2 + S2).sum(1), a[:, 0], np.sqrt(t2[:, 0])


def me_ll2(X, S, n_starts=25, iters=300, seed=0):
    """Best two-group mixture (common tau) with measurement error, per row, by
    generalised EM from many random starts. Returns ll and (pi, m1, m2, tau)."""
    B, n = X.shape
    S2 = np.broadcast_to(np.asarray(S, float) ** 2, X.shape)
    rg = np.random.default_rng(seed)
    t20 = np.maximum(X.var(1, keepdims=True) - S2.mean(1, keepdims=True), 1e-8) / 4
    best = np.full(B, -np.inf)
    bpar = np.zeros((B, 4))
    for _ in range(n_starts):
        m1 = np.take_along_axis(X, rg.integers(0, n, B)[:, None], 1).astype(float).copy()
        m2 = np.take_along_axis(X, rg.integers(0, n, B)[:, None], 1).astype(float).copy()
        same = (m1 == m2)[:, 0]
        if same.any():
            m1[same, 0] = X[same].min(1); m2[same, 0] = X[same].max(1)
        t2 = t20.copy(); p = np.full((B, 1), 0.5)
        for _ in range(iters):
            v = t2 + S2
            l1 = np.log(p) + _me_logpdf(X, m1, v)
            l2 = np.log(1 - p) + _me_logpdf(X, m2, v)
            w = np.exp(l1 - np.logaddexp(l1, l2))
            p = np.clip(w.mean(1, keepdims=True), 1e-4, 1 - 1e-4)
            u = 1 / v
            m1 = (w * u * X).sum(1, keepdims=True) / np.maximum((w * u).sum(1, keepdims=True), 1e-300)
            m2 = ((1 - w) * u * X).sum(1, keepdims=True) / np.maximum(((1 - w) * u).sum(1, keepdims=True), 1e-300)
            t2 = _tau2_step(w * (X - m1) ** 2 + (1 - w) * (X - m2) ** 2, S2, t2)
        v = t2 + S2
        ll = np.logaddexp(np.log(p) + _me_logpdf(X, m1, v),
                          np.log(1 - p) + _me_logpdf(X, m2, v)).sum(1)
        better = ll > best
        best = np.where(better, ll, best)
        cur = np.column_stack([p[:, 0], m1[:, 0], m2[:, 0], np.sqrt(t2[:, 0])])
        bpar[better] = cur[better]
    return best, bpar


def me_posterior_low(X, S, par):
    """Per-locus probability of the lower (more negative) group under me_ll2's fit."""
    S2 = np.broadcast_to(np.asarray(S, float) ** 2, X.shape)
    p, m1, m2, tau = (par[:, i:i + 1] for i in range(4))
    v = tau ** 2 + S2
    l1 = np.log(p) + _me_logpdf(X, m1, v)
    l2 = np.log(1 - p) + _me_logpdf(X, m2, v)
    post1 = np.exp(l1 - np.logaddexp(l1, l2))
    return np.where(m1 <= m2, post1, 1 - post1)


# --------------------------------------------------------------------- analysis
def analyse(dat, pi_pos, name, nboot, n_starts, out_rows, per_locus, mirror_subset=False):
    """mirror_subset=True reproduces the earlier approach of calibrating the mirror on the
    mirror-class loci only. That excludes the very loci that would form the second
    component, so the offsets of the held-out asymmetric loci are measured against a line
    they had no part in defining. Reported here only for comparison."""
    fit_dat = (dat[dat.Dose_response_class.str.lower() == 'mirror']
               if mirror_subset and 'Dose_response_class' in dat.columns else dat)
    base = fit_mirror(fit_dat, pi_pos)
    rho = base['rho']
    D = dat.Deletion_beta.to_numpy(float)
    U = dat.Duplication_beta.to_numpy(float)
    seD = dat.Deletion_SE.to_numpy(float)
    seU = dat.Duplication_SE.to_numpy(float)
    k = math.sqrt(1 + rho * rho)
    x = (rho * D + U) / k                      # signed perpendicular offset
    se = np.sqrt(rho ** 2 * seD ** 2 + seU ** 2) / k
    n = len(x)

    _, mu, sd = ll1(x[None, :])
    mu = float(mu[0]); sd = float(sd[0])
    lat = math.sqrt(max(sd ** 2 - float(np.mean(se ** 2)), 1e-8))

    rg = np.random.default_rng(SEED + 7)
    A = rg.normal(mu, sd, size=(nboot, n))
    Bm = rg.normal(mu, lat, size=(nboot, n)) + rg.normal(0, se, size=(nboot, n))

    cls = (dat['Dose_response_class'] if 'Dose_response_class' in dat.columns
           else pd.Series(['unclassified'] * n)).to_numpy()
    store = {}
    for tag, eq, df in [('free_var', False, 3), ('equal_var', True, 2)]:
        obs, _, _, par = lrt(x[None, :], n_starts, seed=SEED, equal_var=eq)
        obs = float(obs[0]); par = par[0]
        la, _, _, _ = lrt(A, n_starts, seed=SEED + 1, equal_var=eq)
        lb, _, _, _ = lrt(Bm, n_starts, seed=SEED + 2, equal_var=eq)
        pA = (1 + int((la >= obs).sum())) / (nboot + 1)
        pB = (1 + int((lb >= obs).sum())) / (nboot + 1)
        pi_, m1, s1, m2, s2 = par
        if m1 > m2:                            # order components by mean
            pi_, m1, s1, m2, s2 = 1 - pi_, m2, s2, m1, s1
        out_rows.append(dict(
            trait=name, mirror_fit=('mirror-class only' if mirror_subset else 'full data'),
            alternative=tag, extra_params=df, n_loci=n,
            rho_full_fit=rho, mirror_p=base['p'], mirror_m=base['m'],
            offset_mean=mu, offset_sd=sd, mean_locus_SE=float(np.mean(se)), latent_sd=lat,
            LRT=obs, chi2_p_nominal=float(chi2.sf(obs, df)),
            minor_weight=float(min(pi_, 1 - pi_)),
            comp1_mean=m1, comp1_sd=s1, comp2_mean=m2, comp2_sd=s2,
            separation=abs(m2 - m1),
            p_bootstrap_nullA=pA, p_bootstrap_nullB=pB,
            null_median_LRT=float(np.median(la)), null_p95_LRT=float(np.quantile(la, .95)),
            null_frac_near_zero=float((la < 1e-6).mean()),
            n_boot=nboot, n_starts=n_starts))
        store[tag] = dict(obs=obs, par=(pi_, m1, s1, m2, s2),
                          nullA=la, nullB=lb, pA=pA, pB=pB)

    per_locus.append(pd.DataFrame(dict(
        trait=name, CNV_ID=dat.CNV_ID, Dose_response_class=cls,
        Gene_count=dat.Gene_count, Deletion_beta=D, Duplication_beta=U,
        perp_offset=x, perp_offset_SE=se)))
    return dict(name=name, x=x, se=se, rho=rho, mu=mu, sd=sd, cls=cls, fits=store,
                mirror_subset=mirror_subset)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nboot', type=int, default=5000)
    ap.add_argument('--starts', type=int, default=25)
    a = ap.parse_args()
    root = HERE.parent
    H = pd.read_csv(root / 'CNV_height_BMI_corrected_directional_models_v2' /
                    'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
    M = pd.read_csv(root / 'CNV_height_BMI_corrected_directional_models_v2' /
                    'bmi' / 'inputs' / 'sacks_bmi_complete_reciprocal_loci_exact.csv')
    rows, per = [], []
    res = [analyse(H, 0.30, 'Height', a.nboot, a.starts, rows, per),
           analyse(M, 11 / 14, 'BMI', a.nboot, a.starts, rows, per),
           # the mirror-class-only calibration: the mirror line is defined by the mirror
           # loci, so the asymmetric loci are measured against a line they did not shape
           analyse(H, 0.30, 'Height', a.nboot, a.starts, rows, per, mirror_subset=True),
           analyse(M, 11 / 14, 'BMI', a.nboot, a.starts, rows, per, mirror_subset=True)]
    tab = pd.DataFrame(rows)
    (HERE / 'outputs').mkdir(exist_ok=True)
    tab.to_csv(HERE / 'outputs' / 'bimodality_test.csv', index=False)
    pd.concat(per).to_csv(HERE / 'outputs' / 'perp_offsets_per_locus.csv', index=False)
    def tag_of(r):
        return f'{r["name"]}_{"mirrorclass" if r["mirror_subset"] else "fulldata"}'
    np.savez_compressed(
        HERE / 'outputs' / 'null_LRT_draws.npz',
        **{f'{tag_of(r)}_{t}_{k}': r['fits'][t][k]
           for r in res for t in ('free_var', 'equal_var') for k in ('nullA', 'nullB')})
    (HERE / 'outputs' / '_res.pkl').write_bytes(pickle.dumps(res))
    pd.set_option('display.width', 250)
    print(tab.T.to_string())
