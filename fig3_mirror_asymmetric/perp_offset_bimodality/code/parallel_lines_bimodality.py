#!/usr/bin/env python
"""Bimodality as one line vs two parallel mirror lines, with measurement error.

Each locus is summarised by its vertical distance from a mirror line of slope -rho,
    r_i = DUP_i + rho * DEL_i,
and r_i is modelled with two variance components:
    r_i ~ N(a, tau^2 + v_i),   v_i = seDUP_i^2 + rho^2 * seDEL_i^2
tau^2 is the true scatter of loci about their line (shared by all loci); v_i is the
locus's own sampling variance, known from its SEs, so noisy loci count for less.

  one line   (a, rho, tau)
  two lines  (a1, a2, pi, rho, tau): a mixture of two parallel lines, same slope and tau

Two versions:
  joint      rho is estimated together with the lines from ALL loci, with no class labels
  fixed rho  rho held at the mirror-class calibration used by perp_offset_bimodality.py,
             i.e. the current test with only the measurement error added

LR = 2 (logL two lines - logL one line). p is a parametric bootstrap: DUP is redrawn
around the fitted one-line model (each locus with its own SE and its observed DEL), both
models are refitted the same way, and the p is the share of draws with LR >= observed.

Observed data only (no model simulations). Run from the folder above code/:
    python code/parallel_lines_bimodality.py --nboot 2000
Writes outputs/parallel_lines/{parallel_lines_test.csv, parallel_lines_per_locus.csv}.
"""
from __future__ import annotations
import argparse, math
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import norm

HERE = Path(__file__).resolve().parent.parent
V2 = HERE.parent / 'CNV_height_BMI_corrected_directional_models_v2'
SEED = 20261013
P_CUTOFF = 0.05
LOG2PI = math.log(2 * math.pi)


def _norm_logpdf(x, m, var):
    return -0.5 * (LOG2PI + np.log(var) + (x - m) ** 2 / var)


def nll_one(th, D, U, sD2, sU2, rho_fixed):
    a, ltau = th[0], th[1]
    rho = rho_fixed if rho_fixed is not None else math.exp(th[2])
    var = math.exp(2 * ltau) + sU2 + rho * rho * sD2
    return -float(_norm_logpdf(U + rho * D, a, var).sum())


def nll_two(th, D, U, sD2, sU2, rho_fixed):
    a1, a2, lpi, ltau = th[:4]
    rho = rho_fixed if rho_fixed is not None else math.exp(th[4])
    pi_ = expit(lpi)
    r = U + rho * D
    var = math.exp(2 * ltau) + sU2 + rho * rho * sD2
    l = np.logaddexp(math.log(max(pi_, 1e-300)) + _norm_logpdf(r, a1, var),
                     math.log(max(1 - pi_, 1e-300)) + _norm_logpdf(r, a2, var))
    return -float(l.sum())


def _min(f, x0):
    r = minimize(f, x0, method='Nelder-Mead',
                 options={'maxiter': 6000, 'xatol': 1e-9, 'fatol': 1e-10})
    return r


def fit(D, U, seD, seU, rho_fixed=None, rho_start=0.35):
    """Best one-line and two-line fits. Two-line starts split the sorted residuals at the
    six widest gaps and at the 30/50/70% points (each side keeping >= 2 loci)."""
    sD2, sU2 = seD ** 2, seU ** 2
    free = rho_fixed is None
    rho0 = rho_start if free else rho_fixed
    r0 = U + rho0 * D
    x1 = [r0.mean(), math.log(max(r0.std(), 1e-3))] + ([math.log(rho0)] if free else [])
    f1 = lambda th: nll_one(th, D, U, sD2, sU2, rho_fixed)
    best1 = _min(f1, np.array(x1))
    for d in (-.5, .5):                          # a couple of rho restarts
        if free:
            alt = _min(f1, np.array([r0.mean(), x1[1], math.log(rho0) + d]))
            if alt.fun < best1.fun: best1 = alt
    rho1 = math.exp(best1.x[2]) if free else rho_fixed

    f2 = lambda th: nll_two(th, D, U, sD2, sU2, rho_fixed)
    best2 = None
    for rho_s in ({rho0, rho1} if free else {rho_fixed}):
        r = np.sort(U + rho_s * D); n = len(r)
        gaps = np.argsort(np.diff(r))[::-1] + 1
        splits = [k for k in gaps if 2 <= k <= n - 2][:6] + [int(q * n) for q in (.3, .5, .7)]
        for k in dict.fromkeys(splits):
            lo, hi = r[:k], r[k:]
            sd = math.sqrt(max(((lo - lo.mean()) ** 2).sum() + ((hi - hi.mean()) ** 2).sum(), 1e-8) / n)
            x = [lo.mean(), hi.mean(), math.log(k / (n - k)), math.log(max(sd, 1e-3))]
            if free: x.append(math.log(rho_s))
            res = _min(f2, np.array(x))
            if best2 is None or res.fun < best2.fun: best2 = res
    return best1, best2


def unpack(b1, b2, free, rho_fixed):
    a, tau = b1.x[0], math.exp(b1.x[1])
    rho1 = math.exp(b1.x[2]) if free else rho_fixed
    a1, a2, pi_, tau2 = b2.x[0], b2.x[1], expit(b2.x[2]), math.exp(b2.x[3])
    rho2 = math.exp(b2.x[4]) if free else rho_fixed
    if a1 > a2: a1, a2, pi_ = a2, a1, 1 - pi_   # line 1 = lower (displaced) line
    return dict(one_a=a, one_tau=tau, one_rho=rho1, two_a_low=a1, two_a_high=a2,
                two_weight_low=pi_, two_tau=tau2, two_rho=rho2,
                separation=a2 - a1, LR=max(0.0, 2 * (b1.fun - b2.fun)))   # two lines nest one line


def _boot_one(args):
    seed, D, seD, seU, a, rho, tau, rho_fixed = args
    rg = np.random.default_rng(seed)
    U = a - rho * D + rg.normal(0, np.sqrt(tau ** 2 + seU ** 2 + rho ** 2 * seD ** 2))
    b1, b2 = fit(D, U, seD, seU, rho_fixed=rho_fixed, rho_start=rho)
    return max(0.0, 2 * (b1.fun - b2.fun))


def run(dat, label, rho_fixed, nboot, pool):
    D = dat.Deletion_beta.to_numpy(float); U = dat.Duplication_beta.to_numpy(float)
    seD = dat.Deletion_SE.to_numpy(float); seU = dat.Duplication_SE.to_numpy(float)
    free = rho_fixed is None
    b1, b2 = fit(D, U, seD, seU, rho_fixed=rho_fixed)
    res = unpack(b1, b2, free, rho_fixed)
    jobs = [(SEED + i, D, seD, seU, res['one_a'], res['one_rho'], res['one_tau'], rho_fixed)
            for i in range(nboot)]
    null = np.array(pool.map(_boot_one, jobs, chunksize=25))
    res.update(version=label, rho='joint (free)' if free else f'fixed {rho_fixed:.3f}',
               n_loci=len(dat), p_bootstrap=(1 + int((null >= res['LR']).sum())) / (nboot + 1),
               null_median_LR=float(np.median(null)), null_p95_LR=float(np.quantile(null, .95)),
               n_boot=nboot)
    # posterior probability of the lower (displaced) line, per locus
    rho2 = res['two_rho']; r = U + rho2 * D
    var = res['two_tau'] ** 2 + seU ** 2 + rho2 ** 2 * seD ** 2
    l1 = np.log(res['two_weight_low']) + _norm_logpdf(r, res['two_a_low'], var)
    l2 = np.log(1 - res['two_weight_low']) + _norm_logpdf(r, res['two_a_high'], var)
    post = np.exp(l1 - np.logaddexp(l1, l2))
    cls = dat.Dose_response_class.str.lower().to_numpy()
    res['agree_with_classes'] = int(((post > .5) == (cls == 'asymmetric')).sum())
    per = pd.DataFrame(dict(version=label, rho=res['rho'], CNV_ID=dat.CNV_ID, cls=cls,
                            vertical_offset=r, offset_SE=np.sqrt(seU ** 2 + rho2 ** 2 * seD ** 2),
                            P_lower_line=post))
    return res, per


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nboot', type=int, default=2000)
    a = ap.parse_args()
    h = pd.read_csv(V2 / 'height' / 'inputs' / 'sacks_height_all_loci_exact.csv')
    pD = 2 * norm.sf(np.abs(h.Deletion_beta / h.Deletion_SE))
    pU = 2 * norm.sf(np.abs(h.Duplication_beta / h.Duplication_SE))
    kept = h[(pD < P_CUTOFF) | (pU < P_CUTOFF)].reset_index(drop=True)
    # the mirror-class calibrations used by the current test, for the fixed-rho version
    rho_cal = {}
    for lab, path in [('all 33 loci', HERE / 'outputs' / 'bimodality_test.csv'),
                      ('no-effect loci excluded', HERE / 'outputs' / 'excl_noeffect_p05' / 'bimodality_test.csv')]:
        t = pd.read_csv(path)
        rho_cal[lab] = float(t[(t.trait == 'Height') & (t.mirror_fit == 'mirror-class only')].rho_full_fit.iloc[0])
    rows, pers = [], []
    with Pool() as pool:
        for lab, dat in [('no-effect loci excluded', kept), ('all 33 loci', h)]:
            for rf in (None, rho_cal[lab]):
                res, per = run(dat, lab, rf, a.nboot, pool)
                rows.append(res); pers.append(per)
                print(f"{lab:<24s} rho {res['rho']:<13s} n={res['n_loci']}  LR={res['LR']:.2f}  "
                      f"p={res['p_bootstrap']:.4f}  rho_hat={res['two_rho']:.3f}  "
                      f"lines {res['two_a_low']:+.3f}/{res['two_a_high']:+.3f} "
                      f"(weight low {res['two_weight_low']:.2f})  tau={res['two_tau']:.3f} "
                      f"(one-line tau {res['one_tau']:.3f})  classes agree {res['agree_with_classes']}/{res['n_loci']}",
                      flush=True)
    out = HERE / 'outputs' / 'parallel_lines'
    out.mkdir(parents=True, exist_ok=True)
    cols = ['version', 'rho', 'n_loci', 'LR', 'p_bootstrap', 'null_median_LR', 'null_p95_LR',
            'two_rho', 'two_a_low', 'two_a_high', 'separation', 'two_weight_low', 'two_tau',
            'one_rho', 'one_a', 'one_tau', 'agree_with_classes', 'n_boot']
    pd.DataFrame(rows)[cols].to_csv(out / 'parallel_lines_test.csv', index=False)
    pd.concat(pers).to_csv(out / 'parallel_lines_per_locus.csv', index=False)
