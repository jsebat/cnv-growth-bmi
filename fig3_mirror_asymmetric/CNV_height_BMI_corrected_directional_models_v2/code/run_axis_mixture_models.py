#!/usr/bin/env python
"""Axis-aligned mixture models of reciprocal-CNV asymmetry.

The shipped models displace every affected locus by the SAME vector (b*cos(theta),
b*sin(theta)) with theta confined to (5, 85) degrees, so both components are always
non-zero: every locus moves diagonally.

This script adds the alternative geometry: each affected locus is displaced along
ONE axis only -- deletion-side with probability w, duplication-side with probability
1-w. The cloud still drifts diagonally in aggregate, but no individual locus does.
Mechanistically this is deletion-specific (haploinsufficiency) vs duplication-specific
(triplosensitivity) dosage sensitivity, rather than one gene class pushing the trait
the same way in both dosage states.

Two variants, matching the two occupancy structures already in use:

  axis_sparse      one rare hit displaces the locus, on one axis   (q, w, b_DEL, b_DUP)
  axis_continuous  many small events, a fraction w on each axis    (q, w, b_DEL, b_DUP)

The discriminating signature is the DEL-DUP covariance. A common diagonal shift adds
+vv*bD*bU to the covariance; axis-aligned events add nothing, because no single event
moves both coordinates.

All fits are conditional on the mirror component already calibrated in
<trait>/outputs/<trait>_mirror_parameters.csv, exactly as the shipped comparison is,
so conditional AIC/BIC/LOO are directly comparable to the published three models.

Usage:
  python code/run_axis_mixture_models.py --project-dir .
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

SEED = 20260726

ap = argparse.ArgumentParser()
ap.add_argument('--project-dir', default=str(Path(__file__).resolve().parent.parent))
ap.add_argument('--out-subdir', default='outputs')
ap.add_argument('--nstarts', type=int, default=14)
ap.add_argument('--skip-loo', action='store_true')
ap.add_argument('--max-q-continuous', type=float, default=0.8,
                help='Upper bound on the per-gene event frequency for the continuous models '
                     '(shipped value 0.8). Raise to test whether q is pinned at its bound.')
ap.add_argument('--max-q-sparse', type=float, default=0.2,
                help='Upper bound on q for the sparse single-hit models (shipped value 0.2).')
ap.add_argument('--tag', default='', help='Suffix for output filenames, e.g. "_maxq099".')
ap.add_argument('--angle-max', type=float, default=175.,
                help='Upper angular bound for the unconstrained-angle models (shipped models '
                     'are fixed at 85). Above 90 the deletion and duplication components of the '
                     'shift take opposite signs.')
args = ap.parse_args()
PROJECT = Path(args.project_dir).resolve()

# ---------------------------------------------------------------- mirror layer
def bvlogpdf(x, y, m1, m2, v1, v2, cov):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m1 = np.asarray(m1, float); m2 = np.asarray(m2, float)
    v1 = np.asarray(v1, float); v2 = np.asarray(v2, float); cov = np.asarray(cov, float)
    det = v1 * v2 - cov * cov
    if np.any(det <= 1e-14) or np.any(~np.isfinite(det)):
        return None
    dx = x - m1; dy = y - m2
    quad = (v2 * dx * dx + v1 * dy * dy - 2 * cov * dx * dy) / det
    return -np.log(2 * np.pi) - .5 * np.log(det) - .5 * quad

def mirror_components(dat, z):
    D = dat.Deletion_beta.to_numpy(float); U = dat.Duplication_beta.to_numpy(float)
    G = dat.Gene_count.to_numpy(float)
    seD = dat.Deletion_SE.to_numpy(float); seU = dat.Duplication_SE.to_numpy(float)
    p, m, rho, tD, tU, mus = z['p'], z['m'], z['rho'], z['tauD'], z['tauU'], z['sign_mean']
    mean_gene = p * mus * m
    var_gene = max(p * 2 * m * m - (p * mus * m) ** 2, 1e-12)
    v = G * var_gene
    md = G * mean_gene; mu = -rho * G * mean_gene
    vd = v + seD ** 2 + tD ** 2; vu = rho * rho * v + seU ** 2 + tU ** 2; cv = -rho * v
    return D, U, G, md, mu, vd, vu, cv

# ------------------------------------------------------------ likelihood cores
def nll_common(x, cmp_, trait_sign, model, max_q):
    """Shipped geometry: one shared shift vector at an estimated angle.

    The shipped models confine the angle to 5-85 degrees, forcing both components of
    the shift to have the same sign, so every affected locus moves diagonally into
    one quadrant. The 'open_' variants widen the bound to 5-175 degrees. That keeps
    the duplication component pointing toward the asymmetric direction but lets the
    deletion component take either sign.

    Above 90 degrees the fitted displacement attenuates the deletion effect toward
    the null while driving the duplication effect past the null in the same trait
    direction. That is BUFFERING in the sense of Sacks et al.: the deletion side is
    buffered rather than aggravated, and the resulting locus is asymmetric because
    the duplication side carries the effect, not because both dosage states are
    damaging. The shipped 5-85 window cannot represent it at any parameter value.
    """
    D, U, G, md, mu, vd, vu, cv = cmp_
    amin, amax = ANGLE[model]
    q = max_q * expit(x[0]); b = float(np.exp(np.clip(x[1], -50, 50)))
    theta = np.deg2rad(amin + (amax - amin) * expit(x[2]))
    bD = b * np.cos(theta); bU = b * np.sin(theta)
    if BASE[model] == 'continuous':
        vv = G * (2 * q - q * q)
        ll = bvlogpdf(D, U, md + trait_sign * G * q * bD, mu + trait_sign * G * q * bU,
                      vd + vv * bD * bD, vu + vv * bU * bU, cv + vv * bD * bU)
        return 1e100 if ll is None else float(-ll.sum())
    p0 = (1 - q) ** G
    l0 = bvlogpdf(D, U, md, mu, vd, vu, cv)
    l1 = bvlogpdf(D, U, md + trait_sign * bD, mu + trait_sign * bU, vd, vu, cv)
    if l0 is None or l1 is None: return 1e100
    ll0 = np.log(np.maximum(p0, 1e-300)) + l0
    ll1 = np.log(np.maximum(1 - p0, 1e-300)) + l1
    return float(-logsumexp(np.column_stack([ll0, ll1]), axis=1).sum())

def nll_aggregate(x, cmp_, trait_sign, model, max_q):
    """Reparameterised continuous model: only the identified quantity is estimated.

    The (q, b, theta) continuous model is not identified for BMI -- q runs to whatever
    upper bound it is given while b shrinks to hold the product q*b fixed. Only the
    aggregate per-gene shift is estimated by the data. This model drops q entirely and
    estimates that aggregate directly, as a deterministic per-gene drift with no added
    dispersion: 2 parameters instead of 3.

    Note this collapses the diagonal/axis-aligned distinction. For a continuous
    aggregate the geometry only enters through the added second moment (a shared
    diagonal shift contributes +vv*bD*bU to the covariance, axis-aligned events
    contribute nothing). With the variance term dropped as unidentified, the two
    geometries become the same two-parameter mean shift.
    """
    D, U, G, md, mu, vd, vu, cv = cmp_
    aD = float(np.exp(np.clip(x[0], -50, 50)))
    aU = float(np.exp(np.clip(x[1], -50, 50)))
    ll = bvlogpdf(D, U, md + trait_sign * G * aD, mu + trait_sign * G * aU, vd, vu, cv)
    if ll is None: return 1e100
    s = float(-ll.sum())
    return s if np.isfinite(s) else 1e100

def occupancy(G, q, gsat=None):
    """Probability that a CNV of G genes carries at least one dosage-sensitive gene.

    With gsat=None this is the shipped 1-(1-q)^G, which keeps rising across the whole
    observed gene-count range and never flattens. With a finite gsat the effective
    number of eligible genes saturates, Geff = gsat*G/(gsat+G), so occupancy
    approaches 1-(1-q)^gsat and its derivative falls to zero in large CNVs.

    Motivation: the buffered model fits the reciprocal plane best but implies that
    duplication effect DECREASES with gene count, because a large duplication
    displacement multiplied by steadily rising occupancy overwhelms the positive
    trend contributed by the mirror component. Saturating occupancy removes that
    coupling without touching the displacement geometry.
    """
    G = np.asarray(G, float)
    Geff = G if gsat is None else gsat * G / (gsat + G)
    return 1 - (1 - q) ** Geff


def unpack_sat(x, max_q, amin, amax):
    q = float(max_q * expit(x[0]))
    b = float(np.exp(np.clip(x[1], -50, 50)))
    theta = float(np.deg2rad(amin + (amax - amin) * expit(x[2])))
    gsat = float(np.exp(np.clip(x[3], -50, 50)))
    return q, b, theta, gsat


def nll_sat(x, cmp_, trait_sign, model, max_q):
    """Buffered sparse geometry with saturating occupancy (q, b, theta, gsat)."""
    D, U, G, md, mu, vd, vu, cv = cmp_
    amin, amax = ANGLE[model]
    q, b, theta, gsat = unpack_sat(x, max_q, amin, amax)
    if not np.isfinite(gsat) or gsat <= 0: return 1e100
    bD = b * np.cos(theta); bU = b * np.sin(theta)
    p1 = occupancy(G, q, gsat)
    p0 = 1 - p1
    l0 = bvlogpdf(D, U, md, mu, vd, vu, cv)
    l1 = bvlogpdf(D, U, md + trait_sign * bD, mu + trait_sign * bU, vd, vu, cv)
    if l0 is None or l1 is None: return 1e100
    ll0 = np.log(np.maximum(p0, 1e-300)) + l0
    ll1 = np.log(np.maximum(p1, 1e-300)) + l1
    s = float(-logsumexp(np.column_stack([ll0, ll1]), axis=1).sum())
    return s if np.isfinite(s) else 1e100


CS_SHARED = {'cont_single': False, 'cont_single_shared': True}


def unpack_cs(x, max_qc, max_qs, amin, amax, shared):
    """(q_cont, b_cont, theta_cont, q_hit, b_hit, theta_hit); theta_hit tied if shared."""
    ang = lambda t: float(np.deg2rad(amin + (amax - amin) * expit(t)))
    qc = float(max_qc * expit(x[0]))
    bc = float(np.exp(np.clip(x[1], -50, 50)))
    thc = ang(x[2])
    qs = float(max_qs * expit(x[3]))
    bs = float(np.exp(np.clip(x[4], -50, 50)))
    ths = thc if shared else ang(x[5])
    return qc, bc, thc, qs, bs, ths


def nll_cont_single(x, cmp_, trait_sign, model, bounds):
    """Three components: mirror + a weak continuous shift + a rare single hit.

    The two directional layers are superposed, not alternatives. Every locus carries
    the continuous displacement G*q_c*b_c and its added dispersion; on top of that a
    locus is hit with probability 1-(1-q_s)^G and, if hit, displaced by the fixed
    vector b_s. Setting q_c*b_c = 0 recovers the single-hit model and setting q_s = 0
    recovers the continuous model, so both are nested inside this one and the
    likelihood can only improve -- the question is whether it improves enough to pay
    for the extra parameters.

    Identifiability warning: q_c and b_c enter the mean only through the product
    q_c*b_c, and are separated only by the variance term G(2q_c-q_c^2)*b_c^2. With 33
    loci that separation is weak, which is the same pathology the aggregate model was
    written to document. The product is what the data constrain.
    """
    D, U, G, md, mu, vd, vu, cv = cmp_
    amin, amax = ANGLE[model]
    max_qc, max_qs = bounds
    qc, bc, thc, qs, bs, ths = unpack_cs(x, max_qc, max_qs, amin, amax, CS_SHARED[model])
    bcD = bc * np.cos(thc); bcU = bc * np.sin(thc)
    bsD = bs * np.cos(ths); bsU = bs * np.sin(ths)
    vv = G * (2 * qc - qc * qc)
    # continuous layer: shifts the mean of BOTH mixture components and inflates cov
    mdc = md + trait_sign * G * qc * bcD
    muc = mu + trait_sign * G * qc * bcU
    vdc = vd + vv * bcD * bcD
    vuc = vu + vv * bcU * bcU
    cvc = cv + vv * bcD * bcU
    p0 = (1 - qs) ** G
    l0 = bvlogpdf(D, U, mdc, muc, vdc, vuc, cvc)
    l1 = bvlogpdf(D, U, mdc + trait_sign * bsD, muc + trait_sign * bsU, vdc, vuc, cvc)
    if l0 is None or l1 is None: return 1e100
    ll0 = np.log(np.maximum(p0, 1e-300)) + l0
    ll1 = np.log(np.maximum(1 - p0, 1e-300)) + l1
    out = float(-logsumexp(np.column_stack([ll0, ll1]), axis=1).sum())
    return out if np.isfinite(out) else 1e100


def unpack_axis(x, max_q):
    q = max_q * expit(x[0])
    w = float(expit(x[1]))
    bD = float(np.exp(np.clip(x[2], -50, 50)))
    bU = float(np.exp(np.clip(x[3], -50, 50)))
    return float(q), w, bD, bU

def nll_axis(x, cmp_, trait_sign, model, max_q):
    """Axis-aligned geometry: each event moves DEL only or DUP only."""
    D, U, G, md, mu, vd, vu, cv = cmp_
    q, w, bD, bU = unpack_axis(x, max_q)
    if model == 'axis_continuous':
        qD = q * w; qU = q * (1 - w)
        vvD = G * (2 * qD - qD * qD); vvU = G * (2 * qU - qU * qU)
        # No added cross-covariance: no single event moves both coordinates.
        ll = bvlogpdf(D, U,
                      md + trait_sign * G * qD * bD,
                      mu + trait_sign * G * qU * bU,
                      vd + vvD * bD * bD, vu + vvU * bU * bU, cv)
        return 1e100 if ll is None else float(-ll.sum())
    # axis_sparse: three-component mixture per locus
    p0 = (1 - q) ** G
    lnone = bvlogpdf(D, U, md, mu, vd, vu, cv)
    ldel = bvlogpdf(D, U, md + trait_sign * bD, mu, vd, vu, cv)
    ldup = bvlogpdf(D, U, md, mu + trait_sign * bU, vd, vu, cv)
    if lnone is None or ldel is None or ldup is None: return 1e100
    a0 = np.log(np.maximum(p0, 1e-300)) + lnone
    aD = np.log(np.maximum((1 - p0) * w, 1e-300)) + ldel
    aU = np.log(np.maximum((1 - p0) * (1 - w), 1e-300)) + ldup
    s = float(-logsumexp(np.column_stack([a0, aD, aU]), axis=1).sum())
    return s if np.isfinite(s) else 1e100

# ------------------------------------------------------------------- fitting
# Naming. Every model below is conditional on the same fitted mirror layer, so the
# prefix is always "Mirror + ". Two things vary independently and the names say which:
#
#   OCCUPANCY  how many genes in a CNV carry the directional effect
#              "single hit"       rare per-gene frequency q, one hit displaces the
#                                 locus by a fixed vector; most loci are unaffected
#              "continuous shift" every gene contributes, so the displacement grows
#                                 with gene count
#   DIRECTION  where the displacement points, as an angle theta from the DEL axis
#              "constrained angle"  theta in (5, 85) deg  -- the shipped bound, which
#                                   forces the DEL and DUP components to share a sign
#              "free angle"         theta in (5, 175) deg -- lets them take opposite
#                                   signs, i.e. duplication-driven with the deletion
#                                   effect pulled back toward the mirror expectation
#                                   ("buffering")
#              "axis-aligned"       each hit moves along one axis only
#              "duplication axis only"  theta pinned at 90 deg
#
# The earlier labels called one geometry "sparse single hit" and the other
# "buffering", which read as two different mechanisms. They are not: both are
# single-hit mixtures with the same occupancy machinery and they differ only in the
# angular bound.
LABELS = {
    'mirror': 'Mirror alone',
    'continuous': 'Mirror + continuous shift, constrained angle',
    'single': 'Mirror + single hit, constrained angle',
    'axis_sparse': 'Mirror + single hit, axis-aligned',
    'axis_continuous': 'Mirror + continuous shift, axis-aligned',
    'aggregate': 'Mirror + aggregate drift (occupancy dropped)',
    'open_continuous': 'Mirror + continuous shift, free angle',
    'open_single': 'Mirror + single hit, free angle (buffering)',
    'sat_single': 'Mirror + single hit, free angle, saturating occupancy',
    'dup_only': 'Mirror + single hit, duplication axis only',
    'cont_single': 'Mirror + continuous + single hit (3 components)',
    'cont_single_shared': 'Mirror + continuous + single hit, shared angle',
}
MAXQ = {'continuous': args.max_q_continuous, 'single': args.max_q_sparse,
        'axis_continuous': args.max_q_continuous, 'axis_sparse': args.max_q_sparse,
        'open_continuous': args.max_q_continuous, 'open_single': args.max_q_sparse,
        'sat_single': args.max_q_sparse, 'dup_only': args.max_q_sparse,
        'cont_single': args.max_q_sparse, 'cont_single_shared': args.max_q_sparse}
# (continuous bound, single-hit bound) for the two-layer models
MAXQ_CS = (args.max_q_continuous, args.max_q_sparse)
# Which likelihood each key uses, and the angular window it may search.
BASE = {'continuous': 'continuous', 'single': 'single',
        'open_continuous': 'continuous', 'open_single': 'single',
        'sat_single': 'sat_single', 'dup_only': 'single'}
ANGLE = {'continuous': (5., 85.), 'single': (5., 85.),
         'open_continuous': (5., args.angle_max), 'open_single': (5., args.angle_max),
         'sat_single': (5., args.angle_max),
         # angle pinned at 90 degrees: the displacement is pure duplication,
         # so cos(theta) = 0 and the deletion component is exactly zero
         'dup_only': (90., 90.),
         'cont_single': (5., args.angle_max), 'cont_single_shared': (5., args.angle_max)}

def _cs_out(model, x, n):
    amin, amax = ANGLE[model]
    qc, bc, thc, qs, bs, ths = unpack_cs(x, *MAXQ_CS, amin, amax, CS_SHARED[model])
    return dict(q=qs, w_deletion_share=np.nan,
                shift_angle_deg=float(np.rad2deg(ths)), shift_magnitude=bs,
                shift_DEL=float(bs * np.cos(ths)), shift_DUP=float(bs * np.sin(ths)),
                cont_q=qc, cont_angle_deg=float(np.rad2deg(thc)), cont_magnitude=bc,
                cont_shift_DEL=float(bc * np.cos(thc)), cont_shift_DUP=float(bc * np.sin(thc)),
                cont_per_gene_DEL=float(qc * bc * np.cos(thc)),
                cont_per_gene_DUP=float(qc * bc * np.sin(thc)))

def fit_model(dat, base, trait_sign, model, nstarts):
    cmp_ = mirror_components(dat, base)
    n = len(dat)
    if model == 'mirror':
        ll = bvlogpdf(*[cmp_[0], cmp_[1]], *cmp_[3:])
        nllv = float(-ll.sum())
        return dict(model=LABELS[model], q=0., w=np.nan, shift_DEL=0., shift_DUP=0.,
                    negative_log_likelihood=nllv, n_parameters_added=0,
                    conditional_AIC=2 * nllv, conditional_BIC=2 * nllv, raw_par=None)
    if model == 'aggregate':
        x0 = np.array([math.log(2e-4), math.log(2e-3)])
        best = None
        for i in range(nstarts):
            rg = np.random.default_rng(SEED + 211 * i + 13)
            r = minimize(lambda x: nll_aggregate(x, cmp_, trait_sign, model, None),
                         x0 + rg.normal(0, [1.5, 1.5]), method='Nelder-Mead',
                         options={'maxiter': 12000, 'xatol': 1e-9, 'fatol': 1e-9})
            if best is None or r.fun < best.fun: best = r
        aD = float(np.exp(best.x[0])); aU = float(np.exp(best.x[1])); k = 2
        return dict(model=LABELS[model], q=np.nan, w_deletion_share=np.nan,
                    shift_DEL=aD, shift_DUP=aU,
                    negative_log_likelihood=float(best.fun), n_parameters_added=k,
                    conditional_AIC=float(2 * k + 2 * best.fun),
                    conditional_BIC=float(k * np.log(n) + 2 * best.fun), raw_par=best.x.tolist())
    max_q = MAXQ[model]
    if model == 'dup_only':
        # Observed class separation falls almost entirely on the duplication axis
        # (DEL p = 0.65, DUP p = 1e-06), so this fixes the deletion displacement at
        # zero and estimates only (q, b_DUP): two parameters instead of three.
        q0 = .018
        x0 = np.array([math.log(q0 / (max_q - q0)), math.log(.28)])
        best = None
        for i in range(nstarts):
            rg = np.random.default_rng(SEED + 211 * i + 41)
            r = minimize(lambda z: nll_common(np.array([z[0], z[1], 0.0]), cmp_,
                                              trait_sign, model, max_q),
                         x0 + rg.normal(0, [1.5, 1.2]), method='Nelder-Mead',
                         options={'maxiter': 12000, 'xatol': 1e-9, 'fatol': 1e-9})
            if best is None or r.fun < best.fun: best = r
        q = float(max_q * expit(best.x[0])); b = float(np.exp(best.x[1])); k = 2
        return dict(model=LABELS[model], q=q, w_deletion_share=np.nan,
                    shift_angle_deg=90.0, shift_magnitude=b,
                    shift_DEL=0.0, shift_DUP=b,
                    negative_log_likelihood=float(best.fun), n_parameters_added=k,
                    conditional_AIC=float(2 * k + 2 * best.fun),
                    conditional_BIC=float(k * np.log(n) + 2 * best.fun),
                    max_q=max_q, raw_par=[float(best.x[0]), float(best.x[1]), 0.0])
    if model == 'sat_single':
        amin, amax = ANGLE[model]
        q0 = .02; frac = (min(120., amax - 1.) - amin) / (amax - amin)
        x0 = np.array([math.log(q0 / (max_q - q0)), math.log(.30),
                       math.log(frac / (1 - frac)), math.log(25.)])
        best = None
        for i in range(nstarts):
            rg = np.random.default_rng(SEED + 211 * i + 31)
            r = minimize(lambda x: nll_sat(x, cmp_, trait_sign, model, max_q),
                         x0 + rg.normal(0, [1.5, 1.2, 1.6, 1.2]), method='Nelder-Mead',
                         options={'maxiter': 16000, 'xatol': 1e-9, 'fatol': 1e-9})
            if best is None or r.fun < best.fun: best = r
        q, b, theta, gsat = unpack_sat(best.x, max_q, amin, amax); k = 4
        return dict(model=LABELS[model], q=q, w_deletion_share=np.nan,
                    shift_angle_deg=float(np.rad2deg(theta)), shift_magnitude=b,
                    shift_DEL=float(b * np.cos(theta)), shift_DUP=float(b * np.sin(theta)),
                    gene_saturation=gsat,
                    negative_log_likelihood=float(best.fun), n_parameters_added=k,
                    conditional_AIC=float(2 * k + 2 * best.fun),
                    conditional_BIC=float(k * np.log(n) + 2 * best.fun),
                    max_q=max_q, raw_par=best.x.tolist())
    if model in CS_SHARED:
        shared = CS_SHARED[model]
        amin, amax = ANGLE[model]
        fr = lambda deg: (deg - amin) / (amax - amin)
        # start from the two separately fitted optima
        x0 = [math.log(.076 / (MAXQ_CS[0] - .076)), math.log(.085),
              math.log(fr(119.8) / (1 - fr(119.8))),
              math.log(.019 / (MAXQ_CS[1] - .019)), math.log(.40)]
        scales = [1.5, 1.2, 1.6, 1.5, 1.2]
        if not shared:
            x0.append(math.log(fr(119.3) / (1 - fr(119.3)))); scales.append(1.6)
        x0 = np.array(x0); k = len(x0)
        best = None
        for i in range(max(nstarts, 24)):
            rg = np.random.default_rng(SEED + 211 * i + 53 + (3 if shared else 0))
            start = x0 if i == 0 else x0 + rg.normal(0, scales)
            r = minimize(lambda x: nll_cont_single(x, cmp_, trait_sign, model, MAXQ_CS),
                         start, method='Nelder-Mead',
                         options={'maxiter': 24000, 'maxfev': 24000,
                                  'xatol': 1e-9, 'fatol': 1e-9})
            if best is None or r.fun < best.fun: best = r
        out = dict(model=LABELS[model], negative_log_likelihood=float(best.fun),
                   n_parameters_added=k, conditional_AIC=float(2 * k + 2 * best.fun),
                   conditional_BIC=float(k * np.log(n) + 2 * best.fun),
                   raw_par=best.x.tolist())
        out.update(_cs_out(model, best.x, n))
        return out
    axis = model.startswith('axis')
    if axis:
        q0 = .05 if model == 'axis_sparse' else .1
        x0 = np.array([math.log(q0 / (max_q - q0)), 0.0, math.log(.15), math.log(.15)])
        scales = [1.5, 1.5, 1.2, 1.2]
        fn = nll_axis; k = 4
    else:
        q0 = .08 if model == 'continuous' else .018
        b0 = .055 if model == 'continuous' else .25
        frac = (70 - 5) / 80.
        x0 = np.array([math.log(q0 / (max_q - q0)), math.log(b0), math.log(frac / (1 - frac))])
        scales = [1.5, 1.2, 1.6]
        fn = nll_common; k = 3
    best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 211 * i + (7 if axis else (1 if BASE[model] == 'continuous' else 2) + (17 if model.startswith('open') else 0)))
        r = minimize(lambda x: fn(x, cmp_, trait_sign, model, max_q),
                     x0 + rg.normal(0, scales), method='Nelder-Mead',
                     options={'maxiter': 12000, 'xatol': 1e-9, 'fatol': 1e-9})
        if best is None or r.fun < best.fun: best = r
    out = dict(model=LABELS[model], negative_log_likelihood=float(best.fun),
               n_parameters_added=k, conditional_AIC=float(2 * k + 2 * best.fun),
               conditional_BIC=float(k * np.log(n) + 2 * best.fun), raw_par=best.x.tolist())
    if axis:
        q, w, bD, bU = unpack_axis(best.x, max_q)
        out.update(q=q, w_deletion_share=w, shift_DEL=bD, shift_DUP=bU)
    else:
        q = max_q * expit(best.x[0]); b = float(np.exp(best.x[1]))
        amin, amax = ANGLE[model]
        theta = np.deg2rad(amin + (amax - amin) * expit(best.x[2]))
        out.update(q=float(q), w_deletion_share=np.nan, shift_angle_deg=float(np.rad2deg(theta)),
                   shift_magnitude=b, shift_DEL=float(b * np.cos(theta)), shift_DUP=float(b * np.sin(theta)))
    return out

def logpred_row(row, base, trait_sign, model, pars):
    model = BASE.get(model, model)   # open_* share the likelihood of their base model
    one = pd.DataFrame([row])
    D, U, G, md, mu, vd, vu, cv = mirror_components(one, base)
    f = lambda a: float(np.asarray(a)[0])
    D, U, G, md, mu, vd, vu, cv = map(f, (D, U, G, md, mu, vd, vu, cv))
    if model == 'mirror':
        return float(bvlogpdf([D], [U], [md], [mu], [vd], [vu], [cv])[0])
    if model == 'sat_single':
        p1 = float(occupancy(G, pars['q'], pars['gene_saturation']))
        bD = pars['shift_DEL']; bU = pars['shift_DUP']
        l0 = math.log(max(1 - p1, 1e-300)) + float(bvlogpdf([D], [U], [md], [mu], [vd], [vu], [cv])[0])
        l1 = math.log(max(p1, 1e-300)) + float(bvlogpdf([D], [U], [md + trait_sign * bD],
                                                        [mu + trait_sign * bU], [vd], [vu], [cv])[0])
        return float(logsumexp([l0, l1]))
    if model == 'aggregate':
        aD = pars['shift_DEL']; aU = pars['shift_DUP']
        return float(bvlogpdf([D], [U], [md + trait_sign * G * aD], [mu + trait_sign * G * aU],
                              [vd], [vu], [cv])[0])
    if model in CS_SHARED:
        qc = pars['cont_q']; bcD = pars['cont_shift_DEL']; bcU = pars['cont_shift_DUP']
        qs = pars['q']; bsD = pars['shift_DEL']; bsU = pars['shift_DUP']
        vv = G * (2 * qc - qc * qc)
        mdc = md + trait_sign * G * qc * bcD; muc = mu + trait_sign * G * qc * bcU
        vdc = vd + vv * bcD * bcD; vuc = vu + vv * bcU * bcU; cvc = cv + vv * bcD * bcU
        p0 = (1 - qs) ** G
        l0 = math.log(max(p0, 1e-300)) + float(bvlogpdf([D], [U], [mdc], [muc], [vdc], [vuc], [cvc])[0])
        l1 = math.log(max(1 - p0, 1e-300)) + float(bvlogpdf(
            [D], [U], [mdc + trait_sign * bsD], [muc + trait_sign * bsU], [vdc], [vuc], [cvc])[0])
        return float(logsumexp([l0, l1]))
    q = pars['q']; bD = pars['shift_DEL']; bU = pars['shift_DUP']
    if model == 'continuous':
        vv = G * (2 * q - q * q)
        return float(bvlogpdf([D], [U], [md + trait_sign * G * q * bD], [mu + trait_sign * G * q * bU],
                              [vd + vv * bD * bD], [vu + vv * bU * bU], [cv + vv * bD * bU])[0])
    if model == 'axis_continuous':
        w = pars['w_deletion_share']; qD = q * w; qU = q * (1 - w)
        vvD = G * (2 * qD - qD * qD); vvU = G * (2 * qU - qU * qU)
        return float(bvlogpdf([D], [U], [md + trait_sign * G * qD * bD], [mu + trait_sign * G * qU * bU],
                              [vd + vvD * bD * bD], [vu + vvU * bU * bU], [cv])[0])
    p0 = (1 - q) ** G
    if model == 'single':
        l0 = math.log(max(p0, 1e-300)) + float(bvlogpdf([D], [U], [md], [mu], [vd], [vu], [cv])[0])
        l1 = math.log(max(1 - p0, 1e-300)) + float(bvlogpdf([D], [U], [md + trait_sign * bD],
                                                            [mu + trait_sign * bU], [vd], [vu], [cv])[0])
        return float(logsumexp([l0, l1]))
    w = pars['w_deletion_share']
    a0 = math.log(max(p0, 1e-300)) + float(bvlogpdf([D], [U], [md], [mu], [vd], [vu], [cv])[0])
    aD = math.log(max((1 - p0) * w, 1e-300)) + float(bvlogpdf([D], [U], [md + trait_sign * bD], [mu], [vd], [vu], [cv])[0])
    aU = math.log(max((1 - p0) * (1 - w), 1e-300)) + float(bvlogpdf([D], [U], [md], [mu + trait_sign * bU], [vd], [vu], [cv])[0])
    return float(logsumexp([a0, aD, aU]))

def refit_fast(dat, base, trait_sign, model, xstart):
    cmp_ = mirror_components(dat, base)
    if model == 'sat_single':
        max_q = MAXQ[model]; amin, amax = ANGLE[model]
        r = minimize(lambda x: nll_sat(x, cmp_, trait_sign, model, max_q), np.asarray(xstart),
                     method='L-BFGS-B', bounds=[(-12, 12), (-9, 2), (-8, 8), (-2, 8)],
                     options={'maxiter': 1200})
        if not r.success:
            r = minimize(lambda x: nll_sat(x, cmp_, trait_sign, model, max_q),
                         np.asarray(xstart), method='Nelder-Mead', options={'maxiter': 5000})
        q, b, theta, gsat = unpack_sat(r.x, max_q, amin, amax)
        return dict(q=q, shift_DEL=float(b * np.cos(theta)), shift_DUP=float(b * np.sin(theta)),
                    gene_saturation=gsat)
    if model == 'aggregate':
        r = minimize(lambda x: nll_aggregate(x, cmp_, trait_sign, model, None),
                     np.asarray(xstart), method='L-BFGS-B', bounds=[(-20, 2), (-20, 2)],
                     options={'maxiter': 1200})
        if not r.success:
            r = minimize(lambda x: nll_aggregate(x, cmp_, trait_sign, model, None),
                         np.asarray(xstart), method='Nelder-Mead', options={'maxiter': 4000})
        return dict(shift_DEL=float(np.exp(r.x[0])), shift_DUP=float(np.exp(r.x[1])))
    if model in CS_SHARED:
        nb = 5 if CS_SHARED[model] else 6
        bd = [(-12, 12), (-9, 2), (-8, 8), (-12, 12), (-9, 2), (-8, 8)][:nb]
        r = minimize(lambda x: nll_cont_single(x, cmp_, trait_sign, model, MAXQ_CS),
                     np.asarray(xstart), method='L-BFGS-B', bounds=bd,
                     options={'maxiter': 2000})
        if not r.success:
            r = minimize(lambda x: nll_cont_single(x, cmp_, trait_sign, model, MAXQ_CS),
                         np.asarray(xstart), method='Nelder-Mead', options={'maxiter': 8000})
        return _cs_out(model, r.x, len(dat))
    max_q = MAXQ[model]
    fn = nll_axis if model.startswith('axis') else nll_common
    nb = 4 if model.startswith('axis') else 3
    bounds = [(-12, 12), (-12, 12), (-9, 2), (-9, 2)][:nb] if model.startswith('axis') \
        else [(-12, 12), (-9, 2), (-8, 8)]
    r = minimize(lambda x: fn(x, cmp_, trait_sign, model, max_q), np.asarray(xstart),
                 method='L-BFGS-B', bounds=bounds, options={'maxiter': 1200})
    if not r.success:
        r = minimize(lambda x: fn(x, cmp_, trait_sign, model, max_q), np.asarray(xstart),
                     method='Nelder-Mead', options={'maxiter': 4000})
    if model.startswith('axis'):
        q, w, bD, bU = unpack_axis(r.x, max_q)
        return dict(q=q, w_deletion_share=w, shift_DEL=bD, shift_DUP=bU)
    amin, amax = ANGLE[model]
    q = max_q * expit(r.x[0]); b = float(np.exp(r.x[1]))
    theta = np.deg2rad(amin + (amax - amin) * expit(r.x[2]))
    return dict(q=float(q), shift_DEL=float(b * np.cos(theta)), shift_DUP=float(b * np.sin(theta)))

# ------------------------------------------------- posterior predictive checks
STAT_NAMES = ['DEL_mean', 'DUP_mean', 'DEL_DUP_correlation', 'DEL_min', 'DEL_max',
              'DUP_min', 'DUP_max', 'DEL_negative_fraction', 'DUP_negative_fraction']

def empirical(dat):
    D = dat.Deletion_beta.to_numpy(float); U = dat.Duplication_beta.to_numpy(float)
    return np.array([D.mean(), U.mean(), np.corrcoef(D, U)[0, 1], D.min(), D.max(),
                     U.min(), U.max(), np.mean(D < 0), np.mean(U < 0)])

def simulate_model(dat, base, trait_sign, model, pars, nsim, seed):
    model = BASE.get(model, model)
    """Vectorised replicate simulation under a fitted model.

    Equivalent in distribution to the shipped per-draw loop, but samples all
    replicates at once via a per-locus Cholesky factor, so the RNG stream differs
    and PPC values match the shipped ones only to Monte Carlo error.

    The discriminating structure lives in the covariance: a shared diagonal shift
    inflates cov by +vv*outer([bD,bU]) and so pulls the DEL-DUP correlation toward
    zero; axis-aligned events leave cov untouched.
    """
    rg = np.random.default_rng(seed)
    D, U, G, md, mu, vd, vu, cv = mirror_components(dat, base)
    n = len(dat)
    mD = np.tile(md, (nsim, 1)); mU = np.tile(mu, (nsim, 1))
    VD = np.tile(vd, (nsim, 1)); VU = np.tile(vu, (nsim, 1)); CV = np.tile(cv, (nsim, 1))
    s = trait_sign
    if model in CS_SHARED:
        qc = pars['cont_q']; bcD = pars['cont_shift_DEL']; bcU = pars['cont_shift_DUP']
        vv = G * (2 * qc - qc * qc)
        mD += s * G * qc * bcD; mU += s * G * qc * bcU
        VD += vv * bcD * bcD; VU += vv * bcU * bcU; CV += vv * bcD * bcU
        hit = rg.random((nsim, n)) < (1 - (1 - pars['q']) ** G)
        mD += s * pars['shift_DEL'] * hit; mU += s * pars['shift_DUP'] * hit
    elif model == 'continuous':
        q = pars['q']; bD = pars['shift_DEL']; bU = pars['shift_DUP']; vv = G * (2 * q - q * q)
        mD += s * G * q * bD; mU += s * G * q * bU
        VD += vv * bD * bD; VU += vv * bU * bU; CV += vv * bD * bU
    elif model == 'axis_continuous':
        q = pars['q']; w = pars['w_deletion_share']; bD = pars['shift_DEL']; bU = pars['shift_DUP']
        qD = q * w; qU = q * (1 - w)
        mD += s * G * qD * bD; mU += s * G * qU * bU
        VD += G * (2 * qD - qD * qD) * bD * bD; VU += G * (2 * qU - qU * qU) * bU * bU
        # no covariance term: no single event moves both coordinates
    elif model == 'aggregate':
        mD += s * G * pars['shift_DEL']; mU += s * G * pars['shift_DUP']
    elif model == 'sat_single':
        hit = rg.random((nsim, n)) < occupancy(G, pars['q'], pars['gene_saturation'])
        mD += s * pars['shift_DEL'] * hit; mU += s * pars['shift_DUP'] * hit
    elif model == 'single':
        hit = rg.random((nsim, n)) < (1 - (1 - pars['q']) ** G)
        mD += s * pars['shift_DEL'] * hit; mU += s * pars['shift_DUP'] * hit
    elif model == 'axis_sparse':
        hit = rg.random((nsim, n)) < (1 - (1 - pars['q']) ** G)
        on_del = rg.random((nsim, n)) < pars['w_deletion_share']
        mD += s * pars['shift_DEL'] * (hit & on_del)
        mU += s * pars['shift_DUP'] * (hit & ~on_del)
    a = np.sqrt(VD); b = CV / a; c = np.sqrt(np.maximum(VU - b * b, 1e-300))
    z1 = rg.standard_normal((nsim, n)); z2 = rg.standard_normal((nsim, n))
    ds = mD + a * z1; us = mU + b * z1 + c * z2
    dm = ds.mean(1); um = us.mean(1)
    dc = ds - dm[:, None]; uc = us - um[:, None]
    corr = (dc * uc).sum(1) / np.sqrt((dc * dc).sum(1) * (uc * uc).sum(1))
    stats = np.column_stack([dm, um, corr, ds.min(1), ds.max(1), us.min(1), us.max(1),
                             (ds < 0).mean(1), (us < 0).mean(1)])
    return stats, ds, us

def ppc_table(dat, base, trait_sign, fits, models, nsim, trait):
    obs = empirical(dat); rows = []; clouds = {}
    for i, key in enumerate(models):
        st, ds, us = simulate_model(dat, base, trait_sign, key, fits[key], nsim,
                                    SEED + 500 * i + (0 if trait == 'height' else 50))
        clouds[key] = (ds[:150], us[:150])
        for j, nm in enumerate(STAT_NAMES):
            rows.append({'model': LABELS[key], 'statistic': nm, 'observed': obs[j],
                         'simulated_mean': st[:, j].mean(),
                         'simulated_2.5%': np.quantile(st[:, j], .025),
                         'simulated_97.5%': np.quantile(st[:, j], .975),
                         'posterior_predictive_CDF': float(np.mean(st[:, j] <= obs[j]))})
    return pd.DataFrame(rows), clouds

# ------------------------------------------------------------------- burden
def burden_values(base, trait_sign, model, pars, Gvals):
    model = BASE.get(model, model)
    """Expected per-affected-gene DEL and DUP effect, propagated from the fitted
    gene-level model. Extends the shipped burden_values() to the axis-aligned and
    aggregate geometries.

    For the axis models only a fraction of events lands on each axis, so each
    coordinate carries its own share (w for DEL, 1-w for DUP). That is the whole
    point: a diagonal shift moves both burdens together, an axis-aligned mixture
    splits its total effect between them.
    """
    p, m, rho, mus = base['p'], base['m'], base['rho'], base['sign_mean']
    d0 = p * mus * m; u0 = -rho * d0
    G = np.asarray(Gvals, float)
    if model == 'mirror':
        return np.repeat(d0, len(G)), np.repeat(u0, len(G))
    if model == 'continuous':
        return (np.repeat(d0 + trait_sign * pars['q'] * pars['shift_DEL'], len(G)),
                np.repeat(u0 + trait_sign * pars['q'] * pars['shift_DUP'], len(G)))
    if model == 'axis_continuous':
        w = pars['w_deletion_share']
        return (np.repeat(d0 + trait_sign * pars['q'] * w * pars['shift_DEL'], len(G)),
                np.repeat(u0 + trait_sign * pars['q'] * (1 - w) * pars['shift_DUP'], len(G)))
    if model == 'aggregate':
        # Mean displacement is G*a, so the per-gene contribution is a itself.
        return (np.repeat(d0 + trait_sign * pars['shift_DEL'], len(G)),
                np.repeat(u0 + trait_sign * pars['shift_DUP'], len(G)))
    if model in CS_SHARED:
        # continuous layer contributes q_c*b_c per gene; the single hit contributes
        # b_s spread over the G genes of a locus that is hit
        occ = 1 - (1 - pars['q']) ** G
        return (d0 + trait_sign * (pars['cont_per_gene_DEL'] + pars['shift_DEL'] * occ / G),
                u0 + trait_sign * (pars['cont_per_gene_DUP'] + pars['shift_DUP'] * occ / G))
    if model == 'sat_single':
        occ = occupancy(G, pars['q'], pars['gene_saturation'])
        return (d0 + trait_sign * pars['shift_DEL'] * occ / G,
                u0 + trait_sign * pars['shift_DUP'] * occ / G)
    occ = 1 - (1 - pars['q']) ** G
    if model == 'axis_sparse':
        w = pars['w_deletion_share']
        return (d0 + trait_sign * w * pars['shift_DEL'] * occ / G,
                u0 + trait_sign * (1 - w) * pars['shift_DUP'] * occ / G)
    return (d0 + trait_sign * pars['shift_DEL'] * occ / G,
            u0 + trait_sign * pars['shift_DUP'] * occ / G)

def burden_table(dat, base, sign, fits, models, target, units, trait_sd_native):
    rows = []; G = dat.Gene_count.to_numpy(float)
    for key in models:
        d, u = burden_values(base, sign, key, fits[key], G)
        dm = float(d.mean()); um = float(u.mean())
        d_nat = dm * trait_sd_native; u_nat = um * trait_sd_native
        rows.append({
            'model': LABELS[key],
            'DEL_effect_per_affected_gene_SD': dm,
            'DUP_effect_per_affected_gene_SD': um,
            'DEL_effect_per_affected_gene_native': d_nat,
            'DUP_effect_per_affected_gene_native': u_nat,
            'Mace_DEL_target_per_Mb': target, 'Mace_units': units,
            'effective_dosage_genes_per_Mb_to_match_DEL': target / d_nat if d_nat != 0 else np.nan,
            'predicted_DUP_effect_per_Mb_at_matched_density': um / dm * target if dm != 0 else np.nan,
            'DUP_burden_direction': ('adverse' if (um / dm * target) * (1 if sign > 0 else -1) > 0
                                     else 'protective') if dm != 0 else 'undefined',
        })
    return pd.DataFrame(rows)

BURDEN = {'height': (-0.41, 'cm/Mb', 7.99), 'bmi': (0.14, 'kg/m2/Mb', 4.7)}

# ---------------------------------------------------------------------- driver
MODELS = ['mirror', 'continuous', 'single', 'axis_continuous', 'axis_sparse',
          'cont_single', 'cont_single_shared',
          'open_continuous', 'open_single', 'sat_single', 'dup_only', 'aggregate']
TRAITS = [
    ('height', 'sacks_height_all_loci_exact.csv', 'height_mirror_parameters.csv', -1),
    ('bmi', 'sacks_bmi_complete_reciprocal_loci_exact.csv', 'bmi_mirror_parameters.csv', +1),
]

all_tables = {}
for trait, loci_file, par_file, sign in TRAITS:
    dat = pd.read_csv(PROJECT / trait / 'inputs' / loci_file)
    base = pd.read_csv(PROJECT / trait / 'outputs' / par_file).iloc[0].to_dict()
    print(f'\n=== {trait.upper()}  n={len(dat)} loci   mirror: '
          f"p={base['p']:.4f} m={base['m']:.4f} rho={base['rho']:.4f} ===", flush=True)

    fits = {m: fit_model(dat, base, sign, m, args.nstarts) for m in MODELS}

    # Validation: the two shipped models must be recovered by this re-implementation.
    shipped = pd.read_csv(PROJECT / trait / 'outputs' / f'{trait}_single_locus_model_comparison.csv')
    for key in ['continuous', 'single']:
        want = shipped.loc[shipped.model == LABELS[key], 'negative_log_likelihood']
        if len(want):
            d = abs(float(want.iloc[0]) - fits[key]['negative_log_likelihood'])
            flag = 'OK' if d < 1e-3 else 'MISMATCH'
            print(f'  [validate] {LABELS[key]:<40s} dNLL vs shipped = {d:.2e}  {flag}')

    rows = []
    for m in MODELS:
        f = fits[m]
        rows.append({k: f.get(k, np.nan) for k in
                     ['model', 'q', 'w_deletion_share', 'gene_saturation',
                      'shift_DEL', 'shift_DUP', 'shift_angle_deg',
                      'cont_q', 'cont_angle_deg', 'cont_magnitude',
                      'cont_per_gene_DEL', 'cont_per_gene_DUP',
                      'negative_log_likelihood',
                      'n_parameters_added', 'conditional_AIC', 'conditional_BIC']})
    tbl = pd.DataFrame(rows)

    if not args.skip_loo:
        loo_rows = []
        for m in MODELS:
            tot = 0.
            for i in range(len(dat)):
                if m == 'mirror':
                    pars = None
                else:
                    pars = refit_fast(dat.drop(dat.index[i]), base, sign, m, fits[m]['raw_par'])
                v = logpred_row(dat.iloc[i], base, sign, m, pars)
                tot += v
                loo_rows.append({'model': LABELS[m], 'held_out_locus': dat.iloc[i].get('CNV_ID', str(i)),
                                 'held_out_log_predictive_density': v})
            print(f'  [loo] {LABELS[m]:<45s} {tot:9.4f}', flush=True)
        loo = pd.DataFrame(loo_rows)
        sums = loo.groupby('model', as_index=False).held_out_log_predictive_density.sum() \
                  .rename(columns={'held_out_log_predictive_density': 'conditional_LOO_log_predictive_density'})
        tbl = tbl.merge(sums, on='model', how='left')
        loo.to_csv(PROJECT / trait / args.out_subdir / f'{trait}_axis_mixture_conditional_LOO{args.tag}.csv', index=False)

    tbl['delta_conditional_AIC'] = tbl.conditional_AIC - tbl.conditional_AIC.min()
    tbl['delta_conditional_BIC'] = tbl.conditional_BIC - tbl.conditional_BIC.min()
    if 'conditional_LOO_log_predictive_density' in tbl:
        tbl['delta_conditional_LOO'] = (tbl.conditional_LOO_log_predictive_density.max()
                                        - tbl.conditional_LOO_log_predictive_density)
    out = PROJECT / trait / args.out_subdir / f'{trait}_axis_mixture_model_comparison{args.tag}.csv'
    tbl.to_csv(out, index=False)
    all_tables[trait] = tbl

    show = ['model', 'q', 'w_deletion_share', 'shift_DEL', 'shift_DUP',
            'conditional_AIC', 'delta_conditional_AIC', 'delta_conditional_BIC']
    if 'delta_conditional_LOO' in tbl: show.append('delta_conditional_LOO')
    pd.set_option('display.width', 240, 'display.max_colwidth', 48)
    print(tbl[show].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'  wrote {out}')

    # ---- genome-wide burden propagation, all six models ----
    target, units, sd_native = BURDEN[trait]
    bur = burden_table(dat, base, sign, fits, MODELS, target, units, sd_native)
    bout = PROJECT / trait / args.out_subdir / f'{trait}_axis_mixture_burden_predictions{args.tag}.csv'
    bur.to_csv(bout, index=False)
    print(f'\n  --- burden propagation (Mace target {target} {units}) ---')
    print(bur[['model', 'DEL_effect_per_affected_gene_SD', 'DUP_effect_per_affected_gene_SD',
               'effective_dosage_genes_per_Mb_to_match_DEL',
               'predicted_DUP_effect_per_Mb_at_matched_density', 'DUP_burden_direction']]
          .to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'  wrote {bout}')

    # ---- posterior predictive checks, all six models ----
    ppc, clouds = ppc_table(dat, base, sign, fits, MODELS, 1200, trait)
    pout = PROJECT / trait / args.out_subdir / f'{trait}_axis_mixture_posterior_predictive_checks{args.tag}.csv'
    ppc.to_csv(pout, index=False)
    key_stats = ['DEL_mean', 'DUP_mean', 'DEL_DUP_correlation']
    piv = ppc[ppc.statistic.isin(key_stats)].pivot(index='model', columns='statistic',
                                                   values='posterior_predictive_CDF')
    piv = piv.reindex([LABELS[m] for m in MODELS])[key_stats]
    print('\n  --- posterior predictive CDF (0.5 ideal; <0.025 or >0.975 = failure) ---')
    print(piv.to_string(float_format=lambda v: f'{v:.3f}'))
    obs_corr = float(np.corrcoef(dat.Deletion_beta, dat.Duplication_beta)[0, 1])
    sim_corr = ppc[ppc.statistic == 'DEL_DUP_correlation'].set_index('model').simulated_mean
    print(f'\n  observed DEL-DUP correlation = {obs_corr:.3f}; simulated means:')
    for m in MODELS:
        print(f'    {LABELS[m]:<50s} {sim_corr[LABELS[m]]:+.3f}')
    print(f'  wrote {pout}')
    np.savez(PROJECT / trait / args.out_subdir / f'{trait}_ppc_clouds{args.tag}.npz',
             **{f'{k}_{ax}': v[i] for k, v in clouds.items() for i, ax in enumerate(['DEL', 'DUP'])})

print('\nDone.')


