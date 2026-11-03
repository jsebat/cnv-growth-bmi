# Sparse mirror hits (JS 2026-11-03): the 2 x 2 frame (axis: along the mirror line / perpendicular; source: diffuse / sparse hits).
# Adds cell 3, mirror hits: a Poisson number of hits per CNV (rate lambda per gene, or per constrained gene with s_het >= 0.03),
# each moving the locus by h along its trait's mirror direction (more copies -> taller for height; for BMI the fitted mirror
# direction, deletion raising BMI). Fitted for both traits on top of each trait's components:
#   height: mirror base + TWAS + asymmetric single hit (lower-left), with/without mirror hits
#   BMI:    mirror base + TWAS + continuous asymmetry driven by constrained genes (upper-right), with/without mirror hits
# plus mirror hits alone. Same base, likelihood and conditional AIC as models_twas_burden.py / run_corrected_analysis.py.
# Output: ../py_outputs/models_mirror_hits.csv, ../py_outputs/models_mirror_hits_posterior.csv
import math
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import poisson
src = open(__file__.replace("models_mirror_hits.py", "models_twas_burden.py"), encoding="utf-8").read().split("rows = []")[0]
exec(src)
K = 12
def fit(dat, base, sign, twas=True, asym=None, cont=None, mhit=None, nstarts=20, return_post=False):
    """asym: 'single' -> perpendicular single-hit mixture (P = 1-(1-q)^G); cont: column driving continuous asymmetry;
    mhit: column driving Poisson mirror hits (None = no mirror hits)."""
    D, U, G, md0, mu0, vd, vu, cv = components(dat, base); rho = base["rho"]; sm = np.sign(base["sign_mean"])
    Z = dat.Z.to_numpy(float); Zs = Z / Z.std()
    md0 = md0.copy(); mu0 = mu0.copy()
    mdir = sm * np.array([1.0, -rho]) / math.hypot(1, rho)          # direction of a growing mirror dose response
    lay = (["c"] if twas else []) + (["qc", "bc", "tc"] if cont else []) + (["qs", "bs", "ts"] if asym else []) + (["lh", "h"] if mhit else [])
    Xc = dat[cont].to_numpy(float) if cont else None; Xm = dat[mhit].to_numpy(float) if mhit else None
    def un(x):
        p = dict(zip(lay, x)); o = {}
        if twas: o["c"] = p["c"]
        if cont: o["qc"] = .8 * expit(p["qc"]); o["bc"] = np.exp(p["bc"]); o["tc"] = np.deg2rad(5 + 80 * expit(p["tc"]))
        if asym: o["qs"] = .2 * expit(p["qs"]); o["bs"] = np.exp(p["bs"]); o["ts"] = np.deg2rad(5 + 80 * expit(p["ts"]))
        if mhit: o["lh"] = np.exp(p["lh"]); o["h"] = np.exp(p["h"])
        return o
    def comps(o):
        md = md0.copy(); mu = mu0.copy(); vdd = vd.copy(); vuu = vu.copy(); cvv = cv.copy()
        if twas: md = md - o["c"] * Zs; mu = mu + rho * o["c"] * Zs
        if cont:
            bD, bU = o["bc"] * np.cos(o["tc"]), o["bc"] * np.sin(o["tc"]); q = o["qc"]; vv = Xc * (2 * q - q * q)
            md = md + sign * Xc * q * bD; mu = mu + sign * Xc * q * bU; vdd = vdd + vv * bD * bD; vuu = vuu + vv * bU * bU; cvv = cvv + vv * bD * bU
        cols, lps, labels = [], [], []
        asym_states = [(0, None)] if not asym else [(0, None), (1, None)]
        ks = np.arange(K + 1) if mhit else np.array([0])
        for a, _ in asym_states:
            for k in ks:
                m1, m2 = md.copy(), mu.copy()
                if a: m1 = m1 + sign * o["bs"] * np.cos(o["ts"]); m2 = m2 + sign * o["bs"] * np.sin(o["ts"])
                if mhit: m1 = m1 + k * o["h"] * mdir[0]; m2 = m2 + k * o["h"] * mdir[1]
                l = bvlogpdf(D, U, m1, m2, vdd, vuu, cvv)
                if l is None: return None
                lp = np.zeros(len(D))
                if asym: p0 = (1 - o["qs"]) ** G; lp = lp + (np.log(np.maximum(1 - p0, 1e-300)) if a else np.log(np.maximum(p0, 1e-300)))
                if mhit: lp = lp + poisson.logpmf(k, o["lh"] * Xm)
                cols.append(l + lp); labels.append((a, k))
        return np.column_stack(cols), labels
    def nll(x):
        r = comps(un(x))
        return 1e100 if r is None else float(-logsumexp(r[0], axis=1).sum())
    st = {"c": 0.0, "qc": math.log(.08 / .72), "bc": math.log(.03), "tc": 0.0, "qs": math.log(.018 / .182), "bs": math.log(.27), "ts": 2.0, "lh": math.log(.01), "h": math.log(.3)}
    x0 = np.array([st[n] for n in lay]); best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 577 * i)
        r = minimize(nll, x0 + (rg.normal(0, 1.2, len(x0)) if i else 0), method="Nelder-Mead", options={"maxiter": 30000, "xatol": 1e-8, "fatol": 1e-9})
        if best is None or r.fun < best.fun: best = r
    o = un(best.x); out = dict(nll=float(best.fun), k=len(lay), pars={k: (float(np.rad2deg(v)) if k in ("tc", "ts") else float(v)) for k, v in o.items()})
    if return_post and mhit:
        L, labels = comps(o); P = np.exp(L - logsumexp(L, axis=1, keepdims=True)); kk = np.array([l[1] for l in labels])
        out["post"] = pd.DataFrame(dict(CNV_ID=dat.CNV_ID.values, genes=G, constrained=dat.B.values, prior_mirror_hits=o["lh"] * Xm,
                                        posterior_mirror_hits=P @ kk, P_ge1=P[:, kk >= 1].sum(1)))
    return out
rows, posts = [], []
for trait, dat0, bdat, pi, sign, tw, extra in [("height", height, height_mirror, .30, -1, "height", dict(asym="single")),
                                               ("BMI", bmi, bmi, 11 / 14, 1, "bmi", dict(cont="B"))]:
    dat = add_cov(dat0, tw); base = fit_mirror(add_cov(bdat, tw), pi)
    for name, kw in [("background + TWAS", {}), ("background + TWAS + mirror hits (genes)", dict(mhit="Gene_count")),
                     ("background + TWAS + mirror hits (constrained)", dict(mhit="B")),
                     (f"background + TWAS + {'asymmetric single hit' if trait == 'height' else 'continuous asymmetry (constrained)'}", extra),
                     (f"background + TWAS + {'asymmetric single hit' if trait == 'height' else 'continuous asymmetry (constrained)'} + mirror hits (genes)", dict(extra, mhit="Gene_count")),
                     (f"background + TWAS + {'asymmetric single hit' if trait == 'height' else 'continuous asymmetry (constrained)'} + mirror hits (constrained)", dict(extra, mhit="B"))]:
        f = fit(dat, base, sign, return_post=True, **kw)
        rows.append(dict(trait=trait, model=name, k=f["k"], nll=f["nll"], cAIC=2 * f["k"] + 2 * f["nll"], cBIC=f["k"] * math.log(len(dat)) + 2 * f["nll"], **{f"par_{k}": v for k, v in f["pars"].items()}))
        if "post" in f: posts.append(f["post"].assign(trait=trait, model=name))
        print(trait, name, round(f["nll"], 3), flush=True)
R = pd.DataFrame(rows); R["dAIC"] = R.groupby("trait").cAIC.transform(lambda v: v - v.min()); R["dBIC"] = R.groupby("trait").cBIC.transform(lambda v: v - v.min())
R.to_csv(PKG / "py_outputs/models_mirror_hits.csv", index=False); PO = pd.concat(posts); PO.to_csv(PKG / "py_outputs/models_mirror_hits_posterior.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 90)
print(R.sort_values(["trait", "cAIC"]).round(3).to_string(index=False))
