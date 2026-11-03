# Full grid of two-component architectures (JS 2026-11-03), both traits, same base/likelihood as run_corrected_analysis.py.
# Mirror part:     M0 background (gene-level, as Fig. 2); MT + TWAS along the mirror direction; MTH + TWAS + Poisson mirror hits
#                  (rate per constrained gene, s_het >= 0.03).
# Asymmetry part (trait-signed perpendicular direction: lower-left height, upper-right BMI):
#                  none; continuous(genes); continuous(constrained); single hit (P = 1-(1-q)^G); multi-hit(genes);
#                  multi-hit(constrained) (Poisson number of perpendicular hits); single hit + continuous(constrained).
# 3 x 7 = 21 models per trait, conditional AIC/BIC with the mirror base fixed. Output: ../py_outputs/models_grid.csv
import math, itertools
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import poisson
src = open(__file__.replace("models_grid.py", "models_twas_burden.py"), encoding="utf-8").read().split("rows = []")[0]
exec(src)
K = 10
def fit(dat, base, sign, twas, mhit, asym, nstarts=16):
    D, U, G, md0, mu0, vd, vu, cv = components(dat, base); rho = base["rho"]; sm = np.sign(base["sign_mean"])
    Z = dat.Z.to_numpy(float); Zs = Z / Z.std(); B = dat.B.to_numpy(float); Gc = dat.Gene_count.to_numpy(float)
    mdir = sm * np.array([1.0, -rho]) / math.hypot(1, rho)
    cont = {"cont_G": Gc, "cont_B": B, "single_cont_B": B}.get(asym); multi = {"multi_G": Gc, "multi_B": B}.get(asym)
    single = asym in ("single", "single_cont_B")
    lay = (["c"] if twas else []) + (["lh", "h"] if mhit else []) + (["qc", "bc", "tc"] if cont is not None else []) \
          + (["qs", "bs", "ts"] if single else []) + (["la", "ba", "ta"] if multi is not None else [])
    def un(x):
        p = dict(zip(lay, x)); o = dict(p)
        for k in ("lh", "h", "bc", "bs", "la", "ba"):
            if k in o: o[k] = np.exp(p[k])
        if "qc" in o: o["qc"] = .8 * expit(p["qc"])
        if "qs" in o: o["qs"] = .2 * expit(p["qs"])
        for k in ("tc", "ts", "ta"):
            if k in o: o[k] = np.deg2rad(5 + 80 * expit(p[k]))
        return o
    def nll(x):
        o = un(x); md = md0.copy(); mu = mu0.copy(); vdd = vd.copy(); vuu = vu.copy(); cvv = cv.copy()
        if twas: md = md - o["c"] * Zs; mu = mu + rho * o["c"] * Zs
        if cont is not None:
            bD, bU = o["bc"] * np.cos(o["tc"]), o["bc"] * np.sin(o["tc"]); q = o["qc"]; vv = cont * (2 * q - q * q)
            md = md + sign * cont * q * bD; mu = mu + sign * cont * q * bU; vdd = vdd + vv * bD * bD; vuu = vuu + vv * bU * bU; cvv = cvv + vv * bD * bU
        mk = np.arange(K + 1) if mhit else np.array([0])
        ak = np.arange(K + 1) if multi is not None else (np.array([0, 1]) if single else np.array([0]))
        if single: sD, sU = o["bs"] * np.cos(o["ts"]), o["bs"] * np.sin(o["ts"]); p0 = (1 - o["qs"]) ** G
        if multi is not None: aD, aU = o["ba"] * np.cos(o["ta"]), o["ba"] * np.sin(o["ta"])
        cols = []
        for k, a in itertools.product(mk, ak):
            m1 = md + (k * o["h"] * mdir[0] if mhit else 0); m2 = mu + (k * o["h"] * mdir[1] if mhit else 0); lp = np.zeros(len(D))
            if mhit: lp = lp + poisson.logpmf(k, o["lh"] * B)
            if single:
                m1 = m1 + a * sign * sD; m2 = m2 + a * sign * sU; lp = lp + (np.log(np.maximum(1 - p0, 1e-300)) if a else np.log(np.maximum(p0, 1e-300)))
            if multi is not None:
                m1 = m1 + a * sign * aD; m2 = m2 + a * sign * aU; lp = lp + poisson.logpmf(a, o["la"] * multi)
            l = bvlogpdf(D, U, m1, m2, vdd, vuu, cvv)
            if l is None: return 1e100
            cols.append(l + lp)
        return float(-logsumexp(np.column_stack(cols), axis=1).sum())
    if not lay: return dict(nll=nll(np.array([])), k=0, pars={})
    st = {"c": 0.0, "lh": math.log(.01), "h": math.log(.3), "qc": math.log(.08 / .72), "bc": math.log(.03), "tc": 0.0,
          "qs": math.log(.018 / .182), "bs": math.log(.27), "ts": 2.0, "la": math.log(.03), "ba": math.log(.1), "ta": 0.0}
    x0 = np.array([st[n] for n in lay]); best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 733 * i)
        r = minimize(nll, x0 + (rg.normal(0, 1.2, len(x0)) if i else 0), method="Nelder-Mead", options={"maxiter": 30000, "xatol": 1e-8, "fatol": 1e-9})
        if best is None or r.fun < best.fun: best = r
    o = un(best.x)
    return dict(nll=float(best.fun), k=len(lay), pars={k: (float(np.rad2deg(v)) if k in ("tc", "ts", "ta") else float(v)) for k, v in o.items()})
MIRROR = [("background", False, False), ("background + TWAS", True, False), ("background + TWAS + mirror hits", True, True)]
ASYM = ["none", "cont_G", "cont_B", "single", "multi_G", "multi_B", "single_cont_B"]
rows = []
for trait, dat0, bdat, pi, sign, tw in [("height", height, height_mirror, .30, -1, "height"), ("BMI", bmi, bmi, 11 / 14, 1, "bmi")]:
    dat = add_cov(dat0, tw); base = fit_mirror(add_cov(bdat, tw), pi)
    for (mname, twas, mhit), asym in itertools.product(MIRROR, ASYM):
        f = fit(dat, base, sign, twas, mhit, asym)
        rows.append(dict(trait=trait, mirror=mname, asymmetry=asym, k=f["k"], nll=f["nll"], cAIC=2 * f["k"] + 2 * f["nll"],
                         cBIC=f["k"] * math.log(len(dat)) + 2 * f["nll"], **{f"par_{k}": v for k, v in f["pars"].items()}))
        print(trait, mname, asym, round(f["nll"], 3), flush=True)
R = pd.DataFrame(rows); R["dAIC"] = R.groupby("trait").cAIC.transform(lambda v: v - v.min()); R["dBIC"] = R.groupby("trait").cBIC.transform(lambda v: v - v.min())
R.to_csv(PKG / ("py_outputs/models_grid.csv" if TWAS_SOURCE == "mancuso" else f"py_outputs/models_grid_{TWAS_SOURCE}.csv"), index=False)
for t in ["height", "BMI"]:
    print(f"\n{t}: delta AIC (rows mirror part, columns asymmetry part)")
    print(R[R.trait == t].pivot(index="mirror", columns="asymmetry", values="dAIC").loc[[m[0] for m in MIRROR], ASYM].round(1).to_string())
    print(f"{t}: delta BIC"); print(R[R.trait == t].pivot(index="mirror", columns="asymmetry", values="dBIC").loc[[m[0] for m in MIRROR], ASYM].round(1).to_string())
