# Directional dose-response models with TWAS and constrained-gene burden (JS 2026-11-03).
# Extends run_corrected_analysis.py (Fig. 2K-O): same data, same mirror base (gene-level random effects with mean and variance
# proportional to gene count, measurement error, mirror slope rho), same trait-signed perpendicular shift (lower-left height,
# upper-right BMI) and same conditional-AIC comparison. Added terms:
#   T  TWAS: mean shift along the trait's mirror direction, c x summed TWAS z (DEL -c*Z, DUP +rho*c*Z); 1 parameter.
#   C(X) continuous directional shift proportional to X (X = gene count G, as in the original, or constrained-gene count B,
#        GeneBayes s_het >= 0.03): mean += sign*X*q*(bD,bU), variance += X*(2q - q^2)*b b'; 3 parameters (q, b, angle).
#   S  sparse single hit: mixture with P(hit) = 1 - (1 - q)^G, shift sign*(bD,bU); 3 parameters.
#   S + C(B): both, separate magnitudes and angles; 6 parameters.
# Each directional model is fitted with and without T. Summed TWAS z per locus: contiguous run of TWAS genes spanning the CNV
# (Mancuso et al. 2017, mean z over panels), all loci of the model (overlapping intervals included), 0 if no TWAS gene.
#   H(X) multi-hit: Poisson number of hits with rate lambda x X (X = G or B), each hit a fixed shift sign*(bD,bU); 3 parameters.
#        Between single hit (rare, large) and continuous (many, tiny; Gaussian approximation).
# Output: ../py_outputs/models_twas_burden.csv
import math
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
SEED = 20260726
PKG = Path(__file__).resolve().parent.parent
CONC = PKG.parent.parent / "CNV_ASD_growth_analysis" / "concordance"
ROOT = PKG.parent.parent
height = pd.read_csv(PKG / "height/inputs/sacks_height_all_loci_exact.csv")
height_mirror = height[height.Dose_response_class.str.lower().eq("mirror")].copy().reset_index(drop=True)
bmi = pd.read_csv(PKG / "bmi/inputs/sacks_bmi_complete_reciprocal_loci_exact.csv")

# ---------------- locus covariates: summed TWAS z and constrained-gene count ----------------
key = lambda s: str(s).replace("-", "_").replace(".", "").lower()
genes = pd.read_csv(ROOT / "pathway_analysis/Final/cnv_genes_annotated.csv")
gb = pd.read_csv(CONC / "GeneBayes_shet_Zeng2024.tsv", sep="\t"); hg = pd.read_csv(CONC / "hgnc_map.tsv", sep="\t")
gb["sym"] = gb.hgnc.map(dict(zip(hg.hgnc_id, hg.symbol))); shet = gb.dropna(subset=["sym"]).groupby("sym").post_mean.max()
gl = {key(c): sub.Gene.unique() for c, sub in genes.groupby("CNV_ID")}
import os, glob as _glob
TWAS_SOURCE = os.environ.get("TWAS_SOURCE", "mancuso")   # "mancuso" (Fig. 2E) or "brain" (PsychENCODE cortex, Tan et al. 2024 GWAS)
def twas_sums(trait):
    if TWAS_SOURCE == "brain":
        fs = sorted(_glob.glob(str(CONC / "psychencode_twas" / "raw" / f"{trait.upper()}_TAN" / "PECfull.*.dat")))
        t = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
        for c in ("TWAS.Z", "P0", "CHR"): t[c] = pd.to_numeric(t[c], errors="coerce")
        t = t.dropna(subset=["TWAS.Z", "P0", "CHR"]); e2s = dict(zip(hg.ensembl_gene_id, hg.symbol))
        t["GENE"] = t.ID.str.split(".").str[0].map(e2s); t = t.dropna(subset=["GENE"])
        tw = t.rename(columns={"P0": "TxS"})
    else:
        tw = pd.read_csv(CONC / f"twas30/{trait}.txt", sep=r"\s+")
    chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum) if TWAS_SOURCE != "brain" else g["chr"].astype(int); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
    cs = np.concatenate([[0], np.cumsum(g.z.values)]); pos_of = dict(zip(g.GENE, g.index)); out = {}
    for k, gs in gl.items():
        ix = sorted(pos_of[x] for x in gs if x in pos_of)
        out[k] = cs[ix[-1] + 1] - cs[ix[0]] if ix else 0.0
    return out
def add_cov(dat, trait):
    d = dat.copy(); zs = twas_sums(trait); idcol = "CNV_ID"
    d["Z"] = [zs.get(key(c), np.nan) for c in d[idcol]]
    d["B"] = [int(np.nansum(np.array([shet.get(x, np.nan) for x in gl.get(key(c), [])]) >= 0.03)) for c in d[idcol]]
    miss = d[d.Z.isna() | ~d[idcol].map(lambda c: key(c) in gl)]
    if len(miss): print(f"[{trait}] loci without gene list:", miss[idcol].tolist())
    d["Z"] = d.Z.fillna(0.0)
    return d

# ---------------- base mirror model (verbatim from run_corrected_analysis.py) ----------------
def bvlogpdf(x, y, m1, m2, v1, v2, cov):
    det = v1 * v2 - cov * cov
    if np.any(det <= 1e-14) or np.any(~np.isfinite(det)): return None
    dx = x - m1; dy = y - m2
    return -np.log(2 * np.pi) - .5 * np.log(det) - .5 * (v2 * dx * dx + v1 * dy * dy - 2 * cov * dx * dy) / det
def fit_mirror(dat, pi_positive, nstarts=10):
    D = dat.Deletion_beta.to_numpy(float); U = dat.Duplication_beta.to_numpy(float)
    G = dat.Gene_count.to_numpy(float); seD = dat.Deletion_SE.to_numpy(float); seU = dat.Duplication_SE.to_numpy(float); mus = 2 * pi_positive - 1
    def unpack(x): return dict(p=float(expit(x[0])), m=float(np.exp(x[1])), rho=float(np.exp(x[2])), tauD=float(np.exp(x[3])), tauU=float(np.exp(x[4])), sign_mean=mus)
    def nll(x):
        z = unpack(x); p, m, rho, tD, tU = z["p"], z["m"], z["rho"], z["tauD"], z["tauU"]
        mg = p * mus * m; vg = p * 2 * m * m - (p * mus * m) ** 2
        if vg <= 1e-12 or not np.isfinite(vg): return 1e100
        v = G * vg; ll = bvlogpdf(D, U, G * mg, -rho * G * mg, v + seD ** 2 + tD ** 2, rho * rho * v + seU ** 2 + tU ** 2, -rho * v)
        return 1e100 if ll is None else float(-ll.sum())
    p0 = .4 if pi_positive < .5 else .03; m0 = max(.02, abs(D.mean()) / (max(abs(mus), .1) * G.mean() * p0))
    rho0 = max(.15, min(1.0, abs(np.polyfit(D, U, 1)[0])))
    x0 = np.array([math.log(p0 / (1 - p0)), math.log(m0), math.log(rho0), math.log(.02), math.log(.05)]); best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 107 * i + int(pi_positive * 1000))
        r = minimize(nll, x0 + rg.normal(0, [1.4, 1, .6, 1, 1]), method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-8, "fatol": 1e-8})
        if best is None or r.fun < best.fun: best = r
    return unpack(best.x)
def components(dat, z):
    D = dat.Deletion_beta.to_numpy(float); U = dat.Duplication_beta.to_numpy(float); G = dat.Gene_count.to_numpy(float)
    seD = dat.Deletion_SE.to_numpy(float); seU = dat.Duplication_SE.to_numpy(float)
    mg = z["p"] * z["sign_mean"] * z["m"]; vg = max(z["p"] * 2 * z["m"] ** 2 - mg ** 2, 1e-12); v = G * vg
    return D, U, G, G * mg, -z["rho"] * G * mg, v + seD ** 2 + z["tauD"] ** 2, z["rho"] ** 2 * v + seU ** 2 + z["tauU"] ** 2, -z["rho"] * v

# ---------------- extended directional models ----------------
def fit_model(dat, base, sign, cont=None, single=False, twas=False, multi=None, nstarts=16):
    D, U, G, md0, mu0, vd, vu, cv = components(dat, base); rho = base["rho"]
    Z = dat.Z.to_numpy(float); X = None if cont is None else dat[cont].to_numpy(float)
    Zs = Z / Z.std() if Z.std() > 0 else Z
    M = None if multi is None else dat[multi].to_numpy(float)
    layout = (["c"] if twas else []) + (["qc", "bc", "tc"] if cont else []) + (["qs", "bs", "ts"] if single else []) + (["lm", "bm", "tm"] if multi else [])
    def unpack(x):
        p = dict(zip(layout, x)); o = {}
        if twas: o["c"] = p["c"]
        if cont: o["qc"] = .8 * expit(p["qc"]); o["bc"] = np.exp(p["bc"]); o["tc"] = np.deg2rad(5 + 80 * expit(p["tc"]))
        if single: o["qs"] = .2 * expit(p["qs"]); o["bs"] = np.exp(p["bs"]); o["ts"] = np.deg2rad(5 + 80 * expit(p["ts"]))
        if multi: o["lm"] = np.exp(p["lm"]); o["bm"] = np.exp(p["bm"]); o["tm"] = np.deg2rad(5 + 80 * expit(p["tm"]))
        return o
    def nll(x):
        o = unpack(x); md = md0.copy(); mu = mu0.copy(); vdd = vd.copy(); vuu = vu.copy(); cvv = cv.copy()
        if twas: md = md - o["c"] * Zs; mu = mu + rho * o["c"] * Zs
        if cont:
            bD, bU = o["bc"] * np.cos(o["tc"]), o["bc"] * np.sin(o["tc"]); q = o["qc"]; vv = X * (2 * q - q * q)
            md = md + sign * X * q * bD; mu = mu + sign * X * q * bU; vdd = vdd + vv * bD * bD; vuu = vuu + vv * bU * bU; cvv = cvv + vv * bD * bU
        if multi:   # Poisson number of hits, rate lm per gene (or per constrained gene); each hit shifts sign*(bD,bU)
            bD, bU = o["bm"] * np.cos(o["tm"]), o["bm"] * np.sin(o["tm"]); lam = o["lm"] * M; ks = np.arange(0, 25)
            from scipy.stats import poisson
            lp = poisson.logpmf(ks[None, :], lam[:, None]); cols = []
            for kk in ks:
                l = bvlogpdf(D, U, md + sign * kk * bD, mu + sign * kk * bU, vdd, vuu, cvv)
                if l is None: return 1e100
                cols.append(l)
            return float(-logsumexp(lp + np.column_stack(cols), axis=1).sum())
        if not single:
            ll = bvlogpdf(D, U, md, mu, vdd, vuu, cvv); return 1e100 if ll is None else float(-ll.sum())
        bD, bU = o["bs"] * np.cos(o["ts"]), o["bs"] * np.sin(o["ts"]); p0 = (1 - o["qs"]) ** G
        l0 = bvlogpdf(D, U, md, mu, vdd, vuu, cvv); l1 = bvlogpdf(D, U, md + sign * bD, mu + sign * bU, vdd, vuu, cvv)
        if l0 is None or l1 is None: return 1e100
        return float(-logsumexp(np.column_stack([np.log(np.maximum(p0, 1e-300)) + l0, np.log(np.maximum(1 - p0, 1e-300)) + l1]), axis=1).sum())
    if not layout:
        ll = bvlogpdf(D, U, md0, mu0, vd, vu, cv); return dict(nll=float(-ll.sum()), k=0, pars={})
    x0 = []
    for nm in layout:
        x0.append({"lm": math.log(.02), "bm": math.log(.08), "tm": 0.0, "c": 0.0, "qc": math.log(.08 / .72), "bc": math.log(.055), "tc": 0.0, "qs": math.log(.018 / .182), "bs": math.log(.25), "ts": 0.0}[nm])
    x0 = np.array(x0); best = None
    for i in range(nstarts):
        rg = np.random.default_rng(SEED + 313 * i)
        r = minimize(nll, x0 + (rg.normal(0, 1.2, len(x0)) if i else 0), method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-8, "fatol": 1e-9})
        if best is None or r.fun < best.fun: best = r
    o = unpack(best.x)
    pars = {k: (float(np.rad2deg(v)) if k in ("tc", "ts", "tm") else float(v)) for k, v in o.items()}
    return dict(nll=float(best.fun), k=len(layout), pars=pars)

rows = []
for trait, dat0, base_dat, pi, sign, tw in [("height", height, height_mirror, .30, -1, "height"), ("BMI", bmi, bmi, 11 / 14, 1, "bmi")]:
    dat = add_cov(dat0, tw); base_dat = add_cov(base_dat, tw)
    base = fit_mirror(base_dat, pi)
    print(f"\n== {trait}: n={len(dat)} loci; base rho={base['rho']:.3f}; corr(Z, B)={np.corrcoef(dat.Z, dat.B)[0,1]:+.2f}, corr(G, B)={np.corrcoef(dat.Gene_count, dat.B)[0,1]:+.2f}")
    specs = [("Mirror", None, False), ("Mirror + continuous(G)", "Gene_count", False), ("Mirror + continuous(constrained)", "B", False),
             ("Mirror + single hit", None, True), ("Mirror + single hit + continuous(constrained)", "B", True),
             ("Mirror + multi-hit(G)", None, False, "Gene_count"), ("Mirror + multi-hit(constrained)", None, False, "B")]
    for spec in specs:
        name, cont, single = spec[:3]; multi = spec[3] if len(spec) > 3 else None
        for twas in (False, True):
            f = fit_model(dat, base, sign, cont=cont, single=single, twas=twas, multi=multi)
            rows.append(dict(trait=trait, model=name + (" + TWAS" if twas else ""), k=f["k"], nll=f["nll"], cAIC=2 * f["k"] + 2 * f["nll"],
                             cBIC=f["k"] * math.log(len(dat)) + 2 * f["nll"], **{f"par_{k}": v for k, v in f["pars"].items()}))
R = pd.DataFrame(rows)
R["dAIC"] = R.groupby("trait").cAIC.transform(lambda v: v - v.min()); R["dBIC"] = R.groupby("trait").cBIC.transform(lambda v: v - v.min())
(PKG / "py_outputs").mkdir(exist_ok=True); R.to_csv(PKG / "py_outputs/models_twas_burden.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R.sort_values(["trait", "cAIC"]).round(3).to_string(index=False))
