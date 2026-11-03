# Constrained genes, TWAS and CNV dosage effects (Jon 2026-11-03). GeneBayes s_het (Zeng et al. 2024, Supplementary Table 1).
#  1. Genome-wide, per trait (Mancuso et al. 2017 TWAS, height and BMI): do constrained genes have TWAS models, and how strong
#     is their signal? TWAS coverage and mean |z| by s_het tier (< 0.03, 0.03-0.18, >= 0.18).
#  2. Per CNV (the Fig. 2E non-overlapping set): constrained-gene burden (sum of s_het; count of genes with s_het >= 0.03) and the
#     mean TWAS z of constrained vs unconstrained genes.
#  3. Models of the dosage effect, (DUP - DEL)/2: burden alone, summed TWAS z alone, both (the burden term absorbs a uniform lean
#     carried by dosage-sensitive genes; the TWAS term is then the locus-specific part). Least squares on standardised variables
#     and the rank (Spearman) version, with and without Williams. Output: twas_constraint_{genomewide,per_cnv,models}.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
gb = pd.read_csv(HERE + r"\GeneBayes_shet_Zeng2024.tsv", sep="\t")
hg = pd.read_csv(HERE + r"\hgnc_map.tsv", sep="\t")
sym = dict(zip(hg.hgnc_id, hg.symbol)); gb["symbol"] = gb.hgnc.map(sym)
shet = gb.dropna(subset=["symbol"]).groupby("symbol").post_mean.max()
# aliases and previous symbols -> s_het (for older gene names in the TWAS and CNV tables)
alias = {}
for _, r in hg.iterrows():
    if r.symbol in shet.index:
        for col in ("prev_symbol", "alias_symbol"):
            if isinstance(r[col], str):
                for a in r[col].split("|"): alias.setdefault(a, r.symbol)
def s_of(g): return shet.get(g, shet.get(alias.get(g, ""), np.nan))
tier = lambda s: np.where(np.isnan(s), "no s_het", np.where(s >= 0.18, ">=0.18", np.where(s >= 0.03, "0.03-0.18", "<0.03")))
# ---- 1. genome-wide
GW = []
tw = {}
for tr in ["height", "bmi"]:
    t = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+"); z = t.groupby("GENE")["TWAS.Z"].mean(); tw[tr] = z
    zs = pd.DataFrame(dict(z=z)); zs["s"] = [s_of(g) for g in zs.index]; zs["tier"] = tier(zs.s.values)
    univ = pd.DataFrame(dict(s=shet)); univ["tier"] = tier(univ.s.values); univ["has_twas"] = univ.index.isin(z.index) | univ.index.map(lambda g: False)
    for T in ["<0.03", "0.03-0.18", ">=0.18"]:
        u = univ[univ.tier == T]; k = zs[zs.tier == T]
        GW.append(dict(trait=tr, shet_tier=T, genes_with_shet=len(u), frac_with_twas_model=u.has_twas.mean(), n_twas=len(k),
                       mean_abs_z=k.z.abs().mean(), mean_z=k.z.mean(), frac_abs_z_gt_4p5=(k.z.abs() > 4.5).mean()))
    r, p = stats.spearmanr(zs.dropna().s, zs.dropna().z.abs()); GW.append(dict(trait=tr, shet_tier="Spearman |z| vs s_het", mean_abs_z=r, n_twas=int(zs.s.notna().sum()), mean_z=p))
GW = pd.DataFrame(GW); GW.to_csv(HERE + r"\twas_constraint_genomewide.csv", index=False)
pd.set_option("display.width", 250); print(GW.round(4).to_string(index=False))
# ---- 2. per CNV
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
hd = pd.read_csv(HERE + r"\twas_dosage_panel_nonoverlap.csv")[["CNV_ID", "cls", "summed_twas_z", "dosage_effect"]].rename(columns={"summed_twas_z": "sumz_height", "dosage_effect": "dos_height"})
bd = pd.read_csv(HERE + r"\twas_dosage_bmi.csv").rename(columns={"summed_bmi_twas_z": "sumz_bmi", "dosage_effect_bmi": "dos_bmi"})[["CNV_ID", "sumz_bmi", "dos_bmi"]]
P = hd.merge(bd, on="CNV_ID", how="outer")
rows = []
for cid in P.CNV_ID:
    g = genes[genes.CNV_ID == cid].Gene.unique(); s = np.array([s_of(x) for x in g])
    d = dict(CNV_ID=cid, n_genes=len(g), n_with_shet=int(np.sum(~np.isnan(s))), shet_sum=np.nansum(s), n_shet_ge_0p03=int(np.nansum(s >= 0.03)), n_shet_ge_0p18=int(np.nansum(s >= 0.18)))
    for tr in ["height", "bmi"]:
        zz = np.array([tw[tr].get(x, np.nan) for x in g]); con = s >= 0.03
        d[f"mean_z_constrained_{tr}"] = np.nanmean(zz[con]) if np.any(con & ~np.isnan(zz)) else np.nan
        d[f"mean_z_unconstrained_{tr}"] = np.nanmean(zz[~con]) if np.any(~con & ~np.isnan(zz)) else np.nan
        d[f"n_twas_constrained_{tr}"] = int(np.sum(con & ~np.isnan(zz)))
    rows.append(d)
P = P.merge(pd.DataFrame(rows), on="CNV_ID"); P.to_csv(HERE + r"\twas_constraint_per_cnv.csv", index=False)
print(P[["CNV_ID", "n_genes", "n_shet_ge_0p03", "shet_sum", "dos_height", "dos_bmi"]].round(3).to_string(index=False))
# ---- 3. models
def fit(df, y, xs):
    d = df.dropna(subset=[y] + xs); Z = lambda v: (v - v.mean()) / v.std()
    X = np.column_stack([np.ones(len(d))] + [Z(d[x]) for x in xs]); yy = Z(d[y]); b = np.linalg.lstsq(X, yy, rcond=None)[0]
    e = yy - X @ b; n, p = X.shape; s2 = e @ e / (n - p); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X))); pv = 2 * stats.t.sf(np.abs(b / se), n - p)
    R2 = 1 - e @ e / ((yy - yy.mean()) @ (yy - yy.mean()))
    return n, R2, {x: (b[i + 1], pv[i + 1]) for i, x in enumerate(xs)}
M = []
for tr in ["height", "bmi"]:
    for wl, D in [("all", P), ("without Williams", P[P.CNV_ID != "7q11.23_WBS"])]:
        for burden in ["shet_sum", "n_shet_ge_0p03", "n_genes"]:
            for xs in ([burden], [f"sumz_{tr}"], [burden, f"sumz_{tr}"]):
                for kind in ["least squares", "rank"]:
                    DD = D.copy()
                    if kind == "rank":
                        for c in [f"dos_{tr}", burden, f"sumz_{tr}"]: DD[c] = DD[c].rank()
                    n, R2, co = fit(DD, f"dos_{tr}", xs)
                    rec = dict(trait=tr, sample=wl, burden=burden, model=" + ".join(xs), kind=kind, n=n, R2=R2)
                    for x, (b, p) in co.items():
                        k = "burden" if x == burden else "twas"; rec[f"beta_{k}"] = b; rec[f"p_{k}"] = p
                    M.append(rec)
M = pd.DataFrame(M).drop_duplicates(subset=["trait", "sample", "burden", "model", "kind"])
M.to_csv(HERE + r"\twas_constraint_models.csv", index=False)
print(M[(M.kind == "least squares")].round(3).to_string(index=False))
