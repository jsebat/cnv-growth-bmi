# Asymmetric CNVs after removing their duplication-axis offset (Jon 2026-11-02): asymmetric CNVs sit below the mirror line
# (Fig. 1e, 2a). If the single hit adds a roughly constant shift on the duplication axis, then removing that shift should leave
# their duplication effects tracking their genes' TWAS-predicted direction like mirror CNVs.
# Offset: (a) mirror slope fixed (DUP = -rho * DEL + c, rho = 0.351 from the Fig. 3 mirror fit, height_mirror_parameters.csv),
#         c = inverse-variance-weighted mean of DUP + rho*DEL over asymmetric CNVs;
#         (b) free line through the asymmetric CNVs (weighted least squares of DUP on DEL), c = its intercept.
# Adjusted DUP = DUP - c for asymmetric CNVs; mirror CNVs unchanged. Note: a constant offset does not change ranks within the
# asymmetric group, so the within-asymmetric Spearman correlation is unchanged; it changes (i) the pooled all-CNV correlation and
# (ii) per-allele sign concordance. Summed z = contiguous-span TWAS z (as in the circular permutation).
# Pooled correlation tested by the same circular genomic permutation (10,000). Output: twas_asymmetric_offset.csv
import numpy as np, pandas as pd
from scipy import stats

HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
RHO = 0.351140004155472
S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"].set_index("CNV_ID")
by = pd.read_csv(HERE + r"\concordance_twas_by_cnv.csv"); by = by[by.trait == "height"].set_index("CNV_ID")
S["DEL_se"] = np.nan; S["DUP_se"] = np.nan
t3 = pd.read_excel(ROOT + r"\NatureSubmission_04_13_2026\Sacks Supplementary Tables.xlsx", sheet_name="T3 CNV main effects and frequen")
import re; norm = lambda s: re.sub(r"[^0-9a-z]", "", str(s).lower())
t3["allele"] = t3.locus.str.extract(r"_(DEL|DUP)$")[0]; t3["key"] = t3.locus.str.replace(r"_(DEL|DUP)$", "", regex=True).map(norm)
for a in ["DEL", "DUP"]:
    x = t3[t3.allele == a].set_index("key"); S[f"{a}_se"] = [x.seHeight.get(norm(c), np.nan) for c in S.index]
A = S[S.cls == "asymmetric"]
w = 1 / (A.DUP_se ** 2 + (RHO * A.DEL_se) ** 2)
c_fixed = float(np.sum(w * (A.DUP_effect + RHO * A.DEL_effect)) / np.sum(w))
X = np.column_stack([np.ones(len(A)), A.DEL_effect]); W = np.diag(1 / A.DUP_se ** 2)
b = np.linalg.solve(X.T @ W @ X, X.T @ W @ A.DUP_effect.values); c_free, slope_free = float(b[0]), float(b[1])
print(f"offset, mirror slope fixed: c = {c_fixed:.3f} SD; free asymmetric line: intercept {c_free:.3f}, slope {slope_free:.3f}")

rows = []
for lab, c in [("mirror slope fixed", c_fixed), ("free asymmetric line", c_free)]:
    d = S.copy(); d["DUP_adj"] = np.where(d.cls == "asymmetric", d.DUP_effect - c, d.DUP_effect)
    asy = d[d.cls == "asymmetric"]
    conc_raw = int((np.sign(asy.summed_z) == np.sign(asy.DUP_effect)).sum()); conc_adj = int((np.sign(asy.summed_z) == np.sign(asy.DUP_adj)).sum())
    r_all_raw = stats.spearmanr(d.summed_z, d.DUP_effect)[0]; r_all_adj = stats.spearmanr(d.summed_z, d.DUP_adj)[0]
    rp_all_raw = stats.pearsonr(d.summed_z, d.DUP_effect)[0]; rp_all_adj = stats.pearsonr(d.summed_z, d.DUP_adj)[0]
    r_asy = stats.spearmanr(asy.summed_z, asy.DUP_adj)[0]
    rows.append(dict(offset_method=lab, offset_sd=c, n_asymmetric=len(asy), asym_DUP_sign_concordant_raw=conc_raw, asym_DUP_sign_concordant_adjusted=conc_adj,
                     asym_DUP_sign_p_adjusted=stats.binomtest(conc_adj, len(asy), 0.5, alternative="greater").pvalue,
                     spearman_within_asymmetric=r_asy, spearman_all_raw=r_all_raw, spearman_all_adjusted=r_all_adj,
                     pearson_all_raw=rp_all_raw, pearson_all_adjusted=rp_all_adj))
    print(f"\n[{lab}] c = {c:.3f}"); print(asy[["DEL_effect", "DUP_effect", "DUP_adj", "summed_z"]].round(3).to_string())

# circular permutation of the pooled correlation with offset-adjusted duplications (same machinery as twas_circular_permutation.py)
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
tw = pd.read_csv(HERE + r"\twas30\height.txt", sep=r"\s+")
g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
spans = {}
for cid, sub in genes.groupby("CNV_ID"):
    ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
    if ix and cid in S.index: spans[cid] = (ix[0], ix[-1] - ix[0] + 1)
cids = list(spans); order = sorted(cids, key=lambda c: spans[c][0]); clusters, cur, end = [], [], -1
for c in order:
    s0, L = spans[c]
    if cur and s0 <= end: cur.append(c); end = max(end, s0 + L - 1)
    else:
        if cur: clusters.append(cur)
        cur, end = [c], s0 + L - 1
clusters.append(cur)
anchor = {c: min(spans[x][0] for x in cl) for cl in clusters for c in cl}; cl_id = {c: i for i, cl in enumerate(clusters) for c in cl}
rel = np.array([spans[c][0] - anchor[c] for c in cids]); Ls = np.array([spans[c][1] for c in cids]); ci = np.array([cl_id[c] for c in cids]); anc = np.array([anchor[c] for c in cids])
S_obs = np.array([cs[spans[c][0] + spans[c][1]] - cs[spans[c][0]] for c in cids])
rng = np.random.default_rng(11); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; S_perm = cs[st + Ls[None, :]] - cs[st]
Sr = np.apply_along_axis(stats.rankdata, 1, S_perm); Sr = (Sr - Sr.mean(1, keepdims=True)) / Sr.std(1, keepdims=True)
for row in rows:
    c = row["offset_sd"]; y = np.array([S.loc[x, "DUP_effect"] - (c if S.loc[x, "cls"] == "asymmetric" else 0) for x in cids])
    yr = stats.rankdata(y); yr = (yr - yr.mean()) / yr.std(); null = (Sr * yr[None, :]).mean(1)
    r_obs = stats.spearmanr(S_obs, y)[0]; row["spearman_all_adjusted_span"] = r_obs; row["p_circular_all_adjusted"] = (1 + np.sum(null >= r_obs)) / 10001
    y0 = np.array([S.loc[x, "DUP_effect"] for x in cids]); y0r = stats.rankdata(y0); y0r = (y0r - y0r.mean()) / y0r.std()
    row["p_circular_all_raw"] = (1 + np.sum((Sr * y0r[None, :]).mean(1) >= stats.spearmanr(S_obs, y0)[0])) / 10001
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_asymmetric_offset.csv", index=False)
pd.set_option("display.width", 250); print(); print(R.round(4).T.to_string())
