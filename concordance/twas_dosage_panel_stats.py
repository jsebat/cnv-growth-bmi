# Statistics for the single-panel figure (Jon 2026-11-02): per CNV, gene-dosage effect on height = (DUP - DEL)/2 (one value per
# reciprocal CNV, both alleles on one scale: positive = more copies, taller) against the summed TWAS z of its genes (contiguous
# span, as in the circular permutation). Spearman per class (all, mirror, asymmetric) with circular genomic permutation P
# (10,000; locus clusters moved as units). Output: twas_dosage_panel.csv (points), twas_dosage_panel_stats.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"].set_index("CNV_ID")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
tw = pd.read_csv(HERE + r"\twas30\height.txt", sep=r"\s+")
g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
spans = {}
for cid, sub in genes.groupby("CNV_ID"):
    ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
    if ix and cid in S.index and not np.isnan(S.loc[cid, "DEL_effect"]) and not np.isnan(S.loc[cid, "DUP_effect"]): spans[cid] = (ix[0], ix[-1] - ix[0] + 1)
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
x_obs = np.array([cs[spans[c][0] + Ls[i]] - cs[spans[c][0]] for i, c in enumerate(cids)])
y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids]); cls = np.array([S.loc[c, "cls"] for c in cids])
rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Xp = cs[st + Ls[None, :]] - cs[st]
out = []
for grp in ["all", "mirror", "asymmetric"]:
    k = np.ones(len(cids), bool) if grp == "all" else cls == grp
    r = stats.spearmanr(x_obs[k], y[k])[0]
    yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
    null = (Xr * yr[None, :]).mean(1)
    out.append(dict(group=grp, n_cnv=int(k.sum()), spearman_rho=r, p_circular_one_sided=(1 + np.sum(null >= r)) / 10001))
pd.DataFrame(dict(CNV_ID=cids, cls=cls, summed_twas_z=x_obs, DEL=[S.loc[c, "DEL_effect"] for c in cids], DUP=[S.loc[c, "DUP_effect"] for c in cids],
                  dosage_effect=y, n_genes_twas=Ls)).to_csv(HERE + r"\twas_dosage_panel.csv", index=False)
R = pd.DataFrame(out)
# Theil-Sen lines (median pairwise slope; robust to single points, consistent with the rank correlation) for the plot
for grp in ["mirror", "asymmetric"]:
    k = cls == grp; ts = stats.theilslopes(y[k], x_obs[k])
    R.loc[R.group == grp, ["theilsen_slope", "theilsen_intercept", "theilsen_slope_lo", "theilsen_slope_hi"]] = [ts[0], ts[1], ts[2], ts[3]]
# leave-one-out (Jon 2026-11-02): drop each CNV in turn, recompute rho and circular P within the class
loo = []
def rho_p(k):
    r = stats.spearmanr(x_obs[k], y[k])[0]
    yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
    return r, (1 + np.sum((Xr * yr[None, :]).mean(1) >= r)) / 10001
for grp in ["mirror", "asymmetric"]:
    idx = np.where(cls == grp)[0]
    for j in idx:
        k = (cls == grp); k[j] = False
        r, p = rho_p(k); loo.append(dict(group=grp, dropped=cids[j], spearman_rho=r, p_circular_one_sided=p))
L = pd.DataFrame(loo); L.to_csv(HERE + r"\twas_dosage_panel_loo.csv", index=False)
for grp in ["mirror", "asymmetric"]:
    l = L[L.group == grp]
    R.loc[R.group == grp, ["loo_rho_min", "loo_rho_max", "loo_p_max", "loo_p_min"]] = [l.spearman_rho.min(), l.spearman_rho.max(), l.p_circular_one_sided.max(), l.p_circular_one_sided.min()]
w = L[(L.group == "asymmetric") & (L.dropped == "7q11.23_WBS")].iloc[0]
R.loc[R.group == "asymmetric", ["rho_without_williams", "p_without_williams"]] = [w.spearman_rho, w.p_circular_one_sided]
k = (cls == "asymmetric") & (np.array(cids) != "7q11.23_WBS"); ts = stats.theilslopes(y[k], x_obs[k])
R.loc[R.group == "asymmetric", ["theilsen_slope_no_williams", "theilsen_intercept_no_williams"]] = [ts[0], ts[1]]
# pooled line through all CNVs except Williams (Jon 2026-11-02: without Williams mirror and asymmetric CNVs share one line;
# class adds nothing, twas_same_line_test.py)
k = np.array(cids) != "7q11.23_WBS"; r, p = rho_p(k); ts = stats.theilslopes(y[k], x_obs[k])
R = pd.concat([R, pd.DataFrame([dict(group="all without Williams", n_cnv=int(k.sum()), spearman_rho=r, p_circular_one_sided=p,
                                     theilsen_slope=ts[0], theilsen_intercept=ts[1], theilsen_slope_lo=ts[2], theilsen_slope_hi=ts[3])])], ignore_index=True)
R.to_csv(HERE + r"\twas_dosage_panel_stats.csv", index=False); pd.set_option("display.width", 250); print(R.round(4).T.to_string())
print(L.sort_values(["group", "spearman_rho"]).round(3).to_string(index=False))
