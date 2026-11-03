# Is the Williams-Beuren rescue due to constrained genes or to gene count? (Jon 2026-11-03). Height, all 23 CNVs of Fig. 2E.
# For each covariate (summed s_het, genes with s_het >= 0.03, all genes, unconstrained genes): least-squares joint model of the
# dosage effect on covariate + summed TWAS z (standardised), the Williams standardised residual, and the partial Spearman of
# dosage with summed TWAS z given the covariate with circular-permutation P (as twas_constraint_perm.py).
# Output: twas_williams_covariate.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
P = pd.read_csv(HERE + r"\twas_constraint_per_cnv.csv").set_index("CNV_ID"); P["n_unconstrained"] = P.n_with_shet - P.n_shet_ge_0p03
S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"].set_index("CNV_ID")
tw = pd.read_csv(HERE + r"\twas30\height.txt", sep=r"\s+")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
spans = {}
for cid, sub in genes.groupby("CNV_ID"):
    if cid in DROP or cid not in P.index: continue
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
x = np.array([cs[spans[c][0] + Ls[i]] - cs[spans[c][0]] for i, c in enumerate(cids)])
y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids])
rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Xp = cs[st + Ls[None, :]] - cs[st]
rk = lambda v: stats.rankdata(v)
def resid(v, b):
    X = np.column_stack([np.ones(len(b)), b]); return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
Z = lambda v: (v - v.mean()) / v.std()
j = cids.index("7q11.23_WBS"); yy = Z(y); out = []
XT = np.column_stack([np.ones(len(y)), Z(x)]); eT = yy - XT @ np.linalg.lstsq(XT, yy, rcond=None)[0]
out.append(dict(covariate="none (TWAS only)", williams_std_residual=eT[j] / np.sqrt(eT @ eT / (len(y) - 2)), R2=1 - eT @ eT / (yy @ yy),
                partial_spearman_twas=stats.spearmanr(x, y)[0]))
Xr = np.apply_along_axis(rk, 1, Xp)
for cov in ["shet_sum", "n_shet_ge_0p03", "n_genes", "n_unconstrained"]:
    b = P.loc[cids, cov].values.astype(float)
    XX = np.column_stack([np.ones(len(y)), Z(b), Z(x)]); co = np.linalg.lstsq(XX, yy, rcond=None)[0]; e = yy - XX @ co
    n, p = XX.shape; s2 = e @ e / (n - p); se = np.sqrt(np.diag(s2 * np.linalg.inv(XX.T @ XX))); pv = 2 * stats.t.sf(np.abs(co / se), n - p)
    br = rk(b); ry = resid(rk(y), br); obs = np.corrcoef(ry, resid(rk(x), br))[0, 1]
    null = np.array([np.corrcoef(ry, resid(v, br))[0, 1] for v in Xr])
    out.append(dict(covariate=cov, williams_value=P.loc["7q11.23_WBS", cov], williams_rank_of_23=int(rk(b)[j]), beta_cov=co[1], p_cov=pv[1],
                    beta_twas=co[2], p_twas=pv[2], R2=1 - e @ e / (yy @ yy), williams_std_residual=e[j] / np.sqrt(s2),
                    partial_spearman_twas=obs, p_circular=(1 + np.sum(null >= obs)) / 10001))
R = pd.DataFrame(out); R.to_csv(HERE + r"\twas_williams_covariate.csv", index=False)
pd.set_option("display.width", 250); print(R.round(3).to_string(index=False))
