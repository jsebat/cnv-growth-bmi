# Joint model of the CNV dosage effect on constrained-gene burden and summed TWAS z, with circular-permutation P for the TWAS
# term (Jon 2026-11-03). Same CNVs, spans and permutations as Fig. 2E (twas_dosage_panel_stats.py / twas_dosage_bmi.py).
# Statistic: partial Spearman correlation of dosage effect with summed TWAS z given burden (ranks residualised on the burden rank);
# the null keeps each CNV's burden and dosage effect and moves the TWAS z (locus clusters rotated as units, 10,000 times).
# Burden = sum of GeneBayes s_het over the CNV's genes (twas_constraint_per_cnv.csv). Also: least-squares joint fit (standardised),
# Williams residual after burden, and the burden-only correlation. Output: twas_constraint_perm.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
bur = pd.read_csv(HERE + r"\twas_constraint_per_cnv.csv").set_index("CNV_ID")
rk = lambda v: stats.rankdata(v)
def resid(v, b):   # residual of v on b (with intercept)
    X = np.column_stack([np.ones(len(b)), b]); return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
out = []
for trait in ["height", "bmi"]:
    S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == trait].set_index("CNV_ID")
    chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
    tw = pd.read_csv(HERE + rf"\twas30\{trait}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
    N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
    spans = {}
    for cid, sub in genes.groupby("CNV_ID"):
        if cid in DROP or cid not in bur.index: continue
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
    y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids]); b = bur.loc[cids, "shet_sum"].values
    rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
    st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Xp = cs[st + Ls[None, :]] - cs[st]
    for lab, k in [("all", np.ones(len(cids), bool)), ("without Williams", np.array(cids) != "7q11.23_WBS")]:
        yr, br, xr = rk(y[k]), rk(b[k]), rk(x[k])
        ry = resid(yr, br); rx = resid(xr, br); obs = np.corrcoef(ry, rx)[0, 1]
        Xr = np.apply_along_axis(rk, 1, Xp[:, k]); null = np.array([np.corrcoef(ry, resid(v, br))[0, 1] for v in Xr])
        rho_tw = stats.spearmanr(x[k], y[k]); rho_b = stats.spearmanr(b[k], y[k])
        Z = lambda v: (v - v.mean()) / v.std(); XX = np.column_stack([np.ones(k.sum()), Z(b[k]), Z(x[k])]); yy = Z(y[k])
        co = np.linalg.lstsq(XX, yy, rcond=None)[0]; e = yy - XX @ co; n, p = XX.shape; s2 = e @ e / (n - p); se = np.sqrt(np.diag(s2 * np.linalg.inv(XX.T @ XX)))
        pv = 2 * stats.t.sf(np.abs(co / se), n - p); R2 = 1 - e @ e / (yy @ yy)
        XT = np.column_stack([np.ones(k.sum()), Z(x[k])]); eT = yy - XT @ np.linalg.lstsq(XT, yy, rcond=None)[0]; R2T = 1 - eT @ eT / (yy @ yy)
        rec = dict(trait=trait, sample=lab, n_cnv=int(k.sum()), spearman_twas=rho_tw[0], spearman_burden=rho_b[0], p_spearman_burden=rho_b[1],
                   partial_spearman_twas_given_burden=obs, p_circular_one_sided=(1 + np.sum(null >= obs)) / 10001,
                   ls_beta_burden=co[1], ls_p_burden=pv[1], ls_beta_twas=co[2], ls_p_twas=pv[2], R2_twas_only=R2T, R2_burden_plus_twas=R2)
        if lab == "all" and "7q11.23_WBS" in cids:
            j = cids.index("7q11.23_WBS"); full = yy - XX @ co; rec["williams_std_residual_joint"] = full[j] / np.sqrt(s2)
            rec["williams_std_residual_twas_only"] = eT[j] / np.sqrt(eT @ eT / (n - 2))
        out.append(rec)
R = pd.DataFrame(out); R.to_csv(HERE + r"\twas_constraint_perm.csv", index=False)
pd.set_option("display.width", 250); print(R.round(4).T.to_string())
