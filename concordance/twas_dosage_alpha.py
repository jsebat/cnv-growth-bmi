# Fig. 2E test with the TWAS effect size ALPHA (Mancuso et al. 2017; mean over panels) summed over the CNV genes, in place of the
# TWAS z (Jon 2026-11-03: z depends on GWAS N and heritability). Otherwise identical to twas_dosage_bmi.py.
# BMI version of Fig. 2E (Jon 2026-11-03): per CNV, BMI dosage effect (DUP - DEL)/2 against the summed BMI TWAS z of its genes
# (Mancuso et al. 2017 BMI TWAS; contiguous span, as in twas_dosage_panel_stats.py), on the same non-overlapping CNV set as Fig. 2E
# (spanning intervals kept, 22q11.2 A-D). Spearman with circular genomic permutation P (10,000; locus clusters moved as units),
# all CNVs and without Williams. BMI has no asymmetric class (all CNVs mirror for BMI), so only the pooled tests are primary.
# Output: twas_dosage_alpha_stats.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
out = []
for STAT, trait in [(st, tr) for st in ["ALPHA", "TWAS.Z"] for tr in ["height", "bmi"]]:   # height rerun as a check against Fig. 2E
    S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == trait].set_index("CNV_ID")
    chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
    tw = pd.read_csv(HERE + rf"\twas30\{trait}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=(STAT, "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
    N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
    spans = {}
    for cid, sub in genes.groupby("CNV_ID"):
        if cid in DROP: continue
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
    y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids])
    rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
    st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Xp = cs[st + Ls[None, :]] - cs[st]
    def rho_p(k):
        r = stats.spearmanr(x_obs[k], y[k])[0]
        yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
        Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
        null = (Xr * yr[None, :]).mean(1)
        return r, (1 + np.sum(null >= r)) / 10001, (1 + np.sum(np.abs(null) >= abs(r))) / 10001
    for lab, k in [("all", np.ones(len(cids), bool)), ("all without Williams", np.array(cids) != "7q11.23_WBS")]:
        r, p1, p2 = rho_p(k); out.append(dict(statistic=STAT, trait=trait, group=lab, n_cnv=int(k.sum()), spearman_rho=r, p_circular_one_sided=p1, p_circular_two_sided=p2))
R = pd.DataFrame(out); R.to_csv(HERE + r"\twas_dosage_alpha_stats.csv", index=False); print(R.round(4).to_string(index=False))
