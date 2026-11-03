# Do asymmetric CNVs carry genes whose higher expression lowers height? Mean summed TWAS z (contiguous span) in asymmetric vs mirror
# CNVs; null from the circular genomic permutation (locus clusters moved as units, classes fixed), 10,000 permutations; also per
# gene (mean z per gene) to remove the effect of interval size. JS 2026-11-02. Output: twas_class_mean_test.csv
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
asym = np.array([S.loc[c, "cls"] == "asymmetric" for c in cids])
S_obs = np.array([cs[spans[c][0] + Ls[i]] - cs[spans[c][0]] for i, c in enumerate(cids)])
rng = np.random.default_rng(5); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Sp = cs[st + Ls[None, :]] - cs[st]
rows = []
for lab, obs, perm in [("summed z", S_obs, Sp), ("mean z per gene", S_obs / Ls, Sp / Ls[None, :])]:
    d_obs = obs[asym].mean() - obs[~asym].mean(); d_null = perm[:, asym].mean(1) - perm[:, ~asym].mean(1)
    rows.append(dict(statistic=lab, mean_asymmetric=obs[asym].mean(), mean_mirror=obs[~asym].mean(), difference=d_obs,
                     p_circular_one_sided_lower=(1 + np.sum(d_null <= d_obs)) / 10001, p_mannwhitney_two_sided=stats.mannwhitneyu(obs[asym], obs[~asym]).pvalue,
                     n_asym_negative=int((obs[asym] < 0).sum()), n_asym=int(asym.sum()), n_mirror_negative=int((obs[~asym] < 0).sum()), n_mirror=int((~asym).sum())))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_class_mean_test.csv", index=False); print(R.round(4).T.to_string())
