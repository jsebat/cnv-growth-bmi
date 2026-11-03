# Circular genomic permutation null for the CNV-TWAS concordance (Jon 2026-11-02): how often does a random placement of the
# CNVs in the genome, with the same number of genes, give a correlation as strong as observed?
# Genes with TWAS data are ordered along the genome (chromosome, start; TWAS30 coordinates) and treated as a circle. Each CNV is
# represented by the contiguous run of ordered genes spanning its TWAS genes. Overlapping/nested intervals (e.g. 22q11.21 A-D and
# its subregions) form one locus cluster that is shifted as a unit, preserving gene counts, nesting and the local correlation of
# neighbouring genes' TWAS z. Each permutation draws an independent random offset per cluster, recomputes the summed z of every
# CNV, and correlates it with the observed CNV effects (fixed). One-sided P in the expected direction (DEL -, DUP +).
# Output: concordance_twas_circular_permutation.csv
import numpy as np, pandas as pd
from scipy import stats

HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
by = pd.read_csv(HERE + r"\concordance_twas_by_cnv.csv")
NPERM = 10000; rng = np.random.default_rng(7)
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
rows = []
for tr in ["height", "bmi"]:
    tw = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
    N = len(g); zc = np.concatenate([g.z.values, g.z.values]); cs = np.concatenate([[0], np.cumsum(zc)])
    pos_of = dict(zip(g.GENE, g.index))
    eff = by[by.trait == tr].set_index("CNV_ID")
    # each CNV -> (start index, span) on the ordered gene list
    spans = {}
    for cid, sub in genes.groupby("CNV_ID"):
        ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
        if ix: spans[cid] = (ix[0], ix[-1] - ix[0] + 1)
    cids = [c for c in eff.index if c in spans]
    # locus clusters: CNVs whose spans overlap move together
    order = sorted(cids, key=lambda c: spans[c][0]); clusters, cur, end = [], [], -1
    for c in order:
        s, L = spans[c]
        if cur and s <= end: cur.append(c); end = max(end, s + L - 1)
        else:
            if cur: clusters.append(cur)
            cur, end = [c], s + L - 1
    clusters.append(cur)
    anchor = {c: min(spans[x][0] for x in cl) for cl in clusters for c in cl}
    S_of = lambda start, L: cs[start + L] - cs[start]
    S_obs = np.array([S_of(*spans[c]) for c in cids])
    # permutations: per cluster offset k; CNV start -> (anchor + k + (s - anchor)) mod N
    cl_id = {c: i for i, cl in enumerate(clusters) for c in cl}
    rel = np.array([spans[c][0] - anchor[c] for c in cids]); Ls = np.array([spans[c][1] for c in cids]); ci = np.array([cl_id[c] for c in cids])
    anc = np.array([anchor[c] for c in cids])
    K = rng.integers(0, N, size=(NPERM, len(clusters)))
    starts = (anc[None, :] + K[:, ci] + rel[None, :]) % N
    S_perm = cs[starts + Ls[None, :]] - cs[starts]
    cls = eff.loc[cids, "cls"].values
    for gname in (["all", "mirror", "asymmetric"] if tr == "height" else ["all"]):   # no asymmetric CNVs for BMI
        keep = np.ones(len(cids), bool) if gname == "all" else cls == gname
        for a, es in [("DEL", -1), ("DUP", +1)]:
            y = eff.loc[cids, f"{a}_{tr}"].values; ok = keep & ~np.isnan(y)
            r_obs = stats.spearmanr(S_obs[ok], y[ok])[0]
            yr = stats.rankdata(y[ok]); yr = (yr - yr.mean()) / yr.std()
            Sr = np.apply_along_axis(stats.rankdata, 1, S_perm[:, ok]); Sr = (Sr - Sr.mean(1, keepdims=True)) / Sr.std(1, keepdims=True)
            null = (Sr * yr[None, :]).mean(1)
            p = (1 + np.sum(es * null >= es * r_obs)) / (NPERM + 1)
            rows.append(dict(trait=tr, group=gname, allele=a, n_cnv=int(ok.sum()), n_clusters=len(clusters), spearman_rho=r_obs,
                             null_mean=null.mean(), null_sd=null.std(), p_circular_one_sided=p))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\concordance_twas_circular_permutation.csv", index=False)
pd.set_option("display.width", 200); print(R.round(4).to_string(index=False))
