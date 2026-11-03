# Fig. 2E test with brain TWAS (Jon 2026-11-03): FUSION on PsychENCODE cortex weights (PECfull; ConGA-GWAS install), GWAS from
# Tan et al. 2024 (population estimates; BMI and height from the same study, full imputed coverage). Run on Expanse
# (psychencode_twas/twas_cnv.sbatch); per-chromosome .dat files copied to psychencode_twas/raw/<TRAIT>_TAN/.
# Per CNV: dosage effect (DUP - DEL)/2 against summed brain TWAS z over the contiguous run of genes spanning the CNV (genes ordered
# by position; Ensembl IDs mapped to the CNV gene symbols through HGNC). Same non-overlapping CNVs and circular permutation (10,000;
# locus clusters as units) as Fig. 2E; plain Spearman and partial Spearman given constrained-gene burden (summed s_het).
# Also the genome-wide signal (mean |z|, genes with |z| > 4.5) for comparison with the Mancuso et al. TWAS.
# Output: twas_psychencode_stats.csv, twas_psychencode_per_cnv.csv
import numpy as np, pandas as pd, glob, sys
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
bur = pd.read_csv(HERE + r"\twas_constraint_per_cnv.csv").set_index("CNV_ID")
hg = pd.read_csv(HERE + r"\hgnc_map.tsv", sep="\t")
ens2sym = dict(zip(hg.ensembl_gene_id, hg.symbol))
alias = {}
for _, r in hg.iterrows():
    for col in ("prev_symbol", "alias_symbol"):
        if isinstance(r[col], str):
            for a in r[col].split("|"): alias.setdefault(a, r.symbol)
rk = lambda v: stats.rankdata(v)
def resid(v, b):
    X = np.column_stack([np.ones(len(b)), b]); return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
out, per = [], []
for trait, cnvtrait in [("BMI_TAN", "bmi"), ("HEIGHT_TAN", "height")]:
    fs = sorted(glob.glob(HERE + rf"\psychencode_twas\raw\{trait}\PECfull.*.dat"))
    if len(fs) < 22: print(trait, "incomplete:", len(fs), "chromosomes"); continue
    t = pd.concat([pd.read_csv(f, sep="\t") for f in fs]); t = t.dropna(subset=["TWAS.Z"])
    t["ens"] = t.ID.str.split(".").str[0]; t["sym"] = t.ens.map(ens2sym)
    g = t.groupby("ens").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("P0", "first"), sym=("sym", "first")).reset_index()
    g = g.sort_values(["chr", "pos"]).reset_index(drop=True)
    out.append(dict(trait=trait, group="genome-wide", n_cnv=len(g), spearman_rho=g.z.abs().mean(), partial_rho_given_burden=int((g.z.abs() > 4.5).sum()),
                    p_circular_one_sided=g.z.mean(), note="genome-wide: n_cnv = genes with z; spearman_rho = mean |z|; partial = genes |z|>4.5; p = mean z"))
    N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))])
    pos_of = {}
    for i, s in enumerate(g.sym):
        if isinstance(s, str): pos_of.setdefault(s, i)
    S = pd.read_csv(HERE + r"\twas_per_cnv_sign_test.csv"); S = S[S.trait == cnvtrait].set_index("CNV_ID")
    spans = {}
    for cid, sub in genes.groupby("CNV_ID"):
        if cid in DROP or cid not in bur.index: continue
        ix = sorted(pos_of[x] if x in pos_of else pos_of[alias[x]] for x in sub.Gene if x in pos_of or alias.get(x) in pos_of)
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
    for c, xi, yi, Li in zip(cids, x, y, Ls): per.append(dict(trait=trait, CNV_ID=c, summed_brain_twas_z=xi, dosage_effect=yi, n_genes_span=Li, shet_sum=bur.loc[c, "shet_sum"]))
    for lab, k in [("all", np.ones(len(cids), bool)), ("without Williams", np.array(cids) != "7q11.23_WBS")]:
        r = stats.spearmanr(x[k], y[k])[0]
        yr = rk(y[k]); yz = (yr - yr.mean()) / yr.std(); Xr = np.apply_along_axis(rk, 1, Xp[:, k]); Xz = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
        p = (1 + np.sum((Xz * yz[None, :]).mean(1) >= r)) / 10001
        br = rk(b[k]); ry = resid(yr, br); obs = np.corrcoef(ry, resid(rk(x[k]), br))[0, 1]
        null = np.array([np.corrcoef(ry, resid(v, br))[0, 1] for v in Xr]); pp = (1 + np.sum(null >= obs)) / 10001
        out.append(dict(trait=trait, group=lab, n_cnv=int(k.sum()), spearman_rho=r, p_circular_one_sided=p, partial_rho_given_burden=obs, p_partial_circular=pp))
R = pd.DataFrame(out); R.to_csv(HERE + r"\twas_psychencode_stats.csv", index=False); pd.DataFrame(per).to_csv(HERE + r"\twas_psychencode_per_cnv.csv", index=False)
pd.set_option("display.width", 250); print(R.round(4).to_string(index=False))
