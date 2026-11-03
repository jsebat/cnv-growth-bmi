# Per-CNV sign test and FDR-significant TWAS genes inside each CNV (JS 2026-11-02).
# Gene membership: TWAS genes on the contiguous run spanning each CNV's listed genes (as in the circular permutation; captures
# genes listed under older names in TWAS30).
# 1. Sign test per reciprocal CNV and trait: expected gene direction = sign(DUP effect) and -sign(DEL effect); a gene is
#    concordant with the duplication if sign(z) = sign(DUP effect), with the deletion if sign(z) = -sign(DEL effect).
#    One-sided binomial (fraction concordant > 0.5); genes within a CNV are correlated (LD, co-regulation), so descriptive.
# 2. FDR: Benjamini-Hochberg over all gene-tissue TWAS tests genome-wide per trait; a gene is significant if any panel has
#    q < 0.05; reported with its best panel, z, and concordance with the CNV's deletion and duplication effects.
# Outputs: twas_per_cnv_sign_test.csv, twas_fdr_genes_in_cnvs.csv
import numpy as np, pandas as pd
from scipy import stats

HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
by = pd.read_csv(HERE + r"\concordance_twas_by_cnv.csv")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
sign_rows, fdr_rows = [], []
for tr in ["height", "bmi"]:
    tw = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+")
    p = tw["TWAS.P"].values; o = np.argsort(p); m = len(p)
    q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; tw["q"] = np.minimum(q, 1)
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first"), q_min=("q", "min"), n_panel=("q", "size")).reset_index()
    best = tw.loc[tw.groupby("GENE")["q"].idxmin(), ["GENE", "STUDY", "TWAS.Z", "TWAS.P"]].rename(columns={"STUDY": "best_panel", "TWAS.Z": "z_best", "TWAS.P": "p_best"})
    g = g.merge(best, on="GENE"); g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
    pos_of = dict(zip(g.GENE, g.index)); eff = by[by.trait == tr].set_index("CNV_ID")
    for cid, sub in genes.groupby("CNV_ID"):
        ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
        if not ix or cid not in eff.index: continue
        d = g.iloc[ix[0]: ix[-1] + 1]; e = eff.loc[cid]
        row = dict(trait=tr, CNV_ID=cid, cls=e.cls, n_genes_twas=len(d), DEL_effect=e[f"DEL_{tr}"], DEL_p=e[f"DEL_{tr}_p"], DUP_effect=e[f"DUP_{tr}"], DUP_p=e[f"DUP_{tr}_p"],
                   n_z_positive=int((d.z > 0).sum()), summed_z=d.z.sum())
        for a, sgn in [("DEL", -1), ("DUP", +1)]:
            ef = e[f"{a}_{tr}"]
            if np.isnan(ef): continue
            k = int((np.sign(d.z) == sgn * np.sign(ef)).sum())
            row[f"{a}_n_concordant"] = k; row[f"{a}_frac_concordant"] = k / len(d)
            row[f"{a}_p_sign"] = stats.binomtest(k, len(d), 0.5, alternative="greater").pvalue
        sign_rows.append(row)
        for _, r in d[d.q_min < 0.05].iterrows():
            fdr_rows.append(dict(trait=tr, CNV_ID=cid, cls=e.cls, gene=r.GENE, best_panel=r.best_panel, z_best=r.z_best, p_best=r.p_best, q_min=r.q_min,
                                 z_mean_over_panels=r.z, n_panels=r.n_panel, DEL_effect=e[f"DEL_{tr}"], DUP_effect=e[f"DUP_{tr}"],
                                 concordant_with_DEL=np.sign(r.z_best) == -np.sign(e[f"DEL_{tr}"]), concordant_with_DUP=np.sign(r.z_best) == np.sign(e[f"DUP_{tr}"])))
S = pd.DataFrame(sign_rows); F = pd.DataFrame(fdr_rows)
S.to_csv(HERE + r"\twas_per_cnv_sign_test.csv", index=False); F.to_csv(HERE + r"\twas_fdr_genes_in_cnvs.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
cols = ["trait", "CNV_ID", "cls", "n_genes_twas", "DEL_effect", "DUP_effect", "DEL_n_concordant", "DEL_p_sign", "DUP_n_concordant", "DUP_p_sign"]
print(S[cols].round(3).to_string(index=False))
print(); print(F[["trait", "CNV_ID", "cls", "gene", "best_panel", "z_best", "q_min", "concordant_with_DEL", "concordant_with_DUP"]].round(4).to_string(index=False))
for tr in ["height", "bmi"]:
    f = F[F.trait == tr].drop_duplicates(["gene"]) if len(F) else F
    print(tr, "FDR genes in CNVs:", len(F[F.trait == tr]), "gene-CNV pairs,", F[F.trait == tr].gene.nunique() if len(F) else 0, "genes")
