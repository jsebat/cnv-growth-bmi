# Is higher predicted expression generally associated with greater height (or BMI)? JS 2026-11-02 (Jon's question after the
# CNV-TWAS concordance result). Gene-level TWAS z (mean over expression panels, TWAS30 / Mancuso et al. 2017):
#   1. genome-wide mean z and fraction of genes with z > 0, all genes and transcriptome-wide significant genes;
#      SE from a block jackknife over 10-Mb genomic blocks (genes are correlated through LD and co-regulation)
#   2. per CNV: does the summed z just track gene count? Spearman of S with gene count, and partial correlation of the CNV
#      effects with S adjusting for gene count; also with the mean z per gene (S / genes with TWAS) instead of the sum
# Outputs: twas_genomewide_direction.csv, concordance_twas_genecount_adjusted.csv
import numpy as np, pandas as pd
from scipy import stats

HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
out = []
for tr in ["height", "bmi"]:
    tw = pd.read_csv(HERE + rf"\twas30\{tr}.txt", sep=r"\s+")
    thr = 0.05 / len(tw)
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), p_best=("TWAS.P", "min"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["block"] = g["chr"].astype(str) + ":" + (g["pos"] // 10_000_000).astype(str)
    for subset, d in [("all genes", g), ("TWAS-significant genes", g[g.p_best < thr])]:
        blocks = d.block.unique(); B = len(blocks)
        est = lambda x: (x.z.mean(), (x.z > 0).mean())
        m0, f0 = est(d)
        jk = np.array([est(d[d.block != b]) for b in blocks])
        se = np.sqrt((B - 1) / B * ((jk - jk.mean(0)) ** 2).sum(0))
        zt = (m0 / se[0], (f0 - 0.5) / se[1])
        out.append(dict(trait=tr, genes=subset, n_genes=len(d), n_blocks=B, mean_z=m0, se_mean_z=se[0], p_mean=2 * stats.norm.sf(abs(zt[0])),
                        frac_positive=f0, se_frac=se[1], p_frac=2 * stats.norm.sf(abs(zt[1]))))
G = pd.DataFrame(out); G.to_csv(HERE + r"\twas_genomewide_direction.csv", index=False)
pd.set_option("display.width", 220); print(G.round(4).to_string(index=False))

# ---- gene count check on the CNV concordance
by = pd.read_csv(HERE + r"\concordance_twas_by_cnv.csv")
def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r = np.corrcoef(ex, ey)[0, 1]; n = len(x); t = r * np.sqrt((n - 3) / (1 - r ** 2)); return r, 2 * stats.t.sf(abs(t), n - 3)
rows = []
for tr in ["height", "bmi"]:
    d0 = by[by.trait == tr].dropna(subset=["S"]).copy(); d0["zbar"] = d0.S / d0.n_twas
    for gname in (["all", "mirror", "asymmetric"] if tr == "height" else ["all"]):   # no asymmetric CNVs for BMI
        d = d0 if gname == "all" else d0[d0.cls == gname]
        r_sn = stats.spearmanr(d.S, d.n_genes)[0]
        for a in ["DEL", "DUP"]:
            y = d[f"{a}_{tr}"]; ok = y.notna()
            r_eff_n = stats.spearmanr(d.n_genes[ok], y[ok])[0]
            pr, pp = partial_spearman(d.S[ok].values, y[ok].values, d.n_genes[ok].values)
            rm, pm = stats.spearmanr(d.zbar[ok], y[ok])
            rows.append(dict(trait=tr, group=gname, allele=a, n_cnv=int(ok.sum()), rho_S_vs_genecount=r_sn, rho_effect_vs_genecount=r_eff_n,
                             partial_rho_S_given_genecount=pr, p_partial=pp, rho_mean_z_per_gene=rm, p_mean_z=pm))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\concordance_twas_genecount_adjusted.csv", index=False)
print(); print(R.round(3).to_string(index=False))
