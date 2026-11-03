# Can a ~2-fold higher gene count by itself push the summed TWAS z of asymmetric CNVs down? (Jon 2026-11-02)
# 1. Expectation from the genome-wide mean z per gene (TWAS30 height): expected sum = n_genes x mean z.
# 2. Circular permutation null (each CNV moved with its own gene count): null distribution of the class difference in summed z
#    (asymmetric mean - mirror mean); its mean is the shift attributable to gene count alone, its SD the spread.
# 3. Per asymmetric CNV: observed summed z vs random genomic windows with the same number of genes (z score and percentile).
# Output: twas_genecount_null.csv, twas_genecount_null_per_cnv.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]; ns = {}; exec(src, ns)
x_obs, cls, cids, Xp, Ls, g = ns["x_obs"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["Ls"], ns["g"]
asym = cls == "asymmetric"; mz = g.z.mean()
d_obs = x_obs[asym].mean() - x_obs[~asym].mean(); d_null = Xp[:, asym].mean(1) - Xp[:, ~asym].mean(1)
R = pd.DataFrame([dict(genome_mean_z_per_gene=mz, median_genes_mirror=np.median(Ls[~asym]), median_genes_asym=np.median(Ls[asym]),
    mean_genes_mirror=Ls[~asym].mean(), mean_genes_asym=Ls[asym].mean(),
    expected_difference_from_gene_count=mz * (Ls[asym].mean() - Ls[~asym].mean()),
    null_mean_difference_circular=d_null.mean(), null_sd_difference=d_null.std(), observed_difference=d_obs,
    observed_z_vs_null=(d_obs - d_null.mean()) / d_null.std(), p_one_sided=(1 + np.sum(d_null <= d_obs)) / (len(d_null) + 1))])
R.to_csv(HERE + r"\twas_genecount_null.csv", index=False); print(R.round(4).T.to_string())
per = pd.DataFrame(dict(CNV_ID=cids, cls=cls, n_genes=Ls, summed_z=x_obs, null_mean=Xp.mean(0), null_sd=Xp.std(0),
                        z_vs_same_size_windows=(x_obs - Xp.mean(0)) / Xp.std(0), percentile=(Xp < x_obs[None, :]).mean(0)))
per.to_csv(HERE + r"\twas_genecount_null_per_cnv.csv", index=False)
print(per[per.cls == "asymmetric"].round(3).sort_values("percentile").to_string(index=False))
