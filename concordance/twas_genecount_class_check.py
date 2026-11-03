# Gene count and the class difference in summed TWAS z (Jon 2026-11-02). Asymmetric CNVs contain more genes (Fig. 3), and a sum over
# more genes has a wider spread, so the class difference and the pooled correlations are re-checked per gene and adjusted for count.
#  1. Gene count by class (genes with TWAS data on the CNV span).
#  2. Class difference in mean z per gene, and in summed z adjusted for gene count (regression S ~ class + n), circular null.
#  3. Pooled dosage correlation without Williams: summed z, mean z per gene, and partial Spearman on gene count.
# Output: twas_genecount_class_check.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]; ns = {}; exec(src, ns)
x_obs, y, cls, cids, Xp, Ls = ns["x_obs"], ns["y"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["Ls"]
asym = cls == "asymmetric"; rows = []
rows.append(dict(test="genes per CNV (TWAS span)", mirror=np.median(Ls[~asym]), asymmetric=np.median(Ls[asym]), p=stats.mannwhitneyu(Ls[asym], Ls[~asym]).pvalue, note="medians; Mann-Whitney"))
def resid_on_n(v):
    X = np.column_stack([np.ones(len(v)), Ls]); b = np.linalg.lstsq(X, v, rcond=None)[0]; return v - X @ b
for lab, f in [("mean z per gene", lambda S: S / Ls), ("summed z adjusted for gene count", resid_on_n)]:
    obs = f(x_obs); d_obs = obs[asym].mean() - obs[~asym].mean()
    null = np.array([ (lambda q: q[asym].mean() - q[~asym].mean())(f(Xp[b])) for b in range(Xp.shape[0])])
    rows.append(dict(test=f"class difference, {lab}", mirror=obs[~asym].mean(), asymmetric=obs[asym].mean(), p=(1 + np.sum(null <= d_obs)) / (len(null) + 1), note="circular permutation, one-sided (asymmetric lower)"))
k = cids != "7q11.23_WBS"
def circ_rho(xv_obs, xv_perm, kk):
    r = stats.spearmanr(xv_obs[kk], y[kk])[0]; yr = stats.rankdata(y[kk]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, xv_perm[:, kk]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
    return r, (1 + np.sum((Xr * yr[None, :]).mean(1) >= r)) / (xv_perm.shape[0] + 1)
for lab, xo, xp in [("summed z", x_obs, Xp), ("mean z per gene", x_obs / Ls, Xp / Ls[None, :])]:
    r, p = circ_rho(xo, xp, k); rows.append(dict(test=f"pooled dosage correlation without Williams, {lab}", mirror=np.nan, asymmetric=np.nan, rho=r, p=p, note="circular permutation"))
rx, ry, rn = (stats.rankdata(v[k]) for v in (x_obs, y, Ls.astype(float)))
ex = rx - np.polyval(np.polyfit(rn, rx, 1), rn); ey = ry - np.polyval(np.polyfit(rn, ry, 1), rn); pr = np.corrcoef(ex, ey)[0, 1]
rows.append(dict(test="pooled dosage correlation without Williams, partial Spearman on gene count", rho=pr, p=2 * stats.t.sf(abs(pr * np.sqrt((k.sum() - 3) / (1 - pr ** 2))), k.sum() - 3), note="t test"))
rows.append(dict(test="gene count vs dosage effect without Williams", rho=stats.spearmanr(Ls[k], y[k])[0], p=stats.spearmanr(Ls[k], y[k])[1], note="Spearman"))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_genecount_class_check.csv", index=False); pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 80); print(R.round(4).to_string(index=False))
