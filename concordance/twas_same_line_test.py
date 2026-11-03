# Do mirror and asymmetric CNVs fall on the same line once Williams (7q11.23) is dropped? (Jon 2026-11-02)
# Dosage effect y = (DUP - DEL)/2 against summed TWAS z (contiguous span), per CNV. Nested least-squares models with and without
# Williams: (1) one line; (2) + class offset; (3) + class-specific slope. F tests, and HC3-robust P for the class terms.
# Also the pooled Spearman correlation (all CNVs on one line) with the circular permutation null, without Williams, raw and after
# removing the asymmetric offset estimated in model (2).
# Output: twas_same_line_test.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
d = pd.read_csv(HERE + r"\twas_dosage_panel.csv"); d["asym"] = (d.cls == "asymmetric").astype(float)

def ols(X, yv):
    b, *_ = np.linalg.lstsq(X, yv, rcond=None); e = yv - X @ b; return b, e, float(e @ e)
def ftest(rss0, df0, rss1, df1, n):     # df = number of parameters
    F = ((rss0 - rss1) / (df1 - df0)) / (rss1 / (n - df1)); return stats.f.sf(F, df1 - df0, n - df1)
def hc3_p(X, yv, j):
    b, e, _ = ols(X, yv); XtXi = np.linalg.inv(X.T @ X); h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    V = XtXi @ (X.T * (e / (1 - h)) ** 2) @ X @ XtXi; t = b[j] / np.sqrt(V[j, j]); return 2 * stats.t.sf(abs(t), len(yv) - X.shape[1])

rows = []
for lab, dd in [("all CNVs", d), ("without Williams", d[d.CNV_ID != "7q11.23_WBS"])]:
    n = len(dd); x = dd.summed_twas_z.values; a = dd.asym.values; yv = dd.dosage_effect.values; one = np.ones(n)
    X1 = np.column_stack([one, x]); X2 = np.column_stack([one, x, a]); X3 = np.column_stack([one, x, a, x * a])
    b1, _, r1 = ols(X1, yv); b2, _, r2 = ols(X2, yv); b3, _, r3 = ols(X3, yv); tss = float(((yv - yv.mean()) ** 2).sum())
    rows.append(dict(sample=lab, n=n, slope_mirror=b3[1], slope_asym=b3[1] + b3[3], slope_difference=b3[3],
                     p_slope_difference=ftest(r2, 3, r3, 4, n), p_slope_difference_HC3=hc3_p(X3, yv, 3),
                     common_slope=b2[1], offset_asym_common_slope=b2[2], p_offset=ftest(r1, 2, r2, 3, n), p_offset_HC3=hc3_p(X2, yv, 2),
                     p_any_class_difference=ftest(r1, 2, r3, 4, n), r2_one_line=1 - r1 / tss, r2_offset=1 - r2 / tss, r2_separate=1 - r3 / tss))
R = pd.DataFrame(rows)
# pooled Spearman + circular permutation without Williams, raw and offset-removed (set-up part of the panel script: x_obs, y, cls, cids, Xp)
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]
ns = {}; exec(src, ns)
x_obs, y, cls, cids, Xp = ns["x_obs"], ns["y"], ns["cls"], np.array(ns["cids"]), ns["Xp"]
k = cids != "7q11.23_WBS"; off = R.loc[R["sample"] == "without Williams", "offset_asym_common_slope"].item()
def circ(yv, kk):
    r = stats.spearmanr(x_obs[kk], yv[kk])[0]; yr = stats.rankdata(yv[kk]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, kk]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
    return r, (1 + np.sum((Xr * yr[None, :]).mean(1) >= r)) / 10001
y_adj = np.where(cls == "asymmetric", y - off, y)
for lab, yv in [("pooled raw", y), ("pooled offset removed", y_adj)]:
    r, p = circ(yv, k); R.loc[R["sample"] == "without Williams", f"{lab} rho"] = r; R.loc[R["sample"] == "without Williams", f"{lab} p_circular"] = p
R.to_csv(HERE + r"\twas_same_line_test.csv", index=False); pd.set_option("display.width", 200); print(R.round(4).T.to_string())
