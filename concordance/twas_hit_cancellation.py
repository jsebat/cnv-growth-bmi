# Does the dosage effect (DUP - DEL)/2 cancel the single-hit shift? (Jon 2026-11-02)
# A shift that moves deletion and duplication effects in the same direction cancels in DUP - DEL and is kept in DUP + DEL.
#  1. Common component (DUP + DEL)/2: asymmetric vs mirror CNVs (the single hit should make it negative in asymmetric CNVs).
#  2. Dosage component (DUP - DEL)/2: asymmetric vs mirror CNVs after accounting for gene content (offset of asymmetric CNVs from
#     the line fitted to all CNVs; least squares with summed TWAS z as covariate), with and without Williams.
#  3. Seed check: pooled dosage correlation (all CNVs except Williams), circular permutation with a different seed.
# For the two candidate panel sets (22q11.2 A-D; Fig. 3 loci with 22q11.2 A-B and B-D). Output: twas_hit_cancellation.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]
ns = {}; exec(src.replace("default_rng(2026)", "default_rng(99)"), ns)
x, y, cls, cids, Xp, S = ns["x_obs"], ns["y"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["S"]
common = np.array([(S.loc[c, "DUP_effect"] + S.loc[c, "DEL_effect"]) / 2 for c in cids])
SUB = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3"]
SETS = {"22q11.2 A-D": SUB + ["22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"],
        "Fig. 3 loci (22q11.2 A-B, B-D)": SUB + ["22q11.21_A-D", "22q11.21_B_C", "22q11.21_C_D"]}
def ols_p(X, v, j):
    b = np.linalg.lstsq(X, v, rcond=None)[0]; e = v - X @ b; n, p = X.shape; s2 = e @ e / (n - p)
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[j, j]); return b[j], 2 * stats.t.sf(abs(b[j] / se), n - p)
rows = []
for setname, drop in SETS.items():
    keep = ~np.isin(cids, drop)
    for wlab, kk in [("all", keep), ("without Williams", keep & (cids != "7q11.23_WBS"))]:
        a = (cls == "asymmetric")[kk]; m = ~a
        c_a, c_m = common[kk][a], common[kk][m]
        rows.append(dict(cnv_set=setname, sample=wlab, component="common (DUP + DEL)/2 [single hit kept]", mean_asym=c_a.mean(), mean_mirror=c_m.mean(),
                         difference=c_a.mean() - c_m.mean(), p=stats.mannwhitneyu(c_a, c_m).pvalue, test="Mann-Whitney", n_asym=int(a.sum()), n_mirror=int(m.sum())))
        X = np.column_stack([np.ones(kk.sum()), x[kk], a.astype(float)]); b, p = ols_p(X, y[kk], 2)
        rows.append(dict(cnv_set=setname, sample=wlab, component="dosage (DUP - DEL)/2 [single hit cancels], offset given gene content", mean_asym=y[kk][a].mean(),
                         mean_mirror=y[kk][m].mean(), difference=b, p=p, test="least squares, asymmetric offset adjusted for summed TWAS z", n_asym=int(a.sum()), n_mirror=int(m.sum())))
    k = keep & (cids != "7q11.23_WBS")
    r = stats.spearmanr(x[k], y[k])[0]; yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
    rows.append(dict(cnv_set=setname, sample="without Williams", component="seed check (seed 99): pooled dosage rho", difference=r,
                     p=(1 + np.sum((Xr * yr[None, :]).mean(1) >= r)) / (Xp.shape[0] + 1), test="circular permutation"))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_hit_cancellation.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(R[["cnv_set", "sample", "component", "mean_asym", "mean_mirror", "difference", "p"]].round(4).to_string(index=False))
