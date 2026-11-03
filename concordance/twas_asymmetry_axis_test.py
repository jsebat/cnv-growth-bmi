# Does gene content predict the departure from a mirror response? Perpendicular offset from the mirror line (Fig. 2A) vs summed
# TWAS z, and for comparison the position along the mirror line. (JS 2026-11-02)
# Same circular permutation as the dosage-axis test; all CNVs, mirror, asymmetric; with and without Williams.
# Output: twas_asymmetry_axis_test.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]; ns = {}; exec(src, ns)
x_obs, cls, cids, Xp, S = ns["x_obs"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["S"]
RHO = 0.351140004155472   # Fig. 3 mirror fit: DUP = -RHO * DEL; perpendicular offset from the mirror line (as in Fig. 2A)
asy_y = np.array([(S.loc[c, "DUP_effect"] + RHO * S.loc[c, "DEL_effect"]) / np.sqrt(1 + RHO ** 2) for c in cids])
along = np.array([(-S.loc[c, "DEL_effect"] + RHO * S.loc[c, "DUP_effect"]) / np.sqrt(1 + RHO ** 2) for c in cids])   # position along the mirror line
def circ(yv, kk, sided):
    r = stats.spearmanr(x_obs[kk], yv[kk])[0]; yr = stats.rankdata(yv[kk]); yr = (yr - yr.mean()) / yr.std()
    Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, kk]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True); nl = (Xr * yr[None, :]).mean(1)
    return r, (1 + np.sum(np.abs(nl) >= abs(r))) / 10001
rows = []
for lab, base in [("all CNVs", np.ones(len(cids), bool)), ("without Williams", cids != "7q11.23_WBS")]:
    for grp in ["all", "mirror", "asymmetric"]:
        k = base & (np.ones(len(cids), bool) if grp == "all" else cls == grp)
        for ax, yv in [("perpendicular offset (asymmetry)", asy_y), ("along mirror line (dosage)", along)]:
            r, p = circ(yv, k, 2); rows.append(dict(axis=ax, sample=lab, group=grp, n=int(k.sum()), spearman_rho=r, p_circular_two_sided=p))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_asymmetry_axis_test.csv", index=False); print(R.round(4).to_string(index=False))
