# Non-overlapping CNV set (Jon 2026-11-02): drop intervals that span other analysed intervals so that every CNV is an independent
# observation. Primary set (Jon): 22q11.2 keep A-B and B-D (drop A-D, B-C, C-D); elsewhere drop the spanning interval:
# 16p11.2 BP1-BP3, 1q21.1 TAR BP1-BP3, 16p13.11 BP1-BP3. Sensitivity: finest 22q11.2 split (A-B, B-C, C-D; drop B-D).
# Re-runs: class difference in summed z (Welch t, Mann-Whitney, size-matched circular null), per gene; pooled and per-class dosage
# correlation (circular P), with and without Williams; one line vs separate class lines.
# Output: twas_nonoverlapping.csv
import numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]; ns = {}; exec(src, ns)
x_all, y_all, cls_all, cids_all, Xp_all, L_all = ns["x_obs"], ns["y"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["Ls"]
# Jon 2026-11-02 (revised): primary set keeps the spanning intervals and drops the subregions inside them:
# 22q11.2 A-D (drop A-B, B-C, C-D, B-D), 16p11.2 BP1-BP3 (drop BP2-BP3), 1q21.1 TAR BP1-BP3 and 16p13.11 BP1-BP3 (drop BP1-BP2, BP2-BP3).
# The earlier subregion sets are kept as sensitivity analyses.
SPAN = ["16p11.2_BP1-BP3", "1q21.1_TAR_BP1-BP3", "16p13.11_BP1-BP3"]
SUB = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3"]
SETS = {"primary (spanning intervals; 22q A-D)": SUB + ["22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"],
        "sensitivity (subregions; 22q A-B, B-D)": SPAN + ["22q11.21_A-D", "22q11.21_B_C", "22q11.21_C_D"],
        "sensitivity (subregions; 22q A-B, B-C, C-D)": SPAN + ["22q11.21_A-D", "22q11.21_B-D"]}
#  verify_numbers.py (2026-11-02) found 22q11.2 A-D is not one of Fig. 3's 33 loci (it came from the pathway-analysis gene file);
#  within Fig. 3's loci the spanning-interval rule keeps 22q11.2 A-B and B-D.
SETS["Fig. 3 loci, spanning intervals (22q A-B, B-D)"] = SUB + ["22q11.21_A-D", "22q11.21_B_C", "22q11.21_C_D"]
PANEL_SET = "primary (spanning intervals; 22q A-D)"
def ftest(r0, p0, r1, p1, n): return stats.f.sf(((r0 - r1) / (p1 - p0)) / (r1 / (n - p1)), p1 - p0, n - p1)
def rss(X, v): b = np.linalg.lstsq(X, v, rcond=None)[0]; e = v - X @ b; return float(e @ e)
rows = []
for setname, drop in SETS.items():
    keep = ~np.isin(cids_all, drop)
    x, y, cls, Xp, L, cids = x_all[keep], y_all[keep], cls_all[keep], Xp_all[:, keep], L_all[keep], cids_all[keep]
    a = cls == "asymmetric"
    # class difference in summed z and per gene
    for lab, xv, xp in [("summed z", x, Xp), ("mean z per gene", x / L, Xp / L[None, :])]:
        t = stats.ttest_ind(xv[a], xv[~a], equal_var=False); d_obs = xv[a].mean() - xv[~a].mean(); dn = xp[:, a].mean(1) - xp[:, ~a].mean(1)
        rows.append(dict(cnv_set=setname, test=f"class difference, {lab}", n_mirror=int((~a).sum()), n_asym=int(a.sum()), estimate=d_obs,
                         p_welch=t.pvalue, p_mannwhitney=stats.mannwhitneyu(xv[a], xv[~a]).pvalue, p_size_matched=(1 + np.sum(dn <= d_obs)) / (len(dn) + 1)))
    # dosage correlations
    def circ(k):
        r = stats.spearmanr(x[k], y[k])[0]; yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
        Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); Xr = (Xr - Xr.mean(1, keepdims=True)) / Xr.std(1, keepdims=True)
        return r, (1 + np.sum((Xr * yr[None, :]).mean(1) >= r)) / (Xp.shape[0] + 1)
    noW = cids != "7q11.23_WBS"
    for lab, k in [("dosage rho, all", np.ones(len(x), bool)), ("dosage rho, all without Williams", noW), ("dosage rho, mirror", ~a),
                   ("dosage rho, asymmetric", a), ("dosage rho, asymmetric without Williams", a & noW)]:
        r, p = circ(k); rows.append(dict(cnv_set=setname, test=lab, n_mirror=int((~a & k).sum()), n_asym=int((a & k).sum()), estimate=r, p_size_matched=p))
    # one line or two, without Williams
    k = noW; n = int(k.sum()); one = np.ones(n); xa = x[k]; aa = a[k].astype(float); yy = y[k]
    r1 = rss(np.column_stack([one, xa]), yy); r3 = rss(np.column_stack([one, xa, aa, xa * aa]), yy)
    rows.append(dict(cnv_set=setname, test="separate class lines vs one line, without Williams (F test)", n_mirror=int((~a & k).sum()), n_asym=int((a & k).sum()), p_welch=ftest(r1, 2, r3, 4, n)))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_nonoverlapping.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70); print(R.round(4).to_string(index=False))

# ---- panel inputs for the primary non-overlapping set (Theil-Sen line through all CNVs except Williams)
drop = SETS[PANEL_SET]; keep = ~np.isin(cids_all, drop)
P = pd.read_csv(HERE + r"\twas_dosage_panel.csv"); P = P[~P.CNV_ID.isin(drop)]; P.to_csv(HERE + r"\twas_dosage_panel_nonoverlap.csv", index=False)
k = (P.CNV_ID != "7q11.23_WBS").values; ts = stats.theilslopes(P.dosage_effect.values[k], P.summed_twas_z.values[k])
pr = R[R.cnv_set == PANEL_SET].set_index("test")
pd.DataFrame([dict(n=int(k.sum()), rho_all_noW=pr.loc["dosage rho, all without Williams", "estimate"], p_all_noW=pr.loc["dosage rho, all without Williams", "p_size_matched"],
                   rho_mirror=pr.loc["dosage rho, mirror", "estimate"], p_mirror=pr.loc["dosage rho, mirror", "p_size_matched"],
                   rho_asym_noW=pr.loc["dosage rho, asymmetric without Williams", "estimate"], p_asym_noW=pr.loc["dosage rho, asymmetric without Williams", "p_size_matched"],
                   p_separate_lines=pr.loc["separate class lines vs one line, without Williams (F test)", "p_welch"], theilsen_slope=ts[0], theilsen_intercept=ts[1])]
             ).to_csv(HERE + r"\twas_dosage_panel_nonoverlap_stats.csv", index=False)
print("panel inputs (non-overlapping) written")
