# Panel b (Jon 2026-11-02): mean summed TWAS z of the genes in mirror and in asymmetric CNVs, each against its size-matched null
# (circular permutation: each CNV moved to a random genomic position with its own number of genes; 10,000 draws). Same 23
# non-overlapping CNVs as panel a (spanning intervals kept; twas_nonoverlapping.py PANEL_SET). Also the class difference.
# Outputs: twas_class_null_draws.csv (null means per class), twas_class_null_stats.csv
import numpy as np, pandas as pd
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
src = open(HERE + r"\twas_dosage_panel_stats.py", encoding="utf-8").read().split("out = []")[0]; ns = {}; exec(src, ns)
x, cls, cids, Xp, L = ns["x_obs"], ns["cls"], np.array(ns["cids"]), ns["Xp"], ns["Ls"]
keep_ids = set(pd.read_csv(HERE + r"\twas_dosage_panel_nonoverlap.csv").CNV_ID)
k = np.array([c in keep_ids for c in cids])
rows, draws = [], {}
for grp in ["mirror", "asymmetric"]:
    m = k & (cls == grp); obs = x[m].mean(); null = Xp[:, m].mean(1); draws[grp] = null
    rows.append(dict(group=grp, n_cnv=int(m.sum()), median_genes=float(np.median(L[m])), observed_mean_summed_z=obs, null_mean=null.mean(), null_sd=null.std(),
                     null_q025=np.quantile(null, 0.025), null_q975=np.quantile(null, 0.975),
                     p_lower_one_sided=(1 + np.sum(null <= obs)) / (len(null) + 1), p_two_sided=min(1, 2 * min((1 + np.sum(null <= obs)), (1 + np.sum(null >= obs))) / (len(null) + 1))))
d_obs = rows[1]["observed_mean_summed_z"] - rows[0]["observed_mean_summed_z"]; d_null = draws["asymmetric"] - draws["mirror"]
rows.append(dict(group="asymmetric minus mirror", observed_mean_summed_z=d_obs, null_mean=d_null.mean(), null_sd=d_null.std(),
                 null_q025=np.quantile(d_null, 0.025), null_q975=np.quantile(d_null, 0.975), p_lower_one_sided=(1 + np.sum(d_null <= d_obs)) / (len(d_null) + 1),
                 p_two_sided=min(1, 2 * min(1 + np.sum(d_null <= d_obs), 1 + np.sum(d_null >= d_obs)) / (len(d_null) + 1))))
R = pd.DataFrame(rows); R.to_csv(HERE + r"\twas_class_null_stats.csv", index=False)
pd.DataFrame(draws).to_csv(HERE + r"\twas_class_null_draws.csv", index=False)
pd.set_option("display.width", 220); print(R.round(4).to_string(index=False))

# ---- per-CNV values for the strip plot (Jon's sketch): summed z standardised against each CNV's own size-matched null,
#      (observed - null mean) / null SD; class means of the standardised value tested against the same circular draws
mu, sd = Xp.mean(0), Xp.std(0); zs = (x - mu) / sd; Zp = (Xp - mu[None, :]) / sd[None, :]
pts = pd.DataFrame(dict(CNV_ID=cids, cls=cls, n_genes=L, summed_z=x, size_matched_z=zs))[k]
pts.to_csv(HERE + r"\twas_class_strip_points.csv", index=False)
srows = []
for grp in ["mirror", "asymmetric"]:
    m = k & (cls == grp); o = zs[m].mean(); nl = Zp[:, m].mean(1)
    srows.append(dict(group=grp, n_cnv=int(m.sum()), median=float(np.median(zs[m])), mean=o, p_two_sided=min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))))
ma, mm = k & (cls == "asymmetric"), k & (cls == "mirror")
o = zs[ma].mean() - zs[mm].mean(); nl = Zp[:, ma].mean(1) - Zp[:, mm].mean(1)
srows.append(dict(group="asymmetric minus mirror", mean=o, p_two_sided=min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))))
S2 = pd.DataFrame(srows); S2.to_csv(HERE + r"\twas_class_strip_stats.csv", index=False); print(); print(S2.round(4).to_string(index=False))

# ---- both candidate CNV sets for the Fig. 3 build (make_height_figure_v3.py --twas-set): all 32 CNVs' size-matched z, and the
#      class tests for each set
SUB = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3"]
SETS = {"ad": SUB + ["22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"],
        "fig3": SUB + ["22q11.21_A-D", "22q11.21_B_C", "22q11.21_C_D"]}
pd.DataFrame(dict(CNV_ID=cids, cls=cls, n_genes=L, summed_z=x, size_matched_z=zs)).to_csv(HERE + r"\twas_class_strip_points_all.csv", index=False)
out = []
for sname, drop in SETS.items():
    kk = ~np.isin(cids, drop)
    for grp in ["mirror", "asymmetric"]:
        m = kk & (cls == grp); o = zs[m].mean(); nl = Zp[:, m].mean(1)
        out.append(dict(twas_set=sname, group=grp, n_cnv=int(m.sum()), median=float(np.median(zs[m])), mean=o,
                        p_two_sided=min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))))
    ma, mm = kk & (cls == "asymmetric"), kk & (cls == "mirror"); o = zs[ma].mean() - zs[mm].mean(); nl = Zp[:, ma].mean(1) - Zp[:, mm].mean(1)
    out.append(dict(twas_set=sname, group="asymmetric minus mirror", mean=o, p_two_sided=min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))))
pd.DataFrame(out).to_csv(HERE + r"\twas_class_strip_stats_sets.csv", index=False); print(pd.DataFrame(out).round(4).to_string(index=False))
