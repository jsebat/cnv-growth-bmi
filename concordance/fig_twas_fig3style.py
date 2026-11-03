# TWAS concordance panels in the Fig. 3 style (Jon 2026-11-02: "if it's going into Fig. 3 it should use that format, as per the
# github"): matplotlib, palette and type sizes from fig3_mirror_asymmetric/.../code/main_figure.py (INK, LIGHT, BLUE, RED; mirror
# circles, asymmetric squares; bold 7-pt titles, 6.4-pt axis labels, 5.6-pt ticks, 10-pt panel letters).
#  a  CNV dosage effect on height, (DUP - DEL)/2, vs summed TWAS z of its genes; Theil-Sen line through all CNVs except Williams
#     (open symbol). Stats: twas_dosage_panel_nonoverlap_stats.csv.
#  b  per CNV, summed TWAS z standardised against its own size-matched null (random genomic windows with the same number of
#     genes; circular permutation), mirror vs asymmetric; bar = median. Stats: twas_class_strip_stats.csv.
# Same 23 non-overlapping CNVs in both panels (spanning intervals kept). Outputs: fig_twas_fig3style.{png,pdf}
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; RED = '#b03a2e'; BLUE = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fp = lambda p: f"{p:.2g}" if p >= 0.0001 else f"{p:.1e}"

d = pd.read_csv(HERE + r"\twas_dosage_panel_nonoverlap.csv"); ns = pd.read_csv(HERE + r"\twas_dosage_panel_nonoverlap_stats.csv").iloc[0]
pts = pd.read_csv(HERE + r"\twas_class_strip_points.csv"); ss = pd.read_csv(HERE + r"\twas_class_strip_stats.csv").set_index("group")
W = d.CNV_ID == "7q11.23_WBS"

fig, (ax, bx) = plt.subplots(1, 2, figsize=(5.4, 2.55), gridspec_kw=dict(width_ratios=[1.65, 1], wspace=.42), constrained_layout=False)
fig.subplots_adjust(left=.125, right=.955, top=.82, bottom=.21)

# ---- a
for cls, c, mk in [("mirror", BLUE, "o"), ("asymmetric", RED, "s")]:
    s = (d.cls == cls) & ~W
    ax.scatter(d.summed_twas_z[s], d.dosage_effect[s], s=14, color=c, marker=mk, edgecolor="white", lw=.4, zorder=4, label=cls)
ax.scatter(d.summed_twas_z[W], d.dosage_effect[W], s=16, marker="s", facecolor="white", edgecolor=RED, lw=.9, zorder=4)
xr = np.array([d.summed_twas_z[~W].min(), d.summed_twas_z[~W].max()])
ax.plot(xr, ns.theilsen_intercept + ns.theilsen_slope * xr, color=INK, lw=.9, zorder=3)
ax.axhline(0, color=LIGHT, lw=.6, zorder=1); ax.axvline(0, color=LIGHT, lw=.6, zorder=1)
lab = {"16p11.2_BP4_BP5": ("16p11.2 BP4\u2013BP5", 0, .05, "center"), "3q29": ("3q29", 0, -.07, "center"), "1q21.1": ("1q21.1", -.5, 0, "right"),
       "22q11.21_A-D": ("22q11.2 A\u2013D", -.4, .06, "right"), "7q11.23_WBS": ("7q11.23 (Williams)", .5, 0, "left"),
       "16p11.2_BP1-BP3": ("16p11.2 BP1\u2013BP3", 0, .06, "center")}
for cid, (t, dx, dy, ha) in lab.items():
    r = d[d.CNV_ID == cid].iloc[0]
    ax.text(r.summed_twas_z + dx, r.dosage_effect + dy, t, size=5.0, ha=ha, va="center", color=BLUE if r.cls == "mirror" else RED)
ax.text(.02, .98, f"All CNVs except Williams: \u03c1 = {ns.rho_all_noW:.2f}, P = {fp(ns.p_all_noW)}\n"
                  f"mirror \u03c1 = {ns.rho_mirror:.2f}; asymmetric \u03c1 = {ns.rho_asym_noW:.2f} (without Williams)\n"
                  f"separate lines for the two classes: P = {fp(ns.p_separate_lines)}", transform=ax.transAxes, size=5.2, va="top", color=INK, linespacing=1.35)
ax.set_xlabel("summed TWAS z of genes in the CNV\n(common-variant predicted expression; + = higher expression, taller)", size=6.4); ax.set_ylabel("CNV dosage effect on height (SD)\n(DUP \u2212 DEL) / 2", size=6.4)
ax.tick_params(labelsize=5.6, length=2, pad=1.5); ax.set_ylim(-.5, 1.08)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=5.4, loc="lower right", handletextpad=.3, markerscale=1)
ax.set_title("Common-variant predicted expression (TWAS) of the genes\nin each CNV predicts its dosage effect on height", size=7, weight="bold")
fig.text(.004, .985, "a", size=10, weight="bold", va="top")

# ---- b (strip plot, Jon's sketch)
rng = np.random.default_rng(1)
for i, (cls, c, mk, xl) in enumerate([("mirror", BLUE, "o", "mir"), ("asymmetric", RED, "s", "asy")]):
    v = pts.size_matched_z[pts.cls == cls].to_numpy()
    bx.scatter(i + rng.uniform(-.13, .13, len(v)), v, s=12, color=c, marker=mk, edgecolor="white", lw=.4, zorder=4)
    bx.plot([i - .22, i + .22], [np.median(v)] * 2, color=INK, lw=2, zorder=5, solid_capstyle="butt")
    bx.text(i, .95, f"P = {fp(ss.loc[cls, 'p_two_sided'])}", size=5.4, ha="center", va="bottom", color=INK)
bx.axhline(0, color=LIGHT, lw=.6, zorder=1)
bx.text(.5, 1.25, f"asy vs mir: P = {fp(ss.loc['asymmetric minus mirror', 'p_two_sided'])}", size=5.2, ha="center", va="bottom", color=GREY)
bx.set_xticks([0, 1]); bx.set_xticklabels(["mir", "asy"], size=6.4); bx.set_xlim(-.6, 1.6); bx.set_ylim(-1.9, 1.5)
bx.set_ylabel("summed TWAS z vs gene-count-\nmatched null (SD)", size=6.4)
bx.tick_params(labelsize=5.6, length=2, pad=1.5); bx.spines[["top", "right"]].set_visible(False)
bx.set_title("CNV genes vs random windows\nwith equal gene count", size=7, weight="bold")
fig.text(.6, .985, "b", size=10, weight="bold", va="top")

fig.savefig(HERE + r"\fig_twas_fig3style.png", dpi=400); fig.savefig(HERE + r"\fig_twas_fig3style.pdf")
print("written")
