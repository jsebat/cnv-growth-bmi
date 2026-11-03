# Independent check of every number behind the TWAS concordance (Jon 2026-11-02: "double check all your numbers").
# 1. CNV set and classes: my 34 intervals (pathway_analysis cnv_genes_annotated.csv) vs Fig. 3's 33 loci (repo height_loci.csv).
# 2. Effect sizes: values used in the analysis vs Supplementary Table 3 (meta-analysis) vs Fig. 3's height_loci.csv.
# 3. Summed TWAS z: recomputed from the raw TWAS30 file for each CNV (gene-list and contiguous-span versions).
# 4. P-value stability: key circular-permutation tests rerun with a different seed.
import re, numpy as np, pandas as pd
from scipy import stats
HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ROOT = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly"
norm = lambda s: re.sub(r"[^0-9a-z]", "", str(s).lower())
genes = pd.read_csv(ROOT + r"\pathway_analysis\Final\cnv_genes_annotated.csv")
mine = genes.groupby("CNV_ID").Dose_response_class.first()
f3 = pd.read_csv(HERE + r"\fig3_inputs\height_loci.csv").set_index("CNV_ID")
t3 = pd.read_excel(ROOT + r"\NatureSubmission_04_13_2026\Sacks Supplementary Tables.xlsx", sheet_name="T3 CNV main effects and frequen")
t3["allele"] = t3.locus.str.extract(r"_(DEL|DUP)$")[0]; t3["key"] = t3.locus.str.replace(r"_(DEL|DUP)$", "", regex=True).map(norm)
used = pd.read_csv(HERE + r"\concordance_twas_by_cnv.csv"); used = used[used.trait == "height"].set_index("CNV_ID")
print("== 1. CNV sets")
print("in my set, not in Fig. 3:", sorted(set(mine.index) - set(f3.index)))
print("in Fig. 3, not in my set:", sorted(set(f3.index) - set(mine.index)))
cm = [(c, mine[c], f3.loc[c, "Dose_response_class"]) for c in mine.index if c in f3.index and mine[c] != f3.loc[c, "Dose_response_class"]]
print("class disagreements:", cm if cm else "none")
print("\n== 2. Effect sizes (height): used vs Supplementary Table 3 vs Fig. 3 input")
rows = []
for c in mine.index:
    k = norm(c); r = dict(CNV_ID=c)
    for a in ["DEL", "DUP"]:
        x = t3[(t3.key == k) & (t3.allele == a)]
        r[f"{a}_T3"] = x.betaHeight.iloc[0] if len(x) else np.nan
        r[f"{a}_used"] = used.loc[c, f"{a}_height"] if c in used.index else np.nan
        r[f"{a}_fig3"] = f3.loc[c, "Deletion_beta" if a == "DEL" else "Duplication_beta"] if c in f3.index else np.nan
    rows.append(r)
E = pd.DataFrame(rows).set_index("CNV_ID")
for a in ["DEL", "DUP"]:
    d1 = (E[f"{a}_used"] - E[f"{a}_T3"]).abs(); d2 = (E[f"{a}_T3"] - E[f"{a}_fig3"]).abs()
    print(f"{a}: max |used - T3| = {d1.max():.2e}; max |T3 - Fig3| = {d2.max():.2e}; Fig. 3 mismatches > 0.001: {list(E.index[d2 > 1e-3])}")
print(E.loc[[c for c in E.index if c.startswith("22q")]].round(4).to_string())
print("\n== 3. Summed TWAS z recomputed from the raw file")
tw = pd.read_csv(HERE + r"\twas30\height.txt", sep=r"\s+"); gz = tw.groupby("GENE")["TWAS.Z"].mean()
g = tw.groupby("GENE").agg(chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].str.replace("chr", "").replace("X", "23").astype(int); g = g.sort_values(["c", "pos"]).reset_index(drop=True); g["z"] = g.GENE.map(gz)
pos = dict(zip(g.GENE, g.index))
pan = pd.read_csv(HERE + r"\twas_dosage_panel.csv").set_index("CNV_ID")
mx_list, mx_span = 0, 0
for c, sub in genes.groupby("CNV_ID"):
    s_list = gz.reindex(sub.Gene).dropna().sum()
    ix = sorted(pos[x] for x in sub.Gene if x in pos)
    s_span = g.z.iloc[ix[0]: ix[-1] + 1].sum() if ix else np.nan
    if c in used.index and not np.isnan(used.loc[c, "S"]): mx_list = max(mx_list, abs(s_list - used.loc[c, "S"]))
    if c in pan.index: mx_span = max(mx_span, abs(s_span - pan.loc[c, "summed_twas_z"]))
print(f"max |gene-list sum - stored| = {mx_list:.2e}; max |span sum - stored (panel)| = {mx_span:.2e}")
print("dosage effect check (panel): max |(DUP-DEL)/2 - stored| =", f"{(((pan.DUP - pan.DEL) / 2) - pan.dosage_effect).abs().max():.2e}")
