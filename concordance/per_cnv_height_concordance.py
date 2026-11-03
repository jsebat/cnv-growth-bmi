# Per-CNV concordance for height-associated CNV alleles (JS 2026-11-02, Jon): for each reciprocal CNV allele with a height
# association, is the direction predicted by its genes' summed TWAS z (contiguous span, as in the circular test) the observed
# direction? DUP concordant if sign(S) = sign(effect); DEL concordant if sign(S) = -sign(effect). Thresholds: nominal P < 0.05
# and Bonferroni P < 0.05/94 (94 recurrent alleles in the study). Nested intervals are listed separately.
# Output: per_cnv_height_concordance.csv
import numpy as np, pandas as pd
from scipy import stats
S = pd.read_csv("twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"]
rows = []
for _, r in S.iterrows():
    for a, sg in [("DEL", -1), ("DUP", 1)]:
        e, p = r[f"{a}_effect"], r[f"{a}_p"]
        rows.append(dict(CNV_ID=r.CNV_ID, cls=r.cls, allele=a, effect_sd=e, p=p, summed_z=r.summed_z, n_genes_twas=r.n_genes_twas,
                         frac_genes_concordant=r[f"{a}_frac_concordant"], concordant=bool(np.sign(r.summed_z) == sg * np.sign(e))))
D = pd.DataFrame(rows); D.to_csv("per_cnv_height_concordance.csv", index=False)
pd.set_option("display.width", 200)
for lab, thr in [("nominal P < 0.05", 0.05), ("Bonferroni P < 0.05/94", 0.05 / 94)]:
    x = D[D.p < thr]
    print(f"\n== {lab}: {len(x)} CNV alleles")
    print(x.sort_values(["cls", "CNV_ID", "allele"])[["CNV_ID", "cls", "allele", "effect_sd", "p", "summed_z", "frac_genes_concordant", "concordant"]].round(4).to_string(index=False))
    for cls in ["all", "mirror", "asymmetric"]:
        y = x if cls == "all" else x[x.cls == cls]
        if len(y): print(f"  {cls}: {int(y.concordant.sum())}/{len(y)} concordant; binomial P = {stats.binomtest(int(y.concordant.sum()), len(y), 0.5, alternative='greater').pvalue:.3g}")

# ---- Jon 2026-11-02: add prediction strength and summary counts; written to Supplementary Table 29 sheet I
D["height_associated"] = np.where(D.p < 0.05 / 94, "Bonferroni", np.where(D.p < 0.05, "nominal", "no"))
D["gene_content_prediction"] = np.where(D.summed_z.abs() >= 1, "clear (|summed z| >= 1)", "weak (|summed z| < 1)")
D["predicted_direction"] = np.where(D.allele == "DUP", np.sign(D.summed_z), -np.sign(D.summed_z)).astype(int)
D = D[["CNV_ID", "cls", "allele", "effect_sd", "p", "height_associated", "n_genes_twas", "summed_z", "gene_content_prediction",
       "predicted_direction", "frac_genes_concordant", "concordant"]].sort_values(["cls", "CNV_ID", "allele"])
summ = []
x = D[D.p < 0.05]
for pred in ["all", "clear (|summed z| >= 1)", "weak (|summed z| < 1)"]:
    for cls in ["all", "mirror", "asymmetric"]:
        y = x if pred == "all" else x[x.gene_content_prediction == pred]
        y = y if cls == "all" else y[y.cls == cls]
        if len(y):
            k = int(y.concordant.sum())
            summ.append(dict(gene_content_prediction=pred, cnv_class=cls, height_associated_alleles=len(y), concordant=k,
                             fraction=k / len(y), binomial_p_one_sided=stats.binomtest(k, len(y), 0.5, alternative="greater").pvalue,
                             discordant_alleles="; ".join(f"{r.CNV_ID} {r.allele}" for _, r in y[~y.concordant].iterrows())))
Sm = pd.DataFrame(summ); D.to_csv("per_cnv_height_concordance.csv", index=False); Sm.to_csv("per_cnv_height_concordance_summary.csv", index=False)
ST = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\CNV_growth_trajectories\manuscript_materials_2026-10-24\supplementary_tables"
with pd.ExcelWriter(ST + r"\supp_table_S29_cnv_twas_exome_concordance.xlsx", mode="a", if_sheet_exists="replace") as w:
    note = pd.DataFrame({"I. Per-allele height concordance": [
        "Each reciprocal CNV allele: observed height effect and the direction predicted by its genes (duplication: sign of summed TWAS z; deletion: opposite sign).",
        "Summary: height-associated alleles (P < 0.05), split by whether the gene content gives a clear prediction (|summed z| >= 1). One-sided binomial against 0.5; nested intervals listed separately, so alleles are not independent.", ""]})
    note.to_excel(w, sheet_name="I per-allele height concordance", index=False, startrow=0)
    Sm.to_excel(w, sheet_name="I per-allele height concordance", index=False, startrow=5)
    D.to_excel(w, sheet_name="I per-allele height concordance", index=False, startrow=5 + len(Sm) + 3)
print(Sm.round(4).to_string(index=False))
