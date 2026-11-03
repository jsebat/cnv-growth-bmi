# Build Supplementary Table 29 from the concordance outputs (run after the four analysis scripts). JS 2026-11-02.
import pandas as pd

HERE = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\concordance"
ST = r"C:\Users\Jonathan Sebat\sebatlab Dropbox\Manuscripts\PGC_joint_submission\Molly\CNV_ASD_growth_analysis\CNV_growth_trajectories\manuscript_materials_2026-10-24\supplementary_tables"
r = lambda f: pd.read_csv(HERE + "\\" + f)
by = r("concordance_twas_by_cnv.csv")[["trait", "CNV_ID", "cls", "n_genes", "n_twas", "S", "n_sig", "S_sig", "DEL_height", "DEL_height_p", "DUP_height", "DUP_height_p", "DEL_bmi", "DEL_bmi_p", "DUP_bmi", "DUP_bmi_p"]]
by = by.rename(columns={"cls": "class", "n_genes": "genes in interval", "n_twas": "genes with TWAS", "S": "summed TWAS z (mean over panels)", "n_sig": "TWAS-significant genes", "S_sig": "summed z, significant genes"})
T = r("concordance_twas_tests.csv"); T = T[T.score == "S"].drop(columns=["score"]).merge(r("concordance_twas_permutation.csv")[["trait", "group", "allele", "p_perm_one_sided"]], on=["trait", "group", "allele"])
readme = ["Concordance of CNV effects with common-variant (TWAS) and rare-variant (exome) gene effects. See Supplementary Note 4.",
          "CNV classes (mirror, asymmetric) are those of Fig. 3, defined on the height dose-response; for BMI no CNV is asymmetric, so BMI is analysed in all CNVs only.",
          "A: per CNV interval, summed TWAS z of its genes and meta-analysed deletion and duplication effects (SD).",
          "B: correlation and sign-concordance tests (expected sign +1 for duplications, -1 for deletions); gene-label permutation, 10,000 shuffles.",
          "C: exome-wide height (Kosmicki et al.) and BMI (Akbari et al.) genes inside CNV intervals, paired with the deletion effect.",
          "D: genome-wide direction of TWAS z (block jackknife over 10-Mb blocks).",
          "E: tests adjusted for gene count, and with mean z per gene.",
          "F: circular genomic permutation (locus clusters moved as units, gene count preserved), 10,000 permutations.",
          "G: per-CNV sign test (genes on the contiguous run spanning the CNV).",
          "H: TWAS genes at 5% FDR inside each CNV, with concordance of their direction with the CNV deletion and duplication effects.",
          "I: per reciprocal CNV allele, observed height effect versus the direction predicted by its genes; summary of height-associated alleles (P < 0.05) by CNV class and by whether the gene content gives a clear prediction (|summed z| >= 1).",
          "J: asymmetric CNVs after removing their duplication-axis offset from the mirror line (mirror slope fixed at 0.351, or a free asymmetric line): sign concordance and pooled correlation with circular permutation.",
          "K: mean summed TWAS z of the genes in asymmetric versus mirror CNVs, against the circular permutation null.",
          "L: dosage effect (DUP - DEL)/2 per CNV vs summed TWAS z (the main-figure panel): Spearman and circular permutation by class, leave-one-out, Theil-Sen lines.",
          "M: one line or two? Nested least-squares models (one line; + class offset; + class slope), with and without Williams.",
          "N: position along the mirror line and perpendicular offset from it (Fig. 2A asymmetry) vs summed TWAS z.",
          "O: gene count and the class difference in summed z: expected shift from gene count, circular null, and each CNV against random windows with the same number of genes.",
          "P: the main tests repeated in non-overlapping CNVs: primary 23 CNVs keeping the spanning intervals (22q11.2 A-D, 16p11.2 BP1-BP3, 1q21.1 TAR BP1-BP3, 16p13.11 BP1-BP3); sensitivity with subregions instead (22q11.2 A-B and B-D, or A-B, B-C, C-D), including Welch, Mann-Whitney and size-matched tests of the class difference in summed and per-gene z.",
          "Q: per CNV, summed TWAS z standardised against its own size-matched null (Fig. 3 panel b), with class tests from the circular permutation.",
          "R: single-hit cancellation: common component (DUP + DEL)/2 and dosage component (DUP - DEL)/2 in asymmetric vs mirror CNVs, and a seed check of the pooled correlation. Also verify_numbers.py: all effect sizes match Supplementary Table 3 and the Fig. 3 inputs.",
          "Sources: TWAS30 (Mancuso et al. 2017), Kosmicki et al. 2026 medRxiv Supplementary Table 6, Akbari et al. 2021 Table 1."]
with pd.ExcelWriter(ST + r"\supp_table_S29_cnv_twas_exome_concordance.xlsx") as w:
    pd.DataFrame({"Supplementary Table 29": readme}).to_excel(w, sheet_name="README", index=False)
    by.to_excel(w, sheet_name="A per CNV", index=False); T.to_excel(w, sheet_name="B tests", index=False)
    r("concordance_exome_genes.csv").to_excel(w, sheet_name="C exome genes", index=False)
    r("twas_genomewide_direction.csv").to_excel(w, sheet_name="D genome-wide direction", index=False)
    r("concordance_twas_genecount_adjusted.csv").to_excel(w, sheet_name="E gene-count adjusted", index=False)
    r("concordance_twas_circular_permutation.csv").to_excel(w, sheet_name="F circular permutation", index=False)
    r("twas_per_cnv_sign_test.csv").to_excel(w, sheet_name="G per-CNV sign test", index=False)
    r("twas_fdr_genes_in_cnvs.csv").to_excel(w, sheet_name="H FDR 5% genes in CNVs", index=False)
    Sm = r("per_cnv_height_concordance_summary.csv"); Di = r("per_cnv_height_concordance.csv")
    pd.DataFrame({"I. Per-allele height concordance": ["Summary (height-associated alleles, P < 0.05); one-sided binomial against 0.5; nested intervals listed separately, so alleles are not independent."]}).to_excel(w, sheet_name="I per-allele height", index=False)
    Sm.to_excel(w, sheet_name="I per-allele height", index=False, startrow=3)
    Di.to_excel(w, sheet_name="I per-allele height", index=False, startrow=3 + len(Sm) + 3)
    r("twas_asymmetric_offset.csv").to_excel(w, sheet_name="J asymmetric offset removed", index=False)
    r("twas_class_mean_test.csv").to_excel(w, sheet_name="K asym vs mirror gene z", index=False)
    r("twas_dosage_panel.csv").to_excel(w, sheet_name="L dosage panel", index=False)
    r("twas_dosage_panel_stats.csv").to_excel(w, sheet_name="L dosage panel", index=False, startrow=36)
    r("twas_dosage_panel_loo.csv").to_excel(w, sheet_name="L dosage panel", index=False, startrow=44)
    r("twas_same_line_test.csv").to_excel(w, sheet_name="M one line or two", index=False)
    r("twas_asymmetry_axis_test.csv").to_excel(w, sheet_name="N along vs off mirror line", index=False)
    r("twas_genecount_null.csv").to_excel(w, sheet_name="O gene count null", index=False)
    r("twas_genecount_null_per_cnv.csv").to_excel(w, sheet_name="O gene count null", index=False, startrow=4)
    r("twas_nonoverlapping.csv").to_excel(w, sheet_name="P non-overlapping CNVs", index=False)
    r("twas_class_strip_stats.csv").to_excel(w, sheet_name="Q strip size-matched z", index=False)
    r("twas_class_strip_points.csv").to_excel(w, sheet_name="Q strip size-matched z", index=False, startrow=6)
    r("twas_hit_cancellation.csv").to_excel(w, sheet_name="R single hit cancellation", index=False)
print("S29 rebuilt; B rows:", len(T)); print(T[["trait", "group", "allele", "n_cnv", "spearman_rho", "p_perm_one_sided"]].round(3).to_string(index=False))
