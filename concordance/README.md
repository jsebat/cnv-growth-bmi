# CNV effect concordance with common- and rare-variant gene effects (Supplementary Note 4, Supplementary Table 29)

JS 2026-11-02, reviewer request (TWAS concordance), placed after the single-hit model results (Fig. 3).

Scripts (run in this order, Python 3 with pandas/scipy):
1. `concordance_twas_exome.py` - per-CNV summed TWAS z (TWAS30, Mancuso 2017) vs meta-analysed DEL/DUP effects (Sacks Supp T3),
   all / mirror / asymmetric CNVs (classes from pathway_analysis/Final/cnv_genes_annotated.csv); gene-label permutation;
   exome genes inside CNVs (Kosmicki 2026 Supp T6 height; Akbari 2021 Table 1 BMI).
2. `twas_genomewide_direction.py` - is higher predicted expression generally associated with greater height/BMI (no:
   mean z -0.03, P 0.32); gene-count adjustment.
3. `twas_circular_permutation.py` - circular genomic permutation (locus clusters shifted as units, gene count preserved).

Data: `twas30/` (public, github.com/bogdanlab/TWAS30), `exome/kosmicki_2026_supp_tables.xlsx` (medRxiv supplement).
Per-allele (per_cnv_height_concordance.py): 22/32 height-associated alleles concordant; 20/24 where |summed z| >= 1 (Supp T29 sheet I).
Results: height concordance in mirror CNVs (circular P 0.002 DEL, 0.007 DUP), none in asymmetric CNVs or for BMI.
Rare variants: only 3 exome-wide height genes (XYLT1, CNNM3, POR) and no BMI genes fall in the CNVs; all 3 direction-concordant
with the deletion. A full rare-variant test needs all-gene burden statistics (Genebass; requester-pays GCS bucket).
4. `twas_per_cnv_sign_and_fdr.py` - per-CNV sign test (no single CNV significant; min P 0.06) and TWAS genes at 5% FDR inside each
   CNV with concordance flags for the deletion and the duplication (Supp Table 29 sheets G, H). Height: 13/18 gene-CNV pairs
   concordant in mirror CNVs vs 7/14 in asymmetric; 7q11.23 WBS 4/5 significant genes discordant.
5. `twas_asymmetric_offset.py` - Jon's hypothesis: remove the asymmetric duplication-axis offset (c = -0.31 SD with mirror slope
   0.351; -0.28 free line). Result: concordance gets weaker, not stronger (pooled DUP rho 0.61 -> 0.38-0.46; asym DUP signs 9/12 ->
   7-8/12), because the asymmetric CNVs' own genes lean height-lowering. `twas_class_mean_test.py`: mean summed z asym -3.4 vs mirror
   +0.55 (10/12 asym negative), circular P = 0.09 (per gene 0.16) - a trend only. Supp T29 sheets J, K.
6. Single panel for a main figure: `twas_dosage_panel_stats.py` (per CNV, (DUP-DEL)/2 vs summed z; mirror rho 0.64, circular P 0.004;
   asymmetric 0.22, P 0.26) and `fig_twas_dosage_panel.R` -> fig_twas_dosage_panel.{png,pdf}. Draft text: twas_results_text_and_legend.md.


## BMI dosage analysis (Supplementary Note 4, Supplementary Table 29B-C)

`twas_dosage_bmi.py`: the Fig. 2E test with BMI (dosage effect on BMI vs summed BMI TWAS z, same non-overlapping CNVs, circular
permutation; height rerun as a check). BMI is null (rho -0.10, P 0.66, 22 CNVs); Williams is not an outlier for BMI, and
`twas_dosage_bmi_influence.csv` gives each CNV's robust residual, Cook's distance and leave-one-out rho.
