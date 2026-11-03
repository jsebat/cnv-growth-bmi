# Supplementary Tables 24 (update), 27, 28 and 29 into the working workbook (JS 2026-11-02).
#  T24  UK Biobank: body size at 10 vs adult BMI, menarche, genetic instruments (Supplementary Fig. 5; fig_supp5_combined_tables.csv)
#  T27  UK Biobank Mendelian randomization of childhood body mass on adult height (supp_table_S27_..., from ukbb_mr_childhood_obesity.R)
#  T28  Polygenic scores in the clinical cohorts: reduced form and adjustment for measured BMI (T28a.csv, T28b.csv)
#  T29  CNV-TWAS concordance (Fig. 2E-F; BMI in part C): per CNV and tests (T29a_per_cnv.csv, T29b_tests.csv; full detail in
#       supp_table_S29_cnv_twas_exome_concordance.xlsx)
# Layout as 37_supp_table_agebands.R: title row 1, description row 2, table from row 4. A backup of the workbook is taken first.
suppressMessages({ library(data.table); library(openxlsx) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
MM <- file.path(HERE, "manuscript_materials_2026-10-24"); ST <- file.path(MM, "supplementary_tables")
wbf <- file.path(MM, "Sacks Supplementary Tables v2.xlsx"); wb <- loadWorkbook(wbf)
hs <- createStyle(textDecoration = "bold")
put <- function(sheet, title, desc, blocks) {
  if (sheet %in% names(wb)) removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow = 1); addStyle(wb, sheet, hs, rows = 1, cols = 1)
  writeData(wb, sheet, desc, startRow = 2)
  r <- 4
  for (b in blocks) {
    if (!is.null(b$head)) { writeData(wb, sheet, b$head, startRow = r); addStyle(wb, sheet, hs, rows = r, cols = 1); r <- r + 1 }
    writeData(wb, sheet, b$data, startRow = r, headerStyle = hs); r <- r + nrow(b$data) + 3
  }
  setColWidths(wb, sheet, 1:12, "auto")
}
t24 <- fread(file.path(HERE, "fig_supp5_combined_tables.csv"))
put("T24 UKBB body size at 10", "Supplementary Table 24. UK Biobank: childhood body size, pubertal timing and genetic instruments of adult height in CNV carriers (Supplementary Fig. 5)",
    "Panel a: adult height (SD) per ordered step of recalled body size at age 10 and of adult BMI group in the same individuals. b: age at menarche, difference from non-carriers (years), adjusted for age and ten principal components. c: adult height (SD) per year of later menarche within genotype. d: adult height (SD) per s.d. of polygenic score within genotype, adjusted for the height score and principal components; first_stage = recalled body size at 10 per s.d. of score. P_interaction = difference from non-carriers. Participants of European genetic ancestry.",
    list(list(data = t24)))
t27 <- fread(file.path(ST, "supp_table_S27_ukbb_mr_childhood_obesity_adult_height.csv"))
put("T27 UKBB MR childhood BMI", "Supplementary Table 27. UK Biobank Mendelian randomization of childhood body mass and pubertal timing on adult height in CNV carriers (Supplementary Fig. 5d,e)",
    "For each genotype and polygenic score (childhood obesity, childhood BMI, adult BMI): first-stage associations of the score with recalled body size at age 10 and with BMI, the within-genotype association of the score with adult height (adjusted for the height score), its interaction with genotype, and Wald ratios with delta-method standard errors. 22q11.2 deletion carriers (n = 5) omitted.",
    list(list(head = "A. Childhood obesity, childhood BMI and adult BMI scores: first stage, reduced form and Wald ratios", data = t27),
         list(head = "B. Age at menarche on adult height (women): menarche polygenic score as instrument; first stage in all women (Supplementary Fig. 5e)", data = fread(file.path(HERE, "mr_menarche_height_ukbb.csv"))),
         list(head = "C. Adiposity scores and age at menarche (women): years of menarche per s.d. of score", data = fread(file.path(HERE, "ukbb_mr_puberty_pathway_A.csv"))),
         list(head = "D. Two-step mediation (women): indirect effect of each adiposity score on adult height through age at menarche, and total effect", data = fread(file.path(HERE, "ukbb_mr_puberty_pathway_B.csv"))),
         list(head = "E. Multivariable model: adult height per s.d. of the childhood obesity, menarche and adult BMI scores entered together (height score as covariate)", data = fread(file.path(HERE, "ukbb_mr_puberty_pathway_C.csv"))),
         list(head = "Note", data = data.frame(note = "The menarche (Guan et al. 2025), adult BMI (Yengo 2018) and height (Yengo 2022) GWAS include UK Biobank participants, so these instruments are not independent of the UK Biobank sample; the childhood obesity (Bradfield 2019) and childhood BMI (Vogelezang 2020) GWAS do not."))))
put("T28 Clinical PGS reduced form", "Supplementary Table 28. Polygenic scores for childhood obesity, childhood BMI and adult BMI in the clinical cohorts (G2MH carriers versus SSC controls)",
    "Part A: outcome ~ score x carrier + sex + age + height score, for height by age band, head circumference, SRS and SCQ. Part B: within carriers, the score effect with and without adjustment for measured BMI. All carriers are from G2MH and all controls from SSC, so interactions are also cohort contrasts.",
    list(list(head = "A. Score effects in carriers and controls", data = fread(file.path(ST, "T28a.csv"))),
         list(head = "B. Within carriers, with and without adjustment for measured BMI", data = fread(file.path(ST, "T28b.csv")))))
put("T29 CNV TWAS concordance", "Supplementary Table 29. Concordance of CNV effects on height and BMI with common-variant gene effects (TWAS) and exome-wide genes (Fig. 2E-F, Supplementary Note 4)",
    "Summed TWAS z = sum over the genes in the CNV (contiguous run of TWAS genes) of the gene-level height TWAS z (mean over expression panels; Mancuso et al. 2017). Dosage effect = (DUP - DEL)/2. Fig. 2E-F and the BMI analysis use non-overlapping CNVs (spanning intervals kept: 22q11.2 A-D, 16p11.2 BP1-BP3, 1q21.1 TAR BP1-BP3, 16p13.11 BP1-BP3); Williams (7q11.23) is excluded from the line and correlations. Full detail, sensitivity analyses and BMI in supp_table_S29_cnv_twas_exome_concordance.xlsx.",
    list(list(head = "A. Per CNV", data = fread(file.path(ST, "T29a_per_cnv.csv"))),
         list(head = "B. Tests", data = fread(file.path(ST, "T29b_tests.csv"))),
         list(head = "C. BMI, per CNV (as Fig. 2E): dosage effect on BMI, (DUP - DEL)/2, against summed BMI TWAS z; robust residual from the Theil-Sen line (MAD z), Cook's distance and Spearman rho with the CNV left out", data = fread(file.path(ST, "T29c_bmi_per_cnv.csv"))),
         list(head = "D. Constrained genes and TWAS, genome-wide: share of genes with an expression model and TWAS signal by GeneBayes s_het tier (Zeng et al. 2024); last row per trait, Spearman correlation of |z| with s_het (estimate in mean_abs_z, P in mean_z)", data = fread(file.path(ST, "T29d_constraint_genomewide.csv"))),
         list(head = "E. Per CNV (Fig. 2E set): constrained-gene burden (sum of s_het; genes with s_het >= 0.03 and >= 0.18), mean TWAS z of constrained and unconstrained genes, dosage effects", data = fread(file.path(ST, "T29e_constraint_per_cnv.csv"))),
         list(head = "F. Joint model of the dosage effect on constrained-gene burden and summed TWAS z: partial Spearman of TWAS given burden with circular-permutation P (10,000), least-squares standardised coefficients, R2, and the Williams-Beuren residual with and without burden", data = fread(file.path(ST, "T29f_constraint_models.csv"))),
         list(head = "G. Height, all 23 CNVs: which covariate accounts for the Williams-Beuren outlier? Joint least-squares model of the dosage effect on the covariate and summed TWAS z (standardised), the Williams standardised residual, and the partial Spearman of dosage with summed TWAS z given the covariate (circular permutation). Constrained burden, not gene count, reduces the Williams residual; the rank test is significant with or without any covariate", data = fread(file.path(ST, "T29g_williams_covariate.csv"))),
         list(head = "H. Brain TWAS: FUSION with PsychENCODE cortex weights on BMI and height GWAS (Tan et al. 2024, population estimates); genome-wide rows give genes with z (n_cnv), mean |z| (spearman_rho) and genes with |z| > 4.5 (partial_rho_given_burden); other rows the Fig. 2E test (circular permutation) without and with constrained-gene burden", data = fread(file.path(ST, "T29h_brain_twas.csv"))),
         list(head = "I. Fig. 2E test with the TWAS effect-size estimate (ALPHA; Mancuso et al. 2017) in place of z", data = fread(file.path(ST, "T29i_alpha.csv")))))
put("T30 Extended CNV models", "Supplementary Table 30. Extended two-component models of reciprocal CNV effects on height and BMI (Supplementary Note 5)",
    "Mirror part: gene-level background as in Fig. 2K-O; + summed TWAS z along the mirror direction; + Poisson mirror hits (rate per constrained gene, GeneBayes s_het >= 0.03). Asymmetry part (lower-left for height, upper-right for BMI): none; continuous per gene (cont_G) or per constrained gene (cont_B); single hit; Poisson multi-hit per gene (multi_G) or per constrained gene (multi_B); single hit + continuous per constrained gene. Conditional AIC and BIC with the mirror base fixed; dAIC and dBIC relative to the best model within each trait. Parameters: c TWAS coefficient (per s.d. of summed z); lh, h mirror-hit rate and size; qc, bc, tc continuous term; qs, bs, ts single hit; la, ba, ta multi-hit (angles in degrees).",
    list(list(head = "A. Model grid, Mancuso et al. TWAS (21 architectures per trait)", data = fread(file.path(ST, "T30a_model_grid.csv"))),
         list(head = "B. Posterior number of mirror hits per locus in the models that include them", data = fread(file.path(ST, "T30b_mirror_hit_posteriors.csv")))))
# keep table order: T1 ... T26, T27, T28, T29, T30, notes last
nm <- names(wb); num <- suppressWarnings(as.numeric(sub("^T([0-9]+).*", "\\1", nm))); num[is.na(num)] <- 999
num[nm == "T19b Medication prevalence"] <- 19.5
worksheetOrder(wb) <- order(num)
saveWorkbook(wb, wbf, overwrite = TRUE)
cat("written:", paste(names(wb)[worksheetOrder(wb)], collapse = " | "), "\n")
