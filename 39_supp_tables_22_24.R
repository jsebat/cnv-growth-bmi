# Supplementary Tables 22-24 (JS 2026-10-27), written to manuscript_materials_2026-10-24/supplementary_tables/*.csv and as sheets
# of "Sacks Supplementary Tables v2.xlsx".
#  22  Fig. 6 growth tests with UK Biobank adults (17_main_with_ukbb.R + 38_fig6_numbers.R): sample sizes, carrier mean shift in
#      children and adults, age-adjusted CNV x BMI (difference from controls in the BMI slope), five-group omnibus and
#      within-carrier heterogeneity tests.
#  23  Fig. 4 pathway analysis set: the 102 GO BP 2023 / KEGG 2021 gene sets tested, with the three tests and their corrections
#      (pathway_analysis/Final/results/analysis_set_102_pathways.csv).
#  24  UK Biobank: adult height by recalled comparative body size at age 10 and by adult BMI group, same individuals
#      (ukbb_bodysize/ukbb_bodysize10_vs_adult_trends.csv; Supplementary Fig. 5).
suppressMessages({ library(data.table); library(openxlsx) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
MM <- file.path(HERE, "manuscript_materials_2026-10-24"); ST <- file.path(MM, "supplementary_tables")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
# ---- 22
nums <- fread(file.path(HERE, "fig6_numbers.csv")); bea <- fread(file.path(HERE, "main_with_ukbb_bmi_effect_by_age.csv"))
cxb <- unique(bea[, .(outcome, group, cnv_x_bmi = est_main, cnv_x_bmi_se = se_main, cnv_x_bmi_p = p_main)])
t22 <- rbindlist(lapply(c("Height", "Head circumference"), function(o) rbindlist(lapply(c("No CNV", LOCI), function(g) {
  n <- nums[quantity == "n" & outcome == o & group == (if (g == "No CNV") "Controls (tests)" else g)]
  sc <- nums[quantity == "shift_child" & outcome == o & group == g]; sa <- nums[quantity == "shift_adult" & outcome == o & group == g]
  cx <- cxb[outcome == o & group == g]; om <- nums[quantity == "omnibus_contrast" & outcome == o & group == g]
  data.table(Outcome = o, Group = if (g == "No CNV") "Non-carrier controls" else g, `n clinical` = n$n_clinical, `n UK Biobank` = n$n_ukbb,
             `Shift vs controls, < 18 y (SD)` = if (nrow(sc)) round(sc$est, 3) else NA, `SE, < 18 y` = if (nrow(sc)) round(sc$se, 3) else NA, `P, < 18 y` = if (nrow(sc)) signif(sc$p, 2) else NA,
             `Shift vs controls, >= 18 y (SD)` = if (nrow(sa)) round(sa$est, 3) else NA, `SE, >= 18 y` = if (nrow(sa)) round(sa$se, 3) else NA, `P, >= 18 y` = if (nrow(sa)) signif(sa$p, 2) else NA,
             `BMI slope (controls) or CNV x BMI (carriers), SD per SD` = round(cx$cnv_x_bmi, 3), SE = round(cx$cnv_x_bmi_se, 3), P = signif(cx$cnv_x_bmi_p, 2),
             `CNV x BMI in the five-group model` = if (nrow(om)) round(om$est, 3) else NA, `SE (five-group)` = if (nrow(om)) round(om$se, 3) else NA, `P (five-group)` = if (nrow(om)) signif(om$p, 2) else NA)
}))))
tests <- nums[quantity %in% c("omnibus_5group", "omnibus_5group_band", "heterogeneity"),
              .(Outcome = outcome, Test = fcase(quantity == "omnibus_5group", "BMI slope differs among the five genotype groups (age-adjusted; Wald)",
                                                quantity == "omnibus_5group_band", paste0("Same, within age band ", group, " y (F; groups: ", groups, ")"),
                                                quantity == "heterogeneity", "BMI slope differs among the four CNVs, carriers only (Wald)"),
                Statistic = round(fifelse(is.na(F), chisq, F), 2), df = df, P = signif(p, 2))]
ctl <- fread(file.path(HERE, "main_with_ukbb_control_bmi_tests.csv"))   # controls: BMI effect and its age dependence (17_main_with_ukbb.R)
tests <- rbind(ctl[, .(Outcome = outcome, Test = "Controls: BMI effect varies with age (age x BMI tensor interaction)", Statistic = NA_real_, df = NA_integer_, P = signif(p_age_x_bmi, 2))], tests)
fwrite(t22, file.path(ST, "supp_table_S22_fig6_growth_tests_ukbb.csv")); fwrite(tests, file.path(ST, "supp_table_S22b_fig6_genotype_tests.csv"))
# ---- 23
pw <- fread("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/pathway_analysis/Final/results/analysis_set_102_pathways.csv")
t23 <- pw[, .(Library = library, Pathway = pathway, `Set size` = pathway_size, `Genes in CNV universe` = asymmetric_cnv_pathway_n_in_universe,
              `Asymmetric-class genes` = asymmetric_cnv_overlap_n, `Asymmetric OR` = round(asymmetric_cnv_odds_ratio, 2), `Asymmetric P (one-sided)` = signif(asymmetric_cnv_p, 2),
              `Locus-block permutation P` = signif(cnv_locus_perm_p_asymmetric_greater, 2), `Mirror-class genes` = mirror_cnv_overlap_n, `Mirror OR` = round(mirror_cnv_odds_ratio, 2), `Mirror P (one-sided)` = signif(mirror_cnv_p, 2),
              `Kosmicki genes in set` = kosmicki_overlap_n, `Kosmicki OR` = round(kosmicki_odds_ratio, 2), `Kosmicki P` = signif(kosmicki_p, 2), `Kosmicki q (BH)` = signif(q_C, 2), `Kosmicki Bonferroni` = bonf_C,
              `Height-decreasing / increasing` = paste0(direction_negative_n, " / ", direction_positive_n), `Direction OR` = round(direction_odds_ratio, 2), `Direction P` = signif(direction_fisher_p_one_sided, 2),
              `Convergent (asymmetric)` = convergent_asymmetric, `Convergent (mirror)` = convergent_mirror)][order(`Asymmetric P (one-sided)`)]
fwrite(t23, file.path(ST, "supp_table_S23_fig4_pathway_analysis_set.csv"))
# ---- 24
bs <- fread(file.path(HERE, "ukbb_bodysize", "ukbb_bodysize10_vs_adult_trends.csv"))
t24 <- bs[, .(Group = group, Exposure = measure, n, `Height change per step (SD)` = round(trend, 3), SE = round(se, 3), P = signif(p, 2),
              `Mutually adjusted, per step (SD)` = round(trend_adj, 3), `SE (adjusted)` = round(se_adj, 3), `P (adjusted)` = signif(p_adj, 2), `P vs non-carriers` = signif(p_interaction_vs_noCNV, 2))]
fwrite(t24, file.path(ST, "supp_table_S24_ukbb_bodysize10_adult_height.csv"))
# ---- workbook
wbf <- file.path(MM, "Sacks Supplementary Tables v2.xlsx"); wb <- loadWorkbook(wbf); bold <- createStyle(textDecoration = "bold")
put <- function(sh, title, note, tab, tab2 = NULL) {
  if (sh %in% names(wb)) removeWorksheet(wb, sh); addWorksheet(wb, sh)
  writeData(wb, sh, title, startRow = 1); writeData(wb, sh, note, startRow = 2); writeData(wb, sh, tab, startRow = 4, headerStyle = bold)
  if (!is.null(tab2)) writeData(wb, sh, tab2, startRow = nrow(tab) + 7, headerStyle = bold)
  setColWidths(wb, sh, 1:ncol(tab), "auto")
}
put("T22 Fig6 growth tests UKBB", "Supplementary Table 22. Height and head circumference in CNV carriers versus non-carrier controls, with UK Biobank adults (Fig. 6A-J)",
    "Children: external growth references; adults: Z relative to the non-carrier adults of the same cohort (Box-Cox, lambda from UK Biobank). Shift: carrier - control mean, y ~ CNV + s(age), separately below and above 18 y. CNV x BMI: age-adjusted difference between carriers and controls in the BMI slope, y ~ CNV + s(age, by CNV) + s(age, by BMI) + CNV x BMI (controls row: the control BMI slope, y ~ BMI + s(age)). Five-group model: controls and the four CNVs together, each CNV's departure from the control slope; tests of genotype differences below the table. Height controls: clinical non-carriers plus 50,000 random UK Biobank non-carriers; head circumference: clinical cohorts only.",
    t22, tests)
put("T23 Fig4 pathway analysis set", "Supplementary Table 23. The 102 gene sets tested in Fig. 4A-C", "GO Biological Process 2023 and KEGG 2021 Human sets with at least four genes in the CNV gene universe. Asymmetric and mirror tests: one-sided hypergeometric within the CNV gene universe (Fig. 4A plots the two-sided Fisher P of the same 2 x 2 table); locus-block permutation: 5,000 reassignments of the asymmetric label among the 34 intervals. Kosmicki: enrichment of the 207 exome height genes (176 mapped). Direction: height-decreasing genes among the set's Kosmicki genes relative to all mapped Kosmicki genes. Bonferroni threshold 0.05/102; BH q within the set. Convergent: nominal in all three tests.", t23)
put("T24 UKBB body size at 10", "Supplementary Table 24. UK Biobank: adult height by recalled body size at age 10 and by adult BMI group, in the same individuals (Supplementary Fig. 5)", "European-ancestry participants with height, BMI and a valid answer to field 1687 (comparative body size at age 10: thinner, about average, plumper). Height: Box-Cox Z within sex, adjusted for age, age^2 and ten principal components. Trend: linear change in height Z per ordered step of each exposure (thinner < average < plumper; adult BMI <= 30th, 30th-70th, >= 70th centile); mutually adjusted: both exposures in one model; P vs non-carriers: exposure x CNV interaction.", t24)
ord <- names(wb); num <- suppressWarnings(as.numeric(sub("^T(\\d+).*", "\\1", ord))); num[is.na(num)] <- 999; num[grepl("^T19b", ord)] <- 19.5
worksheetOrder(wb) <- order(num); saveWorkbook(wb, wbf, overwrite = TRUE)
cat("Tables 22-24 written; t22", nrow(t22), "rows + ", nrow(tests), "tests; t23", nrow(t23), "rows; t24", nrow(t24), "rows\n")
print(tests)
