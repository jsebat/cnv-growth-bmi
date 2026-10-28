# Supplementary Tables 25 and 26 for Figure 3 (mirror and asymmetric CNV effects on height). JS 2026-10-28.
#  25  (a) the 33 reciprocal loci: class, size, deletion and duplication effects, perpendicular offset from the mirror line
#          and the probability of belonging to the displaced group (measurement-error-aware mixture);
#      (b) every test in Figure 3: bimodality (A-B), mirror vs asymmetric comparisons (C-F), effect vs gene count (G-H),
#          dose-response model comparison and the three models' simulated bimodality and burden predictions (I-K).
#  26  sensitivity analysis excluding the 9 loci with neither a deletion nor a duplication effect at p < 0.05
#      (tables/TableS_sensitivity_excl_noeffect_p05.csv, from make_supp_table_sensitivity.py).
# Reads only outputs already written by the Figure 3 pipeline (see ../README.md); writes tables/supp_table_S25*.csv,
# tables/supp_table_S26*.csv and, when the growth-analysis workbook is present, sheets T25 and T26 of
# "Sacks Supplementary Tables v2.xlsx".
suppressMessages({ library(data.table); library(openxlsx) })
args <- commandArgs(FALSE); me <- sub("^--file=", "", args[grep("^--file=", args)])
HERE <- normalizePath(file.path(dirname(me), ".."), winslash = "/")          # perp_offset_bimodality/
PROJ <- file.path(HERE, "..", "CNV_height_BMI_corrected_directional_models_v2")
ME <- file.path(HERE, "outputs", "me", "all_loci")

# ---- (a) loci
loci <- fread(file.path(PROJ, "height", "inputs", "sacks_height_all_loci_exact.csv"))
po <- fread(file.path(ME, "perp_offsets_per_locus.csv"))[trait == "Height"]
stopifnot(nrow(loci) == 33, nrow(po) == 33, setequal(loci$CNV_ID, po$CNV_ID))
x <- merge(loci, po[, .(CNV_ID, perp_offset, perp_offset_SE, P_lower_group)], by = "CNV_ID")
x[, `:=`(Deletion_P = 2 * pnorm(-abs(Deletion_beta / Deletion_SE)), Duplication_P = 2 * pnorm(-abs(Duplication_beta / Duplication_SE)))]
setorder(x, -Dose_response_class, perp_offset)
s25a <- x[, .(Locus = CNV_ID, Class = Dose_response_class, `Class posterior probability` = signif(Posterior_class_probability, 3),
              `Length (Mb)` = round(CNV_length_Mb, 3), Genes = Gene_count,
              `Deletion effect (SD)` = round(Deletion_beta, 3), `Deletion SE` = round(Deletion_SE, 3), `Deletion P` = signif(Deletion_P, 2),
              `Duplication effect (SD)` = round(Duplication_beta, 3), `Duplication SE` = round(Duplication_SE, 3), `Duplication P` = signif(Duplication_P, 2),
              `Offset from mirror line (SD)` = round(perp_offset, 3), `Offset SE` = round(perp_offset_SE, 3),
              `P(displaced group)` = round(P_lower_group, 3))]

# ---- (b) tests
row <- function(sec, test, stat, value, p = NA_real_, note = "") data.table(Section = sec, Test = test, Statistic = stat, Value = value, P = p, Note = note)
R <- list()
# A-B bimodality
bm <- fread(file.path(ME, "bimodality_test.csv"))[1]
R[[length(R) + 1]] <- row("A-B. Bimodality of the offset", "Mirror line, slope rho (fitted on the 21 mirror-class loci)", "rho", round(bm$rho_full_fit, 3))
R[[length(R) + 1]] <- row("A-B. Bimodality of the offset", "Measurement-error-aware two- vs one-group mixture (equal variance; main analysis)", "likelihood ratio", round(bm$LRT, 2), bm$p_bootstrap_nullB,
                          sprintf("parametric bootstrap, %d replicates, each locus with its own SE", 5000L))
R[[length(R) + 1]] <- row("A-B. Bimodality of the offset", "Displaced group: weight; mean offsets of the two groups (SD); latent SD", "weight; means; tau",
                          NA, NA, sprintf("weight %.2f; means %.3f and %.3f; tau %.3f; %d of 33 loci assigned to their published class",
                                          bm$minor_weight, bm$comp1_mean, bm$comp2_mean, bm$tau, bm$classes_agree))
bp <- fread(file.path(HERE, "outputs", "bimodality_test.csv"))[trait == "Height"]
for (i in seq_len(nrow(bp))) R[[length(R) + 1]] <- row("A-B. Bimodality of the offset",
  sprintf("Plain mixture, line fitted on %s, %s variance", sub("only", "loci", bp$mirror_fit[i]), sub("_var", "", bp$alternative[i])), "likelihood ratio",
  round(bp$LRT[i], 2), bp$p_bootstrap_nullB[i], sprintf("sensitivity; null A (no measurement error) P = %.4f", bp$p_bootstrap_nullA[i]))
# C-F mirror vs asymmetric
cmp <- function(v, lab) { a <- x[Dose_response_class == "asymmetric", get(v)]; m <- x[Dose_response_class == "mirror", get(v)]
  row("C-F. Asymmetric (n = 12) vs mirror (n = 21) loci", lab, "medians, asymmetric vs mirror", NA, t.test(a, m)$p.value,
      sprintf("medians %.3g vs %.3g; means %.3g vs %.3g; Welch t-test P (Mann-Whitney P = %.2g)", median(a), median(m), mean(a), mean(m), wilcox.test(a, m, exact = FALSE)$p.value)) }
R[[length(R) + 1]] <- cmp("CNV_length_Mb", "CNV length (Mb)")
R[[length(R) + 1]] <- cmp("Gene_count", "Gene count")
R[[length(R) + 1]] <- cmp("Deletion_beta", "Deletion effect (SD)")
R[[length(R) + 1]] <- cmp("Duplication_beta", "Duplication effect (SD)")
# G-H effect vs gene count
for (cl in c("mirror", "asymmetric")) for (v in c("Deletion_beta", "Duplication_beta")) {
  d <- x[Dose_response_class == cl]; ct <- suppressWarnings(cor.test(d[[v]], d$Gene_count, method = "spearman", exact = FALSE))
  R[[length(R) + 1]] <- row("G-H. Effect size vs gene count", sprintf("%s effect, %s loci (n = %d)", sub("_beta", "", v), cl, nrow(d)), "Spearman rho", round(unname(ct$estimate), 3), ct$p.value) }
# I-K model comparison
mc <- fread(file.path(PROJ, "height", "outputs", "height_axis_mixture_model_comparison.csv"))
keep <- c("Mirror alone" = "Mirror only", "Mirror + continuous shift, free angle" = "Mirror + continuous (generalized dysfunction)",
          "Mirror + single hit, free angle (buffering)" = "Mirror + sparse single hit")
for (k in names(keep)) { r <- mc[model == k]; stopifnot(nrow(r) == 1)
  R[[length(R) + 1]] <- row("I-K. Dose-response models fitted to the 33 loci", keep[[k]], "delta conditional AIC", round(r$delta_conditional_AIC, 2), NA,
    if (k == "Mirror alone") sprintf("delta BIC %.2f; delta LOO %.2f", r$delta_conditional_BIC, r$delta_conditional_LOO)
    else sprintf("delta BIC %.2f; delta LOO %.2f; fraction of genes q = %.4f; shift DEL %.3f, DUP %.3f SD at %.1f degrees",
                 r$delta_conditional_BIC, r$delta_conditional_LOO, r$q, r$shift_DEL, r$shift_DUP, r$shift_angle_deg)) }
ex <- fread(file.path(ME, "model_expected_bimodality_v2.csv")); sc <- fread(file.path(ME, "model_criteria_scores.csv"))
lab <- c(mirror = "Mirror only", open_continuous = "Mirror + continuous (generalized dysfunction)", open_single = "Mirror + sparse single hit")
for (k in names(lab)) { e <- ex[key == k]; s <- sc[key == k]
  R[[length(R) + 1]] <- row("I-K. Simulated datasets (2,000 per model)", sprintf("%s: bimodality of the simulated offsets", lab[[k]]), "median P vs unimodal null",
    e$median_p_vs_unimodal, NA, sprintf("median likelihood ratio %.2f; %.1f%% of datasets with P < 0.05; %.1f%% reach the observed ratio", e$median_LRT, 100 * e$frac_p_lt_05, 100 * e$P_LRT_ge_observed))
  R[[length(R) + 1]] <- row("I-K. Simulated datasets (2,000 per model)", sprintf("%s: predicted genome-wide burden effect on height", lab[[k]]), "duplication burden (cm/Mb)",
    round(s$burden_DUP_per_Mb, 4), NA, sprintf("deletion burden %.4f cm/Mb", s$burden_DEL_per_Mb)) }
mace <- fread(file.path(PROJ, "code", "mace_burden_reference.csv"))[trait == "height"]
R[[length(R) + 1]] <- row("I-K. Simulated datasets (2,000 per model)", "Published genome-wide burden (Mace et al. 2017)", "duplication burden (cm/Mb)",
  mace[burden_model == "duplication", beta], mace[burden_model == "duplication", p_value], sprintf("deletion burden %.3f cm/Mb (P = %.2g)", mace[burden_model == "deletion", beta], mace[burden_model == "deletion", p_value]))
s25b <- rbindlist(R)
s25b[, P := signif(P, 2)]

# ---- 26 sensitivity
s26 <- fread(file.path(HERE, "tables", "TableS_sensitivity_excl_noeffect_p05.csv"))

# ---- write
TB <- file.path(HERE, "tables")
fwrite(s25a, file.path(TB, "supp_table_S25a_fig3_loci.csv")); fwrite(s25b, file.path(TB, "supp_table_S25b_fig3_tests.csv")); fwrite(s26, file.path(TB, "supp_table_S26_fig3_sensitivity.csv"))
GROW <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/manuscript_materials_2026-10-24"
if (dir.exists(GROW)) {
  for (f in c("supp_table_S25a_fig3_loci.csv", "supp_table_S25b_fig3_tests.csv", "supp_table_S26_fig3_sensitivity.csv")) file.copy(file.path(TB, f), file.path(GROW, "supplementary_tables", f), overwrite = TRUE)
  wbf <- file.path(GROW, "Sacks Supplementary Tables v2.xlsx"); wb <- loadWorkbook(wbf); b <- createStyle(textDecoration = "bold")
  put <- function(sh, title, note, blocks) {
    if (sh %in% names(wb)) removeWorksheet(wb, sh); addWorksheet(wb, sh)
    writeData(wb, sh, title, startRow = 1); addStyle(wb, sh, b, rows = 1, cols = 1); writeData(wb, sh, note, startRow = 2); r <- 4
    for (bl in blocks) { writeData(wb, sh, bl[[1]], startRow = r); addStyle(wb, sh, b, rows = r, cols = 1)
      writeData(wb, sh, bl[[2]], startRow = r + 1, headerStyle = b); r <- r + nrow(bl[[2]]) + 3 }
    setColWidths(wb, sh, 1:14, "auto") }
  put("T25 Fig3 mirror asymmetric", "Supplementary Table 25. Mirror and asymmetric CNV effects on height (Fig. 3)",
      "Effects are SD of height per allele from the meta-analysis (Supplementary Table 3). Offset = signed perpendicular distance from the mirror line DUP = -rho x DEL fitted on the mirror-class loci. P(displaced group) from the measurement-error-aware two-group mixture. Two-sided P values, unadjusted.",
      list(list("(a) The 33 reciprocal loci", s25a), list("(b) Tests", s25b)))
  put("T26 Fig3 sensitivity", "Supplementary Table 26. Figure 3 analyses repeated without the 9 loci where neither the deletion nor the duplication effect reaches P < 0.05",
      "The mirror layer and all three models are refitted on the remaining 24 loci. The main figure uses all 33 loci.", list(list("", s26)))
  ord <- names(wb); num <- suppressWarnings(as.numeric(sub("^T(\\d+).*", "\\1", ord))); num[is.na(num)] <- 999; num[grepl("^T19b", ord)] <- 19.5
  worksheetOrder(wb) <- order(num); saveWorkbook(wb, wbf, overwrite = TRUE)
  cat("workbook sheets:", paste(names(wb)[order(num)], collapse = " | "), "\n")
}
print(s25a, nrows = 40); print(s25b[, .(Section = substr(Section, 1, 8), Test = substr(Test, 1, 60), Value, P, Note = substr(Note, 1, 70))], nrows = 60)
