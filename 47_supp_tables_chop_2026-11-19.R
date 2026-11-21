# Supplementary Tables 17, 18, 21 and 22 rebuilt on the CHOP / fill-gaps data (JS 2026-11-19, Jon: "Regenerate Supp Tables
# 17/18/21/22 and Supp Fig 6 with the new data"). Same content and layout as the v2 workbook sheets (37_supp_table_agebands.R,
# 39_supp_tables_22_24.R, 32_supp_asd_by_cnv.R), read from the new run outputs:
#   growth (T21, T22): run_chop_fillgaps_2026-11-19/expanse_v9e/ (17_main_with_ukbb.R + 38_fig6_numbers.R on Expanse,
#                      v9e_chop_genotype_bands; ANALYSIS_TABLE = analysis_table_chop_2026-11-19.rds, VISIT_RULE = fill_gaps)
#   autism (T17, T18): run_chop_fillgaps_2026-11-19/asd_v9b_random/ (32_supp_asd_by_cnv.R, ADD_CHOP = TRUE)
# Output: run_chop_fillgaps_2026-11-19/supplementary_tables/*.csv and a NEW workbook
#   manuscript_materials_2026-10-24/Sacks Supplementary Tables v3_chop.xlsx (= v2 with T17, T18, T21, T22 replaced); v2 is untouched.
suppressMessages({ library(data.table); library(openxlsx) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
RUN <- file.path(HERE, "run_chop_fillgaps_2026-11-19"); G <- file.path(RUN, "expanse_v9e"); A <- file.path(RUN, "asd_v9b_random")
ST <- file.path(RUN, "supplementary_tables"); dir.create(ST, showWarnings = FALSE)
MM <- file.path(HERE, "manuscript_materials_2026-10-24")
wbf2 <- file.path(MM, "Sacks Supplementary Tables v2.xlsx"); wbf3 <- file.path(MM, "Sacks Supplementary Tables v3_chop.xlsx")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
SRC <- " Data: clinical cohorts including the CHOP 22q11.2 A-D records (G2MH carriers recruited at CHOP removed), one visit per person chosen to fill sparse age bands (Methods), plus UK Biobank adults for height. Rebuilt 2026-11-19."
bold <- createStyle(textDecoration = "bold")

# ---- T21
x <- fread(file.path(G, "main_with_ukbb_bmi_slope_by_band.csv"))[group != "No CNV"]
x[, age_bin := factor(age_bin, c("2-9", "10-18", "19-39", "40-73"))][, outcome := factor(outcome, c("Height", "Head circumference"))]; setorder(x, outcome, group, age_bin)
t21 <- x[, .(Outcome = as.character(outcome), `CNV group` = group, `Age band (years)` = as.character(age_bin), `Carriers (n)` = n_carriers, `Controls (n)` = n_controls,
             `Carrier BMI slope` = round(slope_carriers, 3), `Carrier SE` = round(se_carriers, 3), `Control BMI slope` = round(slope_controls, 3), `Control SE` = round(se_controls, 3),
             `Carrier - control` = round(diff, 3), `Difference SE` = round(se_diff, 3), P = signif(p_diff, 2))]
gb <- fread(file.path(G, "main_with_ukbb_bmi_genotype_by_band.csv"))
t21b <- gb[, .(Outcome = outcome, `Age band (years)` = age_bin, `Test` = "BMI effect differs among controls and the CNVs (genotype x BMI, F test)", F = round(F, 2), df, P = signif(p, 2), `CNVs in the model` = groups)]
fwrite(t21, file.path(ST, "supp_table_S21_bmi_slope_by_ageband.csv")); fwrite(t21b, file.path(ST, "supp_table_S21b_genotype_by_ageband.csv"))
# JS 2026-11-19 (Jon: "this can probably just be a supplementary table"): 0-2 y band, carriers only (no controls under 2 y), from
# run_chop_fillgaps_2026-11-19/preview_0to2/preview_0to2.R (latest visit under 2 y per person; Z ~ BMI Z + sex + age within group;
# head circumference adds the Simons VIP infant visits of hc_visits.rds)
inf <- fread(file.path(RUN, "preview_0to2", "preview_0to2_bands.csv"))[order(outcome, group)]
t21c <- inf[, .(Outcome = outcome, `CNV group` = group, `Age band (years)` = "0-2", `Carriers (n)` = n, `Median age (years)` = round(median_age, 2),
                `Carrier BMI slope` = round(slope, 3), `Carrier SE` = round(se, 3), P = signif(p, 2))]
fwrite(t21c, file.path(ST, "supp_table_S21c_bmi_slope_0to2_carriers.csv"))

# ---- T22 (as 39_supp_tables_22_24.R)
nums <- fread(file.path(G, "fig6_numbers.csv")); bea <- fread(file.path(G, "main_with_ukbb_bmi_effect_by_age.csv"))
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
ctl <- fread(file.path(G, "main_with_ukbb_control_bmi_tests.csv"))
tests <- rbind(ctl[, .(Outcome = outcome, Test = "Controls: BMI effect varies with age (age x BMI tensor interaction)", Statistic = NA_real_, df = NA_integer_, P = signif(p_age_x_bmi, 2))], tests)
fwrite(t22, file.path(ST, "supp_table_S22_fig6_growth_tests_ukbb.csv")); fwrite(tests, file.path(ST, "supp_table_S22b_fig6_genotype_tests.csv"))

# ---- T17, T18 (as 32_supp_asd_by_cnv.R: the S7 / S8 tables)
t17 <- fread(file.path(A, "supp_table_S7_asd_by_cnv.csv")); t18 <- fread(file.path(A, "supp_table_S8_asd_by_cnv_pgs_adjusted.csv"))
fwrite(t17, file.path(ST, "supp_table_S17_asd_by_cnv.csv")); fwrite(t18, file.path(ST, "supp_table_S18_asd_by_cnv_pgs_adjusted.csv"))

# ---- workbook v3 (copy of v2; four sheets replaced)
wb <- loadWorkbook(wbf2)
put <- function(sh, title, note, tab, tab2 = NULL) {
  old <- if (sh %in% names(wb)) read.xlsx(wb, sh, colNames = FALSE, rows = 1:2, cols = 1)[[1]] else character(0)
  if (sh %in% names(wb)) removeWorksheet(wb, sh); addWorksheet(wb, sh)
  writeData(wb, sh, if (is.null(title)) old[1] else title, startRow = 1); writeData(wb, sh, paste0(if (is.null(note)) (if (length(old) > 1) old[2] else "") else note, SRC), startRow = 2)
  writeData(wb, sh, tab, startRow = 4, headerStyle = bold)
  if (!is.null(tab2)) writeData(wb, sh, tab2, startRow = nrow(tab) + 7, headerStyle = bold)
  setColWidths(wb, sh, 1:ncol(tab), "auto")
}
put("T17 ASD within CNV", "Supplementary Table 17. Height, head circumference and BMI by autism diagnosis within CNV group (Supplementary Fig. 4)",
    "GAM y ~ autism + s(age) + sex (+ cohort, + locus for pooled groups); int_*: autism x carrier interaction against SSC probands vs siblings.", t17)
put("T18 ASD within CNV PGS adj", "Supplementary Table 18. Autism contrasts within CNV group with polygenic scores as covariates (same individuals)",
    "Models as Supplementary Table 17 in the subset with polygenic scores, adding each score and all four.", t18)
put("T21 BMI slope by age band", "Supplementary Table 21. BMI effect on height and head circumference within age bands: carriers versus non-carrier controls",
    "Slopes are SD of the outcome per SD of BMI from one least-squares model per CNV and age band, outcome Z ~ CNV x BMI Z + sex + age, fitted to carriers and controls together: control slope = BMI term, carrier slope = BMI + CNV x BMI. Below: whether the BMI effect differs by genotype within each band (controls and all CNVs with at least 15 carriers in one model).",
    t21, rbind(t21b, data.table(Outcome = "", `Age band (years)` = ""), fill = TRUE))
# 0-2 y section below the genotype tests
wb_r <- nrow(t21) + 7 + nrow(t21b) + 3
writeData(wb, "T21 BMI slope by age band", "Infancy (0-2 years), carriers only: no non-carrier controls under 2 years, so within-carrier slopes (Z ~ BMI Z + sex + age; latest visit under 2 years per person; head circumference includes Simons VIP infant visits). Groups with at least 15 carriers.", startRow = wb_r)
writeData(wb, "T21 BMI slope by age band", t21c, startRow = wb_r + 1, headerStyle = bold)
put("T22 Fig6 growth tests UKBB", "Supplementary Table 22. Height and head circumference in CNV carriers versus non-carrier controls, with UK Biobank adults (main growth figure)",
    "Children: external growth references; adults: Z relative to the non-carrier adults of the same cohort (Box-Cox, lambda from UK Biobank). Shift: carrier - control mean, y ~ CNV + s(age), separately in children and adults. CNV x BMI: age-adjusted difference from controls in the BMI slope.",
    t22, tests)
ord <- names(wb); num <- suppressWarnings(as.numeric(sub("^T(\\d+).*", "\\1", ord))); num[is.na(num)] <- 999; num[grepl("^T19b", ord)] <- 19.5
worksheetOrder(wb) <- order(num); saveWorkbook(wb, wbf3, overwrite = TRUE)
cat("written:", wbf3, "\n T17", nrow(t17), "T18", nrow(t18), "T21", nrow(t21), "+", nrow(t21b), "T22", nrow(t22), "+", nrow(tests), "\n")
