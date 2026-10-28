# Supplementary Table 21 (JS 2026-10-27): BMI effect on height and head circumference within age bands, carriers vs non-carrier
# controls, from the pooled band models of 17_main_with_ukbb.R (band_slopes; main_with_ukbb_bmi_slope_by_band.csv, copied from Expanse). Within
# each band: outcome Z ~ CNV * BMI Z + sex + age on carriers + controls; control slope = BMI, carrier slope = BMI + CNV x BMI, P = CNV x BMI.
# (Before 2026-10-27 the carrier and control slopes were separate fits from 36_check_bmi_effect_by_agebin.R.) Controls: clinical non-carriers plus 50,000 random UK Biobank non-carriers
# for height; clinical non-carriers for head circumference. Bands with < 15 carriers or < 30 controls are omitted.
# Output: supplementary_tables/supp_table_S21_bmi_slope_by_ageband.csv and sheet "T21 BMI slope by age band" in
# manuscript_materials_2026-10-24/Sacks Supplementary Tables v2.xlsx.
suppressMessages({ library(data.table); library(openxlsx) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
MM <- file.path(HERE, "manuscript_materials_2026-10-24")
x <- fread(file.path(HERE, "main_with_ukbb_bmi_slope_by_band.csv"))[group != "No CNV"]   # JS 2026-10-27: pooled-model slopes, the same numbers as the Fig. 6 prime-panel bars (17_main_with_ukbb.R, band_slopes)
x[, age_bin := factor(age_bin, c("2-9", "10-18", "19-39", "40-73"))][, outcome := factor(outcome, c("Height", "Head circumference"))]
setorder(x, outcome, group, age_bin)
t21 <- x[, .(Outcome = as.character(outcome), `CNV group` = group, `Age band (years)` = as.character(age_bin), `Carriers (n)` = n_carriers, `Controls (n)` = n_controls,
             `Carrier BMI slope` = round(slope_carriers, 3), `Carrier SE` = round(se_carriers, 3), `Control BMI slope` = round(slope_controls, 3), `Control SE` = round(se_controls, 3),
             `Carrier - control` = round(diff, 3), `Difference SE` = round(se_diff, 3), P = signif(p_diff, 2))]
fwrite(t21, file.path(MM, "supplementary_tables", "supp_table_S21_bmi_slope_by_ageband.csv"))
wbf <- file.path(MM, "Sacks Supplementary Tables v2.xlsx"); wb <- loadWorkbook(wbf); sh <- "T21 BMI slope by age band"
if (sh %in% names(wb)) removeWorksheet(wb, sh)
addWorksheet(wb, sh)
writeData(wb, sh, "Supplementary Table 21. BMI effect on height and head circumference within age bands: carriers versus non-carrier controls", startRow = 1)
writeData(wb, sh, "Slopes are SD of the outcome per SD of BMI from one least-squares model per CNV and age band, outcome Z ~ CNV x BMI Z + sex + age, fitted to carriers and controls together: control slope = BMI term, carrier slope = BMI + CNV x BMI terms; carrier - control and P = the CNV x BMI term (as plotted in Fig. 6 A'-J'). Height controls: clinical non-carriers plus 50,000 random UK Biobank non-carriers; head-circumference controls: clinical non-carriers. Bands with < 15 carriers or < 30 controls omitted.", startRow = 2)
writeData(wb, sh, t21, startRow = 4, headerStyle = createStyle(textDecoration = "bold")); setColWidths(wb, sh, 1:ncol(t21), "auto")
# keep the T* sheets in numeric order, notes last
ord <- names(wb); num <- suppressWarnings(as.numeric(sub("^T(\\d+).*", "\\1", ord))); num[is.na(num)] <- 999; num[grepl("^T19b", ord)] <- 19.5
worksheetOrder(wb) <- order(num)
saveWorkbook(wb, wbf, overwrite = TRUE)
cat("Supplementary Table 21 written:", nrow(t21), "rows; sheets:", paste(names(wb)[order(num)], collapse = " | "), "\n")
