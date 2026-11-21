# 85_old_vs_corrected_controls_v2.R  (JS / Claude 2026-11-20)
# What it does: side-by-side table of the growth results with the OLD control set (paper run 2026-11-19: SSC autistic probands as
#   the child controls; variant sens_controls_ssc_probands, identical to v9e / v10) and the CORRECTED control set (00_control_set.R:
#   SSC designated siblings + SSC parents + Simons VIP non-carriers without autism; variant main), both run by
#   82_growth_blocks_v2.R on Expanse in v11_unaffected_controls_2026-11-20.
#   Sections: control BMI slopes by age band; each CNV's CNV x BMI by band and all ages; five-group omnibus; pooled "All four
#   CNVs" by band and all ages; genotype-by-band F tests; heterogeneity; sex (CNV x BMI x sex) tests; CNV main effects by band (43).
# Inputs: run_v2_2026-11-20/growth_unaffctl/{main, sens_controls_ssc_probands}/*.csv (aggregates copied from Expanse).
# Output: run_v2_2026-11-20/growth_unaffctl/old_vs_corrected_controls.csv. RUNS LOCALLY.
# null model: none; this script only lines up results (tests are described in 17, 38, 43 and 82).
suppressPackageStartupMessages(library(data.table)); options(width = 230)
G <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/run_v2_2026-11-20/growth_unaffctl"
BL <- c("2-9", "10-18", "19-39", "40-73")
get1 <- function(v) { d <- file.path(G, v); rd <- function(f) fread(file.path(d, f))
  bs <- rd("main_with_ukbb_bmi_slope_by_band.csv"); pb <- rd("pooled_band_slopes.csv"); aa <- rd("cnv_x_bmi_all_ages.csv"); gb <- rd("main_with_ukbb_bmi_genotype_by_band.csv")
  sx <- rd("sex_cnv_bmi_growth.csv"); cl <- rd("band_slopes_family_clustered.csv"); me <- if (file.exists(file.path(d, "main_effects_by_band.csv"))) rd("main_effects_by_band.csv") else NULL
  k <- c("outcome", "group", "age_bin")
  bsc <- merge(bs[group != "No CNV"], cl[, .(outcome, group, age_bin, p_cl = p_diff_family_clustered)], by = k, all.x = TRUE)
  pbc <- pb[age_bin != "all ages" & grepl("primary", model)]
  rbind(
    bs[group == "No CNV", .(section = "1 control BMI slope by band (y ~ BMI + sex + age in controls)", outcome, group = "Non-carrier controls", age_bin, n_carriers = NA_integer_, n_controls, est = slope_controls, se = se_controls, p = p_controls)],
    bsc[, .(section = "2 CNV x BMI by band, per genotype (17 band_slopes)", outcome, group, age_bin, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff, p_family_clustered = p_cl)],
    aa[grepl("est_main of 17", quantity), .(section = "3 CNV x BMI, all ages, per genotype (17 est_main)", outcome, group, age_bin = "all ages", n_carriers, n_controls, est, se, p)],
    aa[grepl("omnibus", quantity) & grepl("differs", quantity), .(section = "4 five-group omnibus, all ages (38)", outcome, group, age_bin = "all ages", est = chisq, p)],
    aa[grepl("heterogeneity", quantity), .(section = "4b heterogeneity across the four CNVs, carriers only, all ages (38)", outcome, group, age_bin = "all ages", est = F, p)],
    pbc[, .(section = "5 pooled CNV x BMI by band, locus-adjusted (82)", outcome, group, age_bin, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff, p_family_clustered = p_diff_family_clustered)],
    pb[age_bin != "all ages" & !grepl("primary", model), .(section = "5b pooled CNV x BMI by band, no locus term (as 35f)", outcome, group, age_bin, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff)],
    aa[grepl("pooled carriers", quantity) & grepl("primary", model), .(section = "5c pooled CNV x BMI, all ages, est_main form + locus (82)", outcome, group, age_bin = "all ages", n_carriers, n_controls, est, se, p)],
    gb[, .(section = "6 BMI effect differs by genotype within band, F test (17 bgt)", outcome, group = "Controls + CNVs", age_bin, est = F, p)],
    sx[, .(section = "7 CNV x BMI x sex (male - female) (82)", outcome, group, age_bin = scope, n_carriers = n_carriers_male + n_carriers_female, n_controls = n_controls_male + n_controls_female,
           est = cnv_x_bmi_x_sex, se = se_cnv_x_bmi_x_sex, p = p_cnv_x_bmi_x_sex, p_family_clustered = p_cnv_x_bmi_x_sex_family_clustered)],
    sx[, .(section = "7b CNV x BMI in males (82)", outcome, group, age_bin = scope, n_carriers = n_carriers_male, n_controls = n_controls_male, est = cnv_x_bmi_male, se = se_cnv_x_bmi_male, p = p_cnv_x_bmi_male)],
    sx[, .(section = "7c CNV x BMI in females (82)", outcome, group, age_bin = scope, n_carriers = n_carriers_female, n_controls = n_controls_female, est = cnv_x_bmi_female, se = se_cnv_x_bmi_female, p = p_cnv_x_bmi_female)],
    sx[, .(section = "7d control BMI slope, males", outcome, group, age_bin = scope, n_controls = n_controls_male, est = slope_controls_male, se = se_controls_male)],
    sx[, .(section = "7e control BMI slope, females", outcome, group, age_bin = scope, n_controls = n_controls_female, est = slope_controls_female, se = se_controls_female)],
    if (!is.null(me)) me[group != "Non-carriers", .(section = "8 CNV main effect by band: carrier - control shift (43, Supplementary Fig. 6)", outcome, group, age_bin = band, n_carriers = n, est = shift_vs_controls, se = shift_se, p = shift_p)],
    if (!is.null(me)) me[group == "Non-carriers", .(section = "8b control mean Z by band (43)", outcome, group = "Non-carrier controls", age_bin = band, n_controls = n, est = mean_z, se)],
    fill = TRUE)
}
o <- get1("sens_controls_ssc_probands"); n <- get1("main"); k <- c("section", "outcome", "group", "age_bin")
W <- merge(o, n, by = k, all = TRUE, suffixes = c("_old", "_corrected"))
W[, n_carriers := fcoalesce(n_carriers_corrected, n_carriers_old)][, c("n_carriers_old", "n_carriers_corrected", "p_family_clustered_old") := NULL]
W[, change := fifelse(is.na(p_old) | is.na(p_corrected), fifelse(xor(is.na(est_old), is.na(est_corrected)), "estimable in one set only", ""),
                      fifelse((p_old < 0.05) & !(p_corrected < 0.05), "LOST: p < 0.05 with old controls only", fifelse(!(p_old < 0.05) & (p_corrected < 0.05), "NEW: p < 0.05 with corrected controls only", "")))]
W[grepl("^(2|3|5|7|8 )", section) & !is.na(est_old) & !is.na(est_corrected) & sign(est_old) != sign(est_corrected) & pmin(p_old, p_corrected, na.rm = TRUE) < 0.05, change := trimws(paste(change, "sign differs"))]
W[, `:=`(outcome = factor(outcome, c("Height", "Head circumference", "BMI")), age_bin = factor(age_bin, c("all ages", BL)))]
setcolorder(W, c(k, "n_carriers", "n_controls_old", "est_old", "se_old", "p_old", "n_controls_corrected", "est_corrected", "se_corrected", "p_corrected", "p_family_clustered_corrected", "change"))
setorder(W, section, outcome, group, age_bin); fwrite(W, file.path(G, "old_vs_corrected_controls.csv"))
sh <- function(x) { x <- copy(x); for (v in names(x)) if (is.numeric(x[[v]]) && !is.integer(x[[v]])) set(x, j = v, value = signif(x[[v]], 2)); x[, section := substr(section, 1, 22)]
  setnames(x, c("n_controls_old", "n_controls_corrected", "est_corrected", "se_corrected", "p_corrected", "p_family_clustered_corrected"), c("nctl_old", "nctl_new", "est_new", "se_new", "p_new", "p_new_cl")); x }
for (s in unique(W$section)) if (!grepl("^7[b-e]|^5b|^8b", s)) { cat("\n==", s, "\n"); print(sh(W[section == s])[, section := NULL], nrows = 200) }
cat("\n== rows that change nominal significance or sign:", W[change != "", .N], "of", W[!is.na(p_old) | !is.na(p_corrected), .N], "tests\n")
print(sh(W[change != ""]), nrows = 200)
