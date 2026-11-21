# 83_growth_sensitivity_summary_v2.R  (JS / Claude 2026-11-20, BMI modifier paper v2)
#
# What it does: collects the aggregate outputs of 82_growth_blocks_v2.R (run on Expanse, one folder per variant, copied to
#   run_v2_2026-11-20/growth/<variant>/) into
#   1. reproduction_check_vs_v9e.csv: main variant vs the v9e run (band slopes, BMI effect by age, genotype-by-band tests,
#      fig6 numbers), maximum absolute difference per file;
#   2. sensitivity_summary.csv: band slopes (per genotype and pooled), all-ages CNV x BMI terms, genotype-by-band F tests and
#      the five-group omnibus, main next to each sensitivity variant;
#   3. supplement tables: supp_table_bmi_slope_by_ageband_v2.csv (update of Supplementary Table 21: per genotype + pooled),
#      supp_table_cnv_x_bmi_all_ages_v2.csv (all-ages terms with the five-group omnibus), supp_table_sex_tests_v2.csv.
# Env (2026-11-20): GROWTH_SUBDIR (growth | growth_unaffctl), MAIN_VARIANT (main_r3 | main), REPRO_VARIANT.
# Inputs: run_v2_2026-11-20/<GROWTH_SUBDIR>/{main_r3 (final layout; statistics identical to main and main_r2), sens_*}/*.csv; run_chop_fillgaps_2026-11-19/expanse_v9e/*.csv (reference).
# Outputs: run_v2_2026-11-20/growth/ (files above). RUNS LOCALLY (aggregate tables only). Date: 2026-11-20.
# null model: none; this script only reshapes results (tests are described in 17, 38 and 82).
suppressPackageStartupMessages(library(data.table))
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
G <- file.path(HERE, "run_v2_2026-11-20", Sys.getenv("GROWTH_SUBDIR", "growth")); REF <- file.path(HERE, "run_chop_fillgaps_2026-11-19/expanse_v9e")
VARS <- c(main = Sys.getenv("MAIN_VARIANT", "main_r3"), visit_latest = "sens_visit_latest", no_chop = "sens_no_chop", chop_overlap_keep = "sens_chop_overlap_keep",
          no_chop_g2mh203_kept = "sens_no_chop_g2mh203_kept", controls_ssc_probands = "sens_controls_ssc_probands")
REPRO <- Sys.getenv("REPRO_VARIANT", VARS[["main"]])   # the variant that must reproduce v9e (old control set): main_r3 in growth/, sens_controls_ssc_probands in growth_unaffctl/
VARS <- VARS[dir.exists(file.path(G, VARS))]
BL <- c("2-9", "10-18", "19-39", "40-73"); LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"); POOL <- "All four CNVs"

# ---- 1. reproduction of v9e
rc <- rbindlist(lapply(c("main_with_ukbb_bmi_slope_by_band.csv", "main_with_ukbb_bmi_effect_by_age.csv", "main_with_ukbb_bmi_genotype_by_band.csv",
                         "main_with_ukbb_height_tests.csv", "main_with_ukbb_hc_tests.csv", "main_with_ukbb_control_bmi_tests.csv", "main_with_ukbb_lifecourse.csv", "fig6_numbers.csv"), function(f) {
  a <- fread(file.path(REF, f)); b <- fread(file.path(G, REPRO, f))
  same_shape <- identical(dim(a), dim(b)) && identical(names(a), names(b))
  num <- names(a)[sapply(a, is.numeric)]
  md <- if (same_shape) max(sapply(num, function(v) suppressWarnings(max(abs(a[[v]] - b[[v]]), na.rm = TRUE))), -Inf) else NA_real_
  chr <- if (same_shape) all(sapply(setdiff(names(a), num), function(v) identical(as.character(a[[v]]), as.character(b[[v]])))) else FALSE
  na_same <- if (same_shape) all(sapply(num, function(v) identical(is.na(a[[v]]), is.na(b[[v]])))) else FALSE
  data.table(file = f, rows_v9e = nrow(a), rows_v10 = nrow(b), same_shape, labels_identical = chr, missing_pattern_identical = na_same, max_abs_diff = md,
             md5_identical = unname(tools::md5sum(file.path(REF, f)) == tools::md5sum(file.path(G, REPRO, f))))
}))
fwrite(rc, file.path(G, "reproduction_check_vs_v9e.csv")); print(rc)

# ---- 2. sensitivity summary
one <- function(v) { d <- file.path(G, VARS[[v]])
  bs <- fread(file.path(d, "main_with_ukbb_bmi_slope_by_band.csv"))[group != "No CNV", .(analysis = "BMI slope by age band, CNV x BMI (17 band_slopes)", outcome, group, age_bin, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff)]
  pb <- fread(file.path(d, "pooled_band_slopes.csv"))
  pbb <- pb[age_bin != "all ages" & grepl("primary", model), .(analysis = "pooled BMI slope by age band, CNV x BMI, locus-adjusted (82)", outcome, group, age_bin, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff)]
  pbh <- pb[age_bin != "all ages" & grepl("primary", model) & !is.na(het_p), .(analysis = "heterogeneity across CNVs within band, carriers only, F test (82)", outcome, group = "four CNVs", age_bin, n_carriers, p = het_p)]
  aa <- fread(file.path(d, "cnv_x_bmi_all_ages.csv"))
  a1 <- aa[grepl("est_main of 17", quantity), .(analysis = "CNV x BMI, all ages, age-adjusted (17 est_main)", outcome, group, age_bin = "all ages", n_carriers, n_controls, est, se, p)]
  a2 <- aa[grepl("pooled carriers", quantity) & grepl("primary", model), .(analysis = "pooled CNV x BMI, all ages, est_main form + locus (82)", outcome, group, age_bin = "all ages", n_carriers, n_controls, est, se, p)]
  a3 <- aa[grepl("omnibus", quantity) & grepl("differs", quantity), .(analysis = "five-group omnibus, all ages (38)", outcome, group, age_bin = "all ages", p)]
  a4 <- aa[grepl("heterogeneity", quantity), .(analysis = "heterogeneity across the four CNVs, all ages, carriers only (38)", outcome, group, age_bin = "all ages", p)]
  gb <- fread(file.path(d, "main_with_ukbb_bmi_genotype_by_band.csv"))[, .(analysis = "BMI effect differs by genotype within band, F test (17 bgt)", outcome, group = "Controls + CNVs", age_bin, p)]
  cbind(variant = v, rbind(bs, pbb, pbh, a1, a2, a3, a4, gb, fill = TRUE))
}
L <- rbindlist(lapply(names(VARS), one), fill = TRUE)
W <- dcast(L, analysis + outcome + group + age_bin ~ variant, value.var = c("n_carriers", "est", "se", "p"))
ord <- c("analysis", "outcome", "group", "age_bin", unlist(lapply(names(VARS), function(v) paste(c("n_carriers", "est", "se", "p"), v, sep = "_"))))
W <- W[, ..ord]
W[, `:=`(outcome = factor(outcome, c("Height", "Head circumference")), age_bin = factor(age_bin, c("all ages", BL)), analysis = factor(analysis, unique(L$analysis)))]
setorder(W, outcome, analysis, group, age_bin)
for (v in setdiff(names(VARS), "main")) {   # flags: nominal significance changed, or the estimate changed sign while nominally significant in either run
  pm <- W$p_main; pv <- W[[paste0("p_", v)]]; em <- W$est_main; ev <- W[[paste0("est_", v)]]
  W[, paste0("flag_", v) := fifelse(is.na(pm) | is.na(pv), fifelse(is.na(pm) & is.na(pv), "", "not estimable in one run"),
                                   fifelse((pm < 0.05) != (pv < 0.05), fifelse(pm < 0.05, "main p < 0.05, variant not", "variant p < 0.05, main not"), "")) ]
  sgn <- !is.na(em) & !is.na(ev) & sign(em) != sign(ev) & (pmin(pm, pv) < 0.05); W[sgn, paste0("flag_", v) := paste(get(paste0("flag_", v)), "sign differs")]
}
fwrite(W, file.path(G, "sensitivity_summary.csv"))
cat("\n== sensitivity: rows where nominal significance differs from main\n")
fl <- grep("^flag_", names(W), value = TRUE)
print(W[rowSums(as.matrix(W[, ..fl]) != "") > 0, c("analysis", "outcome", "group", "age_bin", grep("^(est|p|flag)_", names(W), value = TRUE)), with = FALSE][, analysis := substr(as.character(analysis), 1, 40)], nrows = 200)
cat("\n== sensitivity: rows nominally significant in main, all variants\n")
print(W[!is.na(p_main) & p_main < 0.05, c("analysis", "outcome", "group", "age_bin", grep("^(n_carriers|est|p)_", names(W), value = TRUE)), with = FALSE][, analysis := substr(as.character(analysis), 1, 40)], nrows = 200)

# ---- 3. supplement tables (main variant)
d <- file.path(G, VARS[["main"]])
bs <- fread(file.path(d, "main_with_ukbb_bmi_slope_by_band.csv"))[group != "No CNV"][, model := "outcome Z ~ CNV x BMI Z + sex + age"]
pb <- fread(file.path(d, "pooled_band_slopes.csv"))[age_bin != "all ages" & grepl("primary", model)][, model := "outcome Z ~ locus + BMI Z + CNV x BMI Z + sex + age"]
t21 <- rbind(bs, pb, fill = TRUE)[, `:=`(outcome = factor(outcome, c("Height", "Head circumference")), group = factor(group, c(LOCI, POOL)), age_bin = factor(age_bin, BL))]
setorder(t21, outcome, group, age_bin)
gb <- fread(file.path(d, "main_with_ukbb_bmi_genotype_by_band.csv"))
t21 <- merge(t21, gb[, .(outcome, age_bin, p_gb = p)], by = c("outcome", "age_bin"), all.x = TRUE, sort = FALSE)
t21o <- t21[, .(Outcome = as.character(outcome), `CNV group` = as.character(group), `Age band (years)` = as.character(age_bin), `Carriers (n)` = n_carriers, `Controls (n)` = n_controls,
                `Carrier BMI slope` = round(slope_carriers, 3), `Carrier SE` = round(se_carriers, 3), `Control BMI slope` = round(slope_controls, 3), `Control SE` = round(se_controls, 3),
                `Carrier - control (CNV x BMI)` = round(diff, 3), `Difference SE` = round(se_diff, 3), P = signif(p_diff, 2),
                `P heterogeneity across CNVs, carriers only (pooled rows)` = signif(het_p, 2), `P BMI effect differs by genotype, controls + CNVs (same for every row of the band)` = signif(p_gb, 2), Model = model)]
fwrite(t21o, file.path(G, "supp_table_bmi_slope_by_ageband_v2.csv"))
aa <- fread(file.path(d, "cnv_x_bmi_all_ages.csv"))
fwrite(aa[, .(Quantity = quantity, Outcome = outcome, Group = group, Model = model, `Carriers (n)` = n_carriers, `Controls (n)` = n_controls, Estimate = round(est, 4), SE = round(se, 4),
              `Chi-square` = round(chisq, 2), F = round(F, 2), df, P = signif(p, 3))], file.path(G, "supp_table_cnv_x_bmi_all_ages_v2.csv"))
sx <- fread(file.path(d, "sex_cnv_bmi_growth.csv")); r3 <- function(x) round(x, 3)
fwrite(sx[, .(Outcome = outcome, `CNV group` = group, `Age range` = scope, Method = method, `Male carriers (n)` = n_carriers_male, `Female carriers (n)` = n_carriers_female,
              `Male controls (n)` = n_controls_male, `Female controls (n)` = n_controls_female,
              `BMI slope, male carriers` = r3(slope_carriers_male), SE1 = r3(se_carriers_male), `BMI slope, female carriers` = r3(slope_carriers_female), SE2 = r3(se_carriers_female),
              `BMI slope, male controls` = r3(slope_controls_male), SE3 = r3(se_controls_male), `BMI slope, female controls` = r3(slope_controls_female), SE4 = r3(se_controls_female),
              `CNV x BMI, males` = r3(cnv_x_bmi_male), SE5 = r3(se_cnv_x_bmi_male), `P, males` = signif(p_cnv_x_bmi_male, 2),
              `CNV x BMI, females` = r3(cnv_x_bmi_female), SE6 = r3(se_cnv_x_bmi_female), `P, females` = signif(p_cnv_x_bmi_female, 2),
              `CNV x BMI x sex (male - female)` = r3(cnv_x_bmi_x_sex), SE7 = r3(se_cnv_x_bmi_x_sex), `P, three-way` = signif(p_cnv_x_bmi_x_sex, 2))], file.path(G, "supp_table_sex_tests_v2.csv"))
cat("\n== sex three-way terms, main (", nrow(sx), "tests; nominal p < 0.05:", sx[p_cnv_x_bmi_x_sex < 0.05, .N], "; Bonferroni threshold", signif(0.05 / nrow(sx), 2), ")\n")
print(sx[, .(outcome = substr(outcome, 1, 6), group, scope, n_m = n_carriers_male, n_f = n_carriers_female, cxb_m = round(cnv_x_bmi_male, 2), p_m = signif(p_cnv_x_bmi_male, 2),
             cxb_f = round(cnv_x_bmi_female, 2), p_f = signif(p_cnv_x_bmi_female, 2), three = round(cnv_x_bmi_x_sex, 2), se = round(se_cnv_x_bmi_x_sex, 2), p3 = signif(p_cnv_x_bmi_x_sex, 2))], nrows = 100)
cat("83 done\n")
