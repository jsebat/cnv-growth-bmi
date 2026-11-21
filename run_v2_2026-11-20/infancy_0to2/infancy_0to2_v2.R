# infancy_0to2_v2.R  (2026-11-20, BMI modifier paper v2, supplementary tables package)
#
# What it does: reruns the 0-2 year, carriers-only block of Supplementary Table 21 (BMI effect on head circumference and
#   height in infancy). The v3 version (run_chop_fillgaps_2026-11-19/preview_0to2/preview_0to2.R) read CHOP infants from
#   analysis_table_chop, in which rule 2 of 00_qc_growth.R (height outside 60-230 cm -> height AND BMI set to missing) had
#   removed every valid infant length under 60 cm, and with it the BMI of the CHOP visits before about 4 months.
#   Fix (as in 80_longitudinal_v2.R): CHOP visits under 2 y are read from the workbook through chop_rows() of 00_add_chop.R.
#   The other plausibility rules of 00_qc_growth.R are applied to those rows here (height over 230 cm; BMI outside 10-70;
#   |height Z| > 5; |BMI Z| > 5; head circumference outside 30-70 cm or |Z| > 5); only the 60 cm lower bound is not applied.
#   Everything else is as in the v3 script: one record per person (latest visit under 2 y with the outcome and BMI),
#   scoring by 00_scoring.R (UK90 head circumference, WHO 2006 length and BMI-for-age), |Z| <= 5, groups with >= 15 carriers,
#   head circumference adds the Simons VIP infant visits of hc_visits.rds. 00_qc_growth.R itself is not changed.
# Model (per genotype): Z ~ BMI Z + sex + age, least squares, carriers only (there are no non-carrier controls under 2 y).
#   New in v2 (labelled "new"): all four CNVs pooled with a genotype term, and the heterogeneity test BMI x genotype.
# null model: none; t test of the BMI coefficient within genotype, F test of BMI x genotype for heterogeneity.
# Inputs (participant level, controlled): DATA/analysis_table_chop_unaffctl_2026-11-20.rds (carrier rows are the same as in
#   analysis_table_chop_2026-11-19.rds), DATA/hc_visits.rds, the CHOP workbook (00_add_chop.R).
# Outputs (HERE/run_v2_2026-11-20/infancy_0to2/; aggregates only): infancy_0to2_v2_bands.csv (v2 result),
#   infancy_0to2_old_vs_new.csv (v3 file, v3 reproduced here from the table, v2), infancy_0to2_n_by_cohort.csv,
#   infancy_0to2_qc_counts.csv (visits lost to the 60 cm rule, by cohort), infancy_0to2_v2.log.
suppressPackageStartupMessages({ library(data.table) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
OUT <- file.path(HERE, "run_v2_2026-11-20", "infancy_0to2"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(HERE); source("00_scoring.R"); source("00_add_chop.R")
sink(file.path(OUT, "infancy_0to2_v2.log"), split = TRUE); options(width = 220)
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"); ALL4 <- "All four CNVs"

tab <- as.data.table(readRDS(file.path(DATA, "analysis_table_chop_unaffctl_2026-11-20.rds")))
tab[, sex := as.character(sex)][, locus := as.character(locus)]
cl <- tab[!is.na(sex) & age >= 0 & age < 2, .(id, cohort, locus, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw)]

# CHOP visits under 2 y from the workbook, with the plausibility rules of 00_qc_growth.R except the 60 cm lower bound
cw <- as.data.table(chop_rows())[age >= 0 & age < 2]
cw <- cw[, .(id = paste0("CHOP:", MSassigned_ID), cohort = "CHOP", locus = ifelse(CNV == "22q11.2 A-D DEL", "22q11.2 del", "22q11.2 dup"),
             sex = ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA_character_)), age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw)][!is.na(sex)]
n_lt60 <- cw[!is.na(HEIGHT_raw) & HEIGHT_raw < 60, .N]; n_lt60_bmi <- cw[!is.na(HEIGHT_raw) & HEIGHT_raw < 60 & !is.na(BMI_raw), .N]
n_lt60_hc <- cw[!is.na(HEIGHT_raw) & HEIGHT_raw < 60 & !is.na(BMI_raw) & !is.na(HEADCIRCUM_raw), .N]
cw[!is.na(HEIGHT_raw) & HEIGHT_raw > 230, c("HEIGHT_raw", "BMI_raw") := .(NA_real_, NA_real_)]
cw[!is.na(BMI_raw) & (BMI_raw < 10 | BMI_raw > 70), BMI_raw := NA_real_]
cw[, hz := score_height(HEIGHT_raw, age, sex)][!is.na(hz) & abs(hz) > 5, HEIGHT_raw := NA_real_]
cw[, bz := score_bmi_child(BMI_raw, age, sex)][!is.na(bz) & abs(bz) > 5, BMI_raw := NA_real_]
cw[, cz := score_hc(HEADCIRCUM_raw, age, sex)][!is.na(HEADCIRCUM_raw) & (HEADCIRCUM_raw < 30 | HEADCIRCUM_raw > 70 | (!is.na(cz) & abs(cz) > 5)), HEADCIRCUM_raw := NA_real_]
cw[, c("hz", "bz", "cz") := NULL]
tb <- cl[cohort == "CHOP"]
qc <- data.table(cohort = "CHOP", quantity = c("workbook visits under 2 y with a length under 60 cm", "... of these with a BMI", "... of these with a BMI and a head circumference",
                                              "visits under 2 y with height and BMI: analysis table", "visits under 2 y with height and BMI: workbook (v2)",
                                              "visits under 2 y with head circumference and BMI: analysis table", "visits under 2 y with head circumference and BMI: workbook (v2)",
                                              "carriers with head circumference and BMI under 2 y: analysis table", "carriers with head circumference and BMI under 2 y: workbook (v2)",
                                              "carriers with height and BMI under 2 y: analysis table", "carriers with height and BMI under 2 y: workbook (v2)"),
                 n = c(n_lt60, n_lt60_bmi, n_lt60_hc, tb[!is.na(HEIGHT_raw) & !is.na(BMI_raw), .N], cw[!is.na(HEIGHT_raw) & !is.na(BMI_raw), .N],
                       tb[!is.na(HEADCIRCUM_raw) & !is.na(BMI_raw), .N], cw[!is.na(HEADCIRCUM_raw) & !is.na(BMI_raw), .N],
                       tb[!is.na(HEADCIRCUM_raw) & !is.na(BMI_raw), uniqueN(id)], cw[!is.na(HEADCIRCUM_raw) & !is.na(BMI_raw), uniqueN(id)],
                       tb[!is.na(HEIGHT_raw) & !is.na(BMI_raw), uniqueN(id)], cw[!is.na(HEIGHT_raw) & !is.na(BMI_raw), uniqueN(id)]))
# the other cohorts: how many infant lengths under 60 cm did the rule remove? (pre-QC copy of the pre-CHOP table, 2026-10-24)
pq <- file.path(DATA, "analysis_table.preQC_20261024.rds")
if (file.exists(pq)) { p0 <- as.data.table(readRDS(pq)); if (all(c("HEIGHT_raw", "age", "cohort", "locus") %in% names(p0))) {
  o <- p0[as.character(locus) %in% LOCI & age < 2 & !is.na(HEIGHT_raw) & HEIGHT_raw >= 3 & HEIGHT_raw < 60, .(n = .N), by = cohort]
  if (nrow(o)) qc <- rbind(qc, o[, .(cohort, quantity = "carrier visits under 2 y with a length under 60 cm in the pre-QC table (lost to rule 2; not restored here)", n)])
  else qc <- rbind(qc, data.table(cohort = "G2MH, Cardiff, Simons VIP, Searchlight, SSC", quantity = "carrier visits under 2 y with a length under 60 cm in the pre-QC table", n = 0L)) } }
fwrite(qc, file.path(OUT, "infancy_0to2_qc_counts.csv")); cat("== QC counts\n"); print(qc)

hv <- as.data.table(readRDS(file.path(DATA, "hc_visits.rds")))[age >= 0 & age < 2 & !is.na(hc) & !is.na(bmi)]
hv <- hv[order(id, -age)][, .SD[1], by = id][, .(id = paste0("hcv:", id), cohort = paste0(cohort, " (hc_visits)"), locus, sex = as.character(sex), age, HEADCIRCUM_raw = hc, BMI_raw = bmi)]

run <- function(src, label) {
  res <- list(); cnt <- list()
  for (o in c("Head circumference", "Height")) {
    raw <- if (o == "Height") "HEIGHT_raw" else "HEADCIRCUM_raw"
    x <- src[!is.na(get(raw)) & !is.na(BMI_raw)][order(id, -age)][, .SD[1], by = id]
    if (o == "Head circumference") x <- rbind(x, hv, fill = TRUE)
    x[, y := if (o == "Height") score_height(get(raw), age, sex) else score_hc(get(raw), age, sex)][, BMI := score_bmi_child(BMI_raw, age, sex)]
    x <- x[!is.na(y) & !is.na(BMI) & abs(y) <= 5 & abs(BMI) <= 5 & locus %in% LOCI]
    cnt[[o]] <- x[, .N, by = .(locus, cohort)][, outcome := o]
    for (G in LOCI) { d <- x[locus == G]; if (nrow(d) < 15) next
      cf <- summary(lm(y ~ BMI + sex + age, data = d))$coefficients["BMI", ]
      res[[length(res) + 1]] <- data.table(outcome = o, group = G, age_bin = "0-2", model = "Z ~ BMI Z + sex + age", n = nrow(d), slope = cf[1], se = cf[2], p = cf[4],
                                           median_age = median(d$age), n_under_4_months = sum(d$age < 1 / 3), het_F = NA_real_, het_df = NA_real_, het_p = NA_real_) }
    d <- x[locus %in% x[, .N, by = locus][N >= 15, locus]]   # pooled over the genotypes with >= 15 carriers (new in v2)
    m1 <- lm(y ~ BMI + sex + age + locus, data = d); m2 <- lm(y ~ BMI * locus + sex + age, data = d); a <- anova(m1, m2); cf <- summary(m1)$coefficients["BMI", ]
    res[[length(res) + 1]] <- data.table(outcome = o, group = ALL4, age_bin = "0-2", model = paste0("new: Z ~ BMI Z + sex + age + CNV (", uniqueN(d$locus), " CNVs with >= 15 carriers); heterogeneity = F test of BMI x CNV"),
                                         n = nrow(d), slope = cf[1], se = cf[2], p = cf[4], median_age = median(d$age), n_under_4_months = sum(d$age < 1 / 3), het_F = a$F[2], het_df = a$Df[2], het_p = a$`Pr(>F)`[2])
  }
  list(res = rbindlist(res)[, version := label], cnt = rbindlist(cnt)[, version := label])
}
old <- run(cl, "v3 reproduced: CHOP from the analysis table (infant lengths under 60 cm lost)")
new <- run(rbind(cl[cohort != "CHOP"], cw), "v2: CHOP from the workbook (infant lengths kept)")
fwrite(new$res[, !"version"], file.path(OUT, "infancy_0to2_v2_bands.csv"))
cnt <- rbind(old$cnt, new$cnt); cnt[, N := as.character(N)][as.integer(N) < 5, N := "<5"]   # cohort-level counts: cells below 5 masked
fwrite(dcast(cnt, version + outcome + locus ~ cohort, value.var = "N", fill = "0"), file.path(OUT, "infancy_0to2_n_by_cohort.csv"))
v3 <- fread(file.path(HERE, "run_chop_fillgaps_2026-11-19", "preview_0to2", "preview_0to2_bands.csv"))
cmp <- merge(merge(v3[, .(outcome, group, n_v3_file = n, slope_v3_file = slope, se_v3_file = se, p_v3_file = p)],
                   old$res[, .(outcome, group, n_v3_reproduced = n, slope_v3_reproduced = slope, p_v3_reproduced = p)], by = c("outcome", "group"), all = TRUE),
             new$res[, .(outcome, group, n_v2 = n, slope_v2 = slope, se_v2 = se, p_v2 = p, median_age_v2 = median_age, n_under_4_months_v2 = n_under_4_months, het_p_v2 = het_p)], by = c("outcome", "group"), all = TRUE)
cmp[, reproduced := ifelse(is.na(slope_v3_file), NA, abs(slope_v3_file - slope_v3_reproduced) < 1e-8 & n_v3_file == n_v3_reproduced)]
fwrite(cmp, file.path(OUT, "infancy_0to2_old_vs_new.csv"))
cat("\n== sample by cohort (masked)\n"); print(dcast(cnt, version + outcome + locus ~ cohort, value.var = "N", fill = "0"))
cat("\n== old vs new\n"); print(cmp[, .(outcome, group, n_v3_file, slope_v3 = round(slope_v3_file, 3), p_v3 = signif(p_v3_file, 2), reproduced, n_v2, slope_v2 = round(slope_v2, 3), se_v2 = round(se_v2, 3), p_v2 = signif(p_v2, 2),
                                         median_age_v2 = round(median_age_v2, 2), n_under_4_months_v2, het_p_v2 = signif(het_p_v2, 2))])
sink()
