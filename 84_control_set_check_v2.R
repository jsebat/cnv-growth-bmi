# 84_control_set_check_v2.R  (JS / Claude 2026-11-20)
# What it does: documents the non-carrier control set of the growth analysis before and after the correction of 00_control_set.R.
#   1. Rebuild check: 01_build.R with CONTROL_SET=paper_2026-11-19 must reproduce analysis_table_chop_2026-11-19.rds (run
#      separately; here the rebuilt table, if present, is compared column by column).
#   2. Clinical controls ("No CNV") eligible for the height and the head-circumference models (outcome + BMI measured, age 2-73,
#      sex known; one person counted once), by role, sex and age band (2-9, 10-18, 19-39, 40-73 at the latest eligible visit),
#      old table vs corrected table. Roles come from control_role of the corrected table (same person ids in both tables).
#   3. Carriers: confirms the carrier rows are identical in the two tables.
# Inputs: Jon/Claude/CNV_growth_trajectories_data/analysis_table_chop_2026-11-19.rds, analysis_table_chop_unaffctl_2026-11-20.rds,
#   (optional) analysis_table_chop_oldctl_rebuild.rds. Output: run_v2_2026-11-20/growth_unaffctl/control_n_by_role_sex_band.csv.
# RUNS LOCALLY. null model: none; counts only.
suppressPackageStartupMessages(library(data.table)); options(width = 200)
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")
OUT <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/run_v2_2026-11-20/growth_unaffctl")
old <- readRDS(file.path(DATA, "analysis_table_chop_2026-11-19.rds")); new <- readRDS(file.path(DATA, "analysis_table_chop_unaffctl_2026-11-20.rds"))
rb <- file.path(DATA, "analysis_table_chop_oldctl_rebuild.rds")
if (file.exists(rb)) { r <- readRDS(rb); cm <- names(old); setkey(r, id, age); o2 <- copy(old); setkey(o2, id, age)
  cat("rebuild with CONTROL_SET=paper_2026-11-19 identical to the 2026-11-19 table on its columns:", isTRUE(all.equal(as.data.frame(o2[, ..cm]), as.data.frame(r[, ..cm]), check.attributes = FALSE)), "\n") }
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
a <- old[locus %in% LOCI][order(id, age)]; b <- new[locus %in% LOCI][order(id, age)]; cc <- intersect(names(a), names(b))
cat("carrier rows identical in both tables:", isTRUE(all.equal(as.data.frame(a[, ..cc]), as.data.frame(b[, ..cc]), check.attributes = FALSE)), " (", nrow(a), "rows,", uniqueN(a$id), "carriers )\n")
role <- unique(new[, .(id, control_role, family)])
BANDS <- c(2, 10, 19, 40, 74); BLAB <- c("2-9", "10-18", "19-39", "40-73")
tab <- function(t, lab) rbindlist(lapply(c("HEIGHT_raw", "HEADCIRCUM_raw"), function(raw) {
  x <- t[locus == "No CNV" & !is.na(sex) & age >= 2 & age <= 73 & !is.na(get(raw)) & !is.na(BMI_raw)][order(id, -age)][, .SD[1], by = id]
  x <- merge(x[, .(id, cohort, sex, age)], role, by = "id", all.x = TRUE)[is.na(control_role), control_role := paste(cohort, "non-carrier")]
  x[, band := as.character(cut(age, BANDS, labels = BLAB, right = FALSE))]
  y <- rbind(x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female"), median_age = round(median(age), 1)), by = .(role = control_role, band)],
             x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female"), median_age = round(median(age), 1)), by = .(role = control_role)][, band := "all ages"],
             x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female"), median_age = round(median(age), 1)), by = .(band)][, role := "ALL CLINICAL CONTROLS"],
             x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female"), median_age = round(median(age), 1))][, `:=`(role = "ALL CLINICAL CONTROLS", band = "all ages")], use.names = TRUE)
  cbind(control_set = lab, outcome = if (raw == "HEIGHT_raw") "Height" else "Head circumference", y, n_families = NA_integer_)[role == "ALL CLINICAL CONTROLS" & band == "all ages", n_families := uniqueN(x$family)]
}))
res <- rbind(tab(old, "old (paper run 2026-11-19: SSC probands as child controls)"), tab(new, "corrected (unaffected: SSC designated siblings + parents, SVIP without autism)"))
res[, `:=`(band = factor(band, c("all ages", BLAB)), pct_male = round(100 * n_male / n))]; setorder(res, outcome, control_set, role, band)
fwrite(res, file.path(OUT, "control_n_by_role_sex_band.csv"))
print(dcast(res, outcome + role + band ~ substr(control_set, 1, 9), value.var = c("n", "pct_male"))[order(outcome, role, band)], nrows = 200)
print(res[role == "ALL CLINICAL CONTROLS" & band == "all ages"])
