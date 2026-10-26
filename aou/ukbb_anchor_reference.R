# UKBB adult anchor (runs on Expanse; individual data never leave). Aggregates only; cells with n < 5 suppressed.
suppressPackageStartupMessages(library(data.table))
g <- fread("/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv",
           select = c("f.eid", "sex_combined", "CNV_post_filter..1.yes..0.no.", "duplicate_drop", "withdrawn_consent", "ancestry",
                      "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D", "Height", "BMI_raw", "age"))
setnames(g, c("eid", "sex", "pf", "dup", "wd", "anc", "l16", "l22", "h", "bmi", "age"))
cat("raw rows", nrow(g), "\n"); print(table(pf = g$pf, useNA = "ifany")); print(table(dup = g$dup, useNA = "ifany")); print(table(wd = g$wd, useNA = "ifany")); print(table(anc = g$anc, useNA = "ifany"))
# CNV_post_filter is empty in this file; use duplicate_drop == "no", consent retained, and the largest ancestry group ("blue" = EUR)
g <- g[dup == "no" & wd == "no" & anc == "blue" & !is.na(h) & !is.na(bmi) & !is.na(age) & sex %in% c(1, 2)]
cat("after QC", nrow(g), "\n")
g[, locus := fcase(l16 == 1, "16p11.2 del", l16 == 3, "16p11.2 dup", l22 == 1, "22q11.2 del", l22 == 3, "22q11.2 dup",
                   l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
g <- g[!is.na(locus)]
print(g[, .N, by = locus])
# height and BMI Z: residuals on age + age^2 within sex, standardised, using non-carriers as reference
for (v in c("h", "bmi")) for (s in 1:2) {
  ref <- g[sex == s & locus == "No CNV"]; m <- lm(as.formula(paste(v, "~ age + I(age^2)")), data = ref); sdr <- sd(residuals(m))
  g[sex == s, paste0(v, "Z") := (get(v) - predict(m, newdata = .SD)) / sdr]
}
g[, bmi_grp := fifelse(bmiZ >= qnorm(.7), "high", fifelse(bmiZ <= qnorm(.3), "low", "mid"))]
g[, age_bin := cut(age, c(39, 49, 59, 74), labels = c("40-49", "50-59", "60-73"))]
agg <- function(by) { a <- g[, .(n = .N, mean_hZ = mean(hZ), se_hZ = sd(hZ) / sqrt(.N), mean_age = mean(age), mean_bmiZ = mean(bmiZ)), by = by]; a[n < 5, c("mean_hZ", "se_hZ", "mean_age", "mean_bmiZ") := NA]; a }
out <- rbind(agg(c("locus", "bmi_grp"))[, age_bin := "all"], agg(c("locus", "bmi_grp", "age_bin")), agg("locus")[, `:=`(bmi_grp = "all", age_bin = "all")], fill = TRUE)
fwrite(out[order(locus, bmi_grp, age_bin)], "ukbb_anchor_aggregates.csv")
# tests: within each locus and vs non-carriers (adult BMI-height relation)
tests <- rbindlist(lapply(c("16p11.2 del", "16p11.2 dup", "22q11.2 dup"), function(L) {   # 22q11.2 del: only 5 EUR carriers
  x <- g[locus %in% c(L, "No CNV")][, CNV := as.numeric(locus != "No CNV")]
  m1 <- summary(lm(hZ ~ bmiZ * CNV + factor(sex) + age, data = x))$coefficients
  xt <- x[bmi_grp != "mid"][, hi := as.numeric(bmi_grp == "high")]
  m2 <- summary(lm(hZ ~ hi * CNV + factor(sex) + age, data = xt))$coefficients
  mc <- summary(lm(hZ ~ bmiZ + factor(sex) + age, data = x[CNV == 1]))$coefficients
  data.table(locus = L, n_carriers = x[CNV == 1, .N], n_noncarriers = x[CNV == 0, .N],
             CNV_main = m1["CNV", 1], CNV_main_p = m1["CNV", 4],
             slope_noncarriers = m1["bmiZ", 1], slope_carriers = mc["bmiZ", 1], slope_carriers_p = mc["bmiZ", 4],
             CNVxBMI = m1["bmiZ:CNV", 1], CNVxBMI_se = m1["bmiZ:CNV", 2], CNVxBMI_p = m1["bmiZ:CNV", 4],
             tails_hi_minus_lo_noncarriers = m2["hi", 1], tails_CNVxhigh = m2["hi:CNV", 1], tails_CNVxhigh_p = m2["hi:CNV", 4])
}))
fwrite(tests, "ukbb_anchor_tests.csv"); print(tests)
