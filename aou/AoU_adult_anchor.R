### All of Us: adult height by BMI in 16p11.2 / 22q11.2 CNV carriers ("adult anchor") ###
# For Molly to run in the All of Us Researcher Workbench (Jon Sebat, 2026-10-13).
#
# PURPOSE
#   Adult counterpart to the clinical growth-trajectory figure (16p11.2 BP4-BP5 and 22q11.2 A-D carriers): mean height Z
#   of adult carriers by BMI tail (>= 70th vs <= 30th centile) and age band, and the adult CNV x BMI test on height
#   (are high-BMI adult carriers shorter than low-BMI carriers, beyond the BMI-height relation in non-carriers?).
#   Identical in design to the UK Biobank version (ukbb_anchor.R), so results can be pooled.
#
# INPUT  - the same QC'd genomatrix used for the CNV main-effect and CNV x PRS analyses, one per ancestry, e.g.
#          genomatrix/genomatrix_CNV_EUR.20250715.csv (also AAM, LATNAT). Columns used:
#          person_id, sex_at_birth, Height_raw, Height_age, BMI_raw, BMI_age, X16p11.2_proximal, X22q11.2_BP1_BP4_A_D
#          (locus coding: 1 = deletion, 2 = reference, 3 = duplication).
# OUTPUT - two small CSVs of AGGREGATE statistics only, in ./adult_anchor_output/:
#          AoU_<ANC>_adult_anchor_aggregates.csv, AoU_<ANC>_adult_anchor_tests.csv
#          Every group with fewer than 20 participants is suppressed (count shown as "<20", statistics blank), per the
#          All of Us Data and Statistics Dissemination Policy. Please check before exporting.
#
# HOW TO RUN  - edit the two lines under "SETTINGS", then source this file in an R notebook/terminal. Base R only.

## ---------------- SETTINGS ----------------
ANCESTRY  <- "EUR"                                                   # EUR, AAM or LATNAT
GENOMATRIX <- "../../../genomatrix/genomatrix_CNV_EUR.20250715.csv"  # path to that ancestry's genomatrix
MIN_CELL  <- 20                                                      # AoU small-cell threshold
## -------------------------------------------

dir.create("adult_anchor_output", showWarnings = FALSE)
cols <- c("person_id", "sex_at_birth", "Height_raw", "Height_age", "BMI_raw", "BMI_age", "X16p11.2_proximal", "X22q11.2_BP1_BP4_A_D")
g <- read.csv(GENOMATRIX)[, cols]
names(g) <- c("id", "sex", "h", "h_age", "bmi", "bmi_age", "l16", "l22")
g$sex <- ifelse(g$sex == "Female", 2, ifelse(g$sex == "Male", 1, NA))
g$age <- g$h_age
g <- g[!is.na(g$h) & !is.na(g$bmi) & !is.na(g$age) & !is.na(g$sex) & !is.na(g$l16) & !is.na(g$l22), ]
g$locus <- ifelse(g$l16 == 1, "16p11.2 del", ifelse(g$l16 == 3, "16p11.2 dup", ifelse(g$l22 == 1, "22q11.2 del",
           ifelse(g$l22 == 3, "22q11.2 dup", ifelse(g$l16 == 2 & g$l22 == 2, "No CNV", NA)))))
g <- g[!is.na(g$locus), ]
cat("participants after QC:", nrow(g), "\n"); print(table(g$locus))

# Height and BMI Z: residuals on age + age^2 within sex, scaled by the residual SD of NON-CARRIERS (same as UKBB version)
for (v in c("h", "bmi")) {
  z <- rep(NA_real_, nrow(g))
  for (s in 1:2) {
    ref <- g[g$sex == s & g$locus == "No CNV", ]
    m <- lm(as.formula(paste(v, "~ age + I(age^2)")), data = ref)
    idx <- g$sex == s
    z[idx] <- (g[[v]][idx] - predict(m, newdata = g[idx, ])) / sd(residuals(m))
  }
  g[[paste0(v, "Z")]] <- z
}
g$bmi_grp <- ifelse(g$bmiZ >= qnorm(0.7), "high", ifelse(g$bmiZ <= qnorm(0.3), "low", "mid"))
g$age_bin <- as.character(cut(g$age, c(17, 39, 49, 59, 69, 120), labels = c("18-39", "40-49", "50-59", "60-69", "70+")))

## aggregates: n, mean height Z, SE, mean age, mean BMI Z per locus x BMI group (x age band)
agg <- function(by) {
  key <- interaction(g[by], drop = TRUE, sep = "|")
  out <- do.call(rbind, lapply(split(seq_len(nrow(g)), key), function(i) {
    x <- g[i, ]; n <- nrow(x)
    data.frame(x[1, by, drop = FALSE], n = n, mean_hZ = mean(x$hZ), se_hZ = sd(x$hZ) / sqrt(n),
               mean_age = mean(x$age), mean_bmiZ = mean(x$bmiZ), row.names = NULL)
  }))
  small <- out$n < MIN_CELL
  out[small, c("mean_hZ", "se_hZ", "mean_age", "mean_bmiZ")] <- NA
  out$n <- ifelse(small, paste0("<", MIN_CELL), as.character(out$n))
  out
}
a1 <- agg(c("locus", "bmi_grp")); a1$age_bin <- "all"
a2 <- agg(c("locus", "bmi_grp", "age_bin"))
a3 <- agg("locus"); a3$bmi_grp <- "all"; a3$age_bin <- "all"
aggs <- rbind(a1[, names(a2)], a2, a3[, names(a2)])
aggs$cohort <- paste0("AoU.", ANCESTRY)
write.csv(aggs, file.path("adult_anchor_output", paste0("AoU_", ANCESTRY, "_adult_anchor_aggregates.csv")), row.names = FALSE)

## tests: adult BMI-height relation in carriers vs non-carriers (continuous BMI Z and BMI tails)
tests <- do.call(rbind, lapply(c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"), function(L) {
  x <- g[g$locus %in% c(L, "No CNV"), ]; x$CNV <- as.numeric(x$locus != "No CNV")
  nc <- sum(x$CNV == 1)
  res <- data.frame(cohort = paste0("AoU.", ANCESTRY), locus = L, n_carriers = ifelse(nc < MIN_CELL, paste0("<", MIN_CELL), nc),
                    n_noncarriers = sum(x$CNV == 0))
  if (nc < MIN_CELL) return(res)
  m1 <- summary(lm(hZ ~ bmiZ * CNV + factor(sex) + age, data = x))$coefficients
  mc <- summary(lm(hZ ~ bmiZ + factor(sex) + age, data = x[x$CNV == 1, ]))$coefficients
  xt <- x[x$bmi_grp != "mid", ]; xt$hi <- as.numeric(xt$bmi_grp == "high")
  ok_t <- all(table(xt$CNV, xt$hi) >= MIN_CELL)
  m2 <- if (ok_t) summary(lm(hZ ~ hi * CNV + factor(sex) + age, data = xt))$coefficients else NULL
  cbind(res, CNV_main = m1["CNV", 1], CNV_main_se = m1["CNV", 2], CNV_main_p = m1["CNV", 4],
        slope_noncarriers = m1["bmiZ", 1], slope_carriers = mc["bmiZ", 1], slope_carriers_se = mc["bmiZ", 2], slope_carriers_p = mc["bmiZ", 4],
        CNVxBMI = m1["bmiZ:CNV", 1], CNVxBMI_se = m1["bmiZ:CNV", 2], CNVxBMI_p = m1["bmiZ:CNV", 4],
        tails_CNVxhigh = if (ok_t) m2["hi:CNV", 1] else NA, tails_CNVxhigh_se = if (ok_t) m2["hi:CNV", 2] else NA,
        tails_CNVxhigh_p = if (ok_t) m2["hi:CNV", 4] else NA)
}))
write.csv(tests, file.path("adult_anchor_output", paste0("AoU_", ANCESTRY, "_adult_anchor_tests.csv")), row.names = FALSE)
print(tests)
cat("\nDone. Review adult_anchor_output/ (all cells >= ", MIN_CELL, ") before exporting.\n", sep = "")
