# Numbers for the Fig. 6 text and Supplementary Table 22 (JS 2026-10-27), on the data of the main figure v5 (17_main_with_ukbb.R:
# UK Biobank adults, adults scaled to own-cohort non-carriers). RUNS ON EXPANSE (sources 17).
# 1. Sample sizes behind each test (clinical / UK Biobank, carriers and controls).
# 2. Heterogeneity of the age-adjusted BMI slope across the four CNVs, within carriers: Wald test of locus x BMI in
#    y ~ s(age, by = locus) + locus + BMI + locus:BMI (as in 25_main_figure_v3.R), height and head circumference.
# 3. Carrier mean shift vs controls, children (< 18 y) and adults (>= 18 y), age-adjusted: y ~ CNV + s(age) in each period.
# Output: fig6_numbers.csv (long format: quantity, outcome, group, value columns)
source("17_main_with_ukbb.R")
out <- list()
for (o in c("Height", "Head circumference")) {
  X <- if (o == "Height") H else dh; C <- if (o == "Height") Hctl_test else dh[locus == "No CNV"]
  for (G in LOCI) out[[length(out) + 1]] <- data.table(quantity = "n", outcome = o, group = G, n_clinical = X[locus == G & source == "Clinical", .N], n_ukbb = X[locus == G & source == "UK Biobank", .N])
  out[[length(out) + 1]] <- data.table(quantity = "n", outcome = o, group = "Controls (tests)", n_clinical = C[source == "Clinical", .N], n_ukbb = C[source == "UK Biobank", .N])
  car <- X[locus %in% LOCI][, locus := factor(locus)]
  m <- gam(y ~ s(age, k = 8, by = locus) + locus + BMI + locus:BMI, data = car, method = "REML"); a <- anova(m)$pTerms.table
  out[[length(out) + 1]] <- data.table(quantity = "heterogeneity", outcome = o, group = "4 CNVs", F = a["locus:BMI", "F"], df = a["locus:BMI", "df"], p = a["locus:BMI", "p-value"])
  for (G in LOCI) for (per in c("child", "adult")) {
    d <- rbind(X[locus == G], C)[if (per == "child") age < 18 else age >= 18][, cnv := as.numeric(locus == G)]
    if (d[cnv == 1, .N] < 10) next
    s <- summary(gam(y ~ cnv + s(age, k = 8), data = d, method = "REML"))$p.table["cnv", ]
    out[[length(out) + 1]] <- data.table(quantity = paste0("shift_", per), outcome = o, group = G, n_carriers = d[cnv == 1, .N], est = s[1], se = s[2], p = s[4])
  }
}
# 4. Omnibus test across the five genotype groups (JS 2026-10-27, Jon: "does the effect of BMI differ by genotype?"):
#    controls + four CNVs pooled; y ~ group + s(age, by = group) + s(age, by = BMI) + sum_k (CNV_k x BMI); the four CNV_k x BMI
#    terms are each CNV's departure from the control BMI slope (age-adjusted; the BMI effect may vary with age, shared by all
#    groups). Joint 4-df Wald test (chi-square on Vp). Also within age bands: lm(y ~ group * BMI + sex + age), F-test of
#    group:BMI (4 df, or fewer where a CNV has < 15 carriers in the band).
BANDS <- c(2, 10, 19, 40, 74); BLAB <- c("2-9", "10-18", "19-39", "40-73")
for (o in c("Height", "Head circumference")) {
  X <- if (o == "Height") H else dh; C <- if (o == "Height") Hctl_test else dh[locus == "No CNV"]
  d <- rbind(C[, .(y, age, BMI, sex, group = "Controls")], X[locus %in% LOCI, .(y, age, BMI, sex, group = locus)])[, group := factor(group, c("Controls", LOCI))]
  gk <- paste0("g", seq_along(LOCI), "BMI"); for (k in seq_along(LOCI)) set(d, j = gk[k], value = as.numeric(d$group == LOCI[k]) * d$BMI)
  m <- bam(as.formula(paste("y ~ group + s(age, k = 10, by = group) + s(age, k = 10, by = BMI) +", paste(gk, collapse = " + "))), data = d, method = "REML")
  i <- match(gk, names(coef(m))); b <- coef(m)[i]; V <- m$Vp[i, i]; W <- drop(t(b) %*% solve(V, b))
  out[[length(out) + 1]] <- data.table(quantity = "omnibus_5group", outcome = o, group = "all ages", chisq = W, df = length(gk), p = pchisq(W, length(gk), lower.tail = FALSE))
  for (k in seq_along(LOCI)) out[[length(out) + 1]] <- data.table(quantity = "omnibus_contrast", outcome = o, group = LOCI[k], est = b[k], se = sqrt(V[k, k]), p = 2 * pnorm(-abs(b[k] / sqrt(V[k, k]))))
  for (bb in seq_along(BLAB)) {
    db <- d[age >= BANDS[bb] & age < BANDS[bb + 1]]; keep <- names(which(table(db$group) >= 15)); db <- db[group %in% keep][, group := droplevels(group)]
    if (!"Controls" %in% keep || nlevels(db$group) < 2) next
    m1 <- lm(y ~ group * BMI + sex + age, data = db); m0 <- lm(y ~ group + BMI + sex + age, data = db); a <- anova(m0, m1)
    out[[length(out) + 1]] <- data.table(quantity = "omnibus_5group_band", outcome = o, group = BLAB[bb], F = a$F[2], df = a$Df[2], p = a$`Pr(>F)`[2], groups = paste(setdiff(keep, "Controls"), collapse = "; "))
  }
}
res <- rbindlist(out, fill = TRUE); fwrite(res, file.path(WD, "fig6_numbers.csv")); print(res, nrows = 100)
cat("fig6 numbers written\n")
