# Mean height / head-circumference Z by Fig. 6 age band (2-9, 10-18, 19-39, 40-73), carriers and controls, on the data of
# 17_main_with_ukbb.R (children: external reference; adults: own-cohort non-carriers). Carrier shift = age- and sex-adjusted
# difference from controls within the band (lm y ~ cnv + age + sex). JS 2026-11-02. Output: shift_by_band.csv
source("17_main_with_ukbb.R")
BANDS <- c(2, 10, 19, 40, 74); BLAB <- c("2-9", "10-18", "19-39", "40-73"); out <- list()
for (o in c("Height", "Head circumference")) {
  X <- if (o == "Height") H else dh; C <- if (o == "Height") Hctl_test else dh[locus == "No CNV"]
  for (bb in seq_along(BLAB)) {
    Cb <- C[age >= BANDS[bb] & age < BANDS[bb + 1]]
    out[[length(out) + 1]] <- data.table(outcome = o, band = BLAB[bb], group = "Controls", n = nrow(Cb), n_ukbb = Cb[source == "UK Biobank", .N], mean_z = mean(Cb$y), se = sd(Cb$y) / sqrt(nrow(Cb)))
    for (G in LOCI) {
      Xb <- X[locus == G & age >= BANDS[bb] & age < BANDS[bb + 1]]; if (nrow(Xb) < 10) next
      d <- rbind(Xb, Cb)[, cnv := as.numeric(locus == G)]
      s <- summary(lm(y ~ cnv + age + sex, data = d))$coefficients["cnv", ]
      out[[length(out) + 1]] <- data.table(outcome = o, band = BLAB[bb], group = G, n = nrow(Xb), n_ukbb = Xb[source == "UK Biobank", .N], mean_z = mean(Xb$y), se = sd(Xb$y) / sqrt(nrow(Xb)), shift_vs_controls = s[1], shift_se = s[2], shift_p = s[4])
    }
  }
}
res <- rbindlist(out, fill = TRUE); fwrite(res, file.path(WD, "shift_by_band.csv")); print(res[, .(outcome, band, group, n, n_ukbb, mean_z = round(mean_z, 2), shift = round(shift_vs_controls, 2), p = signif(shift_p, 2))], nrows = 100)
