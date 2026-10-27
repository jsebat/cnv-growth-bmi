# Check of the prime panels (JS 2026-10-27, Jon's question on 16p11.2 del height): is the BMI effect on height in adult
# carriers more negative than in adult controls, as the raw panel suggests, and does the fitted age x BMI surface agree?
# No smoothing: within each age bin (2-9, 10-18, 19-39, 40-73 y), height Z ~ BMI Z + sex + age (linear) in carriers and in
# controls separately, and the carrier - control difference from the pooled model height ~ CNV * BMI + sex + age.
# Also the mean of the fitted surface slope (main_with_ukbb_bmi_effect_by_age.csv) over the same bins, for comparison.
# RUNS ON EXPANSE after 17_main_with_ukbb.R (sources it for H, Hctl_test, dh, bea). Output: check_bmi_effect_by_agebin.csv
source("17_main_with_ukbb.R")
BINS <- c(2, 10, 19, 40, 74); LAB <- c("2-9", "10-18", "19-39", "40-73")
bin_slopes <- function(x_c, ctl, outcome) {
  rbindlist(lapply(seq_along(LAB), function(b) {
    xc <- x_c[age >= BINS[b] & age < BINS[b + 1]]; xk <- ctl[age >= BINS[b] & age < BINS[b + 1]]
    if (nrow(xc) < 15 || nrow(xk) < 30) return(NULL)
    sc <- summary(lm(y ~ BMI + sex + age, data = xc))$coefficients["BMI", ]; sk <- summary(lm(y ~ BMI + sex + age, data = xk))$coefficients["BMI", ]
    p <- rbind(xc[, .(y, BMI, sex, age, cnv = 1)], xk[, .(y, BMI, sex, age, cnv = 0)])
    sd <- summary(lm(y ~ cnv * BMI + sex + age, data = p))$coefficients["cnv:BMI", ]
    data.table(outcome = outcome, age_bin = LAB[b], n_carriers = nrow(xc), n_controls = nrow(xk),
               slope_carriers = sc[1], se_carriers = sc[2], slope_controls = sk[1], se_controls = sk[2],
               diff = sd[1], se_diff = sd[2], p_diff = sd[4])
  }))
}
res <- rbindlist(lapply(LOCI, function(L) {
  rbind(cbind(group = L, bin_slopes(H[locus == L], Hctl_test, "Height")),
        cbind(group = L, bin_slopes(dh[locus == L], dh[locus == "No CNV"], "Head circumference")))
}), fill = TRUE)
# fitted-surface slope averaged over the same bins (from bea: outcome, group, who, age, slope)
surf <- bea[, .(age_bin = LAB[findInterval(age, BINS)], outcome, group, who, slope)][!is.na(age_bin) & !is.na(slope),
            .(surface_slope = mean(slope)), by = .(outcome, group, who, age_bin)]
surf <- dcast(surf, outcome + group + age_bin ~ who, value.var = "surface_slope"); setnames(surf, c("Carriers", "Non-carriers"), c("surface_carriers", "surface_controls"), skip_absent = TRUE)
res <- merge(res, surf, by = c("outcome", "group", "age_bin"), all.x = TRUE)
res[, age_bin := factor(age_bin, LAB)]; setorder(res, outcome, group, age_bin)
fwrite(res, file.path(WD, "check_bmi_effect_by_agebin.csv"))
cat("\n== BMI slope on the outcome by age bin (s.d. per s.d.): carriers vs controls, and the fitted surface in the same bins\n")
print(res[, .(outcome, group, age_bin, n_carriers, n_controls, carriers = sprintf("%+.2f (%.2f)", slope_carriers, se_carriers),
              controls = sprintf("%+.2f (%.2f)", slope_controls, se_controls), diff = sprintf("%+.2f", diff), p_diff = signif(p_diff, 2),
              surf_carriers = round(surface_carriers, 2), surf_controls = round(surface_controls, 2))], nrows = 100)
cat("check written\n")
