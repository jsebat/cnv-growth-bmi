# JS 2026-10-23: what drives the CNV x age x BMI interaction (panel B of fig_main_growth_ukbb)?
# Refit the three-way model used for the "CNV x age x BMI vs controls" test (threeway_p in 17_main_with_ukbb.R) for each
# CNV (height, clinical + UKBB) and evaluate the BMI-on-height slope as a function of age, separately for carriers and
# controls, with the carrier - control difference. Slope(age) = [E(height | BMI = +1 SD) - E(height | BMI = -1 SD)] / 2 at
# that age, from the fitted surface (s(BMI) + ti(age, BMI) + carrier-specific ti). SEs from the coefficient covariance.
# Also the predicted height at BMI = +1 / -1 SD by age for carriers and controls.
# Output: ti_surface_by_age.csv, fig_ti_surface_by_age.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
set.seed(1)
WD <- getwd()
ex <- parse(file.path(WD, "17_main_with_ukbb.R"))
k <- which(sapply(ex, function(e) grepl("Hctl_plot <-", paste(deparse(e), collapse = ""))))[1]
for (i in seq_len(k)) eval(ex[[i]], envir = globalenv())

AGES <- c(3, 5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 40, 50, 60, 70)
fit_one <- function(L) {
  x <- rbind(Hctl_test, H[locus == L])[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- bam(y ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) + ti(age, BMI, k = c(10, 10)) +
             ti(age, BMI, k = c(10, 10), by = is_cnv), data = x, method = "REML")
  p3 <- summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]
  nd <- function(a, b, g) data.frame(age = a, BMI = b, is_cnv = factor(g, levels = c("No CNV", "CNV")))
  V <- vcov(m); B <- coef(m)
  rows <- rbindlist(lapply(AGES, function(a) {
    Xc1 <- predict(m, nd(a, 1, "CNV"), type = "lpmatrix");  Xc0 <- predict(m, nd(a, -1, "CNV"), type = "lpmatrix")
    Xn1 <- predict(m, nd(a, 1, "No CNV"), type = "lpmatrix"); Xn0 <- predict(m, nd(a, -1, "No CNV"), type = "lpmatrix")
    sc <- (Xc1 - Xc0) / 2; sn <- (Xn1 - Xn0) / 2; sd_ <- sc - sn
    est <- function(X) c(drop(X %*% B), sqrt(drop(X %*% V %*% t(X))))
    ec <- est(sc); en <- est(sn); ed <- est(sd_)
    # carrier counts near this age (within +/- 2 y) by BMI tail, for context
    nc <- x[is_cnv == "CNV" & abs(age - a) <= 2]
    data.table(locus = L, age = a, slope_carriers = ec[1], se_carriers = ec[2], slope_controls = en[1], se_controls = en[2],
               diff = ed[1], se_diff = ed[2], p_diff = 2 * pnorm(-abs(ed[1] / ed[2])),
               pred_carrier_hiBMI = drop(Xc1 %*% B), pred_carrier_loBMI = drop(Xc0 %*% B),
               pred_control_hiBMI = drop(Xn1 %*% B), pred_control_loBMI = drop(Xn0 %*% B),
               n_carriers_pm2y = nrow(nc), n_carriers_hi_pm2y = nc[BMI >= qnorm(.7), .N], n_carriers_lo_pm2y = nc[BMI <= qnorm(.3), .N],
               threeway_p = p3)
  }))
  cat(L, "threeway p =", signif(p3, 3), "\n"); rows
}
res <- rbindlist(lapply(LOCI, fit_one))
fwrite(res, file.path(WD, "ti_surface_by_age.csv"))
print(res[, .(locus, age, slope_carriers = round(slope_carriers, 2), slope_controls = round(slope_controls, 2), diff = round(diff, 2), p_diff = signif(p_diff, 2),
              n_hi = n_carriers_hi_pm2y, n_lo = n_carriers_lo_pm2y)], nrows = 100)

# ---- figure: slope of height on BMI by age, carriers vs controls, per locus; and the difference ----
long <- rbind(res[, .(locus, age, group = "Carriers", slope = slope_carriers, se = se_carriers)],
              res[, .(locus, age, group = "Non-carriers", slope = slope_controls, se = se_controls)])
long[, locus := factor(locus, LOCI)]
p1 <- ggplot(long, aes(age, slope, colour = group, fill = group)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), alpha = 0.18, colour = NA) + geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ locus, nrow = 1) + scale_colour_manual(values = c(Carriers = "#4B6EA9", `Non-carriers` = "#8C9A9A")) +
  scale_fill_manual(values = c(Carriers = "#4B6EA9", `Non-carriers` = "#8C9A9A")) +
  labs(x = "Age (years)", y = "Height per +1 SD BMI (SD)", colour = NULL, fill = NULL,
       title = "BMI effect on height by age, from the CNV x age x BMI model (height, clinical + UKBB)",
       subtitle = "Slope = [E(height | BMI +1 SD) - E(height | BMI -1 SD)] / 2 at each age; 95% CI") +
  theme_classic(base_size = 14) + theme(legend.position = "top")
res[, locus := factor(locus, LOCI)]
p2 <- ggplot(res, aes(age, diff)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(aes(ymin = diff - 1.96 * se_diff, ymax = diff + 1.96 * se_diff), alpha = 0.18, fill = "#D7191C") + geom_line(colour = "#D7191C", linewidth = 1) +
  geom_point(aes(size = n_carriers_pm2y), colour = "#D7191C") + scale_size_area(max_size = 5, name = "carriers within 2 y") +
  geom_text(aes(label = ifelse(p_diff < 0.05, "*", "")), vjust = -0.6, size = 6) +
  facet_wrap(~ locus, nrow = 1) +
  labs(x = "Age (years)", y = "Carrier - control slope (SD per SD BMI)", title = "Difference in the BMI effect: carriers minus non-carriers (the term tested by CNV x age x BMI); * p < 0.05 at that age") +
  theme_classic(base_size = 14) + theme(legend.position = "top")
fig <- p1 / p2
ggsave(file.path(WD, "fig_ti_surface_by_age.png"), fig, width = 16, height = 10, dpi = 170)
tryCatch(ggsave(file.path(WD, "fig_ti_surface_by_age.pdf"), fig, width = 16, height = 10, device = cairo_pdf),
         error = function(e) ggsave(file.path(WD, "fig_ti_surface_by_age.pdf"), fig, width = 16, height = 10))
cat("done\n")
