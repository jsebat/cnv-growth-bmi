# Quick look (JS 2026-10-24, Jon): 16p11.2 deletion only, BMI split at the 90th centile (> 90th vs <= 90th) instead of the
# 70/30 tails, because the deletion's low-BMI tail is small. Outcomes: height and head circumference (clinical, ages 2-40,
# scored as in the main figure) and the six phenotypes of 21_pheno_bmi_gam.R (its data build is sourced; tests and figure
# skipped). Per outcome: n per group, age-adjusted group difference (gam y ~ hi90 + s(age)), and age x group (gam with
# s(age, by = group), p of the by-smooth). Runs locally. Output: bmi90_16pdel_tests.csv, fig_16pdel_bmi90.{png,pdf}
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "21_pheno_bmi_gam.R"))   # gives D (phenotypes with BMI Z), gm, bmiZ(), score_height/score_hc, pfx
# ---- height and head circumference for 16p11.2 del (clinical, latest visit with the outcome and BMI, ages 2-40) ----
g16 <- gm[cohort != "Searchlight" & CNV == "16p11.2 proximal DEL" & !is.na(BMI_raw) & age >= 2 & age <= 40]
g16[, id := ifelse(cohort == "SVIP", paste0("SVIP:", sampleID), paste(cohort, MSassigned_ID, sep = ":"))]
ht <- g16[!is.na(HEIGHT_raw)][order(id, -age)][, .SD[1], by = id][, .(id, age, y = score_height(HEIGHT_raw, age, sex), BMI = bmiZ(BMI_raw, age, sex), pheno = "Height")]
hc <- g16[!is.na(HEADCIRCUM_raw)][order(id, -age)][, .SD[1], by = id][, .(id, age, y = score_hc(HEADCIRCUM_raw, age, sex), BMI = bmiZ(BMI_raw, age, sex), pheno = "Head circumference")]
X <- rbind(ht, hc, D[locus == "16p11.2 del", .(id, age, y, BMI, pheno)])[!is.na(y) & !is.na(BMI)]
X[, grp := factor(fifelse(BMI > qnorm(.9), "BMI > 90th centile", "BMI \u2264 90th centile"), levels = c("BMI \u2264 90th centile", "BMI > 90th centile"))]
PH <- c("Height", "Head circumference", "SRS", "SCQ", "RBS-R", "FSIQ", "CNB accuracy", "CNB speed")
tests <- rbindlist(lapply(PH, function(ph) {
  x <- X[pheno == ph]; x[, hi := as.numeric(grp == levels(grp)[2])]
  r <- data.table(pheno = ph, n_le90 = sum(x$hi == 0), n_gt90 = sum(x$hi == 1), n_le30 = sum(x$BMI <= qnorm(.3)), n_ge70 = sum(x$BMI >= qnorm(.7)))
  if (sum(x$hi == 1) < 8 || sum(x$hi == 0) < 8) return(r)
  m1 <- tryCatch(gam(y ~ hi + s(age, k = 8), data = x, method = "REML"), error = function(e) NULL)
  m2 <- tryCatch(gam(y ~ hi + s(age, k = 8) + s(age, k = 6, by = hi), data = x, method = "REML"), error = function(e) NULL)
  if (!is.null(m1)) { co <- summary(m1)$p.table["hi", ]; r[, `:=`(diff_gt90_vs_le90 = co[1], se = co[2], p_diff = co[4])] }
  if (!is.null(m2)) r[, p_age_x_group := summary(m2)$s.table["s(age):hi", "p-value"]]
  r
}))
fwrite(tests, file.path(HERE, "bmi90_16pdel_tests.csv")); print(tests[, .(pheno, n_le90, n_gt90, n_le30, n_ge70, diff = round(diff_gt90_vs_le90, 2), p_diff = signif(p_diff, 2), p_age_x_group = signif(p_age_x_group, 2))])
# ---- figure ----
CL <- c("BMI \u2264 90th centile" = "#4B6EA9", "BMI > 90th centile" = "#F0645A")
YL <- c(Height = "Height (SD)", `Head circumference` = "Head circumference (SD)", SRS = "SRS (Z vs SSC unaffected)", SCQ = "SCQ (Z vs SSC siblings)", `RBS-R` = "RBS-R 29 items (Z)", FSIQ = "FSIQ (Z)", `CNB accuracy` = "CNB accuracy (Z)", `CNB speed` = "CNB speed (Z)")
pan <- function(ph) {
  x <- X[pheno == ph]; r <- tests[pheno == ph]; xmax <- if (ph %in% c("Height", "Head circumference")) 40 else 60
  lab <- gsub("= < ", "< ", sprintf("> 90 vs \u2264 90: %+.2f, p = %s\nage \u00d7 group p = %s", r$diff_gt90_vs_le90, pfx(r$p_diff), pfx(r$p_age_x_group)), fixed = TRUE)
  yl <- range(quantile(x$y, c(0.005, 0.995))); yl <- yl + c(-0.5, 0.9) * diff(yl) * 0.22
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") + geom_point(size = 1.3, alpha = 0.65) +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 6), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = xmax, y = yl[2], label = lab, hjust = 1, vjust = 1, size = 3.1, fontface = if (isTRUE(min(r$p_diff, r$p_age_x_group, na.rm = TRUE) < 0.05)) "bold" else "plain") +
    scale_colour_manual(values = CL) + scale_fill_manual(values = CL) + coord_cartesian(xlim = c(2, xmax), ylim = yl) +
    labs(x = "Age (years)", y = YL[[ph]], colour = NULL, fill = NULL, title = ph, subtitle = sprintf("16p11.2 del: n = %d (\u2264 90th) / %d (> 90th)", r$n_le90, r$n_gt90)) +
    theme_classic(base_size = 12) + theme(legend.position = "inside", legend.position.inside = c(0.25, 0.1), legend.text = element_text(size = 8), legend.background = element_blank(), plot.title = element_text(face = "bold"))
}
fig <- wrap_plots(lapply(PH, pan), ncol = 4) + plot_annotation(title = "16p11.2 deletion carriers: BMI above vs at or below the 90th centile", theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(HERE, "fig_16pdel_bmi90.png"), fig, width = 17, height = 9, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_16pdel_bmi90.pdf"), fig, width = 17, height = 9, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_16pdel_bmi90.pdf"), fig, width = 17, height = 9))
cat("done\n")
