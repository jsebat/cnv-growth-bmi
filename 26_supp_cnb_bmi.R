# Supplementary figure (JS / Claude 2026-10-24): Penn CNB subtests and totals against BMI in CNV carriers, in the format of
# the SRS/SCQ panels of the main figure (24/25): age-adjusted Z against BMI Z, linear fit with 95% CI, binned means, points
# coloured by BMI tail, within-carrier BMI slope and p. Two figures: accuracy (10 subtests + total) and speed (11 subtests
# + total; speed = -RT z, higher = faster). Columns = measures, rows = CNV groups. Carriers only (CNB non-carriers have no
# BMI). Data: Anjali's cnb_normalized.csv (z vs G2MH non-carriers), BMI from the nearest growth visit (<= 3 y) as built in
# 21_pheno_bmi_gam.R (object `cnb`). Outputs: supp_cnb_bmi_slopes.csv, fig_supp_cnb_accuracy_bmi.{png,pdf},
# fig_supp_cnb_speed_bmi.{png,pdf}, fig_supp_cnb_bmi_legend.txt
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
sink(tempfile()); source(file.path(HERE, "21_pheno_bmi_gam.R")); sink()
suppressPackageStartupMessages({ library(mgcv); library(ggplot2); library(patchwork) })
LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
ACCN <- c(PCET = "cnb_pcet_acc2_z", CPT = "cnb_cpt_ptp_z", LNB = "cnb_lnb_mcr_z", ER40 = "cnb_er40_cr_z", PMAT = "cnb_pmat_pc_z",
          VOLT = "cnb_volt_cr_z", CPF = "cnb_cpf_cr_z", MEDF = "cnb_medf_pc_z", ADT = "cnb_adt_pc_z", PLOT = "cnb_plot_pc_z")
SPDN <- c(MPRAXIS = "cnb_mpraxis_rtcr_z", PCET = "cnb_pcet_rtcr_z", CPT = "cnb_cpt_tprt_z", LNB = "cnb_lnb_mrtc_z", ER40 = "cnb_er40_rtcr_z",
          PMAT = "cnb_pmat_rtcr_z", VOLT = "cnb_volt_rtcr_z", CPF = "cnb_cpf_rtcr_z", MEDF = "cnb_medf_rtcr_z", ADT = "cnb_adt_rtcr_z", PLOT = "cnb_plot_rtcr_z")
TESTNAME <- c(MPRAXIS = "Motor praxis", PCET = "Conditional exclusion (PCET)", CPT = "Continuous performance (CPT)", LNB = "Letter n-back (LNB)", ER40 = "Emotion recognition (ER40)",
              PMAT = "Matrix reasoning (PMAT)", VOLT = "Visual object learning (VOLT)", CPF = "Face memory (CPF)", MEDF = "Emotion differentiation (MEDF)", ADT = "Age differentiation (ADT)", PLOT = "Line orientation (PLOT)")
cc <- cnb[locus %in% LOC & !is.na(BMI) & cnb_age >= 4 & cnb_age <= 65]
long <- function(cols, dom) rbindlist(lapply(names(cols), function(nm) cc[!is.na(get(cols[[nm]])), .(id = GUID, locus, age = cnb_age, BMI, y = if (dom == "speed") -get(cols[[nm]]) else get(cols[[nm]]), domain = dom, measure = nm)]))
L <- rbind(long(ACCN, "accuracy"), long(SPDN, "speed"),
           cc[!is.na(accuracy), .(id = GUID, locus, age = cnb_age, BMI, y = accuracy, domain = "accuracy", measure = "Total")],
           cc[!is.na(speed), .(id = GUID, locus, age = cnb_age, BMI, y = speed, domain = "speed", measure = "Total")])
L[, y_adj := { m <- gam(y ~ s(age, k = 8) + BMI, method = "REML"); y - predict(m, type = "terms")[, "s(age)"] }, by = .(domain, measure, locus)]
L[, tail3 := factor(fifelse(BMI >= qnorm(.7), "BMI \u2265 70th centile", fifelse(BMI <= qnorm(.3), "BMI \u2264 30th centile", "BMI 30th-70th centile")), c("BMI \u2265 70th centile", "BMI 30th-70th centile", "BMI \u2264 30th centile"))]
cat("== n per domain x measure x locus\n"); print(dcast(L[, .N, by = .(domain, measure, locus)], domain + measure ~ locus, value.var = "N", fill = 0), nrows = 50)

# ---------------------------------------------------------------- slopes and heterogeneity
slope <- function(x) tryCatch({ m <- gam(y ~ s(age, k = 8) + BMI, data = x, method = "REML"); co <- summary(m)$p.table["BMI", ]; as.data.table(list(n = nrow(x), slope = unname(co[1]), se = unname(co[2]), p = unname(co[4]))) },
                              error = function(e) as.data.table(list(n = nrow(x), slope = NA_real_, se = NA_real_, p = NA_real_)))
het <- function(x) tryCatch({ x <- copy(x)[, locus := factor(locus)]; m <- gam(y ~ s(age, k = 8, by = locus) + locus + BMI + locus:BMI, data = x, method = "REML"); a <- anova(m)$pTerms.table
                               as.data.table(list(het_F = unname(a["locus:BMI", "F"]), het_p = unname(a["locus:BMI", "p-value"]))) }, error = function(e) as.data.table(list(het_F = NA_real_, het_p = NA_real_)))
MEAS <- unique(L[, .(domain, measure)])
# JS 2026-10-24: pooled rows (locus as a covariate): overgrowth = 16p del + 22q dup, undergrowth = 16p dup + 22q del, all
POOLS <- list("Overgrowth CNVs" = c("16p11.2 del", "22q11.2 dup"), "Undergrowth CNVs" = c("16p11.2 dup", "22q11.2 del"), "All CNVs" = LOC)
ROWS <- c(LOC, names(POOLS)); FIG_ROWS <- LOC   # JS 2026-10-24: pooled rows tried in the figure and reverted; pooled slopes stay in the CSV
members <- function(g) if (g %in% names(POOLS)) POOLS[[g]] else g
slope_pool <- function(x) tryCatch({ m <- gam(y ~ s(age, k = 8) + locus + BMI, data = copy(x)[, locus := factor(locus)], method = "REML"); co <- summary(m)$p.table["BMI", ]; as.data.table(list(n = nrow(x), slope = unname(co[1]), se = unname(co[2]), p = unname(co[4]))) },
                                   error = function(e) as.data.table(list(n = nrow(x), slope = NA_real_, se = NA_real_, p = NA_real_)))
R <- rbindlist(lapply(seq_len(nrow(MEAS)), function(i) { d <- MEAS$domain[i]; mm <- MEAS$measure[i]; x <- L[domain == d & measure == mm]
  rbind(rbindlist(lapply(LOC, function(g) cbind(data.table(domain = d, measure = mm, group = g), slope(x[locus == g])))),
        rbindlist(lapply(names(POOLS), function(g) cbind(data.table(domain = d, measure = mm, group = g), slope_pool(x[locus %in% POOLS[[g]]]), if (g == "All CNVs") het(x) else NULL)), fill = TRUE), fill = TRUE) }))
fwrite(R, file.path(HERE, "supp_cnb_bmi_slopes.csv"))
cat("\n== BMI slopes (SD per SD); pooled row carries the across-locus heterogeneity p\n"); print(R[, .(domain, measure, group, n, slope = round(slope, 2), p = signif(p, 2), het_p = signif(het_p, 2))], nrows = 200)

# ---------------------------------------------------------------- figures
COLS3 <- c("BMI \u2265 70th centile" = "#F0645A", "BMI 30th-70th centile" = "#B8BEC4", "BMI \u2264 30th centile" = "#00A5AD")
BR <- c(-Inf, -1, 0, 1, 2, Inf); XL <- c(-3.2, 4.2)
bins <- function(x) x[, .(bmi = mean(BMI), m = mean(y_adj), se = sd(y_adj) / sqrt(.N), n = .N), by = .(bin = cut(BMI, BR))][n >= 5]
fitline <- function(x) { m <- lm(y_adj ~ BMI, data = x); nd <- data.frame(BMI = seq(max(XL[1], min(x$BMI)), min(XL[2], max(x$BMI)), length.out = 50)); cbind(nd, as.data.frame(predict(m, nd, se.fit = TRUE)[1:2])) }
pan <- function(d, mm, g, yl, top, left) {
  x <- L[domain == d & measure == mm & locus %in% members(g)]; r <- R[domain == d & measure == mm & group == g]
  fl <- fitline(x); bb <- bins(x)
  lab1 <- sub("= -0.00", "= +0.00", sprintf("slope = %+.2f", r$slope), fixed = TRUE); lab2 <- gsub("= < ", "< ", sprintf("p = %s", pfx(r$p)), fixed = TRUE); bold <- isTRUE(r$p < 0.05)
  ggplot(x, aes(BMI, y_adj)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(aes(colour = tail3), size = 1.1, alpha = 0.65) +
    geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), fill = "grey30", alpha = 0.18) + geom_line(data = fl, aes(BMI, fit), colour = "grey10", linewidth = 0.9) +
    geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), colour = "grey10", size = 0.4, linewidth = 0.5, shape = 18) +
    annotate("text", x = XL[2], y = yl[2], label = lab1, hjust = 1, vjust = 1, size = 2.9) +
    annotate("text", x = XL[2], y = yl[2] - 0.085 * diff(yl), label = lab2, hjust = 1, vjust = 1, size = 2.9, fontface = if (bold) "bold" else "plain") +
    scale_colour_manual(values = COLS3, drop = FALSE) + coord_cartesian(xlim = XL, ylim = yl) +
    labs(x = "BMI (Z)", y = if (left) sprintf("%s (Z vs G2MH non-carriers)\nage-adjusted", if (d == "accuracy") "Accuracy" else "Speed (higher = faster)") else NULL, colour = NULL,
         title = if (top) (if (mm == "Total") sprintf("%s total", tools::toTitleCase(d)) else TESTNAME[[mm]]) else NULL, subtitle = sprintf("%s (n = %d)", g, nrow(x))) +
    theme_classic(base_size = 11) + theme(legend.text = element_text(size = 6.5), legend.position = if (top && mm == names(if (d == "accuracy") ACCN else SPDN)[1]) "inside" else "none", legend.position.inside = c(0.33, 0.1),
                                          legend.background = element_blank(), legend.key.height = unit(6, "pt"), plot.title = element_text(face = "bold", size = 9, hjust = 0.5), plot.subtitle = element_text(size = 9.5))
}
make_fig <- function(d) {
  ms <- c(names(if (d == "accuracy") ACCN else SPDN), "Total"); x <- L[domain == d]; r <- quantile(x$y_adj, c(0.005, 0.995)); yl <- r + c(-0.08, 0.30) * diff(r)
  plots <- unlist(lapply(seq_along(ms), function(j) lapply(seq_along(FIG_ROWS), function(i) pan(d, ms[j], FIG_ROWS[i], yl, top = i == 1, left = j == 1) + if (i == 1) labs(tag = LETTERS[j]) else labs(tag = NULL))), recursive = FALSE)
  wrap_plots(plots, ncol = length(ms), byrow = FALSE) & theme(plot.tag = element_text(size = 14, face = "bold"))
}
for (d in c("accuracy", "speed")) {
  f <- make_fig(d); nc <- length(c(names(if (d == "accuracy") ACCN else SPDN), "Total")); w <- 2.1 * nc + 1
  ggsave(file.path(HERE, sprintf("fig_supp_cnb_%s_bmi.png", d)), f, width = w, height = 14, dpi = 130, limitsize = FALSE)
  tryCatch(ggsave(file.path(HERE, sprintf("fig_supp_cnb_%s_bmi.pdf", d)), f, width = w, height = 14, device = cairo_pdf, limitsize = FALSE), error = function(e) ggsave(file.path(HERE, sprintf("fig_supp_cnb_%s_bmi.pdf", d)), f, width = w, height = 14, limitsize = FALSE))
}

# ---------------------------------------------------------------- legend
nr <- function(d, mm) paste(sapply(LOC, function(g) R[domain == d & measure == mm & group == g, n]), collapse = "/")
hetline <- function(d) paste(sapply(c(names(if (d == "accuracy") ACCN else SPDN), "Total"), function(mm) sprintf("%s %s", if (mm == "Total") "total" else mm, pfx(R[domain == d & measure == mm & group == "All CNVs", het_p]))), collapse = ", ")
leg <- c(
"Supplementary Figure. Body mass and Penn Computerized Neurocognitive Battery (CNB) performance in CNV carriers.",
"Two figures, accuracy (10 subtests and the accuracy total) and speed (11 subtests and the speed total). Columns: CNB measures; rows: 16p11.2 deletion, 16p11.2 duplication, 22q11.2 deletion and 22q11.2 duplication carriers (G2MH). Scores are z relative to G2MH non-carriers (Anjali's cnb_normalized), speed as -z of the median correct response time (higher = faster); totals are the mean of subtest z (at least 5 subtests). Each panel plots the score, adjusted for age within group (age smooth removed), against BMI Z (child WHO; adult centile vs pooled non-carrier adults), with BMI measured within 3 years of the CNB. Points: BMI at or above the 70th centile (red), 30th-70th (grey), at or below the 30th (teal); black line and band: linear fit with 95% CI; diamonds: mean +/- SE in BMI bins (< -1, -1 to 0, 0 to 1, 1 to 2, > 2 SD). Printed: the BMI slope (SD per SD of BMI) from a GAM y ~ s(age) + BMI and its p-value (bold: p < 0.05). G2MH non-carriers with CNB have no BMI, so no control row or carrier-vs-control contrast is shown.",
sprintf("n (16p del/16p dup/22q del/22q dup): accuracy total %s; speed total %s.", nr("accuracy", "Total"), nr("speed", "Total")),
sprintf("Heterogeneity of the BMI slope across the four CNV groups (Wald test of locus x BMI, within carriers): accuracy %s; speed %s.", hetline("accuracy"), hetline("speed")),
"All slopes, including pooled estimates (overgrowth, undergrowth, all CNVs; locus as a covariate), are in supp_cnb_bmi_slopes.csv.")
writeLines(leg, file.path(HERE, "fig_supp_cnb_bmi_legend.txt"))
cat("done\n")

# ---------------------------------------------------------------- JS 2026-10-24: accuracy total only, one panel per CNV group
# The across-CNV heterogeneity of the BMI slope is nominal for the accuracy total (and emotion differentiation) but not for
# the other subtests or any speed measure, so the total gets its own supplementary figure with the heterogeneity test in the legend.
xt <- L[domain == "accuracy" & measure == "Total"]; rt <- quantile(xt$y_adj, c(0.005, 0.995)); ylt <- rt + c(-0.08, 0.30) * diff(rt)
pt <- lapply(seq_along(LOC), function(i) pan("accuracy", "Total", LOC[i], ylt, top = FALSE, left = i == 1) + labs(tag = LETTERS[i]) +
               theme(plot.subtitle = element_text(size = 11), legend.position = if (i == 1) "inside" else "none", legend.position.inside = c(0.3, 0.1)))
ft <- wrap_plots(pt, nrow = 1) + plot_annotation(title = "Penn CNB accuracy total by BMI in CNV carriers", theme = theme(plot.title = element_text(face = "bold", size = 13))) & theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(HERE, "fig_supp_cnb_accuracy_total_bmi.png"), ft, width = 12, height = 4.2, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_supp_cnb_accuracy_total_bmi.pdf"), ft, width = 12, height = 4.2, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_supp_cnb_accuracy_total_bmi.pdf"), ft, width = 12, height = 4.2))
RT <- R[domain == "accuracy" & measure == "Total"]
sl <- function(g) { r <- RT[group == g]; sprintf("%s %+.2f (SE %.2f, p = %s, n = %d)", g, r$slope, r$se, pfx(r$p), r$n) }
legt <- c(
"Supplementary Figure. Body mass and Penn CNB accuracy in CNV carriers.",
sprintf("Penn Computerized Neurocognitive Battery accuracy total (mean z across subtests relative to G2MH non-carriers; at least 5 subtests), adjusted for age within group (age smooth removed), against BMI Z (child WHO; adult centile vs pooled non-carrier adults; BMI within 3 years of the CNB), in carriers of (A) 16p11.2 deletion, (B) 16p11.2 duplication, (C) 22q11.2 deletion and (D) 22q11.2 duplication (G2MH). Points: BMI at or above the 70th centile (red), 30th-70th (grey), at or below the 30th (teal); black line and band: linear fit with 95%% CI; diamonds: mean +/- SE in BMI bins (< -1, -1 to 0, 0 to 1, 1 to 2, > 2 SD). Printed: the BMI slope (SD of accuracy per SD of BMI) from a GAM y ~ s(age) + BMI and its p-value (bold: p < 0.05). G2MH non-carriers with CNB have no BMI, so no control panel or carrier-vs-control contrast is shown."),
sprintf("Slopes: %s; %s; %s; %s. Pooled with locus as a covariate: all carriers %+.2f (p = %s), overgrowth CNVs (16p11.2 del + 22q11.2 dup) %+.2f (p = %s), undergrowth CNVs (16p11.2 dup + 22q11.2 del) %+.2f (p = %s).",
        sl("16p11.2 del"), sl("16p11.2 dup"), sl("22q11.2 del"), sl("22q11.2 dup"), RT[group == "All CNVs", slope], pfx(RT[group == "All CNVs", p]), RT[group == "Overgrowth CNVs", slope], pfx(RT[group == "Overgrowth CNVs", p]), RT[group == "Undergrowth CNVs", slope], pfx(RT[group == "Undergrowth CNVs", p])),
{ hp <- RT[group == "All CNVs", het_p]; sprintf("Heterogeneity of the BMI slope across the four CNV groups (Wald test of the locus x BMI interaction in a GAM with a separate age smooth per locus): F = %.2f on 3 df, p = %s%s. Before plausibility QC of height and BMI (00_qc_growth.R) this test gave p = 0.027, driven by one 16p11.2 duplication carrier with a recorded BMI of 9.2 kg/m2; after QC no accuracy subtest or speed measure shows significant heterogeneity (fig_supp_cnb_bmi_legend.txt).",
        RT[group == "All CNVs", het_F], pfx(hp), if (isTRUE(hp < 0.05)) ", driven by the negative slope in 16p11.2 duplication carriers" else "; the negative slope in 16p11.2 duplication carriers is suggestive only") })
writeLines(legt, file.path(HERE, "fig_supp_cnb_accuracy_total_legend.txt"))
cat("accuracy-total figure written\n")
