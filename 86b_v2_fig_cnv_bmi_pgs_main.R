# 86b_v2_fig_cnv_bmi_pgs_main.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 3: MAIN figure)
# What it does: main figure of section 3 (CNV x BMI x autism polygenic score -> autistic traits; exploratory), four panels:
#   a  SRS against BMI in genotyped carriers below vs above the median autism score (partial residuals: cohort x CNV, age, sex
#      removed; fitted lines with 95% bands), SSC unaffected relatives as a thin grey reference line
#      Printed tests (2026-11-21): BMI x score within carriers (cohort-adjusted model) and the three-way term vs relatives
#   b  BMI slope on SRS by fifth of the autism score (95% CI) with the continuous interaction line, carriers and relatives
#   c  forest of the three-way term for SRS: pooled, per CNV, per cohort, leave-one-cohort-out, leave-one-CNV-out, and SCQ pooled
#   d  effect of BMI on SRS / SCQ: observational estimate next to Mendelian-randomisation estimates (2SLS with three adiposity
#      scores, Wald ratio per score), all genotyped carriers and SSC unaffected relatives
# Inputs : 86_v2_cnv_bmi_pgs_stats.R outputs (cache pgs_stats_2026-11-20.RData) and v2_base() objects. No new statistics here.
# Outputs: autism/fig_cnv_bmi_pgs_main.{png,pdf}, autism/fig_cnv_bmi_pgs_main_plotted.csv.
# null model: none here (tests and the permutation null are in 86_v2_cnv_bmi_pgs_stats.R).
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
v2_base(); load(file.path("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data", "v2_autism_cache", "pgs_stats_2026-11-20.RData"))
BS <- 16; TS <- 4.6; C_LOW <- "#4B6EA9"; C_HIGH <- "#D9822B"
A3 <- "three-way CNV x BMI x autism score vs SSC relatives"
# ---------------------------------------------------------------- a
X3 <- XG[pheno == "SRS" & is_cnv == 1]; med <- median(X3$pgs); X3[, half := fifelse(pgs < med, "lo", "hi")]
m3 <- lm(y ~ cg + poly(age, 2) + male + BMI * pgs, X3); Xm <- model.matrix(m3); nu <- grep("poly\\(age|male|^cg", colnames(Xm)); cf <- coef(m3); cf[is.na(cf)] <- 0
ct <- as.vector(Xm[, nu, drop = FALSE] %*% cf[nu]); X3[, ya := y - ct + mean(ct)]
sl <- X3[, { r <- rob(lm(y ~ cg + poly(age, 2) + male + BMI, .SD), "BMI"); .(n = .N, slope = r[["b"]], se = r[["se"]], p = r[["p"]]) }, by = half]
lab <- c(lo = sprintf("Carriers, autism score below median (n = %d): slope = %+.2f", sl[half == "lo", n], sl[half == "lo", slope]), hi = sprintf("Carriers, autism score above median (n = %d): slope = %+.2f", sl[half == "hi", n], sl[half == "hi", slope]))
X3[, hl := factor(lab[half], lab)]
f3 <- X3[, { f <- lm(ya ~ BMI); g <- seq(quantile(BMI, 0.02), quantile(BMI, 0.98), length.out = 40); pr <- predict(f, data.frame(BMI = g), se.fit = TRUE); .(BMI = g, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }, by = hl]
rel <- RELf(XG[pheno == "SRS"]); mr0 <- lm(y ~ poly(age, 2) + male + BMI, rel); rl <- data.table(BMI = c(-2.8, 3.8))[, fit := mean(rel$y) + coef(mr0)[["BMI"]] * (BMI - mean(rel$BMI))]
HC2 <- setNames(c(C_LOW, C_HIGH), lab)
# tests printed on panel a (JS 2026-11-21): continuous score, from 86_v2_cnv_bmi_pgs_stats.R (RB = threeway_pgs_robustness.csv)
ra2 <- RB[trait == "SRS" & analysis == "within carriers BMI x autism score (cohort x BMI and cohort x score terms)" & subset == "all carriers"]
ra3 <- RB[trait == "SRS" & analysis == A3 & subset == "all carriers"]; stopifnot(nrow(ra2) == 1, nrow(ra3) == 1)
pa <- ggplot(X3, aes(BMI, ya, colour = hl)) + geom_point(size = 1.5, alpha = 0.45) +
  geom_line(data = rl, aes(BMI, fit), inherit.aes = FALSE, colour = C_NON, linewidth = 0.9) +
  geom_ribbon(data = f3, aes(BMI, ymin = lo, ymax = hi, fill = hl), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = f3, aes(BMI, fit), linewidth = 1.5) +
  annotate("label", x = 3.8, y = rl$fit[2] - 0.3, label = sprintf("unaffected relatives (n = %d): slope = %+.2f", nrow(rel), coef(mr0)[["BMI"]]), hjust = 1, vjust = 1, size = TS, colour = "grey30", fill = "white", alpha = 0.85, linewidth = 0) +
  annotate("label", x = -2.75, y = 5.75, hjust = 0, vjust = 1, size = TS, lineheight = 1.05, colour = "grey10", fill = "white", alpha = 0.85, linewidth = 0, fontface = "bold",
           label = sprintf("BMI \u00d7 autism score in carriers = %+.3f, %s\nCNV \u00d7 BMI \u00d7 autism score (carriers vs relatives) = %+.3f, %s", ra2$b, plab(ra2$p), ra3$b, plab(ra3$p))) +
  scale_colour_manual(values = HC2, name = NULL) + scale_fill_manual(values = HC2, guide = "none") + coord_cartesian(xlim = c(-2.8, 3.8), ylim = c(-1.2, 5.8)) + guides(colour = guide_legend(ncol = 1, override.aes = list(alpha = 1, size = 3))) +
  labs(x = "BMI (SD)", y = "Social Responsiveness Scale (SD)", tag = "a", title = "BMI and autistic traits in carriers, by autism polygenic score") + th_v2(BS) + theme(legend.position = "bottom", legend.text = element_text(size = 14.5))
# ---------------------------------------------------------------- b
GC <- c("CNV carriers" = C_CAR, "SSC unaffected relatives" = C_NON)
fq <- FQ[trait == "SRS"]; fl <- FL[trait == "SRS" & abs(score) <= 2]; r3 <- RB[trait == "SRS" & analysis == A3 & subset == "all carriers"]
pb <- ggplot(fq, aes(score_mean, slope, colour = group)) + geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_ribbon(data = fl, aes(score, ymin = slope - 1.96 * se, ymax = slope + 1.96 * se, fill = group), inherit.aes = FALSE, alpha = 0.16) + geom_line(data = fl, aes(score, slope), linewidth = 1.2) +
  geom_pointrange(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), size = 0.9, linewidth = 0.9, position = position_dodge(width = 0.08)) +
  annotate("text", x = 2, y = 0.5, hjust = 1, vjust = 1, size = TS + 0.2, lineheight = 1.05, fontface = face(r3$p),
           label = sprintf("CNV \u00d7 BMI \u00d7 autism score = %+.3f\n%s (permutation %s)", r3$b, plab(r3$p), plab(r3$perm_p))) +
  scale_colour_manual(values = GC, name = NULL, labels = c("CNV carriers" = sprintf("CNV carriers (n = %d)", sum(fq[group == "CNV carriers", n])), "SSC unaffected relatives" = sprintf("Unaffected relatives (n = %d)", sum(fq[group != "CNV carriers", n])))) + scale_fill_manual(values = GC, guide = "none") +
  coord_cartesian(ylim = c(-0.25, 0.5)) + labs(x = "Autism polygenic score (SD)", y = "BMI slope on Social Responsiveness Scale\n(SD per SD of BMI)", tag = "b", title = "The BMI slope falls as the autism polygenic score rises",
       subtitle = "Points: fifths of the score (95% CI); lines: continuous interaction model") + th_v2(BS) + theme(legend.position = "bottom", legend.text = element_text(size = 14.5))
# ---------------------------------------------------------------- c
F <- rbind(RB[trait == "SRS" & analysis == A3 & !is.na(b)], RB[trait == "SCQ" & analysis == A3 & subset == "all carriers"])
F[, sec := fcase(subset == "all carriers" & trait == "SRS", "Pooled", grepl("^CNV:", subset), "By CNV", grepl("^cohort:", subset), "By cohort", grepl("^leave out cohort", subset), "Leave one\ncohort out", grepl("^leave out CNV", subset), "Leave one\nCNV out", default = "SCQ")]
F[, lab := sub("^CNV: |^cohort: |^leave out cohort: |^leave out CNV: ", "", subset)][sec == "Pooled", lab := "All carriers"][sec == "SCQ", lab := "All carriers (SCQ)"][grepl("^Leave", sec), lab := paste("without", lab)]
F[, lab := sprintf("%s  (n = %d)", lab, n_carriers)][, sec := factor(sec, c("Pooled", "By CNV", "By cohort", "Leave one\ncohort out", "Leave one\nCNV out", "SCQ"))]
F[, lab := factor(lab, rev(unique(lab)))][, ptxt := ifelse(!is.na(perm_p), sprintf("%s (permutation %s)", plab(p), pfmt(perm_p)), plab(p))]
xr <- c(min(F$b - 1.96 * F$se), max(F$b + 1.96 * F$se)); xt <- xr[2] + 0.04 * diff(xr)
hq <- RB[trait == "SRS" & analysis == A3 & grepl("^heterogeneity", subset), p]
pc <- ggplot(F, aes(b, lab, colour = sec == "SCQ")) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", linewidth = 0.9) +
  geom_point(aes(shape = p < 0.05), size = 3.6) + geom_text(aes(x = xt, label = ptxt, fontface = face(p)), hjust = 0, size = TS - 0.4, colour = "grey15") +
  facet_grid(sec ~ ., scales = "free_y", space = "free_y") + scale_colour_manual(values = c(`FALSE` = C_CAR, `TRUE` = "grey45"), guide = "none") + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  coord_cartesian(xlim = c(xr[1], xr[2] + 0.66 * diff(xr)), clip = "off") + scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), labels = c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2")) +
  labs(x = "CNV \u00d7 BMI \u00d7 autism score on Social Responsiveness Scale (95% CI)", y = NULL, tag = "c", title = "Three-way interaction, pooled and by subgroup",
       subtitle = sprintf("Filled: p < 0.05. Heterogeneity across the four CNVs %s.\nGrey: Social Communication Questionnaire (SCQ)", plab(hq))) +
  th_v2(BS) + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, size = 13.5), panel.spacing.y = unit(0.45, "lines"))
# ---------------------------------------------------------------- d
ST <- c("CNV carriers: G2MH + Simons VIP + SSC", "SSC unaffected relatives")
EST <- c("Observational (least squares)", "MR, three adiposity scores (2SLS)", "MR, childhood obesity score alone (Wald ratio)", "MR, childhood BMI score alone (Wald ratio)", "MR, adult BMI score alone (Wald ratio)")
ESTL <- c("Observational\n(least squares)", "Mendelian randomisation,\nthree scores together", "MR, childhood obesity\nscore alone", "MR, childhood BMI\nscore alone", "MR, adult BMI\nscore alone")
M <- MR[stratum %in% ST & estimator %in% EST & !(trait == "SCQ" & stratum == ST[2])]
M[, fac := factor(paste0(fifelse(stratum == ST[1], "Carriers, ", "Relatives, "), trait), c("Carriers, SRS", "Carriers, SCQ", "Relatives, SRS"))][, est := factor(estimator, rev(EST), rev(ESTL))]
tx <- M[estimator == EST[2], .(fac, n, lab = sprintf("n = %d\ninstrument strength\nF = %.1f\n%s,\n%s\ndetectable effect\n(80%% power) = %.2f", n, first_stage_F, ifelse(sargan_p < 0.05, "instruments disagree", "instruments do not disagree"), plab(sargan_p), mde80))]
XLd <- c(-0.45, 1.0)
pd <- ggplot(M, aes(b, est, colour = estimator == EST[1])) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_errorbar(aes(xmin = pmax(ci_lo, XLd[1]), xmax = pmin(ci_hi, XLd[2])), width = 0, orientation = "y", linewidth = 0.9) + geom_point(size = 3.6) +
  geom_errorbar(data = M[estimator == EST[2]], aes(xmin = ar_lo, xmax = ar_hi), width = 0.3, orientation = "y", linewidth = 0.5, colour = "grey20", position = position_nudge(y = -0.28)) +
  geom_text(data = tx, aes(x = XLd[1], y = 0.35, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = TS - 0.5, lineheight = 0.95, colour = "grey15") +
  facet_wrap(~fac, nrow = 1) + scale_colour_manual(values = c(`TRUE` = "grey10", `FALSE` = C_HI), guide = "none") + scale_y_discrete(expand = expansion(add = c(3.6, 0.6))) + coord_cartesian(xlim = XLd) + scale_x_continuous(breaks = c(-0.4, 0, 0.4, 0.8)) +
  labs(x = "Effect of BMI on the trait (SD per SD of BMI, 95% CI)", y = NULL, tag = "d", title = "Mendelian randomisation is inconclusive in carriers",
       subtitle = "Adiposity polygenic scores as instruments for BMI.\nThin bracket: weak-instrument-robust (Anderson-Rubin) interval") +
  th_v2(BS) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 15), panel.spacing.x = unit(1.2, "lines"))
fig <- wrap_elements(full = pa | pb) / wrap_elements(full = (pc | pd) + plot_layout(widths = c(0.85, 1))) + plot_layout(heights = c(1, 1.2))
save_v2(fig, "fig_cnv_bmi_pgs_main", 20, 15)
fwrite(rbind(sl[, .(panel = "a", series = paste("carriers, autism score", fifelse(half == "lo", "below", "above"), "median"), x = NA_real_, n, estimate = slope, se, p)], data.table(panel = "a", series = "SSC unaffected relatives", x = NA_real_, n = nrow(rel), estimate = coef(mr0)[["BMI"]], se = NA_real_, p = NA_real_),
             fq[, .(panel = "b", series = paste(group, "fifth", bin), x = score_mean, n, estimate = slope, se, p = NA_real_)], F[, .(panel = "c", series = paste(trait, gsub("\n", " ", sec), lab), x = NA_real_, n = n_carriers, estimate = b, se, p)],
             M[, .(panel = "d", series = paste(fac, estimator), x = NA_real_, n, estimate = b, se, p)]), file.path(V2OUT, "fig_cnv_bmi_pgs_main_plotted.csv"))
cat("86b done\n")
