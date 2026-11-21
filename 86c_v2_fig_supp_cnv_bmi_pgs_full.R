# 86c_v2_fig_supp_cnv_bmi_pgs_full.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 3: SUPPLEMENT)
# What it does: supplementary figure with everything behind the main polygenic-score figure, from the old figures
#   (66_fig_supp_interactions_main.R and 49_fig_supp_bmi_srs_scq_pgs.R; statistics unchanged, re-run and verified identical) plus the
#   new analyses of 86_v2_cnv_bmi_pgs_stats.R:
#   a  SRS against BMI, carriers and SSC unaffected relatives (two-way CNV x BMI; partial residuals) [66 a]
#   b  two-way and three-way interaction estimates for SRS and SCQ, pooled and by sex, CNV, cohort, autism diagnosis [66 c]
#   c  BMI slope in G2MH carriers after adjustment for polygenic scores (growth, adiposity, neurodevelopmental, all ten) [49 a]
#   d  each score's association with the trait next to the value expected if it acted only through BMI [49 b]
#   e  BMI x score interaction for each of the ten scores, within G2MH carriers [49 d]
#   f  Mendelian randomisation in full: six strata x two traits; 2SLS, Wald ratios, leave-one-instrument-out, Anderson-Rubin sets
#   g  NEW: autism diagnosis in carriers, BMI x autism score odds ratio (logistic), pooled, per CNV, per cohort
# Inputs : caches written by v2_base(), 86a and 86. Outputs: autism/fig_supp_cnv_bmi_pgs_full.{png,pdf} and _plotted.csv.
# null model: none here (tests are in 49 / 62 / 66 and 86_v2_cnv_bmi_pgs_stats.R).
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
v2_base(); CD <- file.path("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data", "v2_autism_cache")
load(file.path(CD, "pgs_stats_2026-11-20.RData")); load(file.path(CD, "pgs49_2026-11-20.RData"))
BS <- 14; OC <- c(SRS = C_CAR, SCQ = "#D9822B"); pd <- position_dodge(width = 0.6)
thS <- th_v2(BS) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12.5), legend.position = "top", legend.text = element_text(size = 12.5))
# ---------------------------------------------------------------- a (66 a)
C2 <- C[pheno == "SRS"]; z <- rbind(C2[is_cnv == 1], RELf(C2))
m <- lm(y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI, z); Xm <- model.matrix(m); nu <- grep("poly\\(age|male|^cg", colnames(Xm)); cf <- coef(m); cf[is.na(cf)] <- 0
ct <- as.vector(Xm[, nu, drop = FALSE] %*% cf[nu]); z[, ya := y - ct + ave(ct, is_cnv)][, grp := factor(fifelse(is_cnv == 1, "CNV carriers", "SSC unaffected relatives"), c("CNV carriers", "SSC unaffected relatives"))]
fl <- z[, { f <- lm(ya ~ BMI); g <- seq(quantile(BMI, 0.02), quantile(BMI, 0.98), length.out = 40); pr <- predict(f, data.frame(BMI = g), se.fit = TRUE); .(BMI = g, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }, by = grp]
T2 <- "Two-way: CNV x BMI"; T3 <- "Three-way: CNV x BMI x autism score"; r2 <- R66[trait == "SRS" & section == "Pooled" & type == T2]
GC <- c("CNV carriers" = C_CAR, "SSC unaffected relatives" = C_NON)
pa <- ggplot(z, aes(BMI, ya, colour = grp)) + geom_point(data = z[is_cnv == 0], size = 0.5, alpha = 0.08) + geom_point(data = z[is_cnv == 1], size = 1, alpha = 0.45) +
  geom_ribbon(data = fl, aes(BMI, ymin = lo, ymax = hi, fill = grp), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = fl, aes(BMI, fit), linewidth = 1.2) +
  scale_colour_manual(values = GC, name = NULL, labels = c("CNV carriers" = sprintf("CNV carriers (n = %d)", sum(z$is_cnv == 1)), "SSC unaffected relatives" = sprintf("Unaffected relatives (n = %d)", sum(z$is_cnv == 0)))) + scale_fill_manual(values = GC, guide = "none") + coord_cartesian(xlim = c(-3, 4.2), ylim = c(-2.6, 6.2)) +
  annotate("text", x = -3, y = 6.2, hjust = 0, vjust = 1, size = 4.3, fontface = face(r2$p), label = sprintf("CNV \u00d7 BMI = %+.3f, %s", r2$b, plab(r2$p))) + guides(colour = guide_legend(ncol = 1, override.aes = list(alpha = 1, size = 2.5))) +
  labs(x = "BMI (SD)", y = "Social Responsiveness Scale (SD)", tag = "a", title = "CNV \u00d7 BMI on SRS") + thS + theme(legend.position = "bottom")
# ---------------------------------------------------------------- b (66 c)
F <- copy(R66[!is.na(b)]); ordl <- c("All carriers", "Female", "Male", "Sex difference in carriers (male - female)", LOC4, "G2MH", "Simons VIP", "Simons Searchlight", "No autism diagnosis", "Autism diagnosis")
F[, label := factor(label, rev(ordl), sub("Sex difference in carriers (male - female)", "Male - female (within carriers)", rev(ordl), fixed = TRUE))][, section := factor(section, c("Pooled", "Sex", "CNV", "Cohort", "Autism diagnosis"), c("Pooled", "Sex", "CNV", "Cohort", "Autism\ndiagnosis"))]
F[, trait := factor(trait, c("SRS", "SCQ"))][, type := factor(type, c(T2, T3), c("Two-way: CNV \u00d7 BMI", "Three-way: CNV \u00d7 BMI \u00d7 autism score"))]
pb <- ggplot(F, aes(b, label, colour = trait, group = trait)) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd, linewidth = 0.7) +
  geom_point(aes(shape = p < 0.05), size = 2.6, position = pd) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("p \u2265 0.05", "p < 0.05"), name = NULL) +
  geom_text(aes(x = Inf, label = paste0("n = ", n)), position = pd, hjust = 1.05, size = 3.3, show.legend = FALSE) + facet_grid(section ~ type, scales = "free", space = "free_y") + scale_colour_manual(values = OC, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) + labs(x = "Interaction estimate (95% CI); n = carriers", y = NULL, tag = "b", title = "Interaction estimates by stratum (reference: SSC unaffected relatives, same sex for the sex strata)") +
  thS + theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(0.45, "lines"))
# ---------------------------------------------------------------- c, d, e (49 a, b, d)
A <- copy(A49)[, set := factor(set, levels = rev(unique(set)))][, outcome := factor(outcome, c("SRS", "SCQ"))]
pc <- ggplot(A, aes(b, set, colour = outcome, group = outcome)) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd, linewidth = 0.7) + geom_point(size = 2.6, position = pd) +
  scale_colour_manual(values = OC, name = NULL) + labs(x = "BMI slope (SD of trait per SD of BMI, 95% CI)", y = "Adjusted for", tag = "c", title = "BMI slope under score adjustment",
  subtitle = sprintf("G2MH carriers (SRS n = %d, SCQ n = %d);\nlast two rows: carriers with measured height", A[outcome == "SRS"]$n[1], A[outcome == "SCQ"]$n[1])) + thS
B <- copy(B49)[, lab := factor(LAB[score], levels = rev(LAB))][, outcome := factor(outcome, c("SRS", "SCQ"))]
pdd <- ggplot(B, aes(b, lab)) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se, colour = outcome), width = 0, orientation = "y", linewidth = 0.7) + geom_point(aes(colour = outcome), size = 2.6) +
  geom_point(aes(expected, lab), shape = 4, size = 2.8, stroke = 1) + facet_wrap(~outcome, nrow = 1) + scale_colour_manual(values = OC, guide = "none") +
  labs(x = "SD of trait per SD of polygenic score (95% CI)", y = NULL, tag = "d", title = "Polygenic scores and autistic traits (G2MH carriers)", subtitle = "Dot: observed. Cross: expected if the score acted on the trait only through BMI") + thS
D <- copy(D49)[!grepl("^half_|_adj$", score)][, lab := factor(LAB[score], levels = rev(LAB))][, outcome := factor(outcome, c("SRS", "SCQ"))]
pe <- ggplot(D, aes(b, lab, colour = outcome)) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", linewidth = 0.7) + geom_point(aes(shape = p < 0.05), size = 2.6) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") + facet_wrap(~outcome, nrow = 1) + scale_colour_manual(values = OC, guide = "none") +
  labs(x = "BMI \u00d7 score interaction (95% CI)", y = NULL, tag = "e", title = "BMI \u00d7 score, ten scores (G2MH carriers)", subtitle = "Filled: p < 0.05 (no correction for the ten scores)") + thS
# ---------------------------------------------------------------- f (Mendelian randomisation in full)
STL <- c("CNV carriers: G2MH + Simons VIP + SSC" = "Carriers:\nall genotyped", "CNV carriers: G2MH" = "Carriers:\nG2MH", "CNV carriers: Simons VIP" = "Carriers:\nSimons VIP", "CNV carriers: G2MH (model of 49)" = "Carriers: G2MH,\nseven score covariates",
         "SSC unaffected relatives" = "Unaffected\nrelatives", "SSC autistic probands, no CNV" = "Autistic,\nno CNV")
M <- copy(MR)[, st := factor(STL[stratum], STL)][, trait := factor(trait, c("SRS", "SCQ"))]
el <- unique(M$estimator); M[, est := factor(estimator, rev(el), rev(sub("^MR, ", "", sub(" \\(Wald ratio\\)| \\(2SLS\\)", "", sub("Observational \\(least squares\\)", "observational", el)))))]
XLf <- c(-1, 1.5); tx <- M[grepl("^MR, three", estimator), .(st, trait, lab = sprintf("n = %d; F = %.1f\ninstruments disagree %s", n, first_stage_F, plab(sargan_p)))]
M[, ar_empty := grepl("^empty", ar_note)]
pf <- ggplot(M, aes(b, est, colour = estimator == el[1])) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_errorbar(aes(xmin = pmax(ar_lo, XLf[1]), xmax = pmin(ar_hi, XLf[2])), width = 0.35, orientation = "y", linewidth = 0.45, colour = "grey25", position = position_nudge(y = -0.3), na.rm = TRUE) +
  geom_text(data = M[ar_empty == TRUE], aes(x = XLf[2], label = "robust set empty"), hjust = 1, vjust = 1.6, size = 2.9, colour = "grey25") +
  geom_errorbar(aes(xmin = pmax(ci_lo, XLf[1]), xmax = pmin(ci_hi, XLf[2])), width = 0, orientation = "y", linewidth = 0.75) + geom_point(size = 2.5) +
  geom_text(data = tx, aes(x = XLf[1], y = 0.15, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2, lineheight = 0.95) +
  facet_grid(trait ~ st) + scale_colour_manual(values = c(`TRUE` = "grey10", `FALSE` = C_HI), guide = "none") + scale_y_discrete(expand = expansion(add = c(1.9, 0.6))) + coord_cartesian(xlim = XLf) +
  labs(x = "Effect of BMI on the trait (SD per SD of BMI; 95% CI, clipped at -1 and 1.5); thin bracket below each estimate: Anderson-Rubin weak-instrument-robust 95% set", y = NULL, tag = "f",
       title = "Mendelian randomisation of BMI on autistic traits: all strata, leave-one-instrument-out and weak-instrument-robust sets",
       subtitle = "Instruments: polygenic scores for childhood obesity, childhood BMI and adult BMI. F = joint first-stage F; 'instruments disagree' = Sargan over-identification test") + thS + theme(strip.text.y = element_text(angle = 0))
# ---------------------------------------------------------------- g (diagnosis)
G <- DX[term == "BMI x autism score" & !is.na(logOR) & !grepl("leave out", group)]
G[, smp := fifelse(grepl("^carriers of the SRS", sample), "SRS / SCQ sample (BMI at the questionnaire)", "Growth table (one visit per person)")][, lab := sub("All four CNVs, cohort: ", "", group)][, lab := sub("^SVIP$", "Simons VIP", lab)]
G[, sec := factor(fcase(group == "All four CNVs", "Pooled", group %in% LOC4, "CNV", default = "Cohort"), c("Pooled", "CNV", "Cohort"))][, lab := factor(sprintf("%s", lab), rev(c("All four CNVs", LOC4, "G2MH", "Simons VIP", "Cardiff")))]
SC <- c("SRS / SCQ sample (BMI at the questionnaire)" = C_CAR, "Growth table (one visit per person)" = "#7A7A7A"); hq <- DX[grepl("^heterogeneity", group) & grepl("^carriers of the SRS", sample), p]
pg <- ggplot(G, aes(OR, lab, colour = smp, group = smp)) + geom_vline(xintercept = 1, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = OR_lo, xmax = pmin(OR_hi, 6)), width = 0, orientation = "y", position = pd, linewidth = 0.7) +
  geom_point(aes(shape = p < 0.05), size = 2.6, position = pd) + geom_text(aes(x = 9.5, label = sprintf("%d / %d", n_asd, n)), position = pd, hjust = 1, size = 3.3, show.legend = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") + scale_colour_manual(values = SC, name = NULL) + scale_x_log10(breaks = c(0.5, 1, 2, 4)) + coord_cartesian(xlim = c(0.45, 9.5)) + facet_grid(sec ~ ., scales = "free_y", space = "free_y") +
  guides(colour = guide_legend(ncol = 1)) + labs(x = "Odds ratio for BMI \u00d7 autism score (95% CI); right: diagnosed / all", y = NULL, tag = "g", title = "Autism diagnosis in carriers: BMI \u00d7 autism score",
       subtitle = sprintf("Logistic regression; heterogeneity across CNVs %s", plab(hq))) + thS + theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")
fig <- wrap_elements(full = (pa | pb) + plot_layout(widths = c(0.42, 1))) / wrap_elements(full = (pc | pdd | pe) + plot_layout(widths = c(0.8, 1, 1))) / wrap_elements(full = pf) / wrap_elements(full = (pg | plot_spacer()) + plot_layout(widths = c(1, 0.9))) + plot_layout(heights = c(1.25, 0.8, 1.25, 0.85))
save_v2(fig, "fig_supp_cnv_bmi_pgs_full", 20, 32)
fwrite(rbind(F[, .(panel = "b", outcome = trait, term = paste(section, label, type), n, estimate = b, se, p)], A[, .(panel = "c", outcome, term = set, n, estimate = b, se, p)], B[, .(panel = "d", outcome, term = lab, n, estimate = b, se, p)],
             D[, .(panel = "e", outcome, term = lab, n, estimate = b, se, p)], M[, .(panel = "f", outcome = trait, term = paste(stratum, estimator), n, estimate = b, se, p)], G[, .(panel = "g", outcome = "autism diagnosis (log OR)", term = paste(smp, group), n, estimate = logOR, se, p)], fill = TRUE),
       file.path(V2OUT, "fig_supp_cnv_bmi_pgs_full_plotted.csv"))
cat("86c done\n")
