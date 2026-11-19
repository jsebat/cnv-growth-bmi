# SRS / SCQ panels of the main figure on the Box-Cox scale (JS / Claude 2026-11-14; Jon: Box-Cox is the scaling Peter Visscher
# recommended, report P both ways per CNV, add the pooled result to the figure).
#  - Raw totals are Box-Cox transformed (raw + 1; lambda by profile likelihood in the unaffected relatives of the same form,
#    raw + 1 ~ age + age^2 + sex, as in 50_srs_scq_scale_sensitivity.R) and z-scored against those relatives; this replaces the
#    linear z of 21_pheno_bmi_gam.R for SRS and SCQ only. Box-Cox of the linear Z is the same transform up to the shift constant.
#  - 24_cognition_bmi_slopes.R is then rerun unchanged on that scale (slopes, carrier vs relatives contrasts, panels).
#  - Ordinal check, invariant to any monotone rescaling: cumulative-logit model of the raw total (<= 40 ordered levels), carrier vs
#    relatives difference in the BMI effect (cnvBMI), per SRS form and inverse-variance combined (functions of script 50).
#  - Panels: each CNV prints its slope and "p vs controls: Box-Cox = ..; ordinal = .."; the control panel adds the four CNVs
#    pooled (locus-adjusted carrier fit line + slope + both P values).
# Output: cog_panels_boxcox.RData (no identifiers; for 35b_main_figure_v6_boxcox.R on Expanse), cognition_bmi_slopes_boxcox.csv,
#         fig_cognition_bmi_slopes_boxcox.png (local check).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))
suppressPackageStartupMessages(library(ordinal))
# ---------------------------------------------------------------- Box-Cox scale for SRS and SCQ in D
rs <- list(child = ssc[form == "child" & role == "sibling", srs], adult = ssc[form == "adult" & role == "parent", srs])
D[pheno == "SRS", form2 := fifelse(abs(y - zref(raw, rs$child)) < 1e-8, "child", fifelse(abs(y - zref(raw, rs$adult)) < 1e-8, "adult", NA_character_))]
D[pheno == "SCQ", form2 := "scq"]
stopifnot(!anyNA(D[pheno %in% c("SRS", "SCQ"), form2]))
D[, y_lin := y]
LAM <- list()
for (f in c("child", "adult", "scq")) {
  ph <- if (f == "scq") "SCQ" else "SRS"
  r <- D[pheno == ph & form2 == f & unaffected == TRUE & locus == "No CNV"]; r[, `:=`(yy = raw + 1, a2 = age^2)]
  b <- MASS::boxcox(yy ~ age + a2 + sex, data = r, lambda = seq(-2, 2, 0.01), plotit = FALSE); L <- b$x[which.max(b$y)]; LAM[[f]] <- L
  bc <- function(x) if (abs(L) < 1e-8) log(x + 1) else ((x + 1)^L - 1) / L
  D[pheno == ph & form2 == f, y := zref(bc(raw), bc(r$raw))]
}
cat("Box-Cox lambda:", paste(names(LAM), round(unlist(LAM), 2), collapse = "; "), "\n")
# ---------------------------------------------------------------- 24 on the Box-Cox scale
source(file.path(HERE, "24_cognition_bmi_slopes.R"))
X[, form2 := D$form2[match(paste(pheno, cohort, sex, age, raw, BMI), paste(D$pheno, D$cohort, D$sex, D$age, D$raw, D$BMI))]]
stopifnot(!anyNA(X[pheno %in% c("SRS", "SCQ"), form2]))
# ---------------------------------------------------------------- ordinal P (carriers vs unaffected relatives), per CNV and pooled
ord1 <- function(z) tryCatch({
  z <- copy(z); br <- unique(quantile(z$raw, seq(0, 1, length.out = 41), type = 1)); z[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)]
  z[, is_cnv := as.numeric(grp != "Unaffected relatives")][, cnvBMI := is_cnv * BMI]
  cl <- if (uniqueN(z$grp[z$grp != "Unaffected relatives"]) > 1) "grp + " else ""
  m <- clm(as.formula(paste("ry ~", cl, "is_cnv * poly(age, 2) + BMI + cnvBMI")), data = z); co <- coef(summary(m)); c(co["cnvBMI", 1], co["cnvBMI", 2])
}, error = function(e) { message("clm: ", conditionMessage(e)); rep(NA_real_, 2) })
ord_p <- function(z) { r <- t(sapply(split(z, z$form2), function(d) if (sum(d$grp != "Unaffected relatives") >= 15) ord1(d) else rep(NA_real_, 2)))
  r <- r[!is.na(r[, 1]) & !is.na(r[, 2]), , drop = FALSE]; if (!nrow(r)) return(c(NA, NA, NA)); w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); se <- sqrt(1 / sum(w)); c(b, se, 2 * pnorm(-abs(b / se))) }
R[, `:=`(ord_logodds = NA_real_, ord_se = NA_real_, ord_p = NA_real_)]
for (ph in c("SRS", "SCQ")) for (g in c(LOC, "All CNVs")) {
  xc <- X[pheno == ph & (if (g == "All CNVs") grp %in% LOC else grp == g)]; ctl <- X[pheno == ph & grp == "Unaffected relatives"]
  if (nrow(xc) < 15) next
  o <- ord_p(rbind(xc, ctl)); R[pheno == ph & group == g & variant == "all carriers", `:=`(ord_logodds = o[1], ord_se = o[2], ord_p = o[3])]
}
fwrite(R, file.path(HERE, "cognition_bmi_slopes_boxcox.csv"))
cat("\n== Box-Cox scale: carrier slopes, vs relatives (Box-Cox and ordinal)\n")
print(R[pheno %in% c("SRS", "SCQ") & variant %in% c("reference", "all carriers"), .(pheno, group, n, slope = round(slope, 3), p = signif(p, 2), vs_rel = round(diff, 3), p_bc = signif(diff_p, 2), ord = round(ord_logodds, 3), p_ord = signif(ord_p, 2))], nrows = 50)
# ---------------------------------------------------------------- panels (Box-Cox labels, both P values, pooled result)
YL <- c(SRS = "SRS total (Box-Cox Z vs SSC unaffected)\nage-adjusted", SCQ = "SCQ lifetime (Box-Cox Z vs SSC siblings)\nage-adjusted", `CNB accuracy` = YL[["CNB accuracy"]], `CNB speed` = YL[["CNB speed"]])
pp2 <- function(pb, po) gsub("= < ", "< ", sprintf("vs controls: Box-Cox p = %s\nordinal p = %s", ptxt(pb), ptxt(po)), fixed = TRUE)
POOLC <- "#7B241C"
pooled_line <- function(ph) { x4 <- X[pheno == ph & grp %in% LOC]; m <- lm(y_adj ~ BMI + grp, data = x4)
  g <- seq(max(XL[1], quantile(x4$BMI, 0.02)), min(XL[2], quantile(x4$BMI, 0.98)), length.out = 50)
  data.frame(BMI = g, fit = rowMeans(sapply(LOC, function(L) predict(m, data.frame(BMI = g, grp = L))))) }
panel_cnv <- function(ph, g) {
  x <- X[pheno == ph & grp == g]; yl <- ylim_of(ph, TRUE); x[, tail3 := factor(tail3, names(COLS3))]
  xu <- X[pheno == ph & grp == "Unaffected relatives"]; has_ctl <- nrow(xu) > 0
  fl <- fitline(x); bb <- bins(x); fu <- if (has_ctl) fitline(xu) else NULL
  r <- R[pheno == ph & group == g & variant == "all carriers"]
  lab1 <- if (nrow(r)) sub("= -0.00", "= +0.00", sprintf("slope = %+.2f", r$slope), fixed = TRUE) else "n too small"
  lab2 <- if (!nrow(r)) "" else if (has_ctl) pp2(r$diff_p, r$ord_p) else gsub("= < ", "< ", sprintf("p-value = %s", ptxt(r$p)), fixed = TRUE)
  bold <- nrow(r) && isTRUE((if (has_ctl) r$diff_p else r$p) < 0.05)
  ggplot(x, aes(BMI, y_adj)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    { if (has_ctl) geom_line(data = fu, aes(BMI, fit), colour = REFC[["Unaffected relatives"]], linetype = 2, linewidth = 0.8) } +
    geom_point(aes(colour = tail3), size = 1.3, alpha = 0.65) +
    geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), fill = "grey30", alpha = 0.18) + geom_line(data = fl, aes(BMI, fit), colour = "grey10", linewidth = 1) +
    geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), colour = "grey10", size = 0.5, linewidth = 0.6, shape = 18) +
    annotate("text", x = XL[2], y = yl[2], label = lab1, hjust = 1, vjust = 1, size = 3.1) +
    annotate("text", x = XL[2], y = yl[2] - 0.075 * diff(yl), label = lab2, hjust = 1, vjust = 1, size = 2.9, lineheight = 0.9, fontface = if (bold) "bold" else "plain") +
    scale_colour_manual(values = COLS3, drop = FALSE) + coord_cartesian(xlim = XL, ylim = yl) +
    labs(x = "BMI (Z)", y = YL[[ph]], colour = NULL, subtitle = sprintf("%s (n = %d)", g, nrow(x))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 7), legend.position = "inside", legend.position.inside = c(0.3, 0.1), legend.background = element_blank(), legend.key.height = unit(7, "pt"))
}
panel_ctl <- function(ph, probands = TRUE, show_title = TRUE) {
  G <- if (probands) names(REFC) else "Unaffected relatives"
  x <- X[pheno == ph & grp %in% G]; yl <- ylim_of(ph, !probands)
  fl <- rbindlist(lapply(G, function(g) cbind(grp = g, fitline(x[grp == g])))); bb <- rbindlist(lapply(G, function(g) cbind(grp = g, bins(x[grp == g]))))
  r <- R[pheno == ph & variant == "reference"]; ru <- r[group == "Unaffected relatives"]; rp <- r[group == "SSC ASD probands"]
  ra <- R[pheno == ph & group == "All CNVs" & variant == "all carriers"]; pl <- pooled_line(ph)
  lab <- sprintf("controls: slope = %+.2f\nASD (no CNV): slope = %+.2f\nvs controls p = %s", ru$slope, rp$slope, ptxt(rp$diff_p))
  labp <- sprintf("4 CNVs pooled (red): slope = %+.2f\n%s", ra$slope, pp2(ra$diff_p, ra$ord_p))
  boldp <- isTRUE(ra$diff_p < 0.05)
  ggplot(x, aes(BMI, y_adj, colour = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 0.5, alpha = 0.12) + geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit, fill = grp), alpha = 0.25, colour = NA) +
    geom_line(data = fl, aes(BMI, fit), linewidth = 1) + geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), size = 0.45, linewidth = 0.6, shape = 18) +
    geom_line(data = pl, aes(BMI, fit), inherit.aes = FALSE, colour = POOLC, linewidth = 1.3) +
    annotate("text", x = XL[2], y = yl[2], label = lab, hjust = 1, vjust = 1, size = 2.8, lineheight = 0.9) +
    annotate("text", x = XL[2], y = yl[2] - 0.17 * diff(yl), label = labp, hjust = 1, vjust = 1, size = 2.8, lineheight = 0.9, colour = POOLC, fontface = if (boldp) "bold" else "plain") +
    scale_colour_manual(values = REFC, labels = c("Unaffected relatives" = "controls (siblings and parents)", "SSC ASD probands" = "ASD")) + scale_fill_manual(values = REFC, guide = "none") + coord_cartesian(xlim = XL, ylim = yl) +
    labs(x = "BMI (Z)", y = YL[[ph]], colour = NULL, title = if (show_title) ph else NULL, subtitle = sprintf("Unaffected (n = %d)\nASD, no CNV (n = %d)", sum(x$grp == "Unaffected relatives"), sum(x$grp == "SSC ASD probands"))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.position.inside = c(0.4, 0.2), legend.background = element_blank(), legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5, size = 16)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.4)))
}
# local check: the two columns as they will appear in the main figure
plots <- unlist(lapply(c("SRS", "SCQ"), function(ph) c(list(panel_ctl(ph, probands = TRUE, show_title = TRUE)), lapply(LOC, function(g) panel_cnv(ph, g)))), recursive = FALSE)
ggsave(file.path(HERE, "fig_cognition_bmi_slopes_boxcox.png"), wrap_plots(plots, ncol = 2, byrow = FALSE), width = 9, height = 20, dpi = 130, bg = "white")
X2 <- copy(X); if ("id" %in% names(X2)) X2[, id := NULL]; X <- X2
objs <- c("X", "R", "panel_ctl", "panel_cnv", "fitline", "bins", "ylim_of", "ptxt", "REFC", "COLS3", "XL", "YL", "pp2", "POOLC", "pooled_line", "LAM")
save(list = objs, file = file.path(HERE, "cog_panels_boxcox.RData"))
cat("cog_panels_boxcox.RData written\n")
