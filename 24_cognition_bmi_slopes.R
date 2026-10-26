# BMI effects on social-communication traits (SRS, SCQ), continuous BMI, no age dimension (JS / Claude 2026-10-24).
# Follows the main-figure conventions (rows: non-carriers, 16p del, 16p dup, 22q del, 22q dup; summary panel at the
# bottom; BMI-tail colours) but plots the trait against BMI Z instead of age, because the carrier effects are not
# age-dependent and the 70/30 tails are sparse. Non-carrier row shows unaffected SSC relatives AND SSC non-carrier
# probands (ASD), to test that the BMI-trait association is CNV-specific rather than an ASD phenomenon.
# Trait values are age-adjusted within group (age smooth removed). Fits: gam y ~ s(age) + BMI (slope = SD/SD);
# carriers vs unaffected: bam y ~ cnv + s(age, by = cnv) + s(age, by = BMI) + cnv:BMI (as in 23_unaffected_controls.R).
# Carrier slopes also for non-ASD carriers and with ASD as a covariate. Data build from 21_pheno_bmi_gam.R.
# Outputs: cognition_bmi_slopes.csv, fig_cognition_bmi_slopes.{png,pdf}
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
if (!exists("D")) { sink(tempfile()); source(file.path(HERE, "21_pheno_bmi_gam.R")); sink() }   # D already present when sourced after 23
suppressPackageStartupMessages({ library(mgcv); library(ggplot2); library(patchwork) })
t0 <- Sys.time()
PH <- c("SRS", "SCQ", "CNB accuracy", "CNB speed")   # JS 2026-10-24: CNB added (carriers only; the CNB non-carriers have no BMI)
PH_FIG <- c("SRS", "SCQ")                               # standalone figure of this script
LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
X <- D[pheno %in% PH & !is.na(BMI) & !is.na(y)]
X[, grp := fcase(locus != "No CNV", locus,
                 cohort == "SSC" & role == "proband", "SSC ASD probands",
                 cohort == "SSC" & role %in% c("sibling", "parent"), "Unaffected relatives",
                 default = NA_character_)]
X <- X[!is.na(grp)][, asd_f := factor(asd)]
X[, tail3 := fifelse(BMI >= qnorm(.7), "BMI \u2265 70th centile", fifelse(BMI <= qnorm(.3), "BMI \u2264 30th centile", "BMI 30th-70th centile"))]
# age-adjust within group x trait (remove the centred age smooth)
X[, y_adj := { m <- gam(y ~ s(age, k = 8) + BMI, method = "REML"); y - predict(m, type = "terms")[, "s(age)"] }, by = .(pheno, grp)]
cat("== sample\n"); print(dcast(X[, .N, by = .(pheno, grp)], grp ~ pheno, value.var = "N", fill = 0))

# ---------------------------------------------------------------- slopes
slope <- function(x, extra = "") tryCatch({
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI", extra)), data = x, method = "REML"); co <- summary(m)$p.table["BMI", ]
  as.data.table(list(n = nrow(x), slope = unname(co[1]), se = unname(co[2]), p = unname(co[4]))) }, error = function(e) as.data.table(list(n = nrow(x), slope = NA_real_, se = NA_real_, p = NA_real_)))
vs_unaff <- function(ctl, xc, extra = "") tryCatch({
  z <- rbind(ctl, xc)[, is_cnv := factor(ifelse(grp == "Unaffected relatives", "ctl", "CNV"), levels = c("ctl", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  m <- bam(as.formula(paste("y ~ is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI", extra)), data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]
  as.data.table(list(diff = unname(co[1]), diff_se = unname(co[2]), diff_p = unname(co[4]))) }, error = function(e) as.data.table(list(diff = NA_real_, diff_se = NA_real_, diff_p = NA_real_)))
rows <- list()
for (ph in PH) {
  xx <- X[pheno == ph]; ctl <- xx[grp == "Unaffected relatives"]
  rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = "SSC siblings", variant = "reference"), slope(xx[grp == "Unaffected relatives" & role == "sibling"]))
  if (ph == "SRS") rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = "SSC parents", variant = "reference"), slope(xx[grp == "Unaffected relatives" & role == "parent"]))
  rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = "Unaffected relatives", variant = "reference"), slope(ctl))
  pr <- xx[grp == "SSC ASD probands"]
  rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = "SSC ASD probands", variant = "reference"), slope(pr), vs_unaff(ctl, pr))
  for (g in c(LOC, "All CNVs")) {
    xc <- if (g == "All CNVs") xx[grp %in% LOC] else xx[grp == g]; ex <- if (g == "All CNVs") "+ locus" else ""
    if (nrow(xc) < 15) next
    rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = g, variant = "all carriers"), slope(xc, ex), vs_unaff(ctl, xc, if (g == "All CNVs") "+ locus" else ""))
    xn <- xc[asd %in% 0L]
    if (nrow(xn) >= 15) rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = g, variant = "non-ASD carriers"), slope(xn, ex), vs_unaff(ctl, xn, if (g == "All CNVs") "+ locus" else ""))
    xa <- xc[!is.na(asd)]
    if (nrow(xa) >= 15 && uniqueN(xa$asd) == 2) rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = g, variant = "ASD-adjusted"), slope(xa, paste(ex, "+ asd_f")))
  }
}
R <- rbindlist(rows, fill = TRUE)
R[, `:=`(ci_lo = slope - 1.96 * se, ci_hi = slope + 1.96 * se)]
fwrite(R, file.path(HERE, "cognition_bmi_slopes.csv"))
cat("\n== BMI slopes (SD of trait per SD of BMI)\n"); print(R[, .(pheno, group, variant, n, slope = round(slope, 3), se = round(se, 3), p = signif(p, 2), vs_unaff = round(diff, 3), vs_unaff_p = signif(diff_p, 2))], nrows = 100)

# ---------------------------------------------------------------- figure
COLS3 <- c("BMI \u2265 70th centile" = "#F0645A", "BMI 30th-70th centile" = "#B8BEC4", "BMI \u2264 30th centile" = "#00A5AD")
REFC <- c("Unaffected relatives" = "#8C9A9A", "SSC ASD probands" = "#222222")
YL <- c(SRS = "SRS total (Z vs SSC unaffected)\nage-adjusted", SCQ = "SCQ lifetime (Z vs SSC siblings)\nage-adjusted", `CNB accuracy` = "CNB accuracy (Z vs G2MH non-carriers)\nage-adjusted", `CNB speed` = "CNB speed (Z, higher = faster)\nage-adjusted")
BR <- c(-Inf, -1, 0, 1, 2, Inf)
XL <- c(-3.2, 4.2)
bins <- function(x) x[, .(bmi = mean(BMI), m = mean(y_adj), se = sd(y_adj) / sqrt(.N), n = .N), by = .(bin = cut(BMI, BR))][n >= 5]
fitline <- function(x) { m <- lm(y_adj ~ BMI, data = x); nd <- data.frame(BMI = seq(max(XL[1], min(x$BMI)), min(XL[2], max(x$BMI)), length.out = 50)); cbind(nd, as.data.frame(predict(m, nd, se.fit = TRUE)[1:2])) }
ylim_of <- function(ph, car = FALSE) { x <- X[pheno == ph & grp != "Unaffected relatives" & (!car | grp %in% LOC)]; r <- quantile(x$y_adj, c(0.005, 0.995)); r + c(-0.08, 0.30) * diff(r) }
ptxt <- function(p) if (is.na(p)) "n/a" else pfx(p)
panel_ctl <- function(ph, probands = TRUE, show_title = TRUE) {   # probands = FALSE: unaffected relatives only (main-figure version)
  G <- if (probands) names(REFC) else "Unaffected relatives"
  x <- X[pheno == ph & grp %in% G]; yl <- ylim_of(ph, !probands)
  fl <- rbindlist(lapply(G, function(g) cbind(grp = g, fitline(x[grp == g]))))
  bb <- rbindlist(lapply(G, function(g) cbind(grp = g, bins(x[grp == g]))))
  r <- R[pheno == ph & variant == "reference"]
  ru <- r[group == "Unaffected relatives"]; rp <- r[group == "SSC ASD probands"]
  lab <- if (probands) sprintf("controls: slope = %+.2f\nASD (no CNV): slope = %+.2f", ru$slope, rp$slope) else gsub("= < ", "< ", sprintf("slope = %+.2f, p = %s", ru$slope, ptxt(ru$p)), fixed = TRUE)
  lab3 <- if (probands) gsub("= < ", "< ", sprintf("p-value vs controls = %s", ptxt(rp$diff_p)), fixed = TRUE) else ""
  bold3 <- probands && isTRUE(rp$diff_p < 0.05)
  sub <- if (probands) sprintf("Unaffected (n = %d)\nASD, no CNV (n = %d)", sum(x$grp == "Unaffected relatives"), sum(x$grp == "SSC ASD probands")) else sprintf("Unaffected (n = %d)", nrow(x))
  ggplot(x, aes(BMI, y_adj, colour = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 0.5, alpha = 0.12) + geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit, fill = grp), alpha = 0.25, colour = NA) +
    geom_line(data = fl, aes(BMI, fit), linewidth = 1) + geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), size = 0.45, linewidth = 0.6, shape = 18) +
    annotate("text", x = XL[2], y = yl[2], label = lab, hjust = 1, vjust = 1, size = 3.1) +
    annotate("text", x = XL[2], y = yl[2] - 0.15 * diff(yl), label = lab3, hjust = 1, vjust = 1, size = 3.1, fontface = if (bold3) "bold" else "plain") +
    scale_colour_manual(values = REFC, labels = c("Unaffected relatives" = "controls (siblings and parents)", "SSC ASD probands" = "ASD")) + scale_fill_manual(values = REFC, guide = "none") + coord_cartesian(xlim = XL, ylim = yl) +
    labs(x = "BMI (Z)", y = YL[[ph]], colour = NULL, title = if (show_title) ph else NULL, subtitle = sub) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 8), legend.position = if (probands) "inside" else "none", legend.position.inside = c(0.4, 0.28), legend.background = element_blank(), legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5, size = 16)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.4)))
}
panel_cnv <- function(ph, g) {
  x <- X[pheno == ph & grp == g]; yl <- ylim_of(ph, TRUE); x[, tail3 := factor(tail3, names(COLS3))]
  xu <- X[pheno == ph & grp == "Unaffected relatives"]; has_ctl <- nrow(xu) > 0
  fl <- fitline(x); bb <- bins(x); fu <- if (has_ctl) fitline(xu) else NULL
  r <- R[pheno == ph & group == g & variant == "all carriers"]; rn <- R[pheno == ph & group == g & variant == "non-ASD carriers"]
  # JS 2026-10-24: print the carrier slope and the p-value of its difference from controls only (the within-carrier p and the
  # underpowered non-ASD slope are in cognition_bmi_slopes.csv)
  lab1 <- if (nrow(r)) sub("= -0.00", "= +0.00", sprintf("slope = %+.2f", r$slope), fixed = TRUE) else "n too small"
  lab2 <- if (!nrow(r)) "" else if (has_ctl) gsub("= < ", "< ", sprintf("p-value vs controls = %s", ptxt(r$diff_p)), fixed = TRUE) else gsub("= < ", "< ", sprintf("p-value = %s", ptxt(r$p)), fixed = TRUE)   # no controls: within-carrier slope p
  lab3 <- ""
  bold1 <- FALSE; bold <- nrow(r) && isTRUE((if (has_ctl) r$diff_p else r$p) < 0.05)   # bold only if p < 0.05
  ggplot(x, aes(BMI, y_adj)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    { if (has_ctl) geom_line(data = fu, aes(BMI, fit), colour = REFC[["Unaffected relatives"]], linetype = 2, linewidth = 0.8) } +
    geom_point(aes(colour = tail3), size = 1.3, alpha = 0.65) +
    geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), fill = "grey30", alpha = 0.18) + geom_line(data = fl, aes(BMI, fit), colour = "grey10", linewidth = 1) +
    geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), colour = "grey10", size = 0.5, linewidth = 0.6, shape = 18) +
    annotate("text", x = XL[2], y = yl[2], label = lab1, hjust = 1, vjust = 1, size = 3.1, fontface = if (bold1) "bold" else "plain") +
    annotate("text", x = XL[2], y = yl[2] - 0.075 * diff(yl), label = lab2, hjust = 1, vjust = 1, size = 3.1, fontface = if (bold) "bold" else "plain") +
    annotate("text", x = XL[2], y = yl[2] - 0.15 * diff(yl), label = lab3, hjust = 1, vjust = 1, size = 2.8, colour = "grey30") +
    scale_colour_manual(values = COLS3, drop = FALSE) + coord_cartesian(xlim = XL, ylim = yl) +
    labs(x = "BMI (Z)", y = YL[[ph]], colour = NULL, subtitle = sprintf("%s (n = %d)", g, nrow(x))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 7), legend.position = "inside", legend.position.inside = c(0.3, 0.1), legend.background = element_blank(), legend.key.height = unit(7, "pt"))
}
panel_forest <- function(ph) {
  r <- R[pheno == ph & !is.na(slope)]
  ord <- c("SSC siblings", "SSC parents", "SSC ASD probands", LOC, "All CNVs")
  r[, group := factor(group, levels = rev(ord))][, variant := factor(variant, levels = c("reference", "all carriers", "non-ASD carriers", "ASD-adjusted"))]
  r <- r[!is.na(group)]
  r[, star := fifelse(!is.na(diff_p) & diff_p < 0.05 & variant != "ASD-adjusted", "*", "")]
  SH <- c(reference = 15, `all carriers` = 16, `non-ASD carriers` = 1, `ASD-adjusted` = 17); CV <- c(reference = "#8C9A9A", `all carriers` = "grey10", `non-ASD carriers` = "grey10", `ASD-adjusted` = "#4B6EA9")
  ggplot(r, aes(slope, group, colour = variant, shape = variant)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_hline(yintercept = uniqueN(r[variant != "reference", group]) + 0.5, colour = "grey80") +
    geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.65), size = 0.5, linewidth = 0.6) +
    geom_text(aes(x = ci_hi, label = star), position = position_dodge(width = 0.65), hjust = -0.6, size = 5, show.legend = FALSE) +
    scale_shape_manual(values = SH) + scale_colour_manual(values = CV) +
    labs(x = sprintf("BMI slope on %s (SD per SD of BMI, 95%% CI)", ph), y = NULL, colour = NULL, shape = NULL, subtitle = "BMI slope by group (* vs unaffected p < 0.05)") +
    theme_classic(base_size = 13) + theme(legend.position = "inside", legend.position.inside = c(0.85, 0.86), legend.background = element_blank(), legend.text = element_text(size = 8), legend.key.height = unit(9, "pt"))
}
hdr <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 8, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
if (toupper(Sys.getenv("COG_DRAW_FIG", "TRUE")) != "FALSE") {   # 25_main_figure_v3.R sources this file for its panels and skips the figure
plots <- unlist(lapply(PH_FIG, function(ph) c(list(panel_ctl(ph)), lapply(LOC, function(g) panel_cnv(ph, g)), list(panel_forest(ph)))), recursive = FALSE)
fig <- wrap_plots(plots, ncol = 2, byrow = FALSE, heights = c(1, 1, 1, 1, 1, 1.05)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(HERE, "fig_cognition_bmi_slopes.png"), fig, width = 11.5, height = 25, dpi = 130)
tryCatch(ggsave(file.path(HERE, "fig_cognition_bmi_slopes.pdf"), fig, width = 11.5, height = 25, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_cognition_bmi_slopes.pdf"), fig, width = 11.5, height = 25))
}
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
