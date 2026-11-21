# 84_v2_block_srs_scq.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 1)
# What it does: block for the head / autistic-traits figure: SRS (tags j-n) and SCQ (tags o-s) against BMI Z. First panel of each
#   row = SSC unaffected relatives, SSC autistic probands without a CNV, and the four CNVs pooled; then one panel per CNV, with the
#   carrier BMI slope and the CNV x BMI P (Box-Cox scale and ordinal check). Statistics are those of 51b_cog_panels_boxcox_combined.R,
#   re-run from its own code on the current data by v2_base() (83_v2_autism_common.R); nothing is copied from the old CSV, and the
#   re-run is compared with the old CSV (columns prev_* / moved).
# Inputs : DATA/combined_srs_scq.rds (G2MH + Simons VIP + Simons Searchlight + SSC carriers), SSC relatives and probands (21 / 24).
# Outputs: run_v2_2026-11-20/autism/block_srs_scq.{png,pdf}, block_srs_scq.csv (every printed slope and test),
#          block_srs_scq_binned_means.csv (plotted diamonds), block_srs_scq_n_by_cohort.csv.
# Models (51b / 24): carrier slope = gam y ~ s(age, k = 8) + BMI (+ cohort x CNV); CNV x BMI = bam y ~ carrier + s(age, by = carrier)
#   + s(age, by = BMI) + carrier:BMI (+ cohort x CNV), reference SSC unaffected relatives; ordinal = cumulative logit of the raw
#   total (40 levels), per SRS form, inverse-variance combined. Pooled = four CNVs with cohort x CNV indicators.
# null model: none; Wald tests of the carrier:BMI coefficient (bam, Box-Cox scale) and of the same term in the cumulative-logit model.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
v2_base()
PHS <- c("SRS", "SCQ"); REL <- "Unaffected relatives"; PRO <- "SSC ASD probands"
X <- XP[pheno %in% PHS]; R <- RP[pheno %in% PHS & variant %in% c("reference", "all carriers")]
# ---------------------------------------------------------------- table: re-run vs the 2026-11-14 CSV
old <- fread(file.path(HERE, "cognition_bmi_slopes_boxcox_combined.csv"))[pheno %in% PHS & variant %in% c("reference", "all carriers")]
T <- merge(R, old[, .(pheno, group, variant, prev_n = n, prev_slope = slope, prev_diff = diff, prev_diff_p = diff_p, prev_ord_p = ord_p)], by = c("pheno", "group", "variant"), all.x = TRUE, sort = FALSE)
T[, moved := !(n == prev_n & abs(slope - prev_slope) < 1e-6 & (is.na(diff_p) | abs(diff_p - prev_diff_p) < 1e-6) & (is.na(ord_p) | abs(ord_p - prev_ord_p) < 1e-3))]
setnames(T, c("diff", "diff_se", "diff_p"), c("cnv_x_bmi", "cnv_x_bmi_se", "cnv_x_bmi_p_boxcox"))
T[, source := "51b models re-run 2026-11-20; carriers G2MH + Simons VIP + Simons Searchlight + SSC; reference SSC unaffected relatives"]
fwrite(T, file.path(V2OUT, "block_srs_scq.csv"))
nc <- dcast(X[grp %in% LOC4, .N, by = .(pheno, grp, cohort, sex)], pheno + grp + sex ~ cohort, value.var = "N", fill = 0); fwrite(nc, file.path(V2OUT, "block_srs_scq_n_by_cohort.csv"))
options(width = 220); print(T[, .(pheno, group, n, slope = round(slope, 3), cxb = round(cnv_x_bmi, 3), p_bc = signif(cnv_x_bmi_p_boxcox, 3), p_ord = signif(ord_p, 3), moved)]); print(nc)
# ---------------------------------------------------------------- panels
XLm <- c(-3.2, 4.2); BR <- c(-Inf, -1, 0, 1, 2, Inf)
bins <- function(x) x[, .(bmi = mean(BMI), m = mean(y_adj), se = sd(y_adj) / sqrt(.N), n = .N), by = .(bin = cut(BMI, BR))][n >= 5]
fitline <- function(x) { m <- lm(y_adj ~ BMI, data = x); nd <- data.frame(BMI = seq(max(XLm[1], min(x$BMI)), min(XLm[2], max(x$BMI)), length.out = 50)); cbind(nd, as.data.frame(predict(m, nd, se.fit = TRUE)[1:2])) }
ylim_of <- function(ph) { r <- quantile(X[pheno == ph & grp != REL, y_adj], c(0.005, 0.995)); r + c(-0.06, 0.62) * diff(r) }
pooled_line <- function(ph) { x4 <- X[pheno == ph & grp %in% LOC4]; m <- lm(y_adj ~ BMI + grp, data = x4); g <- seq(max(XLm[1], quantile(x4$BMI, 0.02)), min(XLm[2], quantile(x4$BMI, 0.98)), length.out = 50)
  data.frame(BMI = g, fit = rowMeans(sapply(LOC4, function(L) predict(m, data.frame(BMI = g, grp = L))))) }
YLB <- c(SRS = "SRS total (Box-Cox Z, age-adjusted)", SCQ = "SCQ lifetime (Box-Cox Z, age-adjusted)")
COLS3 <- setNames(c(C_LO, C_MID, C_HI), BAND3)
GC <- c("Unaffected relatives (SSC)" = C_NON, "Autistic, no CNV (SSC probands)" = "#222222", "Four CNVs pooled" = C_CAR)
TS <- 4.1; binrows <- list()
panel_ctl <- function(ph, tag) {
  x <- X[pheno == ph & grp %in% c(REL, PRO)][, g := factor(fifelse(grp == REL, names(GC)[1], names(GC)[2]), names(GC))]; yl <- ylim_of(ph)
  fl <- rbindlist(lapply(levels(x$g)[1:2], function(k) cbind(g = k, fitline(x[g == k]))))[, g := factor(g, names(GC))]
  bb <- rbindlist(lapply(levels(x$g)[1:2], function(k) cbind(g = k, bins(x[g == k]))))[, g := factor(g, names(GC))]
  pl <- as.data.table(pooled_line(ph))[, g := factor(names(GC)[3], names(GC))]
  binrows[[length(binrows) + 1]] <<- cbind(pheno = ph, panel = tag, bb[, .(group = as.character(g), bin, bmi, m, se, n)])
  ru <- R[pheno == ph & group == REL]; rp <- R[pheno == ph & group == PRO]; ra <- R[pheno == ph & group == "All CNVs"]
  l1 <- sprintf("relatives: slope = %+.2f\nautistic, no CNV: slope = %+.2f\n   vs relatives %s", ru$slope, rp$slope, plab(rp$diff_p))
  l2 <- sprintf("four CNVs pooled: slope = %+.2f\nCNV \u00d7 BMI %s (Box-Cox)\n%s (ordinal)", ra$slope, plab(ra$diff_p), plab(ra$ord_p))
  ggplot(x, aes(BMI, y_adj, colour = g)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") + geom_point(size = 0.5, alpha = 0.12, show.legend = FALSE) +
    geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit, fill = g), alpha = 0.25, colour = NA, show.legend = FALSE) +
    geom_line(data = fl, aes(BMI, fit), linewidth = 1.1) + geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), size = 0.5, linewidth = 0.6, shape = 18, show.legend = FALSE) +
    geom_line(data = pl, aes(BMI, fit), linewidth = 1.5) +
    annotate("text", x = XLm[1], y = yl[2], label = l1, hjust = 0, vjust = 1, size = TS, lineheight = 0.92) +
    annotate("label", x = XLm[1], y = yl[2] - 0.2 * diff(yl), label = l2, hjust = 0, vjust = 1, size = TS, lineheight = 0.92, colour = C_CAR, fontface = face(ra$diff_p), fill = "white", alpha = 0.85, linewidth = 0, label.padding = unit(0.12, "lines")) +
    scale_colour_manual(values = GC, drop = FALSE, name = NULL) + scale_fill_manual(values = GC, guide = "none") + coord_cartesian(xlim = XLm, ylim = yl) +
    labs(x = "BMI (Z)", y = YLB[[ph]], tag = tag, title = ph, subtitle = sprintf("relatives n = %d; autistic, no CNV n = %d;\nfour CNVs n = %d", sum(x$grp == REL), sum(x$grp == PRO), ra$n)) + th_v2()
}
panel_cnv <- function(ph, g, tag) {
  x <- X[pheno == ph & grp == g][, band := band3(BMI)]; yl <- ylim_of(ph); fl <- fitline(x); bb <- bins(x); fu <- fitline(X[pheno == ph & grp == REL])
  binrows[[length(binrows) + 1]] <<- cbind(pheno = ph, panel = tag, bb[, .(group = g, bin, bmi, m, se, n)])
  r <- R[pheno == ph & group == g]
  l1 <- sprintf("slope = %+.2f", r$slope); l2 <- sprintf("CNV × BMI %s (Box-Cox)", plab(r$diff_p)); l3 <- sprintf("%s (ordinal)", plab(r$ord_p))
  ggplot(x, aes(BMI, y_adj)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_line(data = fu, aes(BMI, fit), colour = C_NON, linetype = 2, linewidth = 0.9) + geom_point(aes(colour = band), size = 1.4, alpha = 0.7) +
    geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), fill = "grey30", alpha = 0.18) + geom_line(data = fl, aes(BMI, fit), colour = "grey10", linewidth = 1.1) +
    geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), colour = "grey10", size = 0.55, linewidth = 0.6, shape = 18) +
    annotate("text", x = XLm[1], y = yl[2], label = l1, hjust = 0, vjust = 1, size = TS) +
    annotate("text", x = XLm[1], y = yl[2] - 0.07 * diff(yl), label = l2, hjust = 0, vjust = 1, size = TS, fontface = face(r$diff_p)) +
    annotate("text", x = XLm[1], y = yl[2] - 0.14 * diff(yl), label = l3, hjust = 0, vjust = 1, size = TS, fontface = face(r$ord_p)) +
    scale_colour_manual(values = COLS3, drop = FALSE, name = NULL) + coord_cartesian(xlim = XLm, ylim = yl) + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "BMI (Z)", y = NULL, tag = tag, title = g, subtitle = sprintf("n = %d; dashed: relatives", nrow(x))) + th_v2()
}
TAGS <- list(SRS = letters[10:14], SCQ = letters[15:19])
pl <- unlist(lapply(PHS, function(ph) c(list(panel_ctl(ph, TAGS[[ph]][1])), lapply(seq_along(LOC4), function(i) panel_cnv(ph, LOC4[i], TAGS[[ph]][i + 1])))), recursive = FALSE)
fig <- wrap_plots(pl, ncol = 5, byrow = TRUE) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 13))
save_v2(fig, "block_srs_scq", 20, 11)
fwrite(rbindlist(binrows), file.path(V2OUT, "block_srs_scq_binned_means.csv"))
cat("84 done\n")
