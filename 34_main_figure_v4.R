# Main figure v4 (JS / Claude 2026-10-27): v3 layout plus the "prime" panels of the UKBB figure (17_main_with_ukbb.R v2):
# beside each height and head-circumference trajectory panel, the BMI effect on the outcome as a continuous function of
# age, slope(age) = [E(y | BMI = +1 SD, age) - E(y | BMI = -1 SD, age)] / 2, read from the fitted age x BMI surface of the
# carriers-vs-unaffected model (y ~ CNV + s(age) + s(BMI) + s(age, by = CNV) + ti(age, BMI) + ti(age, BMI, by = CNV), REML)
# for carriers and for controls, with 95% CI; the control row uses the within-group model. Clinical cohorts only, same
# data, scaling and tests as v3 (23_unaffected_controls.R, 24_cognition_bmi_slopes.R): the trajectory panel now prints
# the age-adjusted CNV x BMI p and the prime panel the CNV x age x BMI p (both from 23's tests, as reported in the text).
# Columns: Height | Height' | Head circumference | HC' | SRS | SCQ. Output: fig_main_v4.{png,pdf}, fig_main_v4_legend.txt
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))      # A, res, pick(), LOCI, COLS, pfx, YLAB
source(file.path(HERE, "24_cognition_bmi_slopes.R"))     # X, R, panel_ctl(), panel_cnv()
suppressPackageStartupMessages(library(mgcv))
CTL <- "Unaffected controls"; GR1 <- c(CTL, LOCI)
XMAX <- ceiling(quantile(A[pheno %in% c("Height", "Head circumference") & !is.na(age), age], 0.995))   # JS 2026-10-27: axis ends at the data (Jon: the 65-y axis wasted space)

# ---------------------------------------------------------------- trajectory panel: as 23's panel(), with the CNV x BMI line only
traj4 <- function(ph, Gn) {
  x <- A[pheno == ph & locus == Gn & !is.na(tail)][, grp := factor(tail, names(COLS))]; ctl <- Gn == CTL
  r <- pick(ph, Gn, "all carriers vs unaffected"); pv <- if (ctl) NA_real_ else r$cnv_x_bmi_p
  lab <- if (ctl) "" else gsub("= < ", "< ", sprintf("CNV × BMI p = %s", pfx(pv)), fixed = TRUE)
  yl <- range(quantile(A[pheno == ph, y], c(0.005, 0.995))); yl <- yl + c(-0.6, 0.9) * diff(yl) * 0.2
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = if (ctl) 0.5 else 1.2, alpha = if (ctl) 0.15 else 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 6), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = XMAX, y = yl[2], label = lab, hjust = 1, vjust = 1, size = 3, fontface = if (isTRUE(pv < 0.05)) "bold" else "plain") +
    scale_colour_manual(values = COLS, breaks = names(COLS)) + scale_fill_manual(values = COLS, breaks = names(COLS)) + coord_cartesian(xlim = c(2, XMAX), ylim = yl) +
    labs(x = "Age (years)", y = YLAB[[ph]], colour = NULL, fill = NULL, subtitle = sprintf("%s (n = %d)", if (ctl) "Unaffected (SSC sibs/parents, SVIP)" else Gn, nrow(x))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 7), legend.position = "inside", legend.position.inside = c(0.22, 0.1),
                                          legend.background = element_blank(), legend.key.height = unit(7, "pt"))
}

# ---------------------------------------------------------------- prime panel: BMI effect by age (from 17_main_with_ukbb.R v2)
GRID <- seq(2.5, 60, by = 0.5)
grid_of <- function(a) GRID[GRID >= quantile(a, 0.025) & GRID <= quantile(a, 0.975)]   # JS 2026-10-27: draw each curve only over the central 95% of that group's ages (the tails beyond the data support bent downward: smoother edge effect, Jon)
support_of <- function(a_grid, ages, w = 3, nmin = 5) sapply(a_grid, function(a) sum(abs(ages - a) <= w) >= nmin)   # JS 2026-10-27: carrier curve drawn only where >= 5 carriers lie within 3 y (removes the dip across the 18-30 y gap of the clinical carriers)
K_AGE <- 10; K_BMI <- 6; K_TI <- c(6, 5); GAMMA <- 1.5   # JS 2026-10-27: same smoother settings as 17_main_with_ukbb.R
slope_rows <- function(m, grid, groups) {
  V <- vcov(m); B <- coef(m)
  rbindlist(lapply(groups, function(g) {
    nd <- function(b) { d <- data.frame(age = grid, BMI = b); if (!is.null(g)) d$is_cnv <- factor(g, levels = c("No CNV", "CNV")); d }
    Xm <- (predict(m, nd(1), type = "lpmatrix") - predict(m, nd(-1), type = "lpmatrix")) / 2
    data.table(who = if (is.null(g) || g == "No CNV") "Non-carriers" else "Carriers", age = grid, slope = drop(Xm %*% B), se = sqrt(rowSums((Xm %*% V) * Xm)))
  }))
}
bmi_effect_by_age <- function(ph, Gn) {
  ctl <- A[pheno == ph & locus == CTL & !is.na(y) & !is.na(BMI) & !is.na(age), .(y, age, BMI, locus = "No CNV")]
  if (Gn == CTL) {
    m <- bam(y ~ s(age, k = K_AGE) + s(BMI, k = K_BMI) + ti(age, BMI, k = K_TI), data = ctl, method = "REML", gamma = GAMMA)
    return(slope_rows(m, grid_of(ctl$age), list(NULL)))
  }
  xc <- A[pheno == ph & locus == Gn & !is.na(y) & !is.na(BMI) & !is.na(age), .(y, age, BMI, locus = Gn)]
  x <- rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- bam(y ~ is_cnv + s(age, k = K_AGE) + s(BMI, k = K_BMI) + s(age, k = K_AGE, by = is_cnv) + ti(age, BMI, k = K_TI) + ti(age, BMI, k = K_TI, by = is_cnv), data = x, method = "REML", gamma = GAMMA)
  r <- slope_rows(m, grid_of(xc$age), list("No CNV", "CNV")); r[who == "Carriers" & !support_of(age, xc$age), `:=`(slope = NA_real_, se = NA_real_)]; r
}
BEACOLS <- c("Non-carriers" = "#8C9A9A", "Carriers" = "#4B6EA9")
eff_panel <- function(ph, Gn, ylab, ylim = c(-1, 1.2)) {
  r <- bmi_effect_by_age(ph, Gn)[, who := factor(who, names(BEACOLS))]; ctl <- Gn == CTL
  pv <- if (ctl) NA_real_ else pick(ph, Gn, "all carriers vs unaffected")$cnv_x_age_x_bmi_p   # the three-way test of 23, as reported in the text
  lab <- if (ctl) "" else gsub("= < ", "< ", sprintf("CNV × age × BMI p = %s", pfx(pv)), fixed = TRUE)
  carriers <- A[pheno == ph & locus == Gn & !is.na(tail)]
  g <- ggplot(r, aes(age, slope, colour = who, fill = who)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_ribbon(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.9)
  if (!ctl) g <- g + geom_rug(data = carriers, aes(x = age), inherit.aes = FALSE, sides = "b", alpha = 0.35, length = unit(0.05, "npc"), colour = BEACOLS[["Carriers"]])
  g + annotate("text", x = XMAX, y = ylim[2], label = lab, hjust = 1, vjust = 1, size = 3, fontface = if (isTRUE(pv < 0.05)) "bold" else "plain") +
    scale_colour_manual(values = BEACOLS, breaks = if (ctl) "Non-carriers" else names(BEACOLS)) + scale_fill_manual(values = BEACOLS, breaks = if (ctl) "Non-carriers" else names(BEACOLS)) +
    coord_cartesian(xlim = c(2, XMAX), ylim = ylim) +
    labs(x = "Age (years)", y = ylab, colour = NULL, fill = NULL, subtitle = "BMI effect by age (95% CI)") +
    theme_classic(base_size = 13) + theme(legend.position = "inside", legend.position.inside = c(0.25, 0.14), legend.text = element_text(size = 7),
                                          legend.background = element_blank(), legend.key.height = unit(7, "pt"), plot.subtitle = element_text(size = 10))
}

# ---------------------------------------------------------------- assemble: Height | Height' | HC | HC' | SRS | SCQ, rows = controls + 4 CNVs
cog <- function(ph) c(list(panel_ctl(ph, probands = TRUE, show_title = FALSE)), lapply(LOCI, function(g) panel_cnv(ph, g)))
cols <- list(lapply(GR1, function(Gn) traj4("Height", Gn)), lapply(GR1, function(Gn) eff_panel("Height", Gn, "Height per +1 SD BMI (SD)")),
             lapply(GR1, function(Gn) traj4("Head circumference", Gn)), lapply(GR1, function(Gn) eff_panel("Head circumference", Gn, "Head circ. per +1 SD BMI (SD)")),
             cog("SRS"), cog("SCQ"))
tagsets <- list(LETTERS[1:5], paste0(LETTERS[1:5], "'"), LETTERS[6:10], paste0(LETTERS[6:10], "'"), LETTERS[11:15], LETTERS[16:20])
plots <- unlist(lapply(seq_along(cols), function(j) Map(function(p, tg) p + labs(tag = tg), cols[[j]], tagsets[[j]])), recursive = FALSE)
hdr <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 9, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
W <- c(2, 1.3, 2, 1.3, 1, 1)
hdr_row <- (hdr("Height") | hdr("") | hdr("Head circumference") | hdr("") | hdr("SRS") | hdr("SCQ")) + plot_layout(widths = W)
body <- wrap_plots(plots, ncol = 6, byrow = FALSE, widths = W) & theme(plot.tag = element_text(size = 18, face = "bold"))
fig <- hdr_row / body + plot_layout(heights = c(0.04, 1))
ggsave(file.path(HERE, "fig_main_v4.png"), fig, width = 22, height = 20.5, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_main_v4.pdf"), fig, width = 22, height = 20.5, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_main_v4.pdf"), fig, width = 22, height = 20.5))

# ---------------------------------------------------------------- legend: v3 text plus the prime panels
leg3 <- readLines(file.path(HERE, "fig_main_v3_legend.txt"), warn = FALSE)
leg4 <- c(sub("^Figure\\.", "Figure (v4).", leg3[1]),
  "Columns: height (A-E) with the BMI effect on height by age (A'-E'), head circumference (F-J) with the BMI effect by age (F'-J'), SRS (K-O) and SCQ (P-T). Rows as in v3: unaffected non-carrier controls (top), then carriers of 16p11.2 deletion, 16p11.2 duplication, 22q11.2 deletion and 22q11.2 duplication.",
  leg3[-(1:2)],
  "A'-E', F'-J': the effect of BMI on the outcome as a continuous function of age, half the difference in the predicted outcome between BMI = +1 SD and BMI = -1 SD at each age (SD of the outcome per SD of BMI), with 95% confidence bands, read from the fitted age x BMI surface of the model that compares carriers with unaffected controls (y ~ CNV + s(age) + s(BMI) + s(age, by CNV) + ti(age, BMI) + ti(age, BMI, by CNV), REML): carriers in blue, controls in grey; so each panel shows where in the age range the carriers' BMI effect departs from the controls'. The control row shows the within-group surface (y ~ s(age) + s(BMI) + ti(age, BMI)). Tick marks along the bottom mark the carriers' ages. Curves are drawn over the central 95% of each group's ages, and the carrier curve only at ages with at least five carriers within three years (gaps elsewhere); the tensor-product surfaces use basis sizes 10 (age), 6 (BMI) and 6 x 5 (age x BMI) with the smoothing penalty inflated by gamma = 1.5, after the k = 15/10/10 x 10 surfaces proved too wiggly and bent at the ends of the age range. The trajectory panels print the age-adjusted CNV x BMI p and the prime panels the CNV x age x BMI p, both from the tests reported in the text (23_unaffected_controls.R); the CNV x age x BMI term is sensitive to extreme observations and is reported as exploratory.")
writeLines(leg4, file.path(HERE, "fig_main_v4_legend.txt"))
cat("fig_main_v4 written\n")
