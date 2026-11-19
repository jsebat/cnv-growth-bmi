# Grant figure, compact (JS / Claude 2026-11-15, Jon's sketch): left column = head circumference by BMI tail for controls and the four CNVs
#   (35e panels a-e; pooled high - low BMI offset vs controls in a); right column = pooled results only:
#   f  BMI effect on head circumference by age band, four CNVs pooled vs non-carrier controls (band_slopes() of 17, pooled carriers;
#      printed p = CNV x BMI within band)
#   g  SRS vs BMI: controls, ASD without CNV, four CNVs pooled (as panel_ctl of 51b = Fig. 5k, rebuilt here with only the pooled slope and its
#      Box-Cox p vs controls printed; Jon 2026-11-15: "the caption should show just the first 2 lines of red text")
#   h  BMI slope on SRS (SD per SD BMI, 95% CI): controls, ASD no CNV, four CNVs pooled; p vs controls (cognition_bmi_slopes_boxcox_combined.csv, R)
#   Jon 2026-11-15: keys, captions, axis labels and titles enlarged for a proposal page (BIG below scales every in-panel text layer).
# RUNS ON EXPANSE in main_with_ukbb/v8_grant. Output: fig_grant_pooled.{png,pdf}, grant_hc_pooled_band_slopes.csv
source("17_main_with_ukbb.R")
load(file.path(WD, "cog_panels_boxcox_combined.RData"))
LOCI_C <- LOCI; LOC <- LOCI; BR <- c(-Inf, -1, 0, 1, 2, Inf)
BIG <- 1.8                                                         # scale for text layers inside panels (p values, n labels)
bump <- function(p, f = BIG) { for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomText") && !is.null(p$layers[[i]]$aes_params$size))
  p$layers[[i]]$aes_params$size <- p$layers[[i]]$aes_params$size * f; p }
TXT <- theme(axis.title = element_text(size = 21), axis.text = element_text(size = 17), legend.text = element_text(size = 16),
             legend.key.size = unit(20, "pt"), plot.subtitle = element_text(face = "bold", size = 21))
# ---- a-e: HC trajectories
po <- offset_vs_ctl(dh[locus %in% LOCI], dh[locus == "No CNV"])
# Jon 2026-11-15: no p values in panels a-e (the CNV x BMI text layers of 17 are removed; the pooled offset test is printed in the log only)
notext <- function(p) { p$layers <- Filter(function(l) !inherits(l$geom, "GeomText"), p$layers); p }
tr <- lapply(GROUPS, function(G) notext(v2[[paste("HC", G)]] + labs(tag = NULL)) + TXT + theme(legend.position.inside = c(0.72, 0.9)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.6)), fill = guide_legend(override.aes = list(alpha = 0.3))))
cat(sprintf("pooled HC high-low BMI offset vs controls: %+.3f SD, p = %.2e (not printed in the figure)\n", po[["diff"]], po[["p"]]))
tr[[1]] <- tr[[1]] + theme(legend.position.inside = c(0.75, 0.12))
# ---- f: pooled HC BMI effect by age band
bp <- band_slopes(dh[locus %in% LOCI], dh[locus == "No CNV"]); fwrite(bp, file.path(WD, "grant_hc_pooled_band_slopes.csv")); print(bp)
lv <- BLAB[BLAB %in% bp$age_bin]
pf <- band_panel(bp, "HC per +1 SD BMI (SD)", FALSE, yl_band(bp), lv)
pf$layers[[length(pf$layers)]]$aes_params$size <- 5.4           # p labels
pf <- pf + scale_x_discrete(drop = FALSE, labels = pf$scales$get_scales("x")$labels, expand = expansion(add = 0.75))
pf <- pf + labs(subtitle = "Head circumference:\nBMI effect by age (95% CI)") +
  scale_fill_manual(values = BEACOLS, labels = c("Non-carriers" = "Controls", "Carriers" = "CNV carriers")) + TXT +
  theme(axis.text.x = element_text(size = 15), legend.key.size = unit(18, "pt"))
# ---- g: SRS scatter, pooled (rebuilt; only the pooled slope + Box-Cox p printed)
ph_ <- "SRS"; G <- names(REFC); x <- X[pheno == ph_ & grp %in% G]; yl <- ylim_of(ph_, FALSE)
fl <- rbindlist(lapply(G, function(g) cbind(grp = g, fitline(x[grp == g])))); bb <- rbindlist(lapply(G, function(g) cbind(grp = g, bins(x[grp == g]))))
ra <- R[pheno == ph_ & group == "All CNVs" & variant == "all carriers"]; pl <- pooled_line(ph_)
labp <- gsub("= < ", "< ", sprintf("4 CNVs pooled (red):\nslope = %+.2f\nvs controls: p = %s", ra$slope, ptxt(ra$diff_p)), fixed = TRUE)
pg <- ggplot(x, aes(BMI, y_adj, colour = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 0.5, alpha = 0.12) + geom_ribbon(data = fl, aes(BMI, fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit, fill = grp), alpha = 0.25, colour = NA) +
  geom_line(data = fl, aes(BMI, fit), linewidth = 1.2) + geom_pointrange(data = bb, aes(bmi, m, ymin = m - se, ymax = m + se), size = 0.5, linewidth = 0.7, shape = 18) +
  geom_line(data = pl, aes(BMI, fit), inherit.aes = FALSE, colour = POOLC, linewidth = 1.8) +
  annotate("text", x = XL[2], y = yl[2], label = labp, hjust = 1, vjust = 1, size = 5.8, lineheight = 0.95, colour = POOLC, fontface = if (isTRUE(ra$diff_p < 0.05)) "bold" else "plain") +
  scale_colour_manual(values = REFC, labels = c("Unaffected relatives" = "controls", "SSC ASD probands" = "ASD, no CNV")) + scale_fill_manual(values = REFC, guide = "none") +
  coord_cartesian(xlim = XL, ylim = yl) +
  labs(x = "BMI (Z)", y = "SRS (Z, age-adjusted)", colour = NULL, subtitle = "SRS vs BMI") +
  theme_classic(base_size = 15) + TXT + theme(plot.subtitle = element_text(face = "bold", size = 17), legend.position = "inside", legend.position.inside = c(0.3, 0.72),
  legend.background = element_blank()) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.8, linewidth = 1.6)))
# ---- h: SRS BMI slopes, pooled
rs <- R[pheno == "SRS" & ((group %in% c("Unaffected relatives", "SSC ASD probands") & variant == "reference") | (group == "All CNVs" & variant == "all carriers"))]
rs[, lab := factor(c("Unaffected relatives" = "Controls", "SSC ASD probands" = "ASD,\nno CNV", "All CNVs" = "4 CNVs\npooled")[group], c("Controls", "ASD,\nno CNV", "4 CNVs\npooled"))]
rs[, col := c("Unaffected relatives" = "#8C9A9A", "SSC ASD probands" = "#1a1a1a", "All CNVs" = POOLC)[group]]
rs[, ptx := ifelse(is.na(diff_p), "", gsub("= < ", "< ", sprintf("p = %s", pfx(diff_p)), fixed = TRUE))]
print(rs[, .(group, n, slope, se, diff_p)])
ph <- ggplot(rs, aes(lab, slope, fill = col)) + geom_hline(yintercept = 0, colour = "grey40") + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), width = 0.2, colour = "grey20") +
  geom_text(aes(y = slope + 1.96 * se + 0.012, label = ptx, fontface = ifelse(!is.na(diff_p) & diff_p < 0.05, "bold", "plain")), size = 6.5, vjust = 0) +
  scale_fill_identity() + coord_cartesian(ylim = c(-0.02, max(rs$slope + 1.96 * rs$se) + 0.05)) +
  labs(x = NULL, y = "SRS per +1 SD BMI", subtitle = "SRS: BMI effect (95% CI)") +
  theme_classic(base_size = 15) + TXT + theme(axis.text.x = element_text(size = 16))
# ---- layout
P <- c(Map(function(p, tg) p + labs(tag = tg), tr, letters[1:5]), list(pf + labs(tag = "f"), pg + labs(tag = "g"), ph + labs(tag = "h")))
des <- "AF\nBF\nCG\nDG\nEH"
hdr5 <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 11, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
hdr_row <- (hdr5("Head circumference") | hdr5("4 CNVs pooled")) + plot_layout(widths = c(2, 1.6))
body <- wrap_plots(P, design = des, widths = c(2, 1.6)) & theme(plot.tag = element_text(size = 24, face = "bold"))
fig <- hdr_row / body + plot_layout(heights = c(0.04, 1))
ggsave(file.path(WD, "fig_grant_pooled.png"), fig, width = 14, height = 22, dpi = 150)
tryCatch(ggsave(file.path(WD, "fig_grant_pooled.pdf"), fig, width = 14, height = 22, device = cairo_pdf), error = function(e) ggsave(file.path(WD, "fig_grant_pooled.pdf"), fig, width = 14, height = 22))
cat("fig_grant_pooled written\n")
