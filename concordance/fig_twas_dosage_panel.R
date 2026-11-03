# Single panel (Jon 2026-11-02): gene-dosage effect of each reciprocal CNV on height, (DUP - DEL)/2, against the summed TWAS z of
# its genes (higher z = genes whose higher predicted expression goes with greater height). Mirror CNVs blue, asymmetric red (as in
# Fig. 3). One Theil-Sen line through all CNVs except Williams (7q11.23, open symbol): without Williams, mirror and asymmetric CNVs
# share one line (class adds nothing; twas_same_line_test.py). Spearman rho and circular-permutation P from
# twas_dosage_panel_stats.py. Output: fig_twas_dosage_panel.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/concordance"
# Jon 2026-11-02 (revised): non-overlapping CNVs only, keeping the spanning intervals (22q11.2 A-D, 16p11.2 BP1-BP3, 1q21.1 TAR
# BP1-BP3, 16p13.11 BP1-BP3) and dropping the subregions inside them, so every point is independent (twas_nonoverlapping.py)
d <- fread(file.path(HERE, "twas_dosage_panel_nonoverlap.csv")); ns <- fread(file.path(HERE, "twas_dosage_panel_nonoverlap_stats.csv"))
d[, cls := factor(cls, c("mirror", "asymmetric"), c("Mirror", "Asymmetric"))][, williams := CNV_ID == "7q11.23_WBS"]
COL <- c(Mirror = "#2C6FB7", Asymmetric = "#C8312F")
fp <- function(p) formatC(p, format = "fg", digits = 2)
txt1 <- sprintf("All CNVs except Williams: \u03c1 = %.2f, P = %s", ns$rho_all_noW, fp(ns$p_all_noW))
txt2 <- sprintf("Mirror \u03c1 = %.2f; asymmetric \u03c1 = %.2f (without Williams)", ns$rho_mirror, ns$rho_asym_noW)
txt3 <- sprintf("Separate lines for the two classes: P = %s", fp(ns$p_separate_lines))
xr <- range(d[williams == FALSE, summed_twas_z])
LN <- data.table(x0 = xr[1], x1 = xr[2], y0 = ns$theilsen_intercept + ns$theilsen_slope * xr[1], y1 = ns$theilsen_intercept + ns$theilsen_slope * xr[2])
LBL <- c("16p11.2_BP4_BP5" = "16p11.2 BP4\u2013BP5", "3q29" = "3q29", "1q21.1" = "1q21.1", "22q11.21_A-D" = "22q11.2 A\u2013D",
         "7q11.23_WBS" = "7q11.23 (Williams)", "16p11.2_BP1-BP3" = "16p11.2 BP1\u2013BP3")
show <- d[CNV_ID %in% names(LBL)][, lab := LBL[CNV_ID]]
POS <- data.table(CNV_ID = c("1q21.1", "16p11.2_BP4_BP5", "3q29", "22q11.21_A-D", "7q11.23_WBS", "16p11.2_BP1-BP3", "2q11.2"),
                  dx = c(-0.5, 0, 0, -0.35, 0.45, 0, 0), dy = c(0, 0.045, -0.05, 0.07, 0, 0.05, -0.05), hj = c(1, 0.5, 0.5, 1, 0, 0.5, 0.5))
show <- merge(show, POS, by = "CNV_ID")
p <- ggplot(d, aes(summed_twas_z, dosage_effect, colour = cls)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_segment(data = LN, aes(x = x0, xend = x1, y = y0, yend = y1), colour = "grey25", linewidth = 0.7, inherit.aes = FALSE) +
  geom_point(data = d[williams == FALSE], size = 2.2, alpha = 0.9) +
  geom_point(data = d[williams == TRUE], size = 2.2, shape = 21, stroke = 1, fill = "white") +
  geom_text(data = show, aes(x = summed_twas_z + dx, y = dosage_effect + dy, label = lab, hjust = hj), size = 2.5, vjust = 0.5, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = txt1, hjust = -0.04, vjust = 1.4, size = 2.9, colour = "grey15") +
  annotate("text", x = -Inf, y = Inf, label = txt2, hjust = -0.04, vjust = 3.3, size = 2.5, colour = "grey30") +
  annotate("text", x = -Inf, y = Inf, label = txt3, hjust = -0.04, vjust = 5.0, size = 2.5, colour = "grey30") +
  annotate("text", x = Inf, y = -Inf, label = "Mirror", hjust = 1.1, vjust = -3.0, size = 2.9, colour = COL[["Mirror"]]) +
  annotate("text", x = Inf, y = -Inf, label = "Asymmetric", hjust = 1.1, vjust = -1.4, size = 2.9, colour = COL[["Asymmetric"]]) +
  coord_cartesian(ylim = c(-0.45, 0.9)) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(x = "Summed TWAS z of genes in the CNV\n(higher predicted expression \u2192 taller)",
       y = "CNV dosage effect on height (SD)\n(duplication \u2212 deletion) / 2") +
  theme_classic(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(HERE, "fig_twas_dosage_panel.png"), p, width = 4.2, height = 3.9, dpi = 300, bg = "white")
ggsave(file.path(HERE, "fig_twas_dosage_panel.pdf"), p, width = 4.2, height = 3.9, device = cairo_pdf)
cat("panel written\n")
