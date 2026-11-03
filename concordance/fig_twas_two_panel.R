# Two-panel TWAS concordance figure (Jon 2026-11-02).
#  a  CNV dosage effect on height vs summed TWAS z (fig_twas_dosage_panel.R, sourced)
#  b  mean summed TWAS z of mirror and asymmetric CNVs vs their size-matched null (circular permutation, 10,000 draws;
#     twas_class_null_panel.py); grey violin = null distribution of the class mean, bar = 95% interval, diamond = observed.
# Same 23 non-overlapping CNVs in both panels. Output: fig_twas_concordance_2panel.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/concordance"
source(file.path(HERE, "fig_twas_dosage_panel.R"))   # defines p (panel a) and COL, fp
pa <- p
nd <- melt(fread(file.path(HERE, "twas_class_null_draws.csv")), measure.vars = c("mirror", "asymmetric"), variable.name = "cls", value.name = "null_mean")
st <- fread(file.path(HERE, "twas_class_null_stats.csv"))
lv <- c(mirror = "Mirror", asymmetric = "Asymmetric")
nd[, cls := factor(lv[as.character(cls)], c("Mirror", "Asymmetric"))]
so <- st[group %in% names(lv)][, cls := factor(lv[group], c("Mirror", "Asymmetric"))]
so[, ptxt := sprintf("P = %s", fp(p_two_sided))][, ntxt := sprintf("%d CNVs\nmedian %d genes", n_cnv, as.integer(median_genes))]
dif <- st[group == "asymmetric minus mirror"]
pb <- ggplot(nd, aes(cls, null_mean)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_violin(fill = "grey88", colour = "grey60", linewidth = 0.3, width = 0.8, scale = "width", trim = TRUE) +
  geom_errorbar(data = so, aes(x = cls, ymin = null_q025, ymax = null_q975), inherit.aes = FALSE, width = 0.12, colour = "grey40", linewidth = 0.4) +
  geom_point(data = so, aes(x = cls, y = observed_mean_summed_z, fill = cls), inherit.aes = FALSE, shape = 23, size = 3.4, colour = "grey15") +
  geom_text(data = so, aes(x = cls, y = 8.2, label = ptxt), inherit.aes = FALSE, size = 2.8) +
  annotate("text", x = 1.5, y = 9.7, label = sprintf("Asymmetric vs mirror: P = %s", fp(dif$p_two_sided)), size = 2.6, colour = "grey30") +
  scale_x_discrete(labels = XL) +
  scale_fill_manual(values = COL, guide = "none") +
  coord_cartesian(ylim = c(-8, 10)) +
  labs(x = NULL, y = "Mean summed TWAS z per CNV\n(observed vs size-matched null)") +
  theme_classic(base_size = 10)
fig <- (pa | pb) + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 12))
ggsave(file.path(HERE, "fig_twas_concordance_2panel.png"), fig, width = 7.2, height = 3.9, dpi = 300, bg = "white")
ggsave(file.path(HERE, "fig_twas_concordance_2panel.pdf"), fig, width = 7.2, height = 3.9, device = cairo_pdf)
cat("two-panel figure written\n")
