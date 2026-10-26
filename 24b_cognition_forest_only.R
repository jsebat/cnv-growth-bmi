# Forest-plot-only version of the SRS/SCQ supplementary figure (JS / Claude 2026-10-26): the scatter rows of fig_cognition_bmi_slopes
# duplicate main-figure panels K-T, so the manuscript's Supplementary Fig. 2 keeps only the two forest plots. Reads cognition_bmi_slopes.csv
# (written by 24_cognition_bmi_slopes.R), so run after 24. Output: fig_cognition_bmi_forest.{png,pdf}.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
R <- fread(file.path(HERE, "cognition_bmi_slopes.csv"))
LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
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
    labs(x = sprintf("BMI slope on %s (SD per SD of BMI, 95%% CI)", ph), y = NULL, colour = NULL, shape = NULL, title = ph, subtitle = "BMI slope by group (* vs unaffected p < 0.05)") +
    theme_classic(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "inside", legend.position.inside = c(0.84, 0.86), legend.background = element_blank(), legend.text = element_text(size = 8), legend.key.height = unit(9, "pt"))
}
fig <- (panel_forest("SRS") | panel_forest("SCQ")) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(HERE, "fig_cognition_bmi_forest.png"), fig, width = 11.5, height = 4.6, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_cognition_bmi_forest.pdf"), fig, width = 11.5, height = 4.6, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_cognition_bmi_forest.pdf"), fig, width = 11.5, height = 4.6))
print(R[variant != "reference" | group %in% c("SSC siblings", "SSC parents", "SSC ASD probands"), .(pheno, group, variant, n, slope = round(slope, 2), ci_lo = round(ci_lo, 2), ci_hi = round(ci_hi, 2), diff_p = signif(diff_p, 2))], nrows = 80)
cat("done\n")
