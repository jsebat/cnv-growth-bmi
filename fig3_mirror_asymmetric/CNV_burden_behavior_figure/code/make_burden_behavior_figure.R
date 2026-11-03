#!/usr/bin/env Rscript
required <- c("ggplot2", "dplyr", "readr", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
library(ggplot2); library(dplyr); library(readr); library(patchwork)
set.seed(20260728)
root <- normalizePath(".")
dir.create(file.path(root, "figures_R"), showWarnings = FALSE)
height_b <- read_csv(file.path(root, "inputs", "height_burden_predictions.csv"), show_col_types = FALSE)
bmi_b <- read_csv(file.path(root, "inputs", "bmi_burden_predictions.csv"), show_col_types = FALSE)
lit <- read_csv(file.path(root, "inputs", "literature_burden_targets.csv"), show_col_types = FALSE)
models <- c("Mirror alone","Mirror + continuous","Mirror + single hit")
make_clouds <- function(df, trait) {
  bind_rows(lapply(seq_len(nrow(df)), function(i) {
    r <- df[i,]
    dens0 <- r$effective_dosage_genes_per_Mb_to_match_DEL
    dens <- rlnorm(800, meanlog = log(dens0) - 0.5 * log(1 + 0.35^2), sdlog = sqrt(log(1 + 0.35^2)))
    dens <- pmax(0.2, pmin(dens, dens0 * 2.2))
    del_coef <- dens * r$DEL_effect_per_affected_gene_native * rnorm(length(dens), 1, 0.06)
    dup_coef <- dens * r$DUP_effect_per_affected_gene_native * rnorm(length(dens), 1, 0.10)
    tibble(trait = trait, model = r$model, density = dens, DEL_per_Mb = del_coef, DUP_per_Mb = dup_coef,
           matched_DEL = dens0 * r$DEL_effect_per_affected_gene_native,
           matched_DUP = dens0 * r$DUP_effect_per_affected_gene_native)
  }))
}
clouds <- bind_rows(make_clouds(height_b, "height"), make_clouds(bmi_b, "BMI"))
write_csv(clouds, file.path(root, "outputs", "simulated_burden_clouds_R.csv"))
lims <- list(height = list(x = c(-0.52, 0.06), y = c(-0.22, 0.18)), BMI = list(x = c(-0.02, 0.22), y = c(-0.10, 0.08)))
model_colors <- c("Mirror alone" = "grey40", "Mirror + continuous" = "seagreen4", "Mirror + single hit" = "slateblue3")
base_theme <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))
make_obs <- function(trait) {
  x <- lit %>% filter(trait == !!trait, copy_state == "deletion") %>% pull(estimate) %>% .[1]
  xlim <- lims[[trait]]$x; ylim <- lims[[trait]]$y
  p <- ggplot() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + coord_cartesian(xlim = xlim, ylim = ylim) + labs(x = "DEL burden effect (native units per Mb)", y = "DUP burden effect (native units per Mb)") + base_theme
  if (trait == "height") {
    p <- p + annotate("rect", xmin = xlim[1-1], xmax = 0, ymin = ylim[1-1], ymax = 0, alpha = 0.18, fill = "lightsteelblue") + annotate("segment", x = x, xend = x, y = 0.005, yend = ylim[1] * 0.75, arrow = arrow(length = unit(0.15, "cm"))) + annotate("point", x = x, y = 0, size = 2.5) + annotate("text", x = x, y = 0.012, label = "Macé DEL
−0.41 cm/Mb", size = 2.5, vjust = 0) + annotate("text", x = xlim[1] + 0.02, y = ylim[2] - 0.01, label = "Observed burden pattern
DEL < 0 (Macé)
DUP < 0 (Auwerx direction)", hjust = 0, vjust = 1, size = 3)
  } else {
    p <- p + annotate("rect", xmin = 0, xmax = xlim[2], ymin = 0, ymax = ylim[2], alpha = 0.22, fill = "mistyrose") + annotate("segment", x = x, xend = x, y = -0.003, yend = ylim[2] * 0.75, arrow = arrow(length = unit(0.15, "cm"))) + annotate("point", x = x, y = 0, size = 2.5) + annotate("text", x = x, y = -0.009, label = "Macé DEL
+0.14 kg/m²/Mb", size = 2.5, vjust = 1) + annotate("text", x = xlim[1] + 0.02, y = ylim[2] - 0.01, label = "Observed burden pattern
DEL > 0 (Macé)
DUP > 0 (Auwerx direction)", hjust = 0, vjust = 1, size = 3)
  }
  p
}
make_mod <- function(trait, model) {
  d <- clouds %>% filter(trait == !!trait, model == !!model)
  qgrid <- quantile(d$density, probs = seq(0.05, 0.95, length.out = 18))
  mids <- bind_rows(lapply(qgrid, function(q) { s <- d %>% slice_min(abs(density - q), n = 60, with_ties = FALSE); tibble(DEL_per_Mb = mean(s$DEL_per_Mb), DUP_per_Mb = mean(s$DUP_per_Mb)) }))
  ggplot(d, aes(DEL_per_Mb, DUP_per_Mb)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_point(alpha = 0.16, size = 0.6, color = model_colors[model]) + geom_path(data = mids, linewidth = 1.1, color = model_colors[model]) + geom_point(aes(x = matched_DEL, y = matched_DUP), data = d %>% slice(1), inherit.aes = FALSE, shape = 4, stroke = 1.5, size = 3.2, color = "darkorange") + coord_cartesian(xlim = lims[[trait]]$x, ylim = lims[[trait]]$y) + labs(x = "DEL burden effect (native units per Mb)", y = ifelse(model == models[1], "DUP burden effect (native units per Mb)", "")) + base_theme
}
fig <- ((make_obs("height") + plot_annotation(title = "A  Observed height burden evidence")) | (make_mod("height", models[1]) + plot_annotation(title = "B  Mirror
Predicts DUP burden in the wrong direction")) | (make_mod("height", models[2]) + plot_annotation(title = "C  Continuous
Shifts DUP burden below zero")) | (make_mod("height", models[3]) + plot_annotation(title = "D  Single hit
Also shifts DUP burden below zero"))) / ((make_obs("BMI") + plot_annotation(title = "E  Observed BMI burden evidence")) | (make_mod("BMI", models[1]) + plot_annotation(title = "F  Mirror
Predicts DUP burden in the wrong direction")) | (make_mod("BMI", models[2]) + plot_annotation(title = "G  Continuous
Shifts DUP burden above zero")) | (make_mod("BMI", models[3]) + plot_annotation(title = "H  Single hit
Also shifts DUP burden above zero"))) + plot_annotation(title = "Genome-wide CNV burden predictions under mirror, continuous, and sparse single-hit models")
ggsave(file.path(root, "figures_R", "Figure_burden_model_behavior.pdf"), fig, width = 18, height = 9, units = "in")
ggsave(file.path(root, "figures_R", "Figure_burden_model_behavior.png"), fig, width = 18, height = 9, units = "in", dpi = 300)
