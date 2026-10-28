#!/usr/bin/env Rscript
required <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)
set.seed(20260728)
root <- normalizePath(".")
dir.create(file.path(root, "figures_R"), showWarnings = FALSE)
dir.create(file.path(root, "outputs_R"), showWarnings = FALSE)
height <- read_csv(file.path(root, "inputs", "height_loci.csv"), show_col_types = FALSE)
bmi <- read_csv(file.path(root, "inputs", "bmi_loci.csv"), show_col_types = FALSE)
pars <- read_csv(file.path(root, "inputs", "model_parameters.csv"), show_col_types = FALSE)
MIRROR <- "Mirror alone"; CONT <- "Mirror + continuous directional shift"; SINGLE <- "Mirror + sparse single hit"
get_par <- function(trait, model) pars %>% filter(Trait == trait, Model == model) %>% slice(1)
simulate_trait <- function(trait, model, n = 1000, load = 1, seed = 1, include_sampling = TRUE) {
  set.seed(seed)
  dat <- if (trait == "Height") height else bmi
  idx <- sample(seq_len(nrow(dat)), n, replace = TRUE)
  G <- dat$Gene_count[idx]; seD <- dat$Deletion_SE[idx]; seU <- dat$Duplication_SE[idx]
  z <- get_par(trait, model)
  mu_gene <- z$p * z$sign_mean * z$m
  var_gene <- z$p * (2 * z$m^2) - (z$p * z$sign_mean * z$m)^2
  mean_D <- G * mu_gene
  var_D <- pmax(1e-9, G * var_gene + z$tauD^2)
  D <- rnorm(n, mean_D, sqrt(var_D))
  U <- -z$rho * D + rnorm(n, 0, z$tauU)
  if (model == CONT) {
    b <- z$shift_magnitude * load; theta <- z$shift_angle_deg * pi / 180
    shift_mean <- G * z$q * b
    shift_var <- pmax(1e-9, G * (z$q * 2 * b^2 - (z$q * b)^2))
    total <- pmax(0, rnorm(n, shift_mean, sqrt(shift_var)))
    D <- D + z$trait_sign * total * cos(theta)
    U <- U + z$trait_sign * total * sin(theta)
  } else if (model == SINGLE) {
    q_use <- pmin(0.99, z$q * load)
    hit <- runif(n) < (1 - (1 - q_use)^G)
    D <- D + hit * z$trait_sign * z$shift_DEL
    U <- U + hit * z$trait_sign * z$shift_DUP
  }
  if (include_sampling) { D <- D + rnorm(n, 0, seD); U <- U + rnorm(n, 0, seU) }
  tibble(Trait = trait, Model = model, DEL = D, DUP = U, Gene_count = G, SE_DEL = seD, SE_DUP = seU)
}
directional_residual <- function(dat, trait) { z <- get_par(trait, MIRROR); z$trait_sign * (dat$DUP + z$rho * dat$DEL) / sqrt(1 + z$rho^2) }
loads <- tibble(Load = factor(c("Low", "Fitted", "High"), levels = c("Low", "Fitted", "High")), multiplier = c(0.25,1,2))
cloud <- bind_rows(lapply(c("Height", "BMI"), function(trait) bind_rows(lapply(c(CONT, SINGLE), function(model) bind_rows(lapply(seq_len(nrow(loads)), function(i) simulate_trait(trait, model, 900, loads$multiplier[i], seed = 20260728 + match(trait,c("Height","BMI"))*1000 + match(model,c(CONT,SINGLE))*100 + i, include_sampling = FALSE) %>% mutate(Load = loads$Load[i], Load_multiplier = loads$multiplier[i])))))))
write_csv(cloud, file.path(root, "outputs_R", "simulated_clouds_R.csv"))
axis_limits <- list(Height = list(x = c(-1.45, 0.65), y = c(-0.75, 0.65)), BMI = list(x = c(-0.85, 1.35), y = c(-0.95, 0.85)))
base_theme <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold", size = 10), strip.text = element_text(size = 8), legend.title = element_blank())
observed_plot <- function(trait) { dat <- if (trait == "Height") height else bmi; z <- get_par(trait, MIRROR); lim <- axis_limits[[trait]]; p <- ggplot(dat, aes(Deletion_beta, Duplication_beta)) + geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) + geom_abline(slope = -z$rho, intercept = 0, linetype = "dashed") + geom_errorbar(aes(ymin = Duplication_beta - Duplication_SE, ymax = Duplication_beta + Duplication_SE), width = 0, alpha = 0.45) + geom_errorbarh(aes(xmin = Deletion_beta - Deletion_SE, xmax = Deletion_beta + Deletion_SE), height = 0, alpha = 0.45) + coord_cartesian(xlim = lim$x, ylim = lim$y) + labs(title = paste("Observed", tolower(trait), "effects"), x = "DEL effect (SD)", y = "DUP effect (SD)") + base_theme; if (trait == "Height") p <- p + geom_point(aes(shape = Dose_response_class), size = 2) else p <- p + geom_point(size = 2); p }
cloud_plot <- function(trait, model, title) { z <- get_par(trait, MIRROR); lim <- axis_limits[[trait]]; d <- cloud %>% filter(Trait == trait, Model == model) %>% mutate(Load_label = factor(paste0(Load, "
(", Load_multiplier, "x)"), levels = paste0(loads$Load, "
(", loads$multiplier, "x)"))); means <- d %>% group_by(Load_label) %>% summarize(DEL = mean(DEL), DUP = mean(DUP), .groups = "drop"); ggplot(d, aes(DEL, DUP)) + geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) + geom_abline(slope = -z$rho, intercept = 0, linetype = "dashed") + geom_point(alpha = 0.12, size = 0.55) + geom_point(data = means, shape = 4, size = 3, stroke = 1.1) + facet_wrap(~Load_label, nrow = 1) + coord_cartesian(xlim = lim$x, ylim = lim$y) + labs(title = title, x = "DEL effect (SD)", y = "DUP effect (SD)") + base_theme + theme(legend.position = "none") }
density_plot <- function(trait) { dat <- if (trait == "Height") height else bmi; obs <- tibble(DEL = dat$Deletion_beta, DUP = dat$Duplication_beta); d <- bind_rows(tibble(Model = "Observed Sacks loci", Residual = directional_residual(obs, trait)), bind_rows(lapply(c(MIRROR, CONT, SINGLE), function(model) { x <- simulate_trait(trait, model, 15000, 1, seed = 20265728 + match(trait, c("Height", "BMI")) * 100 + match(model, c(MIRROR, CONT, SINGLE))); tibble(Model = model, Residual = directional_residual(x, trait)) }))); write_csv(d, file.path(root, "outputs_R", paste0(tolower(trait), "_directional_residuals_R.csv"))); ggplot(d, aes(Residual, linetype = Model)) + geom_density(linewidth = 0.9, adjust = 1) + geom_vline(xintercept = 0, linetype = "dotted") + geom_rug(data = d %>% filter(Model == "Observed Sacks loci"), sides = "b") + labs(title = "Best-fit displacement distributions", x = paste0("Directional distance from mirror line (SD)
Positive = ", ifelse(trait == "Height", "lower-left", "upper-right")), y = "Density") + base_theme }
height_row <- observed_plot("Height") | cloud_plot("Height", CONT, "Continuous: increasing small-effect load") | cloud_plot("Height", SINGLE, "Single hit: increasing hit frequency") | density_plot("Height")
bmi_row <- observed_plot("BMI") | cloud_plot("BMI", CONT, "Continuous: increasing small-effect load") | cloud_plot("BMI", SINGLE, "Single hit: increasing hit frequency") | density_plot("BMI")
full <- (height_row / bmi_row) + plot_layout(widths = c(1.15,2.45,2.45,1.45), guides = "collect") + plot_annotation(title = "Continuous and sparse single-hit models generate distinct patterns of reciprocal CNV effects", subtitle = "BMI mirror clouds use a moment-matched aggregate approximation to avoid zero-inflated artifacts from weakly identifiable gene-level p and m.", tag_levels = "A")
ggsave(file.path(root, "figures_R", "Figure_main_continuous_vs_single_hit.pdf"), full, width = 20, height = 10, units = "in")
ggsave(file.path(root, "figures_R", "Figure_main_continuous_vs_single_hit.png"), full, width = 20, height = 10, units = "in", dpi = 300)
