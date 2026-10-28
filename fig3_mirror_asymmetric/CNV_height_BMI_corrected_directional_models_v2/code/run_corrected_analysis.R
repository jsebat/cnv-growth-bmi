#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", flag)
  args[hit + 1]
}
project_dir <- normalizePath(get_arg("--project-dir", "."), mustWork = FALSE)
run_loo <- !"--skip-loo" %in% args
formats <- strsplit(get_arg("--formats", "html,word,pdf"), ",", fixed = TRUE)[[1]]
seed <- as.integer(get_arg("--seed", "20260726"))

source(file.path(project_dir, "code", "corrected_model_functions.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

read_csv_base <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
dir.create(file.path(project_dir, "R_outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "R_outputs", "height"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "R_outputs", "bmi"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "R_outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

height <- read_csv_base(file.path(project_dir, "height", "inputs", "sacks_height_all_loci_exact.csv"))
height_mirror <- height[tolower(height$Dose_response_class) == "mirror", , drop = FALSE]
bmi <- read_csv_base(file.path(project_dir, "bmi", "inputs", "sacks_bmi_complete_reciprocal_loci_exact.csv"))

message("Fitting corrected height models...")
H <- fit_trait_models(height, height_mirror, pi_positive = 0.30,
                      trait_sign = -1, seed = seed, run_loo = run_loo)
message("Fitting corrected BMI models...")
B <- fit_trait_models(bmi, bmi, pi_positive = 11/14,
                      trait_sign = +1, seed = seed + 1000, run_loo = run_loo)

write.csv(as.data.frame(H$base[setdiff(names(H$base), "raw_par")]),
          file.path(project_dir, "R_outputs", "height", "height_mirror_parameters.csv"), row.names = FALSE)
write.csv(H$comparison,
          file.path(project_dir, "R_outputs", "height", "height_single_locus_model_comparison.csv"), row.names = FALSE)
write.csv(as.data.frame(B$base[setdiff(names(B$base), "raw_par")]),
          file.path(project_dir, "R_outputs", "bmi", "bmi_mirror_parameters.csv"), row.names = FALSE)
write.csv(B$comparison,
          file.path(project_dir, "R_outputs", "bmi", "bmi_single_locus_model_comparison.csv"), row.names = FALSE)
if (run_loo) {
  write.csv(H$loo$detail, file.path(project_dir, "R_outputs", "height", "height_conditional_LOO.csv"), row.names = FALSE)
  write.csv(B$loo$detail, file.path(project_dir, "R_outputs", "bmi", "bmi_conditional_LOO.csv"), row.names = FALSE)
}

predictions <- function(dat, fit_obj, trait_sign, prefix, outdir) {
  # File suffixes must match the Python reference outputs: the sparse model is
  # keyed "single" internally but written as "single_hit".
  suffix <- c(mirror = "mirror", continuous = "continuous", single = "single_hit")
  for (model in c("mirror", "continuous", "single")) {
    z <- model_locus_expectations(dat, fit_obj$base, trait_sign, model, fit_obj$fits[[model]])
    write.csv(z, file.path(outdir, paste0(prefix, "_locus_predictions_", suffix[[model]], ".csv")), row.names = FALSE)
  }
}
predictions(height, H, -1, "height", file.path(project_dir, "R_outputs", "height"))
predictions(bmi, B, +1, "bmi", file.path(project_dir, "R_outputs", "bmi"))

burden_table <- function(dat, fit_obj, trait_sign, target, native_sd, units) {
  labels <- c(mirror = "Mirror alone", continuous = "Mirror + continuous", single = "Mirror + single hit")
  do.call(rbind, lapply(names(labels), function(model) {
    z <- burden_values(fit_obj$base, trait_sign, model, fit_obj$fits[[model]], dat$Gene_count)
    d <- mean(z$DEL_SD_per_gene); u <- mean(z$DUP_SD_per_gene)
    d_native <- d * native_sd; u_native <- u * native_sd
    data.frame(
      model = labels[[model]], DEL_effect_per_affected_gene_SD = d,
      DUP_effect_per_affected_gene_SD = u,
      DEL_effect_per_affected_gene_native = d_native,
      DUP_effect_per_affected_gene_native = u_native,
      Mace_DEL_target_per_Mb = target, Mace_units = units,
      effective_dosage_genes_per_Mb_to_match_DEL = target / d_native,
      predicted_DUP_effect_per_Mb_at_matched_density = target * u / d
    )
  }))
}
Hbur <- burden_table(height, H, -1, -0.41, 7.99, "cm/Mb")
Bbur <- burden_table(bmi, B, +1, 0.14, 4.7, "kg/m2/Mb")
write.csv(Hbur, file.path(project_dir, "R_outputs", "height", "height_burden_predictions.csv"), row.names = FALSE)
write.csv(Bbur, file.path(project_dir, "R_outputs", "bmi", "bmi_burden_predictions.csv"), row.names = FALSE)

# Height independent evidence: size and within-class gradient.
hm <- height[tolower(height$Dose_response_class) == "mirror", ]
ha <- height[tolower(height$Dose_response_class) == "asymmetric", ]
size_rows <- lapply(c("CNV_length_Mb", "Gene_count"), function(v) {
  trans <- if (v == "CNV_length_Mb") log else log1p
  data.frame(
    variable = v, asymmetric_mean = mean(ha[[v]]), mirror_mean = mean(hm[[v]]),
    asymmetric_median = median(ha[[v]]), mirror_median = median(hm[[v]]),
    Welch_p = t.test(trans(ha[[v]]), trans(hm[[v]]))$p.value,
    Mann_Whitney_p = wilcox.test(ha[[v]], hm[[v]], exact = FALSE)$p.value
  )
})
write.csv(do.call(rbind, size_rows), file.path(project_dir, "R_outputs", "height", "height_size_gene_tests.csv"), row.names = FALSE)

hres <- model_locus_expectations(height, H$base, -1, "mirror", H$fits$mirror)
grad <- list(); kk <- 1
for (subset_name in c("all", "mirror", "asymmetric")) {
  z <- if (subset_name == "all") hres else hres[tolower(hres$Dose_response_class) == subset_name, ]
  for (v in c("CNV_length_Mb", "Gene_count")) {
    tst <- cor.test(z$directional_shift_score, log1p(z[[v]]), method = "spearman", exact = FALSE)
    grad[[kk]] <- data.frame(subset = subset_name, variable = v,
                             Spearman_rho = unname(tst$estimate), p_value = tst$p.value)
    kk <- kk + 1
  }
}
write.csv(do.call(rbind, grad), file.path(project_dir, "R_outputs", "height", "height_within_class_gradient_tests.csv"), row.names = FALSE)

# Core figures regenerated in R.
model_cols <- c("Mirror alone" = "grey45",
                "Mirror + continuous directional shift" = "#4f8a5b",
                "Mirror + sparse single hit" = "#7563a8")
make_scatter <- function(dat, fit_obj, trait_sign, trait, out_path) {
  allpred <- bind_rows(lapply(c("mirror", "continuous", "single"), function(m) {
    z <- model_locus_expectations(dat, fit_obj$base, trait_sign, m, fit_obj$fits[[m]])
    z$model <- fit_obj$fits[[m]]$model
    z
  }))
  p <- ggplot(allpred, aes(Deletion_beta, Duplication_beta)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_segment(aes(x = model_expected_DEL, y = model_expected_DUP,
                     xend = Deletion_beta, yend = Duplication_beta), colour = "grey85", linewidth = 0.25) +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_point(aes(x = model_expected_DEL, y = model_expected_DUP, colour = model), size = 1.7) +
    facet_wrap(~model, nrow = 1) + scale_colour_manual(values = model_cols) +
    labs(title = paste0(trait, ": corrected single-locus model comparison"),
         x = "Deletion effect (SD)", y = "Duplication effect (SD)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")
  ggsave(out_path, p, width = 11, height = 4.1, dpi = 300)
}
make_scatter(height, H, -1, "Height", file.path(project_dir, "R_outputs", "figures", "R_height_single_locus_comparison.png"))
make_scatter(bmi, B, +1, "BMI", file.path(project_dir, "R_outputs", "figures", "R_BMI_single_locus_comparison.png"))

make_burden_plot <- function(tab, trait, out_path) {
  z <- tab %>% select(model, DEL_effect_per_affected_gene_SD, DUP_effect_per_affected_gene_SD) %>%
    pivot_longer(-model, names_to = "copy_state", values_to = "effect")
  p <- ggplot(z, aes(model, effect, fill = copy_state)) + geom_col(position = "dodge") +
    geom_hline(yintercept = 0) + scale_fill_manual(values = c("#b74d45", "#286090"),
      labels = c("DEL", "DUP")) +
    labs(title = paste0(trait, ": propagation to genome-wide burden"), x = NULL,
         y = "Expected effect per affected gene (SD)", fill = NULL) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(out_path, p, width = 7, height = 4.5, dpi = 300)
}
make_burden_plot(Hbur, "Height", file.path(project_dir, "R_outputs", "figures", "R_height_burden.png"))
make_burden_plot(Bbur, "BMI", file.path(project_dir, "R_outputs", "figures", "R_BMI_burden.png"))

# Render dynamic report when rmarkdown is available.
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  out_map <- c(html = "html_document", word = "word_document", pdf = "pdf_document")
  for (fmt in formats) {
    if (!fmt %in% names(out_map)) warning("Unknown format: ", fmt) else {
      try(rmarkdown::render(
        file.path(project_dir, "code", "CNV_corrected_models_report.Rmd"),
        output_format = out_map[[fmt]],
        output_file = paste0("CNV_corrected_models_R_report.",
                            c(html = "html", word = "docx", pdf = "pdf")[[fmt]]),
        output_dir = file.path(project_dir, "R_outputs"),
        params = list(project_dir = project_dir), envir = new.env(parent = globalenv()), quiet = TRUE
      ), silent = FALSE)
    }
  }
}
message("Completed. Outputs: ", file.path(project_dir, "R_outputs"))
