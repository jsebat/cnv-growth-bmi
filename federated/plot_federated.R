# plot_federated.R  (2026-11-20)
# What: figures from the output of meta_combine.R. No individual data are needed or shown.
#   a. Trajectory panels in the style of Fig. 1a-e / 2a-e of the paper (one panel for controls and one per genotype):
#      the curve of each BMI tail (>= 70th centile red, <= 30th centile teal) with its pointwise 95% band, and the
#      binned means as points sized by n. One figure for the pooled curves and one per site.
#   b. Forest of the CNV x BMI term by age band (rows) and genotype (columns): each site and the pooled estimate.
#   c. Forest of the longitudinal BMI coefficient (change in Z to the next visit per SD of BMI Z), by age window.
# Inputs: the folder written by meta_combine.R (--meta).
# Outputs (in --out): fig_trajectories_<outcome>_<scope>.png/.pdf, fig_forest_bands_<outcome>.png/.pdf,
#   fig_forest_longitudinal.png/.pdf, and <figure>_plotted.csv with every plotted number.
# null model: none; the figures show estimates and 95% intervals from meta_combine.R (z test of the pooled estimate,
#   printed as "p = ...", bold when < 0.05).
# Usage: Rscript plot_federated.R --meta meta_out --out figures
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
.fed_dir <- function() { if (exists("FEDERATED_DIR")) return(FEDERATED_DIR)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else getwd() }
source(file.path(.fed_dir(), "00_federated_common.R"))
POOLED <- "Pooled (fixed effects)"
TLAB <- c("BMI >= 70th centile" = "BMI ≥ 70th centile", "BMI <= 30th centile" = "BMI ≤ 30th centile")
TH <- theme_classic(base_size = 15) + theme(plot.tag = element_text(size = 22, face = "bold"), legend.background = element_blank(), plot.subtitle = element_text(size = 12))
save_fig <- function(fig, f, w, h) { ggsave(paste0(f, ".png"), fig, width = w, height = h, dpi = 200, bg = "white")
  tryCatch(ggsave(paste0(f, ".pdf"), fig, width = w, height = h, device = cairo_pdf), error = function(e) ggsave(paste0(f, ".pdf"), fig, width = w, height = h)) }
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

plot_federated <- function(meta, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rdm <- function(f) { p <- file.path(meta, f); if (file.exists(p)) fread(p, na.strings = c("", "NA")) else NULL }
  CV <- rdm("meta_trajectory_curves.csv"); BN <- rdm("meta_trajectory_bins.csv"); AE <- rdm("all_estimates.csv"); written <- character(0)
  # ---------------------------------------------------------------- a. trajectories
  if (!is.null(CV) && !is.null(BN)) for (o in unique(CV$outcome)) for (sc in unique(CV$scope)) {
    cv <- CV[outcome == o & scope == sc & bmi_tail %in% TAILS]; bn <- BN[outcome == o & scope == sc]; if (!nrow(cv)) next
    cv[, tail := factor(TLAB[bmi_tail], TLAB)]; bn[, tail := factor(TLAB[bmi_tail], TLAB)][, x := ifelse(is.na(mean_age), age_bin_mid, mean_age)]
    yr <- range(c(cv$fit, bn$mean_z, quantile(c(cv$lo, cv$hi), c(0.05, 0.95), na.rm = TRUE)), na.rm = TRUE); yr <- c(min(yr[1], -0.5), max(yr[2], 0.5)); yr <- yr + c(-0.15, 0.3) * diff(yr); xmax <- min(73, max(20, ceiling(max(cv$age))))
    smax <- max(bn$n, na.rm = TRUE)
    panels <- lapply(seq_along(GROUPS5), function(i) { G <- GROUPS5[i]; c1 <- cv[group == G]; b1 <- bn[group == G]
      e <- if (is.null(AE) || G == CTL) NULL else AE[site == sc & analysis == "cross-sectional all ages" & outcome == o & group == G & term == "CNV x BMI" & status == "ok"]
      lab <- if (!is.null(e) && nrow(e)) sprintf("CNV × BMI %s", plab(e$p[1])) else ""; bold <- !is.null(e) && nrow(e) && isTRUE(e$p[1] < 0.05)
      nt <- unique(c1[, .(tail, n)]); ns <- if (nrow(c1)) max(c1$n_sites) else 0
      g <- ggplot() + geom_hline(yintercept = 0, linetype = 2, colour = "grey40")
      if (nrow(c1)) g <- g + geom_ribbon(data = c1, aes(age, ymin = lo, ymax = hi, fill = tail), alpha = 0.2) + geom_line(data = c1, aes(age, fit, colour = tail), linewidth = 1)
      if (nrow(b1)) g <- g + geom_point(data = b1, aes(x, mean_z, colour = tail, size = n), alpha = 0.75)
      g + annotate("text", x = xmax, y = yr[1], label = lab, hjust = 1, vjust = 0, size = 3.8, fontface = if (bold) "bold" else "plain") +
        scale_colour_manual(values = setNames(unname(TAIL_COLS), TLAB), breaks = TLAB, drop = FALSE, name = NULL) + scale_fill_manual(values = setNames(unname(TAIL_COLS), TLAB), breaks = TLAB, drop = FALSE, guide = "none") +
        scale_size_area(max_size = 5, limits = c(0, smax), guide = "none") + coord_cartesian(xlim = c(2, xmax), ylim = yr) +
        labs(x = "Age (years)", y = if (i == 1) paste0(o, " (SD)") else NULL, tag = letters[i],
             subtitle = sprintf("%s\nn = %s (BMI tails)%s", if (G == CTL) "Non-carrier controls" else G, if (nrow(nt)) format(sum(nt$n), big.mark = ",") else "0",
                                if (sc == POOLED && ns > 0) sprintf("; %d site%s", ns, if (ns > 1) "s" else "") else "")) +
        TH + theme(legend.position = if (i == 1) "inside" else "none", legend.position.inside = c(0.62, 0.93), legend.text = element_text(size = 9), legend.key.height = unit(9, "pt")) })
    fig <- wrap_plots(panels, nrow = 1) + plot_annotation(title = sprintf("%s by BMI tail: %s", o, sc), subtitle = "Curves: fixed natural-spline basis of log(age), 95% band. Points: binned means, area = n. No individual values.",
                                                         theme = theme(plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 12)))
    f <- file.path(out, sprintf("fig_trajectories_%s_%s", safe(o), safe(sc))); save_fig(fig, f, 20, 5.4)
    fwrite(rbind(cbind(what = "curve", cv[, .(scope, outcome, group, bmi_tail, age, fit, se, lo, hi, n, n_sites)]), cbind(what = "binned mean", bn[, .(scope, outcome, group, bmi_tail, age = x, fit = mean_z, se, n)]), fill = TRUE), paste0(f, "_plotted.csv"))
    written <- c(written, paste0(f, ".png"))
  }
  # ---------------------------------------------------------------- forests
  forest <- function(d, title, xlab, f, row_levels, w = 20, h = 9) { if (!nrow(d)) return(invisible(NULL))
    sites <- c(sort(setdiff(unique(d$site), POOLED)), POOLED)
    d <- d[status == "ok" & is.finite(estimate) & is.finite(se)]; if (!nrow(d)) return(invisible(NULL))
    d[, `:=`(site = factor(site, rev(sites)), group = factor(group, intersect(c(CTL, LOCI, POOL), unique(group))), stratum = factor(stratum, intersect(row_levels, unique(stratum))), pooled = site == POOLED)]
    d[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]; xr <- range(c(d$lo, d$hi, 0)); xr <- xr + c(-0.05, 0.45) * diff(xr)
    d[, lab := ifelse(pooled, plab(p), "")][, nlab := ifelse(is.na(n_carriers), "", ifelse(grepl("^>=", n_carriers), "n masked", paste0("n = ", n_carriers)))]
    g <- ggplot(d, aes(estimate, site, colour = pooled)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8) + geom_point(aes(shape = pooled, size = pooled)) +
      geom_text(aes(x = xr[2], label = ifelse(pooled, lab, nlab), fontface = ifelse(pooled & !is.na(p) & p < 0.05, "bold", "plain")), hjust = 1, size = 3.3, colour = "grey15") +
      scale_colour_manual(values = c(`FALSE` = COL_CTL, `TRUE` = COL_CAR), guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") + scale_size_manual(values = c(`FALSE` = 2.6, `TRUE` = 4.6), guide = "none") +
      facet_grid(stratum ~ group, drop = TRUE) + coord_cartesian(xlim = xr) + labs(x = xlab, y = NULL, title = title, subtitle = "Each site and the fixed-effects pooled estimate, 95% interval. p: pooled estimate.") +
      theme_classic(base_size = 15) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 13), strip.text.y = element_text(angle = 0), panel.border = element_rect(fill = NA, colour = "grey70"),
                                             plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 12))
    save_fig(g, f, w, h); fwrite(d[, .(site, group, stratum, estimate, se, lo, hi, p, n_carriers, n_controls, n_pairs)], paste0(f, "_plotted.csv")); written <<- c(written, paste0(f, ".png")) }
  if (!is.null(AE)) {
    for (o in intersect(OUTCOMES, unique(AE$outcome))) {
      d <- AE[analysis == "cross-sectional band" & outcome == o & term == "CNV x BMI" & (group %in% LOCI | (group == POOL & grepl("(primary)", model, fixed = TRUE)))]
      forest(copy(d), sprintf("CNV × BMI on %s, by age band", tolower(o)), sprintf("Carrier minus control BMI slope (SD of %s per SD of BMI)", tolower(o)), file.path(out, paste0("fig_forest_bands_", safe(o))), paste0(BLAB), h = 2.2 + 1.9 * uniqueN(d$stratum))
    }
    d <- AE[analysis == "longitudinal" & term == "BMI Z at the visit"]
    if (nrow(d)) { d[, stratum := paste0(ifelse(outcome == "Height", "Height", "Head circ."), "\n", stratum)]
      lv <- c(paste0("Height\n", c("All ages", WIN_H)), paste0("Head circ.\n", c("All ages", WIN_C)))
      forest(copy(d), "Longitudinal: BMI Z at a visit and the change in Z to the next visit", "Change in Z to the next visit per SD of BMI Z", file.path(out, "fig_forest_longitudinal"), lv, h = 2.2 + 1.35 * uniqueN(d$stratum)) }
  }
  message("figures written: ", length(written)); invisible(written)
}

if (!exists("FEDERATED_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); opt <- function(k) { i <- match(k, a); if (is.na(i)) NULL else a[i + 1] }
  if (is.null(opt("--meta")) || is.null(opt("--out"))) stop("usage: Rscript plot_federated.R --meta meta_out --out figures")
  plot_federated(opt("--meta"), opt("--out"))
}
