# Decision figure (JS 2026-10-13): within-group age x BMI interaction, controls alone vs each CNV, both GAM
# smoothing-parameter methods, plus the carrier-vs-control 3-way test (both methods).
#
# Data: one record per person (latest visit with the outcome and BMI), ages 2-40, same table as the main figure.
#   Controls = non-carriers (SSC + SVIP; SSC "No CNV" includes probands, siblings and parents).
# Within-group test (continuous BMI):  y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6));  p of ti(age,BMI)
#   = does the BMI effect on the outcome change with age, within that group.
# 3-way (Anjali's test): y ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) + ti(age, BMI) +
#   ti(age, BMI, by = is_cnv);  p of ti(age,BMI):is_cnvCNV = does the carriers' age x BMI surface differ from controls'.
# Each fitted with method = "GCV.Cp" (mgcv default, Anjali) and "REML".
# Curves: GAM by BMI stratum (>= 70th vs <= 30th centile), as in the main figure.
# Output: within_vs_controls_tests.csv, fig_within_vs_controls.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
GROUPS <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
OUT <- c(HEIGHT = "Height", HEADCIRCUM = "Head circumference")
pf <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))

within_p <- function(x, meth) {
  m <- tryCatch(gam(y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6)), data = x, method = meth), error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m)$s.table["ti(age,BMI)", "p-value"]
}
threeway_p <- function(x, meth) {
  x <- copy(x)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- gam(y ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) + ti(age, BMI, k = c(10, 10)) +
             ti(age, BMI, k = c(10, 10), by = is_cnv), data = x, method = meth)
  summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]
}

res <- list(); plots <- list()
for (oc in names(OUT)) {
  d <- tab[!is.na(tab[[oc]]) & !is.na(BMI) & age >= 2 & age <= 40][order(id, -age)][, .SD[1], by = id]
  d[, y := d[[oc]]]
  for (G in GROUPS) {
    x <- d[locus == G]
    r <- data.table(outcome = OUT[[oc]], group = G, n = nrow(x),
                    within_GCV = within_p(x, "GCV.Cp"), within_REML = within_p(x, "REML"),
                    threeway_GCV = if (G == "No CNV") NA_real_ else threeway_p(d[locus %in% c(G, "No CNV")], "GCV.Cp"),
                    threeway_REML = if (G == "No CNV") NA_real_ else threeway_p(d[locus %in% c(G, "No CNV")], "REML"))
    res[[length(res) + 1]] <- r; cat(oc, G, "done\n")
    lv <- c("BMI \u2265 70th centile", "BMI \u2264 30th centile")
    x[, grp := factor(ifelse(BMI >= qnorm(.7), lv[1], ifelse(BMI <= qnorm(.3), lv[2], NA)), levels = lv)]
    xs <- x[!is.na(grp)]
    lab <- sprintf("within-group age \u00d7 BMI p = %s (GCV) / %s (REML)%s", pf(r$within_GCV), pf(r$within_REML),
                   if (G == "No CNV") "" else sprintf("\nCNV \u00d7 age \u00d7 BMI vs controls p = %s (GCV) / %s (REML)", pf(r$threeway_GCV), pf(r$threeway_REML)))
    plots[[paste(oc, G)]] <- ggplot(xs, aes(age, y, colour = grp, fill = grp)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_point(size = if (G == "No CNV") 0.5 else 1.1, alpha = if (G == "No CNV") 0.2 else 0.55, shape = 16) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
      annotate("text", x = 40, y = -4.6, label = lab, hjust = 1, vjust = 0, size = 3.6) +
      scale_colour_manual(values = c("#F0645A", "#00A5AD")) + scale_fill_manual(values = c("#F0645A", "#00A5AD")) +
      coord_cartesian(xlim = c(2, 40), ylim = c(-4.6, 4.6)) +
      labs(x = "Age (years)", y = paste(OUT[[oc]], "(SD)"), colour = NULL, fill = NULL,
           title = if (G == GROUPS[1]) OUT[[oc]] else NULL,
           subtitle = sprintf("%s (n = %d)", if (G == "No CNV") "Non-carrier controls (SSC + SVIP)" else G, nrow(x))) +
      theme_classic(base_size = 15) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.position.inside = c(0.8, 0.93),
                                            legend.background = element_blank(), legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
  }
}
res <- rbindlist(res); fwrite(res, file.path(HERE, "within_vs_controls_tests.csv"))
print(res[, .(outcome, group, n, within_GCV = signif(within_GCV, 2), within_REML = signif(within_REML, 2),
              threeway_GCV = signif(threeway_GCV, 2), threeway_REML = signif(threeway_REML, 2))])
ord <- c(paste("HEIGHT", GROUPS), paste("HEADCIRCUM", GROUPS))
fig <- wrap_plots(plots[ord], ncol = 2, byrow = FALSE) +
  plot_annotation(theme = theme(plot.title = element_text(size = 14), plot.caption = element_text(size = 9)),
                  title = "Age-dependent BMI effects within each group: non-carrier controls vs CNV carriers",
                  caption = paste("Within-group age \u00d7 BMI: y ~ s(age) + s(BMI) + ti(age, BMI), continuous BMI, p of ti(age, BMI). CNV \u00d7 age \u00d7 BMI: carriers vs controls (Anjali's tensor-product test).",
                                  "\nGCV = mgcv default smoothing selection (Anjali); REML = restricted maximum likelihood. Curves: BMI \u2265 70th vs \u2264 30th centile. One record per person, ages 2-40."))
ggsave(file.path(HERE, "fig_within_vs_controls.png"), fig, width = 13, height = 21, dpi = 170)
ggsave(file.path(HERE, "fig_within_vs_controls.pdf"), fig, width = 13, height = 21, device = cairo_pdf)
