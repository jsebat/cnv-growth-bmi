# Main-figure panels H and I (head circumference by age and BMI, 16p11.2 del and dup) restricted to ages < 20 (MS 2026-10-13).
# Same as 02b_cross_sectional_split.R / 14_within_vs_controls.R except the age window is 2 <= age < 20 instead of 2-40:
#   visits at age >= 20 are dropped BEFORE taking each person's latest visit, so adults contribute their latest
#   childhood/adolescent visit if they have one (set LATEST_OVERALL <- TRUE to instead keep only people whose latest visit is < 20).
# Controls for the CNV x age x BMI test are restricted to the same window.
# All p-values are recomputed in the < 20 subset (within-group age x BMI, CNV x age x BMI vs controls, high - low BMI offset).
# BMI split (plotted groups and offset test): DEL <= 90th vs > 90th centile, DUP <= 30th vs >= 70th; the age x BMI and vs-controls tests use continuous BMI.
# Output: hc_16p_under20_tests.csv, fig_hc_16p_under20.png

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
AGE_MAX <- 20
LATEST_OVERALL <- FALSE
oc <- "HEADCIRCUM"; OUTLAB <- "Head circumference"
PANELS <- c(H = "16p11.2 del", I = "16p11.2 dup")
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))

d <- tab[!is.na(get(oc)) & !is.na(BMI) & !is.na(age) & age >= 2]
if (!LATEST_OVERALL) d <- d[age < AGE_MAX]
d <- d[order(id, -age)][, .SD[1], by = id]
if (LATEST_OVERALL) d <- d[age < AGE_MAX]
d[, y := get(oc)]

# tests (as in 14_within_vs_controls.R and 02b test_2way)
within_p <- function(x, meth) {
  m <- tryCatch(gam(y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6)), data = x, method = meth), error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m)$s.table["ti(age,BMI)", "p-value"]
}
threeway_p <- function(x, meth) {
  x <- copy(x)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- tryCatch(gam(y ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) + ti(age, BMI, k = c(10, 10)) +
                      ti(age, BMI, k = c(10, 10), by = is_cnv), data = x, method = meth), error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]
}
offset_test <- function(dc) {
  dc <- copy(dc); dc[, grp := as.ordered(bmi_grp)]
  m <- gam(y ~ grp + s(age, k = 10) + s(age, k = 10, by = grp), data = dc, method = "REML")
  c(offset = unname(summary(m)$p.coeff[2]), offset_p = unname(summary(m)$p.pv[2]))
}

# BMI split per locus (centile cuts): low group <= lo, high group > hi (if hi == lo) or >= hi
SPLIT <- list("16p11.2 del" = c(lo = 90, hi = 90), "16p11.2 dup" = c(lo = 30, hi = 70))   # MS 2026-10-13: DEL lowest 90% vs top 10%; DUP 30/70 tails
# vs-controls GAMs use continuous BMI and are slow: cache them by age window so split/plot changes rerun in seconds
CACHE <- file.path(OUT, sprintf("hc_16p_threeway_cache_age%d%s.rds", AGE_MAX, if (LATEST_OVERALL) "_latestoverall" else ""))
tw <- if (file.exists(CACHE)) readRDS(CACHE) else list()
res <- list(); plots <- list()
for (P in names(PANELS)) {
  L <- PANELS[[P]]
  lo <- SPLIT[[L]][["lo"]]; hi <- SPLIT[[L]][["hi"]]; gt <- if (hi == lo) ">" else "\u2265"
  lv <- c(sprintf("BMI \u2264 %dth centile", lo), sprintf("BMI %s %dth centile", gt, hi))
  x <- d[locus == L]
  dc <- copy(x)[, bmi_grp := factor(ifelse(if (hi == lo) BMI > qnorm(hi / 100) else BMI >= qnorm(hi / 100), lv[2], ifelse(BMI <= qnorm(lo / 100), lv[1], NA)), levels = lv)][!is.na(bmi_grp)]
  o <- offset_test(dc)
  if (is.null(tw[[L]])) { tw[[L]] <- c(GCV = threeway_p(d[locus %in% c(L, "No CNV")], "GCV.Cp"), REML = threeway_p(d[locus %in% c(L, "No CNV")], "REML")); saveRDS(tw, CACHE) }
  r <- data.table(outcome = OUTLAB, locus = L, split = sprintf("<=%d / %s%d", lo, gt, hi), age_window = sprintf("2-<%d", AGE_MAX), n_carriers_all_bmi = nrow(x), n_plotted = nrow(dc),
                  n_high_bmi = dc[bmi_grp == lv[2], .N], n_low_bmi = dc[bmi_grp == lv[1], .N], n_controls = d[locus == "No CNV", .N],
                  within_GCV = within_p(x, "GCV.Cp"), within_REML = within_p(x, "REML"),
                  threeway_GCV = tw[[L]][["GCV"]],
                  threeway_REML = tw[[L]][["REML"]],
                  high_minus_low_offset = o[["offset"]], p_offset = o[["offset_p"]])
  res[[P]] <- r
  lab <- paste0(sprintf("age \u00d7 BMI p = %s (GCV) / %s (REML)", pfx(r$within_GCV), pfx(r$within_REML)),
                sprintf("\nvs controls p = %s (GCV) / %s (REML)", pfx(r$threeway_GCV), pfx(r$threeway_REML)),
                sprintf("\nhigh \u2212 low BMI offset %+.2f, p = %s", r$high_minus_low_offset, pfx(r$p_offset)))
  plots[[P]] <- ggplot(dc, aes(age, y, colour = bmi_grp, fill = bmi_grp)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 1.1, alpha = 0.55, shape = 16) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = AGE_MAX, y = -4.6, label = lab, hjust = 1, vjust = 0, size = 3.7) +
    scale_colour_manual(values = c("#00A5AD", "#F0645A")) + scale_fill_manual(values = c("#00A5AD", "#F0645A")) +
    guides(colour = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(2, AGE_MAX), ylim = c(-4.6, 4.6)) +
    labs(x = "Age (years)", y = paste(OUTLAB, "(SD)"), colour = NULL, fill = NULL,
         subtitle = sprintf("%s, ages 2-%d (n = %d)", L, AGE_MAX - 1, nrow(dc))) +
    theme_classic(base_size = 15) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.position.inside = c(0.75, 0.93),
                                          legend.background = element_blank(), legend.key.height = unit(8, "pt"))
}
res <- rbindlist(res)
fwrite(res, file.path(OUT, "hc_16p_under20_tests.csv"))
print(res, digits = 3)

fig <- wrap_plots(plots, ncol = 1) + plot_annotation(tag_levels = list(names(PANELS))) & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(OUT, "fig_hc_16p_under20.png"), fig, width = 6, height = 8.2, dpi = 200)
