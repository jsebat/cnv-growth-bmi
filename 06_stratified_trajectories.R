# BMI-stratified growth trajectories using every visit (longitudinal + single-visit carriers), age 2-40 y.
# Split rule (JS 2026-10-12): obese loci (16p11.2 del) -> BMI > 90th centile vs the rest;
#   loci that are not obese on average -> tails, BMI <= 30th vs >= 70th centile (middle 40% dropped).
# A person's stratum = their mean BMI Z across visits (so a child does not switch curve between visits).
# Model per locus x outcome (mgcv, REML):
#   Z ~ strat + sex + cohort + s(age, k = 10) + s(age, by = strat_ord, k = 10) + s(id, bs = "re")
#   difference-smooth p = do the two age curves differ in shape; offset p = are they shifted overall.
# Visits cleaned as in 00_helpers.R (duplicates, carried-forward heights, impossible changes, > 2.5 SD jumps).
# Z scale changes at 18 y (childsds growth references below, within-sample Box-Cox residual above) -> dotted line.
# Output: fig_stratified_trajectories.{png,pdf}, stratified_trajectory_tests.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
source(file.path(HERE, "00_helpers.R"))
SPLIT <- c("16p11.2 del" = "30_70", "16p11.2 dup" = "30_70", "22q11.2 del" = "30_70", "22q11.2 dup" = "30_70")   # JS 2026-10-12: 30/70 tails for every locus (was > 90th vs rest for 16p11.2 del)
OUT <- c(HEIGHT = "Height", HEADCIRCUM = "Head circumference")
COL <- c("#00A5AD", "#F0645A")
SPAGHETTI <- Sys.getenv("SPAGHETTI", "0") == "1"   # faint lines joining repeat visits (busy; supplementary version)
SUF <- if (SPAGHETTI) "_spaghetti" else ""

prep <- function(oc, L) {
  d <- tab[locus == L & age >= 2 & age <= 40 & !is.na(tab[[oc]]) & !is.na(BMI)][order(id, age)]
  d <- clean_visits(d, oc); d <- unique(d, by = c("id", "age"))
  d[, y := d[[oc]]]; d[, bmi_mean := mean(BMI), by = id]
  if (SPLIT[[L]] == "p90") {
    lv <- c("BMI \u2264 90th centile", "BMI > 90th centile")
    d[, strat := factor(ifelse(bmi_mean > qnorm(.9), lv[2], lv[1]), levels = lv)]
  } else {
    lv <- c("BMI \u2264 30th centile", "BMI \u2265 70th centile")
    d[, strat := factor(ifelse(bmi_mean >= qnorm(.7), lv[2], ifelse(bmi_mean <= qnorm(.3), lv[1], NA)), levels = lv)]
    d <- d[!is.na(strat)]
  }
  d[, `:=`(id = factor(id), cohort = factor(cohort), strat_ord = as.ordered(strat), grp = as.integer(strat))]
  d
}

res <- list(); plots <- list()
for (oc in names(OUT)) for (L in LOCI) {
  d <- prep(oc, L)
  cterm <- if (nlevels(d$cohort) > 1) "+ cohort" else ""
  m <- gam(as.formula(sprintf("y ~ strat_ord + sex %s + s(age, k = 10) + s(age, by = strat_ord, k = 10) + s(id, bs = 're')", cterm)),
           data = d, method = "REML")
  s <- summary(m)
  r <- data.table(outcome = OUT[oc], locus = L, split = SPLIT[[L]], persons = uniqueN(d$id), visits = nrow(d),
                  persons_multi = d[, .N, by = id][N > 1, .N], persons_high = d[grp == 2, uniqueN(id)],
                  offset = s$p.coeff[2], p_offset = s$p.pv[2], p_curve_shape = s$s.table[grep("strat_ord", rownames(s$s.table)), "p-value"])
  res[[length(res) + 1]] <- r
  # population curves: average over sex and cohorts present, random effect excluded
  nd <- CJ(age = seq(2, 40, 0.25), strat = levels(d$strat), sex = levels(droplevels(d$sex)), cohort = levels(d$cohort))
  rng <- d[, .(lo = min(age), hi = max(age)), by = .(strat = as.character(strat))]
  nd <- merge(nd, rng, by = "strat")[age >= lo & age <= hi][, c("lo", "hi") := NULL]
  nd[, `:=`(strat_ord = factor(strat, levels = levels(d$strat), ordered = TRUE), id = d$id[1], cohort = factor(cohort, levels(d$cohort)))]
  pr <- predict(m, nd, exclude = "s(id)", se.fit = TRUE); nd[, `:=`(fit = pr$fit, se = pr$se.fit)]
  cur <- nd[, .(fit = mean(fit), se = mean(se)), by = .(strat, age)][, strat := factor(strat, levels(d$strat))]
  lab <- sprintf("high \u2212 low offset %+.2f, p = %s\ncurve-shape difference p = %s\n%d persons (%d with repeat visits)",
                 r$offset, format.pval(r$p_offset, 2, eps = 1e-4), format.pval(r$p_curve_shape, 2, eps = 1e-4), r$persons, r$persons_multi)
  plots[[paste(oc, L)]] <- ggplot() +
    annotate("rect", xmin = 8, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey93") +
    geom_vline(xintercept = 18, linetype = 3, colour = "grey50") + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    {if (SPAGHETTI) geom_line(data = d[id %in% d[, .N, by = id][N > 1, id]], aes(age, y, group = id, colour = strat), alpha = 0.25, linewidth = 0.3)} +
    geom_point(data = d, aes(age, y, colour = strat), size = 0.7, alpha = if (SPAGHETTI) 0.45 else 0.3, shape = 16) +
    geom_ribbon(data = cur, aes(age, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se, fill = strat), alpha = 0.25) +
    geom_line(data = cur, aes(age, fit, colour = strat), linewidth = 1) +
    annotate("text", x = 40, y = -4.5, label = lab, hjust = 1, vjust = 0, size = 2.5) +
    scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
    coord_cartesian(xlim = c(2, 40), ylim = c(-4.5, 4.5)) +
    labs(x = "Age (years)", y = paste(OUT[oc], "(SD)"), colour = NULL, fill = NULL,
         title = if (L == LOCI[1]) OUT[oc] else NULL, subtitle = L) +
    theme_classic(base_size = 10) + theme(legend.position = "inside", legend.position.inside = c(0.75, 0.93),
                                           legend.background = element_blank(), legend.key.height = unit(8, "pt"),
                                           plot.title = element_text(face = "bold", hjust = 0.5))
}
res <- rbindlist(res); fwrite(res, file.path(HERE, "stratified_trajectory_tests.csv")); print(res, digits = 3)
ord <- c(rbind(paste("HEIGHT", LOCI), paste("HEADCIRCUM", LOCI)))
fig <- wrap_plots(plots[ord], ncol = 2) +
  plot_annotation(title = "Growth trajectories of 16p11.2 and 22q11.2 carriers by BMI (all visits; G2MH, Cardiff, SVIP, SSC, Searchlight)",
                  subtitle = paste("16p11.2 del: BMI > 90th centile vs rest; other loci: BMI \u2265 70th vs \u2264 30th centile (middle dropped). Stratum = person's mean BMI Z.",
                                   if (SPAGHETTI) "
Faint lines join a child's repeat visits." else "
", "Bold curves = mixed GAM (person random intercept) \u00b1 95% CI. Grey = 8-13 y; dotted = 18 y (Z scale changes)."))
ggsave(file.path(HERE, paste0("fig_stratified_trajectories", SUF, ".png")), fig, width = 11, height = 13, dpi = 200)
ggsave(file.path(HERE, paste0("fig_stratified_trajectories", SUF, ".pdf")), fig, width = 11, height = 13, device = cairo_pdf)
