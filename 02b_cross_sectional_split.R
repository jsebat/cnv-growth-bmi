# Cross-sectional trajectories: height and head circumference Z by age, BMI >= 90th centile vs rest, per CNV.
# One record per person (latest visit with the outcome and BMI), ages 2-40 y.
# Tests per CNV x outcome:
#   3-way  CNV x age x BMI   -- Anjali's tensor-product GAM (continuous BMI, SSC/SVIP/G2MH non-carrier controls):
#          Z ~ is_cnv + s(age) + s(BMI) + s(age, by = is_cnv) + ti(age, BMI) + ti(age, BMI, by = is_cnv)
#   2-way  age x BMI stratum within carriers -- ordered-factor difference smooth (does the age curve differ)
#          plus the parametric stratum offset (is the high-BMI curve shifted overall)
# 02b: same figure with the locus-specific BMI split (16p11.2 del > 90th vs rest; other loci >= 70th vs <= 30th, middle dropped),
# reusing the CNV x age x BMI p-values from cross_sectional_tests.csv (continuous BMI, so the split does not change them).
# Output: fig_cross_sectional_split.{png,pdf}, cross_sectional_tests_split.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
Z90 <- qnorm(0.9)
OUT <- c(HEIGHT = "Height", HEADCIRCUM = "Head circumference")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")

SPLIT <- c("16p11.2 del" = "30_70", "16p11.2 dup" = "30_70", "22q11.2 del" = "30_70", "22q11.2 dup" = "30_70", "No CNV" = "30_70")   # JS 2026-10-12: 30/70 tails for every locus (was > 90th vs rest for 16p11.2 del)
P3 <- fread(file.path(HERE, "cross_sectional_tests.csv"))
W <- fread(file.path(HERE, "within_vs_controls_tests.csv"))
# UK Biobank adult anchor (EUR; aggregates from ukbb_anchor.R on Expanse; cells n < 5 suppressed), height only
UK <- fread(file.path(HERE, "ukbb_anchor_aggregates.csv"))[age_bin != "all" & bmi_grp %in% c("high", "low") & !is.na(mean_hZ)]
UK[, bmi_grp := factor(ifelse(bmi_grp == "high", "BMI ≥ 70th centile", "BMI ≤ 30th centile"), levels = c("BMI ≤ 30th centile", "BMI ≥ 70th centile"))]   # within-group and 3-way p, GCV and REML (14_within_vs_controls.R)
GROUPS <- c("No CNV", LOCI)   # JS 2026-10-13: control row on top
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))
latest_rec <- function(outcome) {
  d <- tab[!is.na(get(outcome)) & !is.na(BMI) & !is.na(age) & age >= 2 & age <= 40]
  d <- d[order(id, -age)][, .SD[1], by = id]
  d[, y := get(outcome)]
  d
}
split_grp <- function(dc, L) {
  if (SPLIT[[L]] == "p90") { lv <- c("BMI \u2264 90th centile", "BMI > 90th centile"); dc[, bmi_grp := factor(ifelse(BMI > Z90, lv[2], lv[1]), levels = lv)] }
  else { lv <- c("BMI \u2264 30th centile", "BMI \u2265 70th centile")
         dc[, bmi_grp := factor(ifelse(BMI >= qnorm(.7), lv[2], ifelse(BMI <= qnorm(.3), lv[1], NA)), levels = lv)] }
  dc[!is.na(bmi_grp)]
}

test_3way <- function(d) {
  d <- copy(d); d[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- gam(y ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) +
             ti(age, BMI, k = c(10, 10)) + ti(age, BMI, k = c(10, 10), by = is_cnv), data = d, method = "REML")
  st <- summary(m)$s.table
  st[grep("is_cnvCNV", rownames(st))[grep("ti\\(", rownames(st)[grep("is_cnvCNV", rownames(st))])], "p-value"]
}

test_2way <- function(d) {
  d <- copy(d); d[, grp := as.ordered(bmi_grp)]
  m <- gam(y ~ grp + s(age, k = 10) + s(age, k = 10, by = grp), data = d, method = "REML")
  s <- summary(m)
  c(curve_p = s$s.table[grep("grp", rownames(s$s.table)), "p-value"],
    offset = unname(s$p.coeff[2]), offset_p = unname(s$p.pv[2]))
}

res <- list(); plots <- list()
for (oc in names(OUT)) {
  d <- latest_rec(oc)
  for (L in GROUPS) {
    dc <- split_grp(d[locus == L], L)
    p3 <- if (L == "No CNV") NA_real_ else P3[outcome == OUT[[oc]] & locus == L, p_CNVxAgexBMI]
    w <- W[outcome == OUT[[oc]] & group == L]
    p2 <- test_2way(dc)
    kids <- dc[age >= 6 & age < 13]
    diff_kids <- if (kids[, uniqueN(bmi_grp)] == 2) t.test(y ~ bmi_grp, data = kids) else NULL
    res[[length(res) + 1]] <- data.table(outcome = OUT[oc], locus = L, split = SPLIT[[L]], n_carriers = nrow(dc),
      n_high_bmi = dc[bmi_grp == levels(bmi_grp)[2], .N], n_controls = d[locus == "No CNV", .N],
      p_CNVxAgexBMI = p3, p_age_curve_differs = p2[["curve_p"]], high_minus_low_offset = p2[["offset"]], p_offset = p2[["offset_p"]],
      diff_6_13y = if (!is.null(diff_kids)) unname(diff(diff_kids$estimate)) else NA, p_6_13y = if (!is.null(diff_kids)) diff_kids$p.value else NA,
      n_6_13y = nrow(kids))
    r <- res[[length(res)]]
    lab <- paste0(sprintf("age \u00d7 BMI p = %s (GCV) / %s (REML)", pfx(w$within_GCV), pfx(w$within_REML)),
                  if (L != "No CNV") sprintf("\nvs controls p = %s (GCV) / %s (REML)", pfx(w$threeway_GCV), pfx(w$threeway_REML)) else "",
                  sprintf("\nhigh \u2212 low BMI offset %+.2f, p = %s", r$high_minus_low_offset, pfx(r$p_offset)))
    plots[[paste(oc, L)]] <- ggplot(dc, aes(age, y, colour = bmi_grp, fill = bmi_grp)) +
      # (6-13 y grey band removed, JS 2026-10-12)
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_point(size = if (L == "No CNV") 0.5 else 1.1, alpha = if (L == "No CNV") 0.2 else 0.55, shape = 16) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
      {if (oc == "HEIGHT" && nrow(UK[locus == L]) > 0) list(
        geom_vline(xintercept = 42, linetype = 3, colour = "grey50"),
        annotate("text", x = 43, y = 4.5, label = "UK Biobank adults", hjust = 0, vjust = 1, size = 3.4, colour = "grey30"),
        geom_pointrange(data = UK[locus == L], aes(x = mean_age, y = mean_hZ, ymin = mean_hZ - 1.96 * se_hZ, ymax = mean_hZ + 1.96 * se_hZ, colour = bmi_grp),
                        inherit.aes = FALSE, shape = 15, size = 0.45, linewidth = 0.7, position = position_dodge(2.5), show.legend = FALSE))} +
      annotate("text", x = if (oc == "HEIGHT") 72 else 40, y = -4.6, label = lab, hjust = 1, vjust = 0, size = 3.7) +
      scale_colour_manual(values = c("#00A5AD", "#F0645A")) + scale_fill_manual(values = c("#00A5AD", "#F0645A")) +
      guides(colour = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +   # >= 70th listed first (JS)
      coord_cartesian(xlim = if (oc == "HEIGHT") c(2, 72) else c(2, 40), ylim = c(-4.6, 4.6)) +
      labs(x = "Age (years)", y = paste(OUT[oc], "(SD)"), colour = NULL, fill = NULL,
           title = if (L == GROUPS[1]) OUT[oc] else NULL, subtitle = sprintf("%s (n = %d)", if (L == "No CNV") "Non-carrier controls (SSC + SVIP)" else L, nrow(dc))) +
      theme_classic(base_size = 15) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.position.inside = if (oc == "HEIGHT") c(0.3, 0.93) else c(0.75, 0.93), legend.background = element_blank(), legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
  }
}
res <- rbindlist(res)
fwrite(res, file.path(HERE, "cross_sectional_tests_split.csv"))
print(res, digits = 3)

# third column: head circumference by ASD, < 18 y (side_panel_hc_asd.csv from 09_side_panels.R)
asd <- fread(file.path(HERE, "side_panel_hc_asd.csv"))
pf <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))
for (L in LOCI) {
  a <- asd[locus == L]
  b <- data.table(grp = factor(c("No ASD", "ASD"), c("No ASD", "ASD")), m = c(a$mean_noASD, a$mean_ASD), se = c(a$se_noASD, a$se_ASD), n = c(a$n_noASD, a$n_ASD))
  plots[[paste("ASD", L)]] <- ggplot(b, aes(grp, m, fill = grp)) + geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(width = 0.65) + geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2) +
    geom_text(aes(y = -4.4, label = paste0("n = ", n)), size = 2.7) +
    annotate("text", x = 1.5, y = 4.3, size = 2.7, vjust = 1,
             label = sprintf("Δ %+.2f, p = %s
adjusted %+.2f, p = %s", a$diff, pf(a$p_welch), a$beta_cohort_adj, pf(a$p_cohort_adj))) +
    scale_fill_manual(values = c("#8C9A9A", "#F28E2B"), guide = "none") + coord_cartesian(ylim = c(-4.6, 4.6)) +
    labs(x = NULL, y = "Head circumference (SD)", title = if (L == LOCI[1]) "HC by ASD, ages 2-40" else NULL, subtitle = L) +
    theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
ov <- asd[grepl("Overgrowth", locus)]; un <- asd[grepl("Undergrowth", locus)]
# JS 2026-10-12: ASD column dropped from the main figure (now panel C of the supplement, fig_supp_bars)
ord <- c(rbind(paste("HEIGHT", LOCI), paste("HEADCIRCUM", LOCI)))
fig <- wrap_plots(plots[ord], ncol = 2) +
  plot_annotation(theme = theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 11), plot.caption = element_text(size = 9)),
                  title = "Height and head circumference by age and BMI in 16p11.2 and 22q11.2 CNV carriers",
                  subtitle = "BMI ≥ 70th vs ≤ 30th centile (middle 40% dropped), all loci.
One record per person (latest visit), G2MH + Cardiff + SVIP + SSC + Searchlight (height only). Curves: GAM \u00b1 95% CI.",
                  caption = paste("CNV \u00d7 age \u00d7 BMI: tensor-product GAM with continuous BMI vs SSC/SVIP/G2MH non-carrier controls (Anjali's test).
",
                                  "Offset: high- minus low-BMI shift of the carrier age curve. Z: < 18 y childsds WHO/UK90; \u2265 18 y within-sample Box-Cox."))
ggsave(file.path(HERE, "fig_cross_sectional_split.png"), fig, width = 12, height = 14, dpi = 200)
ggsave(file.path(HERE, "fig_cross_sectional_split.pdf"), fig, width = 12, height = 14, device = cairo_pdf)

# ---------------------------------------------------------------- main figure (JS 2026-10-12)
# Trajectory panels (rows = loci; left = height, right = head circumference) plus a bottom row with the pooled
# ASD bars (overgrowth / undergrowth, locus-standardised) from supp_bars_results.csv (10_supplement.R).
# Panels lettered down each column (A-E left, F-J right). No title/caption: legend in fig_main_growth_legend.txt.
sb <- fread(file.path(HERE, "supp_bars_results.csv"))
pf2 <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))
asd_bar <- function(panel_prefix, ylab) {
  r <- sb[startsWith(panel, panel_prefix) & grepl("growth|Non-carriers", group)]
  r <- r[order(group != "Non-carriers")]   # non-carriers first (JS)
  r[group == "Non-carriers", group := "Non-carriers\n(SSC sib-pairs)"]
  r[grepl("growth", group), group := sub("\n.*", "\nCNVs", group)][, group := factor(group, unique(group))]
  d <- rbind(r[, .(group, grp = "No ASD", m = mean_0, se = se_0, n = n_0)], r[, .(group, grp = "ASD", m = mean_1, se = se_1, n = n_1)])
  d[, grp := factor(grp, c("No ASD", "ASD"))]
  r[, lab := sprintf("Δ %+.2f, p = %s\nadj. %+.2f\np = %s", diff, pf2(p_welch), beta_adj, pf2(p_adj))][, bold := p_adj < 0.05]
  r[within_family %in% TRUE, `:=`(lab = sprintf("sib-pair\nΔ %+.2f\np = %s", diff, pf2(p_welch)), bold = p_welch < 0.05)]
  yr <- range(c(d$m - d$se, d$m + d$se)); top <- yr[2] + 0.75 * diff(yr); bot <- yr[1] - 0.18 * diff(yr)
  ggplot(d, aes(group, m, fill = grp)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(y = bot, label = n), position = position_dodge(0.8), size = 3.3) +
    geom_text(data = r, aes(group, top, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 3.4, vjust = 1, lineheight = 0.9) +
    scale_fill_manual(values = c("#8C9A9A", "#F28E2B")) + coord_cartesian(ylim = c(bot, top)) +
    labs(x = NULL, y = ylab, fill = NULL, subtitle = "ASD vs no ASD") +
    theme_classic(base_size = 15) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
                                           legend.text = element_text(size = 9), legend.key.size = unit(10, "pt"))   # no grey panel background (JS)
}
plots[["ASD HEIGHT"]] <- asd_bar("B.", "Height Z")
plots[["ASD HEADCIRCUM"]] <- asd_bar("C.", "Head circumference Z")
ordm <- c(paste("HEIGHT", GROUPS), "ASD HEIGHT", paste("HEADCIRCUM", GROUPS), "ASD HEADCIRCUM")   # column-wise -> tags run down each column
figm <- wrap_plots(plots[ordm], ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(HERE, "fig_main_growth.png"), figm, width = 12, height = 20.5, dpi = 200)
ggsave(file.path(HERE, "fig_main_growth.pdf"), figm, width = 12, height = 20.5, device = cairo_pdf)
