# Cross-sectional trajectories: height and head circumference Z by age, BMI >= 90th centile vs rest, per CNV.
# One record per person (latest visit with the outcome and BMI), ages 2-40 y.
# Tests per CNV x outcome:
#   3-way  CNV x age x BMI   -- Anjali's tensor-product GAM (continuous BMI, SSC/SVIP/G2MH non-carrier controls):
#          Z ~ is_cnv + s(age) + s(BMI) + s(age, by = is_cnv) + ti(age, BMI) + ti(age, BMI, by = is_cnv)
#   2-way  age x BMI stratum within carriers -- ordered-factor difference smooth (does the age curve differ)
#          plus the parametric stratum offset (is the high-BMI curve shifted overall)
# Output: fig_cross_sectional.{png,pdf}, cross_sectional_tests.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
Z90 <- qnorm(0.9)
OUT <- c(HEIGHT = "Height", HEADCIRCUM = "Head circumference")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")

latest_rec <- function(outcome) {
  d <- tab[!is.na(get(outcome)) & !is.na(BMI) & !is.na(age) & age >= 2 & age <= 40]
  d <- d[order(id, -age)][, .SD[1], by = id]
  d[, y := get(outcome)]
  d[, bmi_grp := factor(ifelse(BMI > Z90, "BMI > 90th centile", "BMI \u2264 90th centile"),
                        levels = c("BMI \u2264 90th centile", "BMI > 90th centile"))]
  d
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
  for (L in LOCI) {
    dc <- d[locus == L]
    p3 <- test_3way(d[locus %in% c(L, "No CNV")])
    p2 <- test_2way(dc)
    kids <- dc[age >= 6 & age < 13]
    diff_kids <- if (kids[, uniqueN(bmi_grp)] == 2) t.test(y ~ bmi_grp, data = kids) else NULL
    res[[length(res) + 1]] <- data.table(outcome = OUT[oc], locus = L, n_carriers = nrow(dc),
      n_high_bmi = dc[bmi_grp == levels(bmi_grp)[2], .N], n_controls = d[locus == "No CNV", .N],
      p_CNVxAgexBMI = p3, p_age_curve_differs = p2[["curve_p"]], high_minus_low_offset = p2[["offset"]], p_offset = p2[["offset_p"]],
      diff_6_13y = if (!is.null(diff_kids)) unname(diff(diff_kids$estimate)) else NA, p_6_13y = if (!is.null(diff_kids)) diff_kids$p.value else NA,
      n_6_13y = nrow(kids))
    r <- res[[length(res)]]
    lab <- sprintf("CNV \u00d7 age \u00d7 BMI p = %s\nhigh \u2212 low BMI offset %+.2f, p = %s",
                   format.pval(r$p_CNVxAgexBMI, digits = 2, eps = 1e-4), r$high_minus_low_offset, format.pval(r$p_offset, digits = 2, eps = 1e-4))
    plots[[paste(oc, L)]] <- ggplot(dc, aes(age, y, colour = bmi_grp, fill = bmi_grp)) +
      annotate("rect", xmin = 6, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey92") +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_point(size = 1.1, alpha = 0.55, shape = 16) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
      annotate("text", x = 40, y = -4.6, label = lab, hjust = 1, vjust = 0, size = 2.7) +
      scale_colour_manual(values = c("#00A5AD", "#F0645A"), drop = FALSE) +
      scale_fill_manual(values = c("#00A5AD", "#F0645A"), drop = FALSE) +
      coord_cartesian(xlim = c(2, 40), ylim = c(-4.6, 4.6)) +
      labs(x = "Age (years)", y = paste(OUT[oc], "(SD)"), colour = NULL, fill = NULL,
           title = if (L == LOCI[1]) OUT[oc] else NULL, subtitle = sprintf("%s (n = %d)", L, nrow(dc))) +
      theme_classic(base_size = 10) + theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))
  }
}
res <- rbindlist(res)
fwrite(res, file.path(HERE, "cross_sectional_tests.csv"))
print(res, digits = 3)

ord <- c(rbind(paste("HEIGHT", LOCI), paste("HEADCIRCUM", LOCI)))
fig <- wrap_plots(plots[ord], ncol = 2, guides = "collect") +
  plot_annotation(title = "Height and head circumference by age and BMI in 16p11.2 and 22q11.2 CNV carriers",
                  subtitle = "One record per person (latest visit), G2MH + Cardiff + SVIP + SSC + Searchlight (height only). Grey band = 6-13 y. Curves: GAM \u00b1 95% CI.",
                  caption = paste("CNV \u00d7 age \u00d7 BMI: tensor-product GAM with continuous BMI vs SSC/SVIP/G2MH non-carrier controls (Anjali's test).",
                                  "Offset: high- minus low-BMI shift of the carrier age curve. Z: < 18 y childsds WHO/UK90; \u2265 18 y within-sample Box-Cox.")) &
  theme(legend.position = "top")
ggsave(file.path(HERE, "fig_cross_sectional.png"), fig, width = 11, height = 13, dpi = 200)
ggsave(file.path(HERE, "fig_cross_sectional.pdf"), fig, width = 11, height = 13)
