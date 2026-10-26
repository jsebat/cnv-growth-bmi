# Pairwise differences in growth Z-scores by age, per CNV (MS 2026-10-15).
# Question: does head circumference relative to BMI (and height) change across development, and differently by CNV?
# People: one record per person -- latest visit before age 18 with height, BMI AND head circumference all measured (ages 2-17; plots and stats).
# Rows (y = difference in Z): height - BMI, HC - BMI, height - HC. Columns: non-carriers, then the four CNVs.
# Trendline: GAM s(age, k = 6) +/- 95% CI. CNV panels also show the non-carrier curve (dashed grey) for reference.
# Stats per panel (REML):
#   age trend p   -- y ~ s(age): does the difference change with age within the group
#   vs controls   -- y ~ grp + s(age) + s(age, by = grp), grp ordered (non-carrier reference):
#                    shape p = difference smooth (does the age curve differ from controls'), offset = overall shift vs controls
# Output: fig_pairwise_z_differences.png, pairwise_z_differences_tests.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
AGE_MAX <- 18   # MS 2026-10-15: under 18 only (Z reference switches to within-sample at 18); visits >= 18 dropped before taking latest
GROUPS <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
GLAB <- c("No CNV" = "Non-carrier controls", "16p11.2 del" = "16p11.2 BP4-BP5 DEL", "16p11.2 dup" = "16p11.2 BP4-BP5 DUP",
          "22q11.2 del" = "22q11.2 A-D DEL", "22q11.2 dup" = "22q11.2 A-D DUP")
COMP <- c(ht_bmi = "Height \u2212 BMI (SD)", hc_bmi = "HC \u2212 BMI (SD)", ht_hc = "Height \u2212 HC (SD)")
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))

d <- tab[!is.na(HEIGHT) & !is.na(BMI) & !is.na(HEADCIRCUM) & !is.na(age) & age >= 2 & age < AGE_MAX][order(id, -age)][, .SD[1], by = id]
d[, `:=`(ht_bmi = HEIGHT - BMI, hc_bmi = HEADCIRCUM - BMI, ht_hc = HEIGHT - HEADCIRCUM)]
long <- melt(d[, .(id, cohort, locus, age, ht_bmi, hc_bmi, ht_hc)], id.vars = c("id", "cohort", "locus", "age"), variable.name = "comp", value.name = "y")
long[, comp := factor(COMP[as.character(comp)], levels = COMP)][, locus_lab := factor(GLAB[as.character(locus)], levels = GLAB)]

# tests
res <- list()
for (cc in levels(long$comp)) for (L in GROUPS) {
  x <- long[comp == cc & locus == L]
  m0 <- gam(y ~ s(age, k = 6), data = x, method = "REML")
  r <- data.table(comp = cc, locus = L, n = nrow(x), mean_diff = mean(x$y), age_trend_p = summary(m0)$s.table[1, "p-value"],
                  shape_vs_ctrl_p = NA_real_, offset_vs_ctrl = NA_real_, offset_vs_ctrl_p = NA_real_)
  if (L != "No CNV") {
    xx <- copy(long[comp == cc & locus %in% c("No CNV", L)])[, grp := factor(ifelse(locus == "No CNV", "ctrl", "cnv"), levels = c("ctrl", "cnv"), ordered = TRUE)]
    m <- gam(y ~ grp + s(age, k = 6) + s(age, k = 6, by = grp), data = xx, method = "REML"); s <- summary(m)
    r[, `:=`(shape_vs_ctrl_p = s$s.table[grep("grp", rownames(s$s.table)), "p-value"], offset_vs_ctrl = unname(s$p.coeff[2]), offset_vs_ctrl_p = unname(s$p.pv[2]))]
  }
  res[[length(res) + 1]] <- r
}
res <- rbindlist(res)
fwrite(res, file.path(OUT, "pairwise_z_differences_tests.csv"))
print(res, digits = 3)

res[, comp := factor(comp, levels = COMP)][, locus_lab := factor(GLAB[locus], levels = GLAB)]
res[, lab := paste0(sprintf("age trend p = %s", pfx(age_trend_p)),
                    ifelse(locus == "No CNV", "", sprintf("\nvs controls: shape p = %s\noffset %+.2f, p = %s", pfx(shape_vs_ctrl_p), offset_vs_ctrl, pfx(offset_vs_ctrl_p))))]
res[, lab := gsub("p = < ", "p < ", lab, fixed = TRUE)]
res[, nlab := sprintf("n = %d", n)]
ctrl <- long[locus == "No CNV", .(id, cohort, age, y, comp)]   # reference curve, drawn in every column
ctrl_all <- rbindlist(lapply(levels(long$locus_lab)[-1], function(g) copy(ctrl)[, locus_lab := factor(g, levels = GLAB)]))

fig <- ggplot(long, aes(age, y)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(aes(size = locus == "No CNV", alpha = locus == "No CNV"), colour = "#4C6A92", shape = 16, show.legend = FALSE) +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 0.5)) + scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 0.15)) +
  geom_smooth(data = ctrl_all, method = "gam", formula = y ~ s(x, k = 6), se = FALSE, colour = "grey45", linetype = "22", linewidth = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), colour = "#1F3A5F", fill = "#1F3A5F", linewidth = 0.9, alpha = 0.2) +
  geom_text(data = res, aes(x = AGE_MAX, y = -6.4, label = lab), hjust = 1, vjust = 0, size = 3.1, lineheight = 0.9, inherit.aes = FALSE) +
  geom_text(data = res, aes(x = 2.5, y = 6.4, label = nlab), hjust = 0, vjust = 1, size = 3.1, inherit.aes = FALSE) +
  facet_grid(comp ~ locus_lab, switch = "y") +
  coord_cartesian(xlim = c(2, AGE_MAX), ylim = c(-6.4, 6.4)) +
  labs(x = "Age (years)", y = NULL,
       title = "Pairwise differences in growth Z-scores by age",
       subtitle = paste0("One record per person (latest visit with height, BMI and HC), ages 2-", AGE_MAX - 1, ". Solid: GAM \u00b1 95% CI; dashed grey: non-carrier curve.")) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(10, "pt"), plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_pairwise_z_differences.png"), fig, width = 18, height = 10.5, dpi = 200)
