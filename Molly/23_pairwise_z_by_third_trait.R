# Experiment (MS 2026-10-15): pairwise Z differences by age (as in 22_pairwise_z_differences.R), with each panel split by the
# sign of the THIRD trait: height - BMI by HC Z, HC - BMI by height Z, height - HC by BMI Z (> 0 vs <= 0).
# People: one record per person, latest visit before 18 with height, BMI and HC (ages 2-17).
# Lines: GAM s(age, k = 6) per group. Stats (REML): y ~ grp + s(age) + s(age, by = grp), grp ordered (third Z <= 0 reference):
#   offset = overall shift of the third-Z > 0 line, shape p = difference smooth.
# Output: fig_pairwise_z_by_third_trait.png, pairwise_z_by_third_trait_tests.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
AGE_MAX <- 18
GROUPS <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
GLAB <- c("No CNV" = "Non-carrier controls", "16p11.2 del" = "16p11.2 BP4-BP5 DEL", "16p11.2 dup" = "16p11.2 BP4-BP5 DUP",
          "22q11.2 del" = "22q11.2 A-D DEL", "22q11.2 dup" = "22q11.2 A-D DUP")
COMP <- c(ht_bmi = "Height \u2212 BMI (SD)", hc_bmi = "HC \u2212 BMI (SD)", ht_hc = "Height \u2212 HC (SD)")
THIRD <- c(ht_bmi = "HC", hc_bmi = "Height", ht_hc = "BMI")   # trait not in the difference
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))

d <- tab[!is.na(HEIGHT) & !is.na(BMI) & !is.na(HEADCIRCUM) & !is.na(age) & age >= 2 & age < AGE_MAX][order(id, -age)][, .SD[1], by = id]
long <- rbind(d[, .(id, locus, age, comp = "ht_bmi", y = HEIGHT - BMI, z3 = HEADCIRCUM)],
              d[, .(id, locus, age, comp = "hc_bmi", y = HEADCIRCUM - BMI, z3 = HEIGHT)],
              d[, .(id, locus, age, comp = "ht_hc", y = HEIGHT - HEADCIRCUM, z3 = BMI)])
long[, grp := factor(ifelse(z3 > 0, "third trait Z > 0", "third trait Z \u2264 0"), levels = c("third trait Z \u2264 0", "third trait Z > 0"))]
long[, comp_lab := factor(sprintf("%s\nsplit by %s Z", COMP[comp], THIRD[comp]), levels = sprintf("%s\nsplit by %s Z", COMP, THIRD))]
long[, locus_lab := factor(GLAB[as.character(locus)], levels = GLAB)]

res <- long[, {
  x <- copy(.SD)[, g := factor(grp, levels = levels(grp), ordered = TRUE)]
  m <- gam(y ~ g + s(age, k = 6) + s(age, k = 6, by = g), data = x, method = "REML"); s <- summary(m)
  .(n_low = sum(grp == levels(grp)[1]), n_high = sum(grp == levels(grp)[2]),
    offset = unname(s$p.coeff[2]), offset_p = unname(s$p.pv[2]), shape_p = s$s.table[grep("^s\\(age\\):g", rownames(s$s.table)), "p-value"])
}, by = .(comp, comp_lab, locus, locus_lab)]
fwrite(res[, !"comp_lab"], file.path(OUT, "pairwise_z_by_third_trait_tests.csv"))
print(res[, !c("comp_lab", "locus_lab")], digits = 3)
res[, lab := gsub("p = < ", "p < ", sprintf("Z>0 \u2212 Z\u22640 offset %+.2f, p = %s\nshape p = %s", offset, pfx(offset_p), pfx(shape_p)), fixed = TRUE)]
res[, nlab := sprintf("n = %d | %d", n_low, n_high)]

cols <- c("#00A5AD", "#F0645A")
fig <- ggplot(long, aes(age, y, colour = grp, fill = grp)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(aes(size = locus == "No CNV", alpha = locus == "No CNV"), shape = 16) +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 0.5), guide = "none") + scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 0.2), guide = "none") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
  geom_text(data = res, aes(x = AGE_MAX, y = -6.4, label = lab), hjust = 1, vjust = 0, size = 3.1, lineheight = 0.9, inherit.aes = FALSE) +
  geom_text(data = res, aes(x = 2.5, y = 6.4, label = nlab), hjust = 0, vjust = 1, size = 3.1, inherit.aes = FALSE) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_grid(comp_lab ~ locus_lab, switch = "y") +
  coord_cartesian(xlim = c(2, AGE_MAX), ylim = c(-6.4, 6.4)) +
  labs(x = "Age (years)", y = NULL, colour = NULL, fill = NULL,
       title = "Pairwise differences in growth Z-scores by age, split by the third trait",
       subtitle = sprintf("One record per person (latest visit with height, BMI and HC), ages 2-%d. Row label says which trait sets the colour. n = third Z \u2264 0 | > 0.", AGE_MAX - 1)) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(10, "pt"), plot.title = element_text(face = "bold"), legend.position = "top", legend.justification = "left")
ggsave(file.path(OUT, "fig_pairwise_z_by_third_trait.png"), fig, width = 18, height = 11, dpi = 200)
