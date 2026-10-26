# ASHG poster: standalone height panels (non-carrier controls, 16p11.2 del, 22q11.2 del) -- main-figure panels A, B and D as single files (MS 2026-10-14).
# Same data, split, curves, stats and styling as the height column of fig_main_growth (02b_cross_sectional_split.R):
#   one record per person (latest visit with height and BMI), ages 2-40; BMI >= 70th vs <= 30th centile (middle dropped);
#   GAM s(age, k = 6) per BMI group. UK Biobank adult anchors left out (MS 2026-10-14).
# Stats: age x BMI and vs-controls p from within_vs_controls_tests.csv (14_within_vs_controls.R, continuous BMI);
#   high - low BMI offset from the ordered-factor GAM y ~ grp + s(age) + s(age, by = grp), REML (02b test_2way).
# Output: ASHGposter_height_{nocnv,16pdel,22qdel}.png

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
W <- fread(file.path(HERE, "within_vs_controls_tests.csv"))[outcome == "Height"]
lv <- c("BMI \u2264 30th centile", "BMI \u2265 70th centile")
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))

d <- tab[!is.na(HEIGHT) & !is.na(BMI) & !is.na(age) & age >= 2 & age <= 40][order(id, -age)][, .SD[1], by = id]
d[, y := HEIGHT]

test_2way <- function(dc) {
  dc <- copy(dc); dc[, grp := as.ordered(bmi_grp)]
  m <- gam(y ~ grp + s(age, k = 10) + s(age, k = 10, by = grp), data = dc, method = "REML")
  s <- summary(m)
  c(curve_p = s$s.table[grep("grp", rownames(s$s.table)), "p-value"], offset = unname(s$p.coeff[2]), offset_p = unname(s$p.pv[2]))
}

PANELS <- list(c("No CNV", "nocnv"), c("16p11.2 del", "16pdel"), c("22q11.2 del", "22qdel"))
POSTER_NAME <- c("16p11.2 del" = "16p11.2 BP4-BP5 DEL", "22q11.2 del" = "22q11.2 A-D DEL")   # poster labels
for (P in PANELS) {
  L <- P[1]
  dc <- d[locus == L][, bmi_grp := factor(ifelse(BMI >= qnorm(.7), lv[2], ifelse(BMI <= qnorm(.3), lv[1], NA)), levels = lv)][!is.na(bmi_grp)]
  w <- W[group == L]; p2 <- test_2way(dc)
  cat(L, ": n =", nrow(dc), " offset", round(p2[["offset"]], 2), "p", signif(p2[["offset_p"]], 2), " curve p", signif(p2[["curve_p"]], 2), "\n")
  lab <- paste0(sprintf("age \u00d7 BMI p = %s (GCV) / %s (REML)", pfx(w$within_GCV), pfx(w$within_REML)),
                if (L != "No CNV") sprintf("\nvs controls p = %s (GCV) / %s (REML)", pfx(w$threeway_GCV), pfx(w$threeway_REML)) else "",
                sprintf("\nhigh \u2212 low BMI offset %+.2f, p = %s", p2[["offset"]], pfx(p2[["offset_p"]])))
  plot <- ggplot(dc, aes(age, y, colour = bmi_grp, fill = bmi_grp)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = if (L == "No CNV") 0.5 else 1.1, alpha = if (L == "No CNV") 0.2 else 0.55, shape = 16) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = 40, y = -4.6, label = lab, hjust = 1, vjust = 0, size = 3.7) +
    scale_colour_manual(values = c("#00A5AD", "#F0645A")) + scale_fill_manual(values = c("#00A5AD", "#F0645A")) +
    guides(colour = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    coord_cartesian(xlim = c(2, 40), ylim = c(-4.6, 4.6)) +
    labs(x = "Age (years)", y = "Height (SD)", colour = NULL, fill = NULL,
         subtitle = sprintf("%s (n = %d)", if (L == "No CNV") "Non-carrier controls (SSC + SVIP)" else POSTER_NAME[[L]], nrow(dc))) +
    theme_classic(base_size = 15) + theme(legend.text = element_text(size = 3.7 * .pt), legend.position = "inside",   # legend text = p-value text size
                                          legend.position.inside = c((40 - 0.1) / 41.8, 1), legend.justification = c(1, 1), legend.margin = margin(0, 0, 0, 0),   # right edge at x = 40 (stats text), top of panel
                                          legend.background = element_blank(), legend.key.height = unit(11, "pt"), plot.subtitle = element_text(hjust = 0.5),
                                          plot.margin = margin(1, 5.5, 1, 5.5))   # trim top/bottom whitespace (default 5.5 pt)
  ggsave(file.path(OUT, sprintf("ASHGposter_height_%s.png", P[2])), plot, width = 7.5, height = 3 - 9 / 72, dpi = 200)   # was 6 x 3.6; height reduced by the 9 pt of margin removed so the panel itself is unchanged
}
