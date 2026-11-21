# 81_fig_longitudinal_v2.R  (2026-11-20, BMI modifier paper v2, work package A)
#
# What it does: draws the two longitudinal figure blocks from the models of 80_longitudinal_v2.R (no new statistics).
#   block_hc_longitudinal     (20 x 5.5 in; tags g, h, i): g design schematic; h children by genotype and cohort with visit pairs
#                             and median age; i change in head-circumference Z to the next visit against BMI Z at the visit
#                             (covariate-adjusted partial residuals, mean and 95% CI in fifths of BMI, fitted slope and P).
#   block_height_longitudinal (20 x 5.5 in; tags g, h): g change in height Z to the next visit per SD of BMI Z by age window;
#                             h children by genotype, cohort and window.
# Inputs: DATA/longitudinal_v2_panels.rds (participant-level panels written by 80), run_v2_2026-11-20/longitudinal/results_*.csv.
# Outputs (run_v2_2026-11-20/longitudinal/): block_hc_longitudinal.{png,pdf,csv}, block_height_longitudinal.{png,pdf,csv}.
# null model: none; the P values shown are those of 80_longitudinal_v2.R (Wald tests of the BMI coefficients in the mixed models;
#   likelihood-ratio tests for BMI x genotype heterogeneity and BMI x window).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
OUT <- file.path(HERE, "run_v2_2026-11-20", "longitudinal")
P <- readRDS(file.path(DATA, "longitudinal_v2_panels.rds")); IC <- P$IC; IH <- P$IH; PR <- P$PRC; RC <- P$RC; RH <- P$RH
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"); ALL <- "All four CNVs"; G5 <- c(LOCI, ALL)
CARR <- "#4B6EA9"; INK <- "grey15"
COHC <- c("Simons VIP" = "#4B6EA9", "Searchlight" = "#5BB381", "Cardiff" = "#E0A458", "CHOP" = "#7B5EA7")
cohlab <- function(x) factor(ifelse(x == "SVIP", "Simons VIP", x), names(COHC))
fp <- function(p) sapply(p, function(v) if (is.na(v)) "NA" else if (v < 0.001) formatC(v, format = "e", digits = 0) else if (v < 0.01) sprintf("%.4f", v) else if (v < 0.1) sprintf("%.3f", v) else sprintf("%.2f", v))
face <- function(p) ifelse(!is.na(p) & p < 0.05, "bold", "plain")
TH <- theme_classic(base_size = 15) + theme(plot.tag = element_text(face = "bold", size = 22), plot.tag.position = c(0, 1), plot.title = element_text(size = 15, face = "bold"),
                                             plot.subtitle = element_text(size = 13), strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold"),
                                             legend.position = "top", legend.title = element_blank(), plot.margin = margin(8, 10, 6, 14))

# ================================================================ head circumference block
# g: design schematic
bx <- function(x0, x1, y0, y1, fill, col) annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = col, linewidth = 0.7)
pg <- ggplot() + xlim(0, 100) + ylim(0, 100) +
  annotate("segment", x = 3, xend = 98, y = 47, yend = 47, colour = "grey45", linewidth = 0.9, arrow = arrow(length = unit(0.13, "in"), type = "closed")) +
  annotate("text", x = 98, y = 41, label = "age", hjust = 1, size = 4.2, colour = "grey35") +
  annotate("point", x = c(25, 75), y = 47, size = 4.5, colour = INK) +
  bx(1, 48, 58, 96, "#e9eef7", CARR) +
  annotate("text", x = 24.5, y = 89, label = "Visit 1", fontface = "bold", size = 4.8, colour = INK) +
  annotate("text", x = 24.5, y = 72, label = "BMI Z\nhead\ncircumference Z", size = 4.1, colour = INK, lineheight = 0.92) +
  bx(52, 99, 58, 96, "#fdf0e3", "#b8610f") +
  annotate("text", x = 75.5, y = 89, label = "Visit 2 (next)", fontface = "bold", size = 4.8, colour = INK) +
  annotate("text", x = 75.5, y = 72, label = "head\ncircumference Z", size = 4.1, colour = INK, lineheight = 0.92) +
  annotate("text", x = 50, y = 28, label = "Outcome = change in\nhead-circumference Z\n(visit 2 minus visit 1)", fontface = "bold", size = 4.3, colour = INK, lineheight = 0.92) +
  annotate("text", x = 50, y = 6, label = "Predictor: BMI Z at visit 1\nAdjusted for head Z at visit 1, interval,\nage, sex, cohort; child random intercept", size = 3.5, colour = "grey30", lineheight = 0.92) +
  labs(tag = "g", title = "Next-visit design") + theme_void(base_size = 15) +
  theme(plot.tag = element_text(face = "bold", size = 22), plot.tag.position = c(0, 1), plot.title = element_text(size = 15, face = "bold", hjust = 0.5), plot.margin = margin(8, 6, 6, 14))

# h: children by genotype and cohort
sh <- IC[, .(children = uniqueN(id), pairs = .N), by = .(genotype = factor(locus, LOCI), cohort = cohlab(cohort))]
st <- IC[, .(children = uniqueN(id), pairs = .N, median_age = median(age0)), by = .(genotype = factor(locus, LOCI))]
st[, lab := sprintf("n = %d\n%s pairs\nage %.1f y", children, format(pairs, big.mark = ",", trim = TRUE), median_age)]
ph <- ggplot(sh, aes(genotype, children, fill = cohort)) + geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  geom_text(data = st, aes(genotype, children, label = lab), inherit.aes = FALSE, vjust = -0.15, size = 4, lineheight = 0.92, colour = INK) +
  scale_fill_manual(values = COHC, drop = TRUE) + scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
  labs(tag = "h", x = NULL, y = "Children with repeated visits",
       title = sprintf("%d children, %s visit pairs", uniqueN(IC$id), format(nrow(IC), big.mark = ",")), subtitle = "n = children; age = median age at first visit of a pair") + TH +
  theme(axis.text.x = element_text(size = 13, colour = INK, lineheight = 0.85), legend.position = "inside", legend.position.inside = c(0.2, 0.84), legend.background = element_blank())

# i: change in head-circumference Z against BMI Z, partial residuals in fifths of BMI
PR[, panel := factor(panel, G5)]
PR[, bin := cut(B0, unique(quantile(B0, 0:5 / 5)), include.lowest = TRUE, labels = FALSE), by = panel]
bn <- PR[, .(bmi = mean(B0), y = mean(pr), se = sd(pr) / sqrt(.N), n_pairs = .N), by = .(panel, bin)][order(panel, bin)]
bn[, `:=`(lo = y - 1.96 * se, hi = y + 1.96 * se)]
ft <- merge(unique(PR[, .(panel, slope, se, p)]), PR[, .(children = uniqueN(id), pairs = .N, x0 = quantile(B0, 0.02), x1 = quantile(B0, 0.98)), by = panel], by = "panel")
ft[, lab := sprintf("slope = %+.3f\np = %s\nn = %d (%s pairs)", slope, fp(p), children, format(pairs, big.mark = ",", trim = TRUE))]
ln <- ft[, .(x = c(x0, x1), y = slope * c(x0, x1)), by = panel]
hetp <- RC[genotype == ALL & scope == "All ages", p_heterogeneity]
ytop <- max(bn$hi) * 1.9; ybot <- min(bn$lo) * 1.1
pi_ <- ggplot(bn, aes(bmi, y)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_line(data = ln, aes(x, y), colour = CARR, linewidth = 1.1) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = CARR, linewidth = 0.7) + geom_point(size = 2.8, colour = CARR) +
  geom_text(data = ft, aes(x = -Inf, y = Inf, label = lab, fontface = face(p)), hjust = -0.06, vjust = 1.15, size = 4.3, lineheight = 0.92, colour = INK) +
  facet_wrap(~panel, nrow = 1) + coord_cartesian(ylim = c(ybot, ytop)) +
  labs(tag = "i", x = "BMI Z at the visit", y = "Change in head-circumference Z\nto the next visit (adjusted)",
       title = "BMI at a visit and head growth to the next visit",
       subtitle = sprintf("Means and 95%% CI in fifths of BMI; line = fitted slope per SD of BMI.  Heterogeneity across genotypes p = %s", fp(hetp))) + TH + theme(panel.spacing = unit(1, "lines"))
blk <- pg + ph + pi_ + plot_layout(widths = c(4.1, 4.4, 11.5))
ggsave(file.path(OUT, "block_hc_longitudinal.png"), blk, width = 20, height = 5.5, dpi = 200, bg = "white")
ggsave(file.path(OUT, "block_hc_longitudinal.pdf"), blk, width = 20, height = 5.5, bg = "white")
fwrite(rbind(cbind(panel_tag = "h", sh[order(genotype, cohort)]), cbind(panel_tag = "h_total", st[, .(genotype, children, pairs, median_age)]),
             cbind(panel_tag = "i_bins", bn[, .(genotype = panel, bmi_fifth = bin, mean_bmi_z = bmi, mean_change_adjusted = y, se, lo, hi, pairs = n_pairs)]),
             cbind(panel_tag = "i_fit", ft[, .(genotype = panel, slope, se, p, children, pairs)]),
             data.table(panel_tag = "i_heterogeneity", genotype = "BMI x genotype (LRT)", p = hetp), fill = TRUE), file.path(OUT, "block_hc_longitudinal.csv"))

# ================================================================ height block
WIN_H <- c("2-8 y", "8-13 y", "13-18 y")
eh <- RH[scope != "All ages"][, `:=`(genotype = factor(genotype, G5), scope = factor(scope, WIN_H))]
xw <- merge(RH[scope == "All ages", .(genotype, p_bmi_x_window, children_all = children, pairs_all = pairs)], eh[, .(top = max(hi, na.rm = TRUE)), by = genotype], by = "genotype")
xw[, `:=`(genotype = factor(genotype, G5), lab = sprintf("BMI x age window\np = %s\nn = %d (%s pairs)", fp(p_bmi_x_window), children_all, format(pairs_all, big.mark = ",", trim = TRUE)))]
eh[, plab := ifelse(is.na(beta), sprintf("n = %d\n(too few)", children), sprintf("p = %s\nn = %d", fp(p), children))]
yl <- c(min(eh$lo, na.rm = TRUE) - 0.02, max(eh$hi, na.rm = TRUE) + 0.38)
pg2 <- ggplot(eh, aes(scope, beta)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = CARR, linewidth = 0.8, na.rm = TRUE) + geom_point(size = 3.4, colour = CARR, na.rm = TRUE) +
  geom_text(aes(y = ifelse(is.na(hi), 0.02, hi), label = plab, fontface = face(p)), vjust = -0.25, size = 3.6, lineheight = 0.9, colour = INK) + scale_x_discrete(expand = expansion(add = 0.75)) +
  geom_text(data = xw, aes(x = 3.7, y = Inf, label = lab, fontface = face(p_bmi_x_window)), hjust = 1, vjust = 1.1, size = 4.1, lineheight = 0.92, colour = INK) +
  facet_wrap(~genotype, nrow = 1) + coord_cartesian(ylim = yl) +
  labs(tag = "g", x = "Age at the visit", y = "Change in height Z to the next visit\nper SD of BMI Z (95% CI)", title = "BMI at a visit and height gain to the next visit, by age",
       subtitle = "One model per genotype (age window x BMI); all four CNVs adjusted for genotype; n = children") + TH + theme(panel.spacing = unit(1, "lines"), axis.text.x = element_text(size = 12.5, colour = INK))
shh <- IH[, .(children = uniqueN(id), pairs = .N), by = .(genotype = factor(locus, LOCI), window = factor(win, WIN_H), cohort = cohlab(cohort))]
sth <- IH[, .(children = uniqueN(id), pairs = .N), by = .(genotype = factor(locus, LOCI), window = factor(win, WIN_H))]
ph2 <- ggplot(shh, aes(genotype, children, fill = cohort)) + geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(data = sth, aes(genotype, children, label = sprintf("%d\n(%s)", children, format(pairs, big.mark = ",", trim = TRUE))), inherit.aes = FALSE, vjust = -0.15, size = 3.7, lineheight = 0.9, colour = INK) +
  facet_wrap(~window, nrow = 1) + scale_fill_manual(values = COHC, drop = TRUE) + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_x_discrete(labels = function(x) sub("11.2 ", "\n", x)) +
  labs(tag = "h", x = NULL, y = "Children (visit pairs)", title = sprintf("%d children, %s visit pairs", uniqueN(IH$id), format(nrow(IH), big.mark = ","))) + TH +
  theme(axis.text.x = element_text(size = 11.5, colour = INK, lineheight = 0.85), panel.spacing = unit(0.8, "lines"))
blk2 <- pg2 + ph2 + plot_layout(widths = c(13.6, 6.4))
ggsave(file.path(OUT, "block_height_longitudinal.png"), blk2, width = 20, height = 5.5, dpi = 200, bg = "white")
ggsave(file.path(OUT, "block_height_longitudinal.pdf"), blk2, width = 20, height = 5.5, bg = "white")
fwrite(rbind(cbind(panel_tag = "g", eh[, .(genotype, window = scope, cohorts, children, pairs, median_age0, beta, se, lo, hi, p)]),
             cbind(panel_tag = "g_all_ages", RH[scope == "All ages", .(genotype, window = scope, cohorts, children, pairs, median_age0, beta, se, lo, hi, p, p_bmi_x_window, p_heterogeneity)]),
             cbind(panel_tag = "h", shh[order(genotype, window, cohort)]), fill = TRUE), file.path(OUT, "block_height_longitudinal.csv"))
cat("figure blocks written to", OUT, "\n"); print(ft[, .(panel, slope, p, children, pairs)]); print(eh[, .(genotype, scope, beta, p, children)])
