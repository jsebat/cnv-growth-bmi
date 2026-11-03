# Supplementary figure (rows = CNV; columns = height, head circumference, BMI): CNV main effects within the Fig. 6 age bands (2-9, 10-18, 19-39,
# 40-73 y), drawn like the Fig. 6 prime panels. Data of 17_main_with_ukbb.R: children on the external references (WHO height
# and BMI, UK90 head circumference), adults relative to the non-carrier adults of their own cohort (Box-Cox scaling);
# UK Biobank adults in the height and BMI 40-73 band. Bars: mean Z (95% CI) in non-carrier controls (grey) and carriers
# (blue); printed value: P of the age- and sex-adjusted carrier-control difference within the band (lm y ~ cnv + age + sex).
# JS 2026-11-02. Outputs: fig_supp_main_effects_by_band.{png,pdf}, main_effects_by_band.csv
source("17_main_with_ukbb.R")
suppressMessages({ library(ggplot2); library(patchwork) })
BANDS <- c(2, 10, 19, 40, 74); BLAB <- c("2-9", "10-18", "19-39", "40-73"); out <- list()
sets <- list(Height = list(X = H, C = Hctl_test), BMI = list(X = copy(H)[, y := BMI], C = copy(Hctl_test)[, y := BMI]),
             `Head circumference` = list(X = dh, C = dh[locus == "No CNV"]))
for (o in names(sets)) {
  X <- sets[[o]]$X; C <- sets[[o]]$C
  for (bb in seq_along(BLAB)) {
    Cb <- C[age >= BANDS[bb] & age < BANDS[bb + 1]]
    out[[length(out) + 1]] <- data.table(outcome = o, band = BLAB[bb], group = "Non-carriers", n = nrow(Cb), n_ukbb = Cb[source == "UK Biobank", .N], mean_z = mean(Cb$y), se = sd(Cb$y) / sqrt(nrow(Cb)))
    for (G in LOCI) {
      Xb <- X[locus == G & age >= BANDS[bb] & age < BANDS[bb + 1]]; if (nrow(Xb) < 10) next
      d <- rbind(Xb, Cb)[, cnv := as.numeric(locus == G)]
      s <- summary(lm(y ~ cnv + age + sex, data = d))$coefficients["cnv", ]
      out[[length(out) + 1]] <- data.table(outcome = o, band = BLAB[bb], group = G, n = nrow(Xb), n_ukbb = Xb[source == "UK Biobank", .N], mean_z = mean(Xb$y), se = sd(Xb$y) / sqrt(nrow(Xb)), shift_vs_controls = s[1], shift_se = s[2], shift_p = s[4])
    }
  }
}
res <- rbindlist(out, fill = TRUE); fwrite(res, file.path(WD, "main_effects_by_band.csv"))
print(res[, .(outcome, band, group, n, n_ukbb, mean_z = round(mean_z, 2), shift = round(shift_vs_controls, 2), p = signif(shift_p, 2))], nrows = 100)
# ---- figure: rows = outcome, columns = CNV; each panel bars by band, grey non-carriers vs blue carriers (as Fig. 6 A'-J')
ctl <- res[group == "Non-carriers"]
pfun <- function(o, G, ylab) {
  d <- rbind(ctl[outcome == o][, .(band, who = "Non-carriers", mean_z, se, n, p = NA_real_)],
             res[outcome == o & group == G, .(band, who = "Carriers", mean_z, se, n, p = shift_p)])
  d[, band := factor(band, BLAB)][, who := factor(who, c("Non-carriers", "Carriers"))]
  bands_present <- d[who == "Carriers", as.character(band)]
  d <- d[band %in% bands_present]
  nl <- d[, .(lab = paste0("n = ", format(n[who == "Carriers"], big.mark = ",")), ymin = min(mean_z - 1.96 * se)), by = band]
  pl <- d[who == "Carriers" & !is.na(p)][, lab := ifelse(p < 0.001, formatC(p, format = "e", digits = 1), signif(p, 2))][, lab := paste0("p = ", lab)][, bold := p < 0.05]
  pl[, ytxt := pmax(mean_z + 1.96 * se, d[band %in% pl$band, max(mean_z + 1.96 * se), by = band]$V1) + 0.08]
  rng <- range(c(d$mean_z - 1.96 * d$se, d$mean_z + 1.96 * d$se, pl$ytxt + 0.15))
  ggplot(d, aes(band, mean_z, fill = who)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_errorbar(aes(ymin = mean_z - 1.96 * se, ymax = mean_z + 1.96 * se), position = position_dodge(0.7), width = 0.2, linewidth = 0.3) +
    geom_text(data = pl, aes(band, ytxt, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 2.3) +
    scale_fill_manual(values = c(`Non-carriers` = "grey60", Carriers = "#3B6FB6"), name = NULL) +
    scale_x_discrete(labels = setNames(paste0(nl$band, "\n", nl$lab), nl$band)) +
    coord_cartesian(ylim = c(min(-0.3, rng[1] - 0.1), max(0.3, rng[2]))) +
    labs(title = G, subtitle = paste0(o, ": mean Z by age band (95% CI)"), x = "Age band (years)", y = ylab) +
    theme_classic(base_size = 9) + theme(legend.position = "top", legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7.5))
}
ylabs <- c(Height = "Height (SD)", BMI = "BMI (SD)", `Head circumference` = "Head circumference (SD)")
# layout (Jon 2026-11-02): rows = CNV, columns = height, head circumference, BMI
COLS <- c("Height", "Head circumference", "BMI")
plots <- list(); for (G in LOCI) for (o in COLS) plots[[length(plots) + 1]] <- pfun(o, G, ylabs[[o]])
fig <- wrap_plots(plots, ncol = 3, byrow = TRUE) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 11))
ggsave(file.path(WD, "fig_supp_main_effects_by_band.png"), fig, width = 10.5, height = 13, dpi = 250, bg = "white")
ggsave(file.path(WD, "fig_supp_main_effects_by_band.pdf"), fig, width = 10.5, height = 13)
cat("supp main-effects figure written\n")
