# BMI vs SRS / SCQ in CNV carriers, stratified by autism polygenic score (below / above the carriers' median), pooled and per CNV
# (JS / Claude 2026-11-14). G2MH carriers with scores (setup of 47; SRS/SCQ on the Box-Cox scale). Points are partial residuals:
# one reference model per trait, y ~ locus + locus:poly(age, 2) + male + BMI * pgs_asd, with the age and sex terms removed
# (centred), so CNV differences in level remain; the pooled panel also removes the CNV means. Lines: linear fit within each half,
# 95% band. Printed: slope (P) in each half, and the continuous BMI x autism score interaction P (robust SE; pooled model has CNV
# terms, per-CNV models age and sex). Output: fig_bmi_srs_by_autism_pgs.{png,pdf} + _table.csv (aggregate statistics only).
Sys.setenv(SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s47 <- readLines(file.path(HERE, "47_bmi_srs_scq_pgs_causal.R")); eval(parse(text = s47[1:grep("^out <- list", s47)]))
if (!"pgs_asd" %in% names(W)) W <- merge(W, g2[, .(key = GUID, pgs_asd)], by = "key")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
HC <- c("Below median" = "#1F6FB4", "Above median" = "#E07B16"); LT <- c("Below median" = "solid", "Above median" = "22")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
pf <- function(p) ifelse(p < 0.001, sprintf("%.0e", p), sprintf("%.2g", p))
P <- list(); TB <- list()
for (ph in c("SRS", "SCQ")) {
  d <- copy(W[pheno == ph]); med <- median(d$pgs_asd); d[, half := factor(fifelse(pgs_asd < med, "Below median", "Above median"), names(HC))]
  m <- lm(y ~ locus + locus:poly(age, 2) + male + BMI * pgs_asd, d); Xm <- model.matrix(m); nu <- grep("poly\\(age|male", colnames(Xm))
  ct <- as.vector(Xm[, nu, drop = FALSE] %*% coef(m)[nu]); d[, y_adj := y - ct + mean(ct)]
  for (g in c("All four CNVs", LOCI)) {
    x <- if (g == "All four CNVs") copy(d)[, y_adj := y_adj - ave(y_adj, locus) + mean(y_adj)] else d[locus == g]
    fi <- if (g == "All four CNVs") lm(y ~ locus + locus:poly(age, 2) + male + BMI * pgs_asd, x) else lm(y ~ poly(age, 2) + male + BMI * pgs_asd, x)
    ip <- rob(fi, "BMI:pgs_asd")
    sl <- rbindlist(lapply(names(HC), function(h) { xh <- x[half == h]; f <- if (g == "All four CNVs") lm(y ~ locus + locus:poly(age, 2) + male + BMI, xh) else lm(y ~ poly(age, 2) + male + BMI, xh)
      r <- rob(f, "BMI"); data.table(half = h, n = nrow(xh), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }))
    TB[[length(TB) + 1]] <- cbind(data.table(trait = ph, group = g, n = nrow(x), inter_b = ip[["b"]], inter_se = ip[["se"]], inter_p = ip[["p"]]), sl)
    fl <- rbindlist(lapply(names(HC), function(h) { xh <- x[half == h]; f <- lm(y_adj ~ BMI, xh); g2r <- seq(quantile(xh$BMI, 0.03), quantile(xh$BMI, 0.97), length.out = 40)
      pr <- predict(f, data.frame(BMI = g2r), se.fit = TRUE); data.table(half = factor(h, names(HC)), BMI = g2r, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }))
    lab <- sprintf("below: %+.2f (P %s)\nabove: %+.2f (P %s)\nBMI x score P = %s", sl$b[1], pf(sl$p[1]), sl$b[2], pf(sl$p[2]), pf(ip[["p"]]))
    P[[paste(ph, g)]] <- ggplot(x, aes(BMI, y_adj, colour = half)) + geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      geom_point(aes(shape = half), size = 1.4, alpha = 0.55, stroke = 0.3) + geom_ribbon(data = fl, aes(BMI, ymin = lo, ymax = hi, fill = half), inherit.aes = FALSE, alpha = 0.13) +
      geom_line(data = fl, aes(BMI, fit, linetype = half), linewidth = 0.9) +
      annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.04, vjust = 1.15, size = 2.6, colour = "grey15", lineheight = 0.95) +
      scale_colour_manual(values = HC, name = "Autism polygenic score") + scale_fill_manual(values = HC, guide = "none") + scale_linetype_manual(values = LT, name = "Autism polygenic score") +
      scale_shape_manual(values = c(16, 17), name = "Autism polygenic score") + coord_cartesian(ylim = range(W[pheno == ph]$y) + c(-0.1, 0.9)) +
      labs(x = "BMI (z)", y = if (g == "All four CNVs") sprintf("%s (Box-Cox z, age- and sex-adjusted)", ph) else NULL, title = if (ph == "SRS") g else NULL, subtitle = sprintf("n = %d", nrow(x))) +
      theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.subtitle = element_text(size = 8, colour = "grey35"), legend.position = "bottom")
  }
}
TB <- rbindlist(TB); fwrite(TB, file.path(HERE, "fig_bmi_srs_by_autism_pgs_table.csv"))
print(TB[, .(trait, group, n, inter = round(inter_b, 3), p_int = signif(inter_p, 2), half, nh = n, slope = round(b, 3), p = signif(p, 2))])
fig <- wrap_plots(P, nrow = 2, byrow = TRUE) + plot_layout(guides = "collect") + plot_annotation(
  title = "BMI and autistic traits in CNV carriers, by autism polygenic score (G2MH; split at the carriers' median)",
  theme = theme(plot.title = element_text(size = 11, face = "bold"))) & theme(legend.position = "bottom")
ggsave(file.path(HERE, "fig_bmi_srs_by_autism_pgs.png"), fig, width = 15, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(HERE, "fig_bmi_srs_by_autism_pgs.pdf"), fig, width = 15, height = 7.2)
