# Same stratification as 55 (BMI vs SRS / SCQ by autism polygenic score, below / above the group median) in SSC without CNVs:
# autistic probands, designated siblings and (SRS only) parents (JS / Claude 2026-11-14). Autism score from the April 2026 ASD GWAS
# leaving SSC out (ssc_pgs_scores.tsv, residualised on 10 PCs, z within SSC). SRS / SCQ on the Box-Cox scale (00_boxcox_srs_scq.R).
# Points: partial residuals from one reference model per panel, y ~ poly(age, 2) + male [+ SRS form] + BMI * pgs_asd, age, sex and
# form terms removed (centred). Lines: linear fit within each half, 95% band. Printed: slope (P) per half and the continuous
# BMI x autism score interaction P (robust SE). Output: fig_bmi_srs_by_autism_pgs_ssc.{png,pdf} + _table.csv (aggregate only).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(sandwich) })
X[, form2 := D$form2[match(paste(pheno, cohort, sex, age, raw, BMI), paste(D$pheno, D$cohort, D$sex, D$age, D$raw, D$BMI))]]
pg <- fread(file.path(DATA, "ssc_pgs_scores.tsv"))
S <- merge(X[cohort == "SSC" & pheno %in% c("SRS", "SCQ") & grp %in% c("SSC ASD probands", "Unaffected relatives")], pg[, .(individual, pgs_asd)], by.x = "id", by.y = "individual")
S[, panel := fcase(grp == "SSC ASD probands", "Autistic probands", role == "sibling", "Unaffected siblings", role == "parent", "Parents")][, male := as.numeric(sex == "male")]
cat("n by trait x panel:\n"); print(dcast(S[, .N, by = .(pheno, panel)], panel ~ pheno, value.var = "N"))
rob <- function(m, v) { V <- vcovHC(m, type = "HC1"); b <- unname(coef(m)[v]); s <- sqrt(V[v, v]); c(b = b, se = s, p = 2 * pnorm(-abs(b / s))) }
HC <- c("Below median" = "#1F6FB4", "Above median" = "#E07B16"); LT <- c("Below median" = "solid", "Above median" = "22")
pf <- function(p) ifelse(p < 0.001, sprintf("%.0e", p), sprintf("%.2g", p))
PAN <- c("Autistic probands", "Unaffected siblings", "Parents")
P <- list(); TB <- list()
for (ph in c("SRS", "SCQ")) for (g in PAN) {
  x <- copy(S[pheno == ph & panel == g]); if (nrow(x) < 50) next
  med <- median(x$pgs_asd); x[, half := factor(fifelse(pgs_asd < med, "Below median", "Above median"), names(HC))]
  cov <- paste("poly(age, 2) + male", if (uniqueN(x$form2) > 1) "+ form2" else "")
  m <- lm(as.formula(paste("y ~", cov, "+ BMI * pgs_asd")), x); Xm <- model.matrix(m); nu <- grep("poly\\(age|male|form2", colnames(Xm))
  ct <- as.vector(Xm[, nu, drop = FALSE] %*% coef(m)[nu]); x[, y_adj := y - ct + mean(ct)]
  ip <- rob(m, "BMI:pgs_asd")
  sl <- rbindlist(lapply(names(HC), function(h) { xh <- x[half == h]; cvh <- paste("poly(age, 2) + male", if (uniqueN(xh$form2) > 1) "+ form2" else ""); r <- rob(lm(as.formula(paste("y ~", cvh, "+ BMI")), xh), "BMI"); data.table(half = h, n_half = nrow(xh), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }))
  TB[[length(TB) + 1]] <- cbind(data.table(trait = ph, group = g, n = nrow(x), score_main = coef(m)[["pgs_asd"]], inter_b = ip[["b"]], inter_se = ip[["se"]], inter_p = ip[["p"]]), sl)
  fl <- rbindlist(lapply(names(HC), function(h) { xh <- x[half == h]; f <- lm(y_adj ~ BMI, xh); gr <- seq(quantile(xh$BMI, 0.03), quantile(xh$BMI, 0.97), length.out = 40)
    pr <- predict(f, data.frame(BMI = gr), se.fit = TRUE); data.table(half = factor(h, names(HC)), BMI = gr, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }))
  lab <- sprintf("below: %+.2f (P %s)\nabove: %+.2f (P %s)\nBMI x score P = %s", sl$b[1], pf(sl$p[1]), sl$b[2], pf(sl$p[2]), pf(ip[["p"]]))
  yl <- quantile(S[pheno == ph]$y, c(0.005, 0.995)) + c(-0.2, 1.2)
  P[[paste(ph, g)]] <- ggplot(x, aes(BMI, y_adj, colour = half)) + geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(shape = half), size = 0.9, alpha = 0.3, stroke = 0.2) + geom_ribbon(data = fl, aes(BMI, ymin = lo, ymax = hi, fill = half), inherit.aes = FALSE, alpha = 0.18) +
    geom_line(data = fl, aes(BMI, fit, linetype = half), linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.04, vjust = 1.15, size = 2.6, colour = "grey15", lineheight = 0.95) +
    scale_colour_manual(values = HC, name = "Autism polygenic score") + scale_fill_manual(values = HC, guide = "none") + scale_linetype_manual(values = LT, name = "Autism polygenic score") +
    scale_shape_manual(values = c(16, 17), name = "Autism polygenic score") + coord_cartesian(xlim = c(-3.2, 4.2), ylim = yl) +
    labs(x = "BMI (z)", y = if (g == "Autistic probands") sprintf("%s (Box-Cox z)", ph) else NULL, title = if (ph == "SRS") g else NULL, subtitle = sprintf("n = %d", nrow(x))) +
    theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.subtitle = element_text(size = 8, colour = "grey35"), legend.position = "bottom")
}
TB <- rbindlist(TB); fwrite(TB, file.path(HERE, "fig_bmi_srs_by_autism_pgs_ssc_table.csv"))
print(TB[, .(trait, group, n, main = round(score_main, 3), inter = round(inter_b, 3), p_int = signif(inter_p, 2), half, n_half, slope = round(b, 3), p = signif(p, 2))])
lay <- c(P[c("SRS Autistic probands", "SRS Unaffected siblings", "SRS Parents")], P[c("SCQ Autistic probands", "SCQ Unaffected siblings")], list(plot_spacer()))
fig <- wrap_plots(lay, nrow = 2, byrow = TRUE) + plot_layout(guides = "collect") + plot_annotation(
  title = "BMI and autistic traits in the Simons Simplex Collection (no CNV), by autism polygenic score (split at each group's median)", subtitle = "SRS and SCQ: Box-Cox z vs unaffected relatives, age- and sex-adjusted partial residuals",
  theme = theme(plot.title = element_text(size = 11, face = "bold"))) & theme(legend.position = "bottom")
ggsave(file.path(HERE, "fig_bmi_srs_by_autism_pgs_ssc.png"), fig, width = 10.5, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(HERE, "fig_bmi_srs_by_autism_pgs_ssc.pdf"), fig, width = 10.5, height = 7.2)
