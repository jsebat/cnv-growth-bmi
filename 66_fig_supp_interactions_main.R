# Supplementary figure for the manuscript (JS / Claude 2026-11-14): the two-way (CNV x BMI) and three-way (CNV x BMI x autism score)
# interactions on autistic traits in the combined sample (G2MH + Simons VIP + Searchlight + SSC carriers; Box-Cox scale; reference SSC
# unaffected relatives), pooled and in the most informative strata: sex, CNV, cohort and autism diagnosis.
#  a  SRS vs BMI, carriers and SSC relatives (partial residuals: age, sex and carrier cohort/CNV removed from one reference model)
#  b  SRS vs BMI in genotyped carriers, by autism score below / above the median (partial residuals, one reference model)
#  c  forest: two-way and three-way estimates for SRS and SCQ; sex strata compare carriers with same-sex relatives, and the sex
#     difference is tested in the full model (CNV x BMI x sex; CNV x BMI x score x sex)
# Models as 61/62/64 (two-way bam with age smooths, as Fig. 5k-t; three-way lm with robust SE). Output:
# fig_supp_interactions_main.{png,pdf} + _table.csv (aggregate statistics only).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s64c <- readLines(file.path(HERE, "64c_fig_supp_interaction_forest_age_asd.R")); eval(parse(text = s64c[1:grep("^RELf|^SIBf", s64c)[1]]))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
tm <- function(m, parts) { nm <- names(coef(m)); nm[sapply(strsplit(nm, ":"), function(z) setequal(z, parts))][1] }
rows <- list(); addr <- function(ph, sec, lab, typ, r) if (!is.null(r)) rows[[length(rows) + 1]] <<- data.table(trait = ph, section = sec, label = lab, type = typ, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
T2 <- "Two-way: CNV x BMI"; T3 <- "Three-way: CNV x BMI x autism score"
for (ph in c("SRS", "SCQ")) {
  C2 <- C[pheno == ph]; X3 <- X[pheno == ph]; car2 <- C2[is_cnv == 1]; car3 <- X3[is_cnv == 1]; ref2 <- RELf(C2); ref3 <- RELf(X3)
  addr(ph, "Pooled", "All carriers", T2, vs_ref(car2, ref2)); addr(ph, "Pooled", "All carriers", T3, three(car3, ref3))
  for (s in c("female", "male")) { lab <- if (s == "female") "Female" else "Male"
    addr(ph, "Sex", lab, T2, vs_ref(car2[sex == s], ref2[sex == s])); addr(ph, "Sex", lab, T3, three(car3[sex == s], ref3[sex == s])) }
  # sex difference tested within carriers (carrier BMI slope, male - female; the relatives' own sex difference is small and reported
  # separately in sex_check_bmi_slopes.csv), and the three-way analogue (BMI x score x sex within carriers)
  # sex difference = male - female difference of the carrier BMI slopes (two-way) or of the within-carrier BMI x score terms
  # (three-way), each fitted separately by sex (no shared age, cohort or CNV terms); SE from the two independent fits
  dif <- function(f, m) { b <- m[["b"]] - f[["b"]]; se <- sqrt(f[["se"]]^2 + m[["se"]]^2); c(b = b, se = se, p = 2 * pnorm(-abs(b / se))) }
  sf <- slope(car2[sex == "female"]); sm <- slope(car2[sex == "male"])
  addr(ph, "Sex", "Sex difference in carriers (male - female)", T2, c(n = nrow(car2), dif(sf, sm)))
  w3 <- function(d) rob(lm(y ~ cg + poly(age, 2) + BMI * pgs, d), "BMI:pgs")
  addr(ph, "Sex", "Sex difference in carriers (male - female)", T3, c(n = nrow(car3), dif(w3(car3[sex == "female"]), w3(car3[sex == "male"]))))
  for (L in LOC) { if (nrow(car2[locus == L]) >= 15) addr(ph, "CNV", L, T2, vs_ref(car2[locus == L], ref2)); addr(ph, "CNV", L, T3, three(car3[locus == L], ref3)) }
  for (k in c("G2MH", "SVIP", "Searchlight")) { if (nrow(car2[cohort == k]) >= 15) addr(ph, "Cohort", COHL[[k]], T2, vs_ref(car2[cohort == k], ref2)); addr(ph, "Cohort", COHL[[k]], T3, three(car3[cohort == k], ref3)) }
  for (k in c(0L, 1L)) { lab <- if (k == 1L) "Autism diagnosis" else "No autism diagnosis"
    addr(ph, "Autism diagnosis", lab, T2, vs_ref(car2[asd2 %in% k], ref2)); addr(ph, "Autism diagnosis", lab, T3, three(car3[asd2 %in% k], ref3)) }
  cat(ph, "done\n")
}
R <- rbindlist(rows); fwrite(R, file.path(HERE, "fig_supp_interactions_main_table.csv"))
options(width = 200); print(R[, .(trait, section, label, type = substr(type, 1, 9), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
# ---------------------------------------------------------------- panels a, b (SRS)
th <- theme_classic(base_size = 10.5) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "bottom", legend.key.width = unit(1.2, "cm"))
C2 <- C[pheno == "SRS"]; z <- rbind(C2[is_cnv == 1], RELf(C2))
m <- lm(y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI, z); Xm <- model.matrix(m); nu <- grep("poly\\(age|male|^cg", colnames(Xm))
cf <- coef(m); cf[is.na(cf)] <- 0; ct <- as.vector(Xm[, nu, drop = FALSE] %*% cf[nu]); z[, ya := y - ct + ave(ct, is_cnv)][, grp := factor(fifelse(is_cnv == 1, "CNV carriers", "SSC unaffected relatives"), c("CNV carriers", "SSC unaffected relatives"))]
fl <- z[, { f <- lm(ya ~ BMI); g <- seq(quantile(BMI, 0.02), quantile(BMI, 0.98), length.out = 40); pr <- predict(f, data.frame(BMI = g), se.fit = TRUE); .(BMI = g, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }, by = grp]
r2 <- R[trait == "SRS" & section == "Pooled" & type == T2]
GC <- c("CNV carriers" = "#7B241C", "SSC unaffected relatives" = "grey45")
pa <- ggplot(z, aes(BMI, ya, colour = grp)) + geom_point(data = z[is_cnv == 0], size = 0.5, alpha = 0.08) + geom_point(data = z[is_cnv == 1], size = 0.9, alpha = 0.45) +
  geom_ribbon(data = fl, aes(BMI, ymin = lo, ymax = hi, fill = grp), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = fl, aes(BMI, fit), linewidth = 1) +
  scale_colour_manual(values = GC, name = NULL) + scale_fill_manual(values = GC, guide = "none") + coord_cartesian(xlim = c(-3, 4.2)) +
  annotate("text", x = -3, y = Inf, hjust = 0, vjust = 1.2, size = 3, label = sprintf("CNV x BMI = %+.3f (P = %.2g)\ncarriers n = %d, relatives n = %d", r2$b, r2$p, sum(z$is_cnv == 1), sum(z$is_cnv == 0))) +
  labs(x = "BMI (z)", y = "SRS (Box-Cox z, adjusted)", title = "a  Two-way: CNV x BMI (SRS)", subtitle = "Partial residuals: age, sex and cohort/CNV removed") + th
X3 <- X[pheno == "SRS" & is_cnv == 1]; med <- median(X3$pgs); X3[, half := factor(fifelse(pgs < med, "Autism score below median", "Autism score above median"), c("Autism score below median", "Autism score above median"))]
m3 <- lm(y ~ cg + poly(age, 2) + male + BMI * pgs, X3); Xm <- model.matrix(m3); nu <- grep("poly\\(age|male|^cg", colnames(Xm))
cf <- coef(m3); cf[is.na(cf)] <- 0; ct <- as.vector(Xm[, nu, drop = FALSE] %*% cf[nu]); X3[, ya := y - ct + mean(ct)]
f3 <- X3[, { f <- lm(ya ~ BMI); g <- seq(quantile(BMI, 0.03), quantile(BMI, 0.97), length.out = 40); pr <- predict(f, data.frame(BMI = g), se.fit = TRUE); .(BMI = g, fit = pr$fit, lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit) }, by = half]
r3 <- R[trait == "SRS" & section == "Pooled" & type == T3]
HC <- c("Autism score below median" = "#1F6FB4", "Autism score above median" = "#E07B16")
pb <- ggplot(X3, aes(BMI, ya, colour = half)) + geom_point(aes(shape = half), size = 1, alpha = 0.45) +
  geom_ribbon(data = f3, aes(BMI, ymin = lo, ymax = hi, fill = half), inherit.aes = FALSE, alpha = 0.18) + geom_line(data = f3, aes(BMI, fit, linetype = half), linewidth = 1) +
  scale_colour_manual(values = HC, name = NULL) + scale_fill_manual(values = HC, guide = "none") + scale_linetype_manual(values = c("solid", "22"), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
  coord_cartesian(xlim = c(-3, 4.2)) + annotate("text", x = -3, y = Inf, hjust = 0, vjust = 1.2, size = 3, label = sprintf("CNV x BMI x score = %+.3f (P = %.2g)\ngenotyped carriers n = %d", r3$b, r3$p, nrow(X3))) +
  labs(x = "BMI (z)", y = "SRS (Box-Cox z, adjusted)", title = "b  Three-way: CNV x BMI x autism score (SRS)", subtitle = "Genotyped carriers (G2MH, Simons VIP, SSC), split at the median autism score") + th
# ---------------------------------------------------------------- panel c forest
F <- R[!is.na(b)]
ordl <- c("All carriers", "Female", "Male", "Sex difference in carriers (male - female)", LOC, "G2MH", "Simons VIP", "Simons Searchlight", "No autism diagnosis", "Autism diagnosis")
F[, label := factor(label, levels = rev(ordl))][, section := factor(section, c("Pooled", "Sex", "CNV", "Cohort", "Autism diagnosis"))][, trait := factor(trait, c("SRS", "SCQ"))][, type := factor(type, c(T2, T3))]
OC <- c(SRS = "#C0392B", SCQ = "#2E86C1"); pd <- position_dodge(width = 0.6)
pc <- ggplot(F, aes(b, label, colour = trait, group = trait)) + geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd, linewidth = 0.55) +
  geom_point(aes(shape = p < 0.05), size = 2, position = pd) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("P >= 0.05", "P < 0.05"), name = NULL) +
  geom_text(aes(x = Inf, label = paste0("n=", n)), position = pd, hjust = 1.05, size = 2.3, show.legend = FALSE) +
  facet_grid(section ~ type, scales = "free", space = "free_y") + scale_colour_manual(values = OC, name = NULL) +
  labs(x = "Interaction estimate (95% CI)", y = NULL, title = "c  Interaction estimates, pooled and by stratum (SRS and SCQ)",
       subtitle = "Reference SSC unaffected relatives (same sex for the sex strata; sex difference tested within carriers); n = carriers; three-way: genotyped carriers") +
  th + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0), strip.text.x = element_text(face = "bold"), panel.spacing.y = unit(0.5, "lines"), legend.position = "top")
fig <- (pa | pb) / pc + plot_layout(heights = c(1, 1.55))
ggsave(file.path(HERE, "fig_supp_interactions_main.png"), fig, width = 12, height = 13, dpi = 250, bg = "white")
ggsave(file.path(HERE, "fig_supp_interactions_main.pdf"), fig, width = 12, height = 13)
