# Polygenic score x BMI interaction on SRS / SCQ in SSC, without CNVs (JS / Claude 2026-11-14). Jon: is the BMI x genotype
# interaction specific to rare variants of large effect, or also seen with polygenic scores? Autism score from the April 2026 ASD
# GWAS leaving SSC out (asd2026_loSSC); childhood obesity, childhood BMI, adult BMI and height scores as in 30a (all residualised on
# 10 PCs and z-scored within SSC; ssc_pgs_scores.tsv). SRS / SCQ on the Box-Cox scale (00_boxcox_srs_scq.R), as in Fig. 5.
# Groups: SSC probands (no CNV) and unaffected relatives (SRS: siblings + parents; SCQ: siblings); both pooled with a group term.
# Models: y ~ s(age) + sex [+ SRS form] + BMI * PGS (continuous PGS, per SD; and top decile vs rest, the closest analogue of a
# carrier/non-carrier contrast); ordinal (cumulative logit, raw total binned <= 40 levels, per form, inverse-variance combined).
# For scale: the pooled CNV x BMI interaction is 0.10 (SRS) and 0.10 (SCQ) s.d. per s.d. BMI (Fig. 5).
# Output: ssc_pgs_x_bmi.csv (aggregate statistics only).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages(library(ordinal))
X[, form2 := D$form2[match(paste(pheno, cohort, sex, age, raw, BMI), paste(D$pheno, D$cohort, D$sex, D$age, D$raw, D$BMI))]]
pg <- fread(file.path(DATA, "ssc_pgs_scores.tsv"))
S <- merge(X[cohort == "SSC" & grp %in% c("SSC ASD probands", "Unaffected relatives") & pheno %in% c("SRS", "SCQ")], pg, by.x = "id", by.y = "individual")
cat("SSC with scores:\n"); print(S[, .N, by = .(pheno, grp)])
PGS <- c(pgs_asd = "Autism (lo-SSC)", pgs_cobes = "Childhood obesity", pgs_cbmi = "Childhood BMI", pgs_bmi = "Adult BMI", pgs_ht = "Height")
fit <- function(d, v, top) tryCatch({
  d <- copy(d); d[, g := if (top) as.numeric(get(v) >= quantile(get(v), 0.9)) else get(v)][, gB := g * BMI]
  ex <- c(if (uniqueN(d$grp) > 1) "grp", if (uniqueN(d$form2) > 1) "form2")
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + sex + BMI + g + gB", if (length(ex)) paste("+", paste(ex, collapse = " + ")) else "")), data = d, method = "REML")
  co <- summary(m)$p.table; c(co["g", 1], co["g", 4], co["BMI", 1], co["gB", 1], co["gB", 2], co["gB", 4]) }, error = function(e) rep(NA_real_, 6))
ord <- function(d, v) { r <- t(sapply(split(d, d$form2), function(z) tryCatch({
    z <- copy(z); br <- unique(quantile(z$raw, seq(0, 1, length.out = 41), type = 1)); z[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)][, `:=`(g = get(v), a2 = age^2)][, gB := g * BMI]
    m <- clm(as.formula(paste("ry ~ age + a2 + sex + BMI + g + gB", if (uniqueN(z$grp) > 1) "+ grp" else "")), data = z); co <- coef(summary(m)); c(co["gB", 1], co["gB", 2]) }, error = function(e) c(NA, NA))))
  r <- r[!is.na(r[, 1]), , drop = FALSE]; if (!nrow(r)) return(c(NA, NA)); w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); c(b, 2 * pnorm(-abs(b / sqrt(1 / sum(w))))) }
out <- list()
for (ph in c("SRS", "SCQ")) for (gr in c("SSC ASD probands", "Unaffected relatives", "both")) for (v in names(PGS)) {
  d <- S[pheno == ph & (gr == "both" | grp == gr)]
  a <- fit(d, v, FALSE); t <- fit(d, v, TRUE); o <- ord(d, v)
  out[[length(out) + 1]] <- data.table(pheno = ph, group = gr, pgs = PGS[[v]], n = nrow(d), pgs_main = a[1], pgs_main_p = a[2], bmi_slope = a[3],
    pgs_x_bmi = a[4], pgs_x_bmi_se = a[5], pgs_x_bmi_p = a[6], top10_x_bmi = t[4], top10_x_bmi_se = t[5], top10_x_bmi_p = t[6], ord_pgs_x_bmi = o[1], ord_p = o[2])
}
R3 <- rbindlist(out); fwrite(R3, file.path(HERE, "ssc_pgs_x_bmi.csv"))
options(width = 230); print(R3[, .(pheno, group = substr(group, 1, 12), pgs, n, main = round(pgs_main, 3), p_main = signif(pgs_main_p, 2), bmi = round(bmi_slope, 3), pxb = round(pgs_x_bmi, 3), se = round(pgs_x_bmi_se, 3), p = signif(pgs_x_bmi_p, 2), top10 = round(top10_x_bmi, 3), p_top = signif(top10_x_bmi_p, 2), p_ord = signif(ord_p, 2))], nrows = 60)
