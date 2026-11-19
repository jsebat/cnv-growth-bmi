# Scale sensitivity of the BMI effects on SRS and SCQ (JS / Claude 2026-11-14). Jon: the biobank analyses were moved to a Box-Cox
# scale to avoid spurious interactions from extreme shifts in height/BMI; could the same artefact affect SRS and SCQ?
# Concern: SRS/SCQ are raw totals z-scored linearly against unaffected relatives (21_pheno_bmi_gam.R zref). The raw totals are
# right-skewed in relatives and carriers sit several SD above them; on a skewed scale whose spread grows with the mean, any
# predictor has a larger slope in the high-scoring group, so CNV x BMI could be a scale artefact.
# Tests (same samples and models as 23_unaffected_controls.R contrast(): y ~ [locus] + is_cnv + s(age, by = is_cnv) +
#   s(age, by = BMI) + cnvBMI; within carriers bmi_main(): y ~ [locus] + BMI + s(age)):
#   S0 original linear z vs relatives
#   S1 log(raw + 1), z vs relatives
#   S2 Box-Cox (lambda by profile likelihood in relatives, raw + 1 ~ age + age^2 + sex; SRS child and adult forms separately), z vs relatives
#   S3 BMI trimmed to |z| <= 3 (original y)
#   S4 ordinal (cumulative logit, ordinal::clm) on the raw total: invariant to ANY monotone transformation of the score.
#      raw ~ [locus] + is_cnv * poly(age, 3) + BMI + BMI:age + cnvBMI (+ separate thresholds by SRS form); log-odds per SD BMI.
# Diagnostics: skewness and mean/SD of the raw totals by group; mean-SD relation across BMI bins in relatives.
# Output: srs_scq_scale_sensitivity.csv (aggregate statistics only).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
Sys.setenv(UNAFF_DRAW_FIG = "FALSE"); source(file.path(HERE, "23_unaffected_controls.R"))
suppressPackageStartupMessages(library(ordinal))
skew <- function(x) { x <- x[!is.na(x)]; mean((x - mean(x))^3) / sd(x)^3 }
# recover the SRS form from the linear z (child form referenced to siblings, adult form to parents)
rs <- list(child = ssc[form == "child" & role == "sibling", srs], adult = ssc[form == "adult" & role == "parent", srs], scq = ssc[role == "sibling", scq])
S <- A[pheno %in% c("SRS", "SCQ")]
S[pheno == "SRS", form := fifelse(abs(y - zref(raw, rs$child)) < 1e-8, "child", fifelse(abs(y - zref(raw, rs$adult)) < 1e-8, "adult", NA_character_))]
S[pheno == "SCQ", form := "scq"]
stopifnot(!anyNA(S$form))
S[, ref := locus == "Unaffected controls"]
cat("n by pheno x form x group:\n"); print(dcast(S[, .N, by = .(pheno, form, g = fifelse(ref, "relatives", "carriers"))], pheno + form ~ g, value.var = "N"))
# ---------------------------------------------------------------- diagnostics
dg <- S[, .(n = .N, mean_raw = mean(raw), sd_raw = sd(raw), skew_raw = skew(raw), skew_log = skew(log(raw + 1))), by = .(pheno, form, group = fifelse(ref, "Unaffected relatives", locus))][order(pheno, form, group)]
cat("\nraw-score distribution:\n"); print(dg[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
mv <- S[ref == TRUE, .(n = .N, mean_raw = mean(raw), sd_raw = sd(raw)), by = .(pheno, form, bmi_bin = cut(BMI, c(-Inf, -1, 0, 1, Inf)))][order(pheno, form, bmi_bin)]
cat("\nrelatives, mean and SD of the raw score by BMI bin (spread growing with the mean = scale artefact risk):\n"); print(mv[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
# ---------------------------------------------------------------- rescaled outcomes
lam <- list()
for (f in unique(S$form)) {
  r <- S[ref == TRUE & form == f]; r[, `:=`(yy = raw + 1, a2 = age^2)]
  b <- MASS::boxcox(yy ~ age + a2 + sex, data = r, lambda = seq(-2, 2, 0.01), plotit = FALSE); lam[[f]] <- b$x[which.max(b$y)]
  bc <- function(x) if (abs(lam[[f]]) < 1e-8) log(x + 1) else ((x + 1)^lam[[f]] - 1) / lam[[f]]
  S[form == f, `:=`(y_log = zref(log(raw + 1), log(r$raw + 1)), y_bc = zref(bc(raw), bc(r$raw)))]
}
cat("\nBox-Cox lambda by form:", paste(names(lam), round(unlist(lam), 2), collapse = "; "), "\n")
cat("skewness in relatives after Box-Cox:\n"); print(S[ref == TRUE, .(skew_lin = round(skew(y), 2), skew_bc = round(skew(y_bc), 2)), by = .(pheno, form)])
# ---------------------------------------------------------------- contrasts
# S4: per form (SRS child, SRS adult, SCQ) on the raw total binned into <= 40 ordered levels (quantiles of the pooled sample; monotone,
# so the fit stays invariant to any monotone transformation), then inverse-variance combined across SRS forms.
ord1 <- function(z, within) tryCatch({
  z <- copy(z); br <- unique(quantile(z$raw, seq(0, 1, length.out = 41), type = 1)); z[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)]
  z[, is_cnv := as.numeric(locus != "Unaffected controls")][, cnvBMI := is_cnv * BMI]
  cl <- if (uniqueN(z$locus[z$locus != "Unaffected controls"]) > 1) "locus + " else ""
  f <- if (within) paste("ry ~", cl, "poly(age, 2) + BMI") else paste("ry ~", cl, "is_cnv * poly(age, 2) + BMI + cnvBMI")
  m <- clm(as.formula(f), data = z); co <- coef(summary(m)); v <- if (within) "BMI" else "cnvBMI"; c(co[v, 1], co[v, 2])
}, error = function(e) { message("clm: ", conditionMessage(e)); rep(NA_real_, 2) })
ord_fit <- function(z, within = FALSE) {
  r <- t(sapply(split(z, z$form), function(d) if (sum(d$locus != "Unaffected controls") >= 20 || within) ord1(d, within) else rep(NA_real_, 2))); r <- r[!is.na(r[, 1]) & !is.na(r[, 2]), , drop = FALSE]
  if (!nrow(r)) return(rep(NA_real_, 3)); w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); se <- sqrt(1 / sum(w)); c(b, se, 2 * pnorm(-abs(b / se))) }
GR <- c(setNames(as.list(LOCI), LOCI), list(`All CNVs` = LOCI))
out <- list()
for (ph in c("SRS", "SCQ")) for (g in names(GR)) {
  xc <- S[pheno == ph & locus %in% GR[[g]]]; ctl <- S[pheno == ph & ref == TRUE]
  for (sc in c("S0 linear z (paper)", "S1 log", "S2 Box-Cox", "S3 BMI trimmed |z|<=3")) {
    v <- switch(substr(sc, 1, 2), S0 = "y", S1 = "y_log", S2 = "y_bc", S3 = "y")
    a <- copy(xc)[, y := get(v)]; b <- copy(ctl)[, y := get(v)]
    if (substr(sc, 1, 2) == "S3") { a <- a[abs(BMI) <= 3]; b <- b[abs(BMI) <= 3] }
    e <- contrast(b, a); w <- bmi_main(copy(a)); wc <- bmi_main(copy(b))
    out[[length(out) + 1]] <- data.table(pheno = ph, group = g, scale = sc, n_cnv = nrow(a), n_rel = nrow(b), carrier_slope = w[1], carrier_slope_p = w[3], relative_slope = wc[1], relative_slope_p = wc[3], cnv_x_bmi = e[3], cnv_x_bmi_se = e[4], cnv_x_bmi_p = e[5])
  }
  o <- ord_fit(rbind(xc, ctl)); ow <- ord_fit(xc, within = TRUE); orr <- ord_fit(ctl, within = TRUE)
  out[[length(out) + 1]] <- data.table(pheno = ph, group = g, scale = "S4 ordinal (log-odds per SD BMI)", n_cnv = nrow(xc), n_rel = nrow(ctl), carrier_slope = ow[1], carrier_slope_p = ow[3], relative_slope = orr[1], relative_slope_p = orr[3], cnv_x_bmi = o[1], cnv_x_bmi_se = o[2], cnv_x_bmi_p = o[3])
  cat(ph, g, "done\n")
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "srs_scq_scale_sensitivity.csv"))
fwrite(dg, file.path(HERE, "srs_scq_scale_diagnostics.csv"))
options(width = 250); print(R[, .(pheno, group, scale, n_cnv, n_rel, carrier = round(carrier_slope, 3), p_car = signif(carrier_slope_p, 2), rel = round(relative_slope, 3), p_rel = signif(relative_slope_p, 2), cxb = round(cnv_x_bmi, 3), se = round(cnv_x_bmi_se, 3), p_cxb = signif(cnv_x_bmi_p, 2))], nrows = 200)
