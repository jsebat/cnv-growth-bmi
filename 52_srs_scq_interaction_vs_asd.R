# CNV x BMI interaction on SRS / SCQ against two reference groups (JS / Claude 2026-11-14; Jon: report the interaction relative to
# controls and relative to ASD?). Box-Cox scale (00_boxcox_srs_scq.R) and ordinal check (scale-free), per CNV and pooled.
#  reference 1: unaffected SSC relatives (as in the figure); reference 2: SSC ASD probands without a CNV, whose scores sit in the
#  same range as carriers, so this contrast does not depend on comparing slopes at opposite ends of a skewed scale.
# Model (as 24 vs_unaff): y ~ [grp] + is_cnv + s(age, by = is_cnv) + s(age, by = BMI) + cnvBMI; cnvBMI = interaction term.
# Output: srs_scq_interaction_vs_controls_and_asd.csv (aggregate statistics only).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages(library(ordinal))
X[, form2 := D$form2[match(paste(pheno, cohort, sex, age, raw, BMI), paste(D$pheno, D$cohort, D$sex, D$age, D$raw, D$BMI))]]
vs_ref <- function(ref, xc) tryCatch({
  z <- rbind(ref, xc); z[, is_cnv := factor(ifelse(grp %in% LOC, "CNV", "ref"), levels = c("ref", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  ex <- if (uniqueN(xc$grp) > 1) "+ grp" else ""
  m <- bam(as.formula(paste("y ~ is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI", ex)), data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]
  c(co[1], co[2], co[4]) }, error = function(e) rep(NA_real_, 3))
ord1 <- function(z) tryCatch({
  z <- copy(z); br <- unique(quantile(z$raw, seq(0, 1, length.out = 41), type = 1)); z[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)]
  z[, is_cnv := as.numeric(grp %in% LOC)][, cnvBMI := is_cnv * BMI]; cl <- if (uniqueN(z$grp[z$grp %in% LOC]) > 1) "grp + " else ""
  m <- clm(as.formula(paste("ry ~", cl, "is_cnv * poly(age, 2) + BMI + cnvBMI")), data = z); co <- coef(summary(m)); c(co["cnvBMI", 1], co["cnvBMI", 2]) }, error = function(e) rep(NA_real_, 2))
ord_p <- function(z) { r <- t(sapply(split(z, z$form2), function(d) if (sum(d$grp %in% LOC) >= 15 && sum(!d$grp %in% LOC) >= 15) ord1(d) else rep(NA_real_, 2)))
  r <- r[!is.na(r[, 1]) & !is.na(r[, 2]), , drop = FALSE]; if (!nrow(r)) return(c(NA, NA, NA)); w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); se <- sqrt(1 / sum(w)); c(b, se, 2 * pnorm(-abs(b / se))) }
out <- list()
for (ph in c("SRS", "SCQ")) for (g in c(LOC, "All CNVs")) for (rf in c("Unaffected relatives", "SSC ASD probands")) {
  xc <- X[pheno == ph & (if (g == "All CNVs") grp %in% LOC else grp == g)]; ref <- X[pheno == ph & grp == rf]
  if (nrow(xc) < 15) next
  e <- vs_ref(ref, xc); o <- ord_p(rbind(ref, xc))
  out[[length(out) + 1]] <- data.table(pheno = ph, group = g, reference = rf, n_cnv = nrow(xc), n_ref = nrow(ref), ref_mean_y = mean(ref$y), cnv_mean_y = mean(xc$y),
    bc_interaction = e[1], bc_se = e[2], bc_p = e[3], ord_interaction = o[1], ord_se = o[2], ord_p = o[3])
}
R2 <- rbindlist(out); fwrite(R2, file.path(HERE, "srs_scq_interaction_vs_controls_and_asd.csv"))
options(width = 220); print(R2[, .(pheno, group, reference = substr(reference, 1, 12), n_cnv, n_ref, ref_y = round(ref_mean_y, 2), cnv_y = round(cnv_mean_y, 2), bc = round(bc_interaction, 3), p_bc = signif(bc_p, 2), ord = round(ord_interaction, 3), p_ord = signif(ord_p, 2))], nrows = 50)
