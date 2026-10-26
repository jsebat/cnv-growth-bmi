# Supplementary figure (JS / Claude 2026-10-26): height, head circumference and BMI in CNV carriers with versus without an autism
# diagnosis, within each CNV group, pooled overgrowth / undergrowth groups, and the non-carrier comparison (SSC probands with ASD
# versus SSC designated siblings). Restores the ASD-by-CNV panels (F/L of the earlier UKBB main figure, 10_supplement.R) on the
# QC'd data build, and adds BMI. Bars: mean +/- SE of the Z score by ASD status; printed: the age- and sex-adjusted ASD effect
# (GAM y ~ asd + s(age) + sex [+ cohort] [+ locus for pooled groups], REML; Delta in SD, p). Ages 2-65, latest visit with the
# outcome and BMI, as in the main figure. Outputs: supp_asd_by_cnv_tests.csv, fig_supp_asd_by_cnv.{png,pdf}, fig_supp_asd_by_cnv_legend.txt
Sys.setenv(UNAFF_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))   # A (growth rows with asd), gmc, ctl_g, mk_growth, LOCI, POOLS
suppressPackageStartupMessages({ library(mgcv); library(ggplot2); library(patchwork) })
# non-carrier comparison: SSC probands (ASD) vs SSC designated siblings (no ASD), same build as the carriers
pro_g <- gmc[cohort == "SSC" & CNV == "No CNV" & grepl("[.]p[0-9]$", sampleID), .(id, cohort, sampleID, GUID, CNV, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd = 1L)][, locus := "Non-carriers (SSC)"]
sib_g <- ctl_g[cohort == "SSC" & grepl("[.]s[0-9]$", sampleID)][, `:=`(asd = 0L, locus = "Non-carriers (SSC)")][, .(id, cohort, sampleID, GUID, CNV, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd, locus)]
NC <- rbind(mk_growth(rbind(pro_g, sib_g), "HEIGHT_raw", "Height"), mk_growth(rbind(pro_g, sib_g), "HEADCIRCUM_raw", "Head circumference"))
X <- rbind(A[locus %in% LOCI & !is.na(asd), .(id, cohort, locus, sex, age, BMI, y, pheno, asd)], NC[, .(id, cohort, locus, sex, age, BMI, y, pheno, asd)])
# BMI as an outcome: one row per individual (from the Height rows, which carry BMI Z for everyone with height and BMI)
XB <- X[pheno == "Height", .(id, cohort, locus, sex, age, BMI, y = BMI, pheno = "BMI", asd)]
X <- rbind(X, XB)[, asd_f := factor(ifelse(asd == 1L, "ASD", "No ASD"), c("No ASD", "ASD"))]
GROUPS <- c(setNames(as.list(LOCI), LOCI), POOLS[c("Overgrowth CNVs", "Undergrowth CNVs")], list(`Non-carriers (SSC)` = "Non-carriers (SSC)"))
PH <- c("Height", "Head circumference", "BMI")
fit <- function(x) tryCatch({
  x <- copy(x)[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))]
  cov <- paste(c(if (nlevels(x$cohort) > 1) "+ cohort", if (nlevels(x$locus) > 1) "+ locus"), collapse = " ")
  m <- gam(as.formula(paste("y ~ asd + s(age, k = 8) + sex", cov)), data = x, method = "REML"); co <- summary(m)$p.table["asd", ]
  data.table(delta = unname(co[1]), se = unname(co[2]), p = unname(co[4])) }, error = function(e) data.table(delta = NA_real_, se = NA_real_, p = NA_real_))
res <- rbindlist(lapply(PH, function(ph) rbindlist(lapply(names(GROUPS), function(g) {
  x <- X[pheno == ph & locus %in% GROUPS[[g]]]; if (sum(x$asd == 1L) < 5 || sum(x$asd == 0L) < 5) return(NULL)
  cbind(data.table(pheno = ph, group = g, n_noasd = sum(x$asd == 0L), n_asd = sum(x$asd == 1L), mean_noasd = x[asd == 0L, mean(y)], se_noasd = x[asd == 0L, sd(y) / sqrt(.N)], mean_asd = x[asd == 1L, mean(y)], se_asd = x[asd == 1L, sd(y) / sqrt(.N)],
                   age_med_noasd = x[asd == 0L, median(age)], age_med_asd = x[asd == 1L, median(age)]), fit(x)) }))))
fwrite(res, file.path(HERE, "supp_asd_by_cnv_tests.csv"))
cat("== ASD vs no ASD, age- and sex-adjusted (SD)\n"); print(res[, .(pheno, group, n_noasd, n_asd, delta = round(delta, 2), se = round(se, 2), p = signif(p, 2))], nrows = 60)
# ---------------------------------------------------------------- figure
COLS <- c("No ASD" = "#8C9A9A", "ASD" = "#4B6EA9")
YL <- c(Height = "Height (SD)", `Head circumference` = "Head circumference (SD)", BMI = "BMI (Z)")
bars <- rbind(res[, .(pheno, group, grp = "No ASD", m = mean_noasd, se = se_noasd, n = n_noasd)], res[, .(pheno, group, grp = "ASD", m = mean_asd, se = se_asd, n = n_asd)])
bars[, grp := factor(grp, names(COLS))][, group := factor(group, names(GROUPS))]
pan <- function(ph, tag) {
  b <- bars[pheno == ph]; r <- res[pheno == ph][, group := factor(group, names(GROUPS))]
  top <- b[, max(m + se)]; rng <- diff(range(c(b$m - b$se, b$m + b$se, 0))); r[, `:=`(ytxt = top + 0.16 * rng, lab = gsub("= < ", "< ", sprintf("\u0394 %+.2f\np = %s", delta, pfx(p)), fixed = TRUE), bold = p < 0.05)]
  ggplot(b, aes(group, m, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(0.8), width = 0.75) + geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(data = b, aes(label = n, y = pmin(m, 0) - 0.06 * rng, group = grp), position = position_dodge(0.8), size = 2.6, colour = "grey30") +
    geom_text(data = r, aes(group, ytxt, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = COLS) + coord_cartesian(ylim = c(min(0, b[, min(m - se)]) - 0.12 * rng, top + 0.32 * rng)) +
    labs(x = NULL, y = YL[[ph]], fill = NULL, title = ph, tag = tag) + theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5), plot.tag = element_text(size = 16, face = "bold"))
}
fig <- (pan("Height", "A") | pan("Head circumference", "B") | pan("BMI", "C")) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & plot_annotation(title = "Height, head circumference and BMI by autism diagnosis within CNV group", theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(HERE, "fig_supp_asd_by_cnv.png"), fig, width = 15, height = 5.2, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_supp_asd_by_cnv.pdf"), fig, width = 15, height = 5.2, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_supp_asd_by_cnv.pdf"), fig, width = 15, height = 5.2))
# ---------------------------------------------------------------- sensitivity: ASD effects with sex and polygenic scores as covariates
# Same PGS as Tables S4-S6 (G2MH WGS carriers: Anjali's scores + autism 2026 full weights; SSC: scored on Expanse, autism lo-SSC),
# residualised on 10 PCs and z-scored within cohort. Models add each PGS, then all four, to y ~ asd + s(age) + sex [+ cohort/locus];
# rows are restricted to individuals with all four scores so that the unadjusted and adjusted estimates use the same people.
AD <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data")
rdw <- function(p) { x <- fread(p); setnames(x, "#IID", "IID"); x[, .(GUID = sub("^[0-9]+-", "", IID), score = SCORE1_AVG)] }
g2 <- Reduce(function(a, b) merge(a, b, by = "GUID"), list(setnames(rdw(file.path(AD, "pgs_scores/scores_wgs/childhood_obesity.sscore")), "score", "cobes"), setnames(rdw(file.path(AD, "pgs_scores/scores_wgs/bmi.sscore")), "score", "bmi"),
  setnames(rdw(file.path(AD, "pgs_scores/scores_wgs/height.sscore")), "score", "ht"), setnames(rdw(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/g2mh_scores_raw/g2mh_wgs.asd2026_full.sscore")), "score", "asdp")))
pcw <- fread(file.path(AD, "pca_results/naive/cohort.sscore")); setnames(pcw, "#IID", "IID"); pcw <- pcw[, c(list(GUID = sub("^[0-9]+-", "", IID)), .SD), .SDcols = paste0("PC", 1:10, "_AVG")]; setnames(pcw, paste0("PC", 1:10, "_AVG"), paste0("PC", 1:10))
g2 <- merge(g2, pcw, by = "GUID"); rz <- function(v, d) { f <- as.formula(paste("v ~", paste0("PC", 1:10, collapse = " + "))); as.numeric(scale(residuals(lm(f, data = cbind(d, v = v))))) }
g2[, `:=`(pgs_cobes = rz(cobes, g2), pgs_bmi = rz(bmi, g2), pgs_ht = rz(ht, g2), pgs_asd = rz(asdp, g2))]
ss <- fread(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/ssc_pgs_scores.tsv"))
PG <- rbind(g2[, .(key = GUID, pgs_cobes, pgs_bmi, pgs_ht, pgs_asd)], ss[, .(key = paste0("SSC:", individual), pgs_cobes, pgs_bmi, pgs_ht, pgs_asd)])
idmap <- unique(rbind(gmc[, .(id, cohort, sampleID, GUID)], ctl_g[, .(id, cohort, sampleID, GUID)]))
XS <- merge(X, idmap, by = c("id", "cohort"), all.x = TRUE)[, key := fifelse(cohort == "SSC", paste0("SSC:", sampleID), fifelse(!is.na(GUID) & GUID != "", GUID, NA_character_))]
XS <- merge(XS, PG, by = "key")[complete.cases(XS[, .(pgs_cobes, pgs_bmi, pgs_ht, pgs_asd)])]
PGSV <- c(pgs_cobes = "childhood obesity PGS", pgs_bmi = "adult BMI PGS", pgs_ht = "height PGS", pgs_asd = "autism PGS (2026)")
fit_adj <- function(x, covs) tryCatch({
  x <- copy(x)[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))]
  cov <- paste(c(if (nlevels(x$cohort) > 1) "+ cohort", if (nlevels(x$locus) > 1) "+ locus", if (length(covs)) paste("+", covs)), collapse = " ")
  m <- gam(as.formula(paste("y ~ asd + s(age, k = 8) + sex", cov)), data = x, method = "REML"); co <- summary(m)$p.table["asd", ]
  data.table(delta = unname(co[1]), se = unname(co[2]), p = unname(co[4])) }, error = function(e) data.table(delta = NA_real_, se = NA_real_, p = NA_real_))
MODELS <- c(list("age + sex (same individuals)" = character(0)), setNames(as.list(names(PGSV)), paste("+", PGSV)), list("+ all four PGS" = names(PGSV)))
sens <- rbindlist(lapply(PH, function(ph) rbindlist(lapply(names(GROUPS), function(g) {
  x <- XS[pheno == ph & locus %in% GROUPS[[g]]]; if (sum(x$asd == 1L) < 5 || sum(x$asd == 0L) < 5) return(NULL)
  rbindlist(lapply(names(MODELS), function(mn) cbind(data.table(pheno = ph, group = g, model = mn, n_noasd = sum(x$asd == 0L), n_asd = sum(x$asd == 1L)), fit_adj(x, MODELS[[mn]])))) }))))
fwrite(sens, file.path(HERE, "supp_asd_by_cnv_sensitivity.csv"))
cat("\n== sensitivity: ASD effect with PGS covariates (PGS-complete subset)\n"); print(sens[group %in% c("Overgrowth CNVs", "Undergrowth CNVs", "Non-carriers (SSC)", "16p11.2 del") & model %in% c("age + sex (same individuals)", "+ all four PGS"), .(pheno, group, model, n_noasd, n_asd, delta = round(delta, 2), se = round(se, 2), p = signif(p, 2))], nrows = 60)
sn <- function(ph, g, mn, col) sens[pheno == ph & group == g & model == mn][[col]]
sens_txt <- sprintf("Sensitivity (supp_asd_by_cnv_sensitivity.csv; G2MH carriers and SSC non-carriers with polygenic scores, same individuals in every row): adding childhood obesity, adult BMI, height and autism PGS (each, then all four) to the age- and sex-adjusted model leaves the ASD effects unchanged. Head circumference, overgrowth CNVs: Delta %+.2f (p = %s) with age and sex, %+.2f (p = %s) with all four PGS (n = %d / %d); BMI, overgrowth: %+.2f (p = %s) vs %+.2f (p = %s); BMI, undergrowth: %+.2f (p = %s) vs %+.2f (p = %s); BMI, non-carriers: %+.2f (p = %s) vs %+.2f (p = %s). Sex is a covariate in every model.",
  sn("Head circumference", "Overgrowth CNVs", "age + sex (same individuals)", "delta"), pfx(sn("Head circumference", "Overgrowth CNVs", "age + sex (same individuals)", "p")), sn("Head circumference", "Overgrowth CNVs", "+ all four PGS", "delta"), pfx(sn("Head circumference", "Overgrowth CNVs", "+ all four PGS", "p")), sn("Head circumference", "Overgrowth CNVs", "+ all four PGS", "n_noasd"), sn("Head circumference", "Overgrowth CNVs", "+ all four PGS", "n_asd"),
  sn("BMI", "Overgrowth CNVs", "age + sex (same individuals)", "delta"), pfx(sn("BMI", "Overgrowth CNVs", "age + sex (same individuals)", "p")), sn("BMI", "Overgrowth CNVs", "+ all four PGS", "delta"), pfx(sn("BMI", "Overgrowth CNVs", "+ all four PGS", "p")),
  sn("BMI", "Undergrowth CNVs", "age + sex (same individuals)", "delta"), pfx(sn("BMI", "Undergrowth CNVs", "age + sex (same individuals)", "p")), sn("BMI", "Undergrowth CNVs", "+ all four PGS", "delta"), pfx(sn("BMI", "Undergrowth CNVs", "+ all four PGS", "p")),
  sn("BMI", "Non-carriers (SSC)", "age + sex (same individuals)", "delta"), pfx(sn("BMI", "Non-carriers (SSC)", "age + sex (same individuals)", "p")), sn("BMI", "Non-carriers (SSC)", "+ all four PGS", "delta"), pfx(sn("BMI", "Non-carriers (SSC)", "+ all four PGS", "p")))
nn <- function(ph, g, col) res[pheno == ph & group == g][[col]]
leg <- c("Supplementary Figure. Height, head circumference and BMI by autism diagnosis within CNV group.",
  sprintf("Mean +/- SE of (A) height Z (WHO), (B) head-circumference Z (UK90) and (C) BMI Z (WHO child; adult centile vs pooled non-carrier adults, mapped to Z) in carriers with (blue) versus without (grey) an autism diagnosis, for each CNV group, the pooled overgrowth (16p11.2 deletion + 22q11.2 duplication) and undergrowth (16p11.2 duplication + 22q11.2 deletion) groups, and, as the non-carrier comparison, SSC probands with ASD versus SSC designated siblings. Numbers under the bars are n. Printed above each group: the age- and sex-adjusted ASD effect (Delta, SD; GAM y ~ ASD + s(age) + sex, plus cohort and, for pooled groups, locus) and its p-value (bold: p < 0.05). ASD status: G2MH subject table, SVIP diagnosis summary; SSC probands = ASD. Data as in the main figure (ages 2-65, latest visit with the outcome and BMI, plausibility QC applied). Carriers with ASD are mostly children and non-ASD carriers include adult relatives; the adjustment for age addresses this, and 20_asd_age_confound.R shows the head-circumference effect holds within children. Head circumference, pooled overgrowth CNVs: Delta %+.2f, p = %s (n = %d / %d); pooled undergrowth: %+.2f, p = %s; non-carriers: %+.2f, p = %s. BMI, pooled overgrowth: %+.2f, p = %s; undergrowth: %+.2f, p = %s; non-carriers: %+.2f, p = %s. Full statistics in supp_asd_by_cnv_tests.csv.",
          nn("Head circumference", "Overgrowth CNVs", "delta"), pfx(nn("Head circumference", "Overgrowth CNVs", "p")), nn("Head circumference", "Overgrowth CNVs", "n_noasd"), nn("Head circumference", "Overgrowth CNVs", "n_asd"),
          nn("Head circumference", "Undergrowth CNVs", "delta"), pfx(nn("Head circumference", "Undergrowth CNVs", "p")), nn("Head circumference", "Non-carriers (SSC)", "delta"), pfx(nn("Head circumference", "Non-carriers (SSC)", "p")),
          nn("BMI", "Overgrowth CNVs", "delta"), pfx(nn("BMI", "Overgrowth CNVs", "p")), nn("BMI", "Undergrowth CNVs", "delta"), pfx(nn("BMI", "Undergrowth CNVs", "p")), nn("BMI", "Non-carriers (SSC)", "delta"), pfx(nn("BMI", "Non-carriers (SSC)", "p"))))
leg <- c(leg, sens_txt)
writeLines(leg, file.path(HERE, "fig_supp_asd_by_cnv_legend.txt"))
cat("done\n")
