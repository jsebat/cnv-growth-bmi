# Gene x BMI (two-way) on SRS / SCQ across all Simons Searchlight v14 genotypes (CNVs and single genes) with > 100 carriers
# (JS / Claude 2026-11-14). Searchlight all_cohorts: srs_sa (school-age, child form) / srs_adult (adult form) totals, scq_life
# summary score, genetic_status (the CNV or gene) and age_at_evaluation; height and weight parent-reported in
# medical_history_version_3.1 / 3.2 (cm or ft/in; kg or lb). BMI = weight / height^2 from the medical history nearest the
# questionnaire (within 3 y), z by 00_scoring.R (WHO < 19 y; adult centile among clinical non-carrier adults); implausible
# BMI (< 10 or > 70) dropped. Scale: Box-Cox with the SSC-relative lambdas, z vs SSC relatives (00_boxcox_srs_scq.R).
# Model (as Fig. 5k-t): genotype carriers vs SSC designated siblings (and vs all SSC relatives), y ~ is_cnv + s(age, by = is_cnv) +
# s(age, by = BMI) + cnvBMI; carrier BMI slope; ordinal check. Genotypes with > 100 carriers with trait and BMI are tested;
# 50-100 reported as secondary. Output: searchlight_gene_x_bmi.csv, searchlight_gene_x_bmi_counts.csv (aggregate).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages(library(ordinal))
num <- function(v) suppressWarnings(as.numeric(v))
SL <- file.path(DB, "G2MH/Simons_Searchlight/Simons_Searchlight_Dataset_v14.0/all_cohorts/measure_data")
rd <- function(f) fread(file.path(SL, f), colClasses = "character")
toyr <- function(a) { a <- num(a); if (max(a, na.rm = TRUE) > 120) a / 12 else a }
ph <- rbind(rd("srs_sa.csv")[, .(id = sfari_id, gs = genetic_status, raw = num(total), age_ph = toyr(age_at_evaluation), form2 = "child", pheno = "SRS")],
            rd("srs_adult.csv")[, .(id = sfari_id, gs = genetic_status, raw = num(total), age_ph = toyr(age_at_evaluation), form2 = "adult", pheno = "SRS")],
            rd("scq_life.csv")[, .(id = sfari_id, gs = genetic_status, raw = num(summary_score), age_ph = toyr(age_at_evaluation), form2 = "scq", pheno = "SCQ")])
ph <- ph[!is.na(raw) & !is.na(age_ph) & gs != "" & !grepl("negative|normal", gs, ignore.case = TRUE)][order(id, pheno, form2 == "adult", -age_ph)][, .SD[1], by = .(id, pheno)]
mh <- rbind(rd("medical_history_version_3.1.csv"), rd("medical_history_version_3.2.csv"), fill = TRUE)
mh[, h := fifelse(!is.na(num(height_cm)), num(height_cm), (num(height_feet) * 12 + fifelse(is.na(num(height_inches)), 0, num(height_inches))) * 2.54)]
mh[, w := fifelse(!is.na(num(weight_kg)), num(weight_kg), num(weight_lbs) * 0.4536)]
mh <- mh[, .(id = sfari_id, sex_m = tolower(sex), age_v = toyr(age_at_evaluation), bmi_raw = w / (h / 100)^2)][!is.na(bmi_raw) & bmi_raw >= 10 & bmi_raw <= 70 & !is.na(age_v)]
x <- merge(ph, mh, by = "id", allow.cartesian = TRUE)[, dage := abs(age_v - age_ph)][order(id, pheno, dage)][, .SD[1], by = .(id, pheno)][dage <= 3]
x[, sex := sex_m][, BMI := bmiZ(bmi_raw, age_ph, sex)][, age := age_ph][, male := as.numeric(sex == "male")]
x <- x[!is.na(BMI) & age >= 4 & age <= 65]
for (f in c("child", "adult", "scq")) { p <- if (f == "scq") "SCQ" else "SRS"; r <- D[pheno == p & form2 == f & unaffected == TRUE & locus == "No CNV", raw]; L <- LAM[[f]]
  bc <- function(v) if (abs(L) < 1e-8) log(v + 1) else ((v + 1)^L - 1) / L; x[form2 == f, y := (bc(raw) - mean(bc(r))) / sd(bc(r))] }
cnt <- x[, .N, by = .(gs, pheno)]; cnt <- dcast(cnt, gs ~ pheno, value.var = "N", fill = 0)[order(-SRS)]
cnt_all <- dcast(ph[, .N, by = .(gs, pheno)], gs ~ pheno, value.var = "N", fill = 0); setnames(cnt_all, c("SCQ", "SRS"), c("SCQ_any", "SRS_any"))
cnt <- merge(cnt_all, cnt, by = "gs", all.x = TRUE)[order(-SRS)]; fwrite(cnt, file.path(HERE, "searchlight_gene_x_bmi_counts.csv"))
cat("carriers with trait (any) and with trait + BMI, top 30:\n"); print(head(cnt, 30))
# ---------------------------------------------------------------- models
SIB <- D[cohort == "SSC" & locus == "No CNV" & unaffected == TRUE & grepl("[.]s[0-9]+$", id)]; REL <- D[cohort == "SSC" & locus == "No CNV" & unaffected == TRUE]
vs <- function(xc, ref) tryCatch({ z <- rbind(xc[, .(y, BMI, age, is_cnv = 1)], ref[, .(y, BMI, age, is_cnv = 0)])[, is_c := factor(ifelse(is_cnv == 1, "CNV", "ref"), c("ref", "CNV"))][, cnvBMI := is_cnv * BMI]
  m <- bam(y ~ is_c + s(age, k = 10, by = is_c) + s(age, k = 8, by = BMI) + cnvBMI, data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]; c(b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) c(b = NA, se = NA, p = NA))
sl <- function(xc) tryCatch({ m <- gam(y ~ s(age, k = 6) + male + BMI, data = xc, method = "REML"); co <- summary(m)$p.table["BMI", ]; c(b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) c(b = NA, se = NA, p = NA))
out <- list()
for (p in c("SRS", "SCQ")) { tab <- x[pheno == p, .N, by = gs][N >= 50]
  for (g in tab$gs) { xc <- x[pheno == p & gs == g]; tier <- if (nrow(xc) > 100) "N > 100" else "N 50-100"
    a <- sl(xc); b <- vs(xc, SIB[pheno == p]); cc <- vs(xc, REL[pheno == p])
    out[[length(out) + 1]] <- data.table(trait = p, genotype = g, tier = tier, n = nrow(xc), pct_child_form = round(100 * mean(xc$form2 != "adult")), median_age = median(xc$age),
      carrier_slope = a[["b"]], carrier_slope_p = a[["p"]], x_bmi_vs_siblings = b[["b"]], se_sib = b[["se"]], p_sib = b[["p"]], x_bmi_vs_relatives = cc[["b"]], p_rel = cc[["p"]]) } }
R <- rbindlist(out)[order(trait, -n)]; fwrite(R, file.path(HERE, "searchlight_gene_x_bmi.csv"))
options(width = 230); print(R[, .(trait, genotype, tier, n, age = round(median_age, 1), slope = round(carrier_slope, 3), p_slope = signif(carrier_slope_p, 2), x_sib = round(x_bmi_vs_siblings, 3), se = round(se_sib, 3), p_sib = signif(p_sib, 2), x_rel = round(x_bmi_vs_relatives, 3), p_rel = signif(p_rel, 2))], nrows = 100)
# ---------------------------------------------------------------- pooled genotypes (JS 2026-11-14: no single gene reaches N > 100 with BMI)
# all Searchlight genotypes other than proximal 16p11.2 pooled, with genotype fixed effects; also CNVs and single genes separately
vsg <- function(xc, ref) tryCatch({ z <- rbind(xc[, .(y, BMI, age, g = gs, is_cnv = 1)], ref[, .(y, BMI, age, g = "ref", is_cnv = 0)])[, is_c := factor(ifelse(is_cnv == 1, "CNV", "ref"), c("ref", "CNV"))][, cnvBMI := is_cnv * BMI]
  m <- bam(y ~ is_c + g + s(age, k = 10, by = is_c) + s(age, k = 8, by = BMI) + cnvBMI, data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]; c(b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) c(b = NA, se = NA, p = NA))
CNVS <- grep("^[0-9X]|deletion|duplication|NRXN1|KANSL1", unique(x$gs), value = TRUE)
out2 <- list()
for (p in c("SRS", "SCQ")) for (set in c("all genotypes except proximal 16p11.2", "other CNVs", "single genes", "proximal 16p11.2 (del + dup)")) {
  xc <- x[pheno == p][switch(set, `all genotypes except proximal 16p11.2` = !gs %in% c("16p11.2 deletion", "16p11.2 duplication"), `other CNVs` = gs %in% CNVS & !gs %in% c("16p11.2 deletion", "16p11.2 duplication"),
                       `single genes` = !gs %in% CNVS, `proximal 16p11.2 (del + dup)` = gs %in% c("16p11.2 deletion", "16p11.2 duplication"))]
  b <- vsg(xc, SIB[pheno == p]); m <- tryCatch(summary(gam(y ~ gs + s(age, k = 6) + male + BMI, data = xc, method = "REML"))$p.table["BMI", ], error = function(e) rep(NA, 4))
  out2[[length(out2) + 1]] <- data.table(trait = p, set = set, n = nrow(xc), n_genotypes = uniqueN(xc$gs), carrier_slope = m[1], carrier_slope_p = m[4], x_bmi_vs_siblings = b[["b"]], se = b[["se"]], p = b[["p"]]) }
R2 <- rbindlist(out2); fwrite(R2, file.path(HERE, "searchlight_gene_x_bmi_pooled.csv"))
print(R2[, .(trait, set, n, n_genotypes, slope = round(carrier_slope, 3), p_slope = signif(carrier_slope_p, 2), x_sib = round(x_bmi_vs_siblings, 3), se = round(se, 3), p = signif(p, 2))])
