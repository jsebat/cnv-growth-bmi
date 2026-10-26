# Sensitivity analysis for the main figure (JS / Claude 2026-10-26): are the BMI-driven effects (CNV x BMI and
# CNV x age x BMI on height and head circumference vs unaffected controls; BMI slopes on SRS / SCQ vs controls; BMI slopes
# on Penn CNB) driven by sex or by polygenic background? The main-figure tests are rerun within strata of
#   1. sex (male / female; carriers and controls both stratified),
#   2. childhood obesity PGS (Bradfield 2019 trans-ancestral; SBayesRC weights from Anjali), below / above median,
#   3. adult BMI PGS (Yengo 2018; same pipeline), below / above median,
# with carriers and SSC controls both stratified (SSC genotypes scored with the same weights: pgs_obesity_ssc on Expanse;
# G2MH WGS carriers from Anjali's pgs_pca.csv). PGS are residualised on 10 ancestry PCs within cohort and z-scored within
# cohort; the median split is within cohort. SVIP / Searchlight / Cardiff carriers have no scores in this pipeline and are
# excluded from PGS strata (they stay in the sex strata). Within each stratum: within-group BMI slope and age x BMI,
# CNV x BMI shift and CNV x age x BMI vs the controls of the same stratum (models of 23/24/28). Heterogeneity between the two
# strata within carriers: BMI x stratum (age-adjusted; locus covariate when pooled). For PGS, also the continuous
# BMI x PGS interaction within carriers and the PGS main effect.
# Inputs: 23 (A, gmc, ctl_g), 24 (X), Anjali pgs_pca.csv, resources/ssc_pgs_scores.tsv (built by 30a_build_ssc_pgs.R),
#         all_cohorts_adsi_pgs.tsv (SSC PCs), SSC IDs Sanders.csv (SS id -> SSC individual).
# Outputs: supp_table_S4_sex_pgs_strata.csv, sheet "S4 sex and PGS strata" in supp_tables_bmi_effects.xlsx, md section,
#          pgs_strata_sample.csv (who is in which stratum, aggregate counts only in the md).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))
source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages({ library(mgcv); library(openxlsx) })
t0 <- Sys.time()
AD <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data")

# ---------------------------------------------------------------- PGS tables
# G2MH (WGS): raw plink2 scores (Anjali, scores_wgs; IID = "<batch>-<GUID>") residualised on the 10 projected ancestry PCs of
# the WGS cohort (pca_results/naive/cohort.sscore) and z-scored, i.e. the same treatment as the SSC scores in 30a. (Anjali's
# pgs_pca.csv holds the same scores already residualised but covers 1500 of the 1762 scored samples.)
rdw <- function(nm) { x <- fread(file.path(AD, "pgs_scores/scores_wgs", paste0(nm, ".sscore"))); setnames(x, "#IID", "IID"); x[, .(GUID = sub("^[0-9]+-", "", IID), score = SCORE1_AVG)] }
# height: Anjali's GIANT 2022 score; autism: April 2026 ASD GWAS, full weights, scored on the G2MH WGS cohort on Expanse
# (pgs_obesity_ssc/scores/g2mh_wgs.asd2026_full.sscore; G2MH is not a GWAS cohort), copied to Jon/Claude data folder.
rda <- function(p) { x <- fread(p); setnames(x, "#IID", "IID"); x[, .(GUID = sub("^[0-9]+-", "", IID), score = SCORE1_AVG)] }
g2 <- Reduce(function(a, b) merge(a, b, by = "GUID"), list(setnames(rdw("childhood_obesity"), "score", "cobes"), setnames(rdw("bmi"), "score", "bmi"), setnames(rdw("bmi_childhood"), "score", "cbmi"),
                                                          setnames(rdw("height"), "score", "ht"), setnames(rda(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/g2mh_scores_raw/g2mh_wgs.asd2026_full.sscore")), "score", "asd")))
pcw <- fread(file.path(AD, "pca_results/naive/cohort.sscore")); setnames(pcw, "#IID", "IID"); pcw <- pcw[, c(list(GUID = sub("^[0-9]+-", "", IID)), .SD), .SDcols = paste0("PC", 1:10, "_AVG")]; setnames(pcw, paste0("PC", 1:10, "_AVG"), paste0("PC", 1:10))
g2 <- merge(g2, pcw, by = "GUID")
rz <- function(v, d) { f <- as.formula(paste("v ~", paste0("PC", 1:10, collapse = " + "))); as.numeric(scale(residuals(lm(f, data = cbind(d, v = v))))) }
g2[, `:=`(pgs_cobes = rz(cobes, g2), pgs_bmi = rz(bmi, g2), pgs_cbmi = rz(cbmi, g2), pgs_ht = rz(ht, g2), pgs_asd = rz(asd, g2))]
g2 <- g2[, .(GUID, pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)]
# SSC: new scores (Expanse pgs_obesity_ssc) residualised on the ADSI PCs, z within SSC, mapped SS id -> SSC individual id
sscf <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/ssc_pgs_scores.tsv")   # built by 30a from the Expanse scores (participant-level, kept in Jon/Claude)
if (!file.exists(sscf)) stop("SSC PGS scores not found: run 30a_build_ssc_pgs.R after the Expanse job pgs_obesity_ssc finishes")
ss <- fread(sscf)   # columns: individual, pgs_cobes, pgs_bmi, pgs_cbmi (already residualised and z-scored by 30a)
cat("PGS: G2MH", nrow(g2), " SSC", nrow(ss), "\n")

# ---------------------------------------------------------------- analysis tables with sex and PGS
idmap <- unique(rbind(gmc[, .(id, cohort, sampleID, GUID)], ctl_g[, .(id, cohort, sampleID, GUID)]))
G <- A[pheno %in% c("Height", "Head circumference"), .(id, locus, sex, age, BMI, y, pheno, asd)]
G <- merge(G, idmap, by = "id", all.x = TRUE)
G[, key := fifelse(cohort == "SSC", paste0("SSC:", sampleID), fifelse(!is.na(GUID) & GUID != "", GUID, NA_character_))]
C <- X[, .(id, locus = fifelse(grp == "Unaffected relatives", "Unaffected controls", grp), sex, age, BMI, y, pheno, asd, cohort)]
C[, key := fifelse(cohort == "SSC", paste0("SSC:", id), id)]
ALL <- rbind(G[, .(key, locus, sex, age, BMI, y, pheno, asd, cohort)], C[, .(key, locus, sex, age, BMI, y, pheno, asd, cohort)])
ALL <- ALL[locus != "SSC ASD probands"]
P <- rbind(g2[, .(key = GUID, pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)], ss[, .(key = paste0("SSC:", individual), pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)])
ALL <- merge(ALL, P, by = "key", all.x = TRUE)
ALL[, carrier := locus != "Unaffected controls"]
cat("rows with PGS by trait x carrier:\n"); print(dcast(ALL[, .(n = .N, with_pgs = sum(!is.na(pgs_cobes))), by = .(pheno, carrier)], pheno ~ carrier, value.var = c("n", "with_pgs")))
# strata: sex; PGS median split within cohort group (carriers = G2MH; controls = SSC)
ALL[, sex_stratum := fifelse(sex %in% c("male", "female"), sex, NA_character_)]
for (v in c("pgs_cobes", "pgs_bmi", "pgs_ht", "pgs_asd")) {   # median within carriers (G2MH) and within controls (SSC), computed over unique individuals
  mc <- median(unique(ALL[carrier == TRUE & !is.na(get(v)), .(key, p = get(v))])$p); mk <- median(unique(ALL[carrier == FALSE & !is.na(get(v)), .(key, p = get(v))])$p)
  thr <- fifelse(ALL$carrier, mc, mk)
  ALL[, (paste0(v, "_stratum")) := fifelse(is.na(get(v)), NA_character_, fifelse(get(v) > thr, "above median", "below median"))]
  cat(sprintf("%s: median carriers %.3f, controls %.3f\n", v, mc, mk))
}
STRATA <- list(sex = list(var = "sex_stratum", levels = c("male", "female"), label = "Sex"),
               childhood_obesity_PGS = list(var = "pgs_cobes_stratum", levels = c("below median", "above median"), label = "Childhood obesity PGS", cont = "pgs_cobes"),
               adult_BMI_PGS = list(var = "pgs_bmi_stratum", levels = c("below median", "above median"), label = "Adult BMI PGS", cont = "pgs_bmi"),
               height_PGS = list(var = "pgs_ht_stratum", levels = c("below median", "above median"), label = "Height PGS", cont = "pgs_ht"),
               autism_PGS = list(var = "pgs_asd_stratum", levels = c("below median", "above median"), label = "Autism PGS (2026)", cont = "pgs_asd"))
fwrite(ALL[carrier == TRUE, .N, by = .(pheno, locus, sex_stratum, pgs_cobes_stratum, pgs_bmi_stratum, pgs_ht_stratum, pgs_asd_stratum)], file.path(HERE, "pgs_strata_sample.csv"))

# ---------------------------------------------------------------- models (as in 28)
within <- function(x) tryCatch({
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)]
  m1 <- gam(as.formula(paste("y ~ BMI + s(age, k = 8)", cov)), data = x, method = "REML"); co <- summary(m1)$p.table["BMI", ]
  m2 <- gam(as.formula(paste("y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6))", cov)), data = x, method = "REML")
  data.table(n = nrow(x), bmi_slope = unname(co[1]), bmi_se = unname(co[2]), bmi_p = unname(co[4]), age_x_bmi_p = summary(m2)$s.table["ti(age,BMI)", "p-value"]) },
  error = function(e) data.table(n = nrow(x), bmi_slope = NA_real_, bmi_se = NA_real_, bmi_p = NA_real_, age_x_bmi_p = NA_real_))
contrast2 <- function(ref, xc) tryCatch({
  cov <- if (uniqueN(xc$locus) > 1) "locus + " else ""
  z <- rbind(ref[, .(key, locus = "ref", age, BMI, y)], xc[, .(key, locus, age, BMI, y)])[, is_cnv := factor(ifelse(locus == "ref", "ref", "CNV"), levels = c("ref", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  z[, locus := factor(ifelse(locus == "ref", levels(factor(xc$locus))[1], locus))]
  m1 <- bam(as.formula(paste("y ~", cov, "is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI")), data = z, method = "fREML"); co <- summary(m1)$p.table
  m2 <- bam(as.formula(paste("y ~", cov, "is_cnv + s(age, k = 10) + s(BMI, k = 6) + s(age, k = 10, by = is_cnv) + ti(age, BMI, k = c(6, 6)) + ti(age, BMI, k = c(6, 6), by = is_cnv)")), data = z, method = "fREML"); st <- summary(m2)$s.table
  data.table(n_ref = nrow(ref), cnv_x_bmi = co["cnvBMI", 1], cnv_x_bmi_se = co["cnvBMI", 2], cnv_x_bmi_p = co["cnvBMI", 4], cnv_x_age_x_bmi_p = st[grep("^ti\\(age,BMI\\):is_cnv", rownames(st))[1], "p-value"]) },
  error = function(e) data.table(n_ref = nrow(ref), cnv_x_bmi = NA_real_, cnv_x_bmi_se = NA_real_, cnv_x_bmi_p = NA_real_, cnv_x_age_x_bmi_p = NA_real_))
het_strata <- function(x, var) tryCatch({   # BMI x stratum within carriers
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)][, s := factor(get(var))]
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI * s", cov)), data = x, method = "REML"); co <- summary(m)$p.table
  r <- grep("^BMI:s", rownames(co)); data.table(bmi_x_stratum = co[r[1], 1], bmi_x_stratum_se = co[r[1], 2], bmi_x_stratum_p = co[r[1], 4]) },
  error = function(e) data.table(bmi_x_stratum = NA_real_, bmi_x_stratum_se = NA_real_, bmi_x_stratum_p = NA_real_))
cont_pgs <- function(x, v) tryCatch({        # continuous PGS main effect and BMI x PGS within carriers
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)][, pgs := get(v)]
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI * pgs", cov)), data = x, method = "REML"); co <- summary(m)$p.table
  data.table(pgs_main = co["pgs", 1], pgs_main_p = co["pgs", 4], bmi_x_pgs = co["BMI:pgs", 1], bmi_x_pgs_se = co["BMI:pgs", 2], bmi_x_pgs_p = co["BMI:pgs", 4]) },
  error = function(e) data.table(pgs_main = NA_real_, pgs_main_p = NA_real_, bmi_x_pgs = NA_real_, bmi_x_pgs_se = NA_real_, bmi_x_pgs_p = NA_real_))

GROUPS <- c(setNames(as.list(LOCI), LOCI), POOLS)
PH <- c("Height", "Head circumference", "SRS", "SCQ", "CNB accuracy", "CNB speed")
rows <- list()
for (sn in names(STRATA)) { S <- STRATA[[sn]]
  for (ph in PH) { d <- ALL[pheno == ph & !is.na(get(S$var))]
    for (g in names(GROUPS)) {
      xc_all <- d[locus %in% GROUPS[[g]]]; if (nrow(xc_all) < 20) next
      hs <- het_strata(xc_all, S$var); cp <- if (!is.null(S$cont)) cont_pgs(xc_all[!is.na(get(S$cont))], S$cont) else NULL
      for (lv in S$levels) {
        xc <- xc_all[get(S$var) == lv]; if (nrow(xc) < 15) next
        ref <- d[locus == "Unaffected controls" & get(S$var) == lv]
        w <- within(xc); ct <- if (nrow(ref) >= 30 && !ph %like% "CNB") contrast2(ref, xc) else data.table(n_ref = nrow(ref))
        rows[[length(rows) + 1]] <- cbind(data.table(stratification = S$label, stratum = lv, pheno = ph, group = g, members = paste(GROUPS[[g]], collapse = " + ")), w, ct, hs, if (!is.null(cp)) cp else NULL) } } }
  cat(sn, "done", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "min\n") }
T4 <- rbindlist(rows, fill = TRUE)
T4[pheno %like% "CNB", c("n_ref", "cnv_x_bmi", "cnv_x_bmi_se", "cnv_x_bmi_p", "cnv_x_age_x_bmi_p") := NA]
T4[, pheno := factor(pheno, PH)][, group := factor(group, names(GROUPS))]; setorder(T4, stratification, pheno, group, stratum)
fwrite(T4, file.path(HERE, "supp_table_S4_sex_pgs_strata.csv"))
F4 <- T4[, .(Stratification = stratification, Stratum = stratum, Trait = as.character(pheno), Group = as.character(group), n = n, `n controls` = ifelse(is.na(n_ref), "", as.character(n_ref)),
             `BMI slope (SE)` = ifelse(is.na(bmi_slope), "", sprintf("%+.3f (%.3f)", bmi_slope, bmi_se)), `BMI slope p` = pfx(bmi_p), `age x BMI p` = pfx(age_x_bmi_p),
             `CNV x BMI vs controls (SE)` = ifelse(is.na(cnv_x_bmi), "", sprintf("%+.3f (%.3f)", cnv_x_bmi, cnv_x_bmi_se)), `CNV x BMI p` = ifelse(is.na(cnv_x_bmi_p), "", pfx(cnv_x_bmi_p)), `CNV x age x BMI p` = ifelse(is.na(cnv_x_age_x_bmi_p), "", pfx(cnv_x_age_x_bmi_p)),
             `BMI x stratum (carriers) p` = pfx(bmi_x_stratum_p), `PGS main effect p` = ifelse(is.na(pgs_main_p), "", pfx(pgs_main_p)), `BMI x PGS (continuous) p` = ifelse(is.na(bmi_x_pgs_p), "", pfx(bmi_x_pgs_p)))]
# ---------------------------------------------------------------- Table S5: the BMI effects with PGS as covariates
# Does adjusting for polygenic background (each PGS and its interaction with BMI; then all four jointly) change the effect
# of BMI: the carrier BMI slope and the CNV x BMI shift versus unaffected controls? Restricted to individuals with all
# PGS (G2MH carriers, SSC controls) so that the unadjusted and adjusted rows use the same people.
PGSV <- c(pgs_cobes = "childhood obesity PGS", pgs_bmi = "adult BMI PGS", pgs_ht = "height PGS", pgs_asd = "autism PGS (2026)")
slope_adj <- function(x, covs) tryCatch({
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)]
  ex <- if (length(covs)) paste(sprintf("+ %s + %s:BMI", covs, covs), collapse = " ") else ""
  m <- gam(as.formula(paste("y ~ BMI + s(age, k = 8)", cov, ex)), data = x, method = "REML"); co <- summary(m)$p.table
  data.table(n = nrow(x), bmi_slope = co["BMI", 1], bmi_se = co["BMI", 2], bmi_p = co["BMI", 4],
             pgs_x_bmi_p = if (length(covs) == 1) co[paste0("BMI:", covs), 4] else NA_real_) },
  error = function(e) data.table(n = nrow(x), bmi_slope = NA_real_, bmi_se = NA_real_, bmi_p = NA_real_, pgs_x_bmi_p = NA_real_))
contrast_adj <- function(ref, xc, covs) tryCatch({
  cov <- if (uniqueN(xc$locus) > 1) "locus + " else ""; keep <- c("key", "age", "BMI", "y", covs)
  z <- rbind(ref[, c(keep), with = FALSE][, locus := "ref"], xc[, c(keep, "locus"), with = FALSE])[, is_cnv := factor(ifelse(locus == "ref", "ref", "CNV"), levels = c("ref", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  z[, locus := factor(ifelse(locus == "ref", levels(factor(xc$locus))[1], locus))]
  ex <- if (length(covs)) paste(sprintf("+ %s + %s:BMI", covs, covs), collapse = " ") else ""
  m1 <- bam(as.formula(paste("y ~", cov, "is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI", ex)), data = z, method = "fREML"); co <- summary(m1)$p.table
  data.table(n_ref = nrow(ref), cnv_x_bmi = co["cnvBMI", 1], cnv_x_bmi_se = co["cnvBMI", 2], cnv_x_bmi_p = co["cnvBMI", 4]) },
  error = function(e) data.table(n_ref = nrow(ref), cnv_x_bmi = NA_real_, cnv_x_bmi_se = NA_real_, cnv_x_bmi_p = NA_real_))
MODELS <- c(list("unadjusted (same individuals)" = character(0), "+ sex" = "sexm"), setNames(as.list(names(PGSV)), paste("+", PGSV)), list("+ all four PGS" = names(PGSV), "+ all four PGS + sex" = c(names(PGSV), "sexm")))
ALL[, sexm := as.numeric(sex == "male")]   # sex as a covariate (main effect + sex x BMI), same treatment as a PGS
PGSV <- c(PGSV, sexm = "sex (male)")
rows5 <- list()
for (ph in PH) { d <- ALL[pheno == ph & sex %in% c("male", "female")]; d <- d[complete.cases(d[, .(pgs_cobes, pgs_bmi, pgs_ht, pgs_asd)])]
  for (g in names(GROUPS)) { xc <- d[locus %in% GROUPS[[g]]]; if (nrow(xc) < 20) next; ref <- d[locus == "Unaffected controls"]
    for (mn in names(MODELS)) { covs <- MODELS[[mn]]
      rows5[[length(rows5) + 1]] <- cbind(data.table(pheno = ph, group = g, model = mn, covariate = if (length(covs) == 1) PGSV[[covs]] else if (length(covs) == 0) "" else paste(PGSV[covs], collapse = " + ")), slope_adj(xc, covs),
                                          if (nrow(ref) >= 30 && !ph %like% "CNB") contrast_adj(ref, xc, covs) else data.table(n_ref = nrow(ref))) } }
  cat("S5", ph, "done", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "min\n") }
T5 <- rbindlist(rows5, fill = TRUE); T5[, pheno := factor(pheno, PH)][, group := factor(group, names(GROUPS))][, model := factor(model, names(MODELS))]; setorder(T5, pheno, group, model)
fwrite(T5, file.path(HERE, "supp_table_S5_pgs_adjusted.csv"))
F5 <- T5[, .(Trait = as.character(pheno), Group = as.character(group), Model = as.character(model), `n carriers` = n, `n controls` = ifelse(is.na(n_ref), "", as.character(n_ref)),
             `Carrier BMI slope (SE)` = ifelse(is.na(bmi_slope), "", sprintf("%+.3f (%.3f)", bmi_slope, bmi_se)), `BMI slope p` = pfx(bmi_p), `PGS x BMI p` = ifelse(is.na(pgs_x_bmi_p), "", pfx(pgs_x_bmi_p)),
             `CNV x BMI vs controls (SE)` = ifelse(is.na(cnv_x_bmi), "", sprintf("%+.3f (%.3f)", cnv_x_bmi, cnv_x_bmi_se)), `CNV x BMI p` = ifelse(is.na(cnv_x_bmi_p), "", pfx(cnv_x_bmi_p)))]
# ---------------------------------------------------------------- Table S6: power to detect modifier x BMI interactions
# For each modifier (sex, each PGS as a continuous z-score) the interaction term is "modifier x BMI" in the within-carrier
# model y ~ s(age) + BMI * modifier [+ locus]. Power is computed from the term's standard error (normal approximation) for an
# interaction as large as the carriers' own BMI effect (i.e. the BMI effect vanishing in one sex, or changing by its full
# size per SD of PGS): power = P(|Z| > 1.96) with Z = BMI slope / SE(interaction). Also the smallest interaction detectable
# with 80% power (2.8 x SE). Pooled "All CNVs" and per-locus rows.
MODV <- c(sexm = "sex (male vs female)", pgs_cobes = "childhood obesity PGS", pgs_bmi = "adult BMI PGS", pgs_ht = "height PGS", pgs_asd = "autism PGS (2026)")
inter <- function(x, v) tryCatch({
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)][, mod := get(v)]
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI * mod", cov)), data = x, method = "REML"); co <- summary(m)$p.table
  data.table(n = nrow(x), bmi_slope = co["BMI", 1], inter_est = co["BMI:mod", 1], inter_se = co["BMI:mod", 2], inter_p = co["BMI:mod", 4]) },
  error = function(e) data.table(n = nrow(x), bmi_slope = NA_real_, inter_est = NA_real_, inter_se = NA_real_, inter_p = NA_real_))
rows6 <- list()
for (ph in PH) { d <- ALL[pheno == ph & carrier == TRUE & sex %in% c("male", "female")]
  for (g in names(GROUPS)) { xc <- d[locus %in% GROUPS[[g]]]; if (nrow(xc) < 20) next
    for (v in names(MODV)) { x <- xc[!is.na(get(v))]; if (nrow(x) < 20) next
      r <- inter(x, v); r[, `:=`(pheno = ph, group = g, modifier = MODV[[v]])]
      # reference BMI effect: the main-figure carrier slope for this trait/group (unadjusted, from the sensitivity subset)
      r[, `:=`(ref_slope = bmi_slope, mde80 = 2.8 * inter_se, power_ref = pnorm(abs(bmi_slope) / inter_se - 1.96) + pnorm(-abs(bmi_slope) / inter_se - 1.96))]
      rows6[[length(rows6) + 1]] <- r } } }
T6 <- rbindlist(rows6)[, .(pheno, group, modifier, n, bmi_slope, inter_est, inter_se, inter_p, mde80, power_ref)]
T6[, pheno := factor(pheno, PH)][, group := factor(group, names(GROUPS))]; setorder(T6, pheno, group, modifier)
fwrite(T6, file.path(HERE, "supp_table_S6_interaction_power.csv"))
F6 <- T6[, .(Trait = as.character(pheno), Group = as.character(group), Modifier = modifier, `n carriers` = n, `BMI slope (reference effect)` = sprintf("%+.3f", bmi_slope),
             `Interaction SE` = sprintf("%.3f", inter_se), `Detectable interaction (80% power)` = sprintf("%.2f", mde80), `Power for interaction = BMI slope` = sprintf("%.0f%%", 100 * power_ref),
             `Observed interaction (SE)` = sprintf("%+.3f (%.3f)", inter_est, inter_se), `Interaction p` = pfx(inter_p))]
cat("\n== S6 power (All CNVs)\n"); print(F6[Group == "All CNVs"], nrows = 60)
cat("max power across all rows:", round(100 * max(T6$power_ref, na.rm = TRUE)), "%;  rows with power >= 80%:", sum(T6$power_ref >= 0.8, na.rm = TRUE), "of", nrow(T6), "\n")
note6 <- "Table S6: power to detect modifier x BMI interactions within carriers. For each modifier (sex; each PGS as a continuous z-score) the interaction term BMI x modifier was fitted in GAM y ~ s(age) + BMI * modifier [+ locus]. Power (normal approximation from the term's SE) is given for an interaction as large as the carriers' own BMI effect on that trait, i.e. the BMI effect disappearing in one sex or changing by its full size per SD of PGS; 'detectable interaction' is the smallest effect detectable with 80% power (2.8 x SE). Observed interaction estimates and p-values are listed for completeness."

note5 <- c("Table S5: the effect of BMI with sex and polygenic scores as covariates (the confounding check). For each trait and group, the carrier BMI slope (GAM y ~ BMI + s(age) [+ locus]) and the CNV x BMI shift versus unaffected controls (model as in Tables S1-S2) are refit (i) unadjusted, (ii) adding sex and its interaction with BMI, (iii) adding one PGS and its interaction with BMI (childhood obesity, adult BMI, height, or autism; the latter from the April 2026 ASD GWAS, full weights for G2MH carriers and leave-SSC-out weights for SSC controls), (iv) adding all four PGS with their BMI interactions, and (v) all four PGS plus sex. All rows use the same individuals (G2MH carriers and SSC controls with all four scores). PGS x BMI p: whether the BMI slope depends on that PGS within carriers.")
cat("\n== S5 (all CNVs)\n"); print(F5[Group == "All CNVs"], nrows = 100)

note <- c("Table S4: sensitivity of the main-figure BMI effects to sex and polygenic background. The main-figure tests are repeated within strata: sex (male / female; carriers and unaffected controls both stratified); childhood obesity PGS (Bradfield et al. 2019 trans-ancestral GWAS, SBayesRC weights), adult BMI PGS (Yengo et al. 2018), height PGS (Yengo et al. 2022) and autism PGS (April 2026 ASD GWAS; full weights for G2MH carriers, leave-SSC-out weights for SSC controls), each below / above the median. PGS were computed with the same weights on the G2MH WGS cohort (carriers) and on the SSC genotype set (controls), residualised on 10 ancestry principal components and z-scored within cohort; the median split is within cohort. Carriers from SVIP, Searchlight and Cardiff have no PGS in this pipeline and appear only in the sex strata.",
          "Columns: within-stratum BMI slope (SD of trait per SD of BMI, GAM y ~ BMI + s(age) [+ locus]) and age x BMI p; CNV x BMI and CNV x age x BMI versus the unaffected controls of the same stratum (models as in Tables S1-S2); BMI x stratum p: whether the within-carrier BMI slope differs between the two strata (age-adjusted, locus covariate when pooled); PGS main effect and BMI x PGS (continuous) p: within carriers, GAM y ~ s(age) + BMI * PGS [+ locus]. CNB has no control group with BMI.")
wb <- loadWorkbook(file.path(HERE, "supp_tables_bmi_effects.xlsx"))
if ("S4 sex and PGS strata" %in% names(wb)) removeWorksheet(wb, "S4 sex and PGS strata"); addWorksheet(wb, "S4 sex and PGS strata"); writeData(wb, "S4 sex and PGS strata", F4)
setColWidths(wb, "S4 sex and PGS strata", cols = 1:15, widths = c(20, 14, 18, 16, 7, 10, 16, 11, 11, 22, 11, 16, 20, 16, 20)); addStyle(wb, "S4 sex and PGS strata", createStyle(textDecoration = "bold"), rows = 1, cols = 1:15)
if ("S5 PGS adjusted" %in% names(wb)) removeWorksheet(wb, "S5 PGS adjusted"); addWorksheet(wb, "S5 PGS adjusted"); writeData(wb, "S5 PGS adjusted", F5)
setColWidths(wb, "S5 PGS adjusted", cols = 1:10, widths = c(18, 16, 30, 10, 10, 22, 11, 12, 24, 11)); addStyle(wb, "S5 PGS adjusted", createStyle(textDecoration = "bold"), rows = 1, cols = 1:10)
if ("S6 interaction power" %in% names(wb)) removeWorksheet(wb, "S6 interaction power"); addWorksheet(wb, "S6 interaction power"); writeData(wb, "S6 interaction power", F6)
setColWidths(wb, "S6 interaction power", cols = 1:10, widths = c(18, 16, 24, 10, 16, 12, 18, 18, 22, 12)); addStyle(wb, "S6 interaction power", createStyle(textDecoration = "bold"), rows = 1, cols = 1:10)
nt <- read.xlsx(wb, "notes"); nt <- nt[!grepl("^Table S4|^Columns: within-stratum|^Table S5|^Table S6", nt$Note), , drop = FALSE]; nt <- rbind(nt, data.frame(Note = c(note, note5, note6))); removeWorksheet(wb, "notes"); addWorksheet(wb, "notes"); writeData(wb, "notes", nt); setColWidths(wb, "notes", cols = 1, widths = 160)
saveWorkbook(wb, file.path(HERE, "supp_tables_bmi_effects.xlsx"), overwrite = TRUE)
mdt <- function(x, title) { h <- names(x); c(sprintf("## %s", title), "", paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", paste(rep("---", length(h)), collapse = "|"), "|"), apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "") }
md <- readLines(file.path(HERE, "supp_tables_bmi_effects.md")); cut <- grep("^- Table S4|^## Table S4", md); if (length(cut)) md <- md[seq_len(cut[1] - 1)]
writeLines(c(md, paste0("- ", note), "", mdt(F4, "Table S4. Main-figure BMI effects within strata of sex, childhood obesity PGS, adult BMI PGS, height PGS and autism PGS"), paste0("- ", note5), "", mdt(F5, "Table S5. BMI effects with sex and polygenic scores as covariates"), paste0("- ", note6), "", mdt(F6, "Table S6. Power to detect modifier x BMI interactions")), file.path(HERE, "supp_tables_bmi_effects.md"))
cat("\n== S4 (key rows: all CNVs and 16p11.2 del)\n"); print(F4[Group %in% c("All CNVs", "16p11.2 del")], nrows = 200)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
