# 80_longitudinal_v2.R  (2026-11-20, BMI modifier paper v2, work package A: longitudinal analyses)
#
# What it does: ONE longitudinal ("next visit") analysis, rebuilt inside this code base, of whether BMI Z at a visit predicts
#   the change in (1) head-circumference Z and (2) height Z by the next visit in CNV-carrier children (16p11.2 BP4-BP5 and
#   22q11.2 A-D deletions and duplications). It replaces the grant-figure code (COMET_ASD_ARPAH/analysis A1, A4, A7, A8), whose
#   numbers are compared with the ones produced here (compare_with_grant_version.csv).
#
# Model (per genotype, adjusted for cohort; then all four CNVs with one intercept per cohort within genotype, because cohort and
#   genotype are partly confounded (JS 2026-11-21); then heterogeneity across genotypes):
#   change in Z to the next visit ~ Z at the visit + interval + sex (+ cohort | + cohort within genotype) + ns(age at the visit, 3)
#                                   + BMI Z at the visit + (1 | child)                      [lmerTest, REML, Satterthwaite P]
#   by age window : the same model with  window + window:BMI  in place of BMI (window = age at the first visit of the pair);
#                   head circumference < 2, 2-5, 5-18 y; height 2-8, 8-13, 13-18 y (windows of 05_age_varying.R / 07_hc_trajectories.R)
#   by sex        : the same model with  BMI + BMI:sex  (exploratory)
#   A subset in which every child has one interval is fitted by lm (no random intercept).
# null model: none; Wald tests (Satterthwaite) of the BMI coefficients, and likelihood-ratio tests of ML fits for
#   BMI x genotype (heterogeneity), BMI x window and (Wald) BMI x sex.
#
# Data (participant level, controlled; nothing participant-level is written under HERE):
#   DATA/analysis_table_chop_2026-11-19.rds : all cohorts, all visits, raw height / BMI / head circumference (01_build.R + 00_add_chop.R,
#        after 00_qc_growth.R). Used for height in every cohort and for CHOP head circumference.
#   DATA/hc_visits.rds : Simons VIP (v10 growth_measurements, every visit from infancy) and Cardiff head circumference (07_hc_trajectories.R).
# Scoring: 00_scoring.R for every cohort (height and BMI: WHO 2006 < 5 y, WHO 2007 5-18 y; head circumference: UK90). Z recomputed
#   from the raw values here, so no cohort depends on Z scores supplied upstream. childsds scoring is run as a sensitivity.
# Cleaning: height = clean_visits() of 00_helpers.R; head circumference = the rule of 07_hc_trajectories.R (same rule the grant
#   code applied to CHOP). Height pairs are >= 0.5 y apart (00_helpers.R intervals()), head-circumference pairs >= 0.1 y (07).
#
# Outputs (HERE/run_v2_2026-11-20/longitudinal/): results_*.csv, sample_*.csv, sensitivity.csv, sex_interaction.csv,
#   scoring_agreement_childsds_vs_00_scoring.csv, duplicate_check_chop_cardiff.csv, age_coverage.csv, cleaning_log.csv,
#   controls_repeat_check.csv, compare_with_grant_version.csv.  Participant-level panels for 81_fig_longitudinal_v2.R go to
#   DATA/longitudinal_v2_panels.rds.

suppressPackageStartupMessages({ library(data.table); library(lmerTest); library(splines) })
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressPackageStartupMessages(library(childsds))
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "00_helpers.R"))    # HERE, DATA, LOCI, clean_visits(); it also loads the pre-CHOP table as `tab`, replaced below
source(file.path(HERE, "00_scoring.R"))    # score_height, score_bmi_child, score_hc
OUT <- file.path(HERE, "run_v2_2026-11-20", "longitudinal"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 230); set.seed(20261120)
tab <- as.data.table(readRDS(file.path(DATA, "analysis_table_chop_2026-11-19.rds")))
tab[, `:=`(sex = as.character(sex), locus = as.character(locus))]
ALL <- "All four CNVs"
WIN_H <- c("2-8 y", "8-13 y", "13-18 y"); WIN_C <- c("<2 y", "2-5 y", "5-18 y")
MINKIDS <- 15   # an estimate from fewer than 15 children is reported as NA

# ---------------------------------------------------------------- scorers
z_00 <- function(v, age, sex, item) switch(item, height = score_height(v, age, sex), bmi = score_bmi_child(v, age, sex), hc = score_hc(v, age, sex))
z_cs <- function(v, age, sex, item) {   # childsds, as in the grant code (A4, A8) and 07_hc_trajectories.R
  z <- rep(NA_real_, length(v)); ok <- !is.na(v) & !is.na(age) & sex %in% c("male", "female")
  if (item == "hc") { z[ok] <- sds(v[ok], age = age[ok], sex = sex[ok], male = "male", female = "female", item = "hc", ref = uk1990.ref); return(z) }
  y <- ok & age < 5; o <- ok & age >= 5
  if (any(y)) z[y] <- sds(v[y], age = age[y], sex = sex[y], male = "male", female = "female", item = item, ref = who.ref)
  if (any(o)) z[o] <- sds(v[o], age = age[o], sex = sex[o], male = "male", female = "female", item = item, ref = who2007.ref)
  z[!is.finite(z)] <- NA_real_
  z
}
# scorer = "00" (00_scoring.R everywhere), "cs_chop" (childsds for CHOP rows only), "cs_all" (childsds everywhere)
zz <- function(v, age, sex, item, cohort, scorer) {
  a <- z_00(v, age, sex, item); if (scorer == "00") return(a)
  b <- suppressWarnings(z_cs(v, age, sex, item)); if (scorer == "cs_all") return(b)
  ifelse(cohort == "CHOP", b, a)
}

# ---------------------------------------------------------------- height visits and intervals
build_height <- function(scorer = "00", log = FALSE) {
  d <- tab[locus %in% LOCI & sex %in% c("male", "female") & age >= 2 & age < 18 & !is.na(HEIGHT_raw) & !is.na(BMI_raw),
           .(id, cohort, locus, sex, age, HEIGHT_raw, BMI_raw)]
  d[, `:=`(HEIGHT = zz(HEIGHT_raw, age, sex, "height", cohort, scorer), BMI = zz(BMI_raw, age, sex, "bmi", cohort, scorer))]
  d <- d[!is.na(HEIGHT) & !is.na(BMI) & abs(HEIGHT) <= 5 & abs(BMI) <= 5][order(id, age)]
  n0 <- d[, .(visits_before = .N, children_before = uniqueN(id)), by = cohort]
  CLEAN_LOG <<- list(); d <- clean_visits(d, "HEIGHT"); d <- unique(d, by = c("id", "age"))
  if (log) { cl <- CLEAN_LOG[[1]][, .N, by = .(cohort, why)]; n1 <- d[, .(visits_kept = .N), by = cohort]
    assign("LOG_H", merge(merge(n0, n1, by = "cohort"), dcast(cl, cohort ~ why, value.var = "N", fill = 0), by = "cohort", all.x = TRUE), envir = .GlobalEnv) }
  iv <- d[, if (.N >= 2) .(age0 = head(age, -1), dt = diff(age), Z0 = head(HEIGHT, -1), dZ = diff(HEIGHT), B0 = head(BMI, -1),
                           raw0 = head(HEIGHT_raw, -1), sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id][dt >= 0.5]
  iv[, win := cut(age0, c(2, 8, 13, 18), right = FALSE, labels = WIN_H)]
  list(V = d, I = iv)
}

# ---------------------------------------------------------------- head-circumference visits and intervals
clean_hc <- function(d) {   # rule of 07_hc_trajectories.R (per child, age order, against the last kept visit)
  keep <- logical(nrow(d)); why <- character(nrow(d))
  for (ix in split(seq_len(nrow(d)), d$id)) { last <- ix[1]; keep[last] <- TRUE
    for (i in ix[-1]) { dt <- d$age[i] - d$age[last]; dh <- d$hc[i] - d$hc[last]; dz <- d$hcz[i] - d$hcz[last]
      w <- if (dt < ifelse(d$age[i] < 2, 0.04, 0.25)) "duplicate_visit"
           else if (abs(dh) < 0.05 && dt >= 0.1 && d$age[i] < 3) "carried_forward"
           else if (dh < -1.5 || (d$age[i] > 1 && dh / dt > 15)) "impossible_change"
           else if (abs(dz) > 2.5) "jump_gt_2.5SD" else ""
      if (w == "") { keep[i] <- TRUE; last <- i } else why[i] <- w } }
  list(d = d[keep], log = data.table(cohort = d$cohort, why = why)[why != "", .N, by = .(cohort, why)])
}
hv0 <- as.data.table(readRDS(file.path(DATA, "hc_visits.rds")))
cat("hc_visits.rds: children with more than one visit, by cohort (G2MH and SSC have one visit each and cannot contribute):\n")
print(hv0[, .N, by = .(id, cohort)][, .(children = .N, with_repeat = sum(N > 1)), by = cohort])
# CHOP head-circumference visits are read from the workbook with chop_rows() of 00_add_chop.R (the reader 01_build.R uses), NOT from
# analysis_table_chop: 00_qc_growth.R rule 2 (height outside 60-230 cm -> height and BMI NA) was written for children from 2 y and
# deletes every valid infant length under 60 cm, i.e. the BMI of nearly all CHOP visits before about 4 months. Simons VIP infancy
# visits (07_hc_trajectories.R) never went through that rule, so reading the workbook treats the cohorts alike. The table-based
# version is kept as a sensitivity analysis (chop_source = "table").
source(file.path(HERE, "00_add_chop.R"))
chop_wb <- as.data.table(chop_rows())[, .(id = paste("CHOP", sampleID, sep = ":"), cohort = "CHOP", locus = ifelse(CNV == "22q11.2 A-D DEL", "22q11.2 del", "22q11.2 dup"),
                                          sex = ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA_character_)), age, hc = HEADCIRCUM_raw, bmi = BMI_raw, height = HEIGHT_raw)][!is.na(hc)]
cat("CHOP head-circumference visits 0.1-18 y with a BMI: workbook", chop_wb[age >= 0.1 & age < 18 & !is.na(bmi), .N], "| analysis_table_chop", tab[cohort == "CHOP" & age >= 0.1 & age < 18 & !is.na(HEADCIRCUM_raw) & !is.na(BMI_raw), .N],
    "| workbook visits with length < 60 cm:", chop_wb[age >= 0.1 & age < 18 & !is.na(bmi) & height < 60, .N], "
")
build_hc <- function(scorer = "00", log = FALSE, chop_source = "workbook") {
  a <- hv0[cohort %in% c("SVIP", "Cardiff") & locus %in% LOCI, .(id = as.character(id), cohort, locus, sex, age, hc, bmi)]
  b <- if (chop_source == "workbook") chop_wb[, .(id, cohort, locus, sex, age, hc, bmi)] else
    tab[cohort == "CHOP" & locus %in% LOCI & !is.na(HEADCIRCUM_raw), .(id, cohort, locus, sex, age, hc = HEADCIRCUM_raw, bmi = BMI_raw)]
  d <- rbind(a, b)[sex %in% c("male", "female") & age >= 0.1 & age < 18]
  d[!is.na(hc) & (hc < 30 | hc > 70), hc := NA_real_]; d[!is.na(bmi) & (bmi < 10 | bmi > 70), bmi := NA_real_]; d <- d[!is.na(hc)]
  d[, `:=`(hcz = zz(hc, age, sex, "hc", cohort, scorer), bmiz = zz(bmi, age, sex, "bmi", cohort, scorer))]
  d <- d[!is.na(hcz) & abs(hcz) <= 5][order(id, age)]; d[!is.na(bmiz) & abs(bmiz) > 5, bmiz := NA_real_]
  n0 <- d[, .(visits_before = .N, children_before = uniqueN(id)), by = cohort]
  cl <- clean_hc(d); d <- cl$d
  if (log) assign("LOG_C", merge(merge(n0, d[, .(visits_kept = .N, kept_with_bmi = sum(!is.na(bmiz))), by = cohort], by = "cohort"),
                                 dcast(cl$log, cohort ~ why, value.var = "N", fill = 0), by = "cohort", all.x = TRUE), envir = .GlobalEnv)
  iv <- d[!is.na(bmiz)][order(id, age)][, if (.N >= 2) .(age0 = head(age, -1), dt = diff(age), Z0 = head(hcz, -1), dZ = diff(hcz), B0 = head(bmiz, -1),
                                                         raw0 = head(hc, -1), sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id][dt >= 0.1]
  iv[, win := cut(age0, c(0, 2, 5, 18), right = FALSE, labels = WIN_C)]
  list(V = d, I = iv)
}
H <- build_height(log = TRUE); C <- build_hc(log = TRUE)
IH <- H$I; IC <- C$I
fwrite(rbind(cbind(outcome = "Height", LOG_H), cbind(outcome = "Head circumference", LOG_C), fill = TRUE), file.path(OUT, "cleaning_log.csv"))

# ---------------------------------------------------------------- models
adj <- function(x, v) if (uniqueN(x[[v]]) > 1) paste("+", v) else ""
# Cohort adjustment (JS 2026-11-21, final rule of the paper): cohort and genotype are partly confounded (CHOP has only 22q11.2
# carriers, Simons VIP and Searchlight only 16p11.2). A model of one genotype carries a cohort term. A model that pools genotypes
# carries one intercept per cohort within genotype (cg = cohort x genotype cell), which removes cohort differences inside a genotype
# without using cohort to explain differences between genotypes. COH_MODE: "cg" (primary), "additive" (cohort + genotype, primary
# until 2026-11-20) or "genotype" (genotype only, tried on 2026-11-21); the last two are sensitivity rows. The three agree within 0.003.
COH_MODE <- "cg"
adjc <- function(x) if (uniqueN(x$locus) == 1) adj(x, "cohort") else switch(COH_MODE,
  cg = if (uniqueN(paste(x$cohort, x$locus)) > 1) "+ cg" else "", additive = paste(adj(x, "cohort"), "+ locus"), genotype = "+ locus")
with_mode <- function(mode, expr) { old <- COH_MODE; COH_MODE <<- mode; on.exit(COH_MODE <<- old); force(expr) }
core <- function(x) paste("dZ ~ Z0 + dt", adj(x, "sex"), adjc(x), "+ ns(age0, 3)")
fitm <- function(f, x, REML = TRUE) suppressMessages(suppressWarnings(
  { x <- as.data.frame(x); x$cg <- factor(paste(x$cohort, x$locus)); if (anyDuplicated(x$id)) lmer(as.formula(paste(f, "+ (1 | id)")), x, REML = REML) else lm(as.formula(f), x) }))
crow <- function(m, pat) { co <- coef(summary(m)); r <- grep(pat, rownames(co)); data.table(term = rownames(co)[r], beta = co[r, 1], se = co[r, 2], p = co[r, ncol(co)]) }
lrt <- function(m0, m1) tryCatch({ a <- suppressMessages(anova(m0, m1)); p <- a[[grep("^Pr", names(a))]]; p[2] }, error = function(e) NA_real_)
NAROW <- data.table(beta = NA_real_, se = NA_real_, p = NA_real_)
sub_g <- function(I, g) if (g == ALL) I else I[locus == g]
desc <- function(x) data.table(children = uniqueN(x$id), pairs = nrow(x), cohorts = paste(sort(unique(x$cohort)), collapse = "+"),
                               median_age0 = round(median(x$age0), 2), median_interval = round(median(x$dt), 2))
# overall BMI coefficient
fit_all <- function(x) { if (uniqueN(x$id) < MINKIDS) return(NAROW)
  tryCatch(crow(fitm(paste(core(x), "+ B0"), x), "^B0$")[, .(beta, se, p)], error = function(e) NAROW) }
# window-specific BMI coefficients from one model, and the BMI x window test
fit_win <- function(x) { x <- copy(x)[!is.na(win)][, win := droplevels(win)]
  d <- x[, .(children = uniqueN(id), pairs = .N, cohorts = paste(sort(unique(cohort)), collapse = "+"), median_age0 = round(median(age0), 2), median_interval = round(median(dt), 2)), by = .(scope = as.character(win))]
  if (uniqueN(x$id) < MINKIDS || nlevels(x$win) < 2) return(list(est = cbind(d, NAROW[rep(1, nrow(d))]), p_x = NA_real_))
  r <- tryCatch({ m <- fitm(paste(core(x), "+ win + win:B0"), x); k <- crow(m, ":B0$"); k[, scope := sub("^win", "", sub(":B0$", "", term))]
    p <- lrt(fitm(paste(core(x), "+ win + B0"), x, REML = FALSE), fitm(paste(core(x), "+ win + win:B0"), x, REML = FALSE)); list(k = k, p = p) }, error = function(e) NULL)
  if (is.null(r)) return(list(est = cbind(d, NAROW[rep(1, nrow(d))]), p_x = NA_real_))
  e <- merge(d, r$k[, .(scope, beta, se, p)], by = "scope", all.x = TRUE); e[children < MINKIDS, c("beta", "se", "p") := NA_real_]
  list(est = e, p_x = r$p) }
het <- function(x) { if (uniqueN(x$locus) < 2) return(NA_real_)
  f <- sub(" \\+ locus", "", core(x)); tryCatch(lrt(fitm(paste(f, "+ locus + B0"), x, REML = FALSE), fitm(paste(f, "+ locus * B0"), x, REML = FALSE)), error = function(e) NA_real_) }
# full set: per genotype and pooled, overall and by window; heterogeneity overall and within window; BMI x window per genotype
run_set <- function(I, outcome, wins) {
  gs <- c(intersect(LOCI, unique(I$locus)), ALL)
  est <- rbindlist(lapply(gs, function(g) { x <- sub_g(I, g); w <- fit_win(x)
    rbind(cbind(scope = "All ages", desc(x), fit_all(x), p_bmi_x_window = w$p_x), cbind(w$est, p_bmi_x_window = NA_real_), fill = TRUE)[, genotype := g] }))
  hp <- rbind(data.table(scope = "All ages", p_heterogeneity = het(I)),
              rbindlist(lapply(wins, function(w) data.table(scope = w, p_heterogeneity = if (nrow(I[win == w])) het(I[win == w]) else NA_real_))))
  est <- merge(est, hp, by = "scope", all.x = TRUE); est[genotype != ALL, p_heterogeneity := NA_real_]
  est[, `:=`(outcome = outcome, lo = beta - 1.96 * se, hi = beta + 1.96 * se, scope = factor(scope, c("All ages", wins)), genotype = factor(genotype, c(LOCI, ALL)))]
  setorder(est, genotype, scope)
  est[, .(outcome, genotype, scope, cohorts, children, pairs, median_age0, median_interval, beta, se, lo, hi, p, p_bmi_x_window, p_heterogeneity)]
}
RH <- run_set(IH, "Height", WIN_H); RC <- run_set(IC, "Head circumference", WIN_C)
cat("\n== PRIMARY, head circumference\n"); print(RC, digits = 3); cat("\n== PRIMARY, height\n"); print(RH, digits = 3)
fwrite(RC, file.path(OUT, "results_head_circumference.csv")); fwrite(RH, file.path(OUT, "results_height.csv"))

# stratified-by-window fits (the form used in the grant code A8: a separate model per window, no age spline), for comparison
strat_win <- function(I, outcome) rbindlist(lapply(c(LOCI, ALL), function(g) rbindlist(lapply(levels(I$win), function(w) { x <- sub_g(I, g)[win == w]
  r <- if (uniqueN(x$id) < MINKIDS) NAROW else tryCatch(crow(fitm(paste("dZ ~ B0 + Z0 + dt", adj(x, "sex"), adjc(x)), x), "^B0$")[, .(beta, se, p)], error = function(e) NAROW)
  cbind(data.table(outcome = outcome, genotype = g, scope = w, children = uniqueN(x$id), pairs = nrow(x)), r) }))))
SW <- rbind(strat_win(IH, "Height"), strat_win(IC, "Head circumference"))
fwrite(SW, file.path(OUT, "results_by_window_separate_models.csv"))

# ---------------------------------------------------------------- sample sizes
samp <- function(I, outcome) { a <- I[, .(children = uniqueN(id), pairs = .N, median_age0 = round(median(age0), 2), q25_age0 = round(quantile(age0, .25), 2), q75_age0 = round(quantile(age0, .75), 2),
                                           median_interval = round(median(dt), 2)), by = .(genotype = locus, cohort)]
  b <- I[, .(cohort = "All cohorts", children = uniqueN(id), pairs = .N, median_age0 = round(median(age0), 2), q25_age0 = round(quantile(age0, .25), 2), q75_age0 = round(quantile(age0, .75), 2),
             median_interval = round(median(dt), 2)), by = .(genotype = locus)]
  k <- I[, .(genotype = ALL, cohort = "All cohorts", children = uniqueN(id), pairs = .N, median_age0 = round(median(age0), 2), q25_age0 = round(quantile(age0, .25), 2), q75_age0 = round(quantile(age0, .75), 2),
             median_interval = round(median(dt), 2))]
  cbind(outcome = outcome, rbind(a, b, k, use.names = TRUE))[order(factor(genotype, c(LOCI, ALL)), cohort)] }
SC <- samp(IC, "Head circumference"); SH <- samp(IH, "Height")
fwrite(SC, file.path(OUT, "sample_head_circumference_by_genotype_cohort.csv")); fwrite(SH, file.path(OUT, "sample_height_by_genotype_cohort.csv"))
sampw <- function(I, outcome) cbind(outcome = outcome, I[, .(children = uniqueN(id), pairs = .N, median_age0 = round(median(age0), 2)), by = .(genotype = locus, cohort, window = win)][order(factor(genotype, LOCI), window, cohort)])
fwrite(rbind(sampw(IH, "Height"), sampw(IC, "Head circumference")), file.path(OUT, "sample_by_genotype_cohort_window.csv"))
cat("\n== samples\n"); print(SC); print(SH)
# what ages the estimates cover
cover <- function(I, outcome) rbindlist(lapply(c(LOCI, ALL), function(g) { x <- sub_g(I, g)
  cbind(data.table(outcome = outcome, genotype = g, pairs = nrow(x), median_age0 = median(x$age0), q25 = quantile(x$age0, .25), q75 = quantile(x$age0, .75), p90 = quantile(x$age0, .9), max_age0 = max(x$age0)),
        as.data.table(as.list(round(100 * prop.table(table(x$win)), 1)))) }), fill = TRUE)
COV <- rbind(cover(IC, "Head circumference"), cover(IH, "Height"), fill = TRUE)
fwrite(COV, file.path(OUT, "age_coverage.csv")); cat("\n== age coverage (% of visit pairs by window of the index visit)\n"); print(COV, digits = 3)

# ---------------------------------------------------------------- sex (exploratory): BMI x sex
sex_set <- function(I, outcome, wins) rbindlist(lapply(c("All ages", wins), function(w) rbindlist(lapply(c(LOCI, ALL), function(g) {
  x <- sub_g(if (w == "All ages") I else I[win == w], g)
  o <- data.table(outcome = outcome, scope = w, genotype = g, children_f = x[sex == "female", uniqueN(id)], children_m = x[sex == "male", uniqueN(id)], pairs_f = x[sex == "female", .N], pairs_m = x[sex == "male", .N])
  if (min(o$children_f, o$children_m) < MINKIDS) return(o)
  tryCatch({ f <- sub(" \\+ sex", "", core(x)); a <- crow(fitm(paste(f, "+ sex + B0 + B0:sex"), x), "sexmale:B0|B0:sexmale"); b <- crow(fitm(paste(f, "+ sex + sex:B0"), x), ":B0$")
    cbind(o, beta_f = b[grep("female", term), beta], se_f = b[grep("female", term), se], p_f = b[grep("female", term), p],
          beta_m = b[grep("sexmale", term), beta], se_m = b[grep("sexmale", term), se], p_m = b[grep("sexmale", term), p],
          diff_m_minus_f = a$beta, se_diff = a$se, p_bmi_x_sex = a$p) }, error = function(e) o) }), fill = TRUE)), fill = TRUE)
SX <- rbind(sex_set(IC, "Head circumference", WIN_C), sex_set(IH, "Height", WIN_H), fill = TRUE)
# heterogeneity of the BMI x sex term across genotypes (all ages)
sx_het <- function(I) { f <- sub(" \\+ sex", "", sub(" \\+ locus", "", core(I)))
  lrt(fitm(paste(f, "+ locus * B0 + sex * B0 + sex:locus"), I, REML = FALSE), fitm(paste(f, "+ locus * B0 * sex"), I, REML = FALSE)) }
SX[genotype == ALL & scope == "All ages", p_heterogeneity_bmi_x_sex := c(sx_het(IC), sx_het(IH))]
fwrite(SX, file.path(OUT, "sex_interaction.csv")); cat("\n== BMI x sex (exploratory)\n"); print(SX, digits = 3)

# ---------------------------------------------------------------- candidate duplicates: CHOP EHR vs Cardiff 22q11.2 (no ID map exists)
# Strict (the rule set for this pass): same sex and genotype, a visit pair with ages within 0.02 y and raw head circumference or
# height within 0.2 cm.  Cardiff ages are recorded in whole years, so a true duplicate could fail the 0.02 y test; a looser rule
# (CHOP age within [-0.5, +1) y of the Cardiff whole-year age, same measurement within 0.2 cm) is therefore also counted, by the
# number of distinct Cardiff visits matched, with a negative control: the same rule between children of OPPOSITE sex, who cannot
# be the same child, which gives the number of matches expected by chance.  Sharper version of the loose rule: height AND weight
# (weight = BMI x height^2) both agree at the same visit (0.2 cm and 0.5 kg), which a chance height match rarely satisfies.
ca <- tab[cohort == "Cardiff" & locus %in% c("22q11.2 del", "22q11.2 dup"), .(cid = id, locus, csex = sex, cage = age, ch = HEIGHT_raw, chc = HEADCIRCUM_raw, cw = BMI_raw * (HEIGHT_raw / 100)^2)]
ch <- tab[cohort == "CHOP", .(pid = id, locus, psex = sex, page = age, ph = HEIGHT_raw, phc = HEADCIRCUM_raw, pw = BMI_raw * (HEIGHT_raw / 100)^2)]
mm <- merge(ca, ch, by = "locus", allow.cartesian = TRUE)
mm[, `:=`(dage = page - cage, mh = !is.na(ch) & !is.na(ph) & abs(ch - ph) <= 0.2, mc = !is.na(chc) & !is.na(phc) & abs(chc - phc) <= 0.2)]
mm <- mm[(mh | mc) & dage >= -0.5 & dage < 1]
mm[, `:=`(same_sex = csex == psex, strict = abs(dage) <= 0.02, hw = mh & !is.na(cw) & !is.na(pw) & abs(cw - pw) <= 0.5)]
nvis <- ca[, .(cardiff_visits = .N), by = cid]
pp <- merge(mm[, .(matched_cardiff_visits = uniqueN(cage), strict_visits = uniqueN(cage[strict]), hw_visits = uniqueN(cage[hw]), both_measures = any(mh & mc)), by = .(cid, pid, same_sex)], nvis, by = "cid")
cnt <- function(x, lab) data.table(rule = lab, same_sex_child_pairs = x[same_sex == TRUE, .N], same_sex_cardiff_children = x[same_sex == TRUE, uniqueN(cid)], same_sex_chop_children = x[same_sex == TRUE, uniqueN(pid)],
                                   opposite_sex_child_pairs_negative_control = x[same_sex == FALSE, .N])
DUP <- rbind(cnt(pp[strict_visits >= 1], "STRICT: >= 1 visit, age within 0.02 y, HC or height within 0.2 cm"),
             cnt(pp[strict_visits >= 2], "STRICT: >= 2 visits"),
             cnt(pp[matched_cardiff_visits >= 1], "LOOSE (whole-year Cardiff age): >= 1 Cardiff visit matched"),
             cnt(pp[matched_cardiff_visits >= 2], "LOOSE: >= 2 Cardiff visits matched"),
             cnt(pp[matched_cardiff_visits >= 2 & matched_cardiff_visits == cardiff_visits], "LOOSE: every Cardiff visit of the child matched (>= 2)"),
             cnt(pp[matched_cardiff_visits >= 3], "LOOSE: >= 3 Cardiff visits matched"),
             cnt(pp[hw_visits >= 1], "LOOSE + WEIGHT: height (0.2 cm) and weight (0.5 kg) agree at >= 1 Cardiff visit"),
             cnt(pp[hw_visits >= 2], "LOOSE + WEIGHT: height and weight agree at >= 2 Cardiff visits"))
DUP[, `:=`(cardiff_22q_children = uniqueN(ca$cid), chop_children = uniqueN(ch$pid))]
fwrite(DUP, file.path(OUT, "duplicate_check_chop_cardiff.csv")); cat("\n== candidate duplicates CHOP vs Cardiff\n"); print(DUP)
dup_strict <- pp[same_sex == TRUE & strict_visits >= 1, unique(cid)]            # Cardiff record removed (the shorter series), CHOP kept
dup_loose <- pp[same_sex == TRUE & matched_cardiff_visits >= 2, unique(cid)]    # upper bound
dup_hw <- pp[same_sex == TRUE & hw_visits >= 2, unique(cid)]                     # height and weight agree at >= 2 visits
cat("loose candidates with >= 3 Cardiff visits matched: age offsets (CHOP - Cardiff) and which measure matched
")
print(mm[paste(cid, pid) %in% pp[same_sex == TRUE & (matched_cardiff_visits >= 3 | hw_visits >= 2), paste(cid, pid)], .(cage, dage = round(dage, 2), ch, ph, cw = round(cw, 1), pw = round(pw, 1), chc, phc), by = .(cid, pid)][order(cid, pid, cage)], nrows = 300)

# ---------------------------------------------------------------- childsds vs 00_scoring agreement (weakness a)
agree <- function(d, v, item, what) { d <- copy(d); d[, `:=`(a = z_00(get(v), age, sex, item), b = suppressWarnings(z_cs(get(v), age, sex, item)))]
  d[, `:=`(only00 = sum(!is.na(a) & is.na(b)), onlycs = sum(is.na(a) & !is.na(b))), by = cohort]
  d[!is.na(a) & !is.na(b), .(measure = what, visits = .N, scored_by_00_scoring_only = only00[1], scored_by_childsds_only = onlycs[1], correlation = cor(a, b), mean_diff_childsds_minus_00 = mean(b - a), sd_diff = sd(b - a), max_abs_diff = max(abs(b - a)),
                             n_abs_diff_gt_0.1 = sum(abs(b - a) > 0.1)), by = cohort] }
AG <- rbind(agree(C$V, "hc", "hc", "Head-circumference Z (UK90), cleaned HC visits"), agree(C$V[!is.na(bmi)], "bmi", "bmi", "BMI Z (WHO/WHO 2007) at HC visits"),
            agree(H$V, "HEIGHT_raw", "height", "Height Z (WHO/WHO 2007), cleaned height visits 2-18 y"), agree(H$V, "BMI_raw", "bmi", "BMI Z (WHO/WHO 2007) at height visits"))
fwrite(AG, file.path(OUT, "scoring_agreement_childsds_vs_00_scoring.csv")); cat("\n== childsds vs 00_scoring\n"); print(AG, digits = 4)

# ---------------------------------------------------------------- sensitivity analyses
first_iv <- function(I) I[order(id, age0)][, .SD[1], by = id]
sens_one <- function(I, outcome, wins, label) { if (!nrow(I)) return(NULL)
  r <- tryCatch(run_set(I, outcome, wins), error = function(e) { cat("sensitivity failed:", label, conditionMessage(e), "\n"); NULL }); if (is.null(r)) return(NULL)
  cbind(sensitivity = label, r) }
Hc <- build_height("cs_chop"); Cc <- build_hc("cs_chop"); Ha <- build_height("cs_all"); Ca <- build_hc("cs_all"); Ct <- build_hc(chop_source = "table")
sens_all <- function(I, Ics, Ica, outcome, wins, Itab = NULL) rbindlist(list(
  sens_one(I, outcome, wins, "PRIMARY"),
  with_mode("additive", sens_one(I, outcome, wins, "Pooled models with separate cohort and genotype terms (primary until 2026-11-20; per-genotype rows identical to PRIMARY)")),
  with_mode("genotype", sens_one(I, outcome, wins, "Pooled models with a genotype term and no cohort term (per-genotype rows identical to PRIMARY)")),
  sens_one(first_iv(I), outcome, wins, "One interval per child (the first)"),
  sens_one(I[dt >= 0.5 & dt <= 3], outcome, wins, "Intervals of 0.5-3 y only"),
  sens_one(I[cohort != "CHOP"], outcome, wins, "CHOP excluded"),
  sens_one(I[cohort == "CHOP"], outcome, wins, "CHOP only (22q11.2)"),
  sens_one(I[!id %in% dup_strict], outcome, wins, sprintf("Candidate duplicates removed, strict rule (%d Cardiff children; chance-level matches)", length(dup_strict))),
  sens_one(I[!id %in% dup_loose], outcome, wins, sprintf("Candidate duplicates removed, loose rule (%d Cardiff children; upper bound)", length(dup_loose))),
  sens_one(I[!id %in% dup_hw], outcome, wins, sprintf("Candidate duplicates removed, height + weight rule (%d Cardiff children)", length(dup_hw))),
  sens_one(Ics, outcome, wins, "childsds scoring for CHOP"),
  sens_one(Ica, outcome, wins, "childsds scoring for every cohort"),
  if (!is.null(Itab)) sens_one(Itab, outcome, wins, "CHOP from analysis_table_chop (infant lengths < 60 cm removed by 00_qc_growth)")), fill = TRUE)
SENS <- rbind(sens_all(IC, Cc$I, Ca$I, "Head circumference", WIN_C, Ct$I), sens_all(IH, Hc$I, Ha$I, "Height", WIN_H))
fwrite(SENS, file.path(OUT, "sensitivity.csv"))
cat("\n== SENSITIVITY: head circumference, all ages\n"); print(dcast(SENS[outcome == "Head circumference" & scope == "All ages"], sensitivity ~ genotype, value.var = "beta"), digits = 3)
print(dcast(SENS[outcome == "Head circumference" & scope == "All ages"], sensitivity ~ genotype, value.var = "p"), digits = 2)
cat("\n== SENSITIVITY: height, 2-8 y\n"); print(dcast(SENS[outcome == "Height" & scope == "2-8 y"], sensitivity ~ genotype, value.var = "beta"), digits = 3)
print(dcast(SENS[outcome == "Height" & scope == "2-8 y"], sensitivity ~ genotype, value.var = "p"), digits = 2)
cat("\n== SENSITIVITY: heterogeneity and BMI x window P\n"); print(SENS[genotype == ALL & scope == "All ages", .(outcome, sensitivity, children, pairs, beta, p, p_heterogeneity, p_bmi_x_window)], digits = 3)

# ---------------------------------------------------------------- non-carrier controls: repeated measurements?
CTL <- tab[locus == "No CNV", .(rows = .N, distinct_ages = uniqueN(age)), by = .(id, cohort)][, .(non_carriers = .N, with_more_than_one_row = sum(rows > 1), with_more_than_one_age = sum(distinct_ages > 1)), by = cohort]
fwrite(CTL, file.path(OUT, "controls_repeat_check.csv")); cat("\n== non-carrier controls with repeated measurements\n"); print(CTL)
# carriers of other cohorts with repeated child visits (why G2MH and CHU Sainte-Justine do not enter)
REP <- tab[locus %in% LOCI & age < 18, .(height_bmi_ages = uniqueN(age[!is.na(HEIGHT_raw) & !is.na(BMI_raw)]), hc_bmi_ages = uniqueN(age[!is.na(HEADCIRCUM_raw) & !is.na(BMI_raw)])), by = .(id, cohort)][
  , .(carrier_children = .N, with_2plus_height_bmi_visits = sum(height_bmi_ages >= 2), with_2plus_hc_bmi_visits = sum(hc_bmi_ages >= 2)), by = cohort][order(cohort)]
fwrite(REP, file.path(OUT, "carrier_repeat_visits_by_cohort.csv")); cat("\n== carrier children with repeated visits in analysis_table_chop (SVIP head circumference comes from hc_visits.rds instead)\n"); print(REP)

# ---------------------------------------------------------------- comparison with the grant-figure version (COMET_ASD_ARPAH, 2026-11-19)
ref <- rbind(data.table(outcome = "Head circumference", scope = "All ages", genotype = c(LOCI, ALL), ref_children = c(91, 46, 331, 45, 513), ref_pairs = c(287, 159, 1396, 228, 2070),
                        ref_beta = c(0.0821, 0.0504, 0.0121, 0.0478, 0.0289), ref_p = c(0.0392, 0.3474, 0.3427, 0.2638, 0.0124), ref_source = "Head_growth_by_genotype_numbers.txt (A7)"),
             data.table(outcome = "Height", scope = rep(WIN_H, 4), genotype = rep(LOCI, each = 3), ref_children = c(61, 60, 28, 20, 25, 12, 367, 357, 243, 48, 44, 26), ref_pairs = c(93, 90, 41, 20, 40, 24, 1077, 840, 556, 132, 113, 66),
                        ref_beta = c(0.1411, -0.0060, 0.0207, 0.0859, 0.1682, NA, 0.0461, -0.0034, -0.0083, 0.1609, 0.0027, 0.0181), ref_p = c(0.032, 0.93, 0.70, 0.63, 0.12, NA, 7.8e-5, 0.72, 0.51, 3.9e-4, 0.94, 0.74),
                        ref_source = "BMI_by_age_window_numbers.txt (A8; separate model per window)"))
new <- rbind(RC[scope == "All ages", .(outcome, scope = as.character(scope), genotype = as.character(genotype), children, pairs, beta, p, v2_model = "primary")],
             SW[outcome == "Height", .(outcome, scope, genotype, children, pairs, beta, p, v2_model = "separate model per window (as A8)")],
             RH[scope != "All ages", .(outcome, scope = as.character(scope), genotype = as.character(genotype), children, pairs, beta, p, v2_model = "primary (one model, window x BMI)")])
CMP <- merge(new, ref, by = c("outcome", "scope", "genotype"), all.x = TRUE)[!is.na(ref_source)][order(outcome, v2_model, factor(genotype, c(LOCI, ALL)), factor(scope, c("All ages", WIN_H)))]
CMP <- rbind(CMP, data.table(outcome = "Head circumference", scope = "All ages", genotype = "Heterogeneity P (BMI x genotype)", p = RC[genotype == ALL & scope == "All ages", p_heterogeneity], ref_p = 0.0601, v2_model = "primary", ref_source = "A7"), fill = TRUE)
fwrite(CMP, file.path(OUT, "compare_with_grant_version.csv")); cat("\n== v2 vs grant version\n"); print(CMP, digits = 3)

# ---------------------------------------------------------------- panels for the figure script (participant level -> DATA)
partial <- function(I) rbindlist(lapply(c(LOCI, ALL), function(g) { x <- copy(sub_g(I, g)); m <- fitm(paste(core(x), "+ B0"), x); co <- coef(summary(m))
  fx <- if (inherits(m, "lm")) predict(m) else predict(m, re.form = NA)
  x[, .(panel = g, id, cohort, locus, sex, age0, B0, pr = dZ - fx + co["B0", 1] * B0, slope = co["B0", 1], se = co["B0", 2], p = co["B0", ncol(co)])] }))
saveRDS(list(IC = IC, IH = IH, PRC = partial(IC), RC = RC, RH = RH, built = Sys.time()), file.path(DATA, "longitudinal_v2_panels.rds"))
cat("\nDone. Outputs in", OUT, "\n")
