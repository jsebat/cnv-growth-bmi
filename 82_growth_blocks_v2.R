# 82_growth_blocks_v2.R  (JS / Claude 2026-11-20, BMI modifier paper v2: height and head circumference become separate figures)
#
# What it does
#   Sources 38_fig6_numbers.R, which sources 17_main_with_ukbb.R (statistics of both UNCHANGED: per-genotype trajectories, band
#   slopes bsl, all-ages CNV x BMI term est_main in bea, genotype-by-band F tests bgt, five-group omnibus in fig6_numbers.csv),
#   and then ADDS (all new, labelled as new in the outputs):
#   a. Pooled "All four CNVs" BMI slopes per age band, height and head circumference (carriers of all four CNVs vs the same
#      controls), least squares as band_slopes() of 17:
#        primary  "locus-adjusted":  y ~ locus (Controls = reference) + BMI + CNV x BMI + sex + age
#        check    "no locus term":   y ~ CNV * BMI + sex + age  (= band_slopes() on the pooled carriers, the model of
#                                    35f_grant_figure_pooled.R; reproduces the earlier pooled head-circumference numbers)
#      plus a heterogeneity F test of locus x BMI among carriers within the band (loci with >= 15 carriers in the band), and
#      the all-ages age-adjusted pooled CNV x BMI term in the model form of est_main (17, bmi_effect_by_age):
#        "est_main form, no locus term":  y ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI) + CNV x BMI
#        "est_main form + locus" (primary): the same with a locus intercept for carriers
#        "locus-specific age curves":     y ~ group + s(age, by = group) + s(age, by = BMI) + CNV x BMI (the omnibus model of
#                                         38_fig6_numbers.R with its four CNV_k x BMI terms collapsed to one)
#   b. Sex tests (Jon: is there a sex difference in the BMI effect of CNVs on height or head circumference?). Per outcome, per
#      genotype and pooled, all ages and per age band: (i) BMI slope in male and female carriers and controls, (ii) CNV x BMI
#      within each sex, (iii) CNV x BMI x sex from one model with every lower-order term and sex-specific age terms.
#      Bands: plain lm (as band_slopes). All ages: bam, REML, est_main form. Cells with < 15 carriers of either sex (or < 30
#      controls of either sex) are left out.
#   c. Sample sizes per outcome x genotype x cohort x age band (and by sex).
#   d. Figure blocks block_height.png / block_hc.png (20 in wide): top row trajectories a-e (traj2 of 17), second row band
#      panels a'-e' (band_panel of 17) and f = pooled "All four CNVs" band slopes (locus-adjusted model).
#   e. (2026-11-20) family-clustered (sandwich) standard errors of the band-level CNV x BMI terms and of the band-level three-way
#      terms, as a check column (band_slopes_family_clustered.csv; needs the family column of tables built with 00_control_set.R).
#   Sensitivity variants are the same script run in another folder with other inputs (env VISIT_RULE, ANALYSIS_TABLE,
#   DROP_COHORT); 83_growth_sensitivity_summary_v2.R collects them.
#
# RUNS ON EXPANSE (UK Biobank adults are read there; individual data stay on Expanse). Run with the working directory = the
#   job's own run folder, which holds frozen copies of 00_scoring.R, 00_select_visit.R, 17_main_with_ukbb.R,
#   38_fig6_numbers.R, growth_refs_lms.csv, supp_bars_results.csv and the participant table.
# Inputs (env): ANALYSIS_TABLE (participant table .rds), VISIT_RULE (fill_gaps | latest), VARIANT (label written to the
#   outputs), DROP_COHORT (optional: cohort removed from the participant table before 17 runs, e.g. CHOP), MAKE_BLOCKS (1/0).
# Outputs (working directory): everything 17 and 38 write, plus pooled_band_slopes.csv, cnv_x_bmi_all_ages.csv,
#   sex_cnv_bmi_growth.csv, n_by_genotype_cohort.csv, band_slopes_family_clustered.csv, analysis_table_job.rds (participant level, stays), block_height.{png,pdf}, block_height_numbers.csv, block_hc.{png,pdf},
#   block_hc_numbers.csv.
# null model: none; parametric tests. Bands: t test of the CNV x BMI (or CNV x BMI x sex) coefficient and F test of
#   locus x BMI in least-squares models; all ages: t test of the parametric coefficient in a bam (REML) model.

suppressPackageStartupMessages(library(data.table))
VARIANT <- Sys.getenv("VARIANT", "main")
# Job-local participant table (written next to the job, never elsewhere): (i) optional DROP_COHORT sensitivity; (ii) JS 2026-11-20:
# only the five analysis groups are passed to 17 (non-carrier controls and the four CNVs). Tables built with 00_control_set.R keep
# flagged non-controls (SSC probands, SVIP non-carriers with autism, G2MH non-carriers) under their own locus labels; 17 would
# otherwise treat them as extra genotype groups in its genotype-by-band test (bgt). A no-op for the 2026-11-19 table.
dc <- Sys.getenv("DROP_COHORT", "")
t0 <- readRDS(Sys.getenv("ANALYSIS_TABLE")); t1 <- t0
if (nzchar(dc)) { t1 <- t1[cohort != dc]
  cat(sprintf("DROP_COHORT = %s: %d rows / %d people removed; %d rows / %d people left\n", dc, nrow(t0) - nrow(t1), uniqueN(t0$id) - uniqueN(t1$id), nrow(t1), uniqueN(t1$id))) }
G5 <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
if (any(!as.character(t1$locus) %in% G5)) { cat("rows not in the five analysis groups, left out of this run (people):\n"); print(t1[!as.character(locus) %in% G5, .(rows = .N, people = uniqueN(id)), by = locus]) }
t1 <- t1[as.character(locus) %in% G5][, locus := factor(as.character(locus), G5)]
f1 <- file.path(getwd(), "analysis_table_job.rds"); saveRDS(t1, f1); Sys.setenv(ANALYSIS_TABLE = f1)
FAM <- if ("family" %in% names(t1)) unique(t1[, .(id, family)]) else NULL; HAVE_FAM <- !is.null(FAM)
cat("family id for clustered standard errors:", if (HAVE_FAM) sprintf("yes (%d families among %d people)", uniqueN(FAM$family), nrow(FAM)) else "not in this table (clustered columns left NA)", "\n")
rm(t0, t1)
cat("VARIANT =", VARIANT, "| VISIT_RULE =", Sys.getenv("VISIT_RULE", "fill_gaps"), "| ANALYSIS_TABLE =", Sys.getenv("ANALYSIS_TABLE"), "\n")
source("38_fig6_numbers.R")   # sources 17_main_with_ukbb.R; statistics unchanged
OUT2 <- WD
POOL <- "All four CNVs"
addfam <- function(x) { x <- copy(x); x[, fam := if (HAVE_FAM) FAM$family[match(id, FAM$id)] else NA_character_]; if (HAVE_FAM) x[is.na(fam), fam := id]; x }   # UK Biobank: own id
DATA2 <- list(Height = list(car = addfam(H[locus %in% LOCI]), ctl = addfam(Hctl_test)), `Head circumference` = list(car = addfam(dh[locus %in% LOCI]), ctl = addfam(dh[locus == "No CNV"])))
# family-clustered (sandwich, CR1) standard error of one coefficient of an lm; P from t with (clusters - 1) df. Check column only:
# SSC siblings and parents of one family, and Simons VIP relatives, are in the same model.
cl_se <- function(m, cluster, term) { if (!HAVE_FAM || anyNA(cluster)) return(c(NA_real_, NA_real_, NA_real_))
  X <- model.matrix(m); ok <- !is.na(coef(m)); X <- X[, ok, drop = FALSE]; u <- resid(m); A <- chol2inv(chol(crossprod(X))); dimnames(A) <- list(colnames(X), colnames(X))
  S <- rowsum(X * u, cluster); Gn <- nrow(S); N <- nrow(X); K <- ncol(X); V <- A %*% crossprod(S) %*% A * (Gn / (Gn - 1)) * ((N - 1) / (N - K)); dimnames(V) <- dimnames(A)
  se <- sqrt(V[term, term]); c(se, 2 * pt(-abs(coef(m)[[term]] / se), Gn - 1), Gn) }
inband <- function(x, b) x[age >= BANDS[b] & age < BANDS[b + 1]]
co4 <- function(m, term) { s <- summary(m); tb <- if (inherits(m, "gam")) s$p.table else s$coefficients
  if (!term %in% rownames(tb)) return(c(NA_real_, NA_real_, NA_real_)); unname(tb[term, c(1, 2, 4)]) }

# ---------------------------------------------------------------- c. sample sizes
nn <- rbindlist(lapply(names(DATA2), function(o) {
  x <- rbind(DATA2[[o]]$car, DATA2[[o]]$ctl)[, cohort := ifelse(source == "UK Biobank", "UK Biobank", sub(":.*", "", id))]
  x[, age_bin := as.character(cut(age, BANDS, labels = BLAB, right = FALSE))]
  a <- x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female")), by = .(group = locus, cohort, age_bin)]
  b <- x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female")), by = .(group = locus, cohort)][, age_bin := "all ages"]
  cc <- x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female")), by = .(group = locus, age_bin)][, cohort := "all cohorts"]
  d <- x[, .(n = .N, n_male = sum(sex == "male"), n_female = sum(sex == "female")), by = .(group = locus)][, `:=`(cohort = "all cohorts", age_bin = "all ages")]
  cbind(variant = VARIANT, outcome = o, rbind(a, b, cc, d, use.names = TRUE))
}))
nn[, `:=`(group = factor(group, GROUPS), age_bin = factor(age_bin, c("all ages", BLAB)))]; setorder(nn, outcome, group, cohort, age_bin)
fwrite(nn, file.path(OUT2, "n_by_genotype_cohort.csv"))
cat("\n== n per outcome x genotype x cohort (all ages; controls = those used in the tests)\n")
print(dcast(nn[age_bin == "all ages"], outcome + group ~ cohort, value.var = "n", fill = 0))

# ---------------------------------------------------------------- a. pooled band slopes
pooled_band <- function(xc, ctl, adjust) rbindlist(lapply(seq_along(BLAB), function(b) {
  k0 <- inband(ctl, b); c1 <- inband(xc, b); if (nrow(c1) < 15 || nrow(k0) < 30) return(NULL)
  x <- rbind(c1[, .(y, BMI, sex, age, cnv = 1, loc = locus, fam)], k0[, .(y, BMI, sex, age, cnv = 0, loc = "Controls", fam)])
  x[, `:=`(loc = relevel(factor(loc), "Controls"), cnvBMI = cnv * BMI)]
  m <- if (adjust) lm(y ~ loc + BMI + cnvBMI + sex + age, data = x) else lm(y ~ cnv + BMI + cnvBMI + sex + age, data = x)
  B <- coef(m); V <- vcov(m); cf <- summary(m)$coefficients
  # heterogeneity across genotypes among the carriers of this band (loci with >= 15 carriers)
  kl <- c1[, .N, by = locus][N >= 15, locus]; ch <- c1[locus %in% kl][, loc := factor(locus)]
  het <- if (length(kl) >= 2) { a <- anova(lm(y ~ loc + BMI + sex + age, data = ch), lm(y ~ loc * BMI + sex + age, data = ch)); c(a$F[2], a$Df[2], a$`Pr(>F)`[2]) } else rep(NA_real_, 3)
  nl <- sapply(LOCI, function(L) c1[locus == L, .N])
  data.table(age_bin = BLAB[b], model = if (adjust) "locus-adjusted (primary)" else "no locus term (band_slopes of 17, as 35f)",
             n_carriers = nrow(c1), n_controls = nrow(k0), n_16p11.2_del = nl[[1]], n_16p11.2_dup = nl[[2]], n_22q11.2_del = nl[[3]], n_22q11.2_dup = nl[[4]],
             slope_controls = B[["BMI"]], se_controls = sqrt(V["BMI", "BMI"]),
             slope_carriers = B[["BMI"]] + B[["cnvBMI"]], se_carriers = sqrt(V["BMI", "BMI"] + V["cnvBMI", "cnvBMI"] + 2 * V["BMI", "cnvBMI"]),
             diff = B[["cnvBMI"]], se_diff = cf["cnvBMI", 2], p_diff = cf["cnvBMI", 4], p_controls = cf["BMI", 4],
             het_F = het[1], het_df = het[2], het_p = het[3], het_loci = paste(kl, collapse = "; "),
             se_diff_family_clustered = cl_se(m, x$fam, "cnvBMI")[1], p_diff_family_clustered = cl_se(m, x$fam, "cnvBMI")[2], n_families = cl_se(m, x$fam, "cnvBMI")[3])
}))
pooled_all <- function(xc, ctl) {
  x <- rbind(ctl[, .(y, BMI, sex, age, locus)], xc[, .(y, BMI, sex, age, locus)])   # controls first, as in bmi_effect_by_age
  x[, `:=`(is_cnv = factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV")), grp = factor(ifelse(locus == "No CNV", "Controls", locus), c("Controls", LOCI)))]
  x[, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  fits <- list(
    "est_main form, no locus term" = bam(y ~ is_cnv + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI, data = x, method = "REML"),
    "est_main form + locus (primary)" = bam(y ~ grp + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI, data = x, method = "REML"),
    "locus-specific age curves (omnibus model of 38, one CNV x BMI term)" = bam(y ~ grp + s(age, k = 10, by = grp) + s(age, k = 10, by = BMI) + cnvBMI, data = x, method = "REML"))
  nl <- sapply(LOCI, function(L) xc[locus == L, .N])
  rbindlist(lapply(names(fits), function(nm) { v <- co4(fits[[nm]], "cnvBMI")
    data.table(age_bin = "all ages", model = nm, n_carriers = nrow(xc), n_controls = nrow(ctl), n_16p11.2_del = nl[[1]], n_16p11.2_dup = nl[[2]],
               n_22q11.2_del = nl[[3]], n_22q11.2_dup = nl[[4]], diff = v[1], se_diff = v[2], p_diff = v[3]) }))
}
pb <- rbindlist(lapply(names(DATA2), function(o) { d <- DATA2[[o]]
  cbind(variant = VARIANT, outcome = o, group = POOL, rbind(pooled_band(d$car, d$ctl, TRUE), pooled_band(d$car, d$ctl, FALSE), pooled_all(d$car, d$ctl), fill = TRUE)) }))
fwrite(pb, file.path(OUT2, "pooled_band_slopes.csv"))
# per-genotype band models of 17 (band_slopes: y ~ cnv * BMI + sex + age) refitted to add the family-clustered check column
bcl <- rbindlist(lapply(names(DATA2), function(o) rbindlist(lapply(LOCI, function(G) rbindlist(lapply(seq_along(BLAB), function(b) {
  c1 <- inband(DATA2[[o]]$car[locus == G], b); k0 <- inband(DATA2[[o]]$ctl, b); if (nrow(c1) < 15 || nrow(k0) < 30) return(NULL)
  x <- rbind(c1[, .(y, BMI, sex, age, cnv = 1, fam)], k0[, .(y, BMI, sex, age, cnv = 0, fam)])[, cnvBMI := cnv * BMI]
  m <- lm(y ~ cnv + BMI + cnvBMI + sex + age, data = x); cf <- summary(m)$coefficients; v <- cl_se(m, x$fam, "cnvBMI")
  data.table(variant = VARIANT, outcome = o, group = G, age_bin = BLAB[b], n_carriers = nrow(c1), n_controls = nrow(k0), diff = cf["cnvBMI", 1], se_diff = cf["cnvBMI", 2], p_diff = cf["cnvBMI", 4],
             se_diff_family_clustered = v[1], p_diff_family_clustered = v[2], n_families = v[3]) }))))))
bcl <- rbind(bcl, pb[age_bin != "all ages" & model == "locus-adjusted (primary)", .(variant, outcome, group, age_bin, n_carriers, n_controls, diff, se_diff, p_diff, se_diff_family_clustered, p_diff_family_clustered, n_families)])
chk <- merge(bcl[group != POOL], bsl[, .(outcome, group, age_bin, d17 = diff, p17 = p_diff)], by = c("outcome", "group", "age_bin"))
cat("refit of the band models equals band_slopes of 17: max |diff| =", chk[, max(abs(diff - d17))], ", max |p - p| =", chk[, max(abs(p_diff - p17))], "\n")
fwrite(bcl, file.path(OUT2, "band_slopes_family_clustered.csv"))
cat("\n== band-level CNV x BMI: model-based vs family-clustered SE\n"); print(bcl[, .(outcome = substr(outcome, 1, 6), group, age_bin, diff = round(diff, 3), se = round(se_diff, 3), p = signif(p_diff, 2), se_cl = round(se_diff_family_clustered, 3), p_cl = signif(p_diff_family_clustered, 2))], nrows = 100)
cat("\n== pooled 'All four CNVs' BMI slopes (carrier = control + CNV x BMI)\n")
print(pb[, .(outcome = substr(outcome, 1, 6), age_bin, model = substr(model, 1, 28), n_car = n_carriers, n_ctl = n_controls, ctl = round(slope_controls, 3), car = round(slope_carriers, 3),
             diff = round(diff, 3), se = round(se_diff, 3), p = signif(p_diff, 2), het_p = signif(het_p, 2))], nrows = 100)

# all-ages CNV x BMI terms: per genotype (est_main of 17), pooled (new), five-group omnibus and contrasts (38)
em <- unique(bea[group != "No CNV", .(outcome, group, est = est_main, se = se_main, p = p_main)])
em[, n_carriers := mapply(function(o, g) DATA2[[o]]$car[locus == g, .N], outcome, group)][, n_controls := sapply(outcome, function(o) nrow(DATA2[[o]]$ctl))]
aa <- rbind(
  em[, .(quantity = "CNV x BMI, all ages, age-adjusted (est_main of 17)", outcome, group, model = "y ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI) + CNV x BMI", n_carriers, n_controls, est, se, p)],
  pb[age_bin == "all ages", .(quantity = "CNV x BMI, all ages, pooled carriers (new)", outcome, group, model, n_carriers, n_controls, est = diff, se = se_diff, p = p_diff)],
  res[quantity == "omnibus_contrast", .(quantity = "CNV_k x BMI in the five-group model (38, omnibus_contrast)", outcome, group, model = "y ~ group + s(age, by = group) + s(age, by = BMI) + sum CNV_k x BMI", est, se, p)],
  res[quantity == "omnibus_5group", .(quantity = "five-group omnibus: BMI effect differs by genotype (38)", outcome, group = "Controls + four CNVs", model = "joint Wald test of the four CNV_k x BMI terms", chisq, df, p)],
  res[quantity == "heterogeneity", .(quantity = "heterogeneity across the four CNVs, carriers only (38)", outcome, group = "four CNVs", model = "y ~ s(age, by = locus) + locus + BMI + locus:BMI, Wald F of locus:BMI", F, df, p)],
  fill = TRUE)
aa <- cbind(variant = VARIANT, aa); fwrite(aa, file.path(OUT2, "cnv_x_bmi_all_ages.csv"))
cat("\n== all-ages CNV x BMI terms\n"); print(aa[, .(quantity = substr(quantity, 1, 45), outcome = substr(outcome, 1, 6), group, model = substr(model, 1, 24), n_carriers, est = round(est, 3), se = round(se, 3), p = signif(p, 3))], nrows = 100)

# ---------------------------------------------------------------- b. sex tests
mk <- function(xc, ctl) {
  x <- rbind(ctl[, .(y, BMI, sex, age, cnv = 0, loc = "Controls", fam)], xc[, .(y, BMI, sex, age, cnv = 1, loc = locus, fam)])   # controls first
  x[, `:=`(loc = relevel(factor(loc), "Controls"), is_cnv = factor(ifelse(cnv == 1, "CNV", "No CNV"), levels = c("No CNV", "CNV")), male = as.numeric(sex == "male"))]
  x[, `:=`(cnvBMI = cnv * BMI, BMImale = BMI * male, cnvBMImale = cnv * BMI * male, agemale = age * male, cell = interaction(is_cnv, sex, drop = TRUE))]
  x
}
ncell <- function(x) c(cm = x[cnv == 1 & male == 1, .N], cf = x[cnv == 1 & male == 0, .N], km = x[cnv == 0 & male == 1, .N], kf = x[cnv == 0 & male == 0, .N])
sex_band <- function(xc, ctl, pooled, b) {
  x <- mk(inband(xc, b), inband(ctl, b)); n <- ncell(x)
  if (min(n[c("cm", "cf")]) < 15 || min(n[c("km", "kf")]) < 30) return(NULL)
  G <- if (pooled) "loc" else "cnv"
  r <- data.table(scope = BLAB[b], method = "lm within band", n_carriers_male = n[["cm"]], n_carriers_female = n[["cf"]], n_controls_male = n[["km"]], n_controls_female = n[["kf"]])
  for (s in c("male", "female")) {
    xs <- x[sex == s]; if (pooled) xs[, loc := droplevels(loc)]
    m <- lm(as.formula(sprintf("y ~ %s + BMI + cnvBMI + age", G)), data = xs); B <- coef(m); V <- vcov(m); cf <- summary(m)$coefficients
    r[, paste0(c("slope_controls_", "se_controls_", "slope_carriers_", "se_carriers_", "cnv_x_bmi_", "se_cnv_x_bmi_", "p_cnv_x_bmi_"), s) :=
        .(B[["BMI"]], sqrt(V["BMI", "BMI"]), B[["BMI"]] + B[["cnvBMI"]], sqrt(V["BMI", "BMI"] + V["cnvBMI", "cnvBMI"] + 2 * V["BMI", "cnvBMI"]), B[["cnvBMI"]], cf["cnvBMI", 2], cf["cnvBMI", 4])]
  }
  # three-way: every lower-order term (CNV or locus, sex, BMI, their two-way terms) and a sex-specific age term
  m3 <- lm(as.formula(sprintf("y ~ %s * male + BMI + BMImale + cnvBMI + cnvBMImale + age + agemale", G)), data = x); v <- co4(m3, "cnvBMImale")
  vc <- if (is.na(v[1])) rep(NA_real_, 3) else cl_se(m3, x$fam, "cnvBMImale")
  r[, `:=`(cnv_x_bmi_x_sex = v[1], se_cnv_x_bmi_x_sex = v[2], p_cnv_x_bmi_x_sex = v[3], se_cnv_x_bmi_x_sex_family_clustered = vc[1], p_cnv_x_bmi_x_sex_family_clustered = vc[2])]; r
}
slope1 <- function(d, loc_term) tryCatch({   # age-adjusted BMI slope within one group and sex (as m0 of the control branch of bmi_effect_by_age: y ~ BMI + s(age))
  d <- copy(d)[, loc := droplevels(loc)]; kk <- min(15, max(4, floor(uniqueN(d$age) / 2)))
  f <- if (loc_term && nlevels(d$loc) > 1) y ~ BMI + loc + s(age, k = kk) else y ~ BMI + s(age, k = kk)
  co4(bam(f, data = d, method = "REML"), "BMI")[1:2] }, error = function(e) c(NA_real_, NA_real_))
sex_all <- function(xc, ctl, pooled) {
  x <- mk(xc, ctl); n <- ncell(x)
  if (min(n[c("cm", "cf")]) < 15 || min(n[c("km", "kf")]) < 30) return(NULL)
  G <- if (pooled) "loc" else "is_cnv"
  r <- data.table(scope = "all ages", method = "bam REML, est_main form", n_carriers_male = n[["cm"]], n_carriers_female = n[["cf"]], n_controls_male = n[["km"]], n_controls_female = n[["kf"]])
  for (s in c("male", "female")) {
    xs <- x[sex == s]
    m <- bam(as.formula(sprintf("y ~ %s + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI", G)), data = xs, method = "REML"); v <- co4(m, "cnvBMI")
    sk <- slope1(xs[cnv == 0], FALSE); sc <- slope1(xs[cnv == 1], pooled)
    r[, paste0(c("slope_controls_", "se_controls_", "slope_carriers_", "se_carriers_", "cnv_x_bmi_", "se_cnv_x_bmi_", "p_cnv_x_bmi_"), s) := .(sk[1], sk[2], sc[1], sc[2], v[1], v[2], v[3])]
  }
  m3 <- bam(as.formula(sprintf("y ~ %s * male + s(age, k = 15, by = cell) + s(age, k = 10, by = BMI) + s(age, k = 10, by = BMImale) + cnvBMI + cnvBMImale", G)), data = x, method = "REML")
  v <- co4(m3, "cnvBMImale"); r[, `:=`(cnv_x_bmi_x_sex = v[1], se_cnv_x_bmi_x_sex = v[2], p_cnv_x_bmi_x_sex = v[3])]; r
}
sx <- rbindlist(lapply(names(DATA2), function(o) rbindlist(lapply(c(LOCI, POOL), function(G) {
  d <- DATA2[[o]]; xc <- if (G == POOL) d$car else d$car[locus == G]; cat("sex tests:", o, G, "\n")
  r <- rbind(sex_all(xc, d$ctl, G == POOL), rbindlist(lapply(seq_along(BLAB), function(b) sex_band(xc, d$ctl, G == POOL, b)), fill = TRUE), fill = TRUE)
  if (is.null(r) || !nrow(r)) return(NULL); cbind(variant = VARIANT, outcome = o, group = G, r)
}), fill = TRUE)), fill = TRUE)
fwrite(sx, file.path(OUT2, "sex_cnv_bmi_growth.csv"))
cat("\n== sex tests (CNV x BMI within sex; three-way = male minus female)\n")
print(sx[, .(outcome = substr(outcome, 1, 6), group, scope, n_m = n_carriers_male, n_f = n_carriers_female, car_m = round(slope_carriers_male, 2), car_f = round(slope_carriers_female, 2),
             ctl_m = round(slope_controls_male, 2), ctl_f = round(slope_controls_female, 2), cxb_m = round(cnv_x_bmi_male, 2), p_m = signif(p_cnv_x_bmi_male, 2),
             cxb_f = round(cnv_x_bmi_female, 2), p_f = signif(p_cnv_x_bmi_female, 2), three = round(cnv_x_bmi_x_sex, 2), p3 = signif(p_cnv_x_bmi_x_sex, 2))], nrows = 200)

# ---------------------------------------------------------------- d. figure blocks
if (Sys.getenv("MAKE_BLOCKS", "1") == "1") {
  TAG <- theme(plot.tag = element_text(size = 22, face = "bold"))
  mk_block <- function(o, key, ylab_eff, YL, LV) {
    pbo <- pb[outcome == o & model == "locus-adjusted (primary)"]
    yl <- range(c(YL, yl_band(pbo)))
    # layout-only adjustments for the 20 in block (panels 4 in and 3.3 in wide; statistics and panel functions of 17 untouched):
    # smaller subtitles, x axis padded so the outer P labels are not clipped, genotype-test P labels of the control panel staggered
    SUBT <- theme(plot.subtitle = element_text(size = 12.5)); SUBB <- theme(plot.subtitle = element_text(size = 10))
    padx <- function(p) { for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomText")) p$layers[[i]]$aes_params$size <- 2.35   # P labels 2.7 -> 2.35 (narrower panels)
      p + scale_x_discrete(drop = FALSE, labels = p$scales$get_scales("x")$labels, expand = expansion(add = 0.62)) + theme(axis.text.x = element_text(size = 7)) }
    tr <- lapply(seq_along(GROUPS), function(i) v2[[paste(key, GROUPS[i])]] + labs(tag = letters[i]) + TAG + SUBT)
    ef <- lapply(seq_along(GROUPS), function(i) { G <- GROUPS[i]; b <- bsl[outcome == o & group == G]
      p <- suppressMessages(padx(band_panel(b, ylab_eff, G == "No CNV", yl, LV)))
      if (G == "No CNV") {
        pg <- merge(bgt[outcome == o, .(age_bin, p)], b[, .(age_bin)], by = "age_bin")[, age_bin := factor(age_bin, LV)][order(age_bin)]
        pg[, `:=`(y = yl[2] - (0.02 + 0.07 * (seq_len(.N) %% 2 == 0)) * diff(yl), lab = gsub("= < ", "< ", sprintf("p = %s", sapply(p, pfx)), fixed = TRUE), face = ifelse(p < 0.05, "bold", "plain"))]
        p <- p + geom_text(data = pg, aes(age_bin, y, label = lab, fontface = face), inherit.aes = FALSE, vjust = 1, size = 2.5) + theme(axis.text.x = element_text(size = 6.6))
      }
      p + labs(tag = paste0(letters[i], "'"), subtitle = if (G == "No CNV") "Non-carrier controls\nBMI effect by age band (95% CI)\np: BMI effect differs by genotype" else paste0(G, "\nBMI effect by age band (95% CI)")) + TAG + SUBB })
    pf <- suppressMessages(padx(band_panel(pbo, ylab_eff, FALSE, yl, LV))) + labs(tag = "f", subtitle = paste0(POOL, ", adjusted for CNV\nBMI effect by age band (95% CI)")) + TAG + SUBB
    fig <- wrap_plots(tr, nrow = 1) / wrap_plots(c(ef, list(pf)), nrow = 1) + plot_layout(heights = c(1.12, 1))
    f <- file.path(OUT2, paste0("block_", if (o == "Height") "height" else "hc"))
    ggsave(paste0(f, ".png"), fig, width = 20, height = 9.8, dpi = 200, bg = "white")
    tryCatch(ggsave(paste0(f, ".pdf"), fig, width = 20, height = 9.8, device = cairo_pdf), error = function(e) ggsave(paste0(f, ".pdf"), fig, width = 20, height = 9.8))
    # every plotted number and test (the points themselves are participant-level and are not exported)
    trn <- rbindlist(lapply(seq_along(GROUPS), function(i) { G <- GROUPS[i]
      x <- if (G == "No CNV") { if (o == "Height") Hctl_plot else dh[locus == G] } else DATA2[[o]]$car[locus == G]
      x <- copy(x)[, grp := tail_of(BMI)][!is.na(grp)]; e <- em[outcome == o & group == G]
      data.table(panel = letters[i], panel_type = "trajectory by BMI tail", group = G, n_in_panel = nrow(x), n_high_bmi = x[grp == names(COLS)[1], .N], n_low_bmi = x[grp == names(COLS)[2], .N],
                 n_clinical = x[source == "Clinical", .N], n_ukbb = x[source == "UK Biobank", .N],
                 cnv_x_bmi_all_ages = if (nrow(e)) e$est else NA_real_, se = if (nrow(e)) e$se else NA_real_, p = if (nrow(e)) e$p else NA_real_) }))
    efn <- rbindlist(lapply(seq_along(GROUPS), function(i) { G <- GROUPS[i]; b <- bsl[outcome == o & group == G]
      r <- b[, .(panel = paste0(letters[i], "'"), panel_type = "BMI slope by age band", group = G, age_bin, n_carriers, n_controls, slope_controls, se_controls, slope_carriers, se_carriers, diff, se_diff, p = p_diff, p_controls)]
      if (G == "No CNV") r <- merge(r, bgt[outcome == o, .(age_bin, p_differs_by_genotype = p, F_genotype = F, df_genotype = df)], by = "age_bin", all.x = TRUE, sort = FALSE); r }), fill = TRUE)
    pfn <- pbo[, .(panel = "f", panel_type = "BMI slope by age band, pooled, locus-adjusted", group = POOL, age_bin, n_carriers, n_controls, slope_controls, se_controls, slope_carriers, se_carriers, diff, se_diff, p = p_diff, p_controls,
                   het_F, het_df, het_p)]
    fwrite(cbind(variant = VARIANT, outcome = o, rbind(trn, efn, pfn, fill = TRUE)), paste0(f, "_numbers.csv"))
    cat("block written:", f, "\n")
  }
  mk_block("Height", "HEIGHT", "Height per +1 SD BMI (SD)", YLH, LVH)
  mk_block("Head circumference", "HC", "Head circ. per +1 SD BMI (SD)", YLC, LVC)
}
cat("82_growth_blocks_v2 done:", VARIANT, "\n")
