# site_export.R  (2026-11-20)
# What: federated replication of the non-genetic analyses of the Sacks 2026 "BMI modifier" paper inside ONE site.
#   The site that holds the data runs this script; only the aggregate CSV files written here leave the site.
#   a. sample sizes by genotype x age band x sex (small cells suppressed)
#   b. cross-sectional CNV x BMI on height and head circumference: BMI slope by age band in carriers and controls
#      (lm, as band_slopes() of 17_main_with_ukbb.R), pooled "All four CNVs" with a locus term (82_growth_blocks_v2.R),
#      the all-ages age-adjusted CNV x BMI term (mgcv bam, REML), genotype-by-band F test, heterogeneity tests
#   c. sex: CNV x BMI within each sex and the CNV x BMI x sex term, per band and all ages (82_growth_blocks_v2.R)
#   d. longitudinal next-visit models for height and head circumference (80_longitudinal_v2.R: same cleaning rules,
#      intervals, age windows, child random intercept), per genotype and pooled, BMI x window and heterogeneity tests
#   e. autism diagnosis, if given: height / head circumference / BMI Z by diagnosis; logistic regression of diagnosis
#      on BMI Z and head-circumference Z; proportion with a diagnosis by BMI band (85_v2_asd_by_cnv.R)
#   f. autistic traits, if SRS or SCQ are given for carriers AND controls: Box-Cox scale with fixed exponents,
#      BMI slope per group, CNV x BMI versus controls, ordinal check (51b_cog_panels_boxcox_combined.R)
#   g. federated trajectories: per outcome x genotype x BMI tail, binned means on a fixed age grid and an unpenalised
#      least-squares curve on a FIXED natural-cubic-spline basis of log(age) (coefficients + covariance), so that a
#      pooled curve can be drawn without any individual value leaving the site
# Inputs: one long table of visits (see README.md) and a JSON config (config_example.json). Needs, in the same folder,
#   00_federated_common.R, 00_scoring.R, 00_select_visit.R and growth_refs_lms.csv.
#   R packages: data.table, mgcv, MASS, jsonlite, splines; lmerTest (longitudinal); ordinal (optional, ordinal check).
# Outputs (aggregates only, in out_dir): sample_n.csv, estimates.csv, tests.csv, autism_by_bmi_band.csv,
#   trajectory_bins.csv, trajectory_spline.csv, trajectory_spline_matrices.csv, scoring_info.csv, run_info.csv.
#   Every row carries the site name, script version, run date, scoring and estimator. No ids, no dates, no individual values.
# Scoring: children = WHO 2006 / 2007 height and BMI, UK 1990 head circumference (00_scoring.R). Adults (height and
#   BMI from 19 y, head circumference from 18 y) = Box-Cox with the FIXED UK Biobank exponents, then residual of
#   lm(y ~ age + age^2) / residual SD in the site's non-carrier adults of the same sex (score_adult_bc of 00_scoring.R).
# null model: none; parametric tests. Bands: t test of the CNV x BMI (or CNV x BMI x sex) coefficient and F tests of
#   locus x BMI / genotype x BMI in least-squares models (plus family-clustered sandwich SE); all ages: t test of the
#   parametric coefficient in a bam (REML) model; longitudinal: Satterthwaite t test of the BMI coefficient (lmerTest)
#   and likelihood-ratio tests (ML) of BMI x window and BMI x genotype; diagnosis: Wald tests (gam, logistic), Cochran Q;
#   traits: Wald test of carrier:BMI (bam) and of the same term in a cumulative-logit model.
# Usage: Rscript site_export.R --config my_config.json [--input F] [--out DIR] [--site LABEL] [--min_cell N]
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(splines); library(MASS); library(jsonlite) })

`%||%` <- function(a, b) if (is.null(a)) b else a
.fed_dir <- function() {
  if (exists("FEDERATED_DIR")) return(FEDERATED_DIR)   # set it before source()-ing this file from another script
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else getwd()
}
HERE <- .fed_dir()
source(file.path(HERE, "00_federated_common.R"))
source(file.path(HERE, "00_scoring.R"))        # child LMS scorers and score_adult_bc() (reads growth_refs_lms.csv from HERE)
source(file.path(HERE, "00_select_visit.R"))   # select_visit(): one visit per person, rule fill_gaps or latest

DEFAULTS <- list(
  site_label = "site", input_file = NULL, out_dir = "federated_out", age_unit = "years",
  columns = list(id = "id", family = NULL, sex = "sex", age = "age", height_cm = "height_cm", weight_kg = "weight_kg", bmi = NULL, hc_cm = "hc_cm",
                 genotype = "genotype", asd = NULL, srs_total = NULL, scq_total = NULL, age_questionnaire = NULL, cohort = NULL),
  sex_values = list(male = c("male", "Male", "MALE", "M", "m", "1"), female = c("female", "Female", "FEMALE", "F", "f", "2")),
  genotype_map = list(`16p11.2 del` = "16p11.2 del", `16p11.2 dup` = "16p11.2 dup", `22q11.2 del` = "22q11.2 del", `22q11.2 dup` = "22q11.2 dup",
                      control = c("control", "No CNV")),
  asd_values = list(yes = c("1", "TRUE", "true", "yes", "Yes"), no = c("0", "FALSE", "false", "no", "No")),
  min_cell = 5, visit_rule = "fill_gaps", visit_seed = 20261119, adult_boxcox = "fixed_ukbb", min_adult_ref = 100,
  min_model_carriers = 15, min_model_controls = 30, min_longitudinal_children = 15, min_curve_n = 20, min_trait_ref = 30,
  srs_adult_form_from_age = 18, bmi_window_years = 3,
  plausible = list(height_cm = c(30, 250), weight_kg = c(1, 400), bmi = c(5, 100), hc_cm = c(20, 80)),
  sections = c("sample", "cross_sectional", "sex", "longitudinal", "autism_diagnosis", "autistic_traits", "trajectories"))
# Box-Cox exponents, fixed (never estimated at a site unless adult_boxcox = "estimate"):
#   adults: estimated once in UK Biobank non-carriers (17_main_with_ukbb.R, run of 2026-11-20); head circumference untransformed
#   SRS / SCQ raw total + 1: estimated once in SSC unaffected relatives (51b_cog_panels_boxcox_combined.R)
ADULT_LAMBDA_UKBB <- list(height = c(male = 0.65, female = 0.60), bmi = c(male = -0.75, female = -1.05), hc = c(male = 1, female = 1))
TRAIT_LAMBDA <- c(`SRS child form` = 0.32, `SRS adult form` = 0.29, SCQ = -0.16)

read_table <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.frame(readRDS(path)))
  as.data.frame(fread(path, data.table = FALSE, check.names = FALSE, na.strings = c("", "NA", "NaN", "nan", "None")))
}

# ---------------------------------------------------------------------------------------------- input
prep_data <- function(raw, cfg) {
  raw <- as.data.frame(raw); cc <- cfg$columns; n0 <- nrow(raw); info <- list(n_rows_input = n0)
  miss <- setdiff(c(cc$id, cc$sex, cc$age, cc$genotype), names(raw)); if (length(miss)) stop("input table: required column(s) not found: ", paste(miss, collapse = ", "))
  num <- function(nm) if (!is.null(nm) && nm %in% names(raw)) suppressWarnings(as.numeric(raw[[nm]])) else rep(NA_real_, n0)
  chr <- function(nm) if (!is.null(nm) && nm %in% names(raw)) as.character(raw[[nm]]) else rep(NA_character_, n0)
  sx <- as.character(raw[[cc$sex]]); sex <- ifelse(sx %in% unlist(cfg$sex_values$male), "male", ifelse(sx %in% unlist(cfg$sex_values$female), "female", NA_character_))
  unit <- switch(cfg$age_unit, years = 1, months = 12, days = 365.25, stop("age_unit must be years, months or days"))
  g <- as.character(raw[[cc$genotype]]); locus <- rep(NA_character_, n0)
  for (k in names(cfg$genotype_map)) { key <- if (k == "control") CTL else k
    if (!key %in% GROUPS5) stop("genotype_map: unknown group '", k, "'; use ", paste(c(LOCI, "control"), collapse = ", "))
    locus[g %in% as.character(unlist(cfg$genotype_map[[k]]))] <- key }
  a <- chr(cc$asd); asd <- ifelse(a %in% unlist(cfg$asd_values$yes), 1L, ifelse(a %in% unlist(cfg$asd_values$no), 0L, NA_integer_))
  d <- data.table(id = as.character(raw[[cc$id]]), family = chr(cc$family), sex = sex, age = num(cc$age) / unit, HEIGHT_raw = num(cc$height_cm), weight = num(cc$weight_kg),
                  BMI_raw = num(cc$bmi), HEADCIRCUM_raw = num(cc$hc_cm), locus = locus, cohort = chr(cc$cohort), asd = asd, srs = num(cc$srs_total), scq = num(cc$scq_total),
                  age_q = num(cc$age_questionnaire) / unit)
  info$has_family <- !all(is.na(d$family)); info$has_cohort <- !all(is.na(d$cohort))
  d[is.na(family), family := id]; d[is.na(cohort), cohort := cfg$site_label]
  for (v in c("HEIGHT_raw", "weight", "HEADCIRCUM_raw")) { m <- median(d[[v]], na.rm = TRUE); hi <- c(HEIGHT_raw = 220, weight = 250, HEADCIRCUM_raw = 70)[[v]]; lo <- c(HEIGHT_raw = 30, weight = 1.5, HEADCIRCUM_raw = 12)[[v]]
    if (is.finite(m) && (m < lo || m > hi)) stop(sprintf("input table: median %s = %.3g; expected cm (height, head circumference) and kg (weight)", v, m)) }
  pr <- cfg$plausible; nimp <- c(height_cm = 0, weight_kg = 0, bmi = 0, hc_cm = 0)
  for (v in names(nimp)) { col <- c(height_cm = "HEIGHT_raw", weight_kg = "weight", bmi = "BMI_raw", hc_cm = "HEADCIRCUM_raw")[[v]]; r <- unlist(pr[[v]])
    bad <- !is.na(d[[col]]) & (d[[col]] < r[1] | d[[col]] > r[2]); nimp[[v]] <- sum(bad); set(d, which(bad), col, NA_real_) }
  d[is.na(BMI_raw) & !is.na(weight) & !is.na(HEIGHT_raw), BMI_raw := weight / (HEIGHT_raw / 100)^2]
  bad <- !is.na(d$BMI_raw) & (d$BMI_raw < pr$bmi[[1]] | d$BMI_raw > pr$bmi[[2]]); nimp[["bmi"]] <- nimp[["bmi"]] + sum(bad); set(d, which(bad), "BMI_raw", NA_real_)
  info$n_values_outside_plausible_range <- nimp
  info$n_rows_genotype_not_in_map <- sum(is.na(d$locus)); info$n_rows_missing_sex_or_age <- sum(!is.na(d$locus) & (is.na(d$sex) | is.na(d$age)))
  d <- d[!is.na(locus) & !is.na(sex) & !is.na(age)]
  inc <- d[, .(k = uniqueN(sex) > 1 || uniqueN(locus) > 1), by = id][k == TRUE, id]; info$n_persons_inconsistent_sex_or_genotype <- length(inc); d <- d[!id %in% inc]
  d[, asd := if (all(is.na(asd))) NA_integer_ else max(asd, na.rm = TRUE), by = id]
  d[, weight := NULL]; d[, locus := factor(locus, GROUPS5)]
  info$n_rows_used <- nrow(d); info$n_persons <- uniqueN(d$id); attr(d, "info") <- info; d
}

# ---------------------------------------------------------------------------------------------- scoring
# adults: score_adult_bc() of 00_scoring.R, called per sex with the exponent fixed (mode "fixed_ukbb") or estimated in the
# site's own non-carrier adults (mode "estimate": grid -2..2; this is what 17_main_with_ukbb.R falls back to for a table
# without UK Biobank, e.g. BMI in its head-circumference table). Head circumference: exponent 1 in both modes.
score_adult_site <- function(v, age, sex, is_ref, trait, cfg) {
  min_age <- if (trait == "hc") 18 else 19; z <- rep(NA_real_, length(v)); info <- list()
  for (s in c("male", "female")) { i <- which(sex == s); if (!length(i)) next
    fixed <- trait == "hc" || cfg$adult_boxcox == "fixed_ukbb"
    zz <- score_adult_bc(v[i], age[i], sex[i], rep("site", length(i)), is_ref[i], min_age = min_age, min_ref = cfg$min_adult_ref,
                         lambda_from = "(none)", lambda_fixed = if (fixed) ADULT_LAMBDA_UKBB[[trait]][[s]] else NULL)
    z[i] <- as.numeric(zz); a <- attr(zz, "info")
    nref <- sum(is_ref[i] %in% TRUE & !is.na(v[i]) & v[i] > 0 & age[i] >= min_age)
    info[[s]] <- data.table(trait = trait, age_group = sprintf(">= %d y", min_age), sex = s, n_reference_adults = nref, n_scored = sum(!is.na(zz) & age[i] >= min_age),
      boxcox_exponent = if (is.null(a)) NA_real_ else a$lambda[1],
      exponent_source = if (is.null(a)) NA_character_ else if (trait == "hc") "none (not transformed)" else if (fixed) "fixed, UK Biobank non-carriers" else "estimated in this site's non-carrier adults",
      note = if (is.null(a)) sprintf("adults not scored: fewer than %d non-carrier adults of this sex with the measurement", cfg$min_adult_ref) else "")
  }
  list(z = z, info = rbindlist(info))
}
score_table <- function(x, rawcol, trait, cfg) {   # x: one row per person; adds y (outcome Z) and BMI (BMI Z); as 17_main_with_ukbb.R
  ref <- x$locus == CTL; sx <- as.character(x$sex)
  zc <- switch(trait, height = score_height(x[[rawcol]], x$age, sx), hc = score_hc(x[[rawcol]], x$age, sx))
  za <- score_adult_site(x[[rawcol]], x$age, sx, ref, trait, cfg); ad <- x$age >= (if (trait == "hc") 18 else 19); zc[ad] <- za$z[ad]
  bc <- score_bmi_child(x$BMI_raw, x$age, sx); ba <- score_adult_site(x$BMI_raw, x$age, sx, ref, "bmi", cfg); ad2 <- x$age >= 19; bc[ad2] <- ba$z[ad2]
  x <- copy(x)[, `:=`(y = zc, BMI = bc)]
  attr(x, "scoring") <- rbind(za$info, ba$info); x
}
SCORING_GROWTH <- function(cfg) sprintf("children: WHO 2006/2007 LMS Z (height, BMI), UK 1990 LMS Z (head circumference); adults (height, BMI >= 19 y; head circumference >= 18 y): Box-Cox (%s), residual of lm(y ~ age + age^2) / SD in the site's non-carrier adults of the same sex",
  if (cfg$adult_boxcox == "fixed_ukbb") "fixed UK Biobank exponents: height 0.65 men, 0.60 women; BMI -0.75 men, -1.05 women; head circumference untransformed" else "exponent estimated in the site's non-carrier adults; head circumference untransformed")
SCORING_CHILD <- "WHO 2006/2007 LMS Z (height, BMI), UK 1990 LMS Z (head circumference); visits before 18 y only"
SCORING_TRAIT <- "SRS / SCQ raw total + 1, Box-Cox with fixed exponents (SRS child form 0.32, SRS adult form 0.29, SCQ -0.16), Z against the site's controls of the same form; BMI Z as in the growth models"

# ---------------------------------------------------------------------------------------------- small cells
fmt_n <- function(n, minc) ifelse(is.na(n), NA_character_, ifelse(n > 0 & n < minc, paste0("<", minc), as.character(n)))
# closure of a table of counts whose LAST row and LAST column are totals: a hidden cell (0 < n < minc) must not be
# recoverable by subtraction, so every line (parts + their total) with exactly one hidden cell gets a second one.
close_table <- function(n, minc) {
  P <- n > 0 & n < minc; S <- P; R <- nrow(n); C <- ncol(n)
  lines <- c(lapply(seq_len(R), function(r) cbind(r, seq_len(C))), lapply(seq_len(C), function(cc) cbind(seq_len(R), cc)))
  repeat { changed <- FALSE
    for (ln in lines) { s <- S[ln]; if (sum(s) != 1) next
      cand <- which(!s & n[ln] > 0); if (!length(cand)) next
      S[ln[cand[which.min(n[ln][cand])], , drop = FALSE]] <- TRUE; changed <- TRUE }
    if (!changed) break }
  list(hidden = P, masked = S & !P)
}
close_line <- function(parts, minc) {   # parts of a published total; returns TRUE for parts whose count must not be shown
  h <- parts > 0 & parts < minc
  if (sum(h) == 1) { cand <- which(!h & parts > 0); if (length(cand)) h[cand[which.min(parts[cand])]] <- TRUE }
  h
}

# ---------------------------------------------------------------------------------------------- output rows
erow <- function(analysis, outcome, group, stratum, term, model, estimator, est = NA_real_, se = NA_real_, p = NA_real_, df = NA_real_, sex = "all",
                 n_carriers = NA, n_controls = NA, n_pairs = NA, n_asd = NA, se_cl = NA_real_, p_cl = NA_real_, n_fam = NA_real_, status = "ok", in_paper = "yes")
  data.table(analysis = analysis, outcome = outcome, group = group, stratum = stratum, sex = sex, term = term, model = model, estimator = estimator,
             estimate = as.numeric(est), se = as.numeric(se), df = as.numeric(df), p = as.numeric(p), se_family_clustered = as.numeric(se_cl), p_family_clustered = as.numeric(p_cl),
             n_families = as.numeric(n_fam), n_carriers = as.character(n_carriers), n_controls = as.character(n_controls), n_pairs = as.character(n_pairs), n_asd = as.character(n_asd),
             status = status, same_model_as_paper = in_paper)
trow <- function(analysis, outcome, group, stratum, test, model, estimator, statistic = NA_real_, df = NA_character_, p = NA_real_, detail = "", status = "ok")
  data.table(analysis = analysis, outcome = outcome, group = group, stratum = stratum, test = test, model = model, estimator = estimator, statistic = as.numeric(statistic),
             df = as.character(df), p = as.numeric(p), detail = detail, status = status)
co4 <- function(m, term) { s <- summary(m); tb <- if (inherits(m, "gam")) s$p.table else s$coefficients
  if (!term %in% rownames(tb)) return(c(NA_real_, NA_real_, NA_real_)); unname(tb[term, c(1, 2, 4)]) }
inband <- function(x, b) x[age >= BANDS[b] & age < BANDS[b + 1]]
cochran <- function(b, se) { f <- fe_meta(b, se); c(Q = f$Q, df = f$Q_df, p = f$Q_p, I2 = f$I2) }

# ============================================================================================== driver
# keep_data = TRUE attaches the scored one-visit-per-person tables to the returned object IN MEMORY (attribute "analysis_data");
# it is used by tests/ only and nothing of it is ever written to disk.
run_site_export <- function(cfg = list(), data = NULL, write = TRUE, keep_data = FALSE) {
  cfg <- modifyList(DEFAULTS, cfg[!grepl("^_", names(cfg))])
  stopifnot(cfg$visit_rule %in% c("fill_gaps", "latest"), cfg$adult_boxcox %in% c("fixed_ukbb", "estimate"))
  minc <- cfg$min_cell; MINC <- max(cfg$min_model_carriers, minc); MINK <- max(cfg$min_model_controls, minc); MINKIDS <- max(cfg$min_longitudinal_children, minc)
  SITE <- cfg$site_label; SEC <- cfg$sections
  D <- prep_data(data %||% read_table(cfg$input_file), cfg); INFO <- attr(D, "info"); HAVE_FAM <- INFO$has_family
  message(sprintf("site %s: %d persons, %d visit rows", SITE, INFO$n_persons, INFO$n_rows_used))
  EST <- list(); TST <- list(); NOTES <- character(0)
  addE <- function(x) if (!is.null(x) && nrow(x)) EST[[length(EST) + 1]] <<- x
  addT <- function(x) if (!is.null(x) && nrow(x)) TST[[length(TST) + 1]] <<- x
  note <- function(...) { m <- sprintf(...); message("  note: ", m); NOTES <<- c(NOTES, m) }
  SG <- SCORING_GROWTH(cfg)
  # family-clustered (sandwich, CR1) SE of one lm coefficient; P from t with (clusters - 1) df (cl_se of 82_growth_blocks_v2.R)
  cl_se <- function(m, cluster, term) { if (!HAVE_FAM || anyNA(cluster)) return(c(NA_real_, NA_real_, NA_real_))
    X <- model.matrix(m); ok <- !is.na(coef(m)); X <- X[, ok, drop = FALSE]; u <- resid(m); A <- chol2inv(chol(crossprod(X))); dimnames(A) <- list(colnames(X), colnames(X))
    S <- rowsum(X * u, cluster); Gn <- nrow(S); N <- nrow(X); K <- ncol(X); V <- A %*% crossprod(S) %*% A * (Gn / (Gn - 1)) * ((N - 1) / (N - K)); dimnames(V) <- dimnames(A)
    se <- sqrt(V[term, term]); c(se, 2 * pt(-abs(coef(m)[[term]] / se), Gn - 1), Gn) }
  sx <- function(x) if (uniqueN(x$sex) > 1) "+ sex" else ""

  # ---------------------------------------------------------------- one visit per person, scoring (as 17_main_with_ukbb.R)
  Sys.setenv(VISIT_SEED = as.character(cfg$visit_seed))
  cl <- D[age >= 2 & age <= 73]
  pick1 <- function(rawcol) { if (!any(!is.na(cl[[rawcol]]) & !is.na(cl$BMI_raw))) return(NULL)
    s <- select_visit(cl, rawcol, rule = cfg$visit_rule); s[, c("band", "rule_changed") := NULL]; s }
  mk_cs <- function(rawcol, trait) { s <- pick1(rawcol); if (is.null(s)) return(NULL)
    x <- score_table(s, rawcol, trait, cfg); sc <- attr(x, "scoring"); x <- x[!is.na(y) & !is.na(BMI)]
    x <- x[, .(id, fam = family, locus = as.character(locus), sex, age, y, BMI, cohort, asd)]
    list(car = x[locus %in% LOCI], ctl = x[locus == CTL], scoring = sc) }
  DATA2 <- list(Height = mk_cs("HEIGHT_raw", "height"), `Head circumference` = mk_cs("HEADCIRCUM_raw", "hc"))
  DATA2 <- DATA2[!vapply(DATA2, is.null, logical(1))]
  SCI <- rbindlist(lapply(names(DATA2), function(o) cbind(table = paste(o, "cross-sectional"), DATA2[[o]]$scoring)))

  # ---------------------------------------------------------------- a. sample
  NST <- list()   # display status of n per outcome | group | band | sex: "hidden" (< min), "masked" (protects a hidden cell) or absent
  SMP <- list()
  for (o in names(DATA2)) { x <- rbind(DATA2[[o]]$car, DATA2[[o]]$ctl)[, bl := as.character(cut(age, BANDS, labels = BLAB, right = FALSE))]
    tabs <- list()
    for (G in GROUPS5) { xg <- x[locus == G]
      n <- matrix(0, length(BLAB) + 1, 3, dimnames = list(c(BLAB, "all ages"), c("male", "female", "all")))
      for (b in BLAB) for (s in c("male", "female")) n[b, s] <- xg[bl == b & sex == s, .N]
      n[, "all"] <- n[, "male"] + n[, "female"]; n["all ages", ] <- colSums(n[BLAB, , drop = FALSE]); tabs[[G]] <- list(n = n, cl = close_table(n, minc)) }
    npool <- Reduce(`+`, lapply(tabs[LOCI], `[[`, "n")); anyh <- Reduce(`|`, lapply(tabs[LOCI], function(t) t$cl$hidden | t$cl$masked))
    tabs[[POOL]] <- list(n = npool, cl = list(hidden = npool > 0 & npool < minc, masked = anyh & !(npool > 0 & npool < minc)))
    for (G in names(tabs)) for (b in rownames(tabs[[G]]$n)) for (s in colnames(tabs[[G]]$n)) { t <- tabs[[G]]
      st <- if (t$cl$hidden[b, s]) "hidden" else if (t$cl$masked[b, s]) "masked" else "shown"; if (st != "shown") NST[[paste(o, G, b, s, sep = "|")]] <- st
      SMP[[length(SMP) + 1]] <- data.table(outcome = o, group = G, age_band = b, sex = s, n = if (st == "hidden") paste0("<", minc) else if (st == "masked") paste0(">=", minc) else as.character(t$n[b, s]),
        status = c(shown = "ok", hidden = sprintf("suppressed: fewer than %d", minc), masked = "count masked to protect a small cell in the same table")[[st]]) }
  }
  SMP <- rbindlist(SMP)
  nshow <- function(n, o, G, b = "all ages", s = "all") { st <- NST[[paste(o, G, b, s, sep = "|")]]
    if (is.null(st)) as.character(n) else if (st == "hidden") paste0("<", minc) else paste0(">=", minc) }
  nlow <- function(n, thr) if (n == 0) "0" else paste0("<", thr)

  # ---------------------------------------------------------------- b. cross-sectional
  if ("cross_sectional" %in% SEC) for (o in names(DATA2)) { car <- DATA2[[o]]$car; ctl <- DATA2[[o]]$ctl; message("  cross-sectional: ", o)
    LM <- "OLS (lm)"; BAM <- "GAM, REML (mgcv bam)"
    for (b in seq_along(BLAB)) { bl <- BLAB[b]; k0 <- inband(ctl, b); cb <- inband(car, b); nk <- nshow(nrow(k0), o, CTL, bl)
      few <- sprintf("not estimated: fewer than %d carriers or %d controls in this band", MINC, MINK)
      if (nrow(k0) >= MINK) addE(tryCatch({ m <- lm(as.formula(paste("y ~ BMI + age", sx(k0))), data = k0); cf <- summary(m)$coefficients["BMI", ]
        erow("cross-sectional band", o, CTL, bl, "BMI slope (controls-only model)", "Z ~ BMI Z + sex + age, controls", LM, cf[1], cf[2], cf[4], m$df.residual, n_controls = nk) },
        error = function(e) erow("cross-sectional band", o, CTL, bl, "BMI slope (controls-only model)", "Z ~ BMI Z + sex + age, controls", LM, n_controls = nk, status = paste("model failed:", conditionMessage(e)))))
      # per genotype: band_slopes() of 17 (y ~ cnv * BMI + sex + age)
      for (G in LOCI) { c1 <- cb[locus == G]; MOD <- "Z ~ CNV * BMI Z + sex + age, carriers of one genotype + controls"
        if (nrow(c1) < MINC || nrow(k0) < MINK) { addE(erow("cross-sectional band", o, G, bl, "CNV x BMI", MOD, LM, n_carriers = if (nrow(c1) < MINC) nlow(nrow(c1), MINC) else nshow(nrow(c1), o, G, bl), n_controls = if (nrow(k0) < MINK) nlow(nrow(k0), MINK) else nk, status = few)); next }
        addE(tryCatch({
          x <- rbind(c1[, .(y, BMI, sex, age, cnv = 1, fam)], k0[, .(y, BMI, sex, age, cnv = 0, fam)])[, cnvBMI := cnv * BMI]
          m <- lm(as.formula(paste("y ~ cnv + BMI + cnvBMI + age", sx(x))), data = x); B <- coef(m); V <- vcov(m); cf <- summary(m)$coefficients; v <- cl_se(m, x$fam, "cnvBMI"); nc <- nshow(nrow(c1), o, G, bl)
          sc <- sqrt(V["BMI", "BMI"] + V["cnvBMI", "cnvBMI"] + 2 * V["BMI", "cnvBMI"]); bcar <- B[["BMI"]] + B[["cnvBMI"]]
          rbind(erow("cross-sectional band", o, G, bl, "BMI slope, controls", MOD, LM, B[["BMI"]], sqrt(V["BMI", "BMI"]), cf["BMI", 4], m$df.residual, n_carriers = nc, n_controls = nk),
                erow("cross-sectional band", o, G, bl, "BMI slope, carriers", MOD, LM, bcar, sc, 2 * pt(-abs(bcar / sc), m$df.residual), m$df.residual, n_carriers = nc, n_controls = nk),
                erow("cross-sectional band", o, G, bl, "CNV x BMI", MOD, LM, B[["cnvBMI"]], cf["cnvBMI", 2], cf["cnvBMI", 4], m$df.residual, n_carriers = nc, n_controls = nk, se_cl = v[1], p_cl = v[2], n_fam = v[3])) },
          error = function(e) erow("cross-sectional band", o, G, bl, "CNV x BMI", MOD, LM, status = paste("model failed:", conditionMessage(e)))))
      }
      # pooled, all four CNVs (pooled_band of 82): locus-adjusted (primary) and without a locus term
      if (nrow(cb) >= MINC && nrow(k0) >= MINK) for (adjust in c(TRUE, FALSE)) addE(tryCatch({
        MOD <- if (adjust) "Z ~ locus (controls = reference) + BMI Z + CNV x BMI Z + sex + age (primary)" else "Z ~ CNV * BMI Z + sex + age, no locus term"
        x <- rbind(cb[, .(y, BMI, sex, age, cnv = 1, loc = locus, fam)], k0[, .(y, BMI, sex, age, cnv = 0, loc = "Controls", fam)])
        x[, `:=`(loc = relevel(factor(loc), "Controls"), cnvBMI = cnv * BMI)]
        m <- lm(as.formula(paste("y ~", if (adjust) "loc" else "cnv", "+ BMI + cnvBMI + age", sx(x))), data = x); B <- coef(m); V <- vcov(m); cf <- summary(m)$coefficients; v <- cl_se(m, x$fam, "cnvBMI"); nc <- nshow(nrow(cb), o, POOL, bl)
        sc <- sqrt(V["BMI", "BMI"] + V["cnvBMI", "cnvBMI"] + 2 * V["BMI", "cnvBMI"]); bcar <- B[["BMI"]] + B[["cnvBMI"]]
        rbind(erow("cross-sectional band", o, POOL, bl, "BMI slope, controls", MOD, LM, B[["BMI"]], sqrt(V["BMI", "BMI"]), cf["BMI", 4], m$df.residual, n_carriers = nc, n_controls = nk),
              erow("cross-sectional band", o, POOL, bl, "BMI slope, carriers", MOD, LM, bcar, sc, 2 * pt(-abs(bcar / sc), m$df.residual), m$df.residual, n_carriers = nc, n_controls = nk),
              erow("cross-sectional band", o, POOL, bl, "CNV x BMI", MOD, LM, B[["cnvBMI"]], cf["cnvBMI", 2], cf["cnvBMI", 4], m$df.residual, n_carriers = nc, n_controls = nk, se_cl = v[1], p_cl = v[2], n_fam = v[3])) },
        error = function(e) erow("cross-sectional band", o, POOL, bl, "CNV x BMI", "pooled", LM, status = paste("model failed:", conditionMessage(e)))))
      # heterogeneity across genotypes among the carriers of the band (genotypes with enough carriers)
      kl <- cb[, .N, by = locus][N >= MINC, locus]
      if (length(kl) >= 2) addT(tryCatch({ ch <- cb[locus %in% kl][, loc := factor(locus)]
        a <- anova(lm(as.formula(paste("y ~ loc + BMI + age", sx(ch))), data = ch), lm(as.formula(paste("y ~ loc * BMI + age", sx(ch))), data = ch))
        trow("cross-sectional band", o, POOL, bl, "heterogeneity of the BMI slope across genotypes, carriers only (F)", "Z ~ locus * BMI Z + sex + age vs Z ~ locus + BMI Z + sex + age", LM, a$F[2], sprintf("%d, %d", a$Df[2], a$Res.Df[2]), a$`Pr(>F)`[2], paste(sort(kl), collapse = "; ")) },
        error = function(e) NULL))
      # genotype-by-band F test, controls + genotypes (bgt of 17)
      if (nrow(k0) >= MINK && length(kl) >= 1) addT(tryCatch({ X <- rbind(cb[locus %in% kl], k0)[, grp := relevel(factor(ifelse(locus == CTL, "Controls", locus)), "Controls")]
        a <- anova(lm(as.formula(paste("y ~ grp + BMI + age", sx(X))), data = X), lm(as.formula(paste("y ~ grp * BMI + age", sx(X))), data = X))
        trow("cross-sectional band", o, "Controls + genotypes", bl, "BMI effect differs by genotype (F)", "Z ~ group * BMI Z + sex + age vs Z ~ group + BMI Z + sex + age", LM, a$F[2], sprintf("%d, %d", a$Df[2], a$Res.Df[2]), a$`Pr(>F)`[2], paste(sort(kl), collapse = "; ")) },
        error = function(e) NULL))
    }
    # ---- all ages, age-adjusted (est_main of 17; pooled forms of 82)
    nk <- nshow(nrow(ctl), o, CTL)
    if (nrow(ctl) >= MINK) addE(tryCatch({ v <- co4(bam(y ~ BMI + s(age, k = 15), data = ctl, method = "REML"), "BMI")
      erow("cross-sectional all ages", o, CTL, "all ages", "BMI slope (controls-only model)", "Z ~ BMI Z + s(age), controls", BAM, v[1], v[2], v[3], n_controls = nk) },
      error = function(e) erow("cross-sectional all ages", o, CTL, "all ages", "BMI slope (controls-only model)", "Z ~ BMI Z + s(age), controls", BAM, n_controls = nk, status = paste("model failed:", conditionMessage(e)))))
    MODA <- "Z ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z"
    for (G in LOCI) { xc <- car[locus == G]
      if (nrow(xc) < MINC || nrow(ctl) < MINK) { addE(erow("cross-sectional all ages", o, G, "all ages", "CNV x BMI", MODA, BAM, n_carriers = if (nrow(xc) < MINC) nlow(nrow(xc), MINC) else nshow(nrow(xc), o, G), n_controls = if (nrow(ctl) < MINK) nlow(nrow(ctl), MINK) else nk,
                                                      status = sprintf("not estimated: fewer than %d carriers or %d controls", MINC, MINK))); next }
      addE(tryCatch({ x <- rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == CTL, "No CNV", "CNV"), levels = c("No CNV", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
        v <- co4(bam(y ~ is_cnv + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI, data = x, method = "REML"), "cnvBMI")
        erow("cross-sectional all ages", o, G, "all ages", "CNV x BMI", MODA, BAM, v[1], v[2], v[3], n_carriers = nshow(nrow(xc), o, G), n_controls = nk) },
        error = function(e) erow("cross-sectional all ages", o, G, "all ages", "CNV x BMI", MODA, BAM, status = paste("model failed:", conditionMessage(e)))))
    }
    if (nrow(car) >= MINC && nrow(ctl) >= MINK) {
      x <- rbind(ctl[, .(y, BMI, sex, age, locus)], car[, .(y, BMI, sex, age, locus)])
      x[, `:=`(is_cnv = factor(ifelse(locus == CTL, "No CNV", "CNV"), levels = c("No CNV", "CNV")), grp = factor(ifelse(locus == CTL, "Controls", locus), c("Controls", intersect(LOCI, unique(car$locus)))))]
      x[, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
      forms <- list("Z ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z, no locus term" = y ~ is_cnv + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI,
                    "Z ~ locus + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z (primary)" = y ~ grp + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI,
                    "Z ~ locus + s(age, by = locus) + s(age, by = BMI Z) + CNV x BMI Z" = y ~ grp + s(age, k = 10, by = grp) + s(age, k = 10, by = BMI) + cnvBMI)
      for (nm in names(forms)) addE(tryCatch({ v <- co4(bam(forms[[nm]], data = x, method = "REML"), "cnvBMI")
        erow("cross-sectional all ages", o, POOL, "all ages", "CNV x BMI", nm, BAM, v[1], v[2], v[3], n_carriers = nshow(nrow(car), o, POOL), n_controls = nk) },
        error = function(e) erow("cross-sectional all ages", o, POOL, "all ages", "CNV x BMI", nm, BAM, status = paste("model failed:", conditionMessage(e)))))
      # heterogeneity across genotypes, carriers only (38_fig6_numbers.R) and Cochran Q on the per-genotype terms
      kl <- car[, .N, by = locus][N >= MINC, locus]
      if (length(kl) >= 2) { addT(tryCatch({ ch <- car[locus %in% kl][, locus := factor(locus)]
          a <- anova(gam(y ~ s(age, k = 8, by = locus) + locus + BMI + locus:BMI, data = ch, method = "REML"))$pTerms.table
          trow("cross-sectional all ages", o, POOL, "all ages", "heterogeneity of the BMI slope across genotypes, carriers only (Wald F)", "Z ~ s(age, by = locus) + locus + BMI Z + locus:BMI Z", "GAM, REML (mgcv gam)", a["locus:BMI", "F"], a["locus:BMI", "df"], a["locus:BMI", "p-value"], paste(sort(kl), collapse = "; ")) },
          error = function(e) NULL))
        e <- rbindlist(EST)[analysis == "cross-sectional all ages" & outcome == o & group %in% LOCI & term == "CNV x BMI" & status == "ok"]
        if (nrow(e) >= 2) { q <- cochran(e$estimate, e$se)
          addT(trow("cross-sectional all ages", o, POOL, "all ages", "heterogeneity of CNV x BMI across genotypes (Cochran Q)", MODA, "Cochran Q on the per-genotype estimates", q[["Q"]], q[["df"]], q[["p"]], paste(e$group, collapse = "; "))) }
      }
    }
  }

  # ---------------------------------------------------------------- c. sex (82_growth_blocks_v2.R)
  if ("sex" %in% SEC) for (o in names(DATA2)) { message("  sex: ", o)
    mk <- function(xc, ctl) { x <- rbind(ctl[, .(y, BMI, sex, age, cnv = 0, loc = "Controls", fam)], xc[, .(y, BMI, sex, age, cnv = 1, loc = locus, fam)])
      x[, `:=`(loc = relevel(factor(loc), "Controls"), is_cnv = factor(ifelse(cnv == 1, "CNV", "No CNV"), levels = c("No CNV", "CNV")), male = as.numeric(sex == "male"))]
      x[, `:=`(cnvBMI = cnv * BMI, BMImale = BMI * male, cnvBMImale = cnv * BMI * male, agemale = age * male, cell = interaction(is_cnv, sex, drop = TRUE))]; x }
    for (G in c(LOCI, POOL)) { pooled <- G == POOL; xc0 <- if (pooled) DATA2[[o]]$car else DATA2[[o]]$car[locus == G]
      for (scope in c("all ages", BLAB)) { b <- match(scope, BLAB); xc <- if (is.na(b)) xc0 else inband(xc0, b); k0 <- if (is.na(b)) DATA2[[o]]$ctl else inband(DATA2[[o]]$ctl, b)
        if (!nrow(xc) || !nrow(k0)) next
        x <- mk(xc, k0); n <- c(cm = x[cnv == 1 & male == 1, .N], cf = x[cnv == 1 & male == 0, .N], km = x[cnv == 0 & male == 1, .N], kf = x[cnv == 0 & male == 0, .N])
        if (min(n[c("cm", "cf")]) < MINC || min(n[c("km", "kf")]) < MINK) next
        an <- "sex"; est <- if (is.na(b)) "GAM, REML (mgcv bam)" else "OLS (lm)"; Gt <- if (pooled) "loc" else if (is.na(b)) "is_cnv" else "cnv"
        nn <- function(s) list(c = nshow(n[[if (s == "male") "cm" else "cf"]], o, G, scope, s), k = nshow(n[[if (s == "male") "km" else "kf"]], o, CTL, scope, s))
        addE(tryCatch({ out <- list()
          for (s in c("male", "female")) { xs <- x[sex == s]; if (pooled) xs[, loc := droplevels(loc)]; N <- nn(s)
            if (is.na(b)) { MOD <- "within sex: Z ~ CNV (or locus) + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z"
              v <- co4(bam(as.formula(sprintf("y ~ %s + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI", Gt)), data = xs, method = "REML"), "cnvBMI")
              out[[length(out) + 1]] <- erow(an, o, G, scope, "CNV x BMI", MOD, est, v[1], v[2], v[3], sex = s, n_carriers = N$c, n_controls = N$k)
            } else { MOD <- "within sex: Z ~ CNV (or locus) + BMI Z + CNV x BMI Z + age"
              m <- lm(as.formula(sprintf("y ~ %s + BMI + cnvBMI + age", Gt)), data = xs); B <- coef(m); V <- vcov(m); cf <- summary(m)$coefficients
              sc <- sqrt(V["BMI", "BMI"] + V["cnvBMI", "cnvBMI"] + 2 * V["BMI", "cnvBMI"]); bc <- B[["BMI"]] + B[["cnvBMI"]]
              out[[length(out) + 1]] <- rbind(erow(an, o, G, scope, "BMI slope, controls", MOD, est, B[["BMI"]], sqrt(V["BMI", "BMI"]), cf["BMI", 4], m$df.residual, sex = s, n_carriers = N$c, n_controls = N$k),
                                              erow(an, o, G, scope, "BMI slope, carriers", MOD, est, bc, sc, 2 * pt(-abs(bc / sc), m$df.residual), m$df.residual, sex = s, n_carriers = N$c, n_controls = N$k),
                                              erow(an, o, G, scope, "CNV x BMI", MOD, est, B[["cnvBMI"]], cf["cnvBMI", 2], cf["cnvBMI", 4], m$df.residual, sex = s, n_carriers = N$c, n_controls = N$k)) }
          }
          nc <- nshow(n[["cm"]] + n[["cf"]], o, G, scope); nk <- nshow(n[["km"]] + n[["kf"]], o, CTL, scope)
          if (is.na(b)) { MOD3 <- "Z ~ CNV (or locus) * male + s(age, by = CNV x sex cell) + s(age, by = BMI Z) + s(age, by = BMI Z x male) + CNV x BMI Z + CNV x BMI Z x male"
            m3 <- bam(as.formula(sprintf("y ~ %s * male + s(age, k = 15, by = cell) + s(age, k = 10, by = BMI) + s(age, k = 10, by = BMImale) + cnvBMI + cnvBMImale", Gt)), data = x, method = "REML"); v <- co4(m3, "cnvBMImale"); vc <- rep(NA_real_, 3)
          } else { MOD3 <- "Z ~ CNV (or locus) * male + BMI Z + BMI Z x male + CNV x BMI Z + CNV x BMI Z x male + age + age x male"
            m3 <- lm(as.formula(sprintf("y ~ %s * male + BMI + BMImale + cnvBMI + cnvBMImale + age + agemale", Gt)), data = x); v <- co4(m3, "cnvBMImale"); vc <- if (is.na(v[1])) rep(NA_real_, 3) else cl_se(m3, x$fam, "cnvBMImale") }
          out[[length(out) + 1]] <- erow(an, o, G, scope, "CNV x BMI x sex (male minus female)", MOD3, est, v[1], v[2], v[3], sex = "male - female", n_carriers = nc, n_controls = nk, se_cl = vc[1], p_cl = vc[2], n_fam = vc[3])
          rbindlist(out) }, error = function(e) erow(an, o, G, scope, "CNV x BMI x sex (male minus female)", "sex models", est, status = paste("model failed:", conditionMessage(e)))))
      } } }

  # ---------------------------------------------------------------- d. longitudinal (80_longitudinal_v2.R)
  LONGINFO <- NULL
  if ("longitudinal" %in% SEC) {
    if (!requireNamespace("lmerTest", quietly = TRUE)) note("longitudinal section skipped: package lmerTest is not installed") else {
      message("  longitudinal")
      # height: clean_visits() of 00_helpers.R (per child, in age order, against the last kept visit)
      clean_height <- function(d) { keep <- logical(nrow(d)); why <- character(nrow(d))
        for (ix in split(seq_len(nrow(d)), d$id)) { last <- ix[1]; keep[last] <- TRUE
          for (i in ix[-1]) { dt <- d$age[i] - d$age[last]; dr <- d$HEIGHT_raw[i] - d$HEIGHT_raw[last]; dz <- d$HEIGHT[i] - d$HEIGHT[last]
            w <- if (dt < 0.25) "a_duplicate_visit" else if (!is.na(dr) && abs(dr) < 0.1 && d$age[i] < 16) "b_carried_forward"
                 else if (!is.na(dr) && (dr < -2 || (d$age[i] > 3 && dr / dt > 20))) "c_impossible_change" else if (abs(dz) > 2.5) "d_jump_gt_2.5SD" else ""
            if (w == "") { keep[i] <- TRUE; last <- i } else why[i] <- w } }
        list(d = d[keep], why = why[why != ""]) }
      # head circumference: rule of 07_hc_trajectories.R
      clean_hc <- function(d) { keep <- logical(nrow(d)); why <- character(nrow(d))
        for (ix in split(seq_len(nrow(d)), d$id)) { last <- ix[1]; keep[last] <- TRUE
          for (i in ix[-1]) { dt <- d$age[i] - d$age[last]; dh <- d$hc[i] - d$hc[last]; dz <- d$hcz[i] - d$hcz[last]
            w <- if (dt < ifelse(d$age[i] < 2, 0.04, 0.25)) "duplicate_visit" else if (abs(dh) < 0.05 && dt >= 0.1 && d$age[i] < 3) "carried_forward"
                 else if (dh < -1.5 || (d$age[i] > 1 && dh / dt > 15)) "impossible_change" else if (abs(dz) > 2.5) "jump_gt_2.5SD" else ""
            if (w == "") { keep[i] <- TRUE; last <- i } else why[i] <- w } }
        list(d = d[keep], why = why[why != ""]) }
      build_height <- function(G5) {
        d <- D[locus %in% G5 & age >= 2 & age < 18 & !is.na(HEIGHT_raw) & !is.na(BMI_raw), .(id, cohort, locus = as.character(locus), sex, age, HEIGHT_raw, BMI_raw)]
        if (!nrow(d)) return(NULL)
        d[, `:=`(HEIGHT = score_height(HEIGHT_raw, age, sex), BMI = score_bmi_child(BMI_raw, age, sex))]
        d <- d[!is.na(HEIGHT) & !is.na(BMI) & abs(HEIGHT) <= 5 & abs(BMI) <= 5][order(id, age)]; n0 <- nrow(d)
        cl <- clean_height(d); d <- unique(cl$d, by = c("id", "age"))
        iv <- d[, if (.N >= 2) .(age0 = head(age, -1), dt = diff(age), Z0 = head(HEIGHT, -1), dZ = diff(HEIGHT), B0 = head(BMI, -1), sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id]
        if (!nrow(iv)) return(NULL); iv <- iv[dt >= 0.5]; iv[, win := cut(age0, WIN_H_BREAKS, right = FALSE, labels = WIN_H)]
        list(I = iv, log = data.table(outcome = "Height", visits_before_cleaning = n0, visits_kept = nrow(d), as.data.table(as.list(table(factor(cl$why, c("a_duplicate_visit", "b_carried_forward", "c_impossible_change", "d_jump_gt_2.5SD")))))))
      }
      build_hcl <- function(G5) {
        d <- D[locus %in% G5 & age >= 0.1 & age < 18, .(id, cohort, locus = as.character(locus), sex, age, hc = HEADCIRCUM_raw, bmi = BMI_raw)]
        d[!is.na(hc) & (hc < 30 | hc > 70), hc := NA_real_]; d[!is.na(bmi) & (bmi < 10 | bmi > 70), bmi := NA_real_]; d <- d[!is.na(hc)]
        if (!nrow(d)) return(NULL)
        d[, `:=`(hcz = score_hc(hc, age, sex), bmiz = score_bmi_child(bmi, age, sex))]
        d <- d[!is.na(hcz) & abs(hcz) <= 5][order(id, age)]; d[!is.na(bmiz) & abs(bmiz) > 5, bmiz := NA_real_]; n0 <- nrow(d)
        cl <- clean_hc(d); d <- cl$d
        iv <- d[!is.na(bmiz)][order(id, age)][, if (.N >= 2) .(age0 = head(age, -1), dt = diff(age), Z0 = head(hcz, -1), dZ = diff(hcz), B0 = head(bmiz, -1), sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id]
        if (!nrow(iv)) return(NULL); iv <- iv[dt >= 0.1]; iv[, win := cut(age0, WIN_C_BREAKS, right = FALSE, labels = WIN_C)]
        list(I = iv, log = data.table(outcome = "Head circumference", visits_before_cleaning = n0, visits_kept = nrow(d), as.data.table(as.list(table(factor(cl$why, c("duplicate_visit", "carried_forward", "impossible_change", "jump_gt_2.5SD")))))))
      }
      adj <- function(x, v) if (uniqueN(x[[v]]) > 1) paste("+", v) else ""
      core <- function(x) paste("dZ ~ Z0 + dt", adj(x, "sex"), adj(x, "cohort"), adj(x, "locus"), "+ ns(age0, 3)")
      fitm <- function(f, x, REML = TRUE) suppressMessages(suppressWarnings(
        if (anyDuplicated(x$id)) lmerTest::lmer(as.formula(paste(f, "+ (1 | id)")), x, REML = REML) else lm(as.formula(f), x)))
      crow <- function(m, pat) { co <- coef(summary(m)); r <- grep(pat, rownames(co)); data.table(term = rownames(co)[r], beta = co[r, 1], se = co[r, 2], p = co[r, ncol(co)]) }
      lrt <- function(m0, m1) tryCatch({ a <- suppressMessages(anova(m0, m1)); p <- a[[grep("^Pr", names(a))]]; p[2] }, error = function(e) NA_real_)
      LMM <- "linear mixed model, REML, child random intercept, Satterthwaite t (lmerTest); lm when every child has one interval"
      MODL <- "change in Z to the next visit ~ Z at the visit + interval + sex (+ cohort) (+ genotype) + ns(age at the visit, 3) + BMI Z at the visit + (1 | child)"
      nkid <- function(n) if (n > 0 && n < minc) paste0("<", minc) else as.character(n)
      run_set <- function(I, o, wins, groups, paper) { hid <- FALSE
        for (g in groups) { x <- if (g == POOL) I[locus %in% LOCI] else I[locus == g]; if (!nrow(x)) next
          kids <- uniqueN(x$id); anyhid <- g == POOL && hid
          show <- function(k, pairs) list(k = if (anyhid) paste0(">=", minc) else nkid(k), pr = if (k > 0 && k < minc) NA else if (anyhid) NA else pairs)
          if (kids > 0 && kids < minc) hid <- TRUE
          # overall BMI coefficient
          s <- s_all <- show(kids, nrow(x))
          addE(if (kids < MINKIDS) erow("longitudinal", o, g, "All ages", "BMI Z at the visit", MODL, LMM, n_carriers = s$k, n_pairs = s$pr, status = sprintf("not estimated: fewer than %d children", MINKIDS), in_paper = paper) else
            tryCatch({ r <- crow(fitm(paste(core(x), "+ B0"), x), "^B0$"); erow("longitudinal", o, g, "All ages", "BMI Z at the visit", MODL, LMM, r$beta, r$se, r$p, n_carriers = s$k, n_pairs = s$pr, in_paper = paper) },
                     error = function(e) erow("longitudinal", o, g, "All ages", "BMI Z at the visit", MODL, LMM, n_carriers = s$k, n_pairs = s$pr, status = paste("model failed:", conditionMessage(e)), in_paper = paper)))
          # window-specific BMI coefficients from one model, BMI x window LRT
          xw <- copy(x)[!is.na(win)][, win := droplevels(win)]; dsc <- xw[, .(children = uniqueN(id), pairs = .N), by = .(scope = as.character(win))]
          MODW <- sub("BMI Z at the visit [+]", "window + window:BMI Z at the visit +", MODL)
          r <- if (uniqueN(xw$id) < MINKIDS || nlevels(xw$win) < 2) NULL else tryCatch({ m <- fitm(paste(core(xw), "+ win + win:B0"), xw); k <- crow(m, ":B0$"); k[, scope := sub("^win", "", sub(":B0$", "", term))]
            list(k = k, p = lrt(fitm(paste(core(xw), "+ win + B0"), xw, REML = FALSE), fitm(paste(core(xw), "+ win + win:B0"), xw, REML = FALSE))) }, error = function(e) NULL)
          for (w in dsc$scope) { dd <- dsc[scope == w]; s <- show(dd$children, dd$pairs); kk <- if (is.null(r)) NULL else r$k[scope == w]
            addE(if (is.null(kk) || !nrow(kk) || dd$children < MINKIDS) erow("longitudinal", o, g, w, "BMI Z at the visit", MODW, LMM, n_carriers = s$k, n_pairs = s$pr, status = sprintf("not estimated: fewer than %d children, or model failed", MINKIDS), in_paper = paper) else
              erow("longitudinal", o, g, w, "BMI Z at the visit", MODW, LMM, kk$beta, kk$se, kk$p, n_carriers = s$k, n_pairs = s$pr, in_paper = paper)) }
          if (!is.null(r)) addT(trow("longitudinal", o, g, "All ages", "BMI x age window (likelihood ratio, ML)", MODW, "likelihood-ratio test of ML fits", NA, nlevels(xw$win) - 1, r$p, paste(levels(xw$win), collapse = "; ")))
          # sex (exploratory): BMI x sex
          if (min(x[sex == "female", uniqueN(id)], x[sex == "male", uniqueN(id)]) >= MINKIDS) addE(tryCatch({ f <- sub(" \\+ sex", "", core(x))
            a <- crow(fitm(paste(f, "+ sex + B0 + B0:sex"), x), "sexmale:B0|B0:sexmale"); bb <- crow(fitm(paste(f, "+ sex + sex:B0"), x), ":B0$")
            MODS <- sub("BMI Z at the visit [+]", "sex:BMI Z at the visit +", MODL); kf <- nkid(x[sex == "female", uniqueN(id)]); km <- nkid(x[sex == "male", uniqueN(id)]); if (anyhid) kf <- km <- paste0(">=", minc)
            rbind(erow("longitudinal sex (exploratory)", o, g, "All ages", "BMI Z at the visit", MODS, LMM, bb[grep("female", term), beta], bb[grep("female", term), se], bb[grep("female", term), p], sex = "female", n_carriers = kf, in_paper = paper),
                  erow("longitudinal sex (exploratory)", o, g, "All ages", "BMI Z at the visit", MODS, LMM, bb[grep("sexmale", term), beta], bb[grep("sexmale", term), se], bb[grep("sexmale", term), p], sex = "male", n_carriers = km, in_paper = paper),
                  erow("longitudinal sex (exploratory)", o, g, "All ages", "BMI x sex (male minus female)", sub("BMI Z at the visit [+]", "BMI Z at the visit + BMI Z:sex +", MODL), LMM, a$beta, a$se, a$p, sex = "male - female", n_carriers = s_all$k, in_paper = paper)) },
            error = function(e) NULL))
        }
        # heterogeneity across genotypes (LRT of BMI x genotype, ML), overall and within window
        het <- function(x) { if (uniqueN(x$locus) < 2) return(NA_real_); f <- sub(" \\+ locus", "", core(x))
          tryCatch(lrt(fitm(paste(f, "+ locus + B0"), x, REML = FALSE), fitm(paste(f, "+ locus * B0"), x, REML = FALSE)), error = function(e) NA_real_) }
        Ic <- I[locus %in% LOCI]
        if (POOL %in% groups && uniqueN(Ic$locus) >= 2) for (w in c("All ages", wins)) { x <- if (w == "All ages") Ic else Ic[win == w]; if (uniqueN(x$id) < MINKIDS) next; p <- het(x)
          addT(trow("longitudinal", o, POOL, w, "heterogeneity of the BMI coefficient across genotypes (likelihood ratio, ML)", sub("BMI Z at the visit [+]", "genotype * BMI Z at the visit +", MODL), "likelihood-ratio test of ML fits", NA, uniqueN(x$locus) - 1, p, paste(sort(unique(x$locus)), collapse = "; "))) }
      }
      Hl <- build_height(GROUPS5); Cl <- build_hcl(GROUPS5); LONGINFO <- rbindlist(list(Hl$log, Cl$log), fill = TRUE)
      if (!is.null(Hl)) { run_set(Hl$I, "Height", WIN_H, c(intersect(LOCI, unique(Hl$I$locus)), POOL), "yes")
        if (Hl$I[locus == CTL, uniqueN(id)] >= MINKIDS) run_set(Hl$I[locus == CTL], "Height", WIN_H, CTL, "no (new: the paper has no control with repeated visits)") }
      if (!is.null(Cl)) { run_set(Cl$I, "Head circumference", WIN_C, c(intersect(LOCI, unique(Cl$I$locus)), POOL), "yes")
        if (Cl$I[locus == CTL, uniqueN(id)] >= MINKIDS) run_set(Cl$I[locus == CTL], "Head circumference", WIN_C, CTL, "no (new: the paper has no control with repeated visits)") }
    } }

  # ---------------------------------------------------------------- e. autism diagnosis (85_v2_asd_by_cnv.R)
  ABB <- NULL
  if ("autism_diagnosis" %in% SEC && any(!is.na(D$asd))) { message("  autism diagnosis")
    PH <- list(Height = if (!is.null(DATA2$Height)) DATA2$Height$car[, .(id, locus, cohort, sex, age, y, asd)],
               `Head circumference` = if (!is.null(DATA2$`Head circumference`)) DATA2$`Head circumference`$car[, .(id, locus, cohort, sex, age, y, asd)],
               BMI = if (!is.null(DATA2$Height)) DATA2$Height$car[, .(id, locus, cohort, sex, age, y = BMI, asd)])
    PH <- lapply(PH[!vapply(PH, is.null, logical(1))], function(x) x[!is.na(asd)])
    GAMe <- "GAM, REML (mgcv gam)"; MODD <- "Z ~ autism diagnosis + s(age) + sex (+ cohort) (+ genotype)"
    for (ph in names(PH)) { hid <- FALSE; got <- list()
      for (g in c(LOCI, POOL)) { x <- if (g == POOL) PH[[ph]] else PH[[ph]][locus == g]; n1 <- sum(x$asd == 1L); n0 <- sum(x$asd == 0L); if (!nrow(x)) next
        small <- (n1 < minc) || (n0 < minc); if (small && g != POOL) hid <- TRUE
        if (small) { addE(erow("autism diagnosis", ph, g, "all ages", "Z difference, autism diagnosis minus none", MODD, GAMe, status = sprintf("suppressed: fewer than %d with or without a diagnosis", minc))); next }
        msk <- g == POOL && hid
        addE(tryCatch({ x <- copy(x)[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))]
          cv <- paste(c(if (nlevels(x$sex) > 1) "+ sex", if (nlevels(x$cohort) > 1) "+ cohort", if (nlevels(x$locus) > 1) "+ locus"), collapse = " ")
          v <- co4(gam(as.formula(paste("y ~ asd + s(age, k = 8)", cv)), data = x, method = "REML"), "asd"); if (g != POOL) got[[g]] <- v
          erow("autism diagnosis", ph, g, "all ages", "Z difference, autism diagnosis minus none", MODD, GAMe, v[1], v[2], v[3], n_carriers = if (msk) paste0(">=", minc) else nrow(x), n_asd = if (msk) paste0(">=", minc) else n1) },
          error = function(e) erow("autism diagnosis", ph, g, "all ages", "Z difference, autism diagnosis minus none", MODD, GAMe, status = paste("model failed:", conditionMessage(e)))))
      }
      if (length(got) >= 2) { q <- cochran(sapply(got, `[`, 1), sapply(got, `[`, 2))
        addT(trow("autism diagnosis", ph, POOL, "all ages", "heterogeneity of the diagnosis difference across genotypes (Cochran Q)", MODD, "Cochran Q on the per-genotype estimates", q[["Q"]], q[["df"]], q[["p"]], paste(names(got), collapse = "; "))) }
    }
    # logistic regression of diagnosis on BMI Z and on head-circumference Z; proportions by band
    PR <- list(`BMI Z` = if (!is.null(PH$BMI)) copy(PH$BMI)[, z := y], `Head-circumference Z` = if (!is.null(PH$`Head circumference`)) copy(PH$`Head circumference`)[, z := y])
    PR <- PR[!vapply(PR, is.null, logical(1))]; MODG <- "autism diagnosis ~ Z + ns(age, 3) + sex + cohort (+ genotype)"; abb <- list()
    for (pr in names(PR)) { hid <- FALSE; got <- list(); cells <- list()
      for (g in c(LOCI, POOL)) { d <- if (g == POOL) PR[[pr]] else PR[[pr]][locus == g]; if (!nrow(d)) next
        keep <- d[, .(v = sum(asd) > 0 & sum(asd) < .N), by = cohort][v == TRUE, cohort]; d2 <- d[cohort %in% keep]; n1 <- sum(d2$asd); n0 <- sum(1 - d2$asd); thr <- max(8, minc)
        if (n1 < thr || n0 < thr) { if (g != POOL && (sum(d$asd) < minc || sum(1 - d$asd) < minc)) hid <- TRUE
          addE(erow("autism diagnosis", pr, g, "all ages", "log odds ratio of autism diagnosis per SD", MODG, "logistic regression (glm)", status = sprintf("not estimated: fewer than %d with or without a diagnosis", thr)))
        } else addE(tryCatch({ tm <- c("z", "ns(age, 3)", if (uniqueN(d2$sex) > 1) "sex", if (uniqueN(d2$cohort) > 1) "cohort", if (uniqueN(d2$locus) > 1) "locus")
          m <- suppressWarnings(glm(as.formula(paste("asd ~", paste(tm, collapse = " + "))), data = d2, family = binomial)); co <- coef(summary(m))["z", ]; if (g != POOL) got[[g]] <- co[1:2]; msk <- g == POOL && hid
          erow("autism diagnosis", pr, g, "all ages", "log odds ratio of autism diagnosis per SD", MODG, "logistic regression (glm)", co[[1]], co[[2]], co[[4]], n_carriers = if (msk) paste0(">=", minc) else nrow(d2), n_asd = if (msk) paste0(">=", minc) else n1) },
          error = function(e) erow("autism diagnosis", pr, g, "all ages", "log odds ratio of autism diagnosis per SD", MODG, "logistic regression (glm)", status = paste("model failed:", conditionMessage(e)))))
        # proportion with a diagnosis by band of the Z score
        d[, band := band3(z)]; tb <- d[, .(n = .N, n_asd = sum(asd)), by = band][data.table(band = factor(BAND3, BAND3)), on = "band"]; tb[is.na(n), `:=`(n = 0L, n_asd = 0L)]
        h <- close_line(tb$n, minc) | close_line(tb$n_asd, minc) | close_line(tb$n - tb$n_asd, minc); if (sum(h) == 1) h <- h | close_line(ifelse(h, 1, tb$n), minc)
        cells[[g]] <- data.table(predictor = pr, group = g, band = BAND3, n = tb$n, n_asd = tb$n_asd, hide = h)
      }
      if (length(got) >= 2) { q <- cochran(sapply(got, `[[`, 1), sapply(got, `[[`, 2))
        addT(trow("autism diagnosis", pr, POOL, "all ages", "heterogeneity of the log odds ratio across genotypes (Cochran Q)", MODG, "Cochran Q on the per-genotype estimates", q[["Q"]], q[["df"]], q[["p"]], paste(names(got), collapse = "; "))) }
      if (length(cells)) { cc <- rbindlist(cells); anyh <- cc[group != POOL, .(a = any(hide)), by = band]; cc[group == POOL, hide := hide | anyh$a[match(band, anyh$band)] %in% TRUE]
        abb[[pr]] <- cc[, .(predictor, group, band, n = ifelse(n > 0 & n < minc, paste0("<", minc), as.character(n)), n_with_diagnosis = ifelse(hide, NA_character_, as.character(n_asd)),
                            proportion = ifelse(hide | n == 0, NA_real_, n_asd / n), status = ifelse(hide, sprintf("suppressed: a count below %d, or protects one", minc), "ok"))] }
    }
    ABB <- rbindlist(abb)
  }

  # ---------------------------------------------------------------- f. autistic traits (51b_cog_panels_boxcox_combined.R)
  if ("autistic_traits" %in% SEC && (any(!is.na(D$srs)) || any(!is.na(D$scq)))) { message("  autistic traits")
    for (ph in c("SRS", "SCQ")) { col <- tolower(ph); if (!any(!is.na(D[[col]]))) next
      P <- D[!is.na(get(col)), .(raw_total = get(col)[1], age_ph = (if (all(is.na(age_q))) age else age_q)[1]), by = .(id, locus = as.character(locus), sex, cohort)]
      V <- D[!is.na(BMI_raw), .(id, age_v = age, BMI_raw)]
      M <- merge(P, V, by = "id", allow.cartesian = TRUE)[, dage := abs(age_v - age_ph)][order(id, dage)][, .SD[1], by = id][dage <= cfg$bmi_window_years]
      if (!nrow(M)) next
      bc <- score_bmi_child(M$BMI_raw, M$age_v, M$sex); ba <- score_adult_site(M$BMI_raw, M$age_v, M$sex, M$locus == CTL, "bmi", cfg); ad <- M$age_v >= 19; bc[ad] <- ba$z[ad]
      M[, `:=`(BMI = bc, age = age_ph)]; M <- M[!is.na(BMI) & !is.na(age)]
      M[, form2 := if (ph == "SCQ") "SCQ" else ifelse(age < cfg$srs_adult_form_from_age, "SRS child form", "SRS adult form")]
      M[, y := NA_real_]
      for (f in unique(M$form2)) { L <- TRAIT_LAMBDA[[f]]; tr <- function(v) if (abs(L) < 1e-8) log(v + 1) else ((v + 1)^L - 1) / L; ref <- M[form2 == f & locus == CTL, raw_total]
        if (length(ref) < cfg$min_trait_ref) { note("%s, %s: fewer than %d controls with this form; rows of this form left out", ph, f, cfg$min_trait_ref); next }
        M[form2 == f, y := (tr(raw_total) - mean(tr(ref))) / sd(tr(ref))] }
      M <- M[!is.na(y)]; ctl <- M[locus == CTL]; if (nrow(ctl) < MINK) { note("%s: fewer than %d controls with the questionnaire and a BMI; section skipped", ph, MINK); next }
      M[, cg := ifelse(locus %in% LOCI, paste(cohort, locus), "ref")]; ctl <- M[locus == CTL]
      GAMe <- "GAM, REML (mgcv gam)"; nk <- as.character(nrow(ctl))
      addE(tryCatch({ v <- co4(gam(y ~ s(age, k = 8) + BMI, data = ctl, method = "REML"), "BMI"); erow("autistic traits", ph, CTL, "all ages", "BMI slope", "trait Z ~ s(age) + BMI Z", GAMe, v[1], v[2], v[3], n_controls = nk) }, error = function(e) NULL))
      hid <- FALSE
      for (g in c(LOCI, POOL)) { xc <- if (g == POOL) M[locus %in% LOCI] else M[locus == g]; if (!nrow(xc)) next
        if (nrow(xc) < MINC) { if (nrow(xc) < minc) hid <- TRUE; addE(erow("autistic traits", ph, g, "all ages", "CNV x BMI", "", "GAM, fREML (mgcv bam)", n_carriers = nlow(nrow(xc), MINC), n_controls = nk, status = sprintf("not estimated: fewer than %d carriers", MINC))); next }
        ex <- if (uniqueN(xc$cg) > 1) "+ cg" else ""; nc <- if (g == POOL && hid) paste0(">=", minc) else as.character(nrow(xc))
        addE(tryCatch({ v <- co4(gam(as.formula(paste("y ~ s(age, k = 8) + BMI", ex)), data = xc, method = "REML"), "BMI")
          erow("autistic traits", ph, g, "all ages", "BMI slope, carriers", "trait Z ~ s(age) + BMI Z (+ cohort x genotype)", GAMe, v[1], v[2], v[3], n_carriers = nc) }, error = function(e) NULL))
        addE(tryCatch({ z <- rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == CTL, "ctl", "CNV"), levels = c("ctl", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
          v <- co4(bam(as.formula(paste("y ~ is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI", ex)), data = z, method = "fREML"), "cnvBMI")
          erow("autistic traits", ph, g, "all ages", "CNV x BMI", "trait Z ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI Z) + CNV x BMI Z (+ cohort x genotype)", "GAM, fREML (mgcv bam)", v[1], v[2], v[3], n_carriers = nc, n_controls = nk) },
          error = function(e) erow("autistic traits", ph, g, "all ages", "CNV x BMI", "", "GAM, fREML (mgcv bam)", status = paste("model failed:", conditionMessage(e)))))
        # ordinal check: cumulative-logit model of the raw total, per form, inverse-variance combined (ord1 / ord_p of 51b)
        if (requireNamespace("ordinal", quietly = TRUE)) addE(tryCatch({
          r <- lapply(split(rbind(xc, ctl), by = "form2"), function(zz) { if (sum(zz$locus != CTL) < MINC) return(NULL)
            tryCatch({ zz <- copy(zz); br <- unique(quantile(zz$raw_total, seq(0, 1, length.out = 41), type = 1)); zz[, ry := factor(findInterval(raw_total, br, rightmost.closed = TRUE), ordered = TRUE)]
              zz[, is_cnv := as.numeric(locus != CTL)][, cnvBMI := is_cnv * BMI]; cl1 <- if (uniqueN(zz$locus[zz$locus != CTL]) > 1) "locus + " else ""
              m <- ordinal::clm(as.formula(paste("ry ~", cl1, "is_cnv * poly(age, 2) + BMI + cnvBMI")), data = zz); co <- coef(summary(m)); c(co["cnvBMI", 1], co["cnvBMI", 2]) }, error = function(e) NULL) })
          r <- do.call(rbind, r[!vapply(r, is.null, logical(1))]); if (is.null(r) || !nrow(r)) stop("no form with enough carriers")
          r <- r[is.finite(r[, 1]) & is.finite(r[, 2]), , drop = FALSE]; f <- fe_meta(r[, 1], r[, 2])
          erow("autistic traits", ph, g, "all ages", "CNV x BMI (ordinal check, log odds)", "ordered raw total (<= 40 levels) ~ (genotype +) is_cnv * poly(age, 2) + BMI Z + CNV x BMI Z; per form, inverse-variance combined", "cumulative logit (ordinal clm)", f$est, f$se, f$p, n_carriers = nc, n_controls = nk) },
          error = function(e) NULL))
      }
    } }
  SCORING_OF <- function(an) ifelse(grepl("^longitudinal", an), SCORING_CHILD, ifelse(an == "autistic traits", SCORING_TRAIT, SG))

  # ---------------------------------------------------------------- g. federated trajectories
  BINS <- list(); SPL <- list(); MAT <- list()
  if ("trajectories" %in% SEC) for (o in names(DATA2)) for (G in GROUPS5) { x <- (if (G == CTL) DATA2[[o]]$ctl else DATA2[[o]]$car[locus == G])[age >= 2 & age <= 73]
    x <- copy(x)[, `:=`(tl = tail_of(BMI), ab = age_bin(age))]
    for (tname in unname(TAILS)) { xt <- x[tl == tname]; n <- nrow(xt)
      bt <- xt[, .(n = .N, mean_z = mean(y), sd_z = sd(y), mean_age = mean(age)), by = ab][data.table(ab = factor(AGE_BIN_LABELS, AGE_BIN_LABELS)), on = "ab"]; bt[is.na(n), n := 0L]
      hide <- bt$n > 0 & bt$n < minc; mask <- close_line(bt$n, minc) & !hide; tot_hidden <- n > 0 & n < minc
      BINS[[length(BINS) + 1]] <- data.table(outcome = o, group = G, bmi_tail = tname, age_bin = AGE_BIN_LABELS, age_bin_mid = AGE_BIN_MID,
        n = ifelse(hide, paste0("<", minc), ifelse(mask, paste0(">=", minc), as.character(bt$n))), mean_z = ifelse(hide | bt$n == 0, NA_real_, bt$mean_z),
        se = ifelse(hide | bt$n < 2, NA_real_, bt$sd_z / sqrt(bt$n)), mean_age = ifelse(hide | bt$n == 0, NA_real_, round(bt$mean_age, 1)),
        status = ifelse(bt$n == 0, "empty", ifelse(hide, sprintf("suppressed: fewer than %d", minc), ifelse(mask, "ok (count masked to protect a small cell)", "ok"))))
      cid <- paste(SITE, o, G, tname, sep = " | "); base <- data.table(curve_id = cid, outcome = o, group = G, bmi_tail = tname, basis = SPLINE_ID,
        knots_age = paste(SPLINE_KNOTS_AGE, collapse = ";"), boundary_knots_age = paste(SPLINE_BOUNDARY_AGE, collapse = ";"))
      if (n < max(cfg$min_curve_n, length(SPLINE_NAMES) + minc)) { SPL[[length(SPL) + 1]] <- cbind(base, n = if (tot_hidden) paste0("<", minc) else as.character(n), status = sprintf("no curve: fewer than %d persons in this tail", max(cfg$min_curve_n, length(SPLINE_NAMES) + minc))); next }
      X <- spline_basis(xt$age); XtX <- crossprod(X); Xty <- crossprod(X, xt$y); pin <- psd_inv(XtX); beta <- drop(pin$inv %*% Xty); rss <- sum((xt$y - drop(X %*% beta))^2); dfr <- n - pin$rank; s2 <- rss / dfr
      ag <- sort(xt$age); k <- minc
      SPL[[length(SPL) + 1]] <- cbind(base, n = as.character(n), rank = pin$rank, df_residual = dfr, sigma2 = s2, rss = rss, age_lo = floor(2 * ag[k]) / 2, age_hi = ceiling(2 * ag[n - k + 1]) / 2, status = "ok")
      V <- s2 * pin$inv; ij <- expand.grid(i = seq_along(SPLINE_NAMES), j = seq_along(SPLINE_NAMES))
      MAT[[length(MAT) + 1]] <- rbind(data.table(curve_id = cid, quantity = "coefficient", row = SPLINE_NAMES, col = "", value = beta),
                                      data.table(curve_id = cid, quantity = "information_Xty", row = SPLINE_NAMES, col = "", value = drop(Xty)),
                                      data.table(curve_id = cid, quantity = "covariance", row = SPLINE_NAMES[ij$i], col = SPLINE_NAMES[ij$j], value = V[cbind(ij$i, ij$j)]),
                                      data.table(curve_id = cid, quantity = "information_XtX", row = SPLINE_NAMES[ij$i], col = SPLINE_NAMES[ij$j], value = XtX[cbind(ij$i, ij$j)]))
    } }

  # ---------------------------------------------------------------- assemble, label, write
  lab <- function(d, scoring_label, est_label) { if (is.null(d) || !nrow(d)) return(NULL); d <- copy(d)
    set(d, j = "estimator", value = est_label); set(d, j = "scoring", value = scoring_label)
    cbind(data.table(site = SITE, script_version = FED_VERSION, run_date = as.character(Sys.Date())), d) }
  E <- rbindlist(EST, fill = TRUE); if (nrow(E)) E[, scoring := SCORING_OF(analysis)]
  Tt <- rbindlist(TST, fill = TRUE); if (nrow(Tt)) Tt[, scoring := SCORING_OF(analysis)]
  hdr <- function(d) if (is.null(d) || !nrow(d)) NULL else cbind(data.table(site = SITE, script_version = FED_VERSION, run_date = as.character(Sys.Date())), d)
  res <- list(
    sample_n = lab(SMP, SG, "count of persons (one visit per person, cross-sectional sample of the outcome)"),
    estimates = hdr(E), tests = hdr(Tt),
    autism_by_bmi_band = lab(ABB, SG, "proportion with an autism diagnosis by band of the Z score"),
    trajectory_bins = lab(rbindlist(BINS), SG, "mean Z in a fixed age bin (one visit per person)"),
    trajectory_spline = lab(rbindlist(SPL, fill = TRUE), SG, "unpenalised least squares on a fixed natural cubic spline basis of log(age)"),
    trajectory_spline_matrices = lab(rbindlist(MAT), SG, "unpenalised least squares on a fixed natural cubic spline basis of log(age)"),
    scoring_info = lab(rbind(SCI, if (!is.null(LONGINFO) && nrow(LONGINFO)) data.table(table = paste(LONGINFO$outcome, "longitudinal"), trait = "cleaning", note = apply(LONGINFO[, -1], 1, function(r) paste(names(LONGINFO)[-1], trimws(r), sep = " = ", collapse = "; "))), fill = TRUE), SG, "scoring and cleaning summary"))
  if (!is.null(res$scoring_info)) { res$scoring_info$n_reference_adults <- fmt_n(res$scoring_info$n_reference_adults, minc); res$scoring_info$n_scored <- fmt_n(res$scoring_info$n_scored, minc) }
  ri <- c(site = SITE, script_version = FED_VERSION, run_date = as.character(Sys.Date()), R = R.version.string, mgcv = as.character(packageVersion("mgcv")),
          visit_rule = cfg$visit_rule, visit_seed = cfg$visit_seed, adult_boxcox = cfg$adult_boxcox, min_cell = minc, min_model_carriers = MINC, min_model_controls = MINK,
          min_longitudinal_children = MINKIDS, min_curve_n = cfg$min_curve_n, min_adult_ref = cfg$min_adult_ref, spline_basis = SPLINE_ID, age_bins = paste(AGE_BIN_LABELS, collapse = ";"),
          family_id_given = HAVE_FAM, cohort_column_given = INFO$has_cohort, autism_diagnosis_given = any(!is.na(D$asd)), srs_given = any(!is.na(D$srs)), scq_given = any(!is.na(D$scq)),
          n_rows_input = INFO$n_rows_input, n_rows_used = INFO$n_rows_used, n_rows_genotype_not_in_map = INFO$n_rows_genotype_not_in_map, n_rows_missing_sex_or_age = INFO$n_rows_missing_sex_or_age,
          n_persons_inconsistent_sex_or_genotype = INFO$n_persons_inconsistent_sex_or_genotype,
          setNames(INFO$n_values_outside_plausible_range, paste0("n_values_outside_plausible_range_", names(INFO$n_values_outside_plausible_range))),
          sections = paste(SEC, collapse = ";"), if (length(NOTES)) setNames(NOTES, paste0("note_", seq_along(NOTES))))
  res$run_info <- data.table(site = SITE, item = names(ri), value = as.character(ri))
  if (write) { dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(res)) if (!is.null(res[[nm]])) fwrite(res[[nm]], file.path(cfg$out_dir, paste0(nm, ".csv")))
    message("wrote aggregate tables to ", normalizePath(cfg$out_dir), "  (look at them before sending; nothing else should leave the site)") }
  if (keep_data) attr(res, "analysis_data") <- DATA2
  invisible(res)
}

if (!exists("FEDERATED_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); opt <- function(k) { i <- match(k, a); if (is.na(i)) NULL else a[i + 1] }
  if (is.null(opt("--config")) && is.null(opt("--input"))) stop("usage: Rscript site_export.R --config my_config.json [--input F] [--out DIR] [--site LABEL] [--min_cell N]")
  cfg <- if (!is.null(opt("--config"))) fromJSON(opt("--config"), simplifyVector = TRUE) else list()
  if (!is.null(opt("--input"))) cfg$input_file <- opt("--input")
  if (!is.null(opt("--out"))) cfg$out_dir <- opt("--out")
  if (!is.null(opt("--site"))) cfg$site_label <- opt("--site")
  if (!is.null(opt("--min_cell"))) cfg$min_cell <- as.numeric(opt("--min_cell"))
  run_site_export(cfg)
}
