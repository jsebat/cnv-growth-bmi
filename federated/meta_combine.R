# meta_combine.R  (2026-11-20)
# What: combines the aggregate tables written by site_export.R at several sites (the lab's own export is one site).
#   1. Every scalar estimate: fixed-effects inverse-variance meta-analysis across sites, with Cochran Q and I2.
#      Sites are strata: nothing is re-estimated, and a suppressed or not-estimated cell is never pooled.
#   2. Heterogeneity across the four genotypes: Cochran Q on the pooled per-genotype estimates.
#   3. Sample sizes and autism-by-BMI-band counts: sums over the sites that show the cell.
#   4. Federated trajectories. Binned means: n-weighted mean of the site means (suppressed or masked bins left out).
#      Curves: every site fitted the same fixed spline basis; the coefficient vectors are pooled by multivariate
#      inverse-variance weighting (generalised least squares with each site's covariance),
#          beta = (sum_i P_i)^-1 sum_i P_i beta_i ,  P_i = XtX_i / sigma2_i  (= inverse of the site covariance),
#      and the pooled curve B(age) beta gets a pointwise 95% band from (sum_i P_i)^-1. A site whose ages do not cover
#      the whole basis has a rank-deficient P_i and simply adds no information there. The pooled curve is drawn only
#      at ages inside the supported range of at least one site. --trajectory_weights common_sigma uses P_i = XtX_i
#      (one residual variance for all sites), which reproduces the least-squares fit to the stacked data exactly.
# Inputs: two or more export folders (--in), each with estimates.csv, tests.csv, sample_n.csv, trajectory_*.csv ...
# Outputs (in --out): meta_estimates.csv, all_estimates.csv (site rows + pooled rows, for the forests),
#   meta_heterogeneity_across_genotypes.csv, site_tests.csv, meta_sample_n.csv, meta_autism_by_bmi_band.csv,
#   meta_trajectory_bins.csv, meta_trajectory_curves.csv, meta_trajectory_coefficients.csv, meta_info.csv.
# null model: none; pooled estimate tested by a z test (estimate / SE); heterogeneity across sites and across
#   genotypes by Cochran Q (chi-square).
# Usage: Rscript meta_combine.R --in siteA_out siteB_out [...] --out meta_out [--trajectory_weights site_covariance|common_sigma]
suppressPackageStartupMessages({ library(data.table); library(splines) })
.fed_dir <- function() { if (exists("FEDERATED_DIR")) return(FEDERATED_DIR)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else getwd() }
source(file.path(.fed_dir(), "00_federated_common.R"))
POOLED <- "Pooled (fixed effects)"

rd <- function(dirs, f) rbindlist(lapply(dirs, function(d) { p <- file.path(d, f); if (!file.exists(p)) return(NULL)
  fread(p, colClasses = list(character = intersect(c("n", "n_carriers", "n_controls", "n_pairs", "n_asd", "n_with_diagnosis", "df", "value_chr"), names(fread(p, nrows = 0)))), na.strings = c("", "NA")) }), fill = TRUE)
numn <- function(x) suppressWarnings(as.numeric(x))   # "<5" and ">=5" become NA: such cells are not summed
sumn <- function(x) { v <- numn(x); if (all(is.na(x))) NA_character_ else if (anyNA(v[!is.na(x)])) paste0(">=", sum(v, na.rm = TRUE)) else as.character(sum(v)) }

meta_combine <- function(dirs, out = NULL, trajectory_weights = "site_covariance") {
  stopifnot(trajectory_weights %in% c("site_covariance", "common_sigma")); res <- list()
  # ---------------------------------------------------------------- 1. scalar estimates
  E <- rd(dirs, "estimates.csv"); if (!nrow(E)) stop("no estimates.csv found in the input folders")
  if (anyDuplicated(unique(E[, .(site, script_version)])$site)) warning("a site appears with more than one script version")
  if (uniqueN(E$script_version) > 1) warning("sites ran different script versions: ", paste(unique(E$script_version), collapse = " | "))
  KEY <- c("analysis", "outcome", "group", "stratum", "sex", "term", "model")
  E[, use := status == "ok" & is.finite(estimate) & is.finite(se) & se > 0]
  M <- E[, { u <- which(use); f <- fe_meta(estimate[u], se[u])
    cl <- if (length(u) && all(is.finite(se_family_clustered[u]))) fe_meta(estimate[u], se_family_clustered[u]) else fe_meta(numeric(0), numeric(0))
    .(n_sites_pooled = f$k, sites_pooled = paste(site[u], collapse = "; "), sites_not_pooled = paste(site[!use], collapse = "; "),
      estimate = f$est, se = f$se, z = f$z, p = f$p, ci_lo = f$est - 1.96 * f$se, ci_hi = f$est + 1.96 * f$se, Q = f$Q, Q_df = f$Q_df, Q_p = f$Q_p, I2 = f$I2,
      estimate_family_clustered_weights = cl$est, se_family_clustered = cl$se, p_family_clustered = cl$p,
      n_carriers = sumn(n_carriers[u]), n_controls = sumn(n_controls[u]), n_pairs = sumn(n_pairs[u]), n_asd = sumn(n_asd[u]),
      estimator = paste(unique(estimator[u]), collapse = " | "), scoring = paste(unique(scoring[u]), collapse = " | "),
      same_model_as_paper = paste(unique(same_model_as_paper), collapse = " | ")) }, by = KEY]
  M[, meta_method := "fixed effects, inverse-variance weights; Cochran Q and I2 across sites"]
  res$meta_estimates <- M
  res$all_estimates <- rbind(E[, c("site", KEY, "estimate", "se", "p", "n_carriers", "n_controls", "n_pairs", "n_asd", "status", "estimator"), with = FALSE],
    M[n_sites_pooled > 0, c(list(site = POOLED), .SD, list(status = "ok")), .SDcols = c(KEY, "estimate", "se", "p", "n_carriers", "n_controls", "n_pairs", "n_asd")][, estimator := "fixed-effects meta-analysis"], fill = TRUE)
  # ---------------------------------------------------------------- 2. heterogeneity across genotypes (pooled estimates)
  res$meta_heterogeneity_across_genotypes <- M[group %in% LOCI & n_sites_pooled > 0, { f <- fe_meta(estimate, se)
    .(n_genotypes = f$k, genotypes = paste(group, collapse = "; "), Q = f$Q, Q_df = f$Q_df, Q_p = f$Q_p, I2 = f$I2, mean_of_genotypes = f$est, se = f$se) }, by = setdiff(KEY, "group")][n_genotypes >= 2]
  Tt <- rd(dirs, "tests.csv"); if (nrow(Tt)) res$site_tests <- Tt
  # ---------------------------------------------------------------- 3. counts
  S <- rd(dirs, "sample_n.csv")
  if (nrow(S)) res$meta_sample_n <- S[, .(n_sites = .N, n_sites_shown = sum(!is.na(numn(n))), n_sum_of_shown = sum(numn(n), na.rm = TRUE),
    note = if (anyNA(numn(n))) "at least this many: one or more sites suppressed or masked the cell" else ""), by = .(outcome, group, age_band, sex)]
  A <- rd(dirs, "autism_by_bmi_band.csv")
  if (nrow(A)) res$meta_autism_by_bmi_band <- A[, { u <- status == "ok" & !is.na(numn(n)) & !is.na(numn(n_with_diagnosis)); n1 <- sum(numn(n[u])); k1 <- sum(numn(n_with_diagnosis[u]))
    .(n_sites_pooled = sum(u), n_sites_suppressed = sum(!u), n = n1, n_with_diagnosis = k1, proportion = if (n1 > 0) k1 / n1 else NA_real_) }, by = .(predictor, group, band)]
  # ---------------------------------------------------------------- 4. trajectories
  B <- rd(dirs, "trajectory_bins.csv")
  if (nrow(B)) { B[, nn := numn(n)]
    pb <- B[, { u <- status == "ok" & !is.na(nn) & is.finite(mean_z); N <- sum(nn[u])
      .(n_sites_pooled = sum(u), n_sites_left_out = sum(status != "empty" & !u), n = N, mean_z = if (N > 0) sum(nn[u] * mean_z[u]) / N else NA_real_,
        se = if (N > 0 && all(is.finite(se[u]))) sqrt(sum(nn[u]^2 * se[u]^2)) / N else NA_real_, mean_age = if (N > 0) sum(nn[u] * mean_age[u]) / N else NA_real_) }, by = .(outcome, group, bmi_tail, age_bin, age_bin_mid)]
    res$meta_trajectory_bins <- rbind(B[status == "ok" & !is.na(nn), .(scope = site, outcome, group, bmi_tail, age_bin, age_bin_mid, n = nn, mean_z, se, mean_age)],
                                      pb[n > 0, .(scope = POOLED, outcome, group, bmi_tail, age_bin, age_bin_mid, n, mean_z, se, mean_age, n_sites_pooled, n_sites_left_out)], fill = TRUE) }
  C <- rd(dirs, "trajectory_spline.csv"); MT <- rd(dirs, "trajectory_spline_matrices.csv")
  if (nrow(C) && nrow(MT)) {
    C <- C[status == "ok"]; p <- length(SPLINE_NAMES)
    if (any(C$basis != SPLINE_ID)) stop("a site used a different spline basis; all sites must run the same package version.\n  expected: ", SPLINE_ID, "\n  found: ", paste(setdiff(unique(C$basis), SPLINE_ID), collapse = " | "))
    getm <- function(cid, q) { m <- MT[curve_id == cid & quantity == q]; if (q %in% c("coefficient", "information_Xty")) return(setNames(m$value, m$row)[SPLINE_NAMES])
      A <- matrix(NA_real_, p, p, dimnames = list(SPLINE_NAMES, SPLINE_NAMES)); A[cbind(m$row, m$col)] <- m$value; A }
    Bg <- spline_basis(CURVE_GRID); curves <- list(); coefs <- list()
    curve_rows <- function(scope, o, G, tl, beta, V, sup, n, k) { ok <- sup
      data.table(scope = scope, outcome = o, group = G, bmi_tail = tl, age = CURVE_GRID[ok], fit = drop(Bg[ok, , drop = FALSE] %*% beta),
                 se = sqrt(pmax(0, rowSums((Bg[ok, , drop = FALSE] %*% V) * Bg[ok, , drop = FALSE]))), n = n, n_sites = k) }
    for (kk in split(C, by = c("outcome", "group", "bmi_tail"), drop = TRUE)) { o <- kk$outcome[1]; G <- kk$group[1]; tl <- kk$bmi_tail[1]
      P <- matrix(0, p, p); rhs <- rep(0, p); st <- list(); anysup <- rep(FALSE, length(CURVE_GRID))
      for (i in seq_len(nrow(kk))) { cid <- kk$curve_id[i]; b <- getm(cid, "coefficient"); X2 <- getm(cid, "information_XtX"); V <- getm(cid, "covariance")
        sup <- CURVE_GRID >= kk$age_lo[i] & CURVE_GRID <= kk$age_hi[i]; anysup <- anysup | sup
        curves[[length(curves) + 1]] <- curve_rows(kk$site[i], o, G, tl, b, V, sup, as.numeric(kk$n[i]), 1L)
        w <- if (trajectory_weights == "site_covariance") 1 / kk$sigma2[i] else 1; xy <- getm(cid, "information_Xty"); if (anyNA(xy)) xy <- drop(X2 %*% b); P <- P + w * X2; rhs <- rhs + w * xy; st[[i]] <- list(b = b, X2 = X2) }
      pin <- psd_inv(P); beta <- drop(pin$inv %*% rhs); names(beta) <- SPLINE_NAMES; N <- sum(as.numeric(kk$n))
      V <- if (trajectory_weights == "site_covariance") pin$inv else {
        rss <- sum(kk$rss) + sum(vapply(st, function(s) drop(t(s$b - beta) %*% s$X2 %*% (s$b - beta)), numeric(1))); (rss / (N - pin$rank)) * pin$inv }
      curves[[length(curves) + 1]] <- curve_rows(POOLED, o, G, tl, beta, V, anysup, N, nrow(kk))
      ij <- expand.grid(i = seq_len(p), j = seq_len(p))
      coefs[[length(coefs) + 1]] <- rbind(data.table(outcome = o, group = G, bmi_tail = tl, quantity = "coefficient", row = SPLINE_NAMES, col = "", value = beta),
                                          data.table(outcome = o, group = G, bmi_tail = tl, quantity = "covariance", row = SPLINE_NAMES[ij$i], col = SPLINE_NAMES[ij$j], value = V[cbind(ij$i, ij$j)]))
    }
    CV <- rbindlist(curves); CV[, `:=`(lo = fit - 1.96 * se, hi = fit + 1.96 * se)]
    # pooled high-minus-low BMI tail difference (independent samples), at ages supported in both tails
    w <- dcast(CV[scope == POOLED], outcome + group + age ~ bmi_tail, value.var = c("fit", "se"))
    hn <- paste0("fit_", TAILS[["high"]]); ln <- paste0("fit_", TAILS[["low"]]); hs <- paste0("se_", TAILS[["high"]]); ls <- paste0("se_", TAILS[["low"]])
    if (all(c(hn, ln) %in% names(w))) { w <- w[!is.na(get(hn)) & !is.na(get(ln))]
      CV <- rbind(CV, w[, .(scope = POOLED, outcome, group, bmi_tail = "high minus low BMI tail", age, fit = get(hn) - get(ln), se = sqrt(get(hs)^2 + get(ls)^2))][, `:=`(lo = fit - 1.96 * se, hi = fit + 1.96 * se)], fill = TRUE) }
    CV[, `:=`(basis = SPLINE_ID, pooling = if (trajectory_weights == "site_covariance") "multivariate inverse-variance (GLS with each site's covariance)" else "one residual variance (equals least squares on the stacked data)")]
    res$meta_trajectory_curves <- CV; res$meta_trajectory_coefficients <- rbindlist(coefs)
  }
  res$meta_info <- data.table(item = c("script_version", "run_date", "n_sites", "sites", "site_script_versions", "trajectory_weights", "spline_basis"),
    value = c(FED_VERSION, as.character(Sys.Date()), uniqueN(E$site), paste(unique(E$site), collapse = "; "), paste(unique(E$script_version), collapse = " | "), trajectory_weights, SPLINE_ID))
  if (!is.null(out)) { dir.create(out, recursive = TRUE, showWarnings = FALSE); for (nm in names(res)) fwrite(res[[nm]], file.path(out, paste0(nm, ".csv"))); message("wrote ", normalizePath(out)) }
  invisible(res)
}

if (!exists("FEDERATED_SOURCE_ONLY") && !interactive()) {
  a <- commandArgs(TRUE); i <- match("--in", a); o <- match("--out", a); tw <- match("--trajectory_weights", a)
  if (is.na(i) || is.na(o)) stop("usage: Rscript meta_combine.R --in siteA_out siteB_out [...] --out meta_out [--trajectory_weights site_covariance|common_sigma]")
  flags <- which(startsWith(a, "--")); end <- min(c(flags[flags > i], length(a) + 1)) - 1
  meta_combine(a[(i + 1):end], a[o + 1], if (is.na(tw)) "site_covariance" else a[tw + 1])
}
