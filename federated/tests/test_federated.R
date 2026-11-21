# test_federated.R  (2026-11-20)
# What: end-to-end test of the federated package on two SYNTHETIC sites (tests/make_synthetic.R): generates the data,
#   runs site_export.R for each site through its command line, then meta_combine.R and plot_federated.R, and checks
#   (1) the simulated effects are recovered within sampling error (CNV x BMI on height at 2-9 y in one genotype, the
#       next-visit effect on head circumference, the age-varying height trajectory, the diagnosis and trait effects);
#   (2) the pooled spline curve from the two sites equals the curve fitted to the two sites' data stacked on the same
#       basis (least squares with one residual variance; weighted least squares with each site's variance);
#   (3) small cells are suppressed and cannot be recovered by subtraction; (4) no output contains an id;
#   (5) the meta-analysis equals a hand calculation; (6) every table is labelled; (7) the figures are written.
# Inputs: none. Outputs: PASS / FAIL lines; everything is written to a temporary folder outside the repository
#   (or to the folder given as first argument).
# null model: none; "recovered" = |estimate - simulated effect| < 4 SE + a small allowance stated in each check.
# Usage: Rscript tests/test_federated.R [work_dir]
a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)][1]), fixed = TRUE)
TESTS <- dirname(normalizePath(f)); MOD <- dirname(TESTS)
work <- commandArgs(TRUE)[1]; if (is.na(work)) work <- file.path(tempdir(), "federated_test")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
if (startsWith(normalizePath(work), normalizePath(MOD))) stop("choose a work folder outside the repository")
RS <- file.path(R.home("bin"), "Rscript"); nfail <- 0
check <- function(ok, what) { ok <- isTRUE(ok); cat(if (ok) "PASS " else "FAIL ", what, "\n"); if (!ok) nfail <<- nfail + 1; invisible(ok) }
run <- function(...) { st <- system2(RS, shQuote(c(...))); if (st != 0) stop("command failed: ", paste(..., collapse = " ")) }

FEDERATED_SOURCE_ONLY <- TRUE; FEDERATED_DIR <- MOD
source(file.path(TESTS, "make_synthetic.R")); source(file.path(MOD, "site_export.R")); source(file.path(MOD, "meta_combine.R"))
syn <- make_all(work, file.path(MOD, "growth_refs_lms.csv")); TR <- syn$truth; FL <- syn$files

# ---- unit checks
n <- rbind(c(3, 12, 15), c(20, 30, 50), c(23, 42, 65)); cl <- close_table(n, 5)
check(cl$hidden[1, 1] && all(rowSums(cl$hidden | cl$masked)[1:2] >= 2) && all(colSums(cl$hidden | cl$masked)[1:2] >= 2), "close_table(): a hidden cell cannot be recovered from its row or column")
f1 <- fe_meta(c(0.2, 0.4), c(0.1, 0.2)); check(abs(f1$est - (0.2 / 0.01 + 0.4 / 0.04) / (1 / 0.01 + 1 / 0.04)) < 1e-12 && abs(f1$Q - ((0.2 - f1$est)^2 / 0.01 + (0.4 - f1$est)^2 / 0.04)) < 1e-12, "fe_meta(): inverse-variance mean and Cochran Q")

# ---- pipeline through the command line
for (s in names(FL)) run(file.path(MOD, "site_export.R"), "--config", FL[[s]]$config)
dirs <- vapply(FL, `[[`, "", "out"); MO <- file.path(work, "meta_out"); FG <- file.path(work, "figures")
run(file.path(MOD, "meta_combine.R"), "--in", dirs, "--out", MO)
run(file.path(MOD, "plot_federated.R"), "--meta", MO, "--out", FG)
M <- fread(file.path(MO, "meta_estimates.csv"), na.strings = c("", "NA")); AE <- fread(file.path(MO, "all_estimates.csv"), na.strings = c("", "NA"))
get1 <- function(an, o, g, st, tm, sx = "all", mod = NULL) { r <- M[analysis == an & outcome == o & group == g & stratum == st & term == tm & sex == sx]; if (!is.null(mod)) r <- r[grepl(mod, model, fixed = TRUE)]; stopifnot(nrow(r) == 1); r }
near <- function(r, truth, allow = 0) r$n_sites_pooled >= 1 && abs(r$estimate - truth) < 4 * r$se + allow
show <- function(r, truth) sprintf("(estimate %.3f, SE %.3f, simulated %.2f, %d site%s)", r$estimate, r$se, truth, r$n_sites_pooled, if (r$n_sites_pooled > 1) "s" else "")

# ---- (1) simulated effects recovered
r <- get1("cross-sectional band", "Height", "16p11.2 del", "2-9", "CNV x BMI"); check(near(r, TR$cnvxbmi_height_16pdel_2_9) && r$n_sites_pooled == 2 && r$p < 1e-4, paste("CNV x BMI on height, 16p11.2 del, 2-9 y", show(r, TR$cnvxbmi_height_16pdel_2_9)))
r <- get1("cross-sectional band", "Height", "16p11.2 del", "10-18", "CNV x BMI"); check(near(r, 0), paste("no CNV x BMI on height, 16p11.2 del, 10-18 y", show(r, 0)))
r <- get1("cross-sectional band", "Height", "22q11.2 del", "2-9", "CNV x BMI"); check(near(r, 0), paste("no CNV x BMI on height, 22q11.2 del, 2-9 y", show(r, 0)))
r <- get1("cross-sectional band", "Height", "No CNV", "2-9", "BMI slope (controls-only model)"); check(near(r, TR$ctl_slope), paste("control BMI slope on height, 2-9 y", show(r, 0.15)))
r <- get1("cross-sectional band", "Height", "All four CNVs", "2-9", "CNV x BMI", mod = "(primary)"); check(r$estimate > 0.03 && r$estimate < 0.36, sprintf("pooled four-CNV CNV x BMI at 2-9 y lies between 0 and the 16p11.2 del effect (%.3f)", r$estimate))
r <- get1("cross-sectional band", "Head circumference", "16p11.2 del", "2-9", "CNV x BMI"); check(near(r, 0), paste("no CNV x BMI on head circumference, 16p11.2 del, 2-9 y", show(r, 0)))
r <- get1("cross-sectional all ages", "Head circumference", "16p11.2 dup", "all ages", "CNV x BMI"); check(near(r, 0), paste("all-ages CNV x BMI (bam), head circumference, 16p11.2 dup", show(r, 0)))
r <- get1("sex", "Height", "16p11.2 del", "2-9", "CNV x BMI x sex (male minus female)", sx = "male - female"); check(near(r, 0), paste("no CNV x BMI x sex, height, 16p11.2 del, 2-9 y", show(r, 0)))
r <- get1("longitudinal", "Head circumference", "22q11.2 dup", "All ages", "BMI Z at the visit"); check(near(r, TR$next_visit_hc_22qdup, 0.01) && r$p < 1e-4, paste("next-visit effect on head circumference, 22q11.2 dup", show(r, 0.12)))
r <- get1("longitudinal", "Head circumference", "16p11.2 del", "All ages", "BMI Z at the visit"); check(near(r, 0, 0.01), paste("no next-visit effect on head circumference, 16p11.2 del", show(r, 0)))
r <- get1("longitudinal", "Height", "22q11.2 del", "All ages", "BMI Z at the visit"); check(near(r, 0, 0.02), paste("no next-visit effect on height, 22q11.2 del", show(r, 0)))
r <- get1("autism diagnosis", "BMI Z", "All four CNVs", "all ages", "log odds ratio of autism diagnosis per SD"); check(near(r, TR$asd_logor_bmi), paste("diagnosis log odds per SD of BMI, pooled genotypes", show(r, 0.4)))
r <- get1("autistic traits", "SRS", "All four CNVs", "all ages", "CNV x BMI"); check(near(r, TR$trait_cnvxbmi, 0.03) && r$n_sites_pooled == 2, paste("CNV x BMI on SRS (Box-Cox scale)", show(r, 0.15)))
r <- get1("autistic traits", "SRS", "All four CNVs", "all ages", "CNV x BMI (ordinal check, log odds)"); check(r$estimate > 0 && r$p < 0.01, sprintf("ordinal check agrees in sign (log odds %.3f, p %.1e)", r$estimate, r$p))
r <- get1("autistic traits", "SCQ", "All four CNVs", "all ages", "CNV x BMI"); check(near(r, TR$trait_cnvxbmi, 0.03) && r$n_sites_pooled == 1, paste("CNV x BMI on SCQ, given by one site only", show(r, 0.15)))

# ---- (1c, 2) trajectories
CV <- fread(file.path(MO, "meta_trajectory_curves.csv")); eb <- dnorm(qnorm(0.7)) / 0.3   # mean BMI Z in a 30% tail of N(0, 1)
tr <- CV[scope == "Pooled (fixed effects)" & outcome == "Height" & group == "22q11.2 del" & age %in% c(4, 8, 14, 30, 55)]
tr[, truth := truth_height_22qdel(age) + TR$ctl_slope * ifelse(bmi_tail == TAILS[["high"]], eb, -eb)]; tr <- tr[bmi_tail %in% TAILS]
check(nrow(tr) == 10 && all(abs(tr$fit - tr$truth) < 4 * tr$se + 0.08), sprintf("age-varying height trajectory of 22q11.2 del recovered at 4, 8, 14, 30, 55 y (max |curve - simulated| = %.3f)", max(abs(tr$fit - tr$truth)))); print(tr[, .(bmi_tail, age, fit = round(fit, 3), se = round(se, 3), truth = round(truth, 3))])
check(tr[age == 4 & bmi_tail == TAILS[["high"]], fit] < tr[age == 30 & bmi_tail == TAILS[["high"]], fit] - 0.2, "pooled curve rises from 4 y to 30 y as simulated")
sb <- CV[scope == "SYNTHETIC_B" & outcome == "Height"]; check(nrow(sb) > 0 && max(sb$age) < 19.6, sprintf("the paediatric site's own curves stop where its data stop (max age %.1f)", max(sb$age)))
check(CV[scope == "Pooled (fixed effects)" & outcome == "Height" & group == "22q11.2 del", max(age)] > 60, "the pooled curve covers the ages supported by at least one site")
# stacked data, same basis: analysis tables rebuilt in memory with the same config (nothing is written)
AD <- lapply(FL, function(x) attr(run_site_export(modifyList(x$cfg, list(sections = "trajectories")), write = FALSE, keep_data = TRUE), "analysis_data"))
SPc <- rbindlist(lapply(dirs, function(d) fread(file.path(d, "trajectory_spline.csv"), colClasses = list(character = "n"))), fill = TRUE)[status == "ok"]
mcs <- meta_combine(dirs, NULL, "common_sigma")$meta_trajectory_curves; Bg <- spline_basis(CURVE_GRID); dmax <- c(ols_fit = 0, ols_se = 0, gls_fit = 0, gls_se = 0); ncmp <- 0
for (o in OUTCOMES) for (G in GROUPS5) for (tl in unname(TAILS)) { k <- SPc[outcome == o & group == G & bmi_tail == tl]; if (nrow(k) < 2) next
  st <- rbindlist(lapply(names(AD), function(s) { d <- AD[[s]][[o]]; x <- (if (G == CTL) d$ctl else d$car[locus == G])[age >= 2 & age <= 73]; x <- x[tail_of(BMI) %in% tl]
    data.table(site = paste0("SYNTHETIC_", s), age = x$age, y = x$y) }))
  st[, w := 1 / k$sigma2[match(site, k$site)]]; X <- spline_basis(st$age); ncmp <- ncmp + 1
  m0 <- lm.fit(X, st$y); V0 <- chol2inv(chol(crossprod(X))) * sum(m0$residuals^2) / (nrow(X) - ncol(X))
  bw <- solve(crossprod(X * sqrt(st$w)), crossprod(X * st$w, st$y)); Vw <- solve(crossprod(X * sqrt(st$w)))
  for (v in list(list("ols", mcs, drop(m0$coefficients), V0), list("gls", CV, drop(bw), Vw))) { pc <- v[[2]][scope == "Pooled (fixed effects)" & outcome == o & group == G & bmi_tail == tl]; ix <- match(pc$age, CURVE_GRID)
    dmax[[paste0(v[[1]], "_fit")]] <- max(dmax[[paste0(v[[1]], "_fit")]], abs(pc$fit - drop(Bg[ix, ] %*% v[[3]])))
    dmax[[paste0(v[[1]], "_se")]] <- max(dmax[[paste0(v[[1]], "_se")]], abs(pc$se - sqrt(rowSums((Bg[ix, ] %*% v[[4]]) * Bg[ix, ])))) } }
check(ncmp >= 12 && dmax[["ols_fit"]] < 1e-6 && dmax[["ols_se"]] < 1e-6, sprintf("pooled curve (one residual variance) = least squares on the stacked data: %d curves, max difference fit %.1e, SE %.1e", ncmp, dmax[["ols_fit"]], dmax[["ols_se"]]))
check(dmax[["gls_fit"]] < 1e-6 && dmax[["gls_se"]] < 1e-6, sprintf("pooled curve (each site's covariance) = weighted least squares on the stacked data: max difference fit %.1e, SE %.1e", dmax[["gls_fit"]], dmax[["gls_se"]]))
BN <- fread(file.path(MO, "meta_trajectory_bins.csv")); b1 <- BN[outcome == "Height" & group == "No CNV" & bmi_tail == TAILS[["high"]] & age_bin == "6-8"]
check(nrow(b1) == 3 && abs(b1[scope == "Pooled (fixed effects)", mean_z] - sum(b1[scope != "Pooled (fixed effects)", n * mean_z]) / sum(b1[scope != "Pooled (fixed effects)", n])) < 1e-10, "binned means are pooled by n-weighting")

# ---- (3) suppression
SB <- fread(file.path(dirs[["B"]], "sample_n.csv"), colClasses = list(character = "n"))
check(SB[outcome == "Height" & group == "22q11.2 dup" & age_band == "10-18" & sex == "all", n] == "<5", "the 3-person cell (site B, 22q11.2 dup, 10-18 y) is written as <5")
x <- SB[outcome == "Height" & group == "22q11.2 dup"]; nonnum <- function(v) is.na(suppressWarnings(as.numeric(v)))
check(all(x[age_band == "10-18", sum(nonnum(n))] >= 2, x[sex == "all", sum(nonnum(n))] >= 2) && nonnum(SB[outcome == "Height" & group == "All four CNVs" & age_band == "10-18" & sex == "all", n]),
      "the small cell cannot be recovered from the sex, age-band or all-genotype totals")
allfiles <- c(list.files(dirs, pattern = "[.]csv$", full.names = TRUE), list.files(MO, pattern = "[.]csv$", full.names = TRUE), list.files(FG, pattern = "[.]csv$", full.names = TRUE))
small <- 0; for (p in list.files(dirs, pattern = "[.]csv$", full.names = TRUE)) { d <- fread(p, colClasses = "character"); for (cn in grep("^n$|^n_", names(d), value = TRUE)) { if (cn %in% c("n_families", "n_reference_adults", "n_scored")) next
  v <- suppressWarnings(as.numeric(d[[cn]])); small <- small + sum(!is.na(v) & v > 0 & v < 5) } }
check(small == 0, sprintf("no count between 1 and 4 in any site output (found %d)", small))
E_B <- fread(file.path(dirs[["B"]], "estimates.csv"), na.strings = c("", "NA"))
check(all(is.na(E_B[analysis == "cross-sectional band" & group == "22q11.2 dup" & stratum == "10-18", estimate])) && nrow(E_B[stratum %in% c("19-39", "40-73") & status == "ok"]) == 0,
      "site B: no estimate for the small cell, and none for adults (too few non-carrier adults to score them)")
TB <- fread(file.path(dirs[["B"]], "trajectory_bins.csv"), colClasses = list(character = "n")); check(all(is.na(TB[grepl("^<", n), mean_z])), "suppressed age bins carry no mean")

# ---- (4) no ids, (6) labels
ids <- unique(unlist(lapply(FL, function(x) { d <- fread(x$data, select = c("person", "fam_id")); c(d$person, d$fam_id) })))
leak <- 0; badcol <- character(0)
for (p in allfiles) { d <- fread(p, colClasses = "character"); badcol <- c(badcol, grep("^(id|person|fam|family|fam_id|subject|participant)$", names(d), value = TRUE, ignore.case = TRUE))
  for (cn in names(d)) leak <- leak + sum(grepl("FAKE", d[[cn]], fixed = TRUE)) + sum(d[[cn]] %in% ids) }
check(leak == 0 && !length(badcol), sprintf("no id in any output (%d files; %d cells with an id; id-like columns: %s)", length(allfiles), leak, if (length(badcol)) paste(badcol, collapse = ",") else "none"))
lab_ok <- all(vapply(c("sample_n", "estimates", "tests", "autism_by_bmi_band", "trajectory_bins", "trajectory_spline", "trajectory_spline_matrices", "scoring_info"), function(nm) {
  d <- fread(file.path(dirs[["A"]], paste0(nm, ".csv")), colClasses = "character"); all(c("site", "script_version", "run_date", "scoring", "estimator") %in% names(d)) && !anyNA(d$site) && all(nzchar(d$estimator)) }, logical(1)))
check(lab_ok, "every table row carries site, script version, date, scoring and estimator")
sp <- fread(file.path(dirs[["A"]], "trajectory_spline.csv")); check(all(sp$knots_age == "6;11;16;30") && all(sp$boundary_knots_age == "2;73"), "the spline knots are written into the output")

# ---- (5) meta-analysis by hand
k <- AE[analysis == "cross-sectional band" & outcome == "Height" & group == "16p11.2 del" & stratum == "2-9" & term == "CNV x BMI"]; s2 <- k[site != "Pooled (fixed effects)"]; w <- 1 / s2$se^2
r <- get1("cross-sectional band", "Height", "16p11.2 del", "2-9", "CNV x BMI")
check(abs(r$estimate - sum(w * s2$estimate) / sum(w)) < 1e-10 && abs(r$se - sqrt(1 / sum(w))) < 1e-10 && abs(r$Q - sum(w * (s2$estimate - r$estimate)^2)) < 1e-8 && r$Q_df == 1, "meta-analysis equals the hand calculation (estimate, SE, Q)")
check(!is.na(r$se_family_clustered) && r$n_carriers == as.character(sum(as.numeric(s2$n_carriers))), "family-clustered SE pooled, carrier counts summed")
H <- fread(file.path(MO, "meta_heterogeneity_across_genotypes.csv")); h <- H[analysis == "cross-sectional band" & outcome == "Height" & stratum == "2-9" & term == "CNV x BMI"]
check(nrow(h) == 1 && h$Q_p < 0.01, sprintf("heterogeneity across genotypes detected where it was simulated (Q p = %.1e)", h$Q_p))

# ---- (7) figures
pn <- list.files(FG, pattern = "[.]png$", full.names = TRUE)
check(length(pn) >= 9 && all(file.size(pn) > 20000) && all(file.exists(sub("[.]png$", ".pdf", pn))), sprintf("figures written (%d PNG, each with a PDF)", length(pn)))
cat(sprintf("\nwork folder: %s\n%s\n", normalizePath(work), if (nfail == 0) "ALL TESTS PASSED" else sprintf("%d TEST(S) FAILED", nfail)))
quit(status = if (nfail == 0) 0 else 1)
