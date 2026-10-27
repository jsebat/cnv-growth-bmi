# One scoring scheme for every subject in every cohort (JS 2026-10-13: "everyone should be scored the same way").
# Raw measurements -> Z with external references, identical for clinical cohorts, SSC siblings and UK Biobank:
#   height : WHO 2006 (< 5 y) / WHO 2007 (5-19 y); adults scored at the age-19 reference (growth complete)
#   head circumference : UK90; adults at the top reference age (18 y)
#   BMI, children (< 19 y): WHO 2006 / WHO 2007
#   BMI, adults (>= 19 y): no external age-specific reference exists, so BMI is placed on the centile scale of ALL
#     non-carrier adults pooled across cohorts (SSC/SVIP parents + UK Biobank), by sex and 5-year age band, and
#     converted to a normal score (qnorm of the mid-rank centile). Same rule for every cohort.
# Reference LMS tables (growth_refs_lms.csv) were exported from childsds 0.9.11 (who.ref, who2007.ref, uk1990.ref;
#   WHO 2007 is BCPE with tau = 2, i.e. LMS). Z = ((y/M)^L - 1) / (L*S), parameters linearly interpolated in age.
#   Base R only, so the scoring is identical locally and on Expanse. sex: "male"/"female".
.REFDIR <- if (exists("HERE")) HERE else getwd()
.LMS <- read.csv(file.path(.REFDIR, "growth_refs_lms.csv"))

.lmsz <- function(y, age, sex, ref, item) {
  z <- rep(NA_real_, length(y))
  for (s in c("male", "female")) {
    t <- .LMS[.LMS$ref == ref & .LMS$item == item & .LMS$sex == s, ]; t <- t[!is.na(t$mu) & !is.na(t$sigma) & !is.na(t$nu), ]; t <- t[order(t$age), ]   # JS 2026-10-21: drop reference ages with missing LMS parameters (UK90 female HC is NA above 17.67 y, which made every adult female HC Z NA)
    i <- which(sex == s & !is.na(y) & !is.na(age) & y > 0)
    if (!length(i)) next
    a <- pmin(pmax(age[i], min(t$age)), max(t$age))
    L <- approx(t$age, t$nu, a)$y; M <- approx(t$age, t$mu, a)$y; S <- approx(t$age, t$sigma, a)$y
    z[i] <- ifelse(abs(L) < 1e-8, log(y[i] / M) / S, ((y[i] / M)^L - 1) / (L * S))
  }
  z[abs(z) > 10] <- NA
  z
}
score_height <- function(h, age, sex) {
  a <- pmin(age, 19); y <- !is.na(a) & a < 5
  ifelse(y, .lmsz(h, a, sex, "who", "height"), .lmsz(h, a, sex, "who2007", "height"))
}
score_hc <- function(hc, age, sex) .lmsz(hc, pmin(age, 18), sex, "uk1990", "hc")
score_bmi_child <- function(b, age, sex) {
  a <- pmin(age, 19); y <- !is.na(a) & a < 5
  ifelse(y, .lmsz(b, a, sex, "who", "bmi"), .lmsz(b, a, sex, "who2007", "bmi"))
}
adult_band <- function(age) cut(age, c(18.99, seq(24, 74, 5), 120), labels = FALSE)
# b/age/sex: everyone to score; ref_*: pooled non-carrier adults (age >= 19)
score_bmi <- function(b, age, sex, ref_b, ref_age, ref_sex) {
  z <- score_bmi_child(b, age, sex); ad <- !is.na(age) & age >= 19
  z[ad] <- NA
  rb <- adult_band(ref_age); xb <- adult_band(age)
  for (s in c("male", "female")) for (k in unique(xb[ad & sex == s])) {
    ref <- sort(ref_b[ref_sex == s & rb == k & !is.na(ref_b)]); i <- which(ad & sex == s & xb == k & !is.na(b))
    if (length(ref) < 50 || !length(i)) next
    pct <- (findInterval(b[i], ref, left.open = TRUE) + findInterval(b[i], ref)) / 2 / length(ref)   # mid-rank centile
    z[i] <- qnorm(pmin(pmax(pct, 0.001), 0.999))
  }
  z
}

# ---------------------------------------------------------------- JS 2026-10-27: adult scaling relative to the subject's own cohort
# Jon's decision (2026-10-27, for the BMI-modifier figure with UK Biobank adults): children stay on the external growth
# references above; every ADULT is scored relative to the NON-CARRIER ADULTS OF THEIR OWN COHORT with Molly's Box-Cox
# scaling (biobank Methods, "Scaling of the response variable"): within cohort and sex, the Box-Cox lambda is chosen by
# profile likelihood on the reference fit lm(y ~ age + age^2) (MASS::boxcox), y is transformed, and Z = residual from the
# reference model on the transformed scale / the reference residual SD. Carriers are scored with the reference model, so
# a cohort's non-carrier adults have mean 0 and SD 1 by construction and carriers are read relative to them.
# "Cohort" is the reference population: all clinical cohorts form one cohort (its non-carrier adults are the SSC parents
# and SVIP non-carrier family members, as in the Fig. 6 control group; G2MH, Cardiff and Searchlight have no measured
# non-carrier adults), and UK Biobank is its own. A cohort x sex cell needs at least min_ref reference adults, else NA.
# Lambda (Jon, 2026-10-27, second decision): ONE transformation per trait and sex, estimated on the non-carrier adults of
# the anchor cohort (lambda_from = "UK Biobank", the largest reference) and applied in every cohort; the residualisation
# (lm(y ~ age + age^2) on the transformed scale) stays within cohort. Reason: estimated per cohort, the clinical height
# lambda ran to 3.6-4.2 (heavy lower tail of the SSC-parent heights) while UKBB gave 0.6, and the CNV x BMI estimates were
# the same to two decimals either way. If the anchor cohort is absent (head circumference, which UK Biobank lacks) the
# caller fixes lambda (lambda_fixed = 1, no transform: adult HC is close to normal); otherwise a cohort without the anchor
# falls back to its own estimate on the default (-2, 2) grid. info$lambda_source says which applied.
# Returns Z with attr "info" (cohort, sex, n_ref, n_scored, lambda, lambda_source, sd_resid). Requires MASS.
score_adult_bc <- function(y, age, sex, cohort, is_ref, min_age = 19, min_ref = 100, lambda_grid = seq(-2, 2, 0.05),
                           lambda_from = "UK Biobank", lambda_fixed = NULL) {
  z <- rep(NA_real_, length(y)); info <- list()
  ad <- !is.na(age) & age >= min_age & !is.na(y) & y > 0 & !is.na(sex) & !is.na(cohort)
  est_lambda <- function(r) { d <- data.frame(y = y[r], age = age[r], age2 = age[r]^2)
    b <- MASS::boxcox(y ~ age + age2, data = d, lambda = lambda_grid, plotit = FALSE); b$x[which.max(b$y)] }
  anchor <- list()
  for (s in c("male", "female")) {
    if (!is.null(lambda_fixed)) { anchor[[s]] <- list(l = lambda_fixed, src = "fixed"); next }
    r0 <- which(ad & cohort == lambda_from & sex == s & is_ref %in% TRUE)
    anchor[[s]] <- if (length(r0) >= min_ref) list(l = est_lambda(r0), src = lambda_from) else list(l = NA_real_, src = "own")
  }
  for (co in unique(cohort[ad])) for (s in c("male", "female")) {
    r <- which(ad & cohort == co & sex == s & is_ref %in% TRUE); i <- which(ad & cohort == co & sex == s)
    if (length(r) < min_ref || !length(i)) next
    l <- anchor[[s]]$l; src <- anchor[[s]]$src; if (is.na(l)) { l <- est_lambda(r); src <- "own" }
    tr <- function(v) if (abs(l) < 1e-8) log(v) else (v^l - 1) / l
    d <- data.frame(y = y[r], age = age[r], age2 = age[r]^2); d$yb <- tr(d$y); m <- lm(yb ~ age + age2, data = d); sdr <- sd(resid(m))
    z[i] <- (tr(y[i]) - predict(m, newdata = data.frame(age = age[i], age2 = age[i]^2))) / sdr
    info[[paste(co, s)]] <- data.frame(cohort = co, sex = s, n_ref = length(r), n_scored = length(i), lambda = l, lambda_source = src, sd_resid = sdr)
  }
  z[!is.na(z) & abs(z) > 10] <- NA
  attr(z, "info") <- do.call(rbind, info); z
}
# Combined scorers: external reference for children, own-cohort Box-Cox for adults (height and BMI from 19 y, HC from 18 y).
score_height_cohort <- function(h, age, sex, cohort, is_ref) {
  z <- score_height(h, age, sex); zb <- score_adult_bc(h, age, sex, cohort, is_ref, min_age = 19)
  ad <- !is.na(age) & age >= 19; z[ad] <- zb[ad]; attr(z, "info") <- attr(zb, "info"); z
}
score_bmi_cohort <- function(b, age, sex, cohort, is_ref) {
  z <- score_bmi_child(b, age, sex); zb <- score_adult_bc(b, age, sex, cohort, is_ref, min_age = 19)
  ad <- !is.na(age) & age >= 19; z[ad] <- zb[ad]; attr(z, "info") <- attr(zb, "info"); z
}
score_hc_cohort <- function(hc, age, sex, cohort, is_ref) {   # no UK Biobank head circumference: lambda fixed at 1 (see above)
  z <- score_hc(hc, age, sex); zb <- score_adult_bc(hc, age, sex, cohort, is_ref, min_age = 18, lambda_fixed = 1)
  ad <- !is.na(age) & age >= 18; z[ad] <- zb[ad]; attr(z, "info") <- attr(zb, "info"); z
}
