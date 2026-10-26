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
