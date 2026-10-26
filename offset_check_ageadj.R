# JS 2026-10-21 diagnostic: why is "offset vs controls" less significant when stratified by age group?
# Reuses the data-building part of 17_main_with_ukbb.R (everything up to Hctl_plot), then compares
#   m0 (v1): y ~ is_cnv * hi + s(age, by = is_cnv)                       BMI-tail effect assumed constant over age
#   m1     : y ~ is_cnv * hi + s(age, by = is_cnv) + s(age, by = hi)     BMI-tail effect allowed to vary with age (shared)
# is_cnvCNV:hi = carrier - non-carrier offset. If m1 kills the signal, the v1 p came from age confounding of the BMI effect.
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(1)
WD <- getwd()
ex <- parse(file.path(WD, "17_main_with_ukbb.R"))
k <- which(sapply(ex, function(e) grepl("Hctl_plot <-", paste(deparse(e), collapse = ""))))[1]
for (i in seq_len(k)) eval(ex[[i]], envir = globalenv())
dh <- latest("HEADCIRCUM_raw")[, .(id, source = "Clinical", locus = as.character(locus), sex = as.character(sex), age,
                                   y = score_hc(HEADCIRCUM_raw, age, as.character(sex)), BMI = sb_(BMI_raw, age, as.character(sex)))][!is.na(y) & !is.na(BMI)]
prep <- function(xc, ctl) {
  x <- rbind(xc, ctl)[, grp := tail_of(BMI)][!is.na(grp)]
  x[, `:=`(hi = as.numeric(grp == names(COLS)[1]), is_cnv = factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV")))]
  x[, hif := factor(hi)]
  x
}
one <- function(x) {
  m0 <- bam(y ~ is_cnv * hi + s(age, k = 10, by = is_cnv), data = x, method = "REML")
  m1 <- bam(y ~ is_cnv * hi + s(age, k = 10, by = is_cnv) + s(age, k = 10, by = hi), data = x, method = "REML")
  c0 <- summary(m0)$p.table; c1 <- summary(m1)$p.table
  # age composition of the carriers' BMI tails and of the controls' tails
  ag <- function(d) sprintf("%.0f%%<13y, %.0f%%>=40y", 100 * mean(d$age < 13), 100 * mean(d$age >= 40))
  data.table(n_cnv_hi = x[is_cnv == "CNV" & hi == 1, .N], n_cnv_lo = x[is_cnv == "CNV" & hi == 0, .N],
             cnv_ages = ag(x[is_cnv == "CNV"]), ctl_ages = ag(x[is_cnv == "No CNV"]),
             v1_ctl_offset = c0["hi", 1], v1_diff = c0["is_cnvCNV:hi", 1], v1_p = c0["is_cnvCNV:hi", 4],
             ageadj_diff = c1["is_cnvCNV:hi", 1], ageadj_se = c1["is_cnvCNV:hi", 2], ageadj_p = c1["is_cnvCNV:hi", 4])
}
res <- rbindlist(lapply(LOCI, function(L) rbind(
  cbind(outcome = "Height", group = L, one(prep(H[locus == L], Hctl_test))),
  cbind(outcome = "Head circumference", group = L, one(prep(dh[locus == L], dh[locus == "No CNV"]))))))
fwrite(res, file.path(WD, "offset_check_ageadj.csv"))
print(res, digits = 3)

# stratified carrier - non-carrier offset difference under coarser age groupings
strat <- function(x, brk, labs) {
  x <- copy(x)[, ag := cut(age, brk, labels = labs, right = FALSE, include.lowest = TRUE)]
  rbindlist(lapply(labs, function(g) {
    d <- x[ag == g]; n <- c(d[is_cnv == "CNV" & hi == 1, .N], d[is_cnv == "CNV" & hi == 0, .N])
    if (min(n) < 5) return(data.table(age_group = g, n_hi = n[1], n_lo = n[2], diff = NA_real_, p = NA_real_))
    m <- tryCatch(gam(y ~ is_cnv * hi + s(age, k = 5, by = is_cnv), data = d, method = "REML"), error = function(e) lm(y ~ is_cnv * hi + age * is_cnv, data = d))
    co <- if (inherits(m, "gam")) summary(m)$p.table else summary(m)$coefficients
    data.table(age_group = g, n_hi = n[1], n_lo = n[2], diff = co["is_cnvCNV:hi", 1], p = co["is_cnvCNV:hi", 4])
  }))
}
G2 <- list(brk = c(2, 18, 74), labs = c("2-17 y", "18-73 y"))
G3 <- list(brk = c(2, 13, 18, 74), labs = c("2-12 y", "13-17 y", "18-73 y"))
st <- rbindlist(lapply(LOCI, function(L) rbind(
  cbind(outcome = "Height", group = L, grouping = "2 groups", strat(prep(H[locus == L], Hctl_test), G2$brk, G2$labs)),
  cbind(outcome = "Height", group = L, grouping = "3 groups", strat(prep(H[locus == L], Hctl_test), G3$brk, G3$labs)),
  cbind(outcome = "Head circumference", group = L, grouping = "2 groups", strat(prep(dh[locus == L], dh[locus == "No CNV"]), G2$brk, G2$labs)),
  cbind(outcome = "Head circumference", group = L, grouping = "3 groups", strat(prep(dh[locus == L], dh[locus == "No CNV"]), G3$brk, G3$labs)))))
fwrite(st, file.path(WD, "offset_check_groupings.csv"))
print(st, digits = 3, nrows = 200)
