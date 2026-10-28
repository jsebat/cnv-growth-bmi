# High- vs low-BMI offset (BMI z >= 70th vs <= 30th centile, as plotted in Fig. 6 A-J) by age group, carriers vs controls
# (JS 2026-10-27, Jon: would this make better prime panels?). Within children (< 18 y) and adults (>= 18 y) separately:
#   offset within a group   : y ~ hi + s(age)                        (hi = 1 for the high-BMI tail, 0 for the low tail)
#   carrier - control       : y ~ cnv * hi + s(age, by = cnv) + cnv  -> coefficient cnv:hi (difference in offsets), P
# Middle BMI (30th-70th centile) excluded, as in the panels. Controls: height = clinical + 50,000 UKBB; HC = clinical.
# RUNS ON EXPANSE (sources 17). Output: tail_offset_by_agegroup.csv
source("17_main_with_ukbb.R")
AG <- list("< 18 y" = c(0, 18), ">= 18 y" = c(18, 200))
prep <- function(d) d[, hi := as.numeric(tail_of(BMI) == names(COLS)[1])][!is.na(hi)]
one_off <- function(d) { if (nrow(d) < 10 || length(unique(d$hi)) < 2) return(c(NA, NA, NA))
  k <- min(10, length(unique(round(d$age))) - 1); s <- summary(gam(y ~ hi + s(age, k = max(3, k)), data = d, method = "REML"))
  c(s$p.coeff["hi"], s$se["hi"], s$p.pv["hi"]) }
run <- function(xc, ctl, outcome, L) rbindlist(lapply(names(AG), function(a) {
  r <- AG[[a]]; c1 <- prep(copy(xc)[age >= r[1] & age < r[2]]); c0 <- prep(copy(ctl)[age >= r[1] & age < r[2]])
  oc <- one_off(c1); ok <- one_off(c0)
  dd <- NA_real_; dse <- NA_real_; dp <- NA_real_
  if (nrow(c1) >= 10 && length(unique(c1$hi)) == 2) {
    x <- rbind(c1[, .(y, age, hi, cnv = 1)], c0[, .(y, age, hi, cnv = 0)])[, cnvf := factor(cnv)]
    s <- summary(gam(y ~ cnvf * hi + s(age, k = 10, by = cnvf), data = x, method = "REML"))
    dd <- s$p.coeff["cnvf1:hi"]; dse <- s$se["cnvf1:hi"]; dp <- s$p.pv["cnvf1:hi"] }
  data.table(outcome = outcome, group = L, age_group = a, n_carriers_hi = sum(c1$hi == 1), n_carriers_lo = sum(c1$hi == 0),
             n_controls_hi = sum(c0$hi == 1), n_controls_lo = sum(c0$hi == 0),
             offset_carriers = oc[1], se_carriers = oc[2], p_carriers = oc[3], offset_controls = ok[1], se_controls = ok[2], p_controls = ok[3],
             delta = dd, se_delta = dse, p_delta = dp)
}))
res <- rbindlist(c(lapply(LOCI, function(L) run(H[locus == L], Hctl_test, "Height", L)),
                   lapply(LOCI, function(L) run(dh[locus == L], dh[locus == "No CNV"], "Head circumference", L))))
fwrite(res, file.path(WD, "tail_offset_by_agegroup.csv"))
cat("\n== high - low BMI tail offset by age group (SD), carriers vs controls\n")
print(res[, .(outcome = substr(outcome, 1, 6), group, age_group, n_car = n_carriers_hi + n_carriers_lo, off_car = round(offset_carriers, 2),
              off_ctl = round(offset_controls, 2), delta = round(delta, 2), se = round(se_delta, 2), p = signif(p_delta, 2))], nrows = 50)
