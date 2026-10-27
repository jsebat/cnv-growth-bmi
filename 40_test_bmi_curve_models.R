# Which model for the prime-panel BMI-by-age curves reproduces the raw within-age-band BMI slopes? (JS 2026-10-27)
# Jon: the 16p11.2 dup C' curve was positive in childhood while panel C and the raw 2-9 y slope (-0.37) are negative.
# Candidates, each fitted separately within a group (carriers of one CNV, or controls):
#   tensor_g1.5 : y ~ s(age, 10) + s(BMI, 6) + ti(age, BMI, 6 x 5), gamma 1.5   (current)
#   tensor_g1   : same, gamma 1
#   vc_k6       : y ~ s(age, 10) + s(age, 6, by = BMI)   varying-coefficient: linear in BMI, slope a smooth function of age
#   vc_k8       : y ~ s(age, 10) + s(age, 8, by = BMI)
#   vc_k6_g1.5  : vc_k6 with gamma 1.5
# Curve = half the predicted difference between BMI +1 and -1 SD; averaged over grid ages within each band (2-9, 10-18,
# 19-39, 40-73) inside the group's central 95% of ages. Raw = lm(y ~ BMI + sex + age) within the band (as in 36).
# RUNS ON EXPANSE (sources 17). Output: test_bmi_curve_models.csv
source("17_main_with_ukbb.R"); TV <- list()
BINS <- c(2, 10, 19, 40, 74); LAB <- c("2-9", "10-18", "19-39", "40-73")
MODELS <- list(
  tensor_g1.5 = list(f = y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 5)), g = 1.5),
  tensor_g1   = list(f = y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 5)), g = 1),
  vc_k6       = list(f = y ~ s(age, k = 10) + s(age, k = 6, by = BMI), g = 1),
  vc_k8       = list(f = y ~ s(age, k = 10) + s(age, k = 8, by = BMI), g = 1),
  vc_k6_g1.5  = list(f = y ~ s(age, k = 10) + s(age, k = 6, by = BMI), g = 1.5),
  vc_k8_g1.5  = list(f = y ~ s(age, k = 10) + s(age, k = 8, by = BMI), g = 1.5),
  vc_k8_g2    = list(f = y ~ s(age, k = 10) + s(age, k = 8, by = BMI), g = 2),
  vc_k8_ps    = list(f = y ~ s(age, k = 10) + s(age, k = 8, by = BMI, bs = "ps", m = c(2, 1)), g = 1))   # P-spline, first-order penalty: shrinks toward a constant slope
curve_bands <- function(d) {
  gr <- grid_of(d$age); band <- cut(gr, BINS, right = FALSE, labels = LAB)
  out <- lapply(names(MODELS), function(nm) {
    m <- bam(MODELS[[nm]]$f, data = d, method = "REML", gamma = MODELS[[nm]]$g)
    X <- (predict(m, data.frame(age = gr, BMI = 1), type = "lpmatrix") - predict(m, data.frame(age = gr, BMI = -1), type = "lpmatrix")) / 2
    sl <- drop(X %*% coef(m)); TV[[paste(nm, unique(d$locus)[1], d$y[1])]] <<- data.table(model = nm, tv = sum(abs(diff(sl))))
    data.table(model = nm, band = band, slope = sl)[, .(curve = mean(slope)), by = .(model, band)]
  })
  dcast(rbindlist(out), band ~ model, value.var = "curve")
}
raw_bands <- function(d) rbindlist(lapply(seq_along(LAB), function(b) {
  x <- d[age >= BINS[b] & age < BINS[b + 1]]; if (nrow(x) < 15) return(NULL)
  cf <- summary(lm(y ~ BMI + sex + age, data = x))$coefficients["BMI", ]
  data.table(band = LAB[b], n = nrow(x), raw = cf[1], raw_se = cf[2])
}))
one <- function(d, outcome, group) { r <- raw_bands(d); cb <- curve_bands(d); cbind(outcome = outcome, group = group, merge(r, cb, by = "band", all.x = TRUE)) }
res <- rbindlist(c(
  list(one(Hctl_test, "Height", "No CNV"), one(dh[locus == "No CNV"], "Head circumference", "No CNV")),
  lapply(LOCI, function(L) one(H[locus == L], "Height", L)),
  lapply(LOCI, function(L) one(dh[locus == L], "Head circumference", L))), fill = TRUE)
fwrite(res, file.path(WD, "test_bmi_curve_models.csv"))
cat("\n== raw within-band BMI slope vs curve average, by model\n"); print(res[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)], nrows = 100)
nm <- names(MODELS)
cat("\n== mean absolute deviation from raw (weighted by 1/se^2) and largest deviation, by model\n")
print(rbindlist(lapply(nm, function(m) res[!is.na(get(m)), .(model = m, wMAD = round(weighted.mean(abs(get(m) - raw), 1 / raw_se^2), 3), maxdev = round(max(abs(get(m) - raw)), 2))])))
cat("
== mean total variation of the curves (wiggliness), by model
"); print(rbindlist(TV)[, .(tv = round(mean(tv), 3)), by = model])

