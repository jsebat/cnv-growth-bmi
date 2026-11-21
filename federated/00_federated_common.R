# 00_federated_common.R  (2026-11-20)
# What: constants and small helpers shared by site_export.R, meta_combine.R and plot_federated.R, so that every site
#   uses the same genotype labels, age bands, age bins, BMI tails and the same FIXED spline basis for the federated
#   trajectories. Do not edit at a site: the pooled curve is only valid when every site uses the same basis.
# Inputs: none. Outputs: none (constants and functions only).
# null model: none; helpers only. fe_meta() = fixed-effects inverse-variance mean with Cochran Q and I2.
FED_VERSION <- "cnv-growth-bmi federated 1.0.0 (2026-11-20)"
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"); CTL <- "No CNV"; POOL <- "All four CNVs"
GROUPS5 <- c(CTL, LOCI)
OUTCOMES <- c("Height", "Head circumference")
# age bands of the cross-sectional models (17_main_with_ukbb.R): [2,10) [10,19) [19,40) [40,74)
BANDS <- c(2, 10, 19, 40, 74); BLAB <- c("2-9", "10-18", "19-39", "40-73")
# longitudinal age windows (80_longitudinal_v2.R); window = age at the first visit of a pair
WIN_H <- c("2-8 y", "8-13 y", "13-18 y"); WIN_H_BREAKS <- c(2, 8, 13, 18)
WIN_C <- c("<2 y", "2-5 y", "5-18 y");    WIN_C_BREAKS <- c(0, 2, 5, 18)
# BMI tails of the trajectory panels: fixed cut-offs on the BMI Z scale (70th and 30th centile of a standard normal)
TAILS <- c(high = "BMI >= 70th centile", low = "BMI <= 30th centile")
TAIL_COLS <- c("BMI >= 70th centile" = "#F0645A", "BMI <= 30th centile" = "#00A5AD")
tail_of <- function(z) factor(ifelse(z >= qnorm(0.7), TAILS[["high"]], ifelse(z <= qnorm(0.3), TAILS[["low"]], NA)), levels = unname(TAILS))
BAND3 <- c("BMI <= 30th centile", "BMI 30th-70th centile", "BMI >= 70th centile")
band3 <- function(z) factor(ifelse(z >= qnorm(0.7), BAND3[3], ifelse(z <= qnorm(0.3), BAND3[1], BAND3[2])), BAND3)
COL_CAR <- "#4B6EA9"; COL_CTL <- "#8C9A9A"; COL_F <- "#C2528B"; COL_M <- "#3C7FB1"

# ---- federated trajectories: common age bins and common spline basis
AGE_BIN_BREAKS <- c(2, 4, 6, 8, 10, 13, 16, 19, 30, 40, 50, 60, 73.0001)
AGE_BIN_LABELS <- c("2-4", "4-6", "6-8", "8-10", "10-13", "13-16", "16-19", "19-30", "30-40", "40-50", "50-60", "60-73")
AGE_BIN_MID <- c(3, 5, 7, 9, 11.5, 14.5, 17.5, 24.5, 35, 45, 55, 66.5)
age_bin <- function(age) factor(as.character(cut(age, AGE_BIN_BREAKS, labels = AGE_BIN_LABELS, right = FALSE)), AGE_BIN_LABELS)
# natural cubic spline of log(age): interior knots and boundary knots in years, FIXED for every site
SPLINE_KNOTS_AGE <- c(6, 11, 16, 30); SPLINE_BOUNDARY_AGE <- c(2, 73)
SPLINE_NAMES <- c("intercept", paste0("ns", seq_len(length(SPLINE_KNOTS_AGE) + 1)))
SPLINE_ID <- sprintf("ns(log(age)); knots at %s y; boundary knots %s y; intercept + %d basis functions",
                     paste(SPLINE_KNOTS_AGE, collapse = ", "), paste(SPLINE_BOUNDARY_AGE, collapse = " and "), length(SPLINE_KNOTS_AGE) + 1)
spline_basis <- function(age) {
  B <- cbind(1, splines::ns(log(age), knots = log(SPLINE_KNOTS_AGE), Boundary.knots = log(SPLINE_BOUNDARY_AGE)))
  colnames(B) <- SPLINE_NAMES; B
}
CURVE_GRID <- c(seq(2, 19, by = 0.25), seq(19.5, 73, by = 0.5))

# ---- small statistics helpers
# fixed-effects inverse-variance meta-analysis; Q, I2 as metafor::rma(method = "FE")
fe_meta <- function(b, se) {
  ok <- is.finite(b) & is.finite(se) & se > 0; b <- b[ok]; se <- se[ok]; k <- length(b)
  if (!k) return(list(k = 0L, est = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, Q = NA_real_, Q_df = NA_real_, Q_p = NA_real_, I2 = NA_real_))
  w <- 1 / se^2; est <- sum(w * b) / sum(w); s <- sqrt(1 / sum(w)); Q <- sum(w * (b - est)^2); df <- k - 1
  list(k = k, est = est, se = s, z = est / s, p = 2 * pnorm(-abs(est / s)), Q = if (k > 1) Q else NA_real_, Q_df = if (k > 1) df else NA_real_,
       Q_p = if (k > 1) pchisq(Q, df, lower.tail = FALSE) else NA_real_, I2 = if (k > 1) max(0, 100 * (Q - df) / Q) else NA_real_)
}
# pseudo-inverse of a symmetric positive semi-definite matrix (a site whose ages do not cover the whole basis has a
# rank-deficient information matrix; directions without information get zero weight)
psd_inv <- function(A, tol = 1e-10) {
  A <- (A + t(A)) / 2; e <- eigen(A, symmetric = TRUE); keep <- e$values > tol * max(e$values, 0)
  if (!any(keep)) return(list(inv = A * 0, rank = 0L))
  U <- e$vectors[, keep, drop = FALSE]; list(inv = U %*% (t(U) / e$values[keep]), rank = sum(keep))
}
pfmt <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))
plab <- function(p) gsub("= < ", "< ", paste("p =", pfmt(p)), fixed = TRUE)
script_dir <- function(default = getwd()) {
  if (exists("FEDERATED_DIR")) return(FEDERATED_DIR)   # set before source()-ing a package script from another script
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(gsub("~+~", " ", f[1], fixed = TRUE))) else default
}
