# NOTE 2026-10-27: era_basis is now a natural cubic spline (smooth through the knots); era_select_knots still searches
# knot positions with piecewise-linear hinges, which is only used to locate the eras.
# BMI effect on an outcome as a piecewise-linear function of age, with knots at developmental breakpoints estimated in the
# CONTROLS and then held fixed for every carrier group (JS 2026-10-27, Jon: "add notches at the inflection points so that the
# key eras of development are fit accordingly"; "the same knots in B and B'"; "the knots are determined based on panel A").
#
# Knot selection (era_select_knots, controls only, per outcome): for K = 1, 2, 3 knots on the grid 3-30 y (step 1 y, knots at
# least 3 y apart), fit y ~ ns(age, 8) + BMI x {1, age, (age - k1)+, ..., (age - kK)+} by least squares and keep the knot set
# with the lowest BIC across all K. Least squares with a fixed-df natural spline for the age main effect keeps the search fast;
# the chosen knots are then used unchanged in the REML fits below and in the trajectory curves (natural cubic splines of age
# with the same interior knots), so panels A-E and A'-E' share one set of eras.
# Model (bam, REML; only the age main effects are smoothed, the BMI-by-age terms are unpenalised):
#   controls only       : y ~ s(age) + BMI x basis(age)
#   carriers vs controls: y ~ CNV + s(age, by CNV) + BMI x basis(age) + CNV x BMI x basis(age)
# slope(age) = d E(y) / d BMI = basis(age) . coefficients (exact, the model is linear in BMI), pointwise SE from Vp.
# Tests (Wald, chi-square on Vp): p_any = all carrier terms (K + 2 df: the carriers' BMI effect differs from the controls' at
# some age); p_age = carrier terms other than the constant (K + 1 df: the difference changes with age; the analogue of the
# tensor CNV x age x BMI test). Era contrasts: carrier - control slope at 3, 8, 14, 30 and 55 y.
ERA_KNOTS <- c(5, 11, 18)   # default only; the figure scripts use the knots selected in the controls
ERA_BOUNDARY <- c(2, 74)   # fixed boundary knots, so every group and every prediction grid uses the same basis
era_basis <- function(a, knots = ERA_KNOTS) { B <- cbind(1, splines::ns(a, knots = knots, Boundary.knots = ERA_BOUNDARY)); colnames(B) <- paste0("b", 0:(ncol(B) - 1)); B }   # JS 2026-10-27 (Jon: "the curve needs some rounding"): BMI effect = natural cubic spline in age with the era knots, instead of the piecewise-linear hinges
era_terms <- function(x, knots) {   # adds BMI x basis columns (b*) and, if is_cnv is present, carrier x BMI x basis columns (c*)
  B <- era_basis(x$age, knots) * x$BMI; for (j in seq_len(ncol(B))) data.table::set(x, j = colnames(B)[j], value = B[, j])
  if ("is_cnv" %in% names(x)) { cc <- as.numeric(x$is_cnv == "CNV"); for (j in seq_len(ncol(B))) data.table::set(x, j = sub("^b", "c", colnames(B)[j]), value = B[, j] * cc) }
  x
}
vp <- function(m, nm) { i <- match(nm, names(coef(m))); m$Vp[i, i, drop = FALSE] }   # Vp has no dimnames
wald_p <- function(m, terms) { b <- coef(m)[terms]; V <- vp(m, terms); s <- tryCatch(drop(t(b) %*% solve(V, b)), error = function(e) NA_real_); pchisq(s, length(terms), lower.tail = FALSE) }
era_convergence_age <- function(ctl, tails = TRUE, df_age = 8, df_gap = 6, from = 5, grid = seq(2, 73, 0.25)) {
  # JS 2026-10-27 (Jon: anchor the deviations with a knot where the high- and low-BMI curves converge): age at which the
  # controls' BMI gap (high vs low tail, or the continuous BMI slope) first changes sign after `from` years; NA if it never does.
  x <- data.table::copy(ctl)[!is.na(y) & !is.na(age) & !is.na(BMI), .(y, age, BMI)]
  if (tails) x <- x[BMI >= qnorm(0.7) | BMI <= qnorm(0.3)][, BMI := as.numeric(BMI >= qnorm(0.7))]
  bk <- range(c(x$age, grid)); A <- splines::ns(x$age, df = df_age, Boundary.knots = bk); G <- splines::ns(x$age, df = df_gap, Boundary.knots = bk)
  f <- .lm.fit(cbind(1, A, x$BMI, x$BMI * G), x$y); b <- f$coefficients; ng <- ncol(G)
  gap <- b[ncol(A) + 2] + drop(predict(G, grid) %*% b[(ncol(A) + 3):(ncol(A) + 2 + ng)])
  s <- sign(gap); i <- which(grid >= from & c(FALSE, diff(s) != 0)); if (!length(i)) return(NA_real_)
  j <- i[1]; grid[j - 1] + (grid[j] - grid[j - 1]) * abs(gap[j - 1]) / (abs(gap[j - 1]) + abs(gap[j]))   # linear interpolation of the crossing
}
era_select_knots <- function(ctl, Kmax = 3, grid = 3:30, min_gap = 3, df_age = 8, tails = FALSE, anchor = NULL) {
  # tails = TRUE: select on the high (>= 70th centile) vs low (<= 30th) BMI gap shown in panel A; anchor: a knot that is always
  # included (the convergence age, rounded to 0.5 y); the other knots are chosen by BIC among sets at least min_gap from it.
  x <- data.table::copy(ctl)[!is.na(y) & !is.na(age) & !is.na(BMI), .(y, age, BMI)]
  if (tails) x <- x[BMI >= qnorm(0.7) | BMI <= qnorm(0.3)][, BMI := as.numeric(BMI >= qnorm(0.7))]
  anchor <- if (length(anchor) && is.finite(anchor)) round(anchor * 2) / 2 else numeric(0)
  A <- splines::ns(x$age, df = df_age); base <- cbind(1, A, x$BMI, x$BMI * x$age)
  bic_of <- function(kn) { X <- cbind(base, if (length(kn)) sapply(kn, function(k) x$BMI * pmax(x$age - k, 0)) else NULL); f <- .lm.fit(X, x$y); n <- nrow(X); n * log(sum(f$residuals^2) / n) + ncol(X) * log(n) }
  cand <- grid[if (length(anchor)) abs(grid - anchor) >= min_gap else rep(TRUE, length(grid))]
  res <- list(data.frame(K = length(anchor), knots = paste(anchor, collapse = ','), bic = bic_of(anchor)))
  for (K in seq_len(Kmax - length(anchor))) {
    if (length(cand) < K) break
    cmb <- utils::combn(cand, K); if (K > 1) cmb <- cmb[, apply(cmb, 2, function(v) all(diff(v) >= min_gap)), drop = FALSE]
    sets <- lapply(seq_len(ncol(cmb)), function(j) sort(c(anchor, cmb[, j])))
    res[[K + 1]] <- data.frame(K = K + length(anchor), knots = sapply(sets, paste, collapse = ','), bic = sapply(sets, bic_of))
  }
  res <- do.call(rbind, res); best <- res[which.min(res$bic), ]
  kn <- if (best$knots == '') numeric(0) else as.numeric(strsplit(best$knots, ',')[[1]])
  attr(kn, 'bic_table') <- do.call(rbind, lapply(split(res, res$K), function(d) d[which.min(d$bic), ])); attr(kn, 'anchor') <- anchor; kn
}
era_fit <- function(xc, ctl, knots, k_age = 10) {   # xc = carriers (NULL for the control row), ctl = controls; both with y, age, BMI
  nb <- length(knots) + 2; bt <- paste0("b", 0:(nb - 1)); ct <- paste0("c", 0:(nb - 1))
  if (is.null(xc)) {
    x <- era_terms(data.table::copy(ctl)[, .(y, age, BMI)], knots)
    f <- as.formula(paste("y ~ s(age, k =", k_age, ") +", paste(bt, collapse = " + ")))
  } else {
    x <- rbind(data.table::copy(ctl)[, .(y, age, BMI, is_cnv = "No CNV")], data.table::copy(xc)[, .(y, age, BMI, is_cnv = "CNV")])[, is_cnv := factor(is_cnv, levels = c("No CNV", "CNV"))]
    x <- era_terms(x, knots)
    f <- as.formula(paste("y ~ is_cnv + s(age, k =", k_age, ", by = is_cnv) +", paste(c(bt, ct), collapse = " + ")))
  }
  m <- mgcv::bam(f, data = x, method = "REML"); m$era_knots <- knots; m
}
era_slopes <- function(m, grid, carriers = TRUE) {   # data.table(who, age, slope, se) for controls and (if fitted) carriers
  E <- era_basis(grid, m$era_knots); bn <- colnames(E); cn <- sub("^b", "c", bn); B <- coef(m)
  out <- list(data.table::data.table(who = "Non-carriers", age = grid, slope = drop(E %*% B[bn]), se = sqrt(rowSums((E %*% vp(m, bn)) * E))))
  if (carriers && all(cn %in% names(B))) {
    Ec <- cbind(E, E); idx <- c(bn, cn)
    out[[2]] <- data.table::data.table(who = "Carriers", age = grid, slope = drop(Ec %*% B[idx]), se = sqrt(rowSums((Ec %*% vp(m, idx)) * Ec)))
  }
  data.table::rbindlist(out)
}
era_tests <- function(m) {   # p_any, p_age, and carrier - control slope differences at era centres
  cn <- sub("^b", "c", colnames(era_basis(1, m$era_knots))); if (!all(cn %in% names(coef(m)))) return(data.table::data.table(p_any = NA_real_, p_age = NA_real_))
  A <- c(3, 8, 14, 30, 55); E <- era_basis(A, m$era_knots); d <- drop(E %*% coef(m)[cn]); s <- sqrt(rowSums((E %*% vp(m, cn)) * E))
  data.table::data.table(p_any = wald_p(m, cn), p_age = wald_p(m, cn[-1]),
    diff_3 = d[1], p_3 = 2 * pnorm(-abs(d[1] / s[1])), diff_8 = d[2], p_8 = 2 * pnorm(-abs(d[2] / s[2])), diff_14 = d[3], p_14 = 2 * pnorm(-abs(d[3] / s[3])),
    diff_30 = d[4], p_30 = 2 * pnorm(-abs(d[4] / s[4])), diff_55 = d[5], p_55 = 2 * pnorm(-abs(d[5] / s[5])))
}
