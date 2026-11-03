#!/usr/bin/env Rscript
##############################################
## Can a two-cluster model fit without tauU?  Molly's continuous mirror + continuous deleterious model
## (deleterious along (1, rho), sigma_d supplies the off-line scatter) PLUS a fixed-size single hit in a
## free direction.  No tauU.  Joint fit on all 34 height loci; rho fixed at log2(3/2) or fitted.
## Run: Rscript code/hybrid_no_tauU.R   (from PGC_asym_model_diff/)
##############################################
src <- readLines("code/factor_test.R")
eval(parse(text = src[seq_len(grep("^## All configurations", src) - 2)]))
ll_hyb <- function(x, rho_free) {
  i <- 0; nx <- function() { i <<- i + 1; x[i] }
  mu_m <- nx(); sm2 <- exp(2 * nx()); rho <- if (rho_free) exp(nx()) else RHO0
  mu_d <- nx(); sd2 <- exp(2 * nx()); th <- nx(); pi_ <- 10^nx(); beta <- nx()
  if (sm2 > 1 || sd2 > 1 || pi_ >= 1 || pi_ < 1e-6 || abs(beta) > 3 || rho > 5) return(-1e10)
  a <- c(1, -rho); b <- c(1, rho); u <- c(cos(th), sin(th))
  mD <- G * (mu_m * a[1] + mu_d * b[1]); mU <- G * (mu_m * a[2] + mu_d * b[2])
  v11 <- G * (sm2 + sd2) + vD; v22 <- G * rho^2 * (sm2 + sd2) + vU; v12 <- G * rho * (sd2 - sm2)
  pG <- 1 - (1 - pi_)^G
  sum(lse(log1p(-pG) + ld2(mD, mU, v11, v22, v12), log(pG) + ld2(mD + beta * u[1], mU + beta * u[2], v11, v22, v12)))
}
for (rf in c(FALSE, TRUE)) {
  f <- function(x) { v <- -ll_hyb(x, rf); if (is.finite(v)) v else 1e10 }
  x0 <- c(-0.003, log(0.05), if (rf) log(0.35), -0.005, log(0.02), 0, log10(0.02), 0.3)
  best <- NULL
  for (s in 1:60) {
    xs <- x0 + rnorm(length(x0), 0, 0.3); xs[length(x0) - 2] <- runif(1, -pi, pi); xs[length(x0) - 1] <- runif(1, -3.5, -0.5)
    o <- list(par = xs, value = f(xs))
    for (r in 1:3) o <- optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-11, maxit = 6000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  x <- best$par; k <- length(x); n <- length(x0)
  rho <- if (rf) exp(x[3]) else RHO0; th <- x[n - 2]; pG <- 1 - (1 - 10^x[n - 1])^Gmed; beta <- x[n]
  vec <- beta * c(cos(th), sin(th)); perp <- sum(vec * c(rho, 1)) / sqrt(1 + rho^2)
  cat(sprintf("rho %s: k=%d logL=%.2f AIC=%.2f (Molly continuous: logL 15.76, AIC -23.51)\n  rho=%.3f sigma_d=%.4f P(hit|G=23)=%.2f shift=(%.3f, %.3f) perp=%.3f\n",
              if (rf) "fitted" else "fixed 0.585", k, -best$value, 2 * k + 2 * best$value, rho, exp(x[if (rf) 6 else 5]), pG, vec[1], vec[2], perp))
}
