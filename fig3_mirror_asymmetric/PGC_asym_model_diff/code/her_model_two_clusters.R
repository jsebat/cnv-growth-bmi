#!/usr/bin/env Rscript
##############################################
## Variations on Molly's model WITHOUT tauU: can any of them give a good-fitting two-cluster pattern?
## Base = her continuous mirror m_g ~ N(mu_m, s_m^2) along a = (1, -rho) + continuous deleterious
## s_g ~ N(mu_d, s_d^2) along u_d.  Added: a CNV-level single hit (P = 1 - (1 - pi)^G) shifting by
## beta * u_h, beta ~ N(beta, tau^2).  No tauU anywhere; off-line scatter can only come from s_d (and tau).
## Toggles: rho fixed log2(3/2) / fitted; u_d = (1, rho) (hers) / free angle; tau = 0 / free.
## pi is profiled (everything else refit) so the fit cost of forcing a minority-hit solution is visible.
## Separation index = |perpendicular shift| / latent perpendicular SD of the un-hit cloud at G = 23
## (mirror spread is along the line, so only s_d contributes); > ~2 = visibly separate clusters.
## Run: Rscript code/her_model_two_clusters.R   (from PGC_asym_model_diff/)
##############################################
src <- readLines("code/factor_test.R")
eval(parse(text = src[seq_len(grep("^## All configurations", src) - 2)]))

unpack <- function(x, cfg) {
  i <- 0; nx <- function() { i <<- i + 1; x[i] }
  p <- list(mu_m = nx(), sm2 = exp(2 * nx()))
  p$rho <- if (cfg$rho == "free") exp(nx()) else RHO0
  p$mu_d <- nx(); p$sd2 <- exp(2 * nx())
  p$phi <- if (cfg$ud == "free") nx() else NA
  p$theta <- nx(); p$beta <- nx()
  p$tau2 <- if (cfg$tau == "free") exp(2 * nx()) else 0
  p
}
npar <- function(cfg) 6 + (cfg$rho == "free") + (cfg$ud == "free") + (cfg$tau == "free")   # + pi (profiled)

ll_model <- function(p, lpi, idx = rep(TRUE, length(G))) {
  if (p$sm2 > 1 || p$sd2 > 1 || p$rho > 5 || abs(p$beta) > 3 || p$tau2 > 4) return(-1e10)
  g <- G[idx]; a <- c(1, -p$rho)
  ud <- if (is.na(p$phi)) c(1, p$rho) else c(cos(p$phi), sin(p$phi)); uh <- c(cos(p$theta), sin(p$theta))
  mD <- g * (p$mu_m * a[1] + p$mu_d * ud[1]); mU <- g * (p$mu_m * a[2] + p$mu_d * ud[2])
  v11 <- g * (p$sm2 * a[1]^2 + p$sd2 * ud[1]^2) + vD[idx]
  v22 <- g * (p$sm2 * a[2]^2 + p$sd2 * ud[2]^2) + vU[idx]
  v12 <- g * (p$sm2 * a[1] * a[2] + p$sd2 * ud[1] * ud[2])
  pG <- 1 - (1 - 10^lpi)^g
  l0 <- log1p(-pG) + ld2(mD, mU, v11, v22, v12, idx)
  l1 <- log(pG) + ld2(mD + p$beta * uh[1], mU + p$beta * uh[2],
                      v11 + p$tau2 * uh[1]^2, v22 + p$tau2 * uh[2]^2, v12 + p$tau2 * uh[1] * uh[2], idx)
  sum(lse(l0, l1))
}

fit_at_pi <- function(cfg, lpi, start = NULL, nstart = 20) {
  f <- function(x) { v <- -ll_model(unpack(x, cfg), lpi); if (is.finite(v)) v else 1e10 }
  x0 <- c(-0.003, log(0.05), if (cfg$rho == "free") log(0.35), -0.005, log(0.02),
          if (cfg$ud == "free") atan2(RHO0, 1), -1.6, 0.28, if (cfg$tau == "free") log(0.05))
  nth <- 6 + (cfg$rho == "free") + (cfg$ud == "free") - 1   # index of theta
  best <- NULL
  for (s in seq_len(nstart)) {
    xs <- if (s == 1 && !is.null(start)) start else x0 + rnorm(length(x0), 0, 0.3)
    if (s > 1) { xs[nth] <- runif(1, -pi, pi); xs[nth + 1] <- runif(1, 0.05, 0.6) }
    o <- list(par = xs, value = f(xs))
    for (r in 1:3) o <- optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-11, maxit = 6000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  best
}

summ <- function(p, lpi) {
  ud <- if (is.na(p$phi)) c(1, p$rho) else c(cos(p$phi), sin(p$phi))
  n <- c(p$rho, 1) / sqrt(1 + p$rho^2)
  vec <- p$beta * c(cos(p$theta), sin(p$theta)); perp <- sum(vec * n)
  sd_perp <- sqrt(Gmed * p$sd2) * abs(sum(ud * n))
  tibble(P_hit = 1 - (1 - 10^lpi)^Gmed, shift_DEL = vec[1], shift_DUP = vec[2], perp = perp,
         sd_perp_unhit = sd_perp, separation = abs(perp) / sd_perp, rho_hat = p$rho,
         sigma_d = sqrt(p$sd2), tau_hat = sqrt(p$tau2),
         ud_angle = atan2(ud[2], ud[1]) * 180 / pi)
}

cfgs <- expand.grid(rho = c("fixed", "free"), ud = c("fixed", "free"), tau = c("zero", "free"), stringsAsFactors = FALSE)
lpis <- seq(-3, -0.4, by = 0.2)
rows <- list()
for (i in seq_len(nrow(cfgs))) {
  cfg <- as.list(cfgs[i, ]); start <- NULL
  for (lp in lpis) {
    b <- fit_at_pi(cfg, lp, start); start <- b$par
    k <- npar(cfg) + 1
    rows[[length(rows) + 1]] <- bind_cols(tibble(!!!cfg, log10_pi = lp, k = k, logL = -b$value, AIC = 2 * k + 2 * b$value),
                                          summ(unpack(b$par, cfg), lp))
  }
  cat(sprintf("[%d/%d] rho %s, u_d %s, tau %s\n", i, nrow(cfgs), cfg$rho, cfg$ud, cfg$tau))
}
prof <- bind_rows(rows)
write_csv(prof, "outputs/her_model_two_clusters_profile.csv")
saveRDS(list(cfgs = cfgs), "outputs/her_model_two_clusters_cfgs.rds")

options(width = 220)
cat("\nReference: Molly continuous (k = 4) logL 15.76, AIC -23.51\n\nBest over pi, per variant:\n")
print(prof %>% group_by(rho, ud, tau) %>% slice_max(logL, n = 1) %>% ungroup() %>%
        transmute(rho, ud, tau, k, logL = round(logL, 2), AIC = round(AIC, 2), P_hit = round(P_hit, 2),
                  dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3), perp = round(perp, 3),
                  sep = round(separation, 2), rho_hat = round(rho_hat, 3), sigma_d = round(sigma_d, 4),
                  ud_ang = round(ud_angle), tau_hat = round(tau_hat, 3)))
cat("\nBest minority-hit fit (P_hit 0.2-0.5), per variant:\n")
print(prof %>% filter(P_hit >= 0.2, P_hit <= 0.5) %>% group_by(rho, ud, tau) %>% slice_max(logL, n = 1) %>% ungroup() %>%
        transmute(rho, ud, tau, k, logL = round(logL, 2), AIC = round(AIC, 2), P_hit = round(P_hit, 2),
                  dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3), perp = round(perp, 3),
                  sep = round(separation, 2), rho_hat = round(rho_hat, 3), sigma_d = round(sigma_d, 4),
                  ud_ang = round(ud_angle), tau_hat = round(tau_hat, 3)))
