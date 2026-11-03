#!/usr/bin/env Rscript
##############################################
## Why does Molly's v4 asymmetric (CNV-hit) model not reproduce the perpendicular,
## clustered single-hit shift of Jon's v2 model?  Height only.
##
## One likelihood (Molly's cnv_hit_deleterious_model.R generalised) with toggles for each
## structural difference between the two models:
##   dir   : "fixed" = shift along (1, rho) (Molly);  "free" = any angle theta (Jon, but unbounded)
##   tau   : "free"  = hit size ~ N(beta, tau^2) (Molly);  "zero" = fixed hit vector (Jon)
##   rho   : "fixed" = log2(3/2) (Molly);  "free" = fitted (Jon: 0.351)
##   tauU  : "zero"  = no off-line residual (Molly);  "free" = extra DUP-only variance (Jon: 0.051)
##   fit   : "joint" = mirror + shift fit together on all loci (Molly)
##           "cond"  = mirror fit on the 21 mirror-class loci, fixed, shift fit on all loci (Jon)
## Each hit model is compared with a continuous deleterious model with the same toggles.
##
## y = (DEL, DUP) ~ N2(G mu_m a + S u, G s_m^2 a a' + v_S u u' + diag(SE_D^2, SE_U^2 + tauU^2))
## a = (1, -rho); u = (1, rho) if dir fixed, (cos theta, sin theta) if free.
## Run: Rscript code/factor_test.R   (from PGC_asym_model_diff/)
##############################################
suppressMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork); library(tidyr)})
set.seed(1)
dat <- read_csv("inputs/height_loci.csv", show_col_types = FALSE)
G <- dat$Gene_count; D <- dat$Deletion_beta; U <- dat$Duplication_beta
vD <- dat$Deletion_SE^2; vU <- dat$Duplication_SE^2
is_mirror <- dat$Dose_response_class == "mirror"
RHO0 <- log2(3 / 2); Gmed <- median(G)

# per-locus bivariate normal log density; all arguments vectors over loci (or scalars)
ld2 <- function(mD, mU, v11, v22, v12, idx = TRUE) {
  det <- v11 * v22 - v12^2
  dx <- D[idx] - mD; dy <- U[idx] - mU
  -log(2 * base::pi) - 0.5 * log(det) - 0.5 * (v22 * dx^2 + v11 * dy^2 - 2 * v12 * dx * dy) / det
}
lse <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }

# unpack parameter vector x into a named list, given the config
pnames <- function(cfg, kind) {
  mir <- if (cfg$fit == "cond") character(0) else
    c("mu_m", "lsm", if (cfg$rho == "free") "lrho", if (cfg$tauU == "free") "ltU")
  dir <- if (cfg$dir == "free") "theta"
  sh <- switch(kind, mirror = character(0),
               hit = c("lpi", "beta", if (cfg$tau == "free") "ltau"),
               cont = c("mu_d", "lsd"))
  c(mir, if (kind != "mirror") dir, sh)
}
build <- function(x, nm, fixed) { p <- fixed; p[nm] <- x; p }

loglik <- function(p, kind, cfg, idx = rep(TRUE, length(G))) {
  g <- G[idx]; rho <- exp(p$lrho); tU2 <- exp(2 * p$ltU); sm2 <- exp(2 * p$lsm)
  if (!is.finite(rho) || rho > 5 || sm2 > 1 || tU2 > 1) return(-1e10)
  u <- if (cfg$dir == "free") c(cos(p$theta), sin(p$theta)) else c(1, rho)
  mD0 <- g * p$mu_m; mU0 <- -rho * g * p$mu_m
  v11 <- g * sm2 + vD[idx]; v22 <- rho^2 * g * sm2 + vU[idx] + tU2; v12 <- -rho * g * sm2
  if (kind == "mirror") return(sum(ld2(mD0, mU0, v11, v22, v12, idx)))
  if (kind == "cont") {
    sd2 <- exp(2 * p$lsd); if (sd2 > 1) return(-1e10)
    mS <- g * p$mu_d; vS <- g * sd2
    return(sum(ld2(mD0 + mS * u[1], mU0 + mS * u[2], v11 + vS * u[1]^2, v22 + vS * u[2]^2, v12 + vS * u[1] * u[2], idx)))
  }
  pi_ <- 10^p$lpi; tau2 <- if (cfg$tau == "free") exp(2 * p$ltau) else 0
  if (pi_ >= 1 || pi_ < 1e-6 || abs(p$beta) > 3 || tau2 > 4) return(-1e10)
  pG <- 1 - (1 - pi_)^g
  l0 <- log1p(-pG) + ld2(mD0, mU0, v11, v22, v12, idx)
  l1 <- log(pG) + ld2(mD0 + p$beta * u[1], mU0 + p$beta * u[2],
                      v11 + tau2 * u[1]^2, v22 + tau2 * u[2]^2, v12 + tau2 * u[1] * u[2], idx)
  sum(lse(l0, l1))
}

defaults <- list(mu_m = -0.003, lsm = log(0.05), lrho = log(RHO0), ltU = -Inf, theta = 0,
                 lpi = log10(0.05), beta = -0.2, ltau = log(0.1), mu_d = -0.005, lsd = log(0.02))

fit_one <- function(kind, cfg, fixed = defaults, idx = rep(TRUE, length(G)), nstart = 40) {
  nm <- pnames(cfg, kind)
  if (cfg$tauU == "free" && cfg$fit == "joint") fixed$ltU <- log(0.03)
  if (length(nm) == 0) return(list(par = fixed, ll = loglik(fixed, kind, cfg, idx), k = 0))
  f <- function(x) { v <- -loglik(build(as.list(x), nm, fixed), kind, cfg, idx); if (is.finite(v)) v else 1e10 }
  best <- NULL
  for (s in seq_len(nstart)) {
    x0 <- unlist(fixed[nm])
    jit <- c(mu_m = 0.01, lsm = 0.5, lrho = 0.3, ltU = 1, theta = pi, lpi = 1, beta = 0.15,
             ltau = 1, mu_d = 0.01, lsd = 0.7)[nm]
    if (s > 1) x0 <- x0 + rnorm(length(nm), 0, jit)
    if ("theta" %in% nm && s > 1) x0["theta"] <- runif(1, -pi, pi)
    if ("lpi" %in% nm && s > 1) x0["lpi"] <- runif(1, -3.5, -0.2)
    o <- list(par = x0, value = f(x0))
    for (r in 1:3) o <- optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-11, maxit = 5000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  list(par = build(as.list(best$par), nm, fixed), ll = -best$value, k = length(nm))
}

# summary of a fitted shift: vector, P(hit) at median G, perpendicular and along-line parts
describe <- function(p, kind, cfg) {
  rho <- exp(p$lrho)
  u <- if (cfg$dir == "free") c(cos(p$theta), sin(p$theta)) else c(1, rho)
  if (kind == "hit") { vec <- p$beta * u; pG <- 1 - (1 - 10^p$lpi)^Gmed; tau <- if (cfg$tau == "free") exp(p$ltau) else 0 }
  else if (kind == "cont") { vec <- Gmed * p$mu_d * u; pG <- NA; tau <- NA }
  else { vec <- c(0, 0); pG <- NA; tau <- NA }
  nrm <- c(rho, 1) / sqrt(1 + rho^2); al <- c(1, -rho) / sqrt(1 + rho^2)
  tibble(rho_hat = rho, tauU_hat = exp(p$ltU), shift_DEL = vec[1], shift_DUP = vec[2],
         perp = sum(vec * nrm), along = sum(vec * al),
         angle_deg = atan2(vec[2], vec[1]) * 180 / pi,
         P_hit_medG = pG, pi_hat = if (kind == "hit") 10^p$lpi else NA, tau_hat = tau)
}

##############################################
## All configurations
##############################################
grid <- expand.grid(dir = c("fixed", "free"), tau = c("free", "zero"), rho = c("fixed", "free"),
                    tauU = c("zero", "free"), fit = c("joint", "cond"), stringsAsFactors = FALSE)
rows <- list(); fits <- list()
for (i in seq_len(nrow(grid))) {
  cfg <- as.list(grid[i, ])
  fixed <- defaults
  kmir <- 0
  if (cfg$fit == "cond") {   # mirror on mirror-class loci only, then held fixed
    mc <- cfg; mc$fit <- "joint"
    m <- fit_one("mirror", mc, idx = is_mirror)
    fixed <- m$par
  }
  mir_all <- fit_one("mirror", cfg, fixed = fixed)
  for (kind in c("hit", "cont")) {
    if (kind == "cont" && cfg$tau == "zero") next   # tau does not apply to the continuous model
    f <- fit_one(kind, cfg, fixed = fixed)
    key <- paste(kind, paste(unlist(cfg), collapse = "_"), sep = ":")
    fits[[key]] <- list(f = f, cfg = cfg, kind = kind)
    rows[[key]] <- bind_cols(tibble(model = kind, !!!cfg, k = f$k, logL = f$ll, AIC = 2 * f$k - 2 * f$ll,
                                    logL_mirror_alone = mir_all$ll, k_mirror_alone = mir_all$k),
                             describe(f$par, kind, cfg))
  }
  cat(sprintf("[%2d/%d] %s\n", i, nrow(grid), paste(unlist(cfg), collapse = " ")))
}
res <- bind_rows(rows)
# delta AIC of each hit model vs the continuous model sharing its dir/rho/tauU/fit settings
cont_ref <- res %>% filter(model == "cont") %>% select(dir, rho, tauU, fit, AIC_cont = AIC)
res <- res %>% left_join(cont_ref, by = c("dir", "rho", "tauU", "fit")) %>%
  mutate(dAIC_hit_minus_cont = ifelse(model == "hit", AIC - AIC_cont, NA)) %>% select(-AIC_cont)
write_csv(res, "outputs/height_factor_grid.csv")
saveRDS(fits, "outputs/height_factor_fits.rds")

print(res %>% filter(model == "hit") %>%
        transmute(fit, dir, tau, rho, tauU, logL = round(logL, 2), dAIC_vs_cont = round(dAIC_hit_minus_cont, 2),
                  P_hit = round(P_hit_medG, 3), dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3),
                  perp = round(perp, 3), angle = round(angle_deg), tau_hat = round(tau_hat, 3), rho_hat = round(rho_hat, 3), tauU_hat = round(tauU_hat, 3)) %>%
        arrange(fit, dir, tau, rho, tauU), n = 50, width = 200)
print(res %>% filter(model == "cont") %>% transmute(fit, dir, rho, tauU, logL = round(logL, 2), dDEL = round(shift_DEL, 3),
                                                    dDUP = round(shift_DUP, 3), angle = round(angle_deg)), n = 50, width = 200)
