#!/usr/bin/env Rscript
##############################################
## Power: with 34 height loci (real gene counts and SEs), how often does AIC pick a two-cluster model
## over the continuous model when clusters are truly there?  And how often when they are not?
##
## Analysis models (fit to every simulated dataset, no tauU in either):
##   CONT (k = 4): Molly's continuous mirror + continuous deleterious along (1, rho), rho = log2(3/2)
##   HIT  (k = 7): CONT + a fixed-size CNV-level single hit in a free direction (P = 1 - (1 - pi)^G)
## Generating scenarios:
##   A  "two clusters, as fitted":   HIT at its MLE on the real data (P(hit|G=23) ~0.44, shift ~0.13 SD)
##   B  "Jon geometry":              Jon-type single hit fitted jointly (rho 0.31, tauU 0.05, P 0.35, shift (0.25, -0.36))
##   C  "no clusters" (null):        CONT at its MLE on the real data
## Each dataset: true effects drawn from the scenario per locus, plus N(0, SE^2) sampling noise.
## Run: Rscript code/power_simulation.R [n_rep]   (from PGC_asym_model_diff/)
##############################################
suppressMessages({library(dplyr); library(readr); library(parallel)})
N_REP <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(N_REP)) N_REP <- 300
dat <- read_csv("inputs/height_loci.csv", show_col_types = FALSE)
G <- dat$Gene_count; vD <- dat$Deletion_SE^2; vU <- dat$Duplication_SE^2
RHO0 <- log2(3 / 2)

ld2 <- function(D, U, mD, mU, v11, v22, v12) {
  det <- v11 * v22 - v12^2; dx <- D - mD; dy <- U - mU
  -log(2 * pi) - 0.5 * log(det) - 0.5 * (v22 * dx^2 + v11 * dy^2 - 2 * v12 * dx * dy) / det
}
lse <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
# continuous part shared by both models: x[1:4] = mu_m, log s_m, mu_d, log s_d
cont_moments <- function(x) {
  sm2 <- exp(2 * x[2]); sd2 <- exp(2 * x[4]); r <- RHO0
  list(mD = G * (x[1] + x[3]), mU = G * r * (x[3] - x[1]),
       v11 = G * (sm2 + sd2) + vD, v22 = G * r^2 * (sm2 + sd2) + vU, v12 = G * r * (sd2 - sm2), ok = sm2 < 1 && sd2 < 1)
}
ll_cont <- function(x, D, U) { m <- cont_moments(x); if (!m$ok) return(-1e10); sum(ld2(D, U, m$mD, m$mU, m$v11, m$v22, m$v12)) }
# x[5:7] = theta, log10 pi, beta
ll_hit <- function(x, D, U) {
  m <- cont_moments(x); pi_ <- 10^x[6]
  if (!m$ok || pi_ >= 1 || pi_ < 1e-6 || abs(x[7]) > 3) return(-1e10)
  pG <- 1 - (1 - pi_)^G; u <- c(cos(x[5]), sin(x[5]))
  sum(lse(log1p(-pG) + ld2(D, U, m$mD, m$mU, m$v11, m$v22, m$v12),
          log(pG) + ld2(D, U, m$mD + x[7] * u[1], m$mU + x[7] * u[2], m$v11, m$v22, m$v12)))
}
nm <- function(f, x0, rounds = 2) {
  o <- list(par = x0, value = f(x0))
  for (r in seq_len(rounds)) o <- optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-10, maxit = 4000))
  o
}
fit_both <- function(D, U, hit_start, nstart = 8) {
  fc <- function(x) { v <- -ll_cont(x, D, U); if (is.finite(v)) v else 1e10 }
  c0 <- nm(fc, c(-0.003, log(0.05), -0.005, log(0.02)), 3)
  fh <- function(x) { v <- -ll_hit(x, D, U); if (is.finite(v)) v else 1e10 }
  starts <- list(c(c0$par, 0, -5, 0.01), hit_start)   # first start = continuous fit (pi ~ 0) so HIT >= CONT
  for (s in seq_len(nstart)) starts[[length(starts) + 1]] <- c(c0$par + rnorm(4, 0, 0.2), runif(1, -pi, pi), runif(1, -3, -0.5), runif(1, 0.05, 0.5))
  best <- NULL
  for (x0 in starts) { o <- nm(fh, x0); if (is.null(best) || o$value < best$value) best <- o }
  pG <- 1 - (1 - 10^best$par[6])^23; vec <- best$par[7] * c(cos(best$par[5]), sin(best$par[5]))
  c(logL_cont = -c0$value, logL_hit = -best$value, P_hit = pG, shift_DEL = vec[1], shift_DUP = vec[2])
}

## generators
real_D <- dat$Deletion_beta; real_U <- dat$Duplication_beta
set.seed(11)
fc <- function(x) -ll_cont(x, real_D, real_U); cont_mle <- nm(fc, c(-0.003, log(0.05), -0.005, log(0.02)), 3)$par
fh <- function(x) { v <- -ll_hit(x, real_D, real_U); if (is.finite(v)) v else 1e10 }
hit_mle <- NULL
for (s in 1:40) { o <- nm(fh, c(cont_mle + rnorm(4, 0, 0.2), runif(1, -pi, pi), runif(1, -3, -0.5), runif(1, 0.05, 0.5)), 3)
                  if (is.null(hit_mle) || o$value < hit_mle$value) hit_mle <- o }
hit_mle <- hit_mle$par
jon <- readRDS("outputs/height_factor_fits.rds")[["hit:free_zero_free_free_joint"]]$f$par
cat(sprintf("Real data: CONT logL %.2f (AIC %.2f), HIT logL %.2f (AIC %.2f)\n", -fc(cont_mle), 8 + 2 * fc(cont_mle), -fh(hit_mle), 14 + 2 * fh(hit_mle)))

gen_true <- function(scn) {
  n <- length(G)
  if (scn == "C" || scn == "A") {
    x <- if (scn == "C") cont_mle else hit_mle
    M <- rnorm(n, G * x[1], sqrt(G) * exp(x[2])); S <- rnorm(n, G * x[3], sqrt(G) * exp(x[4]))
    D <- M + S; U <- RHO0 * (S - M)
    if (scn == "A") { hit <- runif(n) < 1 - (1 - 10^x[6])^G; D <- D + hit * x[7] * cos(x[5]); U <- U + hit * x[7] * sin(x[5]) }
  } else {
    p <- jon; rho <- exp(p$lrho)
    M <- rnorm(n, G * p$mu_m, sqrt(G) * exp(p$lsm)); hit <- runif(n) < 1 - (1 - 10^p$lpi)^G
    D <- M + hit * p$beta * cos(p$theta); U <- -rho * M + rnorm(n, 0, exp(p$ltU)) + hit * p$beta * sin(p$theta)
  }
  list(D = D + rnorm(n, 0, sqrt(vD)), U = U + rnorm(n, 0, sqrt(vU)))
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("G", "vD", "vU", "RHO0", "ld2", "lse", "cont_moments", "ll_cont", "ll_hit", "nm", "fit_both",
                    "gen_true", "cont_mle", "hit_mle", "jon"))
res <- list()
for (scn in c("A", "B", "C")) {
  t0 <- Sys.time()
  out <- parLapply(cl, seq_len(N_REP), function(i, scn) {
    set.seed(1e5 * match(scn, c("A", "B", "C")) + i); y <- gen_true(scn)
    fit_both(y$D, y$U, hit_mle)
  }, scn = scn)
  res[[scn]] <- as_tibble(do.call(rbind, out)) %>% mutate(scenario = scn, rep = row_number())
  cat(sprintf("scenario %s: %d reps in %.1f min\n", scn, N_REP, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
sims <- bind_rows(res) %>% mutate(dAIC = (14 - 2 * logL_hit) - (8 - 2 * logL_cont))   # negative favours HIT
write_csv(sims, "outputs/power_simulation_reps.csv")
summ <- sims %>% group_by(scenario) %>%
  summarise(n = n(), picks_HIT = mean(dAIC < 0), HIT_by_2 = mean(dAIC < -2), HIT_by_4 = mean(dAIC < -4),
            median_dAIC = median(dAIC), q10 = quantile(dAIC, 0.1), q90 = quantile(dAIC, 0.9),
            P_hit_recovered_median = median(P_hit[dAIC < 0]), .groups = "drop") %>%
  mutate(scenario = recode(scenario, A = "A: two clusters as fitted to real data", B = "B: Jon geometry (rho 0.31, tauU, 0.3 SD shift)",
                           C = "C: no clusters (Molly continuous)"))
real_dAIC <- (14 + 2 * fh(hit_mle)) - (8 + 2 * fc(cont_mle))
summ <- summ %>% mutate(real_data_dAIC = real_dAIC)
write_csv(summ, "outputs/power_simulation_summary.csv")
options(width = 200); print(summ)
