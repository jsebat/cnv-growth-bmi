#!/usr/bin/env Rscript
##############################################
## Molly continuous (u_d free) + single hit, no tauU, tau = 0: profile the HIT DIRECTION theta
## (pi on a grid, everything else refit). Asks: what does a DUP-down hit (Jon's geometry) cost?
## theta = -90 deg is pure DUP-down; -180 pure DEL-left. rho fixed 0.35 (not identified when u_d is free).
## Run: Rscript code/her_model_theta_profile.R   (from PGC_asym_model_diff/)
##############################################
src <- readLines("code/her_model_two_clusters.R")
eval(parse(text = src[seq_len(grep("^cfgs <- expand.grid", src) - 1)]))
RHO0 <- 0.35; cfg <- list(rho = "fixed", ud = "free", tau = "zero")
nth <- 6   # theta index in the parameter vector for this cfg
rows <- list()
for (thd in seq(-180, 0, by = 15)) {
  th <- thd * pi / 180; best_th <- NULL
  for (lp in seq(-3, -0.4, by = 0.2)) {
    f <- function(y) { x <- append(y, th, after = nth - 1); v <- -ll_model(unpack(x, cfg), lp); if (is.finite(v)) v else 1e10 }
    bb <- NULL
    for (s in 1:10) {
      y0 <- c(-0.003, log(0.05), -0.005, log(0.03), runif(1, -pi, pi), runif(1, 0.05, 0.5))
      o <- list(par = y0, value = f(y0))
      for (r in 1:3) o <- optim(o$par, f, method = "Nelder-Mead", control = list(reltol = 1e-11, maxit = 6000))
      if (is.null(bb) || o$value < bb$value) bb <- o
    }
    x <- append(bb$par, th, after = nth - 1)
    rows[[length(rows) + 1]] <- bind_cols(tibble(theta_deg = thd, log10_pi = lp, logL = -bb$value), summ(unpack(x, cfg), lp))
  }
  cat(sprintf("theta %d done\n", thd))
}
prof <- bind_rows(rows); write_csv(prof, "outputs/her_model_theta_profile.csv")
options(width = 200)
cat("\nBest over pi at each hit direction (k = 9 incl. theta; Molly continuous logL 15.76, AIC -23.51):\n")
print(prof %>% group_by(theta_deg) %>% slice_max(logL, n = 1) %>% ungroup() %>%
        transmute(theta_deg, logL = round(logL, 2), AIC = round(18 - 2 * logL, 2), P_hit = round(P_hit, 2),
                  dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3), sigma_d = round(sigma_d, 4), ud_ang = round(ud_angle)))
cat("\nAt theta = -90 (DUP-down), logL by P(hit):\n")
print(prof %>% filter(theta_deg == -90) %>% transmute(P_hit = round(P_hit, 2), logL = round(logL, 2), dDUP = round(shift_DUP, 3), sigma_d = round(sigma_d, 4)))
