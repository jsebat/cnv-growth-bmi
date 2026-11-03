#!/usr/bin/env Rscript
##############################################
## Follow-up to her_model_two_clusters.R: in Molly's model + single hit, WITHOUT tauU, profile rho
## (u_d free, tau = 0; pi profiled on a grid, everything else refit) to see what the Jon-like geometry
## (rho ~0.3-0.35, DUP-down hit) costs, then draw latent and noisy clouds at selected fits.
## Run: Rscript code/her_model_rho_profile.R   (from PGC_asym_model_diff/)
##############################################
src <- readLines("code/her_model_two_clusters.R")
eval(parse(text = src[seq_len(grep("^cfgs <- expand.grid", src) - 1)]))

cfg <- list(rho = "fixed", ud = "free", tau = "zero")
rhos <- c(0.25, 0.30, 0.35, 0.40, 0.50, 0.585, 0.65, 0.75, 0.85)
lpis <- seq(-3, -0.4, by = 0.2)
rows <- list(); pars <- list()
for (r in rhos) {
  RHO0 <- r; start <- NULL
  for (lp in lpis) {
    b <- fit_at_pi(cfg, lp, start, nstart = 15); start <- b$par
    rows[[length(rows) + 1]] <- bind_cols(tibble(rho_fixed = r, log10_pi = lp, logL = -b$value),
                                          summ(unpack(b$par, cfg), lp))
    pars[[sprintf("%.3f_%.1f", r, lp)]] <- b$par
  }
  cat(sprintf("rho %.3f done\n", r))
}
prof <- bind_rows(rows)
write_csv(prof, "outputs/her_model_rho_profile.csv"); saveRDS(pars, "outputs/her_model_rho_profile_pars.rds")
options(width = 200)
# k = 8 when rho is counted as fitted (mu_m, s_m, rho, mu_d, s_d, phi, theta, beta) + pi = 9
best <- prof %>% group_by(rho_fixed) %>% slice_max(logL, n = 1) %>% ungroup() %>% mutate(AIC_if_rho_fitted = 2 * 9 - 2 * logL)
cat("\nBest over pi at each fixed rho (Molly continuous: logL 15.76, AIC -23.51):\n")
print(best %>% transmute(rho_fixed, logL = round(logL, 2), AIC9 = round(AIC_if_rho_fitted, 2), P_hit = round(P_hit, 2),
                         dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3), perp = round(perp, 3),
                         sigma_d = round(sigma_d, 4), ud_ang = round(ud_angle)))

## clouds: best fit at rho 0.35, 0.585, and the overall best rho
sim_lat <- function(p, lpi, n = 4000, noise = FALSE, seed = 3) {
  set.seed(seed); i <- sample(seq_along(G), n, TRUE); g <- G[i]
  ud <- c(cos(p$phi), sin(p$phi)); uh <- c(cos(p$theta), sin(p$theta)); a <- c(1, -p$rho)
  M <- rnorm(n, g * p$mu_m, sqrt(g * p$sm2)); S <- rnorm(n, g * p$mu_d, sqrt(g * p$sd2))
  hit <- runif(n) < 1 - (1 - 10^lpi)^g
  del <- M * a[1] + S * ud[1] + hit * p$beta * uh[1]; dup <- M * a[2] + S * ud[2] + hit * p$beta * uh[2]
  if (noise) { del <- del + rnorm(n, 0, sqrt(vD[i])); dup <- dup + rnorm(n, 0, sqrt(vU[i])) }
  tibble(DEL = del, DUP = dup, hit = hit)
}
obs <- dat %>% transmute(DEL = Deletion_beta, DUP = Duplication_beta, class = Dose_response_class)
pick <- unique(c(0.35, 0.585, best$rho_fixed[which.max(best$logL)]))
panels <- list()
for (r in pick) {
  bb <- best %>% filter(rho_fixed == r); RHO0 <- r
  p <- unpack(pars[[sprintf("%.3f_%.1f", r, bb$log10_pi)]], cfg)
  for (nz in c(FALSE, TRUE)) {
    s <- sim_lat(p, bb$log10_pi, noise = nz)
    panels[[length(panels) + 1]] <- ggplot(s, aes(DEL, DUP)) +
      geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
      geom_point(aes(colour = ifelse(hit, "hit", "not hit")), size = 0.35, alpha = 0.25) +
      geom_abline(slope = -r, intercept = 0, linetype = "dashed", colour = "#7b3294") +
      geom_point(data = obs, aes(shape = class), colour = "#ff7f0e", size = 1.6) +
      scale_colour_manual(values = c(hit = "#d62728", `not hit` = "#1f77b4"), name = NULL) +
      scale_shape_manual(values = c(asymmetric = 17, mirror = 16, unclassified = 4), name = NULL) +
      coord_cartesian(xlim = c(-1.4, 0.6), ylim = c(-0.8, 0.6)) +
      labs(title = if (!nz) sprintf("Molly continuous + single hit, no tauU, rho = %.3f", r) else NULL,
           subtitle = sprintf("%s; logL %.1f, P(hit|G=23) %.2f, shift (%.2f, %.2f)", if (nz) "with locus SE noise" else "latent",
                              bb$logL, bb$P_hit, bb$shift_DEL, bb$shift_DUP),
           x = "DEL effect (SD)", y = "DUP effect (SD)") +
      theme_bw(base_size = 8) + theme(plot.title = element_text(size = 8, face = "bold"))
  }
}
fig <- wrap_plots(panels, ncol = 2, byrow = TRUE) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("figures/Figure_her_model_no_tauU_clouds.png", fig, width = 9, height = 4 * length(pick), dpi = 200)
ggsave("figures/Figure_her_model_no_tauU_clouds.pdf", fig, width = 9, height = 4 * length(pick))
