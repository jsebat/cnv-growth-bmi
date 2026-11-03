#!/usr/bin/env Rscript
##############################################
## Common-scale comparison of Molly's and Jon's height models + noise-free simulated clouds.
## Every model: same 34 loci, same likelihood, all parameters (mirror included) fit jointly,
## AIC counts all parameters, plus leave-one-locus-out log predictive density (full refit per locus).
## Uses the likelihood and fitter from factor_test.R (the part above "All configurations").
## Run: Rscript code/loo_and_clouds.R   (from PGC_asym_model_diff/)
##############################################
src <- readLines("code/factor_test.R")
eval(parse(text = src[seq_len(grep("^## All configurations", src) - 2)]))
fits <- readRDS("outputs/height_factor_fits.rds")

key <- function(kind, dir, tau, rho, tauU) paste0(kind, ":", paste(dir, tau, rho, tauU, "joint", sep = "_"))
models <- tribble(
  ~label,                                                       ~key,
  "Molly continuous (dir (1,rho), rho 0.585)",                   key("cont", "fixed", "free", "fixed", "zero"),
  "Molly CNV hit (dir (1,rho), rho 0.585, tau free)",            key("hit",  "fixed", "free", "fixed", "zero"),
  "Continuous, free dir, rho free, tauU free",                   key("cont", "free",  "free", "free",  "free"),
  "Jon single hit, joint (free dir, rho free, tauU, tau = 0)",   key("hit",  "free",  "zero", "free",  "free"),
  "Jon hit but rho fixed 0.585",                                 key("hit",  "free",  "zero", "fixed", "free"),
  "Jon hit but dir fixed (1,rho)",                               key("hit",  "fixed", "zero", "free",  "free"),
  "Jon hit but tau free",                                        key("hit",  "free",  "free", "free",  "free"),
  "Mirror alone, rho 0.585",                                     "mirror:fixed_zero",
  "Mirror alone, rho free, tauU free",                           "mirror:free_free")

get_fit <- function(k) {
  if (startsWith(k, "mirror:")) {
    r <- strsplit(sub("mirror:", "", k), "_")[[1]]
    cfg <- list(dir = "fixed", tau = "free", rho = r[1], tauU = r[2], fit = "joint")
    return(list(f = fit_one("mirror", cfg), cfg = cfg, kind = "mirror"))
  }
  fits[[k]]
}

loo <- function(m) {
  n <- length(G); lp <- numeric(n)
  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    fi <- fit_one(m$kind, m$cfg, fixed = m$f$par, idx = keep, nstart = 6)   # start at full-data MLE
    lp[i] <- loglik(fi$par, m$kind, m$cfg, idx = seq_len(n) == i)
  }
  lp
}

out <- list(); lpm <- list()
for (j in seq_len(nrow(models))) {
  m <- get_fit(models$key[j]); lp <- loo(m); lpm[[models$label[j]]] <- lp
  out[[j]] <- tibble(model = models$label[j], k = m$f$k, logL = m$f$ll, AIC = 2 * m$f$k - 2 * m$f$ll,
                     LOO_lpd = sum(lp)) %>% bind_cols(describe(m$f$par, m$kind, m$cfg))
  cat(sprintf("[%d/%d] %s  LOO = %.2f\n", j, nrow(models), models$label[j], sum(lp)))
}
cmp <- bind_rows(out) %>% mutate(dAIC = AIC - min(AIC), dLOO = max(LOO_lpd) - LOO_lpd) %>% arrange(AIC)
write_csv(cmp, "outputs/height_common_scale_comparison.csv")
write_csv(as_tibble(lpm) %>% mutate(CNV_ID = dat$CNV_ID, class = dat$Dose_response_class, .before = 1),
          "outputs/height_loo_pointwise.csv")
options(width = 220)
print(cmp %>% transmute(model, k, logL = round(logL, 2), AIC = round(AIC, 2), dAIC = round(dAIC, 2),
                        LOO = round(LOO_lpd, 2), dLOO = round(dLOO, 2), P_hit = round(P_hit_medG, 2),
                        dDEL = round(shift_DEL, 3), dDUP = round(shift_DUP, 3), rho = round(rho_hat, 3)))

##############################################
## Noise-free (latent) and noisy simulated clouds at each MLE
##############################################
simulate <- function(m, n = 4000, noise = FALSE, seed = 7) {
  set.seed(seed); p <- m$f$par; cfg <- m$cfg
  i <- sample(seq_along(G), n, TRUE); g <- G[i]; rho <- exp(p$lrho)
  Mv <- rnorm(n, g * p$mu_m, sqrt(g) * exp(p$lsm))
  u <- if (cfg$dir == "free") c(cos(p$theta), sin(p$theta)) else c(1, rho)
  S <- if (m$kind == "hit") {
    hit <- runif(n) < 1 - (1 - 10^p$lpi)^g
    hit * rnorm(n, p$beta, if (cfg$tau == "free") exp(p$ltau) else 0)
  } else if (m$kind == "cont") rnorm(n, g * p$mu_d, sqrt(g) * exp(p$lsd)) else 0
  del <- Mv + S * u[1]; dup <- -rho * Mv + S * u[2] + rnorm(n, 0, exp(p$ltU))
  if (noise) { del <- del + rnorm(n, 0, sqrt(vD[i])); dup <- dup + rnorm(n, 0, sqrt(vU[i])) }
  tibble(DEL = del, DUP = dup, hit = if (m$kind == "hit") hit else NA, rho = rho)
}
show <- models$label[c(2, 5, 6, 7, 4, 3)]
obs <- dat %>% transmute(DEL = Deletion_beta, DUP = Duplication_beta, class = Dose_response_class)
panels <- list()
for (lab in show) for (noise in c(FALSE, TRUE)) {
  m <- get_fit(models$key[models$label == lab]); s <- simulate(m, noise = noise)
  r <- cmp %>% filter(model == lab)
  pp <- ggplot(s, aes(DEL, DUP)) +
    geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_point(aes(colour = if (all(is.na(hit))) "continuous" else ifelse(hit, "hit", "not hit")),
               size = 0.35, alpha = 0.25) +
    geom_abline(slope = -s$rho[1], intercept = 0, linetype = "dashed", colour = "#7b3294") +
    geom_point(data = obs, aes(shape = class), colour = "#ff7f0e", size = 1.6) +
    scale_colour_manual(values = c(continuous = "#1f77b4", hit = "#d62728", `not hit` = "#1f77b4"), name = NULL) +
    scale_shape_manual(values = c(asymmetric = 17, mirror = 16, unclassified = 4), name = NULL) +
    coord_cartesian(xlim = c(-1.4, 0.6), ylim = c(-0.8, 0.6)) +
    labs(title = if (!noise) lab else NULL,
         subtitle = sprintf("%s;  AIC %.1f, LOO %.1f%s", if (noise) "with locus SE noise" else "latent (no SE noise)",
                            r$AIC, r$LOO_lpd, if (!is.na(r$P_hit_medG)) sprintf(", P(hit|G=23) %.2f", r$P_hit_medG) else ""),
         x = "DEL effect (SD)", y = "DUP effect (SD)") +
    theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 8, face = "bold"))
  panels[[paste(lab, noise)]] <- pp
}
fig <- wrap_plots(panels, ncol = 4, byrow = TRUE) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/Figure_height_factor_clouds.png", fig, width = 15, height = 11, dpi = 200)
ggsave("figures/Figure_height_factor_clouds.pdf", fig, width = 15, height = 11)
