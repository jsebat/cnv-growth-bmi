# Corrected bivariate mirror / continuous / single-hit CNV models
# Height asymmetry is constrained to the lower-left; BMI to the upper-right.

inv_logit <- function(x) 1 / (1 + exp(-x))
log_sum_exp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

bvn_logpdf <- function(x, y, m1, m2, v1, v2, cov12) {
  det <- v1 * v2 - cov12^2
  if (any(!is.finite(det)) || any(det <= 1e-14)) return(rep(-Inf, length(x)))
  dx <- x - m1
  dy <- y - m2
  quad <- (v2 * dx^2 + v1 * dy^2 - 2 * cov12 * dx * dy) / det
  -log(2 * pi) - 0.5 * log(det) - 0.5 * quad
}

check_input <- function(dat) {
  req <- c("Deletion_beta", "Duplication_beta", "Deletion_SE", "Duplication_SE", "Gene_count")
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  dat <- dat[complete.cases(dat[, req]), , drop = FALSE]
  dat$Gene_count <- pmax(1, as.numeric(dat$Gene_count))
  dat
}

mirror_components <- function(dat, base) {
  dat <- check_input(dat)
  G <- dat$Gene_count
  mean_gene <- base$p * base$sign_mean * base$m
  var_gene <- base$p * 2 * base$m^2 - (base$p * base$sign_mean * base$m)^2
  var_gene <- max(var_gene, 1e-12)
  latent_var <- G * var_gene
  list(
    D = dat$Deletion_beta,
    U = dat$Duplication_beta,
    G = G,
    mean_D = G * mean_gene,
    mean_U = -base$rho * G * mean_gene,
    var_D = latent_var + dat$Deletion_SE^2 + base$tauD^2,
    var_U = base$rho^2 * latent_var + dat$Duplication_SE^2 + base$tauU^2,
    cov_DU = -base$rho * latent_var
  )
}

fit_mirror <- function(dat, pi_positive, n_starts = 12, seed = 20260726) {
  dat <- check_input(dat)
  D <- dat$Deletion_beta
  U <- dat$Duplication_beta
  G <- dat$Gene_count
  seD <- dat$Deletion_SE
  seU <- dat$Duplication_SE
  sign_mean <- 2 * pi_positive - 1

  unpack <- function(par) list(
    p = inv_logit(par[1]),
    m = exp(par[2]),
    rho = exp(par[3]),
    tauD = exp(par[4]),
    tauU = exp(par[5])
  )
  objective <- function(par) {
    z <- unpack(par)
    mean_gene <- z$p * sign_mean * z$m
    var_gene <- z$p * 2 * z$m^2 - (z$p * sign_mean * z$m)^2
    if (!is.finite(var_gene) || var_gene <= 1e-12) return(1e100)
    latent_var <- G * var_gene
    ll <- bvn_logpdf(
      D, U,
      G * mean_gene,
      -z$rho * G * mean_gene,
      latent_var + seD^2 + z$tauD^2,
      z$rho^2 * latent_var + seU^2 + z$tauU^2,
      -z$rho * latent_var
    )
    if (any(!is.finite(ll))) return(1e100)
    -sum(ll)
  }

  p0 <- if (pi_positive < 0.5) 0.4 else 0.03
  m0 <- max(0.02, abs(mean(D)) / (max(abs(sign_mean), 0.1) * mean(G) * p0))
  slope0 <- tryCatch(abs(coef(lm(U ~ D))[2]), error = function(e) 0.5)
  rho0 <- min(1, max(0.15, slope0))
  start <- c(qlogis(p0), log(m0), log(rho0), log(0.02), log(0.05))

  set.seed(seed)
  best <- NULL
  for (i in seq_len(n_starts)) {
    init <- start + rnorm(5, 0, c(1.4, 1, 0.6, 1, 1))
    fit <- optim(init, objective, method = "Nelder-Mead", control = list(maxit = 12000, reltol = 1e-11))
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  z <- unpack(best$par)
  c(z, list(
    sign_mean = sign_mean,
    pi_positive = pi_positive,
    model = "Mirror alone",
    negative_log_likelihood = best$value,
    logLik = -best$value,
    n_parameters = 5,
    AIC = 10 + 2 * best$value,
    BIC = 5 * log(nrow(dat)) + 2 * best$value,
    raw_par = best$par
  ))
}

fit_directional <- function(dat, base, trait_sign,
                            model = c("mirror", "continuous", "single"),
                            max_q = NULL, n_starts = 16,
                            angle_min = 5, angle_max = 85,
                            seed = 20260726) {
  model <- match.arg(model)
  dat <- check_input(dat)
  cmp <- mirror_components(dat, base)
  if (model == "mirror") {
    nll <- -sum(bvn_logpdf(cmp$D, cmp$U, cmp$mean_D, cmp$mean_U,
                          cmp$var_D, cmp$var_U, cmp$cov_DU))
    return(list(
      model = "Mirror alone", q = 0, shift_magnitude = 0,
      shift_angle_deg = 45, shift_DEL = 0, shift_DUP = 0,
      negative_log_likelihood = nll, n_parameters_added = 0,
      conditional_AIC = 2 * nll, conditional_BIC = 2 * nll,
      max_q = NA_real_
    ))
  }
  if (is.null(max_q)) max_q <- if (model == "continuous") 0.8 else 0.2

  unpack <- function(par) {
    q <- max_q * inv_logit(par[1])
    magnitude <- exp(par[2])
    angle_deg <- angle_min + (angle_max - angle_min) * inv_logit(par[3])
    theta <- angle_deg * pi / 180
    list(q = q, magnitude = magnitude, angle_deg = angle_deg,
         shift_D = magnitude * cos(theta), shift_U = magnitude * sin(theta))
  }
  objective <- function(par) {
    z <- unpack(par)
    if (model == "continuous") {
      extra_var <- cmp$G * (2 * z$q - z$q^2)
      ll <- bvn_logpdf(
        cmp$D, cmp$U,
        cmp$mean_D + trait_sign * cmp$G * z$q * z$shift_D,
        cmp$mean_U + trait_sign * cmp$G * z$q * z$shift_U,
        cmp$var_D + extra_var * z$shift_D^2,
        cmp$var_U + extra_var * z$shift_U^2,
        cmp$cov_DU + extra_var * z$shift_D * z$shift_U
      )
      if (any(!is.finite(ll))) return(1e100)
      return(-sum(ll))
    }
    p0 <- (1 - z$q)^cmp$G
    ll0 <- log(pmax(p0, 1e-300)) +
      bvn_logpdf(cmp$D, cmp$U, cmp$mean_D, cmp$mean_U,
                 cmp$var_D, cmp$var_U, cmp$cov_DU)
    ll1 <- log(pmax(1 - p0, 1e-300)) +
      bvn_logpdf(cmp$D, cmp$U,
                 cmp$mean_D + trait_sign * z$shift_D,
                 cmp$mean_U + trait_sign * z$shift_U,
                 cmp$var_D, cmp$var_U, cmp$cov_DU)
    -sum(log_sum_exp2(ll0, ll1))
  }

  q0 <- if (model == "continuous") 0.08 else 0.018
  magnitude0 <- if (model == "continuous") 0.055 else 0.25
  angle0 <- 70
  start <- c(qlogis(q0 / max_q), log(magnitude0),
             qlogis((angle0 - angle_min) / (angle_max - angle_min)))
  set.seed(seed + ifelse(model == "continuous", 100, 200))
  best <- NULL
  for (i in seq_len(n_starts)) {
    init <- start + rnorm(3, 0, c(1.5, 1.2, 1.6))
    fit <- optim(init, objective, method = "Nelder-Mead",
                 control = list(maxit = 12000, reltol = 1e-11))
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  z <- unpack(best$par)
  label <- if (model == "continuous")
    "Mirror + continuous directional shift" else "Mirror + sparse single hit"
  list(
    model = label, q = z$q, shift_magnitude = z$magnitude,
    shift_angle_deg = z$angle_deg, shift_DEL = z$shift_D,
    shift_DUP = z$shift_U, negative_log_likelihood = best$value,
    n_parameters_added = 3,
    conditional_AIC = 6 + 2 * best$value,
    conditional_BIC = 3 * log(nrow(dat)) + 2 * best$value,
    max_q = max_q, raw_par = best$par
  )
}

model_locus_expectations <- function(dat, base, trait_sign, model, fit) {
  cmp <- mirror_components(dat, base)
  if (model == "mirror") {
    expected_D <- cmp$mean_D
    expected_U <- cmp$mean_U
    # 0, not NA: the mirror model has no single-hit component, matching the
    # Python reference. NA here silently broke anything summing this column.
    hit_prob <- rep(0, nrow(dat))
  } else if (model == "continuous") {
    expected_D <- cmp$mean_D + trait_sign * cmp$G * fit$q * fit$shift_DEL
    expected_U <- cmp$mean_U + trait_sign * cmp$G * fit$q * fit$shift_DUP
    hit_prob <- rep(NA_real_, nrow(dat))
  } else {
    hit_prob <- 1 - (1 - fit$q)^cmp$G
    expected_D <- cmp$mean_D + trait_sign * hit_prob * fit$shift_DEL
    expected_U <- cmp$mean_U + trait_sign * hit_prob * fit$shift_DUP
  }
  cbind(dat,
        model_expected_DEL = expected_D,
        model_expected_DUP = expected_U,
        mirror_expected_DEL = cmp$mean_D,
        mirror_expected_DUP = cmp$mean_U,
        single_hit_probability = hit_prob,
        residual_DEL_from_mirror = dat$Deletion_beta - cmp$mean_D,
        residual_DUP_from_mirror = dat$Duplication_beta - cmp$mean_U,
        directional_shift_score = trait_sign *
          ((dat$Deletion_beta - cmp$mean_D) + (dat$Duplication_beta - cmp$mean_U)) / sqrt(2))
}

burden_values <- function(base, trait_sign, model, fit, gene_counts) {
  d0 <- base$p * base$sign_mean * base$m
  u0 <- -base$rho * d0
  G <- as.numeric(gene_counts)
  if (model == "mirror") {
    d <- rep(d0, length(G)); u <- rep(u0, length(G))
  } else if (model == "continuous") {
    d <- rep(d0 + trait_sign * fit$q * fit$shift_DEL, length(G))
    u <- rep(u0 + trait_sign * fit$q * fit$shift_DUP, length(G))
  } else {
    occupancy <- 1 - (1 - fit$q)^G
    d <- d0 + trait_sign * fit$shift_DEL * occupancy / G
    u <- u0 + trait_sign * fit$shift_DUP * occupancy / G
  }
  data.frame(DEL_SD_per_gene = d, DUP_SD_per_gene = u)
}

conditional_loo <- function(dat, base, trait_sign, fits, seed = 20260726) {
  rows <- list(); k <- 1
  for (model in c("mirror", "continuous", "single")) {
    vals <- numeric(nrow(dat))
    for (i in seq_len(nrow(dat))) {
      train <- dat[-i, , drop = FALSE]
      test <- dat[i, , drop = FALSE]
      fit <- if (model == "mirror") fits[[model]] else
        fit_directional(train, base, trait_sign, model,
                        max_q = ifelse(model == "continuous", 0.8, 0.2),
                        n_starts = 4, seed = seed + i)
      cmp <- mirror_components(test, base)
      if (model == "mirror") {
        vals[i] <- bvn_logpdf(cmp$D, cmp$U, cmp$mean_D, cmp$mean_U,
                              cmp$var_D, cmp$var_U, cmp$cov_DU)
      } else if (model == "continuous") {
        extra_var <- cmp$G * (2 * fit$q - fit$q^2)
        vals[i] <- bvn_logpdf(
          cmp$D, cmp$U,
          cmp$mean_D + trait_sign * cmp$G * fit$q * fit$shift_DEL,
          cmp$mean_U + trait_sign * cmp$G * fit$q * fit$shift_DUP,
          cmp$var_D + extra_var * fit$shift_DEL^2,
          cmp$var_U + extra_var * fit$shift_DUP^2,
          cmp$cov_DU + extra_var * fit$shift_DEL * fit$shift_DUP
        )
      } else {
        p0 <- (1 - fit$q)^cmp$G
        l0 <- log(pmax(p0, 1e-300)) + bvn_logpdf(
          cmp$D, cmp$U, cmp$mean_D, cmp$mean_U, cmp$var_D, cmp$var_U, cmp$cov_DU)
        l1 <- log(pmax(1 - p0, 1e-300)) + bvn_logpdf(
          cmp$D, cmp$U,
          cmp$mean_D + trait_sign * fit$shift_DEL,
          cmp$mean_U + trait_sign * fit$shift_DUP,
          cmp$var_D, cmp$var_U, cmp$cov_DU)
        vals[i] <- log_sum_exp2(l0, l1)
      }
      # Label held-out folds by locus, not by row index: the integer index was
      # meaningless outside this run and made the output impossible to join
      # against the Python reference.
      locus <- if ("CNV_ID" %in% names(dat)) as.character(dat$CNV_ID[i]) else as.character(i)
      rows[[k]] <- data.frame(model = fits[[model]]$model, held_out_locus = locus,
                              held_out_log_predictive_density = vals[i],
                              stringsAsFactors = FALSE); k <- k + 1
    }
  }
  detail <- do.call(rbind, rows)
  summary <- aggregate(held_out_log_predictive_density ~ model, detail, sum)
  names(summary)[2] <- "conditional_LOO_log_predictive_density"
  list(detail = detail, summary = summary)
}

fit_trait_models <- function(dat, mirror_dat, pi_positive, trait_sign,
                             seed = 20260726, run_loo = TRUE) {
  base <- fit_mirror(mirror_dat, pi_positive, seed = seed)
  fits <- list(
    mirror = fit_directional(dat, base, trait_sign, "mirror"),
    continuous = fit_directional(dat, base, trait_sign, "continuous", 0.8, seed = seed),
    single = fit_directional(dat, base, trait_sign, "single", 0.2, seed = seed)
  )
  comparison <- do.call(rbind, lapply(fits, function(x) {
    keep <- c("model", "q", "shift_magnitude", "shift_angle_deg", "shift_DEL", "shift_DUP",
              "negative_log_likelihood", "n_parameters_added", "conditional_AIC", "conditional_BIC",
              "max_q")
    as.data.frame(x[keep], stringsAsFactors = FALSE)
  }))
  rownames(comparison) <- NULL
  loo <- NULL
  if (run_loo) {
    loo <- conditional_loo(dat, base, trait_sign, fits, seed)
    comparison <- merge(comparison, loo$summary, by = "model", all.x = TRUE, sort = FALSE)
  }
  comparison$delta_conditional_AIC <- comparison$conditional_AIC - min(comparison$conditional_AIC)
  comparison$delta_conditional_BIC <- comparison$conditional_BIC - min(comparison$conditional_BIC)
  if (run_loo) comparison$delta_conditional_LOO <-
    max(comparison$conditional_LOO_log_predictive_density) - comparison$conditional_LOO_log_predictive_density
  list(base = base, fits = fits, comparison = comparison, loo = loo)
}
