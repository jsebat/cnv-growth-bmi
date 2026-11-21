# make_synthetic.R  (2026-11-20)
# What: SYNTHETIC test sites for the federated package, with known effects. Nothing here is real: ids are
#   FAKE_<site>_<number>, families FAKEFAM_<site>_<number>, genotypes and measurements are simulated.
#     site A  mixed ages (children and adults), controls with adults, SRS and SCQ, autism diagnosis
#     site B  paediatric clinic: carriers and age-matched controls from infancy to 18 y, few adults, repeated visits,
#             SRS only; its 22q11.2 dup group has exactly 3 persons aged 10-18 (a small cell that must be suppressed)
#   Simulated truth (Z units), returned as TRUTH:
#     height     Z = main(genotype, age) + 0.15 x BMI Z + 0.40 x BMI Z for 16p11.2 del before 10 y (the CNV x BMI effect)
#                22q11.2 del main effect rises with age: -1.2 at 2 y to -0.7 from 19 y (the age-varying trajectory)
#     head circ. Z = main(genotype) + 0.15 x BMI Z; in 22q11.2 dup children each next visit adds 0.12 x BMI Z
#                (the next-visit effect); no next-visit effect in the other genotypes
#     diagnosis  log odds of an autism diagnosis in carriers rise 0.4 per SD of BMI Z
#     SRS, SCQ   trait Z = 0.8 x carrier + (0.05 + 0.15 x carrier) x BMI Z + noise (CNV x BMI = 0.15)
# Inputs: growth_refs_lms.csv of the package (to turn simulated Z scores into cm and kg).
# Outputs (in the folder given as argument; default a temporary folder, never the repository): visits_<site>.tsv,
#   config_<site>.json. make_all() returns the paths and TRUTH.
# null model: none; this script only simulates data.
# Usage: Rscript tests/make_synthetic.R [out_dir]      (also sourced by tests/test_federated.R)
TRUTH <- list(ctl_slope = 0.15, cnvxbmi_height_16pdel_2_9 = 0.40, next_visit_hc_22qdup = 0.12, asd_logor_bmi = 0.4, trait_cnvxbmi = 0.15,
              main_height = c(`No CNV` = 0, `16p11.2 del` = -0.4, `16p11.2 dup` = -0.2, `22q11.2 del` = NA, `22q11.2 dup` = 0),
              main_hc = c(`No CNV` = 0, `16p11.2 del` = 0.7, `16p11.2 dup` = -0.5, `22q11.2 del` = -0.6, `22q11.2 dup` = 0.1),
              mean_bmi = c(`No CNV` = 0, `16p11.2 del` = 0.8, `16p11.2 dup` = -0.4, `22q11.2 del` = 0, `22q11.2 dup` = 0.2))
truth_height_22qdel <- function(age) -1.2 + 0.5 * pmin(1, log(age / 2) / log(19 / 2))
ADULT_L <- list(height = c(male = 0.65, female = 0.60), bmi = c(male = -0.75, female = -1.05))

.inv_lms <- function(z, age, sex, ref, item, LMS) {   # Z -> measurement with the same table and interpolation as the scorer
  y <- rep(NA_real_, length(z))
  for (s in c("male", "female")) { t <- LMS[LMS$ref == ref & LMS$item == item & LMS$sex == s, ]; t <- t[!is.na(t$mu) & !is.na(t$sigma) & !is.na(t$nu), ]; t <- t[order(t$age), ]
    i <- which(sex == s); if (!length(i)) next
    a <- pmin(pmax(age[i], min(t$age)), max(t$age)); L <- approx(t$age, t$nu, a)$y; M <- approx(t$age, t$mu, a)$y; S <- approx(t$age, t$sigma, a)$y
    y[i] <- ifelse(abs(L) < 1e-8, M * exp(S * z[i]), M * pmax(1e-6, 1 + L * S * z[i])^(1 / L)) }
  y
}
.inv_bc <- function(t, L) (1 + L * t)^(1 / L)
.who <- function(z, age, sex, item, LMS) ifelse(age < 5, .inv_lms(z, age, sex, "who", item, LMS), .inv_lms(z, age, sex, "who2007", item, LMS))

# n: named vector of persons per genotype; ages: function(n) giving the age at the first visit; visits: max visits of a child
make_site <- function(label, n, ages, visits, seed, LMS, scq = TRUE, small_cell = FALSE) {
  set.seed(seed); N <- sum(n); g <- rep(names(n), n); car <- g != "No CNV"
  id <- sprintf("FAKE_%s_%06d", label, seq_len(N)); sex <- sample(c("male", "female"), N, TRUE)
  fam <- sprintf("FAKEFAM_%s_%06d", label, seq_len(N)); k <- which(!car); pair <- k[seq_len(2 * (length(k) %/% 4))]   # half of the controls come in pairs
  fam[pair[seq(2, length(pair), 2)]] <- fam[pair[seq(1, length(pair), 2)]]; fe <- rnorm(N, 0, 0.5); fe <- fe[match(fam, fam)]
  a0 <- ages(N)
  if (small_cell) { j <- which(g == "22q11.2 dup"); a0[j] <- runif(length(j), 2.2, 4.5); a0[j[1:3]] <- c(11.5, 13.2, 15.4) }
  b <- rnorm(N, TRUTH$mean_bmi[g], 1); ph <- sqrt(1 - 0.15^2 - 0.25) * rnorm(N) + fe; pc <- rnorm(N, 0, sqrt(1 - 0.15^2))
  asd <- ifelse(car, rbinom(N, 1, plogis(-1 + TRUTH$asd_logor_bmi * b)), 0L)
  tz <- 0.8 * car + (0.05 + TRUTH$trait_cnvxbmi * car) * b + rnorm(N)
  srs <- pmax(0, round(.inv_bc(5.15 + 1.5 * tz, 0.32) - 1)); tz2 <- 0.8 * car + (0.05 + TRUTH$trait_cnvxbmi * car) * b + rnorm(N)
  scqv <- pmax(0, round(.inv_bc(1.42 + 0.46 * pmin(tz2, 8), -0.16) - 1))
  rows <- vector("list", N); nv <- ifelse(a0 >= 18, 1L, sample.int(visits, N, TRUE))
  for (i in seq_len(N)) { m <- nv[i]; age <- a0[i] + c(0, cumsum(runif(m - 1, 0.8, 1.3))); age <- age[age < 73]; m <- length(age); bt <- b[i] + rnorm(m, 0, 0.1)
    mh <- if (g[i] == "22q11.2 del") truth_height_22qdel(pmax(age, 2)) else TRUTH$main_height[[g[i]]]
    hz <- mh + TRUTH$ctl_slope * bt + (g[i] == "16p11.2 del" & age < 10) * TRUTH$cnvxbmi_height_16pdel_2_9 * bt + ph[i] + rnorm(m, 0, 0.08)
    cz <- TRUTH$main_hc[[g[i]]] + TRUTH$ctl_slope * b[i] + pc[i] + c(0, cumsum(rnorm(m - 1, 0, 0.1) + (g[i] == "22q11.2 dup") * TRUTH$next_visit_hc_22qdup * head(bt, -1)))
    rows[[i]] <- data.frame(i = i, age = age, hz = hz, bz = bt, cz = cz) }
  V <- do.call(rbind, rows); s <- sex[V$i]; ad <- V$age >= 19
  h <- .who(V$hz, V$age, s, "height", LMS); bm <- .who(V$bz, V$age, s, "bmi", LMS)
  for (sx in c("male", "female")) { w <- ad & s == sx; mu <- c(male = 176, female = 163)[[sx]]; L <- ADULT_L$height[[sx]]; t0 <- (mu^L - 1) / L
    h[w] <- .inv_bc(t0 + mu^(L - 1) * 6.8 * V$hz[w], L); L <- ADULT_L$bmi[[sx]]; t0 <- (26^L - 1) / L; bm[w] <- .inv_bc(t0 + 26^(L - 1) * 4.2 * V$bz[w], L) }
  hc <- ifelse(V$age < 18, .inv_lms(V$cz, V$age, s, "uk1990", "hc", LMS), ifelse(s == "male", 57, 55.3) + 1.8 * V$cz)
  out <- data.frame(person = id[V$i], fam_id = fam[V$i], sex_mf = ifelse(s == "male", "M", "F"), age_years = round(V$age, 3), height = round(h, 2), weight = round(bm * (h / 100)^2, 3),
                    head_circ = round(hc, 2), cnv = c(`No CNV` = "non-carrier", `16p11.2 del` = "16p DEL", `16p11.2 dup` = "16p DUP", `22q11.2 del` = "22q DEL", `22q11.2 dup` = "22q DUP")[g[V$i]],
                    asd_dx = c("no", "yes")[asd[V$i] + 1], stringsAsFactors = FALSE)
  first <- !duplicated(V$i); kid <- a0[V$i] >= 4 & a0[V$i] < 17
  out$srs <- ifelse(first & kid, srs[V$i], NA); out$scq <- if (scq) ifelse(first & kid, scqv[V$i], NA) else NA
  out$age_at_questionnaire <- ifelse(first & kid, round(V$age + 0.2, 3), NA); out$clinic <- label
  out
}

make_all <- function(out_dir, lms_file) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); LMS <- read.csv(lms_file)
  specs <- list(
    A = list(n = c(`No CNV` = 3200, `16p11.2 del` = 600, `16p11.2 dup` = 260, `22q11.2 del` = 340, `22q11.2 dup` = 160), visits = 2, seed = 11, scq = TRUE, small_cell = FALSE,
             ages = function(n) ifelse(runif(n) < 0.55, runif(n, 2, 18.9), runif(n, 19, 72))),
    B = list(n = c(`No CNV` = 1400, `16p11.2 del` = 400, `16p11.2 dup` = 150, `22q11.2 del` = 300, `22q11.2 dup` = 48), visits = 5, seed = 22, scq = FALSE, small_cell = TRUE,
             ages = function(n) ifelse(runif(n) < 0.97, runif(n, 0.3, 15.5), runif(n, 19, 45))))
  files <- list()
  for (s in names(specs)) { sp <- specs[[s]]; d <- make_site(s, sp$n, sp$ages, sp$visits, sp$seed, LMS, sp$scq, sp$small_cell)
    f <- file.path(out_dir, sprintf("visits_%s.tsv", s)); write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    cfg <- list(site_label = paste0("SYNTHETIC_", s), input_file = f, out_dir = file.path(out_dir, paste0("export_", s)),
      columns = list(id = "person", family = "fam_id", sex = "sex_mf", age = "age_years", height_cm = "height", weight_kg = "weight", hc_cm = "head_circ", genotype = "cnv",
                     asd = "asd_dx", srs_total = "srs", scq_total = if (sp$scq) "scq" else NULL, age_questionnaire = "age_at_questionnaire", cohort = "clinic"),
      sex_values = list(male = "M", female = "F"),
      genotype_map = list(`16p11.2 del` = "16p DEL", `16p11.2 dup` = "16p DUP", `22q11.2 del` = "22q DEL", `22q11.2 dup` = "22q DUP", control = "non-carrier"),
      asd_values = list(yes = "yes", no = "no"), min_cell = 5)
    cf <- file.path(out_dir, sprintf("config_%s.json", s)); writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, null = "null"), cf)
    files[[s]] <- list(data = f, config = cf, out = cfg$out_dir, cfg = cfg) }
  list(files = files, truth = TRUTH)
}

if (sys.nframe() == 0 && !exists("FEDERATED_SOURCE_ONLY")) {
  a <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)][1]), fixed = TRUE)
  out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- file.path(tempdir(), "federated_synthetic")
  r <- make_all(out, file.path(dirname(dirname(normalizePath(f))), "growth_refs_lms.csv")); cat("synthetic sites written to", normalizePath(out), "\n")
}
