# Does BMI act on autistic traits (SRS, SCQ) in CNV carriers, or does the BMI slope reflect genetics shared between body mass and
# growth? (JS / Claude 2026-11-14). Carriers of the four CNVs pooled (G2MH, the only carriers with polygenic scores), within-carrier
# models with CNV fixed effects, CNV-specific quadratic age, sex. Scores (residualised on 10 PCs, z within cohort): adiposity =
# childhood obesity (cobes), childhood BMI (cbmi), adult BMI (bmi); growth = adult height (ht), infant head circumference (hci,
# proxy for brain growth: no intracranial-volume score exists in the lab libraries), age at menarche (aam, Guan 2025 direct effect).
#  T1 first stage: each score -> measured BMI (partial F)
#  T2 reduced form: each score -> SRS / SCQ, alone and with the score's measured-BMI path removed (y ~ score + BMI)
#  T3 observational BMI slope under adjustment: none; growth scores; adiposity scores; all scores; + measured height and head size
#  T4 two-stage least squares, BMI instrumented by the three adiposity scores, growth scores as covariates; robust SE; first-stage F;
#     Sargan over-identification test; per-instrument Wald ratios; expected reduced form under full mediation (first stage x slope)
# Setup (scores, phenotype table U) reused from mr_reduced_form.R up to its analysis block. Output: mr_bmi_srs_scq_causal.csv
# (aggregate statistics only).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
src <- readLines(file.path(HERE, "mr_reduced_form.R")); stop_at <- grep("res <- list()", src, fixed = TRUE)[1] - 1
eval(parse(text = src[1:stop_at]))
suppressPackageStartupMessages(library(sandwich))
# extra G2MH scores: menarche (Jon/Claude data folder) and infant head circumference (Anjali scores_wgs)
MP <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/menarche_pgs")
gm <- fread(file.path(MP, "g2mh_wgs.age_menarche_guan2025.sscore")); setnames(gm, "#IID", "IID"); gm <- gm[, .(GUID = sub("^[0-9]+-", "", IID), aam = SCORE1_AVG)]
hc <- rdw("headcircumference_infant"); setnames(hc, "score", "hci")
ex <- merge(merge(gm, hc, by = "GUID"), pcw, by = "GUID"); ex[, `:=`(pgs_aam = rz(aam, ex), pgs_hci = rz(hci, ex))]
# measured growth of the same carriers: mean height z (children and adults, each scored against its own reference) and child HC z
gro <- ALL[carrier == TRUE & !is.na(key) & pheno %in% c("Height", "Head circumference"), .(z = mean(y, na.rm = TRUE)), by = .(key, pheno)]
gro <- dcast(gro, key ~ pheno, value.var = "z"); setnames(gro, c("Height", "Head circumference"), c("ht_z", "hc_z"))
W <- unique(U[pheno %in% c("SRS", "SCQ") & carrier == TRUE], by = c("key", "pheno"))
W <- merge(W, ex[, .(key = GUID, pgs_aam, pgs_hci)], by = "key")
W <- merge(W, gro, by = "key", all.x = TRUE)
W[, `:=`(locus = factor(locus), male = as.numeric(sex == "male"))]
cat("carriers with all scores:\n"); print(W[, .N, by = .(pheno, locus)][order(pheno, locus)])
cat("score correlations (one row per carrier):\n"); print(round(cor(unique(W[, .(key, pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht, pgs_hci, pgs_aam)])[, -1]), 2))

BASE <- "locus + locus:poly(age, 2) + male"
ADI <- c("pgs_cobes", "pgs_cbmi", "pgs_bmi"); GRO <- c("pgs_ht", "pgs_hci", "pgs_aam"); ALLS <- c(ADI, GRO)
rob <- function(m, v) { V <- vcovHC(m, type = "HC1"); b <- unname(coef(m)[v]); s <- sqrt(V[v, v]); c(b = unname(b), se = s, p = 2 * pnorm(-abs(b / s))) }
fm <- function(y, terms) as.formula(paste(y, "~", paste(c(BASE, terms), collapse = " + ")))
out <- list(); add <- function(...) out[[length(out) + 1]] <<- data.table(...)
for (ph in c("SRS", "SCQ")) {
  d <- W[pheno == ph]
  # T1 first stage (the same carriers)
  for (v in ALLS) { m <- lm(fm("BMI", v), d); r <- rob(m, v); add(outcome = ph, table = "T1 first stage (score -> BMI)", term = v, adjust = "", n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]], F = (r[["b"]] / r[["se"]])^2) }
  m <- lm(fm("BMI", ADI), d); m0 <- lm(fm("BMI", NULL), d); a <- anova(m0, m)
  add(outcome = ph, table = "T1 first stage (score -> BMI)", term = "three adiposity scores jointly", adjust = "", n = nrow(d), b = NA, se = NA, p = a$`Pr(>F)`[2], F = a$F[2])
  # T2 reduced form, alone and given measured BMI
  for (v in ALLS) for (adj in c("", "BMI")) { m <- lm(fm("y", c(v, if (adj == "BMI") "BMI")), d); r <- rob(m, v)
    add(outcome = ph, table = "T2 reduced form (score -> trait)", term = v, adjust = if (adj == "") "none" else "measured BMI", n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]], F = NA) }
  # T3 observational BMI slope under adjustment
  sets <- list(none = NULL, `growth scores (ht, hci, aam)` = GRO, `adiposity scores` = ADI, `all six scores` = ALLS)
  for (s in names(sets)) { m <- lm(fm("y", c("BMI", sets[[s]])), d); r <- rob(m, "BMI")
    add(outcome = ph, table = "T3 BMI slope", term = "BMI", adjust = s, n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]], F = NA) }
  for (g in c("ht_z", "hc_z")) { dd <- d[!is.na(get(g))]; for (tt in list("BMI", c("BMI", g))) { m <- lm(fm("y", tt), dd); r <- rob(m, "BMI")
    add(outcome = ph, table = "T3 BMI slope", term = "BMI", adjust = paste0(if (length(tt) == 1) "none" else paste("measured", g), " (", g, " subset)"), n = nrow(dd), b = r[["b"]], se = r[["se"]], p = r[["p"]], F = NA) } }
  dg <- d[!is.na(ht_z) & !is.na(hc_z)]
  for (s in c("none (growth subset)", "measured height + HC", "measured height + HC + all six scores")) {
    tt <- switch(s, `none (growth subset)` = "BMI", `measured height + HC` = c("BMI", "ht_z", "hc_z"), c("BMI", "ht_z", "hc_z", ALLS))
    m <- lm(fm("y", tt), dg); r <- rob(m, "BMI"); add(outcome = ph, table = "T3 BMI slope", term = "BMI", adjust = s, n = nrow(dg), b = r[["b"]], se = r[["se"]], p = r[["p"]], F = NA) }
  # T4 2SLS: BMI instrumented by the adiposity scores, growth scores as covariates
  Xc <- model.matrix(fm("y", GRO), d); Z <- cbind(Xc, as.matrix(d[, ..ADI])); X <- cbind(Xc, BMI = d$BMI); y <- d$y
  Xh <- Z %*% solve(crossprod(Z), crossprod(Z, X)); A <- solve(crossprod(Xh, X)); b <- A %*% crossprod(Xh, y); e <- as.vector(y - X %*% b)
  V <- A %*% crossprod(Xh * e) %*% t(A) * nrow(X) / (nrow(X) - ncol(X)); k <- ncol(X)
  sarg <- nrow(d) * summary(lm(e ~ Z - 1))$r.squared
  fs <- anova(lm(fm("BMI", GRO), d), lm(fm("BMI", c(GRO, ADI)), d))
  add(outcome = ph, table = "T4 2SLS (BMI <- adiposity scores; growth scores as covariates)", term = "BMI (IV)", adjust = sprintf("Sargan chi2 %.2f df 2 P %.3g", sarg, pchisq(sarg, 2, lower.tail = FALSE)),
      n = nrow(d), b = b[k], se = sqrt(V[k, k]), p = 2 * pnorm(-abs(b[k] / sqrt(V[k, k]))), F = fs$F[2])
  obs <- rob(lm(fm("y", c("BMI", GRO)), d), "BMI")[["b"]]
  for (v in ADI) { f1 <- rob(lm(fm("BMI", c(v, GRO)), d), v); rf <- rob(lm(fm("y", c(v, GRO)), d), v)
    add(outcome = ph, table = "T4 Wald ratio and mediation check", term = v, adjust = sprintf("first stage %.3f; reduced form %.3f; expected under full mediation (first stage x observed slope %.3f) %.3f", f1[["b"]], rf[["b"]], obs, f1[["b"]] * obs),
        n = nrow(d), b = rf[["b"]] / f1[["b"]], se = abs(rf[["se"]] / f1[["b"]]), p = rf[["p"]], F = (f1[["b"]] / f1[["se"]])^2) }
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "mr_bmi_srs_scq_causal.csv"))
options(width = 250); print(R[, .(outcome, table, term, adjust, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2), F = round(F, 1))], nrows = 200)
