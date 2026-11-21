# 86_v2_cnv_bmi_pgs_stats.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 3: statistics)
# What it does (section 3, exploratory): CNV x BMI x autism polygenic score on autistic traits.
#   1. Carries forward the 62 / 64 / 66 estimates (re-run by v2_base(); statistics unchanged) and adds, for the three-way term
#      (SRS primary, SCQ for comparison): leave-one-cohort-out, leave-one-CNV-out, and a permutation P.
#   2. BMI slope on SRS / SCQ by fifth of the autism score, carriers and SSC relatives (for the main figure).
#   3. NEW secondary outcome: autism diagnosis in carriers, logistic diagnosis ~ BMI x autism score + covariates, per CNV, pooled,
#      by cohort.
#   4. Mendelian randomisation of BMI -> SRS / SCQ, re-run and extended: G2MH carriers as in 49 (ten scores); NEW all genotyped
#      carriers (G2MH + Simons VIP + SSC; the three adiposity scores, height and autism scores exist in all three cohorts), SSC
#      unaffected relatives and SSC autistic probands without a CNV as separate strata; Anderson-Rubin weak-instrument-robust
#      intervals, leave-one-instrument-out, minimal detectable effect at 80% power.
# Inputs : v2_base() objects (XG = genotyped SRS / SCQ sample with the April 2026 autism score, XP, XA, S, idmap); 86a cache (G2MH
#   ten scores, W49); Simons VIP scores (Anjali searchlight/pgs/combined_scores.tsv, same SBayesRC weights as G2MH; residualised on
#   10 global PCs, z within the 393 genomes); SSC scores (DATA/ssc_pgs_scores.tsv).
# Outputs: autism/threeway_pgs_robustness.csv, autism/pgs_bmi_slope_by_score_fifth.csv, autism/asd_dx_bmi_x_pgs_logistic.csv,
#   autism/mr_bmi_autistic_traits_v2.csv; participant-level cache DATA/v2_autism_cache/pgs_stats_2026-11-20.RData for 86b / 86c.
# Models : three-way (64 three()): y ~ carrier * (poly(age, 2) + male) + cohort:CNV + carrier * BMI * score, HC1 robust SE;
#   within carriers (62): y ~ cohort:CNV + poly(age, 2) + male + BMI * score (+ cohort x BMI + cohort x score).
# null model: three-way and within-carrier terms: the autism score is permuted within cohort x CNV among carriers and, separately,
#   among SSC relatives (10,000 permutations; statistic = robust t of the interaction; two-sided); breaks any link between the score
#   and BMI x trait while keeping the BMI-trait relation, cohort and CNV composition. All other tests: Wald (HC1), Cochran Q.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
suppressPackageStartupMessages(library(splines))
v2_base(); DBX <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; DATA0 <- file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"); AN <- file.path(DBX, "G2MH/NetworkPaper/AnalysisByAnjali/data")
load(file.path(DATA0, "v2_autism_cache", "pgs49_2026-11-20.RData"))
NPERM <- as.integer(Sys.getenv("V2_NPERM", "10000")); options(width = 230)
# ================================================================ 1. three-way term: robustness
within <- function(car, full = TRUE) { if (nrow(car) < 15) return(NULL); cgt <- if (uniqueN(car$cg) > 1) "cg +" else ""; ct <- if (full && uniqueN(car$cohort) > 1) "+ cohort:BMI + cohort:pgs" else ""
  c(n = nrow(car), rob(lm(as.formula(paste("y ~", cgt, "poly(age, 2) + male + BMI * pgs", ct)), car), "BMI:pgs")) }
rows <- list(); addr <- function(ph, an, sub, r) if (!is.null(r)) rows[[length(rows) + 1]] <<- data.table(trait = ph, analysis = an, subset = sub, n_carriers = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
A3 <- "three-way CNV x BMI x autism score vs SSC relatives"; AW <- "within carriers BMI x autism score (cohort x BMI and cohort x score terms)"; AW0 <- "within carriers BMI x autism score (no cohort interaction terms)"
for (ph in c("SRS", "SCQ")) { car <- XG[pheno == ph & is_cnv == 1]; ref <- RELf(XG[pheno == ph])
  sets <- c(list(`all carriers` = car), setNames(lapply(LOC4, function(L) car[locus == L]), paste("CNV:", LOC4)), setNames(lapply(c("G2MH", "SVIP", "SSC"), function(k) car[cohort == k]), paste("cohort:", c("G2MH", "Simons VIP", "SSC carriers"))),
            setNames(lapply(c("G2MH", "SVIP", "SSC"), function(k) car[cohort != k]), paste("leave out cohort:", c("G2MH", "Simons VIP", "SSC carriers"))), setNames(lapply(LOC4, function(L) car[locus != L]), paste("leave out CNV:", LOC4)))
  for (s in names(sets)) { addr(ph, A3, s, three(sets[[s]], ref)); addr(ph, AW, s, within(sets[[s]])); addr(ph, AW0, s, within(sets[[s]], FALSE)) }
  r <- rob(lm(y ~ poly(age, 2) + male + BMI * pgs, ref), "BMI:pgs"); rows[[length(rows) + 1]] <- data.table(trait = ph, analysis = "BMI x autism score in SSC relatives", subset = "SSC unaffected relatives", n_carriers = nrow(ref), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
RB <- rbindlist(rows)
# permutation P (score permuted within cohort x CNV in carriers and within the relatives)
permfit <- function(z, f0, pcols, B = NPERM, seed = 20261120) {   # f0: model without score terms; pcols: columns multiplied by the score, the last one is tested
  m0 <- lm(f0, z); X0 <- model.matrix(m0)[, !is.na(coef(m0)), drop = FALSE]; M <- sapply(pcols, function(e) eval(parse(text = e), z)); y <- z$y
  stat <- function(p) { Xp <- cbind(X0, M * p); j <- ncol(Xp); XtXi <- chol2inv(chol(crossprod(Xp))); b <- XtXi %*% crossprod(Xp, y); e <- as.vector(y - Xp %*% b); a <- as.vector(Xp %*% XtXi[, j]); b[j] / sqrt(sum(a^2 * e^2)) }
  t0 <- stat(z$pgs); idx <- split(seq_len(nrow(z)), z$cg); set.seed(seed)
  run <- function(B) replicate(B, { p <- z$pgs; for (ix in idx) p[ix] <- p[ix][sample.int(length(ix))]; stat(p) })
  tp <- run(B); k <- sum(abs(tp) >= abs(t0))
  while (k < 5 && length(tp) < 2e5) { tp <- c(tp, run(B * 4)); k <- sum(abs(tp) >= abs(t0)) }   # P at the floor: add replicates
  c(t_obs = t0, perm_p = (1 + k) / (length(tp) + 1), n_perm = length(tp), n_as_extreme = k) }
PM <- rbindlist(lapply(c("SRS", "SCQ"), function(ph) { car <- XG[pheno == ph & is_cnv == 1]; ref <- RELf(XG[pheno == ph]); cols <- c("y", "BMI", "pgs", "age", "male", "cg", "is_cnv")
  z <- rbind(car[, ..cols], ref[, ..cols]); p3 <- permfit(z, y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI, c("rep(1, length(y))", "is_cnv", "BMI", "is_cnv * BMI"))
  pw <- permfit(car[, ..cols], y ~ cg + poly(age, 2) + male + BMI, c("rep(1, length(y))", "BMI"))
  rbind(data.table(trait = ph, analysis = A3, subset = "all carriers", t(p3)), data.table(trait = ph, analysis = AW0, subset = "all carriers", t(pw))) }))
RB <- merge(RB, PM, by = c("trait", "analysis", "subset"), all.x = TRUE, sort = FALSE)
RB[, robust_t := b / se]
for (ph in c("SRS", "SCQ")) for (an in c(A3, AW0)) { r <- RB[trait == ph & analysis == an & grepl("^CNV:", subset)]; q <- cochranQ(r$b, r$se)
  RB <- rbind(RB, data.table(trait = ph, analysis = an, subset = "heterogeneity across the four CNVs (Cochran Q)", p = q[["p"]]), fill = TRUE) }
old <- fread(file.path(HERE, "fig_supp_interactions_main_table.csv")); cat("66 re-run identical to saved table:", isTRUE(all.equal(old$b, R66$b)) && all(old$n == R66$n), "\n")
RB[, note := "Box-Cox scale; carriers G2MH + Simons VIP + SSC with the April 2026 autism score; reference SSC unaffected relatives; HC1 robust SE; exploratory"]
fwrite(RB, file.path(V2OUT, "threeway_pgs_robustness.csv"))
print(RB[, .(trait, analysis = substr(analysis, 1, 30), subset, n_carriers, b = round(b, 3), se = round(se, 3), p = signif(p, 3), perm_p = signif(perm_p, 3), n_perm)], nrows = 200)
# ================================================================ 2. BMI slope by fifth of the autism score
fifth <- function(d, f, grp, ph) { d <- copy(d)[, bin := cut(pgs, quantile(pgs, 0:5 / 5), include.lowest = TRUE, labels = FALSE)]
  d[, { r <- rob(lm(f, .SD), "BMI"); .(score_mean = mean(pgs), n = .N, slope = r[["b"]], se = r[["se"]]) }, by = bin][order(bin)][, `:=`(group = grp, trait = ph)] }
line <- function(d, f, grp, ph) { m <- lm(f, d); V <- vcovHC(m, type = "HC1"); g <- seq(-2.2, 2.2, 0.1); b <- coef(m); s <- sqrt(V["BMI", "BMI"] + g^2 * V["BMI:pgs", "BMI:pgs"] + 2 * g * V["BMI", "BMI:pgs"])
  data.table(group = grp, trait = ph, score = g, slope = b[["BMI"]] + g * b[["BMI:pgs"]], se = s) }
FQ <- rbindlist(lapply(c("SRS", "SCQ"), function(ph) rbind(fifth(XG[pheno == ph & is_cnv == 1], y ~ cg + poly(age, 2) + male + BMI, "CNV carriers", ph), fifth(RELf(XG[pheno == ph]), y ~ poly(age, 2) + male + BMI, "SSC unaffected relatives", ph))))
FL <- rbindlist(lapply(c("SRS", "SCQ"), function(ph) rbind(line(XG[pheno == ph & is_cnv == 1], y ~ cg + poly(age, 2) + male + BMI * pgs, "CNV carriers", ph), line(RELf(XG[pheno == ph]), y ~ poly(age, 2) + male + BMI * pgs, "SSC unaffected relatives", ph))))
fwrite(rbind(FQ[, .(trait, group, kind = "fifth of the autism score", bin, score = score_mean, n, slope, se)], FL[, .(trait, group, kind = "continuous model", bin = NA, score, n = NA, slope, se)]), file.path(V2OUT, "pgs_bmi_slope_by_score_fifth.csv"))
print(FQ[, .(trait, group, bin, score = round(score_mean, 2), n, slope = round(slope, 3), se = round(se, 3))])
# ================================================================ 3. autism diagnosis ~ BMI x autism score (carriers)
P1 <- XG[is_cnv == 1 & !is.na(asd2)][order(id, pheno != "SRS")][, .SD[1], by = id][, .(id, cohort, locus, male, age, BMI, pgs, asd = asd2)]
g <- merge(XA[pheno == "BMI" & locus %in% LOC4], idmap, by = c("id", "cohort"))[, key := fifelse(cohort %in% c("G2MH", "Cardiff"), GUID, sampleID)]
P2 <- merge(g, S, by.x = "key", by.y = "id")[(cohort %in% c("G2MH", "Cardiff") & src == "G2MH") | (cohort == "SVIP" & src == "SVIP") | (cohort == "SSC" & src == "SSC")][, .(id, cohort, locus, male = as.numeric(sex == "male"), age, BMI, pgs, asd)]
dlg <- function(d) tryCatch({ d <- copy(d); keep <- d[, .(v = sum(asd) > 0 & sum(asd) < .N), by = cohort][v == TRUE, cohort]; d <- d[cohort %in% keep]; if (sum(d$asd) < 10 || sum(1 - d$asd) < 10) stop()
  tm <- c("BMI * pgs", "ns(age, 3)", "male", if (uniqueN(d$cohort) > 1) "cohort", if (uniqueN(d$locus) > 1) "locus"); m <- glm(as.formula(paste("asd ~", paste(tm, collapse = " + "))), data = d, family = binomial); co <- coef(summary(m))
  data.table(n = nrow(d), n_asd = sum(d$asd), cohorts_used = paste(sort(unique(d$cohort)), collapse = "+"), term = c("BMI", "autism score", "BMI x autism score"), logOR = co[c("BMI", "pgs", "BMI:pgs"), 1], se = co[c("BMI", "pgs", "BMI:pgs"), 2], p = co[c("BMI", "pgs", "BMI:pgs"), 4]) },
  error = function(e) data.table(n = nrow(d), n_asd = sum(d$asd), term = "BMI x autism score", logOR = NA_real_, se = NA_real_, p = NA_real_))
DX <- rbindlist(lapply(list(list("carriers of the SRS / SCQ sample with the score (BMI at the questionnaire visit)", P1), list("growth-table carriers with the score (fill_gaps visit; G2MH CHOP-site carriers excluded by the CHOP overlap rule)", P2)), function(s) { d <- s[[2]]
  o <- list(cbind(sample = s[[1]], group = "All four CNVs", dlg(d))); for (L in LOC4) o[[length(o) + 1]] <- cbind(sample = s[[1]], group = L, dlg(d[locus == L]))
  for (k in sort(unique(d$cohort))) o[[length(o) + 1]] <- cbind(sample = s[[1]], group = paste("All four CNVs, cohort:", k), dlg(d[cohort == k]))
  for (k in sort(unique(d$cohort))) o[[length(o) + 1]] <- cbind(sample = s[[1]], group = paste("All four CNVs, leave out", k), dlg(d[cohort != k]))
  r <- rbindlist(o, fill = TRUE); q <- r[term == "BMI x autism score" & group %in% LOC4, cochranQ(logOR, se)]
  rbind(r, data.table(sample = s[[1]], group = "heterogeneity across the four CNVs (Cochran Q)", term = "BMI x autism score", p = q[["p"]]), fill = TRUE) }), fill = TRUE)
DX[, `:=`(OR = exp(logOR), OR_lo = exp(logOR - 1.96 * se), OR_hi = exp(logOR + 1.96 * se))][, note := "logistic: diagnosis ~ BMI x autism score + ns(age,3) + sex + cohort + CNV; logOR per SD (x SD); carriers only; exploratory"]
fwrite(DX, file.path(V2OUT, "asd_dx_bmi_x_pgs_logistic.csv"))
print(DX[term == "BMI x autism score" | is.na(term), .(sample = substr(sample, 1, 22), group, n, n_asd, OR = round(OR, 3), lo = round(OR_lo, 2), hi = round(OR_hi, 2), p = signif(p, 3))], nrows = 100)
# ================================================================ 4. Mendelian randomisation
rzp <- function(v, P) as.numeric(scale(residuals(lm(as.formula(paste("v ~", paste(names(P), collapse = " + "))), data = cbind(P, v = v)))))
cs <- fread(file.path(AN, "searchlight/pgs/combined_scores.tsv"), select = c("IID", "childhood_obesity", "bmi_childhood", "bmi", "height"))
pcs <- fread(file.path(AN, "searchlight/pgs/analysis_dataset.tsv"), select = c("IID", paste0("GLOBAL_PC", 1:10)))
meta <- unique(fread(file.path(DBX, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/Searchlight.WGS2.2023_05.sample_metadata.tsv"))[, .(searchlight_id, sfari_iid)])
sv <- merge(merge(cs, pcs, by = "IID"), meta, by.x = "IID", by.y = "searchlight_id"); PCm <- sv[, paste0("GLOBAL_PC", 1:10), with = FALSE]
sv[, `:=`(pgs_cobes = rzp(childhood_obesity, PCm), pgs_cbmi = rzp(bmi_childhood, PCm), pgs_bmi = rzp(bmi, PCm), pgs_ht = rzp(height, PCm))][, id := sub("-", ".", sfari_iid, fixed = TRUE)]
ss <- fread(file.path(DATA0, "ssc_pgs_scores.tsv"))
ADS <- rbind(g2all[, .(id = GUID, src = "G2MH", pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht)], sv[, .(id, src = "SVIP", pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht)], ss[, .(id = individual, src = "SSC", pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht)])
XM <- merge(XG, ADS, by = c("id", "src")); XM <- XM[complete.cases(XM[, .(pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht, pgs)])]
PR <- merge(XP[grp == "SSC ASD probands" & pheno %in% c("SRS", "SCQ"), .(id, pheno, sex, age, BMI, y)], ss[, .(id = individual, pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht, pgs = pgs_asd)], by = "id")[, male := as.numeric(sex == "male")]
INS <- c("pgs_cobes", "pgs_cbmi", "pgs_bmi"); INSL <- c(pgs_cobes = "childhood obesity score", pgs_cbmi = "childhood BMI score", pgs_bmi = "adult BMI score")
mrfit <- function(d, base, cov, ins = INS) {
  d <- d[complete.cases(d[, c("y", "BMI", cov, ins), with = FALSE])]; Xc <- model.matrix(as.formula(paste("~", paste(c(base, cov), collapse = " + "))), d); q <- qr(Xc); kx <- q$rank; n <- nrow(d)
  y <- qr.resid(q, d$y); x <- qr.resid(q, d$BMI); Z <- qr.resid(q, as.matrix(d[, ..ins])); hc <- sqrt(n / (n - kx - 1))
  ols_b <- sum(x * y) / sum(x^2); ols_se <- sqrt(sum(x^2 * (y - ols_b * x)^2)) / sum(x^2) * hc
  iv <- function(Zs) { Zs <- as.matrix(Zs); L <- ncol(Zs); qz <- qr(Zs); xh <- qr.fitted(qz, x); b <- sum(xh * y) / sum(xh * x); e <- y - b * x; se <- sqrt(sum(xh^2 * e^2)) / abs(sum(xh * x)) * hc
    r2 <- sum(xh^2) / sum(x^2); Fs <- (r2 / L) / ((1 - r2) / (n - kx - L)); sg <- if (L > 1) pchisq(n * sum(qr.fitted(qz, e)^2) / sum(e^2), L - 1, lower.tail = FALSE) else NA_real_
    g <- seq(-10, 10, 0.005); Py <- qr.fitted(qz, y); ex <- sum(Py^2) - 2 * g * sum(Py * xh) + g^2 * sum(xh^2); tot <- sum(y^2) - 2 * g * sum(x * y) + g^2 * sum(x^2)   # Anderson-Rubin: are the instruments unrelated to y - b0 * BMI?
    pa <- pf((ex / L) / ((tot - ex) / (n - kx - L)), L, n - kx - L, lower.tail = FALSE); acc <- g[pa > 0.05]
    arn <- if (!length(acc)) "empty (instruments reject every value: they disagree)" else if (min(acc) <= -10 + 1e-9 && max(acc) >= 10 - 1e-9 && length(acc) == length(g)) "unbounded (no information)" else if (min(acc) <= -10 + 1e-9 || max(acc) >= 10 - 1e-9) "unbounded on one or both sides" else if (any(diff(acc) > 0.0051)) "disjoint" else "bounded"
    list(b = b, se = se, F = Fs, sargan_p = sg, ar_lo = if (length(acc)) min(acc) else NA_real_, ar_hi = if (length(acc)) max(acc) else NA_real_, ar_note = arn) }
  a <- iv(Z); o <- list(data.table(estimator = "Observational (least squares)", b = ols_b, se = ols_se, p = 2 * pnorm(-abs(ols_b / ols_se))),
    data.table(estimator = "MR, three adiposity scores (2SLS)", b = a$b, se = a$se, p = 2 * pnorm(-abs(a$b / a$se)), first_stage_F = a$F, sargan_p = a$sargan_p, ar_lo = a$ar_lo, ar_hi = a$ar_hi, ar_note = a$ar_note,
               mde80 = 2.802 * a$se, power_at_observational_slope = pnorm(abs(ols_b) / a$se - 1.96)))
  for (j in seq_along(ins)) { z <- Z[, j]; fs <- sum(z * x) / sum(z^2); fse <- sqrt(sum(z^2 * (x - fs * z)^2)) / sum(z^2) * hc; rf <- sum(z * y) / sum(z^2); rse <- sqrt(sum(z^2 * (y - rf * z)^2)) / sum(z^2) * hc; w <- iv(z)
    o[[length(o) + 1]] <- data.table(estimator = paste0("MR, ", INSL[[ins[j]]], " alone (Wald ratio)"), b = rf / fs, se = abs(rse / fs), p = 2 * pnorm(-abs(rf / rse)), first_stage_F = (fs / fse)^2, ar_lo = w$ar_lo, ar_hi = w$ar_hi, ar_note = w$ar_note) }
  for (j in seq_along(ins)) { w <- iv(Z[, -j, drop = FALSE]); o[[length(o) + 1]] <- data.table(estimator = paste0("MR, leaving out the ", INSL[[ins[j]]], " (2SLS)"), b = w$b, se = w$se, p = 2 * pnorm(-abs(w$b / w$se)), first_stage_F = w$F, sargan_p = w$sargan_p, ar_lo = w$ar_lo, ar_hi = w$ar_hi, ar_note = w$ar_note) }
  cbind(n = n, rbindlist(o, fill = TRUE)) }
COVC <- c("pgs_ht", "pgs"); W49[, cg := locus]
MR <- rbindlist(lapply(c("SRS", "SCQ"), function(ph) { car <- XM[pheno == ph & is_cnv == 1]; rel <- RELf(XM[pheno == ph]); pr <- PR[pheno == ph]; bc <- "cg + locus:poly(age, 2) + male"; o <- list()
  ad <- function(st, spec, r) o[[length(o) + 1]] <<- cbind(trait = ph, stratum = st, covariates = spec, r)
  ad("CNV carriers: G2MH + Simons VIP + SSC", "cohort x CNV, CNV-specific quadratic age, sex, height and autism scores", mrfit(car, bc, COVC))
  ad("CNV carriers: G2MH", "CNV, CNV-specific quadratic age, sex, height and autism scores", mrfit(car[cohort == "G2MH"], bc, COVC))
  ad("CNV carriers: Simons VIP", "CNV, CNV-specific quadratic age, sex, height and autism scores", mrfit(car[cohort == "SVIP"], bc, COVC))
  ad("CNV carriers: G2MH (model of 49)", "CNV, CNV-specific quadratic age, sex, three growth and four neurodevelopmental scores", mrfit(W49[pheno == ph], BASE, c(GRO, ND)))
  ad("SSC unaffected relatives", "quadratic age, sex, height and autism scores", mrfit(rel, "poly(age, 2) + male", COVC))
  ad("SSC autistic probands, no CNV", "quadratic age, sex, height and autism scores", mrfit(pr, "poly(age, 2) + male", COVC)); rbindlist(o, fill = TRUE) }), fill = TRUE)
MR[, `:=`(ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se)][, note := "SD of trait (Box-Cox) per SD of BMI; instruments = polygenic scores for childhood obesity, childhood BMI, adult BMI; HC1 robust SE; Wald-ratio SE = reduced-form SE / first stage; AR = Anderson-Rubin 95% set on a grid of -10..10; mde80 = smallest causal effect detectable with 80% power at two-sided 0.05 (2.80 x SE)"]
c49 <- C49[est == "2SLS, three adiposity scores"]; cat("2SLS of 49 reproduced by mrfit (SRS, SCQ):", round(MR[stratum == "CNV carriers: G2MH (model of 49)" & grepl("^MR, three", estimator), b], 4), "vs", round(c49$b, 4), "; Sargan", signif(MR[stratum == "CNV carriers: G2MH (model of 49)" & grepl("^MR, three", estimator), sargan_p], 3), "vs", signif(c49$sargan_p, 3), "\n")
fwrite(MR, file.path(V2OUT, "mr_bmi_autistic_traits_v2.csv"))
print(MR[, .(trait, stratum = substr(stratum, 1, 30), est = substr(estimator, 1, 44), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2), F = round(first_stage_F, 1), sarg = signif(sargan_p, 2), ar_lo = round(ar_lo, 2), ar_hi = round(ar_hi, 2), ar = substr(ar_note, 1, 14), mde = round(mde80, 2), pw = round(power_at_observational_slope, 2))], nrows = 200)
save(RB, FQ, FL, DX, MR, XM, PR, P1, P2, file = file.path(DATA0, "v2_autism_cache", "pgs_stats_2026-11-20.RData"))
cat("86 done\n")
