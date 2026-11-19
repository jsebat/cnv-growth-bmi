# CNV x BMI x autism polygenic score on SRS / SCQ: formal three-way test and sensitivity analyses (JS / Claude 2026-11-14).
# Jon: is the within-carrier BMI x autism score interaction a three-way CNV x BMI x PGS interaction? Test it formally and probe the
# weak points. Box-Cox scale (00_boxcox_srs_scq.R) unless stated. G2MH carriers (the only carriers with scores) + SSC without CNV.
#  T1 three-way: carriers + SSC unaffected relatives (and, separately, + SSC autistic probands):
#     y ~ is_cnv * (poly(age,2) + male) + locus(carriers) + is_cnv * BMI * pgs_asd ; is_cnv:BMI:pgs_asd = three-way term
#  T2 within carriers, BMI x autism score (continuous, robust SE) under sensitivity specifications:
#     scale (Box-Cox / linear z / ordinal), covariate x BMI and covariate x score terms (Keller 2014: age, sex, CNV, 10 PCs),
#     BMI x the other nine scores, ASD diagnosis (and ASD x BMI), medication (any psychotropic, antipsychotic, weight-gain class,
#     and their x BMI), BMI winsorised at |z| 3, high-influence points removed (Cook's D > 4/n), children only, SRS child form only,
#     leave-one-CNV-out, CNV heterogeneity of the interaction, permutation P (score permuted within CNV, 10,000), selection checks
#     (score-BMI correlation in carriers, by ASD status), interaction in carriers with and without an ASD diagnosis.
#  T3 replication attempt: SSC carriers of the same CNVs (SSC genotypes, lo-SSC autism score), same model.
# Output: cnv_bmi_pgs_threeway.csv (aggregate statistics only).
Sys.setenv(SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s47 <- readLines(file.path(HERE, "47_bmi_srs_scq_pgs_causal.R")); eval(parse(text = s47[1:grep("^out <- list", s47)]))
suppressPackageStartupMessages(library(ordinal))
if (!"pgs_asd" %in% names(W)) W <- merge(W, g2[, .(key = GUID, pgs_asd)], by = "key")
s48 <- readLines(file.path(HERE, "48_bmi_srs_scq_asd_pgs.R")); eval(parse(text = s48[grep("^nd <- Reduce", s48):grep("^nd <- merge", s48)]))
W <- merge(W, nd[, .(key = GUID, pgs_adhd, pgs_edu, pgs_iq)], by = "key")
# raw score, SRS form, linear z, PCs, medication
W <- merge(W, D[cohort == "G2MH" & pheno %in% c("SRS", "SCQ"), .(key = id, pheno, raw, form2, y_lin, asd)], by = c("key", "pheno"), all.x = TRUE)
W <- merge(W, pcw[, c("GUID", paste0("PC", 1:10)), with = FALSE], by.x = "key", by.y = "GUID", all.x = TRUE)
s27 <- readLines(file.path(HERE, "27_medication_sensitivity.R")); eval(parse(text = s27[grep("^m <- fread", s27):grep("^M <- merge\\(M, unique", s27)]))
W <- merge(W, M[, .(key = GUID, psychotropic, antipsychotic, weight_gain, med_known = known)], by = "key", all.x = TRUE)
print(W[, .N, by = pheno]); cat("missing form/raw:", sum(is.na(W$raw)), " PCs:", sum(is.na(W$PC1)), " medication known:", sum(W$med_known %in% TRUE), "\n")
SC10 <- c("pgs_cobes", "pgs_cbmi", "pgs_bmi", "pgs_ht", "pgs_hci", "pgs_aam", "pgs_asd", "pgs_adhd", "pgs_edu", "pgs_iq")
BASEW <- "locus + locus:poly(age, 2) + male"
ix <- function(d, extra = NULL, yv = "y", base = BASEW) { d <- copy(d)[, yy := get(yv)]
  m <- lm(as.formula(paste("yy ~", base, "+ BMI * pgs_asd", if (length(extra)) paste("+", paste(extra, collapse = " + ")) else "")), d); r <- rob(m, "BMI:pgs_asd"); c(n = nrow(d), r) }
ord_ix <- function(d) { r <- t(sapply(split(d, d$form2), function(z) tryCatch({ z <- copy(z); br <- unique(quantile(z$raw, seq(0, 1, length.out = 41), type = 1))
    z[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)][, `:=`(a1 = as.numeric(scale(age)), a2 = as.numeric(scale(age))^2)]
    m <- clm(ry ~ locus + a1 + a2 + male + BMI * pgs_asd, data = z); co <- coef(summary(m)); c(co["BMI:pgs_asd", 1], co["BMI:pgs_asd", 2]) }, error = function(e) c(NA, NA))))
  r <- r[!is.na(r[, 1]), , drop = FALSE]; w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); se <- sqrt(1 / sum(w)); c(n = nrow(d), b = b, se = se, p = 2 * pnorm(-abs(b / se))) }
out <- list(); add <- function(ph, block, spec, r) out[[length(out) + 1]] <<- data.table(trait = ph, block = block, spec = spec, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
# ---------------------------------------------------------------- SSC (no CNV) with scores, for the three-way models and T3
pg <- fread(file.path(DATA, "ssc_pgs_scores.tsv"))
SS <- merge(X[cohort == "SSC" & pheno %in% c("SRS", "SCQ") & grp %in% c("Unaffected relatives", "SSC ASD probands")], pg, by.x = "id", by.y = "individual")[, male := as.numeric(sex == "male")]
for (ph in c("SRS", "SCQ")) {
  d <- W[pheno == ph]
  # ---------------------------------------------------------------- T1 three-way
  for (rf in c("Unaffected relatives", "SSC ASD probands")) {
    z <- rbind(d[, .(y, BMI, pgs_asd, age, male, locus = as.character(locus), is_cnv = 1)], SS[pheno == ph & grp == rf, .(y, BMI, pgs_asd, age, male, locus = "none", is_cnv = 0)])
    m <- lm(y ~ is_cnv * (poly(age, 2) + male) + locus + is_cnv * BMI * pgs_asd, z); r <- rob(m, "is_cnv:BMI:pgs_asd")
    add(ph, "T1 three-way CNV x BMI x autism score", paste("reference:", rf), c(n = nrow(z), r))
    r0 <- rob(m, "BMI:pgs_asd"); add(ph, "T1 three-way CNV x BMI x autism score", paste("BMI x score in reference:", rf), c(n = sum(z$is_cnv == 0), r0))
  }
  # ---------------------------------------------------------------- T2 sensitivity within carriers
  B <- "T2 within carriers: BMI x autism score"
  add(ph, B, "main (Box-Cox)", ix(d))
  add(ph, B, "linear z scale", ix(d, yv = "y_lin"))
  add(ph, B, "ordinal (scale-free)", ord_ix(d))
  pcs <- paste0("PC", 1:10)
  add(ph, B, "+ covariate x BMI and covariate x score (age, sex, CNV, 10 PCs)", ix(d[!is.na(PC1)], c(pcs, paste0(c("poly(age, 2)", "male", "locus", pcs), ":BMI"), paste0(c("poly(age, 2)", "male", "locus", pcs), ":pgs_asd"))))
  add(ph, B, "+ BMI x each of the other nine scores", ix(d, c(setdiff(SC10, "pgs_asd"), paste0(setdiff(SC10, "pgs_asd"), ":BMI"))))
  add(ph, B, "+ ASD diagnosis and ASD x BMI", ix(d[!is.na(asd)], c("factor(asd)", "factor(asd):BMI")))
  dm <- d[med_known %in% TRUE]
  add(ph, B, "medication known (same subset, unadjusted)", ix(dm))
  add(ph, B, "+ psychotropic, antipsychotic, weight-gain class and x BMI", ix(dm, c("psychotropic", "antipsychotic", "weight_gain", "psychotropic:BMI", "antipsychotic:BMI", "weight_gain:BMI")))
  add(ph, B, "no current psychotropic", ix(dm[psychotropic == FALSE]))
  add(ph, B, "BMI winsorised at |z| 3", ix(copy(d)[, BMI := pmin(pmax(BMI, -3), 3)]))
  m0 <- lm(as.formula(paste("y ~", BASEW, "+ BMI * pgs_asd")), d); keep <- cooks.distance(m0) <= 4 / nrow(d)
  add(ph, B, sprintf("Cook's D > 4/n removed (%d points)", sum(!keep)), ix(d[keep]))
  add(ph, B, "children (< 18 y)", ix(d[age < 18]))
  add(ph, B, "adults (>= 18 y)", ix(d[age >= 18]))
  if (ph == "SRS") add(ph, B, "SRS child form only", ix(d[form2 == "child"]))
  add(ph, B, "carriers with ASD diagnosis", ix(d[asd %in% 1L]))
  add(ph, B, "carriers without ASD diagnosis", ix(d[asd %in% 0L]))
  for (L in levels(d$locus)) add(ph, B, paste("leave out", L), ix(d[locus != L]))
  for (L in levels(d$locus)) { dl <- d[locus == L]; add(ph, B, paste("only", L), ix(dl, base = "poly(age, 2) + male")) }
  mh <- lm(y ~ locus * (poly(age, 2) + male) + locus * BMI * pgs_asd, d); mh0 <- lm(y ~ locus * (poly(age, 2) + male) + locus * BMI + locus * pgs_asd + BMI:pgs_asd, d)
  a <- anova(mh0, mh); add(ph, B, "heterogeneity of the interaction across CNVs (3 df F test)", c(n = nrow(d), b = NA, se = NA, p = a$`Pr(>F)`[2]))
  set.seed(20261114); t0 <- coef(summary(m0))["BMI:pgs_asd", 3]
  tp <- replicate(10000, { dd <- copy(d)[, pgs_asd := sample(pgs_asd), by = locus]; coef(summary(lm(as.formula(paste("y ~", BASEW, "+ BMI * pgs_asd")), dd)))["BMI:pgs_asd", 3] })
  add(ph, B, "permutation P (score permuted within CNV, 10,000)", c(n = nrow(d), b = coef(m0)[["BMI:pgs_asd"]], se = NA, p = (1 + sum(abs(tp) >= abs(t0))) / 10001))
  # selection checks
  for (s in list(list("all carriers", d), list("ASD diagnosis", d[asd %in% 1L]), list("no ASD diagnosis", d[asd %in% 0L]))) {
    r <- rob(lm(as.formula(paste("BMI ~", BASEW, "+ pgs_asd")), s[[2]]), "pgs_asd"); add(ph, "T2 selection check: autism score -> BMI", s[[1]], c(n = nrow(s[[2]]), r)) }
  r <- rob(glm(asd ~ locus + pgs_asd + BMI + pgs_asd:BMI, family = binomial, data = d[!is.na(asd)]), "pgs_asd:BMI"); add(ph, "T2 selection check: ASD diagnosis ~ score x BMI (log-odds)", "all carriers", c(n = sum(!is.na(d$asd)), r))
  cat(ph, "done\n")
}
# ---------------------------------------------------------------- T3 replication: SSC carriers of the four CNVs
sc <- gmc[cohort == "SSC" & CNV %in% names(lab4) & !is.na(BMI_raw) & !is.na(age), .(id = sampleID, locus = lab4[CNV], sex, age, BMI_raw)][!duplicated(id)]
sc <- merge(merge(sc, srs_ssc, by = "id", all.x = TRUE), scq_ssc, by = "id", all.x = TRUE)
sc[, BMI := bmiZ(BMI_raw, age, sex)][, male := as.numeric(sex == "male")]
sc <- merge(sc, pg[, .(id = individual, pgs_asd)], by = "id")
cat("SSC carriers with scores:", nrow(sc), " SRS:", sum(!is.na(sc$srs)), " SCQ:", sum(!is.na(sc$scq)), "\n"); print(sc[, .N, by = locus])
for (ph in c("SRS", "SCQ")) {
  if (ph == "SRS") { x <- sc[!is.na(srs) & !is.na(BMI)][, `:=`(raw = srs, f = form)] } else { x <- sc[!is.na(scq) & !is.na(BMI)][, `:=`(raw = scq, f = "scq")] }
  for (f in unique(x$f)) { r <- D[pheno == ph & form2 == f & unaffected == TRUE & locus == "No CNV", raw]; L <- LAM[[f]]
    bc <- function(v) if (abs(L) < 1e-8) log(v + 1) else ((v + 1)^L - 1) / L; x[f == f, y := (bc(raw) - mean(bc(r))) / sd(bc(r))] }
  if (nrow(x) >= 15) { m <- lm(y ~ locus + poly(age, 2) + male + BMI * pgs_asd, x); r <- rob(m, "BMI:pgs_asd"); add(ph, "T3 replication: SSC carriers", "BMI x autism score", c(n = nrow(x), r))
    r2 <- rob(m, "BMI"); add(ph, "T3 replication: SSC carriers", "BMI slope at mean score", c(n = nrow(x), r2)) }
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "cnv_bmi_pgs_threeway.csv"))
options(width = 230); print(R[, .(trait, block = substr(block, 1, 40), spec = substr(spec, 1, 62), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
