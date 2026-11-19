# Add-on to 47 (JS / Claude 2026-11-14, Jon: "should we add the autism PGS?"). Same carriers, same base model. Adds the autism score
# (April 2026 ASD GWAS, already in mr_reduced_form.R as pgs_asd) and the scores genetically correlated with BMI that could confound the
# BMI -> autistic-trait slope: ADHD, educational attainment, IQ (Anjali scores_wgs, residualised on 10 PCs, z within G2MH).
#  A  each score -> BMI and -> trait
#  B  BMI slope: none / + ASD / + ASD ADHD EA IQ / + those + the six growth/adiposity scores of 47
#  C  childhood-BMI score reduced form (the pleiotropic one in 47): alone / + ASD / + neurodevelopmental set
#  D  2SLS of 47 with the neurodevelopmental scores as extra covariates
# Output: mr_bmi_srs_scq_asd.csv (aggregate statistics only).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s47 <- readLines(file.path(HERE, "47_bmi_srs_scq_pgs_causal.R")); eval(parse(text = s47[1:grep("^out <- list", s47)]))
nd <- Reduce(function(a, b) merge(a, b, by = "GUID"), list(setnames(rdw("adhd"), "score", "adhd"), setnames(rdw("edu"), "score", "edu"), setnames(rdw("iq"), "score", "iq")))
nd <- merge(nd, pcw, by = "GUID"); nd[, `:=`(pgs_adhd = rz(adhd, nd), pgs_edu = rz(edu, nd), pgs_iq = rz(iq, nd))]
W <- merge(W, nd[, .(key = GUID, pgs_adhd, pgs_edu, pgs_iq)], by = "key"); if (!"pgs_asd" %in% names(W)) W <- merge(W, g2[, .(key = GUID, pgs_asd)], by = "key")
ND <- c("pgs_asd", "pgs_adhd", "pgs_edu", "pgs_iq")
cat("n:\n"); print(W[, .N, by = pheno])
cat("score correlations:\n"); print(round(cor(unique(W[, c("key", ND, ALLS), with = FALSE])[, -1]), 2)[ND, ])
for (ph in c("SRS", "SCQ")) {
  d <- W[pheno == ph]
  for (v in ND) for (y in c("BMI", "y")) { r <- rob(lm(fm(y, v), d), v)
    add(outcome = ph, table = "A score -> BMI / trait", term = v, adjust = if (y == "BMI") "-> BMI" else "-> trait", n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  sets <- list(none = NULL, `ASD score` = "pgs_asd", `ASD ADHD EA IQ` = ND, `ASD ADHD EA IQ + six scores of 47` = c(ND, ALLS))
  for (s in names(sets)) { r <- rob(lm(fm("y", c("BMI", sets[[s]])), d), "BMI")
    add(outcome = ph, table = "B BMI slope", term = "BMI", adjust = s, n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  for (s in list(NULL, "pgs_asd", ND)) { r <- rob(lm(fm("y", c("pgs_cbmi", s)), d), "pgs_cbmi")
    add(outcome = ph, table = "C childhood-BMI score reduced form", term = "pgs_cbmi", adjust = if (is.null(s)) "none" else paste(s, collapse = " "), n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  COV <- c(GRO, ND); Xc <- model.matrix(fm("y", COV), d); Z <- cbind(Xc, as.matrix(d[, ..ADI])); X <- cbind(Xc, BMI = d$BMI); y <- d$y
  Xh <- Z %*% solve(crossprod(Z), crossprod(Z, X)); A <- solve(crossprod(Xh, X)); b <- A %*% crossprod(Xh, y); e <- as.vector(y - X %*% b)
  V <- A %*% crossprod(Xh * e) %*% t(A) * nrow(X) / (nrow(X) - ncol(X)); k <- ncol(X); sarg <- nrow(d) * summary(lm(e ~ Z - 1))$r.squared
  add(outcome = ph, table = "D 2SLS (+ neurodevelopmental covariates)", term = "BMI (IV)", adjust = sprintf("Sargan chi2 %.2f P %.3g", sarg, pchisq(sarg, 2, lower.tail = FALSE)),
      n = nrow(d), b = b[k], se = sqrt(V[k, k]), p = 2 * pnorm(-abs(b[k] / sqrt(V[k, k]))))
}
R <- rbindlist(out, fill = TRUE); fwrite(R, file.path(HERE, "mr_bmi_srs_scq_asd.csv"))
options(width = 250); print(R[, .(outcome, table, term, adjust, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
