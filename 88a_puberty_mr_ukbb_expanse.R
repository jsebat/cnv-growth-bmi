# 88a_puberty_mr_ukbb_expanse.R  (work package E, BMI modifier paper v2; 2026-11-20). Runs on Expanse; UK Biobank individual
# data stay there, only the aggregate CSVs written here are copied back.
# What it does (all UK Biobank, European-ancestry "blue" cluster, QC'd; same filters as ukbb_bodysize10_vs_adult.R and
# ukbb_mr_puberty_pathway.R of 2026-10-21 / 2026-11-03):
#   J  adult height by recalled body size at age 10 (field 1687) in 16p11.2 del / dup, 22q11.2 dup carriers and non-carriers:
#      cell means with 95% CI and the linear trend per category (model and Box-Cox scale of ukbb_bodysize10_vs_adult.R).
#   K  age at menarche (field 2714) by genotype, and adult height per year of later menarche within genotype (fig_supp5_combined.R).
#   M1 one-sample polygenic-score Mendelian randomisation in NON-CARRIER WOMEN with age at menarche (the same women in every step):
#      childhood-obesity score -> age at menarche; age at menarche -> adult height (menarche score as instrument, Wald ratio and
#      2SLS); direct effect of the childhood-obesity score on height with menarche instrumented (2SLS); total; indirect = a x b.
#      Variants: without the height score; plus the adult-BMI score; childhood-BMI score as exposure; a menarche allele score
#      built from Perry 2014 lead SNPs (a GWAS without UK Biobank) as instrument; Box-Cox height scale; women with a body-size answer.
#   M2 two-sample SNP-level MR: SNP-exposure effects from GWAS without UK Biobank (Perry 2014 menarche, Vogelezang 2020 childhood
#      BMI, Bradfield 2019 childhood obesity EUR), SNP-outcome effects estimated here in the non-carrier women; inverse-variance
#      weighted (multiplicative random effects), weighted median, MR-Egger, and multivariable IVW (childhood BMI + menarche -> height).
#   NC negative / positive control outcomes for the scores: maternal smoking around birth (1787), dark hair (1747), comparative
#      height at age 10 (1697); and the reduced forms in men.
# Inputs: genomatrix_UKBB_2025.csv, Molly's PGS_*_genomewide.tsv, jsebat/pgs_menarche/ukbb PGS, field tsv files,
#         INSTR_DIR/instrument_candidates.csv + dosage_chr*.raw (88b + plink2 --export A).
# Outputs (working directory): panel_j_means.csv, panel_j_trends.csv, panel_k_menarche.csv, mr_onesample.csv, mr_twosample.csv,
#         mr_twosample_snps.csv, mr_controls.csv, mr_sample_counts.csv, mr_score_correlations.csv
# null model: none; regression Wald tests (t / z), delta-method SE for ratios, 2SLS asymptotic SE, bootstrap SE for the weighted median.
suppressPackageStartupMessages({ library(data.table); library(MASS) })
WD <- getwd(); options(width = 250); set.seed(20261120)
INS <- Sys.getenv("INSTR_DIR", "../02_instruments")
GM <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"; MS <- "/expanse/projects/sebat1/msacks/UKBB"
AAM <- "/expanse/projects/sebat1/jsebat/pgs_menarche/ukbb/PGS_guan_2025_ext_agemenarche_cojo_genomewide.tsv"
GROUPS <- c("16p11.2 DEL", "No CNV", "16p11.2 DUP", "22q11.2 DUP")
pcn <- paste0("PC", 1:10); pcs <- paste(pcn, collapse = " + ")
gm0 <- fread(GM, select = c("project_40980_IID", "f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", pcn, "Height", "BMI_raw", "age", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"))
setnames(gm0, c("project_40980_IID", "f.eid", "sex_combined", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"), c("IID", "eid", "sex", "l16", "l22"))
gm0 <- gm0[duplicate_drop == "no" & withdrawn_consent == "no" & ancestry == "blue" & sex %in% 1:2 & !is.na(Height) & !is.na(age) & l16 %in% 1:3 & l22 %in% 1:3]
gm0 <- gm0[complete.cases(gm0[, pcn, with = FALSE])]
gm0[, locus := fcase(l16 == 1, "16p11.2 DEL", l16 == 3, "16p11.2 DUP", l22 == 1, "22q11.2 DEL", l22 == 3, "22q11.2 DUP", l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
gm0 <- gm0[!is.na(locus)]; gm0[, age_sq := age^2]; gm0[, sexf := factor(sex)]
rdf <- function(file, field, valid) { x <- fread(file.path(MS, "comparitiveBodySize", file)); setnames(x, "f.eid", "eid")
  inst <- grep(paste0("^f[.]", field, "[.]"), names(x), value = TRUE); x[, v := NA_real_]; for (k in inst) x[is.na(v) & get(k) %in% valid, v := as.numeric(get(k))]; x[!is.na(v), .(eid, v)] }
bs <- rdf("field1687_comparative_body_size_age10.tsv", 1687, 1:3); bs[, step10 := fcase(v == 1, -1, v == 3, 0, v == 2, 1)]
co <- function(m, term) { s <- summary(m)$coefficients; c(est = s[term, 1], se = s[term, 2], p = s[term, 4]) }
row <- function(..., est, se, p = NULL, df = Inf) { if (is.null(p)) p <- 2 * pt(-abs(est / se), df); q <- qt(0.975, df); data.table(..., estimate = est, se = se, lower95 = est - q * se, upper95 = est + q * se, p = p) }

# =========================================================================================== J: height by body size at 10
dj <- merge(gm0[!is.na(BMI_raw)], bs[, .(eid, step10)], by = "eid")            # sample of ukbb_bodysize10_vs_adult.R
covs <- paste("age + age_sq +", pcs)
bcz <- function(y, dd) { z <- rep(NA_real_, nrow(dd))                           # Box-Cox within sex, as in ukbb_bodysize10_vs_adult.R
  for (s in 1:2) { i <- which(dd$sex == s); x <- copy(dd[i])[, y := y[i]]
    b <- boxcox(as.formula(paste("y ~", covs)), data = x, lambda = seq(-2, 2, 0.05), plotit = FALSE); l <- b$x[which.max(b$y)]
    x[, yb := if (abs(l) < 1e-8) log(y) else (y^l - 1) / l]; r <- resid(lm(as.formula(paste("yb ~", covs)), data = x)); z[i] <- r / sd(r) }
  z }
dj[, hZ := bcz(Height, dj)]; dj[, bZ := bcz(BMI_raw, dj)]
dj[, bmi_pct := frank(bZ, ties.method = "average") / .N][, stepAd := fifelse(bmi_pct <= 0.30, -1, fifelse(bmi_pct >= 0.70, 1, 0))]
dj[, size10 := c("Thinner", "About average", "Plumper")[step10 + 2]]
Jm <- dj[locus %in% GROUPS, { n <- .N; m <- mean(hZ); s <- sd(hZ) / sqrt(n); q <- if (n > 1) qt(0.975, n - 1) else NA_real_
  .(n = n, n_women = sum(sex == 2), mean_height_sd = m, se = s, lower95 = m - q * s, upper95 = m + q * s) }, by = .(group = locus, body_size_at_10 = size10)]
fwrite(Jm[order(group, body_size_at_10)], file.path(WD, "panel_j_means.csv"))
Jt <- rbindlist(lapply(GROUPS, function(G) rbindlist(lapply(list(all = 1:2, women = 2, men = 1), function(sx) { x <- dj[locus == G & sex %in% sx]
  rbindlist(lapply(c(step10 = "Body size at age 10", stepAd = "Adult BMI group"), function(lab) { v <- if (lab == "Body size at age 10") "step10" else "stepAd"
    m <- lm(as.formula(paste("hZ ~", v)), x); b <- co(m, v)
    ip <- if (G == "No CNV") NA_real_ else { xi <- dj[locus %in% c(G, "No CNV") & sex %in% sx][, CNV := as.numeric(locus == G)]; summary(lm(as.formula(paste("hZ ~", v, "* CNV")), xi))$coefficients[paste0(v, ":CNV"), 4] }
    cbind(row(group = G, measure = lab, n = nrow(x), est = b[["est"]], se = b[["se"]], df = m$df.residual), p_interaction_vs_noncarriers = ip) })) }), idcol = "sexes")))
fwrite(Jt, file.path(WD, "panel_j_trends.csv")); cat("\n== J trends ==\n"); print(Jt[sexes == "all"])

# =========================================================================================== main table with scores
gm <- gm0[locus %in% GROUPS]
rp <- function(f, nm) { p <- fread(f)[n_chr == 22]; setnames(p, "PGS", nm); p[, c("IID", nm), with = FALSE] }
d <- merge(gm, rp(file.path(MS, "PRS_childhood_obesity/PGS_CHILDHOOD_OBESITY.TRANS_ANCESTRAL_2019.RESULTS_cojo_genomewide.tsv"), "pgs_cobes"), by = "IID")
d <- merge(d, rp(file.path(MS, "PRS_childhood_obesity/PGS_GCST90002409_childhood_BMI_cojo_genomewide.tsv"), "pgs_cbmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_BMI_Yengo_genomewide.tsv"), "pgs_bmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_Height_genomewide.tsv"), "pgs_ht"), by = "IID")
d <- merge(d, rp(AAM, "pgs_aam"), by = "IID")
d[, hZ := { r <- resid(lm(as.formula(paste("Height ~ sexf + age + age_sq +", pcs)), d)); r / sd(r) }]     # as ukbb_mr_puberty_pathway.R
d[, hZbc := bcz(Height, d)]                                                                                 # Box-Cox within sex (check)

# ---- instrument SNP dosages (plink2 --export A) and the Perry 2014 menarche allele score
C <- fread(file.path(INS, "instrument_candidates.csv"), colClasses = list(character = c("a1", "a2", "ea_perry", "oa_perry", "ea_vogelezang", "oa_vogelezang", "ea_bradfield", "oa_bradfield"), numeric = c("frq_perry", "frq_vogelezang", "frq_bradfield")))
G <- NULL
for (i in 1:22) { f <- file.path(INS, sprintf("dosage_chr%d.raw", i)); if (!file.exists(f)) next
  x <- fread(f); x <- x[, c("IID", setdiff(names(x), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))), with = FALSE]
  x <- x[match(d$IID, IID)]; x[, IID := NULL]; G <- if (is.null(G)) x else cbind(G, x) }
hdr <- names(G); C2 <- data.table(col = hdr, snp = sub("_[ACGT]$", "", hdr), counted = sub("^.*_", "", hdr))
C2 <- merge(C2, C, by = "snp"); C2[, other := fifelse(counted == a1, a2, a1)]
G <- as.matrix(G[, C2$col, with = FALSE]); storage.mode(G) <- "double"
cat("dosage matrix:", nrow(G), "x", ncol(G), "; people with no genotype row:", sum(rowSums(!is.na(G)) == 0), "\n")
have <- rowSums(!is.na(G)) > 0
C2[, ukb_frq := colMeans(G[have, , drop = FALSE], na.rm = TRUE) / 2][, miss := colMeans(is.na(G[have, , drop = FALSE]))]
for (j in seq_len(ncol(G))) { k <- is.na(G[, j]) & have; if (any(k)) G[k, j] <- 2 * C2$ukb_frq[j] }        # mean-impute sporadic missing calls
comp <- function(a) chartr("ACGT", "TGCA", a)
align <- function(ea, oa, beta, frq, ca, ot, uf) {                                                        # effect per counted allele
  s <- fifelse(ea == ca & oa == ot, 1, fifelse(ea == ot & oa == ca, -1, fifelse(comp(ea) == ca & comp(oa) == ot, 1, fifelse(comp(ea) == ot & comp(oa) == ca, -1, NA_real_))))
  pal <- ca == comp(ot); s[pal & pmin(uf, 1 - uf) > 0.42] <- NA_real_                                       # ambiguous palindromic SNPs dropped
  fa <- fifelse(s == 1, frq, 1 - frq); s[!is.na(fa) & abs(fa - uf) > 0.2] <- NA_real_                       # allele-frequency check where the GWAS gives one
  beta * s }
for (t in c("perry", "vogelezang", "bradfield")) C2[, paste0("b_", t) := align(get(paste0("ea_", t)), get(paste0("oa_", t)), get(paste0("beta_", t)), get(paste0("frq_", t)), counted, other, ukb_frq)]
cat("harmonised betas available: Perry", sum(!is.na(C2$b_perry)), " Vogelezang", sum(!is.na(C2$b_vogelezang)), " Bradfield", sum(!is.na(C2$b_bradfield)), " of", nrow(C2), "\n")
# LD pruning on UK Biobank genotypes (20,000 random people): greedy by exposure P, r2 < 0.01 with every kept SNP on the chromosome
ii <- sample(which(have), 20000); R2 <- cor(G[ii, ])^2
prune <- function(sel, pcol, bcol) { x <- C2[get(sel) == TRUE & !is.na(get(bcol))][order(get(pcol))]; keep <- character(0)
  for (k in seq_len(nrow(x))) { same <- keep[C2[match(keep, col), chr] == x$chr[k]]; if (!length(same) || all(R2[x$col[k], same] < 0.01)) keep <- c(keep, x$col[k]) }
  keep }
I_pe <- prune("sel_perry", "p_perry", "b_perry"); I_vo <- prune("sel_vogelezang", "p_vogelezang", "b_vogelezang"); I_br <- prune("sel_bradfield", "p_bradfield", "b_bradfield")
cat("instruments after harmonisation and LD pruning: Perry menarche", length(I_pe), " Vogelezang childhood BMI", length(I_vo), " Bradfield childhood obesity", length(I_br), "\n")
d[, score_perry := as.numeric(G[, I_pe, drop = FALSE] %*% C2[match(I_pe, col), b_perry])]
d[!have, score_perry := NA_real_]
for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi", "pgs_ht", "pgs_aam")) { r <- resid(lm(as.formula(paste(v, "~", pcs)), d)); set(d, j = v, value = r / sd(r)) }
ok <- !is.na(d$score_perry); r <- resid(lm(as.formula(paste("score_perry ~", pcs)), d[ok])); set(d, i = which(ok), j = "score_perry", value = r / sd(r))
fwrite(as.data.table(round(cor(d[ok, .(pgs_cobes, pgs_cbmi, pgs_bmi, pgs_aam, score_perry, pgs_ht)]), 4), keep.rownames = "score"), file.path(WD, "mr_score_correlations.csv"))
d[, rowid := .I]

men <- rdf("field2714_age_at_menarche.tsv", 2714, 8:18); setnames(men, "v", "menarche")
w <- merge(d[sex == 2], men, by = "eid")
w <- merge(w, bs[, .(eid, step10)], by = "eid", all.x = TRUE)
ms <- rdf("field1787_maternal_smoking_around_birth.tsv", 1787, 0:1); setnames(ms, "v", "mat_smoke"); w <- merge(w, ms, by = "eid", all.x = TRUE)
hc <- rdf("field1747_hair_colour.tsv", 1747, 1:5); hc[, dark_hair := as.numeric(v %in% 4:5)]; w <- merge(w, hc[, .(eid, dark_hair)], by = "eid", all.x = TRUE)
h10 <- rdf("field1697_comparative_height_size_age10.tsv", 1697, 1:3); h10[, height10 := fcase(v == 1, -1, v == 3, 0, v == 2, 1)]; w <- merge(w, h10[, .(eid, height10)], by = "eid", all.x = TRUE)
cat("women with menarche:\n"); print(w[, .N, by = locus])

# =========================================================================================== K: menarche by genotype
w[, locus0 := relevel(factor(locus, GROUPS), "No CNV")]
mb <- lm(as.formula(paste("menarche ~ locus0 + age +", pcs)), w)
K1 <- rbindlist(lapply(GROUPS, function(Gp) { x <- w[locus == Gp]
  if (Gp == "No CNV") cbind(data.table(quantity = "menarche difference from non-carriers (years)", group = Gp, n_women = nrow(x), mean_age_at_menarche = mean(x$menarche), sd_age_at_menarche = sd(x$menarche)), data.table(estimate = 0, se = NA_real_, lower95 = NA_real_, upper95 = NA_real_, p = NA_real_))
  else { b <- co(mb, paste0("locus0", Gp)); row(quantity = "menarche difference from non-carriers (years)", group = Gp, n_women = nrow(x), mean_age_at_menarche = mean(x$menarche), sd_age_at_menarche = sd(x$menarche), est = b[["est"]], se = b[["se"]], df = mb$df.residual) } }))
K2 <- rbindlist(lapply(GROUPS, function(Gp) { x <- w[locus == Gp]; m <- lm(hZ ~ menarche, x); b <- co(m, "menarche")
  ip <- if (Gp == "No CNV") NA_real_ else summary(lm(hZ ~ menarche * I(locus == Gp), w[locus %in% c(Gp, "No CNV")]))$coefficients["menarche:I(locus == Gp)TRUE", 4]
  cbind(row(quantity = "adult height (SD) per year of later menarche", group = Gp, n_women = nrow(x), mean_age_at_menarche = mean(x$menarche), sd_age_at_menarche = sd(x$menarche), est = b[["est"]], se = b[["se"]], df = m$df.residual), p_interaction_vs_noncarriers = ip) }))
K <- rbind(K1, K2, fill = TRUE); fwrite(K, file.path(WD, "panel_k_menarche.csv")); cat("\n== K ==\n"); print(K)

# =========================================================================================== M1: one-sample score MR, non-carrier women
S <- w[locus == "No CNV"]
cnt <- rbind(data.table(set = "EUR QC with height, PCs, five scores (both sexes, 4 genotype groups)", n = nrow(d)),
             data.table(set = "women with age at menarche, all 4 genotype groups", n = nrow(w)),
             data.table(set = paste("women with age at menarche:", GROUPS), n = sapply(GROUPS, function(g) w[locus == g, .N])),
             data.table(set = "MR sample: non-carrier women with age at menarche", n = nrow(S)),
             data.table(set = "MR sample with a body-size-at-10 answer", n = S[!is.na(step10), .N]),
             data.table(set = "MR sample with genotypes for the SNP instruments", n = S[!is.na(score_perry), .N]),
             data.table(set = "non-carrier men (reduced forms only)", n = d[sex == 1 & locus == "No CNV", .N]))
fwrite(cnt, file.path(WD, "mr_sample_counts.csv")); print(cnt)
tsls <- function(y, x, z, W) {                 # one endogenous regressor x, one instrument z, exogenous W (with intercept)
  X <- cbind(x, W); Z <- cbind(z, W); Xh <- Z %*% solve(crossprod(Z), crossprod(Z, X)); b <- solve(crossprod(Xh, X), crossprod(Xh, y))
  e <- y - X %*% b; s2 <- sum(e^2) / (nrow(X) - ncol(X)); V <- s2 * solve(crossprod(Xh)); list(b = as.numeric(b), se = sqrt(diag(V)), names = colnames(X)) }
onesample <- function(S, expo, inst, adj, out = "hZ", label) {
  S <- S[complete.cases(S[, c(expo, inst, adj, out, "menarche"), with = FALSE])]; n <- nrow(S)
  base <- paste(c(adj, "age", pcn), collapse = " + ")
  ma <- lm(as.formula(paste("menarche ~", expo, "+", base)), S); a <- co(ma, expo)
  mf <- lm(as.formula(paste("menarche ~", inst, "+", base)), S); f <- co(mf, inst)
  mr <- lm(as.formula(paste(out, "~", inst, "+", base)), S); rf <- co(mr, inst)
  wald <- rf[["est"]] / f[["est"]]; wse <- sqrt(rf[["se"]]^2 / f[["est"]]^2 + rf[["est"]]^2 * f[["se"]]^2 / f[["est"]]^4)
  mt <- lm(as.formula(paste(out, "~", expo, "+", base)), S); tot <- co(mt, expo)
  mo <- lm(as.formula(paste(out, "~ menarche +", base)), S); ob <- co(mo, "menarche")
  W <- cbind(1, as.matrix(S[, c(expo, adj, "age", pcn), with = FALSE])); colnames(W)[1:2] <- c("int", expo)
  mff <- lm(as.formula(paste("menarche ~", inst, "+", expo, "+", base)), S); f2 <- co(mff, inst)
  t2 <- tsls(S[[out]], cbind(menarche = S$menarche), cbind(S[[inst]]), W)
  bM <- t2$b[1]; bMse <- t2$se[1]; dir <- t2$b[3]; dirse <- t2$se[3]      # X = [menarche, intercept, exposure score, ...]
  ind <- a[["est"]] * bM; indse <- sqrt(a[["est"]]^2 * bMse^2 + bM^2 * a[["se"]]^2)
  rbind(row(analysis = label, step = "exposure score -> age at menarche", unit = "years per SD of score", n = n, F_stat = (a[["est"]] / a[["se"]])^2, est = a[["est"]], se = a[["se"]]),
        row(analysis = label, step = "menarche instrument -> age at menarche (first stage)", unit = "years per SD of instrument", n = n, F_stat = (f[["est"]] / f[["se"]])^2, est = f[["est"]], se = f[["se"]]),
        row(analysis = label, step = "age at menarche -> adult height (Wald ratio)", unit = "SD height per year later", n = n, F_stat = (f[["est"]] / f[["se"]])^2, est = wald, se = wse),
        row(analysis = label, step = "age at menarche -> adult height (2SLS, exposure score as covariate)", unit = "SD height per year later", n = n, F_stat = (f2[["est"]] / f2[["se"]])^2, est = bM, se = bMse),
        row(analysis = label, step = "age at menarche -> adult height (observational)", unit = "SD height per year later", n = n, F_stat = NA_real_, est = ob[["est"]], se = ob[["se"]]),
        row(analysis = label, step = "exposure score -> adult height, total", unit = "SD height per SD of score", n = n, F_stat = NA_real_, est = tot[["est"]], se = tot[["se"]]),
        row(analysis = label, step = "exposure score -> adult height, direct (menarche instrumented)", unit = "SD height per SD of score", n = n, F_stat = (f2[["est"]] / f2[["se"]])^2, est = dir, se = dirse),
        row(analysis = label, step = "exposure score -> adult height, through menarche (a x b)", unit = "SD height per SD of score", n = n, F_stat = NA_real_, est = ind, se = indse))
}
M1 <- rbind(
  onesample(S, "pgs_cobes", "pgs_aam", "pgs_ht", label = "main: childhood-obesity score, menarche score (Guan 2025), height score adjusted"),
  onesample(S, "pgs_cobes", "pgs_aam", character(0), label = "without height-score adjustment"),
  onesample(S, "pgs_cobes", "pgs_aam", c("pgs_ht", "pgs_bmi"), label = "plus adult-BMI score adjustment"),
  onesample(S, "pgs_cbmi", "pgs_aam", "pgs_ht", label = "childhood-BMI score (Vogelezang 2020) as exposure"),
  onesample(S, "pgs_bmi", "pgs_aam", "pgs_ht", label = "adult-BMI score (Yengo 2018) as exposure"),
  onesample(S, "pgs_cobes", "score_perry", "pgs_ht", label = "menarche instrument without UK Biobank: Perry 2014 allele score"),
  onesample(S, "pgs_cbmi", "score_perry", "pgs_ht", label = "childhood-BMI score and Perry 2014 allele score (no UK Biobank in either)"),
  onesample(S, "pgs_cobes", "score_perry", character(0), label = "Perry 2014 allele score, without height-score adjustment"),
  onesample(S, "pgs_cobes", "score_perry", c("pgs_ht", "pgs_bmi"), label = "Perry 2014 allele score, plus adult-BMI score adjustment"),
  onesample(S, "pgs_cobes", "pgs_aam", "pgs_ht", out = "hZbc", label = "height on the Box-Cox within-sex scale"),
  onesample(S[!is.na(step10)], "pgs_cobes", "pgs_aam", "pgs_ht", label = "women with a body-size-at-10 answer"),
  onesample(S[!is.na(score_perry)], "pgs_cobes", "pgs_aam", "pgs_ht", label = "women with SNP-instrument genotypes"))
# per category of recalled body size at 10: Wald ratios with the childhood-obesity score as instrument for body size
S10 <- S[!is.na(step10)]; base <- paste(c("pgs_ht", "age", pcn), collapse = " + ")
for (ex in c("pgs_cobes", "pgs_cbmi")) { lab <- paste0("per category of recalled body size at 10 (instrument: ", if (ex == "pgs_cobes") "childhood-obesity" else "childhood-BMI", " score)")
  f10 <- co(lm(as.formula(paste("step10 ~", ex, "+", base)), S10), ex)
  wr <- function(num) { e <- num[["est"]] / f10[["est"]]; c(e, sqrt(num[["se"]]^2 / f10[["est"]]^2 + num[["est"]]^2 * f10[["se"]]^2 / f10[["est"]]^4)) }
  a10 <- wr(co(lm(as.formula(paste("menarche ~", ex, "+", base)), S10), ex)); t10 <- wr(co(lm(as.formula(paste("hZ ~", ex, "+", base)), S10), ex))
  o1 <- co(lm(as.formula(paste("menarche ~ step10 + age +", pcs)), S10), "step10"); o2 <- co(lm(as.formula(paste("hZ ~ step10 + age +", pcs)), S10), "step10")
  Ff <- (f10[["est"]] / f10[["se"]])^2
  M1 <- rbind(M1, row(analysis = lab, step = "score -> body size at 10 (first stage)", unit = "categories per SD of score", n = nrow(S10), F_stat = Ff, est = f10[["est"]], se = f10[["se"]]),
              row(analysis = lab, step = "body size at 10 -> age at menarche (Wald ratio)", unit = "years per category plumper", n = nrow(S10), F_stat = Ff, est = a10[1], se = a10[2]),
              row(analysis = lab, step = "body size at 10 -> adult height, total (Wald ratio)", unit = "SD height per category plumper", n = nrow(S10), F_stat = Ff, est = t10[1], se = t10[2]),
              row(analysis = lab, step = "body size at 10 -> age at menarche (observational)", unit = "years per category plumper", n = nrow(S10), F_stat = NA_real_, est = o1[["est"]], se = o1[["se"]]),
              row(analysis = lab, step = "body size at 10 -> adult height (observational)", unit = "SD height per category plumper", n = nrow(S10), F_stat = NA_real_, est = o2[["est"]], se = o2[["se"]])) }
M1[, R2_first_stage := F_stat / (F_stat + n)]
fwrite(M1, file.path(WD, "mr_onesample.csv")); cat("\n== M1 ==\n"); print(M1[, .(analysis = substr(analysis, 1, 60), step = substr(step, 1, 62), n, F = round(F_stat), est = signif(estimate, 3), se = signif(se, 2), p = signif(p, 2))], nrows = 200)

# =========================================================================================== controls
ctl <- function(D, out, sc, lab, sample) { D <- D[!is.na(get(out)) & !is.na(get(sc))]; m <- lm(as.formula(paste(out, "~", sc, "+ pgs_ht + age +", pcs)), D); b <- co(m, sc)
  row(sample = sample, outcome = lab, score = sc, n = nrow(D), est = b[["est"]], se = b[["se"]], df = m$df.residual) }
SC <- c("pgs_cobes", "pgs_cbmi", "pgs_aam", "score_perry", "pgs_bmi")
Mn <- d[sex == 1 & locus == "No CNV"]
NCt <- rbind(rbindlist(lapply(SC, function(s) rbind(
    ctl(S, "mat_smoke", s, "negative control: mother smoked around birth (0/1)", "non-carrier women (MR sample)"),
    ctl(S, "dark_hair", s, "negative control: dark brown or black hair (0/1)", "non-carrier women (MR sample)"),
    ctl(S, "height10", s, "positive control: comparative height at age 10 (-1 shorter, 0 average, +1 taller)", "non-carrier women (MR sample)"),
    ctl(S, "step10", s, "exposure check: comparative body size at age 10 (-1 thinner, 0 average, +1 plumper)", "non-carrier women (MR sample)"),
    ctl(S, "menarche", s, "age at menarche (years)", "non-carrier women (MR sample)"),
    ctl(S, "hZ", s, "adult height (SD)", "non-carrier women (MR sample)"),
    ctl(Mn, "hZ", s, "adult height (SD)", "non-carrier men")))))
fwrite(NCt, file.path(WD, "mr_controls.csv")); cat("\n== controls ==\n"); print(NCt[, .(sample = substr(sample, 1, 20), outcome = substr(outcome, 1, 50), score, n, est = signif(estimate, 3), p = signif(p, 2))], nrows = 100)

# =========================================================================================== M2: two-sample SNP-level MR
Sg <- S[!is.na(score_perry)]; Gs <- G[Sg$rowid, , drop = FALSE]
Cm <- cbind(1, as.matrix(Sg[, c("age", pcn), with = FALSE])); Q <- qr(Cm); dfres <- nrow(Sg) - ncol(Cm) - 1
snpreg <- function(y) { ry <- qr.resid(Q, y); t(apply(Gs, 2, function(g) { rg <- qr.resid(Q, g); b <- sum(rg * ry) / sum(rg^2); c(b, sqrt(sum((ry - b * rg)^2) / dfres / sum(rg^2))) })) }
om <- snpreg(Sg$menarche); oh <- snpreg(Sg$hZ)
C2[, `:=`(ukb_menarche_b = om[, 1], ukb_menarche_se = om[, 2], ukb_height_b = oh[, 1], ukb_height_se = oh[, 2])]
C2[, `:=`(inst_perry = col %in% I_pe, inst_vogelezang = col %in% I_vo, inst_bradfield = col %in% I_br)]
fwrite(C2[, .(snp, chr, pos, counted_allele = counted, other_allele = other, ukb_freq = ukb_frq, inst_perry, inst_vogelezang, inst_bradfield, b_perry, se_perry, p_perry, b_vogelezang, se_vogelezang, p_vogelezang,
              b_bradfield, se_bradfield, p_bradfield, ukb_menarche_b, ukb_menarche_se, ukb_height_b, ukb_height_se)], file.path(WD, "mr_twosample_snps.csv"))
ivw <- function(bx, by, sy) { w <- 1 / sy^2; b <- sum(w * bx * by) / sum(w * bx^2); k <- length(bx); Qs <- sum(w * (by - b * bx)^2); phi <- max(1, Qs / (k - 1)); se <- sqrt(phi / sum(w * bx^2))
  list(b = b, se = se, p = 2 * pnorm(-abs(b / se)), Q = Qs, Qp = pchisq(Qs, k - 1, lower.tail = FALSE)) }
egger <- function(bx, by, sy) { s <- sign(bx); bx <- bx * s; by <- by * s; m <- summary(lm(by ~ bx, weights = 1 / sy^2)); sc <- min(1, m$sigma); k <- length(bx)
  b <- m$coefficients[2, 1]; se <- m$coefficients[2, 2] / sc; i0 <- m$coefficients[1, 1]; i0se <- m$coefficients[1, 2] / sc
  list(b = b, se = se, p = 2 * pt(-abs(b / se), k - 2), int = i0, int_se = i0se, int_p = 2 * pt(-abs(i0 / i0se), k - 2)) }
wmed0 <- function(bx, by, wt) { r <- by / bx; o <- order(r); r <- r[o]; wt <- wt[o]; cw <- cumsum(wt) - 0.5 * wt; cw <- cw / sum(wt); if (!any(cw < 0.5)) return(r[1]); k <- max(which(cw < 0.5)); if (k == length(r)) return(r[k])
  r[k] + (r[k + 1] - r[k]) * (0.5 - cw[k]) / (cw[k + 1] - cw[k]) }
wmed <- function(bx, sx, by, sy) { wt <- bx^2 / sy^2; b <- wmed0(bx, by, wt); bs <- replicate(2000, wmed0(rnorm(length(bx), bx, sx), rnorm(length(by), by, sy), wt)); se <- sd(bs); list(b = b, se = se, p = 2 * pnorm(-abs(b / se))) }
two <- function(inst, bcol, scol, ocol, label, unit) { x <- C2[col %in% inst]; bx <- x[[bcol]]; sx <- x[[scol]]; by <- x[[paste0(ocol, "_b")]]; sy <- x[[paste0(ocol, "_se")]]
  k <- length(bx); Fm <- mean((bx / sx)^2); a <- ivw(bx, by, sy); e <- egger(bx, by, sy); m <- wmed(bx, sx, by, sy)
  mk <- function(meth, r, ...) cbind(data.table(analysis = label, method = meth, unit = unit, n_snps = k, mean_F = Fm, n_outcome = nrow(Sg), estimate = r$b, se = r$se, lower95 = r$b - 1.96 * r$se, upper95 = r$b + 1.96 * r$se, p = r$p), ...)
  rbind(mk("Inverse-variance weighted", a, data.table(Q = a$Q, Q_p = a$Qp)), mk("Weighted median", m), mk("MR-Egger", e, data.table(egger_intercept = e$int, egger_intercept_se = e$int_se, egger_intercept_p = e$int_p)), fill = TRUE) }
M2 <- rbind(
  two(I_vo, "b_vogelezang", "se_vogelezang", "ukb_menarche", "childhood BMI (Vogelezang 2020) -> age at menarche", "years per SD of childhood BMI"),
  two(I_br, "b_bradfield", "se_bradfield", "ukb_menarche", "childhood obesity (Bradfield 2019, EUR) -> age at menarche", "years per log-odds of childhood obesity"),
  two(I_pe, "b_perry", "se_perry", "ukb_height", "age at menarche (Perry 2014) -> adult height", "SD height per year later"),
  two(I_vo, "b_vogelezang", "se_vogelezang", "ukb_height", "childhood BMI (Vogelezang 2020) -> adult height, total", "SD height per SD of childhood BMI"),
  two(I_br, "b_bradfield", "se_bradfield", "ukb_height", "childhood obesity (Bradfield 2019, EUR) -> adult height, total", "SD height per log-odds of childhood obesity"), fill = TRUE)
# multivariable IVW: height on childhood adiposity and menarche together (union of instruments, pruned jointly)
mv <- function(Ia, bA, labA, unitA) { cand <- C2[col %in% union(Ia, I_pe) & !is.na(get(bA)) & !is.na(b_perry)]; cand[, pmin_ := pmin(get(sub("b_", "p_", bA)), p_perry, na.rm = TRUE)]; cand <- cand[order(pmin_)]; keep <- character(0)
  for (k in seq_len(nrow(cand))) { same <- keep[C2[match(keep, col), chr] == cand$chr[k]]; if (!length(same) || all(R2[cand$col[k], same] < 0.01)) keep <- c(keep, cand$col[k]) }
  x <- C2[col %in% keep]; m <- summary(lm(x$ukb_height_b ~ 0 + x[[bA]] + x$b_perry, weights = 1 / x$ukb_height_se^2)); sc <- min(1, m$sigma); k <- nrow(x)
  est <- m$coefficients[, 1]; se <- m$coefficients[, 2] / sc
  data.table(analysis = paste0("multivariable: ", labA, " + age at menarche (Perry 2014) -> adult height"), method = "Multivariable inverse-variance weighted",
             unit = c(paste0("SD height per ", unitA, ", direct"), "SD height per year later, adjusted for childhood adiposity"), n_snps = k, mean_F = NA_real_, n_outcome = nrow(Sg),
             estimate = est, se = se, lower95 = est - 1.96 * se, upper95 = est + 1.96 * se, p = 2 * pt(-abs(est / se), k - 2)) }
M2 <- rbind(M2, mv(I_vo, "b_vogelezang", "childhood BMI (Vogelezang 2020)", "SD of childhood BMI"), mv(I_br, "b_bradfield", "childhood obesity (Bradfield 2019, EUR)", "log-odds of childhood obesity"), fill = TRUE)
# replication of the external menarche effects in UK Biobank (instrument validity)
x <- C2[col %in% I_pe]; rr <- summary(lm(ukb_menarche_b ~ 0 + b_perry, x, weights = 1 / ukb_menarche_se^2))$coefficients
M2 <- rbind(M2, data.table(analysis = "check: Perry 2014 SNP effects on menarche replicate in UK Biobank women", method = "weighted slope through origin", unit = "UKB years per Perry year", n_snps = nrow(x), mean_F = mean((x$b_perry / x$se_perry)^2), n_outcome = nrow(Sg),
                           estimate = rr[1, 1], se = rr[1, 2], lower95 = rr[1, 1] - 1.96 * rr[1, 2], upper95 = rr[1, 1] + 1.96 * rr[1, 2], p = rr[1, 4]), fill = TRUE)
fwrite(M2, file.path(WD, "mr_twosample.csv")); cat("\n== M2 ==\n"); print(M2[, .(analysis = substr(analysis, 1, 70), method = substr(method, 1, 26), n_snps, F = round(mean_F), est = signif(estimate, 3), se = signif(se, 2), p = signif(p, 2), Qp = signif(Q_p, 2), eggP = signif(egger_intercept_p, 2))], nrows = 100)
cat("done\n")
