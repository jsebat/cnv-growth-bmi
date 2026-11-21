# 92_v2_srs_multimodal.R
# Ported to the BMI modifier paper (JS / Claude 2026-11-21) from COMET_ASD_ARPAH/analysis/A10_srs_multimodal.R (ARPA-H session, same day), unchanged
# except the output folders: aggregates -> run_v2_2026-11-20/multimodal/, participant-level -> DATA/v2_autism_cache/ (not in git).
# The model rules are the paper's: adults Box-Cox, children WHO / UK90; one intercept per cohort within CNV; per CNV + pooled with
# heterogeneity; SSC carriers dropped (ascertained probands); family-grouped cross-validation and family bootstrap.
# COMET-ASD Figure 2 k-l (2026-11-21, Jon): a model of autistic traits (SRS) in CNV carriers built from growth, sex and genetic
# background, and how well it separates children with high and low SRS.  "Incorporate as much of the data as we can (CNV, height,
# BMI, head circumference, sex, PGS, etc.) ... show our ability to predict a high vs a low group."
# Revised the same day after the rule check by the Sacks BMI-paper session (Jon: "Don't change any rules! Tell arpah if they are
# breaking any") and Jon's request for an importance measure ("loss from drop analysis"):
#   (1) adults are Box-Cox scaled (Jon's rule: adults Box-Cox, children against the WHO / UK90 norms), (2) one analysis per CNV with a
#   heterogeneity test, (3) the 17 SSC carriers are dropped (ascertained as autistic probands; cells of 8, 8 and 1 cannot carry their
#   own intercept in cross-validation), (4) families, not children, are resampled, and the R2 range is labelled as the spread over
#   cross-validation repeats, (5) drop analysis.
# Sample: CNV carriers with an SRS in the combined SRS sample of the Sacks BMI paper (cnv-growth-bmi 60/51b; cache
#   DATA/v2_autism_cache/base_2026-11-20.RData, object C): G2MH, Simons VIP, Simons Searchlight (SSC dropped, see 3).
#   Outcome: SRS total, Box-Cox Z against SSC unaffected relatives of the same form.  BMI Z as in C (Sacks paper).
# Growth measures at the visit nearest the SRS (within 3 y):
#   height and head circumference; children (height < 19 y, head circumference < 18 y): WHO 2006/2007 and UK90 Z (00_scoring.R);
#   adults: Box-Cox scaling of 00_scoring.R score_adult_bc() against the non-carrier adults of the clinical cohorts (latest visit of
#   each non-carrier in analysis_table_chop_unaffctl_2026-11-20.rds: SSC parents and siblings, Simons VIP non-carriers), with the
#   exponents of the paper's growth run (17_main_with_ukbb.R in fig2_unaffctl_2026-11-20, from UK Biobank adults: height 0.65 men,
#   0.60 women; head circumference fixed at 1), residualised on age + age^2 within sex.
#   G2MH and Searchlight measures from that table (linked by the cache id map); Simons VIP from its v10 growth_measurements.csv
#   (the growth table uses another Simons VIP id).
# Predictors: genotype x cohort (one intercept per cell), age (natural spline, 3 df), sex, BMI, height, head circumference, ASD PGS
#   (April 2026 autism score of the Sacks paper, standardised in these carriers), BMI x sex, BMI x ASD PGS.  Missing height, head
#   circumference or ASD PGS: value 0 plus a missing indicator, so every carrier is used (a prediction device; the coefficient of a
#   measured predictor is estimated in the carriers who have it).
#   Note: BMI x sex and BMI x ASD PGS here are WITHIN-CARRIER terms, not the carrier-versus-relative three-way terms of the paper.
# Analyses
#   1. Coefficients: each factor alone (+ genotype x cohort + age) and all together, pooled.  Per CNV: the full model within each CNV
#      (cohort intercepts within the CNV), and Cochran's Q across the four CNVs for each term.
#   2. Prediction: 10-fold cross-validation with families kept in one fold, 50 repeats.  Benchmark = genotype x cohort + age; full =
#      benchmark + every factor.  Out-of-sample R2 per repeat (range = spread over the repeats, not a confidence interval).
#   3. Drop analysis (importance): in the same folds, the full model refitted without one factor at a time (a factor goes with its
#      missing indicator and its interactions: BMI = BMI, BMI x sex, BMI x ASD PGS; sex = sex, BMI x sex; ASD PGS = PGS, indicator,
#      BMI x ASD PGS), and without each interaction alone.  Loss = R2(full) - R2(reduced), per repeat; mean, spread over repeats, and
#      share of the full model's gain over the benchmark.
#   4. Separation within genotype x cohort and age: children ranked by what the added factors contribute (full minus benchmark
#      prediction, both out of sample); SRS minus the benchmark prediction, and the share with an autism diagnosis, by fifth; top
#      minus bottom fifth with a family bootstrap 95% CI (2,000 draws).
#   # null model: none.  Coefficients: Wald tests; heterogeneity: Cochran's Q.  Prediction: out-of-sample comparison with the
#   # benchmark.  Exploratory; no correction for multiple testing.
# Outputs: aggregate CSVs in figures/ (srs_multimodal_*.csv, SRS_multimodal_numbers.txt); participant-level predictions in
#   COMET_ASD_ARPAH_data/srs_multimodal_panels.rds (not in git).
suppressPackageStartupMessages({ library(data.table); library(splines) })
REPO <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/COMET_ASD_ARPAH"; PRIV <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/COMET_ASD_ARPAH_data"
MMOUT <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/run_v2_2026-11-20/multimodal"; dir.create(MMOUT, showWarnings = FALSE)
MMPRIV <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/v2_autism_cache"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
SV <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0"
source(file.path(HERE, "00_scoring.R")); set.seed(20261121); options(width = 220)
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
e <- new.env(); load(file.path(DATA, "v2_autism_cache", "base_2026-11-20.RData"), envir = e)
C <- as.data.table(e$C)[is_cnv == 1 & pheno == "SRS" & locus %in% LOCI, .(id, cohort, locus, sex, age, BMI, y, form2, asd = asd2)]
n_ssc <- C[cohort == "SSC", .N]; C <- C[cohort != "SSC"]
pg <- unique(as.data.table(e$XG)[is_cnv == 1 & pheno == "SRS" & !is.na(pgs), .(id, pgs)], by = "id"); C <- merge(C, pg, by = "id", all.x = TRUE)
# other polygenic scores of the paper (86_v2_cnv_bmi_pgs_stats.R object XM: childhood obesity, adult BMI and height scores, residualised
# on ancestry PCs and standardised as in the paper; genotyped G2MH and Simons VIP carriers)
e86 <- new.env(); load(file.path(DATA, "v2_autism_cache", "pgs_stats_2026-11-20.RData"), envir = e86)
op <- unique(as.data.table(e86$XM)[is_cnv == 1 & pheno == "SRS", .(id, pgs_cobes, pgs_bmi, pgs_ht)], by = "id"); C <- merge(C, op, by = "id", all.x = TRUE)
stopifnot(!anyDuplicated(C$id))
# ---- growth measures nearest the SRS visit (within 3 y), with the paper's child / adult scoring
T0 <- as.data.table(readRDS(file.path(DATA, "analysis_table_chop_unaffctl_2026-11-20.rds")))
REF <- T0[locus == "No CNV" & latest %in% TRUE, .(age, sex = as.character(sex), h = HEIGHT_raw, hc = HEADCIRCUM_raw)]
LAMBDA <- list(height = c(male = 0.65, female = 0.60), hc = c(male = 1, female = 1))   # paper's growth run (fig2d.55921993.log)
zgrow <- function(v, age, sex, item) {   # children: external reference; adults: Box-Cox against REF, exponent of the paper
  z <- if (item == "height") score_height(v, age, sex) else score_hc(v, age, sex); amin <- if (item == "height") 19 else 18; info <- list()
  for (s in c("male", "female")) { i <- which(sex == s & !is.na(age) & age >= amin); if (!length(i)) next
    rv <- if (item == "height") REF$h else REF$hc; r <- which(REF$sex == s & REF$age >= amin & !is.na(rv))
    zb <- score_adult_bc(c(rv[r], v[i]), c(REF$age[r], age[i]), rep(s, length(r) + length(i)), rep("Clinical", length(r) + length(i)),
                         c(rep(TRUE, length(r)), rep(FALSE, length(i))), min_age = amin, lambda_fixed = LAMBDA[[item]][[s]])
    z[i] <- zb[length(r) + seq_along(i)]; info[[s]] <- attr(zb, "info") }
  attr(z, "info") <- do.call(rbind, info); z }
T <- T0[, .(tid = id, tage = age, h = HEIGHT_raw, hc = HEADCIRCUM_raw, family)]
im <- as.data.table(e$idmap)
k <- unique(rbind(C[, .(id, tid = paste0(cohort, ":", id))], merge(C[, .(id)], im[, .(id = sampleID, tid = id)], by = "id"), merge(C[, .(id)], im[!is.na(GUID), .(id = GUID, tid = id)], by = "id")))
g1 <- merge(merge(k, T, by = "tid", allow.cartesian = TRUE), C[, .(id, cohort)], by = "id")[cohort != "SVIP"]
sv <- fread(file.path(SV, "growth_measurements.csv"))[, .(id = individual, tage = eval_age_months / 12, h = growth_measurements_raw.height, hc = growth_measurements_raw.head_circumference, family = as.character(family))]
g <- rbind(g1[, .(id, tage, h, hc, family = as.character(family))], sv[id %in% C[cohort == "SVIP", id]])
g <- merge(g, C[, .(id, sage = age, sex)], by = "id")[, d := abs(tage - sage)][d <= 3]
nh <- g[!is.na(h)][order(id, d)][, .SD[1], by = id]; nc <- g[!is.na(hc)][order(id, d)][, .SD[1], by = id]
nh[, hgt := zgrow(h, tage, sex, "height")]; nc[, hcz := zgrow(hc, tage, sex, "hc")]
cat("adult scaling, height:\n"); print(attr(zgrow(nh$h, nh$tage, nh$sex, "height"), "info")); cat("adult scaling, head circumference:\n"); print(attr(zgrow(nc$hc, nc$tage, nc$sex, "hc"), "info"))
fam <- unique(g[!is.na(family) & family != "", .(id, family)], by = "id")
C <- merge(merge(merge(C, nh[, .(id, hgt, hgt_age = tage)], by = "id", all.x = TRUE), nc[, .(id, hcz, hc_age = tage)], by = "id", all.x = TRUE), fam, by = "id", all.x = TRUE)
C[!is.na(hgt) & abs(hgt) > 5, hgt := NA_real_]; C[!is.na(hcz) & abs(hcz) > 5, hcz := NA_real_]
C[is.na(family), family := fifelse(grepl("^[0-9]+\\.", id), sub("\\..*$", "", id), id)]   # Simons VIP / Searchlight ids carry the family before the dot
C[, `:=`(cg = factor(paste(cohort, locus)), male = as.numeric(sex == "male"), pgs_s = (pgs - mean(pgs, na.rm = TRUE)) / sd(pgs, na.rm = TRUE))]
for (v in c("pgs_cobes", "pgs_bmi", "pgs_ht")) C[, (v) := (get(v) - mean(get(v), na.rm = TRUE)) / sd(get(v), na.rm = TRUE)]
C[, opgs_m := as.numeric(is.na(pgs_cobes) | is.na(pgs_bmi) | is.na(pgs_ht))]; for (v in c("pgs_cobes", "pgs_bmi", "pgs_ht")) C[opgs_m == 1, (v) := 0]
for (v in c("hgt", "hcz", "pgs_s")) { C[, paste0(v, "_m") := as.numeric(is.na(get(v)))]; C[is.na(get(v)), (v) := 0] }
C[, `:=`(form2 = factor(form2), locus = factor(locus, LOCI))]
smp <- C[, .(carriers = .N, adults_19plus = sum(age >= 19), adult_form = sum(form2 == "adult"), with_other_PGS = sum(opgs_m == 0), with_height = sum(hgt_m == 0), with_head_circumference = sum(hcz_m == 0), with_ASD_PGS = sum(pgs_s_m == 0),
             with_all_three = sum(hgt_m + hcz_m + pgs_s_m == 0), autism_known = sum(!is.na(asd)), autism_yes = sum(asd %in% 1), families = uniqueN(family)), by = .(cohort, locus)][order(locus, cohort)]
smp <- rbind(smp, C[, .(cohort = "All", locus = "All", carriers = .N, adults_19plus = sum(age >= 19), adult_form = sum(form2 == "adult"), with_other_PGS = sum(opgs_m == 0), with_height = sum(hgt_m == 0), with_head_circumference = sum(hcz_m == 0), with_ASD_PGS = sum(pgs_s_m == 0),
             with_all_three = sum(hgt_m + hcz_m + pgs_s_m == 0), autism_known = sum(!is.na(asd)), autism_yes = sum(asd %in% 1), families = uniqueN(family))])
cat("== sample (SSC carriers dropped:", n_ssc, ")\n"); print(smp)
# ---- 1. coefficients: each factor alone, all together, and per CNV
# benchmark: genotype x cohort, age, and SRS form (child / adult, with its own age slope: the form differs in informant and age)
BASE <- "y ~ cg + ns(age, 3) + form2 + form2:age"
TERMS <- list(Sex = c("male", "male"), BMI = c("BMI", "BMI"), Height = c("hgt + hgt_m", "hgt"), `Head circumference` = c("hcz + hcz_m", "hcz"), `ASD PGS` = c("pgs_s + pgs_s_m", "pgs_s"),
              `Childhood obesity PGS` = c("pgs_cobes + opgs_m", "pgs_cobes"), `Adult BMI PGS` = c("pgs_bmi + opgs_m", "pgs_bmi"), `Height PGS` = c("pgs_ht + opgs_m", "pgs_ht"))
ADD <- c("male", "BMI", "hgt", "hgt_m", "hcz", "hcz_m", "pgs_s", "pgs_s_m", "pgs_cobes", "pgs_bmi", "pgs_ht", "opgs_m", "BMI:male", "BMI:pgs_s", "BMI:locus")
fml <- function(add, base = BASE) as.formula(paste(base, "+", paste(add, collapse = " + ")))
cf1 <- function(m, term) { s <- summary(m)$coefficients; rv <- paste(rev(strsplit(term, ":")[[1]]), collapse = ":"); if (!term %in% rownames(s) && rv %in% rownames(s)) term <- rv
  if (!term %in% rownames(s)) return(data.table(est = NA_real_, se = NA_real_, p = NA_real_)); data.table(est = s[term, 1], se = s[term, 2], p = s[term, 4]) }
alone <- rbindlist(lapply(names(TERMS), function(nm) cbind(factor = nm, model = "alone", cf1(lm(as.formula(paste(BASE, "+", TERMS[[nm]][1])), C), TERMS[[nm]][2]))))
mf <- lm(fml(ADD), C)
LAB <- c(male = "Sex", BMI = "BMI (16p11.2 del)", hgt = "Height", hcz = "Head circumference", pgs_s = "ASD PGS", pgs_cobes = "Childhood obesity PGS", pgs_bmi = "Adult BMI PGS", pgs_ht = "Height PGS",
         `BMI:male` = "BMI x sex (male)", `BMI:pgs_s` = "BMI x ASD PGS", `BMI:locus16p11.2 dup` = "BMI x CNV: 16p11.2 dup vs del", `BMI:locus22q11.2 del` = "BMI x CNV: 22q11.2 del vs 16p11.2 del",
         `BMI:locus22q11.2 dup` = "BMI x CNV: 22q11.2 dup vs 16p11.2 del")
TERMS$BMI <- c("BMI", "BMI")
together <- rbindlist(lapply(names(LAB), function(t) cbind(factor = LAB[[t]], model = "together", cf1(mf, t))))
CO <- rbind(alone, together)[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
CO[, n := ifelse(factor == "Height", sum(C$hgt_m == 0), ifelse(factor == "Head circumference", sum(C$hcz_m == 0), ifelse(grepl("ASD PGS", factor), sum(C$pgs_s_m == 0),
           ifelse(grepl("PGS", factor), sum(C$opgs_m == 0), nrow(C)))))]
# the CNV-specific BMI slopes of the full model (BMI main effect = 16p11.2 del; slope of each CNV = BMI + BMI x CNV)
bb <- coef(mf); okb <- !is.na(bb); bb <- bb[okb]; vv <- vcov(mf)[names(bb), names(bb)]
bl <- function(L) { t <- if (L == "16p11.2 del") "BMI" else c("BMI", paste0("BMI:locus", L)); w <- setNames(rep(0, length(bb)), names(bb)); w[t] <- 1
  est <- sum(w * bb); se <- sqrt(drop(t(w) %*% vv %*% w)); data.table(factor = paste("BMI slope,", L), model = "together", est = est, se = se, p = 2 * pnorm(-abs(est / se))) }
CO <- rbind(CO, rbindlist(lapply(LOCI, bl))[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se, n = as.integer(table(C$locus)[LOCI]))])
cat("CNV x BMI block, F test:", anova(lm(fml(setdiff(ADD, "BMI:locus")), C), mf)$`Pr(>F)`[2], "\n")
cat("\n== pooled coefficients (SD of SRS per unit: sex = male minus female; others per SD)\n"); print(CO, digits = 3)
PER <- rbindlist(lapply(LOCI, function(G) { x <- C[locus == G]; base <- paste(if (uniqueN(x$cohort) > 1) "y ~ cohort + ns(age, 3)" else "y ~ ns(age, 3)", if (uniqueN(x$form2) > 1) "+ form2 + form2:age" else "")
  add <- ADD[sapply(ADD, function(a) { v <- strsplit(a, ":")[[1]]; all(sapply(v, function(z) uniqueN(x[[z]]) > 1)) })]
  m <- lm(fml(add, base), x); lb <- LAB[!grepl("^BMI:locus", names(LAB))]; lb["BMI"] <- "BMI"
  rbindlist(lapply(names(lb), function(t) cbind(locus = G, n = nrow(x), cohorts = paste(sort(unique(x$cohort)), collapse = "+"), factor = lb[[t]], cf1(m, t)))) }))
het <- PER[, { ok <- !is.na(est) & !is.na(se); if (sum(ok) < 2) .(Q = NA_real_, df = NA_integer_, p_het = NA_real_) else { w <- 1 / se[ok]^2; b <- sum(w * est[ok]) / sum(w); q <- sum(w * (est[ok] - b)^2)
  .(Q = q, df = sum(ok) - 1L, p_het = pchisq(q, sum(ok) - 1, lower.tail = FALSE)) } }, by = factor]
PER <- merge(PER, het, by = "factor", sort = FALSE)
cat("\n== per CNV (full model within each CNV) and heterogeneity across CNVs\n"); print(dcast(PER, factor + p_het ~ locus, value.var = "est"), digits = 3); print(PER, digits = 3)
# ---- 2-3. prediction and drop analysis, in the same folds
DROPS <- list(BMI = c("BMI", "BMI:male", "BMI:pgs_s", "BMI:locus"), Sex = c("male", "BMI:male"), Height = c("hgt", "hgt_m"), `Head circumference` = c("hcz", "hcz_m"),
              `ASD PGS` = c("pgs_s", "pgs_s_m", "BMI:pgs_s"), `Adiposity and height PGS` = c("pgs_cobes", "pgs_bmi", "pgs_ht", "opgs_m"),
              `BMI x CNV` = "BMI:locus", `BMI x sex (male)` = "BMI:male", `BMI x ASD PGS` = "BMI:pgs_s")
MODELS <- c(list(benchmark = as.formula(BASE), full = fml(ADD), benchmark_no_form = as.formula("y ~ cg + ns(age, 3)")), lapply(DROPS, function(d) fml(setdiff(ADD, d))))
fams <- unique(C$family); R <- 50
P <- lapply(MODELS, function(z) matrix(NA_real_, nrow(C), R)); r2 <- matrix(NA_real_, R, length(MODELS), dimnames = list(NULL, names(MODELS)))
r2f <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
for (r in seq_len(R)) { fo <- setNames(sample(rep(1:10, length.out = length(fams))), fams)[C$family]
  for (kk in 1:10) { tr <- C[fo != kk]; te <- fo == kk
    for (nm in names(MODELS)) P[[nm]][te, r] <- suppressWarnings(predict(lm(MODELS[[nm]], tr), C[te])) }
  for (nm in names(MODELS)) r2[r, nm] <- r2f(C$y, P[[nm]][, r]) }
C[, `:=`(pred_b = rowMeans(P$benchmark), pred_f = rowMeans(P$full))]
gain <- r2[, "full"] - r2[, "benchmark"]
fg <- r2[, "benchmark"] - r2[, "benchmark_no_form"]
R2 <- data.table(model = c("benchmark without SRS form: genotype x cohort + age", "benchmark: genotype x cohort + age + SRS form (+ form x age)",
                           "full: + sex, BMI, height, head circumference, ASD PGS, adiposity and height PGS, BMI x sex, BMI x ASD PGS, BMI x CNV", "gain (full - benchmark)", "form gain (benchmark - benchmark without form)"),
                 r2 = c(mean(r2[, "benchmark_no_form"]), mean(r2[, "benchmark"]), mean(r2[, "full"]), mean(gain), mean(fg)),
                 lo = c(quantile(r2[, "benchmark_no_form"], .025), quantile(r2[, "benchmark"], .025), quantile(r2[, "full"], .025), quantile(gain, .025), quantile(fg, .025)),
                 hi = c(quantile(r2[, "benchmark_no_form"], .975), quantile(r2[, "benchmark"], .975), quantile(r2[, "full"], .975), quantile(gain, .975), quantile(fg, .975)),
                 range_is = "2.5-97.5% of 50 CV repeats (not a confidence interval)")
DR <- rbindlist(lapply(names(DROPS), function(nm) { l <- r2[, "full"] - r2[, nm]
  data.table(factor = nm, removed_terms = paste(DROPS[[nm]], collapse = ", "), loss_r2 = mean(l), lo = quantile(l, .025), hi = quantile(l, .975), share_of_gain = mean(l) / mean(gain),
             range_is = "2.5-97.5% of 50 CV repeats (not a confidence interval)") }))[order(-loss_r2)]
cat("\n== out-of-sample R2\n"); print(R2, digits = 3); cat("\n== drop analysis: loss of out-of-sample R2 when a factor is removed\n"); print(DR, digits = 3)
# ---- 4. separation within genotype x cohort and age
fifth <- function(x) cut(rank(x, ties.method = "first"), quantile(seq_along(x), 0:5 / 5), include.lowest = TRUE, labels = 1:5)
C[, `:=`(add = pred_f - pred_b, srs_vs_bench = y - pred_b)]; C[, q_a := fifth(add)]
BW <- C[, .(n = .N, added = mean(add), srs_vs_bench = mean(srs_vs_bench), se = sd(srs_vs_bench) / sqrt(.N), autism = mean(asd, na.rm = TRUE), n_autism_known = sum(!is.na(asd))), by = .(fifth = q_a)][order(fifth)]
fi <- split(seq_len(nrow(C)), C$family)
bw <- t(replicate(2000, { i <- unlist(fi[sample(length(fi), replace = TRUE)], use.names = FALSE); x <- C[i]
  c(x[q_a == 5, mean(srs_vs_bench)] - x[q_a == 1, mean(srs_vs_bench)], x[q_a == 5, mean(asd, na.rm = TRUE)] - x[q_a == 1, mean(asd, na.rm = TRUE)]) }))
TW <- data.table(quantity = c("SRS beyond genotype x cohort + age, top minus bottom fifth of the added-factor prediction (SD)", "autism diagnosis, same fifths (share)"), model = "added factors",
                 est = c(BW[fifth == 5, srs_vs_bench] - BW[fifth == 1, srs_vs_bench], BW[fifth == 5, autism] - BW[fifth == 1, autism]),
                 lo = apply(bw, 2, quantile, .025, na.rm = TRUE), hi = apply(bw, 2, quantile, .975, na.rm = TRUE), range_is = "family bootstrap 95% CI (2,000 draws)")
cat("\n== within genotype x cohort and age\n"); print(BW, digits = 3); print(TW, digits = 3)
# ---- 4b. CORE model for the figure (Jon, 2026-11-21: "core factors only (and figure legend can mention what's in the full model)").
# Core = benchmark + sex, BMI, ASD PGS, BMI x sex, BMI x ASD PGS: the blocks whose removal lowered out-of-sample R2 in the full-model
# drop analysis above.  NOTE: the core set was chosen from that same cross-validation, so its R2 gain is optimistic; it is reported
# as exploratory, and the full-model numbers above stay the primary record.  Same folds scheme (new draws), drop analysis within the
# core, separation by fifth of the core added-factor prediction, and the data for the sex and ASD PGS scatter panels.
CORE <- c("male", "BMI", "pgs_s", "pgs_s_m", "BMI:male", "BMI:pgs_s")
CDROPS <- list(BMI = c("BMI", "BMI:male", "BMI:pgs_s"), Sex = c("male", "BMI:male"), `ASD PGS` = c("pgs_s", "pgs_s_m", "BMI:pgs_s"), `BMI x sex` = "BMI:male", `BMI x ASD PGS` = "BMI:pgs_s")
CM <- c(list(benchmark = as.formula(BASE), core = fml(CORE)), lapply(CDROPS, function(d) fml(setdiff(CORE, d))))
PC <- lapply(CM, function(z) matrix(NA_real_, nrow(C), R)); r2c <- matrix(NA_real_, R, length(CM), dimnames = list(NULL, names(CM)))
for (r in seq_len(R)) { fo <- setNames(sample(rep(1:10, length.out = length(fams))), fams)[C$family]
  for (kk in 1:10) { tr <- C[fo != kk]; te <- fo == kk; for (nm in names(CM)) PC[[nm]][te, r] <- suppressWarnings(predict(lm(CM[[nm]], tr), C[te])) }
  for (nm in names(CM)) r2c[r, nm] <- r2f(C$y, PC[[nm]][, r]) }
gc_ <- r2c[, "core"] - r2c[, "benchmark"]
R2C <- data.table(model = c("benchmark: genotype x cohort + age + SRS form", "core: + sex, BMI, ASD PGS, BMI x sex, BMI x ASD PGS", "gain (core - benchmark)"),
                  r2 = c(mean(r2c[, "benchmark"]), mean(r2c[, "core"]), mean(gc_)), lo = c(quantile(r2c[, "benchmark"], .025), quantile(r2c[, "core"], .025), quantile(gc_, .025)),
                  hi = c(quantile(r2c[, "benchmark"], .975), quantile(r2c[, "core"], .975), quantile(gc_, .975)), range_is = "2.5-97.5% of 50 CV repeats (not a confidence interval); core set chosen from the full-model drop analysis")
DRC <- rbindlist(lapply(names(CDROPS), function(nm) { l <- r2c[, "core"] - r2c[, nm]
  data.table(factor = nm, removed_terms = paste(CDROPS[[nm]], collapse = ", "), loss_r2 = mean(l), lo = quantile(l, .025), hi = quantile(l, .975), share_of_gain = mean(l) / mean(gc_),
             range_is = "2.5-97.5% of 50 CV repeats (not a confidence interval)") }))[order(-loss_r2)]
C[, `:=`(predc_b = rowMeans(PC$benchmark), predc = rowMeans(PC$core))][, `:=`(addc = predc - predc_b, srs_vs_benchc = y - predc_b)][, q_c := fifth(addc)]
BWC <- C[, .(n = .N, added = mean(addc), srs_vs_bench = mean(srs_vs_benchc), se = sd(srs_vs_benchc) / sqrt(.N), autism = mean(asd, na.rm = TRUE), n_autism_known = sum(!is.na(asd))), by = .(fifth = q_c)][order(fifth)]
bwc <- t(replicate(2000, { i <- unlist(fi[sample(length(fi), replace = TRUE)], use.names = FALSE); x <- C[i]
  c(x[q_c == 5, mean(srs_vs_benchc)] - x[q_c == 1, mean(srs_vs_benchc)], x[q_c == 5, mean(asd, na.rm = TRUE)] - x[q_c == 1, mean(asd, na.rm = TRUE)]) }))
TWC <- data.table(quantity = c("SRS beyond benchmark, top minus bottom fifth of the core added-factor prediction (SD)", "autism diagnosis, same fifths (share)"), model = "core added factors",
                  est = c(BWC[fifth == 5, srs_vs_bench] - BWC[fifth == 1, srs_vs_bench], BWC[fifth == 5, autism] - BWC[fifth == 1, autism]),
                  lo = apply(bwc, 2, quantile, .025, na.rm = TRUE), hi = apply(bwc, 2, quantile, .975, na.rm = TRUE), range_is = "family bootstrap 95% CI (2,000 draws)")
# scatter data: partial residual of SRS for the BMI terms of the core model (in-sample fit), and the BMI slope by sex and by ASD PGS
mc <- lm(fml(CORE), C); bc <- coef(mc); bc <- bc[!is.na(bc)]; vc <- vcov(mc)[names(bc), names(bc)];   # aliased (collinear) benchmark terms dropped
cat("aliased terms in the core model:", names(coef(mc))[is.na(coef(mc))], "\n"); g <- function(t) if (t %in% names(bc)) t else paste(rev(strsplit(t, ":")[[1]]), collapse = ":")
kBMI <- g("BMI"); kS <- g("BMI:male"); kP <- g("BMI:pgs_s")
C[, bmi_part := BMI * (bc[kBMI] + bc[kS] * male + bc[kP] * pgs_s)][, pr_bmi := y - fitted(mc) + bmi_part]
sl <- function(w) { w <- w[names(bc)]; w[is.na(w)] <- 0; est <- sum(w * bc); se <- sqrt(drop(t(w) %*% vc %*% w)); c(est = est, se = se, p = 2 * pnorm(-abs(est / se))) }
z <- setNames(rep(0, length(bc)), names(bc))
SL <- rbind(data.table(panel = "sex", group = "girls", t(sl(replace(z, kBMI, 1)))), data.table(panel = "sex", group = "boys", t(sl(replace(z, c(kBMI, kS), 1)))),
            # ASD PGS slopes at the sample's share of boys (the sex slopes are at the mean ASD PGS, 0)
            data.table(panel = "ASD PGS", group = "low ASD PGS (-1 SD)", t(sl(replace(replace(replace(z, kBMI, 1), kP, -1), kS, mean(C$male))))),
            data.table(panel = "ASD PGS", group = "high ASD PGS (+1 SD)", t(sl(replace(replace(replace(z, kBMI, 1), kP, 1), kS, mean(C$male))))))
SL[, n := c(sum(C$male == 0), sum(C$male == 1), sum(C$pgs_s_m == 0), sum(C$pgs_s_m == 0))]
SL <- rbind(SL, data.table(panel = c("sex", "ASD PGS"), group = "interaction", est = c(bc[kS], bc[kP]), se = sqrt(c(vc[kS, kS], vc[kP, kP])), p = 2 * pnorm(-abs(c(bc[kS], bc[kP]) / sqrt(c(vc[kS, kS], vc[kP, kP])))), n = c(nrow(C), sum(C$pgs_s_m == 0))))
# overall accuracy and P (Jon: "is there an accuracy and p-value for the overall prediction?")
#   (a) permutation P of the out-of-sample R2 gain: the added factors (sex, BMI, ASD PGS and its indicator, as one block per carrier)
#       are shuffled among carriers of the same genotype x cohort cell, which keeps the benchmark and the cell composition and breaks
#       only the link between the added factors and SRS.  Statistic = mean gain over 5 fixed family-grouped 10-fold splits, the same
#       splits for the observed data and every permutation; 1,000 permutations; P = (1 + #null >= observed) / 1001.
#       # null model: the added factors permuted within genotype x cohort.
#   (b) accuracy for autism diagnosis: AUC of the out-of-sample core added-factor prediction (addc) for a recorded diagnosis, with the
#       Mann-Whitney P; and the out-of-sample correlation of addc with SRS beyond the benchmark.
FOLD5 <- lapply(1:5, function(i) setNames(sample(rep(1:10, length.out = length(fams))), fams)[C$family])
cvgain <- function(D) mean(sapply(FOLD5, function(fo) { pb <- pc <- numeric(nrow(D)); for (kk in 1:10) { tr <- D[fo != kk]; te <- fo == kk
  pb[te] <- suppressWarnings(predict(lm(CM$benchmark, tr), D[te])); pc[te] <- suppressWarnings(predict(lm(CM$core, tr), D[te])) }; r2f(D$y, pc) - r2f(D$y, pb) }))
obs_gain <- cvgain(C); NPERM <- 1000; blk <- c("male", "BMI", "pgs_s", "pgs_s_m")
null_gain <- replicate(NPERM, { D <- copy(C); D[, (blk) := .SD[sample(.N)], by = cg, .SDcols = blk]; cvgain(D) })
pperm <- (1 + sum(null_gain >= obs_gain)) / (NPERM + 1)
ka <- C[!is.na(asd)]; auc <- ka[, { r <- rank(addc); n1 <- sum(asd == 1); n0 <- sum(asd == 0); (sum(r[asd == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }]
pw <- wilcox.test(addc ~ asd, data = ka)$p.value; rr <- cor.test(C$addc, C$srs_vs_benchc)
ACC <- data.table(quantity = c("out-of-sample R2 gain (mean of 5 fixed splits)", "permutation P of the R2 gain (1,000 permutations within genotype x cohort)",
                               "null R2 gain, mean (95% range)", "AUC of the core added-factor prediction for autism diagnosis", "Mann-Whitney P (AUC)",
                               "correlation of the added-factor prediction with SRS beyond the benchmark (out of sample)", "P (correlation)"),
                  est = c(obs_gain, pperm, mean(null_gain), auc, pw, unname(rr$estimate), rr$p.value), lo = c(NA, NA, quantile(null_gain, .025), NA, NA, rr$conf.int[1], NA),
                  hi = c(NA, NA, quantile(null_gain, .975), NA, NA, rr$conf.int[2], NA), n = c(nrow(C), nrow(C), nrow(C), nrow(ka), nrow(ka), nrow(C), nrow(C)))
ACC[grepl("^permutation P", quantity), note := "P < 0.001 (floor of 1,000 permutations); does not account for choosing the core set from the full-model drop analysis"]
fwrite(ACC, file.path(MMOUT, "srs_core_accuracy.csv"))
# core model per CNV (rule: one analysis per CNV, then pooled with a heterogeneity test); Cochran's Q across the four CNVs
LABC <- c(male = "Sex", BMI = "BMI", pgs_s = "ASD PGS", `BMI:male` = "BMI x sex (male)", `BMI:pgs_s` = "BMI x ASD PGS")
PERC <- rbindlist(lapply(LOCI, function(G) { x <- C[locus == G]; base <- paste(if (uniqueN(x$cohort) > 1) "y ~ cohort + ns(age, 3)" else "y ~ ns(age, 3)", if (uniqueN(x$form2) > 1) "+ form2 + form2:age" else "")
  m <- lm(fml(CORE, base), x); rbindlist(lapply(names(LABC), function(t) cbind(locus = G, n = nrow(x), n_with_ASD_PGS = sum(x$pgs_s_m == 0), cohorts = paste(sort(unique(x$cohort)), collapse = "+"), factor = LABC[[t]], cf1(m, t)))) }))
hetc <- PERC[, { ok <- !is.na(est) & !is.na(se); w <- 1 / se[ok]^2; b <- sum(w * est[ok]) / sum(w); q <- sum(w * (est[ok] - b)^2); .(Q = q, df = sum(ok) - 1L, p_het = pchisq(q, sum(ok) - 1, lower.tail = FALSE)) }, by = factor]
PERC <- merge(PERC, hetc, by = "factor", sort = FALSE)
fwrite(PERC, file.path(MMOUT, "srs_core_per_cnv.csv")); cat("\n== CORE model per CNV\n"); print(PERC, digits = 3)
cat("\n== CORE model\n"); print(ACC, digits = 3); print(R2C, digits = 3); print(DRC, digits = 3); print(BWC, digits = 3); print(TWC, digits = 3); print(SL, digits = 3)
fwrite(R2C, file.path(MMOUT, "srs_core_r2.csv")); fwrite(DRC, file.path(MMOUT, "srs_core_drop_analysis.csv")); fwrite(BWC, file.path(MMOUT, "srs_core_within_cell_by_fifth.csv"))
fwrite(TWC, file.path(MMOUT, "srs_core_separation.csv")); fwrite(SL, file.path(MMOUT, "srs_core_bmi_slopes_by_sex_and_pgs.csv"))
saveRDS(C[, .(locus, cohort, male, pgs_s, pgs_s_m, BMI, pr_bmi, q_c, srs_vs_benchc, asd)], file.path(MMPRIV, "srs_core_scatter.rds"))
# ---- 5. sensitivity (G2MH only; medication data exist only in G2MH, so they are not in the pooled model): psychotropic medication.
# Classes exactly as in the paper (27_medication_sensitivity.R; its RX list and weight-gain / weight-loss rules are read from that
# script): current weight-gain medication (antipsychotics, valproate, lithium, mirtazapine, gabapentinoids) and current weight-loss
# medication (stimulants, atomoxetine, topiramate, metformin, GLP-1 agonists).  Medication status is known for carriers with any row in
# the G2MH medication release (an empty row = asked, none).  Same full model within G2MH, with and without the two flags.
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
m <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_medications_data_release_022026.csv"), colClasses = "character")
m[, nm := tolower(trimws(meds_name))][, current := meds_current == "1"]
s27 <- readLines(file.path(HERE, "27_medication_sensitivity.R"), warn = FALSE); i0 <- grep("^RX <- list\\(", s27); i1 <- i0 + which(grepl("anxiolytic", s27[i0:length(s27)]))[1] - 1
eval(parse(text = s27[i0:i1])); for (k in names(RX)) m[, (k) := grepl(RX[[k]], nm)]
eval(parse(text = grep("^m\\[, weight_gain :=", s27, value = TRUE))); eval(parse(text = grep("^m\\[, weight_loss :=", s27, value = TRUE)))
med <- m[, .(wg = as.numeric(any(current & weight_gain)), wl = as.numeric(any(current & weight_loss)), ap = as.numeric(any(current & antipsychotic))), by = .(sampleID = rarecnv_id)]
med <- merge(med, im[cohort == "G2MH" & !is.na(GUID), .(id = GUID, sampleID)], by = "sampleID")[, .(wg = max(wg), wl = max(wl), ap = max(ap)), by = id]
G2 <- merge(C[cohort == "G2MH"], med, by = "id")
MS <- if (nrow(G2) >= 100) { b0 <- lm(fml(ADD), G2); b1 <- lm(fml(c(ADD, "wg", "wl")), G2)
  tt <- c(names(LAB), "wg", "wl"); rbindlist(lapply(tt, function(t) cbind(term = ifelse(t %in% names(LAB), LAB[t], c(wg = "current weight-gain medication", wl = "current weight-loss medication")[t]),
    without_medication = cf1(b0, t)[, .(est0 = est, p0 = p)], with_medication = cf1(b1, t)[, .(est1 = est, p1 = p)])))[, `:=`(n = nrow(G2), on_weight_gain = sum(G2$wg), on_weight_loss = sum(G2$wl), on_antipsychotic = sum(G2$ap))] } else NULL
cat("\n== medication sensitivity, G2MH carriers with medication status\n"); print(MS, digits = 3)
# ---- outputs
F <- MMOUT
if (!is.null(MS)) fwrite(MS, file.path(F, "srs_multimodal_medication_sensitivity_g2mh.csv"))
for (f in c("srs_multimodal_by_fifth.csv")) if (file.exists(file.path(F, f))) file.remove(file.path(F, f))   # benchmark-fifth table of the first version, no longer produced
fwrite(CO, file.path(F, "srs_multimodal_coefficients.csv")); fwrite(PER, file.path(F, "srs_multimodal_per_cnv.csv")); fwrite(DR, file.path(F, "srs_multimodal_drop_analysis.csv"))
fwrite(BW, file.path(F, "srs_multimodal_within_cell_by_fifth.csv")); fwrite(smp, file.path(F, "srs_multimodal_sample.csv"))
fwrite(rbind(R2[, .(quantity = "out-of-sample R2", model, est = r2, lo, hi, range_is)], TW), file.path(F, "srs_multimodal_prediction.csv"))
saveRDS(list(C = C[, .(cohort, locus, sex, age, y, pred_b, pred_f, add, srs_vs_bench, q_a, asd)], built = Sys.time()), file.path(MMPRIV, "srs_multimodal_panels.rds"))
sink(file.path(F, "SRS_multimodal_numbers.txt")); cat("SRS model in CNV carriers (growth, sex, genetic background). Built", format(Sys.time()), "by 92_v2_srs_multimodal.R. Aggregate numbers only.\n")
cat("SSC carriers dropped:", n_ssc, "\n\n== sample\n"); print(smp); cat("\n== pooled coefficients\n"); print(CO, digits = 3); cat("\n== per CNV and heterogeneity\n"); print(PER, digits = 3)
cat("\n== out-of-sample R2\n"); print(R2, digits = 3); cat("\n== drop analysis\n"); print(DR, digits = 3); cat("\n== within genotype x cohort and age\n"); print(BW, digits = 3); print(TW, digits = 3); cat("\n== medication sensitivity (G2MH)\n"); print(MS, digits = 3); sink()
cat("done\n")
