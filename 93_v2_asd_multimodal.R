# 93_v2_asd_multimodal.R
# Ported to the BMI modifier paper (JS / Claude 2026-11-21) from COMET_ASD_ARPAH/analysis/A11_asd_multimodal.R (ARPA-H session, same day), unchanged
# except the output folders: aggregates -> run_v2_2026-11-20/multimodal/, participant-level -> DATA/v2_autism_cache/ (not in git).
# The model rules are the paper's: adults Box-Cox, children WHO / UK90; one intercept per cohort within CNV; per CNV + pooled with
# heterogeneity; SSC carriers dropped (ascertained probands); family-grouped cross-validation and family bootstrap.
# COMET-ASD / BMI paper supplement (2026-11-21; Jon, relayed by the Sacks BMI-paper session: "Can we do the same thing with ASD
# diagnosis (as a supplemental)").  The model of A10_srs_multimodal.R with autism diagnosis as the outcome.
# Sample: the paper's autism-diagnosis sample (cnv-growth-bmi 85_v2_asd_by_cnv.R via 83_v2_autism_common.R; cache object XA, BMI rows:
#   one visit per carrier with BMI Z and a recorded diagnosis): G2MH, Cardiff, Simons VIP, CHU Sainte-Justine, CHOP (22q11.2 A-D).
#   SSC carriers are excluded (ascertained as autistic probands).  Searchlight is not in this sample.  Cohort x CNV cells without both
#   outcomes are dropped (rule of 85); the number removed is reported.
# Outcome: autism diagnosis (0/1).  Logistic regression.
# Predictors (all at the diagnosis-sample visit): benchmark = genotype x cohort (one intercept per cell) + age (natural spline, 3 df);
#   core = benchmark + sex, BMI Z (paper's, as in XA), ASD PGS (with missing indicator), BMI x sex, BMI x ASD PGS;
#   full = core + height Z, head-circumference Z (each with a missing indicator), childhood-obesity / adult-BMI / height PGS (shared
#   missing indicator), BMI x CNV.
#   Height and head circumference: children (height < 19 y, HC < 18 y) WHO 2006/2007 and UK90 Z; adults Box-Cox against the non-carrier
#   adults of the clinical cohorts with the paper's exponents (height 0.65 men / 0.60 women; HC 1), as in A10.  Raw values from
#   analysis_table_chop_unaffctl_2026-11-20.rds at the same visit (id and age); Simons VIP from its v10 growth_measurements.csv, visit
#   within 0.1 y of the diagnosis-sample age.
#   ASD PGS: the April 2026 autism score of the paper (cache S), linked through the cache id map (GUID or sample id); standardised in
#   these carriers.  Other scores: as assembled in 86_v2_cnv_bmi_pgs_stats.R (G2MH: pgs49 cache g2all; Simons VIP: combined_scores.tsv
#   residualised on 10 global PCs), standardised in these carriers.
# Analyses (as A10): coefficients pooled (core and full) and per CNV with Cochran's Q; family-grouped 10-fold cross-validation, 50
#   repeats: out-of-sample AUC and log loss of benchmark, core and full; drop-one importance (loss of AUC and increase in log loss when
#   a block is removed); share with a diagnosis by fifth of the added-factor prediction (core) within genotype x cohort and age, with a
#   family bootstrap CI; genotyped-subset check (carriers with an ASD PGS only, no PGS missing indicator).
#   Range over CV repeats = 2.5-97.5% of 50 repeats, not a confidence interval.
#   # null model: none.  Wald tests; Cochran's Q; out-of-sample comparison with the benchmark.  Exploratory.
# Caveats: CHOP records do not say whether a carrier without a recorded diagnosis was evaluated; how diagnosis is recorded differs by
#   cohort (absorbed by the cohort x CNV intercepts, not removed).
# Outputs: aggregate CSVs figures/asd_multimodal_*.csv and figures/ASD_multimodal_numbers.txt.
suppressPackageStartupMessages({ library(data.table); library(splines) })
REPO <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/COMET_ASD_ARPAH"
MMOUT <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/run_v2_2026-11-20/multimodal"; dir.create(MMOUT, showWarnings = FALSE)
MMPRIV <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/v2_autism_cache"
F <- MMOUT
DBX <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; DATA <- file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data")
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
SV <- file.path(DBX, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0"); AN <- file.path(DBX, "G2MH/NetworkPaper/AnalysisByAnjali/data")
source(file.path(HERE, "00_scoring.R")); set.seed(20261121); options(width = 220)
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
e <- new.env(); load(file.path(DATA, "v2_autism_cache", "base_2026-11-20.RData"), envir = e)
X <- as.data.table(e$XA)[pheno == "BMI" & locus %in% LOCI & !is.na(asd), .(id, cohort, locus, sex = as.character(sex), age, BMI, asd)]
n0 <- X[, .N, by = cohort]; X <- X[cohort != "SSC"]; n_ssc <- n0[cohort == "SSC", N]
cells <- X[, .(n = .N, n_asd = sum(asd)), by = .(cohort, locus)]; bad <- cells[n_asd == 0 | n_asd == n]
X <- X[!paste(cohort, locus) %in% bad[, paste(cohort, locus)]]
cat("SSC carriers excluded:", n_ssc, "| cells without both outcomes dropped:", nrow(bad), "(", sum(bad$n), "carriers:", paste(bad[, paste(cohort, locus)], collapse = "; "), ")\n")
# ---- growth measures at the same visit
T0 <- as.data.table(readRDS(file.path(DATA, "analysis_table_chop_unaffctl_2026-11-20.rds")))
REF <- T0[locus == "No CNV" & latest %in% TRUE, .(age, sex = as.character(sex), h = HEIGHT_raw, hc = HEADCIRCUM_raw)]
LAMBDA <- list(height = c(male = 0.65, female = 0.60), hc = c(male = 1, female = 1))
zgrow <- function(v, age, sex, item) { z <- if (item == "height") score_height(v, age, sex) else score_hc(v, age, sex); amin <- if (item == "height") 19 else 18
  for (s in c("male", "female")) { i <- which(sex == s & !is.na(age) & age >= amin & !is.na(v)); if (!length(i)) next
    rv <- if (item == "height") REF$h else REF$hc; r <- which(REF$sex == s & REF$age >= amin & !is.na(rv))
    zb <- score_adult_bc(c(rv[r], v[i]), c(REF$age[r], age[i]), rep(s, length(r) + length(i)), rep("Clinical", length(r) + length(i)), c(rep(TRUE, length(r)), rep(FALSE, length(i))),
                         min_age = amin, lambda_fixed = LAMBDA[[item]][[s]]); z[i] <- zb[length(r) + seq_along(i)] }
  z }
g1 <- merge(X[cohort != "SVIP", .(id, age)], T0[, .(id, age, h = HEIGHT_raw, hc = HEADCIRCUM_raw, family = as.character(family))], by = c("id", "age"))
svr <- fread(file.path(SV, "growth_measurements.csv"))[, .(sid = individual, tage = eval_age_months / 12, h = growth_measurements_raw.height, hc = growth_measurements_raw.head_circumference, family = as.character(family))]
g2 <- merge(X[cohort == "SVIP", .(id, sid = sub("^SVIP:", "", id), age)], svr, by = "sid", allow.cartesian = TRUE)[abs(tage - age) <= 0.1][order(id, abs(tage - age))][, .SD[1], by = id]
G <- rbind(g1, g2[, .(id, age, h, hc, family)], fill = TRUE)
X <- merge(X, G[, .(id, h, hc, family)], by = "id", all.x = TRUE)
X[, `:=`(hgt = zgrow(h, age, sex, "height"), hcz = zgrow(hc, age, sex, "hc"))]
X[!is.na(hgt) & abs(hgt) > 5, hgt := NA_real_]; X[!is.na(hcz) & abs(hcz) > 5, hcz := NA_real_]
X[is.na(family), family := id]
# ---- polygenic scores
im <- as.data.table(e$idmap); S <- as.data.table(e$S)
k <- merge(X[, .(id)], im[, .(id, sampleID, GUID)], by = "id", all.x = TRUE)
k[, key := fifelse(!is.na(GUID) & GUID %in% S$id, GUID, fifelse(!is.na(sampleID) & sampleID %in% S$id, sampleID, NA_character_))]
k <- unique(k[!is.na(key), .(id, key)], by = "id"); X <- merge(X, merge(k, unique(S[, .(key = id, pgs)], by = "key"), by = "key")[, .(id, pgs, key)], by = "id", all.x = TRUE)
e49 <- new.env(); load(file.path(DATA, "v2_autism_cache", "pgs49_2026-11-20.RData"), envir = e49)
rzp <- function(v, P) as.numeric(scale(residuals(lm(as.formula(paste("v ~", paste(names(P), collapse = " + "))), data = cbind(P, v = v)))))
cs <- fread(file.path(AN, "searchlight/pgs/combined_scores.tsv"), select = c("IID", "childhood_obesity", "bmi", "height"))
pcs <- fread(file.path(AN, "searchlight/pgs/analysis_dataset.tsv"), select = c("IID", paste0("GLOBAL_PC", 1:10)))
meta <- unique(fread(file.path(DBX, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/Searchlight.WGS2.2023_05.sample_metadata.tsv"))[, .(searchlight_id, sfari_iid)])
sv <- merge(merge(cs, pcs, by = "IID"), meta, by.x = "IID", by.y = "searchlight_id"); PCm <- sv[, paste0("GLOBAL_PC", 1:10), with = FALSE]
sv[, `:=`(pgs_cobes = rzp(childhood_obesity, PCm), pgs_bmi = rzp(bmi, PCm), pgs_ht = rzp(height, PCm))][, key := sub("-", ".", sfari_iid, fixed = TRUE)]
ADS <- unique(rbind(as.data.table(e49$g2all)[, .(key = GUID, pgs_cobes, pgs_bmi, pgs_ht)], sv[, .(key, pgs_cobes, pgs_bmi, pgs_ht)]), by = "key")
X <- merge(X, ADS, by = "key", all.x = TRUE)
for (v in c("pgs", "pgs_cobes", "pgs_bmi", "pgs_ht")) X[, (v) := (get(v) - mean(get(v), na.rm = TRUE)) / sd(get(v), na.rm = TRUE)]
X[, opgs_m := as.numeric(is.na(pgs_cobes) | is.na(pgs_bmi) | is.na(pgs_ht))]; for (v in c("pgs_cobes", "pgs_bmi", "pgs_ht")) X[opgs_m == 1, (v) := 0]
setnames(X, "pgs", "pgs_s"); for (v in c("hgt", "hcz", "pgs_s")) { X[, paste0(v, "_m") := as.numeric(is.na(get(v)))]; X[is.na(get(v)), (v) := 0] }
X[, `:=`(cg = factor(paste(cohort, locus)), male = as.numeric(sex == "male"), locus = factor(locus, LOCI))]
stopifnot(!anyDuplicated(X$id))
smp <- X[, .(carriers = .N, autism = sum(asd), adults_19plus = sum(age >= 19), with_height = sum(hgt_m == 0), with_head_circumference = sum(hcz_m == 0), with_ASD_PGS = sum(pgs_s_m == 0),
             with_other_PGS = sum(opgs_m == 0), families = uniqueN(family)), by = .(cohort, locus)][order(locus, cohort)]
smp <- rbind(smp, X[, .(cohort = "All", locus = "All", carriers = .N, autism = sum(asd), adults_19plus = sum(age >= 19), with_height = sum(hgt_m == 0), with_head_circumference = sum(hcz_m == 0),
             with_ASD_PGS = sum(pgs_s_m == 0), with_other_PGS = sum(opgs_m == 0), families = uniqueN(family))])
cat("== sample\n"); print(smp)
# ---- models
BASE <- "asd ~ cg + ns(age, 3)"
CORE <- c("male", "BMI", "pgs_s", "pgs_s_m", "BMI:male", "BMI:pgs_s")
FULL <- c(CORE, "hgt", "hgt_m", "hcz", "hcz_m", "pgs_cobes", "pgs_bmi", "pgs_ht", "opgs_m", "BMI:locus")
fml <- function(add, base = BASE) as.formula(paste(base, "+", paste(add, collapse = " + ")))
fit <- function(f, d) suppressWarnings(glm(f, binomial, d))
cf1 <- function(m, term) { s <- summary(m)$coefficients; rv <- paste(rev(strsplit(term, ":")[[1]]), collapse = ":"); if (!term %in% rownames(s) && rv %in% rownames(s)) term <- rv
  if (!term %in% rownames(s)) return(data.table(logOR = NA_real_, se = NA_real_, p = NA_real_)); data.table(logOR = s[term, 1], se = s[term, 2], p = s[term, 4]) }
LAB <- c(male = "Sex (male)", BMI = "BMI", pgs_s = "ASD PGS", `BMI:male` = "BMI x sex (male)", `BMI:pgs_s` = "BMI x ASD PGS", hgt = "Height", hcz = "Head circumference",
         pgs_cobes = "Childhood obesity PGS", pgs_bmi = "Adult BMI PGS", pgs_ht = "Height PGS")
coefs <- function(add, label, d = X) { m <- fit(fml(add), d); rbindlist(lapply(intersect(names(LAB), sub(":.*", "", add) |> c(add)), function(t) cbind(model = label, factor = LAB[[t]], cf1(m, t)))) }
CO <- rbind(coefs(CORE, "core"), coefs(FULL, "full"))[, `:=`(OR = exp(logOR), lo = exp(logOR - 1.96 * se), hi = exp(logOR + 1.96 * se))]
cat("\n== pooled coefficients (log odds per SD; sex = male vs female)\n"); print(CO, digits = 3)
PER <- rbindlist(lapply(LOCI, function(L) { x <- X[locus == L]; base <- if (uniqueN(x$cohort) > 1) "asd ~ cohort + ns(age, 3)" else "asd ~ ns(age, 3)"
  m <- fit(fml(CORE, base), x); rbindlist(lapply(c("male", "BMI", "pgs_s", "BMI:male", "BMI:pgs_s"), function(t) cbind(locus = L, n = nrow(x), n_asd = sum(x$asd), n_ASD_PGS = sum(x$pgs_s_m == 0),
    cohorts = paste(sort(unique(x$cohort)), collapse = "+"), factor = LAB[[t]], cf1(m, t)))) }))
het <- PER[, { ok <- !is.na(logOR) & !is.na(se) & se < 10; w <- 1 / se[ok]^2; b <- sum(w * logOR[ok]) / sum(w); q <- sum(w * (logOR[ok] - b)^2); .(Q = q, df = sum(ok) - 1L, p_het = pchisq(q, sum(ok) - 1, lower.tail = FALSE)) }, by = factor]
PER <- merge(PER, het, by = "factor", sort = FALSE); cat("\n== core model per CNV\n"); print(PER, digits = 3)
# ---- cross-validation: AUC, log loss, drop-one
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ll <- function(y, p) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -mean(y * log(p) + (1 - y) * log(1 - p)) }
CDROPS <- list(BMI = c("BMI", "BMI:male", "BMI:pgs_s"), Sex = c("male", "BMI:male"), `ASD PGS` = c("pgs_s", "pgs_s_m", "BMI:pgs_s"), `BMI x sex` = "BMI:male", `BMI x ASD PGS` = "BMI:pgs_s")
FDROPS <- list(BMI = c("BMI", "BMI:male", "BMI:pgs_s", "BMI:locus"), Sex = c("male", "BMI:male"), Height = c("hgt", "hgt_m"), `Head circumference` = c("hcz", "hcz_m"), `ASD PGS` = c("pgs_s", "pgs_s_m", "BMI:pgs_s"),
               `Adiposity and height PGS` = c("pgs_cobes", "pgs_bmi", "pgs_ht", "opgs_m"), `BMI x CNV` = "BMI:locus", `BMI x sex` = "BMI:male", `BMI x ASD PGS` = "BMI:pgs_s")
MODELS <- c(list(benchmark = as.formula(BASE), core = fml(CORE), full = fml(FULL)), setNames(lapply(CDROPS, function(d) fml(setdiff(CORE, d))), paste("core -", names(CDROPS))),
            setNames(lapply(FDROPS, function(d) fml(setdiff(FULL, d))), paste("full -", names(FDROPS))))
cvrun <- function(D, MOD, R = 50) { fams <- unique(D$family); P <- lapply(MOD, function(z) matrix(NA_real_, nrow(D), R)); A <- L <- matrix(NA_real_, R, length(MOD), dimnames = list(NULL, names(MOD)))
  for (r in seq_len(R)) { fo <- setNames(sample(rep(1:10, length.out = length(fams))), fams)[D$family]
    for (kk in 1:10) { tr <- D[fo != kk]; te <- fo == kk; for (nm in names(MOD)) P[[nm]][te, r] <- suppressWarnings(predict(fit(MOD[[nm]], tr), D[te], type = "response")) }
    for (nm in names(MOD)) { A[r, nm] <- auc(D$asd, P[[nm]][, r]); L[r, nm] <- ll(D$asd, P[[nm]][, r]) } }
  list(P = P, A = A, L = L) }
CV <- cvrun(X, MODELS)
rng <- function(v) c(mean = mean(v), lo = unname(quantile(v, .025)), hi = unname(quantile(v, .975)))
PERF <- rbindlist(lapply(c("benchmark", "core", "full"), function(nm) data.table(model = nm, metric = c("AUC", "log loss"), rbind(rng(CV$A[, nm]), rng(CV$L[, nm])))))
PERF <- rbind(PERF, rbindlist(lapply(c("core", "full"), function(nm) data.table(model = paste(nm, "- benchmark"), metric = c("AUC gain", "log-loss reduction"),
              rbind(rng(CV$A[, nm] - CV$A[, "benchmark"]), rng(CV$L[, "benchmark"] - CV$L[, nm]))))))[, range_is := "2.5-97.5% of 50 CV repeats (not a confidence interval)"]
DR <- rbindlist(lapply(names(MODELS)[grepl(" - ", names(MODELS))], function(nm) { base <- sub(" - .*", "", nm)
  data.table(model = base, factor = sub(".* - ", "", nm), t(rng(CV$A[, base] - CV$A[, nm])), ll_increase = mean(CV$L[, nm] - CV$L[, base])) }))[order(model, -mean)]
setnames(DR, c("mean", "lo", "hi"), c("auc_loss", "auc_loss_lo", "auc_loss_hi")); DR[, range_is := "2.5-97.5% of 50 CV repeats (not a confidence interval)"]
cat("\n== out-of-sample performance\n"); print(PERF, digits = 3); cat("\n== drop-one importance\n"); print(DR, digits = 3)
# ---- separation within genotype x cohort and age (core)
fifth <- function(x) cut(rank(x, ties.method = "first"), quantile(seq_along(x), 0:5 / 5), include.lowest = TRUE, labels = 1:5)
lg <- function(p) log(p / (1 - p))
X[, add := rowMeans(lg(pmin(pmax(CV$P$core, 1e-6), 1 - 1e-6))) - rowMeans(lg(pmin(pmax(CV$P$benchmark, 1e-6), 1 - 1e-6)))][, q := fifth(add)]
BW <- X[, .(n = .N, autism = mean(asd), expected = mean(rowMeans(CV$P$benchmark)[.I])), by = .(fifth = q)][order(fifth)]
BW <- X[, .(n = .N, autism = mean(asd)), by = .(fifth = q)][order(fifth)]
fi <- split(seq_len(nrow(X)), X$family)
bt <- replicate(2000, { i <- unlist(fi[sample(length(fi), replace = TRUE)], use.names = FALSE); x <- X[i]; x[q == 5, mean(asd)] - x[q == 1, mean(asd)] })
SEP <- data.table(quantity = "autism diagnosis, top minus bottom fifth of the core added-factor prediction (share)", est = BW[fifth == 5, autism] - BW[fifth == 1, autism],
                  lo = quantile(bt, .025), hi = quantile(bt, .975), range_is = "family bootstrap 95% CI (2,000 draws)")
cat("\n== share with a diagnosis by fifth of the core added-factor prediction\n"); print(BW, digits = 3); print(SEP, digits = 3)
# ---- genotyped-subset check (carriers with an ASD PGS; no PGS missing indicator)
XG <- X[pgs_s_m == 0][, cg := droplevels(cg)]
cg_ok <- XG[, .(n = .N, a = sum(asd)), by = cg][a > 0 & a < n, cg]; XG <- XG[cg %in% cg_ok][, cg := droplevels(cg)]
CORE_G <- setdiff(CORE, "pgs_s_m"); FULL_G <- setdiff(FULL, c("pgs_s_m", "opgs_m"))
COG <- rbind(coefs(CORE_G, "core, genotyped subset", XG), coefs(FULL_G[!(FULL_G %in% c("pgs_cobes", "pgs_bmi", "pgs_ht")) | XG[, all(opgs_m == 0)]], "full, genotyped subset", XG))[, OR := exp(logOR)]
CVG <- cvrun(XG, list(benchmark = as.formula(BASE), core = fml(CORE_G)), R = 50)
PG <- data.table(model = c("benchmark", "core", "core - benchmark"), metric = c("AUC", "AUC", "AUC gain"), rbind(rng(CVG$A[, "benchmark"]), rng(CVG$A[, "core"]), rng(CVG$A[, "core"] - CVG$A[, "benchmark"])),
                 n = nrow(XG), n_asd = sum(XG$asd))
cat("\n== genotyped subset (n", nrow(XG), ", autism", sum(XG$asd), ")\n"); print(COG, digits = 3); print(PG, digits = 3)
# ---- outputs
fwrite(smp, file.path(F, "asd_multimodal_sample.csv")); fwrite(CO, file.path(F, "asd_multimodal_coefficients.csv")); fwrite(PER, file.path(F, "asd_multimodal_core_per_cnv.csv"))
fwrite(PERF, file.path(F, "asd_multimodal_performance.csv")); fwrite(DR, file.path(F, "asd_multimodal_drop_analysis.csv")); fwrite(BW, file.path(F, "asd_multimodal_by_fifth.csv"))
fwrite(SEP, file.path(F, "asd_multimodal_separation.csv")); fwrite(rbind(COG, fill = TRUE), file.path(F, "asd_multimodal_genotyped_subset_coefficients.csv")); fwrite(PG, file.path(F, "asd_multimodal_genotyped_subset_performance.csv"))
sink(file.path(F, "ASD_multimodal_numbers.txt")); cat("Autism-diagnosis model in CNV carriers (supplement). Built", format(Sys.time()), "by 93_v2_asd_multimodal.R. Aggregate numbers only.\n")
cat("SSC carriers excluded:", n_ssc, "| cells without both outcomes dropped:", nrow(bad), "(", sum(bad$n), "carriers )\n\n== sample\n"); print(smp)
cat("\n== coefficients\n"); print(CO, digits = 3); cat("\n== core per CNV\n"); print(PER, digits = 3); cat("\n== performance\n"); print(PERF, digits = 3)
cat("\n== drop-one\n"); print(DR, digits = 3); cat("\n== by fifth\n"); print(BW, digits = 3); print(SEP, digits = 3); cat("\n== genotyped subset\n"); print(COG, digits = 3); print(PG, digits = 3); sink()
cat("done\n")
