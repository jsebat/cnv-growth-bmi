# Sensitivity analysis for the main figure (JS / Claude 2026-10-26): are the BMI-driven effects (CNV x BMI and
# CNV x age x BMI on height and head circumference vs unaffected controls; BMI slopes on SRS / SCQ vs controls; BMI slopes
# on Penn CNB) driven by sex or by polygenic background? The main-figure tests are rerun within strata of
#   1. sex (male / female; carriers and controls both stratified),
#   2. childhood obesity PGS (Bradfield 2019 trans-ancestral; SBayesRC weights from Anjali), below / above median,
#   3. adult BMI PGS (Yengo 2018; same pipeline), below / above median,
# with carriers and SSC controls both stratified (SSC genotypes scored with the same weights: pgs_obesity_ssc on Expanse;
# G2MH WGS carriers from Anjali's pgs_pca.csv). PGS are residualised on 10 ancestry PCs within cohort and z-scored within
# cohort; the median split is within cohort. SVIP / Searchlight / Cardiff carriers have no scores in this pipeline and are
# excluded from PGS strata (they stay in the sex strata). Within each stratum: within-group BMI slope and age x BMI,
# CNV x BMI shift and CNV x age x BMI vs the controls of the same stratum (models of 23/24/28). Heterogeneity between the two
# strata within carriers: BMI x stratum (age-adjusted; locus covariate when pooled). For PGS, also the continuous
# BMI x PGS interaction within carriers and the PGS main effect.
# Inputs: 23 (A, gmc, ctl_g), 24 (X), Anjali pgs_pca.csv, resources/ssc_pgs_scores.tsv (built by 30a_build_ssc_pgs.R),
#         all_cohorts_adsi_pgs.tsv (SSC PCs), SSC IDs Sanders.csv (SS id -> SSC individual).
# Outputs: supp_table_S4_sex_pgs_strata.csv, sheet "S4 sex and PGS strata" in supp_tables_bmi_effects.xlsx, md section,
#          pgs_strata_sample.csv (who is in which stratum, aggregate counts only in the md).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))
source(file.path(HERE, "24_cognition_bmi_slopes.R"))
suppressPackageStartupMessages({ library(mgcv); library(openxlsx) })
t0 <- Sys.time()
AD <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data")

# ---------------------------------------------------------------- PGS tables
# G2MH (WGS): raw plink2 scores (Anjali, scores_wgs; IID = "<batch>-<GUID>") residualised on the 10 projected ancestry PCs of
# the WGS cohort (pca_results/naive/cohort.sscore) and z-scored, i.e. the same treatment as the SSC scores in 30a. (Anjali's
# pgs_pca.csv holds the same scores already residualised but covers 1500 of the 1762 scored samples.)
rdw <- function(nm) { x <- fread(file.path(AD, "pgs_scores/scores_wgs", paste0(nm, ".sscore"))); setnames(x, "#IID", "IID"); x[, .(GUID = sub("^[0-9]+-", "", IID), score = SCORE1_AVG)] }
# height: Anjali's GIANT 2022 score; autism: April 2026 ASD GWAS, full weights, scored on the G2MH WGS cohort on Expanse
# (pgs_obesity_ssc/scores/g2mh_wgs.asd2026_full.sscore; G2MH is not a GWAS cohort), copied to Jon/Claude data folder.
rda <- function(p) { x <- fread(p); setnames(x, "#IID", "IID"); x[, .(GUID = sub("^[0-9]+-", "", IID), score = SCORE1_AVG)] }
g2 <- Reduce(function(a, b) merge(a, b, by = "GUID"), list(setnames(rdw("childhood_obesity"), "score", "cobes"), setnames(rdw("bmi"), "score", "bmi"), setnames(rdw("bmi_childhood"), "score", "cbmi"),
                                                          setnames(rdw("height"), "score", "ht"), setnames(rda(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/g2mh_scores_raw/g2mh_wgs.asd2026_full.sscore")), "score", "asd")))
pcw <- fread(file.path(AD, "pca_results/naive/cohort.sscore")); setnames(pcw, "#IID", "IID"); pcw <- pcw[, c(list(GUID = sub("^[0-9]+-", "", IID)), .SD), .SDcols = paste0("PC", 1:10, "_AVG")]; setnames(pcw, paste0("PC", 1:10, "_AVG"), paste0("PC", 1:10))
g2 <- merge(g2, pcw, by = "GUID")
rz <- function(v, d) { f <- as.formula(paste("v ~", paste0("PC", 1:10, collapse = " + "))); as.numeric(scale(residuals(lm(f, data = cbind(d, v = v))))) }
g2[, `:=`(pgs_cobes = rz(cobes, g2), pgs_bmi = rz(bmi, g2), pgs_cbmi = rz(cbmi, g2), pgs_ht = rz(ht, g2), pgs_asd = rz(asd, g2))]
g2 <- g2[, .(GUID, pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)]
# SSC: new scores (Expanse pgs_obesity_ssc) residualised on the ADSI PCs, z within SSC, mapped SS id -> SSC individual id
sscf <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/ssc_pgs_scores.tsv")   # built by 30a from the Expanse scores (participant-level, kept in Jon/Claude)
if (!file.exists(sscf)) stop("SSC PGS scores not found: run 30a_build_ssc_pgs.R after the Expanse job pgs_obesity_ssc finishes")
ss <- fread(sscf)   # columns: individual, pgs_cobes, pgs_bmi, pgs_cbmi (already residualised and z-scored by 30a)
cat("PGS: G2MH", nrow(g2), " SSC", nrow(ss), "\n")

# ---------------------------------------------------------------- analysis tables with sex and PGS
idmap <- unique(rbind(gmc[, .(id, cohort, sampleID, GUID)], ctl_g[, .(id, cohort, sampleID, GUID)]))
G <- A[pheno %in% c("Height", "Head circumference"), .(id, locus, sex, age, BMI, y, pheno, asd)]
G <- merge(G, idmap, by = "id", all.x = TRUE)
G[, key := fifelse(cohort == "SSC", paste0("SSC:", sampleID), fifelse(!is.na(GUID) & GUID != "", GUID, NA_character_))]
C <- X[, .(id, locus = fifelse(grp == "Unaffected relatives", "Unaffected controls", grp), sex, age, BMI, y, pheno, asd, cohort)]
C[, key := fifelse(cohort == "SSC", paste0("SSC:", id), id)]
ALL <- rbind(G[, .(key, locus, sex, age, BMI, y, pheno, asd, cohort)], C[, .(key, locus, sex, age, BMI, y, pheno, asd, cohort)])
ALL <- ALL[locus != "SSC ASD probands"]
P <- rbind(g2[, .(key = GUID, pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)], ss[, .(key = paste0("SSC:", individual), pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, pgs_asd)])
ALL <- merge(ALL, P, by = "key", all.x = TRUE)
ALL[, carrier := locus != "Unaffected controls"]
# ---- first-stage strength of the BMI instruments (JS 2026-11-02, reviewer MR request)
H1 <- unique(ALL[pheno == "Height" & !is.na(pgs_bmi) & !is.na(BMI), .(key, locus, carrier, sex, age, BMI, pgs_cobes, pgs_bmi, pgs_cbmi, pgs_ht, y)])
H1[, grp := fifelse(carrier, "carriers", "controls")][, agegrp := fifelse(age < 18, "child", "adult")]
fs <- function(d, v) { if (nrow(d) < 30) return(data.table(n = nrow(d), r2 = NA_real_, F = NA_real_, beta = NA_real_)); m <- lm(as.formula(paste("BMI ~", v, "+ sex + age")), data = d); m0 <- lm(BMI ~ sex + age, data = d); r2 <- 1 - sum(resid(m)^2) / sum(resid(m0)^2); data.table(n = nrow(d), r2 = round(r2, 4), F = round(anova(m0, m)$F[2], 1), beta = round(coef(m)[v], 3)) }
cat("\nfirst stage (BMI z ~ PGS + sex + age): incremental R2 and F\n")
for (v in c("pgs_cobes", "pgs_bmi", "pgs_cbmi")) { cat("==", v, "\n"); print(H1[, fs(.SD, v), by = .(grp, agegrp)][order(grp, agegrp)]); print(H1[carrier == TRUE, fs(.SD, v), by = .(locus)][order(locus)]) }
cat("\ncarriers with PGS by locus x age group (height rows):\n"); print(H1[carrier == TRUE, .N, by = .(locus, agegrp)][order(locus, agegrp)])
# reduced form: outcome ~ PGS within carriers vs controls, and the PGS x carrier interaction (height, children), per instrument
cat("\nreduced form, children: y ~ PGS*carrier + sex + age\n")
for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi")) { d <- H1[agegrp == "child"]; m <- lm(as.formula(paste("y ~", v, "* carrier + sex + age")), data = d); s <- summary(m)$coefficients; print(round(s[grep(v, rownames(s)), ], 4)) }
cat("\nmeasured BMI, children, same individuals: y ~ BMI*carrier + sex + age\n"); d <- H1[agegrp == "child"]; print(round(summary(lm(y ~ BMI * carrier + sex + age, data = d))$coefficients[c("BMI", "BMI:carrierTRUE"), ], 4))
