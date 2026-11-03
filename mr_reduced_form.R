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
# ---- reduced-form "genetic exposure" analysis (JS 2026-11-02): outcome ~ PGS x carrier + sex + age, childhood obesity PGS vs adult BMI PGS
U <- unique(ALL[!is.na(pgs_cobes) & !is.na(BMI) & sex %in% c("male", "female"), .(key, pheno, locus, carrier, sex, age, BMI, y, pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht)])
U[, grp := fifelse(carrier, "carriers", "controls")]
bands <- list(`Height 2-9` = quote(pheno == "Height" & age >= 2 & age < 10), `Height 10-18` = quote(pheno == "Height" & age >= 10 & age < 19),
              `Height adult` = quote(pheno == "Height" & age >= 19), `Head circumference child` = quote(pheno == "Head circumference" & age < 18),
              `SRS` = quote(pheno == "SRS"), `SCQ` = quote(pheno == "SCQ"))
res <- list()
for (nm in names(bands)) for (G in c("all carriers", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")) for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi")) {
  d <- U[eval(bands[[nm]])]; d <- if (G == "all carriers") d else d[locus == G | !carrier]
  if (d[carrier == TRUE, .N] < 15) next
  f <- as.formula(paste("y ~", v, "* carrier + sex + age + pgs_ht")); m <- lm(f, data = d); s <- summary(m)$coefficients
  mc <- lm(as.formula(paste("y ~", v, "+ sex + age + pgs_ht")), data = d[carrier == TRUE]); sc <- summary(mc)$coefficients[v, ]
  res[[length(res) + 1]] <- data.table(outcome = nm, group = G, pgs = v, n_carriers = d[carrier == TRUE, .N], n_controls = d[carrier == FALSE, .N],
    slope_controls = s[v, 1], slope_carriers = sc[1], se_carriers = sc[2], p_carriers = sc[4], interaction = s[paste0(v, ":carrierTRUE"), 1], int_se = s[paste0(v, ":carrierTRUE"), 2], int_p = s[paste0(v, ":carrierTRUE"), 4])
}
R <- rbindlist(res); fwrite(R, file.path(HERE, "mr_reduced_form_pgs.csv"))
print(R[, .(outcome, group, pgs, n_carriers, ctl = round(slope_controls, 3), car = round(slope_carriers, 3), p_car = signif(p_carriers, 2), int = round(interaction, 3), int_p = signif(int_p, 2))], nrows = 200)
