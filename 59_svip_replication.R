# Replication in Simons VIP (16p11.2) of the BMI effects on SRS / SCQ and of the BMI x autism score interaction (JS / Claude
# 2026-11-14). Jon: the paper combined SVIP and G2MH; the SRS/SCQ analyses (21_pheno_bmi_gam.R) had used G2MH only.
# SVIP phenotypes: Simons_VIP__16p11.2_Dataset_v10.0 srs_parent.total (child form), srs_adult.total (adult form),
# scq_life.summary_score, eval ages in months; carrier status svip_subjects.genetic_status_16p (deletion / duplication; negative =
# non-carrier family member). BMI: SVIP rows of the growth genomatrix (gmc), visit nearest the evaluation within 3 y, z by 00_scoring.
# Scale: the same Box-Cox transform and SSC reference as the main analysis (00_boxcox_srs_scq.R lambdas; z vs SSC relatives).
# Scores: Anjali's SVIP genomes (data/searchlight/pgs/analysis_dataset.tsv; searchlight_id -> sfari_iid via Searchlight.WGS2
# metadata); autism = Grove et al. 2019 (iPSYCH-PGC), the score available for SVIP; residualised on 10 global PCs, z within SVIP.
# For a like-for-like comparison G2MH is re-run with the same Grove 2019 score (Anjali scores_wgs/asd.sscore).
# Tests: (1) BMI slope in SVIP carriers; CNV x BMI vs SVIP non-carrier relatives and vs SSC relatives; (2) BMI x autism score within
# SVIP carriers; (3) G2MH with the Grove score; (4) G2MH + SVIP pooled with a cohort term. Output: svip_replication.csv.
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s57 <- readLines(file.path(HERE, "57_cnv_bmi_pgs_threeway.R")); eval(parse(text = s57[1:grep("^out <- list\\(\\); add", s57)]))
AN <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data"); SV <- file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0")
fx <- function(v) sub(".", "-", v, fixed = TRUE); num <- function(v) suppressWarnings(as.numeric(v))
sm <- fread(file.path(SV, "svip_subjects.csv"))[, .(id = sfari_id, status = genetic_status_16p, sex_sv = sex, rel = relationship_to_iip)]
sp <- fread(file.path(SV, "srs_parent.csv"))[, .(id = individual, raw = num(srs_parent.total), eval_m = num(srs_parent.measure.eval_age_months), form2 = "child", pheno = "SRS")]
sa <- fread(file.path(SV, "srs_adult.csv")); sa <- sa[, .(id = individual, raw = num(get(grep("^srs_adult[.]total$", names(sa), value = TRUE))), eval_m = num(get(grep("eval_age_months", names(sa), value = TRUE)[1])), form2 = "adult", pheno = "SRS")]
sq <- fread(file.path(SV, "scq_life.csv")); sq <- sq[, .(id = individual, raw = num(scq_life.summary_score), eval_m = num(get(grep("eval_age_months", names(sq), value = TRUE)[1])), form2 = "scq", pheno = "SCQ")]
P2 <- rbind(sp, sa, sq)[!is.na(raw)]; P2 <- P2[order(id, pheno, form2 == "adult")][, .SD[1], by = .(id, pheno)]   # one SRS per person, child form preferred
P2 <- merge(P2, sm, by = "id"); P2[, age_ph := eval_m / 12]
cat("SVIP phenotype rows:\n"); print(P2[, .N, by = .(pheno, form2, status)][order(pheno, form2, status)])
# BMI: SVIP genomatrix visits
gs <- gmc[cohort == "SVIP" & !is.na(BMI_raw) & !is.na(age)]
cat("genomatrix SVIP sampleID examples:", head(unique(gs$sampleID), 4), "\n")
gs[, id := sampleID]; if (!any(P2$id %in% gs$id)) gs[, id := sub("-", ".", sampleID, fixed = TRUE)]
mm <- merge(P2, gs[, .(id, age_v = age, BMI_raw, sex_g = sex)], by = "id", allow.cartesian = TRUE)[, dage := abs(age_v - age_ph)][order(id, pheno, dage)][, .SD[1], by = .(id, pheno)][dage <= 3]
mm[, sex := fifelse(!is.na(sex_g), sex_g, tolower(sex_sv))][, BMI := bmiZ(BMI_raw, age_ph, sex)][, male := as.numeric(sex == "male")][, age := age_ph]
for (f in c("child", "adult", "scq")) { ph <- if (f == "scq") "SCQ" else "SRS"; r <- D[pheno == ph & form2 == f & unaffected == TRUE & locus == "No CNV", raw]; L <- LAM[[f]]
  bc <- function(v) if (abs(L) < 1e-8) log(v + 1) else ((v + 1)^L - 1) / L; mm[form2 == f, y := (bc(raw) - mean(bc(r))) / sd(bc(r))] }
mm[, locus := fcase(status == "deletion", "16p11.2 del", status == "duplication", "16p11.2 dup", status == "negative", "Non-carrier relative", default = NA_character_)]
mm <- mm[!is.na(locus) & !is.na(BMI) & !is.na(y) & age >= 4]
cat("SVIP with BMI and trait:\n"); print(mm[, .N, by = .(pheno, locus)][order(pheno, locus)])
# scores (Grove 2019 autism + adiposity), residualised on 10 global PCs within the SVIP genomes
sc <- fread(file.path(AN, "searchlight/pgs/analysis_dataset.tsv"), select = c("IID", paste0("GLOBAL_PC", 1:10), "PGS_asd", "PGS_childhood_obesity", "PGS_bmi_childhood", "PGS_bmi"))
meta <- unique(fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/Searchlight.WGS2.2023_05.sample_metadata.tsv"))[, .(searchlight_id, sfari_iid)])
sc <- merge(sc, meta, by.x = "IID", by.y = "searchlight_id")[, id := sub("-", ".", sfari_iid, fixed = TRUE)]
rzg <- function(v) as.numeric(scale(residuals(lm(as.formula(paste("v ~", paste0("GLOBAL_PC", 1:10, collapse = " + "))), data = cbind(sc, v = v)))))
sc[, `:=`(pgs_asd = rzg(PGS_asd), pgs_cobes = rzg(PGS_childhood_obesity), pgs_cbmi = rzg(PGS_bmi_childhood), pgs_bmi = rzg(PGS_bmi))]
if (!any(mm$id %in% sc$id)) sc[, id := sfari_iid]
mm <- merge(mm, sc[, .(id, pgs_asd, pgs_cobes, pgs_cbmi, pgs_bmi)], by = "id", all.x = TRUE)
cat("SVIP with autism score:\n"); print(mm[!is.na(pgs_asd), .N, by = .(pheno, locus)][order(pheno, locus)])
# ---------------------------------------------------------------- tests
SSR <- X[cohort == "SSC" & grp == "Unaffected relatives"]   # SSC relatives (Box-Cox y), reference for CNV x BMI
for (ph in c("SRS", "SCQ")) {
  car <- mm[pheno == ph & locus %in% c("16p11.2 del", "16p11.2 dup")][, locus := factor(locus)]; rel <- mm[pheno == ph & locus == "Non-carrier relative"]
  B <- "SVIP: BMI effect"
  for (g in c("16p11.2 del", "16p11.2 dup", "both")) { d <- if (g == "both") car else car[locus == g]; if (nrow(d) < 10) next
    m <- lm(as.formula(paste("y ~", if (g == "both") "locus +" else "", "poly(age, 2) + male + BMI")), d); add(ph, B, paste("carrier BMI slope,", g), c(n = nrow(d), rob(m, "BMI"))) }
  if (nrow(rel) >= 10) { add(ph, B, "non-carrier relatives BMI slope", c(n = nrow(rel), rob(lm(y ~ poly(age, 2) + male + BMI, rel), "BMI")))
    z <- rbind(car[, .(y, BMI, age, male, locus = as.character(locus), is_cnv = 1)], rel[, .(y, BMI, age, male, locus = "none", is_cnv = 0)])
    add(ph, B, "CNV x BMI vs SVIP non-carrier relatives", c(n = nrow(z), rob(lm(y ~ is_cnv * (poly(age, 2) + male) + locus + is_cnv * BMI, z), "is_cnv:BMI"))) }
  z <- rbind(car[, .(y, BMI, age, male, locus = as.character(locus), is_cnv = 1)], SSR[pheno == ph, .(y, BMI, age, male = as.numeric(sex == "male"), locus = "none", is_cnv = 0)])
  add(ph, B, "CNV x BMI vs SSC relatives", c(n = nrow(z), rob(lm(y ~ is_cnv * (poly(age, 2) + male) + locus + is_cnv * BMI, z), "is_cnv:BMI")))
  B <- "SVIP: BMI x autism score (Grove 2019)"
  d <- car[!is.na(pgs_asd)]
  if (nrow(d) >= 15) { m <- lm(y ~ locus + poly(age, 2) + male + BMI * pgs_asd, d); add(ph, B, "within SVIP carriers", c(n = nrow(d), rob(m, "BMI:pgs_asd")))
    med <- median(d$pgs_asd); for (h in c("below", "above")) { dh <- d[if (h == "below") pgs_asd < med else pgs_asd >= med]
      add(ph, B, paste("BMI slope, autism score", h, "median"), c(n = nrow(dh), rob(lm(y ~ locus + poly(age, 2) + male + BMI, dh), "BMI"))) }
    dr <- rel[!is.na(pgs_asd)]; if (nrow(dr) >= 15) add(ph, B, "within SVIP non-carrier relatives", c(n = nrow(dr), rob(lm(y ~ poly(age, 2) + male + BMI * pgs_asd, dr), "BMI:pgs_asd"))) }
  # G2MH with the same Grove 2019 score
  B <- "G2MH: BMI x autism score (Grove 2019)"
  gr <- rdw("asd"); gr <- merge(gr, pcw, by = "GUID"); gr[, pgs_grove := rz(score, gr)]
  dg <- merge(W[pheno == ph], gr[, .(key = GUID, pgs_grove)], by = "key")
  m <- lm(as.formula(paste("y ~", BASEW, "+ BMI * pgs_grove")), dg); add(ph, B, "within G2MH carriers", c(n = nrow(dg), rob(m, "BMI:pgs_grove")))
  cat("correlation of the Grove 2019 and April 2026 autism scores in G2MH:", round(cor(dg$pgs_grove, dg$pgs_asd), 2), "\n")
  # pooled G2MH + SVIP carriers, Grove score, cohort term
  B <- "G2MH + SVIP carriers pooled (Grove 2019)"
  pp <- rbind(dg[, .(y, BMI, age, male, locus = as.character(locus), pgs = pgs_grove, cohort = "G2MH")], d[, .(y, BMI, age, male, locus = as.character(locus), pgs = pgs_asd, cohort = "SVIP")])
  m <- lm(y ~ cohort + locus + locus:poly(age, 2) + male + BMI * pgs, pp); add(ph, B, "BMI x autism score, + cohort", c(n = nrow(pp), rob(m, "BMI:pgs")))
  m2 <- lm(y ~ cohort * (BMI * pgs) + locus + locus:poly(age, 2) + male, pp); add(ph, B, "difference between cohorts (cohort x BMI x score)", c(n = nrow(pp), rob(m2, "cohortSVIP:BMI:pgs")))
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "svip_replication.csv"))
options(width = 230); print(R[, .(trait, block, spec, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 100)
