# CNV x BMI and CNV x age x BMI for behavioural / cognitive phenotypes present in both G2MH and SSC (JS 2026-10-24).
# Same GAM machinery as the main growth figure (17_main_with_ukbb.R), SSC non-carriers as controls, no UK Biobank; runs locally.
#
# Phenotypes (both cohorts):
#   SRS   : Social Responsiveness Scale raw total (65 items, 0-3). G2MH: child form (srs_child_*_collateral) or adult form
#           (srs_adult_*_collateral) summed from items (value - 1; 5 = unknown -> NA; 17 reverse-keyed items), >= 60 items,
#           scaled to 65. SSC: srs_parent total (probands, designated siblings; child form) and srs_adult total (parents; adult
#           form). Z within form against SSC unaffected (child form: designated siblings; adult form: parents).
#   SCQ   : Social Communication Questionnaire lifetime total. G2MH scq_total_collateral; SSC scq_life (probands, siblings).
#           Z against SSC designated siblings.
#   RBS-R : Repetitive Behavior Scale-Revised. The G2MH release carries 29 of the 43 items; sum of those 29 items in both
#           cohorts (SSC rbs_r_raw, probands only). Z against SSC probands.
#   FSIQ  : full-scale IQ. G2MH smry_iq_fsiq; SSC ssc_diagnosis_full_scale_iq (probands only). Z = (IQ - 100) / 15.
# Direction: SRS, SCQ, RBS-R higher = more impairment; FSIQ higher = better.
# Age at phenotype: G2MH scales_age (SRS/SCQ/RBS) or smry_iq_age (FSIQ); SSC: age in the growth genomatrix (single visit).
# BMI: raw BMI from the growth genomatrix (G2MH: visit within 3 y of the phenotype age; SSC: the visit), Z by 00_scoring.R
#   (WHO < 19 y; adults: centile among clinical non-carrier adults). Tails >= 70th / <= 30th centile for the plots.
# Controls = all SSC non-carriers with the phenotype (as for head circumference); sensitivity: SSC unaffected only (no probands).
# CNB accuracy / speed (G2MH only): carriers only, since G2MH non-carriers have no BMI; tests = BMI effect and age x BMI within carriers.
# Tests per phenotype x CNV (REML):
#   within    : y ~ s(age) + s(BMI) + ti(age, BMI)                                   p of ti(age,BMI), within a group
#   bmi_main  : y ~ BMI + s(age)                                                     age-adjusted BMI effect within a group
#   pooled carrier groups add a locus main effect; SSC-alone rows give BMI effects in the control set and ASD (probands) vs unaffected
#   cnv_x_bmi : y ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI) + is_cnv:BMI       age-adjusted constant shift (as in the
#               trajectory panels of the main figure)
#   threeway  : y ~ is_cnv + s(age) + s(BMI) + s(age, by = is_cnv) + ti(age, BMI) + ti(age, BMI, by = is_cnv)   p of the by-CNV ti
# Output: pheno_bmi_gam_tests.csv (ASD-adjusted + stratified; the pre-ASD run is kept as pheno_bmi_gam_tests.noASDadj_20261024.csv),
#   pheno_bmi_sample.csv, fig_supp_pheno_bmi_gam.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
DB   <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")
SSC  <- file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9")
source(file.path(HERE, "00_scoring.R"))
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
lab4 <- c("16p11.2 proximal DEL" = "16p11.2 del", "16p11.2 proximal DUP" = "16p11.2 dup", "22q11.2 A-D DEL" = "22q11.2 del", "22q11.2 A-D DUP" = "22q11.2 dup")

# ---------------------------------------------------------------- growth genomatrix: BMI, age, sex, CNV
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
gm <- gm[!(cohort == "G2MH" & Cardiff_mask %in% TRUE)][, CNV := trimws(CNV)][, sex := ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA))]
# JS 2026-10-24: SSC ages. The "age" in Age_individual_v14_parent_age.txt / SSC_Age_of_Subjects (and, for SSC parents, the
# age in Molly's genomatrix) is a snapshot age as of 2012-09, not age at measurement: it overstates visit age by ~3 y
# (median 3.2 y). Probands in the genomatrix carry age_at_ados (= age at evaluation), which is correct. Siblings and
# parents were measured at the proband evaluation visit, so visit age = proband "Age at Evaluation" (ssc_case_HC.csv)
# + (proband birth month - own birth month). This reproduces SSC's own sibling height Z (r = 0.96). See
# README_ages_and_line_endings.md in the SSC folder. Sibling ages -> `sibage`; SSC parent ages in gm are corrected here.
agef <- fread(file.path(SSC, "Age_individual_v14_parent_age.txt"))
setnames(agef, c("id()", "individual_properties.role"), c("id", "role0"))
agef[, bm := as.integer(substr(birth, 1, 4)) * 12L + as.integer(substr(birth, 6, 7))]
ssc_ev <- fread(file.path(SSC, "ssc_case_HC.csv"))[, .(proband = Proband_ID, age_eval = `Age at Evaluation (Non-PHI) :: Proband SSC Height, Weight, Head Circumference`)]
ssc_pro <- merge(agef[role0 == "proband", .(family, pbm = bm, proband = id)], ssc_ev, by = "proband")
visit_age <- function(roles) merge(agef[role0 %in% roles, .(id, family, sex = tolower(sex), obm = bm)], ssc_pro, by = "family")[, .(id, sex, age = (age_eval + pbm - obm) / 12)][!is.na(age)]
sibage <- visit_age("designated-sibling")
parage <- visit_age(c("mother", "father"))
gm[cohort == "SSC" & sampleID %in% parage$id, age := parage[match(sampleID, id), age]]
cat("SSC parent ages in genomatrix replaced by visit age (snapshot age - ~3 y):", sum(gm$cohort == "SSC" & gm$sampleID %in% parage$id), "rows\n")
# JS 2026-10-24: plausibility QC of height / BMI / head circumference before any scoring (rules in 00_qc_growth.R; flagged rows in implausible_height_weight.csv)
source(file.path(HERE, "00_qc_growth.R")); gm <- qc_growth(gm, label = "growth genomatrix")
REFad <- readRDS(file.path(DATA, "analysis_table.rds"))[locus == "No CNV" & age >= 19 & !is.na(BMI_raw) & BMI_raw >= 10 & BMI_raw <= 70]   # JS 2026-10-24: plausibility bound (00_qc_growth.R)[order(id, -age)][, .SD[1], by = id]
bmiZ <- function(b, a, s) score_bmi(b, a, s, REFad$BMI_raw, REFad$age, as.character(REFad$sex))

# ---------------------------------------------------------------- SSC phenotypes (all non-carrier SSC individuals in the genomatrix)
ssc <- gm[cohort == "SSC" & CNV == "No CNV" & !is.na(BMI_raw) & !is.na(age), .(id = sampleID, cohort = "SSC", locus = "No CNV", sex, age, BMI_raw)]
# JS 2026-10-24: SSC designated siblings have no height/weight in the growth genomatrix. Their BMI comes from the SSC
sibcu <- fread(file.path(SSC, "Designated Unaffected Sibling Data", "ssc_commonly_used.csv"))[, .(id = individual, bmi = suppressWarnings(as.numeric(bmi)))][!is.na(bmi)]
sib <- merge(sibcu, sibage, by = "id")
sib <- sib[id %in% gm[cohort == "SSC" & CNV == "No CNV", sampleID] & !(id %in% ssc$id) & !is.na(age), .(id, cohort = "SSC", locus = "No CNV", sex, age, BMI_raw = bmi)]
cat("SSC designated siblings added from SSC phenotype files (non-carriers, with BMI and age):", nrow(sib), "\n")
ssc <- rbind(ssc, sib)
ssc[, role := sub("^[0-9]+[.]", "", id)][, role := fcase(grepl("^p", role), "proband", grepl("^s", role), "sibling", role %in% c("mo", "fa"), "parent", default = "other")]
rd <- function(sub, f) {   # some SSC v15 files have CR-only line endings: normalise before parsing
  txt <- paste(readLines(file.path(SSC, sub, f), warn = FALSE), collapse = "\n"); txt <- gsub("\r\n|\r", "\n", txt); fread(text = txt)
}
srs_child <- rbind(rd("Proband Data", "srs_parent.csv")[, .(id = individual, srs = as.numeric(total))], rd("Designated Unaffected Sibling Data", "srs_parent.csv")[, .(id = individual, srs = as.numeric(total))])
srs_adult <- rbind(rd("Mother Data", "srs_adult.csv")[, .(id = individual, srs = as.numeric(total))], rd("Father Data", "srs_adult.csv")[, .(id = individual, srs = as.numeric(total))],
                   rd("Designated Unaffected Sibling Data", "srs_adult.csv")[, .(id = individual, srs = as.numeric(total))])
srs_child <- srs_child[!is.na(srs)][, form := "child"]; srs_adult <- srs_adult[!is.na(srs)][, form := "adult"]
srs_ssc <- rbind(srs_child, srs_adult)[!duplicated(id)]
scq_ssc <- rbindlist(lapply(c("Proband Data", "Designated Unaffected Sibling Data", "Other Sibling Data"), function(s) rd(s, "scq_life.csv")[, .(id = individual, scq = as.numeric(summary_score))]))[!is.na(scq)][!duplicated(id)]
rbs_raw <- rd("Proband Data", "rbs_r_raw.csv")
# G2MH release carries RBS-R items 3,4,6,7,8,9,11,13,17,18,19,20,21,24,25,26,28,29,31,32,33,34,35,36,37,38,39,41,43 (29 of 43)
RBS_ITEMS <- c(3, 4, 6, 7, 8, 9, 11, 13, 17, 18, 19, 20, 21, 24, 25, 26, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 43)
qcols <- sapply(RBS_ITEMS, function(i) grep(sprintf("^q%02d_", i), names(rbs_raw), value = TRUE)[1])
rbs_raw[, (qcols) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = qcols]
rbs_raw[, n_ans := rowSums(!is.na(.SD)), .SDcols = qcols][, rbs29 := ifelse(n_ans >= 26, rowSums(.SD, na.rm = TRUE) * 29 / n_ans, NA_real_), .SDcols = qcols]
rbs_ssc <- rbs_raw[!is.na(rbs29), .(id = individual, rbs = rbs29)][!duplicated(id)]
iq_ssc <- rd("Proband Data", "ssc_core_descriptive.csv")[, .(id = individual, fsiq = as.numeric(ssc_diagnosis_full_scale_iq))][!is.na(fsiq)]
ssc <- merge(ssc, srs_ssc, by = "id", all.x = TRUE); ssc <- merge(ssc, scq_ssc, by = "id", all.x = TRUE); ssc <- merge(ssc, rbs_ssc, by = "id", all.x = TRUE); ssc <- merge(ssc, iq_ssc, by = "id", all.x = TRUE)
ssc[, `:=`(age_pheno = age, BMI = bmiZ(BMI_raw, age, sex), unaffected = role %in% c("sibling", "parent"))]
cat("SSC non-carriers with BMI:", nrow(ssc), "; SRS", sum(!is.na(ssc$srs)), "(child", sum(ssc$form == "child", na.rm = TRUE), "adult", sum(ssc$form == "adult", na.rm = TRUE), ") SCQ", sum(!is.na(ssc$scq)), " RBS29", sum(!is.na(ssc$rbs)), " FSIQ", sum(!is.na(ssc$fsiq)), "\n")
cat("SSC SRS raw total by role x form (mean):\n"); print(ssc[!is.na(srs), .(n = .N, mean = round(mean(srs), 1), sd = round(sd(srs), 1)), by = .(role, form)][order(form, role)])

# ---------------------------------------------------------------- G2MH phenotypes (release) linked to the genomatrix by GUID
rel <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"))
num <- function(v) suppressWarnings(as.numeric(v))
REV <- c(3, 7, 11, 12, 15, 17, 21, 22, 26, 32, 38, 40, 43, 45, 48, 52, 55)   # SRS-2 reverse-keyed items (school-age and adult forms)
srs_sum <- function(prefix) {
  M <- as.matrix(rel[, paste0(prefix, 1:65, "_collateral"), with = FALSE]); M <- apply(M, 2, num); M[M == 5] <- NA; M <- M - 1
  M[, REV] <- 3 - M[, REV]; n <- rowSums(!is.na(M)); ifelse(n >= 60, rowSums(M, na.rm = TRUE) * 65 / n, NA_real_)
}
rel[, srs_c := srs_sum("srs_child_")][, srs_a := srs_sum("srs_adult_")]
rel[, srs := fifelse(!is.na(srs_c), srs_c, srs_a)][, form := fifelse(!is.na(srs_c), "child", fifelse(!is.na(srs_a), "adult", NA_character_))]
rb <- paste0("rbs_r_", RBS_ITEMS); rel[, (rb) := lapply(.SD, num), .SDcols = rb]
rel[, n_rbs := rowSums(!is.na(.SD)), .SDcols = rb][, rbs := ifelse(n_rbs >= 26, rowSums(.SD, na.rm = TRUE) * 29 / n_rbs, NA_real_), .SDcols = rb]
g2 <- rel[, .(GUID = guid, srs, form, scq = num(scq_total_collateral), rbs, fsiq = num(smry_iq_fsiq), scales_age = num(scales_age), iq_age = num(smry_iq_age), srs2_asd)]
cat("G2MH SRS computed total by form x srs2_asd category (validation: No ASD < Subclinical; blank = not categorised):\n")
print(g2[!is.na(srs), .(n = .N, mean = round(mean(srs), 1), sd = round(sd(srs), 1)), by = .(form, srs2_asd)][order(form, srs2_asd)])
gg <- gm[cohort == "G2MH" & CNV %in% names(lab4) & !is.na(GUID) & GUID != "", .(GUID, locus = lab4[CNV], sex, age, BMI_raw)]
nearest_bmi <- function(ph_age_col) {
  m <- merge(g2[!is.na(get(ph_age_col)), .(GUID, age_pheno = get(ph_age_col))], gg, by = "GUID", allow.cartesian = TRUE)
  m[, dage := abs(age - age_pheno)][order(GUID, dage)][, .SD[1], by = GUID][dage <= 3, .(GUID, age_pheno, locus, sex, BMI_raw)]
}
g2s <- merge(nearest_bmi("scales_age"), g2[, .(GUID, srs, form, scq, rbs)], by = "GUID")[, `:=`(cohort = "G2MH", id = GUID, unaffected = FALSE)]
g2i <- merge(nearest_bmi("iq_age"), g2[, .(GUID, fsiq)], by = "GUID")[, `:=`(cohort = "G2MH", id = GUID, unaffected = FALSE)]
g2s[, BMI := bmiZ(BMI_raw, age_pheno, sex)]; g2i[, BMI := bmiZ(BMI_raw, age_pheno, sex)]
cat("G2MH carriers with BMI near the scales age:", nrow(g2s), " (SRS", sum(!is.na(g2s$srs)), " SCQ", sum(!is.na(g2s$scq)), " RBS29", sum(!is.na(g2s$rbs)), "); near IQ age:", nrow(g2i), "(FSIQ", sum(!is.na(g2i$fsiq)), ")\n")

# ---------------------------------------------------------------- standardised outcomes, one table per phenotype
zref <- function(x, ref) (x - mean(ref, na.rm = TRUE)) / sd(ref, na.rm = TRUE)
mk <- function(ph) {
  if (ph == "SRS") {
    cs <- ssc[!is.na(srs), .(id, cohort, locus, role, unaffected, sex, age = age_pheno, BMI, raw = srs, form)]
    cg <- g2s[!is.na(srs), .(id, cohort, locus, role = "carrier", unaffected, sex, age = age_pheno, BMI, raw = srs, form)]
    d <- rbind(cs, cg)
    d[form == "child", y := zref(raw, ssc[form == "child" & role == "sibling", srs])]
    d[form == "adult", y := zref(raw, ssc[form == "adult" & role == "parent", srs])]
  } else if (ph == "SCQ") {
    d <- rbind(ssc[!is.na(scq), .(id, cohort, locus, role, unaffected, sex, age = age_pheno, BMI, raw = scq)], g2s[!is.na(scq), .(id, cohort, locus, role = "carrier", unaffected, sex, age = age_pheno, BMI, raw = scq)])
    d[, y := zref(raw, ssc[role == "sibling", scq])]
  } else if (ph == "RBS-R") {
    d <- rbind(ssc[!is.na(rbs), .(id, cohort, locus, role, unaffected, sex, age = age_pheno, BMI, raw = rbs)], g2s[!is.na(rbs), .(id, cohort, locus, role = "carrier", unaffected, sex, age = age_pheno, BMI, raw = rbs)])
    d[, y := zref(raw, ssc[role == "proband", rbs])]
  } else {
    d <- rbind(ssc[!is.na(fsiq), .(id, cohort, locus, role, unaffected, sex, age = age_pheno, BMI, raw = fsiq)], g2i[!is.na(fsiq), .(id, cohort, locus, role = "carrier", unaffected, sex, age = age_pheno, BMI, raw = fsiq)])
    d[, y := (raw - 100) / 15]
  }
  d[!is.na(y) & !is.na(BMI) & !is.na(age) & age >= 4 & age <= 65][, pheno := ph]
}
# ---- CNB accuracy / speed (G2MH only; Anjali's z vs G2MH non-carriers, composites as in 19_cnb_by_cnv_bmi.R) ----
# JS 2026-10-24 (Jon): G2MH non-carriers are the controls here. Only 3 of 224 CNB non-carriers have any BMI, so the controls
# enter without BMI (one reference curve); CNV x BMI and CNV x age x BMI vs controls are not estimable for the CNB.
cnb <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/cnb_processing/cnb_normalized.csv"))
CLAB <- c(no_cnv = "No CNV", `16p11.2_Deletion` = "16p11.2 del", `16p11.2_Duplication` = "16p11.2 dup", `22q11.2_Deletion` = "22q11.2 del", `22q11.2_Duplication` = "22q11.2 dup")
cnb <- cnb[final_group_class %in% names(CLAB) & !is.na(cnb_age)][, locus := CLAB[final_group_class]]
ACC <- c("cnb_pcet_acc2_z", "cnb_cpt_ptp_z", "cnb_lnb_mcr_z", "cnb_er40_cr_z", "cnb_pmat_pc_z", "cnb_volt_cr_z", "cnb_cpf_cr_z", "cnb_medf_pc_z", "cnb_adt_pc_z", "cnb_plot_pc_z")
SPD <- c("cnb_mpraxis_rtcr_z", "cnb_pcet_rtcr_z", "cnb_cpt_tprt_z", "cnb_lnb_mrtc_z", "cnb_er40_rtcr_z", "cnb_pmat_rtcr_z", "cnb_volt_rtcr_z", "cnb_cpf_rtcr_z", "cnb_medf_rtcr_z", "cnb_adt_rtcr_z", "cnb_plot_rtcr_z")
cnb[, n_acc := rowSums(!is.na(.SD)), .SDcols = ACC][, n_spd := rowSums(!is.na(.SD)), .SDcols = SPD]
cnb[, accuracy := ifelse(n_acc >= 5, rowMeans(.SD, na.rm = TRUE), NA_real_), .SDcols = ACC][, speed := ifelse(n_spd >= 5, -rowMeans(.SD, na.rm = TRUE), NA_real_), .SDcols = SPD]
cnb[, sex := ifelse(sex_final == "Male", "male", ifelse(sex_final == "Female", "female", NA))]
vb <- merge(cnb[, .(GUID, cnb_age)], gg[, .(GUID, age, BMI_raw)], by = "GUID", allow.cartesian = TRUE)[, dage := abs(age - cnb_age)][order(GUID, dage)][, .SD[1], by = GUID][dage <= 3]
cnb <- merge(cnb, vb[, .(GUID, bmi_age = age, bmi_raw = BMI_raw)], by = "GUID", all.x = TRUE)
relb <- rel[, .(rarecnv_id, rel_age = num(dem_age), h = num(med_height_m) * 100 + num(med_height_cm), w = num(med_weight))][!is.na(h) & !is.na(w) & h >= 60 & h <= 230 & w > 5 & w < 300 & !is.na(rel_age), .(rarecnv_id, rel_age, rel_bmi = w / (h / 100)^2)][rel_bmi >= 10 & rel_bmi <= 70]   # JS 2026-10-24: plausibility bounds as in 00_qc_growth.R
cnb <- merge(cnb, relb, by = "rarecnv_id", all.x = TRUE)
cnb[is.na(bmi_raw) & !is.na(rel_bmi) & abs(rel_age - cnb_age) <= 3, `:=`(bmi_age = rel_age, bmi_raw = rel_bmi)]
cnb[, BMI := bmiZ(bmi_raw, bmi_age, sex)]
mk_cnb <- function(meas, ph) {
  d <- cnb[!is.na(get(meas)), .(id = GUID, cohort = "G2MH", locus, role = ifelse(locus == "No CNV", "control", "carrier"), unaffected = locus == "No CNV", sex, age = cnb_age, BMI, raw = get(meas), y = get(meas), pheno = ph)]
  d[age >= 4 & age <= 65 & (locus == "No CNV" | !is.na(BMI))]
}
PHEN <- c("SRS", "SCQ", "RBS-R", "FSIQ", "CNB accuracy", "CNB speed")
CNB_PH <- c("CNB accuracy", "CNB speed")
D <- rbind(rbindlist(lapply(setdiff(PHEN, CNB_PH), mk), fill = TRUE), mk_cnb("accuracy", "CNB accuracy"), mk_cnb("speed", "CNB speed"), fill = TRUE)
D[, tail := fifelse(BMI >= qnorm(.7), "BMI \u2265 70th centile", fifelse(BMI <= qnorm(.3), "BMI \u2264 30th centile", NA_character_))]
smp <- D[, .(n = .N, n_hi = sum(tail == "BMI \u2265 70th centile", na.rm = TRUE), n_lo = sum(tail == "BMI \u2264 30th centile", na.rm = TRUE), age_med = as.numeric(median(age)), mean_y = mean(y)), by = .(pheno, locus)][order(pheno, locus)]
fwrite(smp, file.path(HERE, "pheno_bmi_sample.csv")); cat("\n== sample\n"); print(smp)

# ---------------------------------------------------------------- ASD status (JS 2026-10-24, Jon): covariate in every model + stratified sensitivity
# SSC: probands = ASD, siblings / parents / other siblings = no ASD. G2MH carriers: the ASD assignment used for panels F/L
# (Jon/Claude/G2MH_CNV_growth/subjects.csv, G2MH psychiatric assessment incl. Cardiff). Unknown -> "unknown" level (rows kept).
s_asd <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd"))
asdmap <- rbind(data.table(cohort = "SSC", id = ssc$id, asd = as.integer(ssc$role == "proband")),
                s_asd[!is.na(asd) & asd %in% c("0", "1"), .(cohort = "G2MH", id = sfari_id, asd = as.integer(asd))])
asdmap <- unique(asdmap, by = c("cohort", "id"))
D <- merge(D, asdmap, by = c("cohort", "id"), all.x = TRUE)
D[, asd_f := factor(fifelse(is.na(asd), "unknown", fifelse(asd == 1, "ASD", "noASD")), levels = c("noASD", "ASD", "unknown"))]
cat("ASD status by cohort x group (D rows):\n"); print(dcast(D[pheno == "SRS" | pheno == "CNB accuracy", .N, by = .(cohort, locus, asd_f)], cohort + locus ~ asd_f, value.var = "N", fill = 0))

RUN_TESTS <- toupper(Sys.getenv("PHENO_RUN_TESTS", "TRUE")) != "FALSE"   # PHENO_RUN_TESTS=FALSE -> skip the ~70 min test block, redraw the figure from the saved table
if (RUN_TESTS) {
# ---------------------------------------------------------------- GAM tests
# Per phenotype, (1) each CNV vs SSC controls; (2) pooled carriers (overgrowth = 16p del + 22q dup, undergrowth = 16p dup + 22q del,
# all four) vs SSC controls with locus as a covariate; (3) SSC alone: BMI effect and age x BMI in the control set, and ASD (probands)
# vs unaffected within SSC; (4) CNB: carriers only. All models carry ASD status as a covariate where it varies; the
# "ASD-stratified" rows compare like with like (ASD carriers vs SSC probands; non-ASD carriers vs SSC unaffected) and ASD vs
# non-ASD within carriers.
POOLS <- list("Overgrowth CNVs" = c("16p11.2 del", "22q11.2 dup"), "Undergrowth CNVs" = c("16p11.2 dup", "22q11.2 del"), "All CNVs" = LOCI)
cov_locus <- function(x) if (uniqueN(x$locus) > 2 || (uniqueN(x$locus) == 2 && !"No CNV" %in% x$locus)) "locus + " else ""
ASD_ADJ <- TRUE   # toggled below: the test block runs once without and once with the ASD covariate (Jon 2026-10-24: SCQ/SRS feed case assignment)
cov_asd <- function(x) if (ASD_ADJ && uniqueN(x$asd_f) > 1) "asd_f + " else ""
covs <- function(x) paste0(cov_locus(x), cov_asd(x))
within_p <- function(x) tryCatch({ f <- as.formula(paste("y ~", covs(x), "s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6))"))
  summary(gam(f, data = x, method = "REML"))$s.table["ti(age,BMI)", "p-value"] }, error = function(e) NA_real_)
bmi_main <- function(x) tryCatch({ f <- as.formula(paste("y ~", covs(x), "BMI + s(age, k = 8)")); co <- summary(gam(f, data = x, method = "REML"))$p.table["BMI", ]; c(co[1], co[2], co[4]) }, error = function(e) rep(NA_real_, 3))
threeway_p <- function(x) tryCatch({   # x has is_cnv (reference level first)
  f <- as.formula(paste("y ~", covs(x), "is_cnv + s(age, k = 10) + s(BMI, k = 6) + s(age, k = 10, by = is_cnv) + ti(age, BMI, k = c(6, 6)) + ti(age, BMI, k = c(6, 6), by = is_cnv)"))
  st <- summary(gam(f, data = x, method = "REML"))$s.table; st[grep("^ti\\(age,BMI\\):is_cnv", rownames(st))[1], "p-value"] }, error = function(e) NA_real_)
cnvbmi <- function(x) tryCatch({   # age-adjusted constant shift of the BMI slope in level 2 of is_cnv vs level 1
  x[, cnvBMI := as.numeric(is_cnv == levels(is_cnv)[2]) * BMI]
  f <- as.formula(paste("y ~", covs(x), "is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI"))
  co <- summary(gam(f, data = x, method = "REML"))$p.table; mn <- grep("^is_cnv", rownames(co))[1]
  c(co["cnvBMI", 1], co["cnvBMI", 2], co["cnvBMI", 4], co[mn, 1], co[mn, 4]) }, error = function(e) rep(NA_real_, 5))
mkz <- function(ctl, xc) rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
contrast <- function(ctl, xc, min_n = 30) {   # carrier-vs-control contrast: constant shift + three-way; NA if either side < min_n
  if (nrow(xc) < min_n || nrow(ctl) < min_n) return(rep(NA_real_, 6))
  z <- mkz(ctl, xc); e <- cnvbmi(copy(z)); c(e[4], e[5], e[1], e[2], e[3], threeway_p(copy(z)))
}
test_group <- function(ph, gname, loci_in) {
  xc <- D[pheno == ph & locus %in% loci_in]; ctl <- D[pheno == ph & locus == "No CNV"]; iscnb <- ph %in% CNB_PH
  base <- data.table(pheno = ph, group = gname, n_cnv = nrow(xc), n_cnv_asd = sum(xc$asd_f == "ASD"), n_cnv_noasd = sum(xc$asd_f == "noASD"),
                     n_cnv_hi = sum(xc$tail == "BMI ≥ 70th centile", na.rm = TRUE), n_cnv_lo = sum(xc$tail == "BMI ≤ 30th centile", na.rm = TRUE),
                     n_ctl = if (iscnb) NA_integer_ else nrow(ctl), n_ctl_unaffected = if (iscnb) NA_integer_ else sum(ctl$unaffected))
  r <- copy(base)[, comparison := if (iscnb) "carriers only (G2MH), ASD-adjusted" else "carriers vs SSC non-carriers, ASD-adjusted"]
  if (nrow(xc) < 30) return(r)
  bm <- bmi_main(copy(xc)); r[, `:=`(carrier_bmi_slope = bm[1], carrier_bmi_se = bm[2], carrier_bmi_p = bm[3], within_age_x_bmi_p = within_p(copy(xc)))]
  out <- r
  # ASD vs non-ASD within carriers (BMI slope in each stratum, and the ASD x BMI shift)
  xa <- xc[asd_f == "ASD"]; xn <- xc[asd_f == "noASD"]
  if (nrow(xa) >= 30) { b <- bmi_main(copy(xa)); out <- rbind(out, copy(base)[, `:=`(comparison = "ASD carriers only: BMI effect", n_cnv = nrow(xa), carrier_bmi_slope = b[1], carrier_bmi_se = b[2], carrier_bmi_p = b[3], within_age_x_bmi_p = within_p(copy(xa)))], fill = TRUE) }
  if (nrow(xn) >= 30) { b <- bmi_main(copy(xn)); out <- rbind(out, copy(base)[, `:=`(comparison = "non-ASD carriers only: BMI effect", n_cnv = nrow(xn), carrier_bmi_slope = b[1], carrier_bmi_se = b[2], carrier_bmi_p = b[3], within_age_x_bmi_p = within_p(copy(xn)))], fill = TRUE) }
  if (nrow(xa) >= 30 && nrow(xn) >= 30) { z <- rbind(xn, xa)[, is_cnv := factor(ifelse(asd_f == "ASD", "ASD", "noASD"), levels = c("noASD", "ASD"))]; e <- cnvbmi(copy(z))
    out <- rbind(out, copy(base)[, `:=`(comparison = "within carriers: ASD vs non-ASD", n_cnv = nrow(xa), n_ctl = nrow(xn), cnv_main_z = e[4], cnv_main_p = e[5], cnv_x_bmi = e[1], cnv_x_bmi_se = e[2], cnv_x_bmi_p = e[3], cnv_x_age_x_bmi_p = threeway_p(copy(z)))], fill = TRUE) }
  if (iscnb) { out[, note := "G2MH controls lack BMI; CNV x BMI vs controls not estimable"]; return(out) }
  if (nrow(ctl) < 30) return(out)
  e <- contrast(ctl, xc); out[1, `:=`(cnv_main_z = e[1], cnv_main_p = e[2], cnv_x_bmi = e[3], cnv_x_bmi_se = e[4], cnv_x_bmi_p = e[5], cnv_x_age_x_bmi_p = e[6])]
  eu <- contrast(ctl[unaffected == TRUE], xc); out[1, `:=`(unaff_cnv_x_bmi = eu[3], unaff_cnv_x_bmi_p = eu[5], unaff_cnv_x_age_x_bmi_p = eu[6])]
  # stratum-matched contrasts
  ea <- contrast(ctl[asd_f == "ASD"], xa); if (any(is.finite(ea))) out <- rbind(out, copy(base)[, `:=`(comparison = "ASD carriers vs SSC probands", n_cnv = nrow(xa), n_ctl = sum(ctl$asd_f == "ASD"), cnv_main_z = ea[1], cnv_main_p = ea[2], cnv_x_bmi = ea[3], cnv_x_bmi_se = ea[4], cnv_x_bmi_p = ea[5], cnv_x_age_x_bmi_p = ea[6])], fill = TRUE)
  en <- contrast(ctl[asd_f == "noASD"], xn); if (any(is.finite(en))) out <- rbind(out, copy(base)[, `:=`(comparison = "non-ASD carriers vs SSC unaffected", n_cnv = nrow(xn), n_ctl = sum(ctl$asd_f == "noASD"), cnv_main_z = en[1], cnv_main_p = en[2], cnv_x_bmi = en[3], cnv_x_bmi_se = en[4], cnv_x_bmi_p = en[5], cnv_x_age_x_bmi_p = en[6])], fill = TRUE)
  cat(ph, gname, "done\n"); out
}
run_tests <- function() rbindlist(fill = TRUE, lapply(PHEN, function(ph) rbind(
  rbindlist(fill = TRUE, lapply(LOCI, function(L) test_group(ph, L, L))),
  rbindlist(fill = TRUE, lapply(names(POOLS), function(P) test_group(ph, P, POOLS[[P]]))))))
ASD_ADJ <- FALSE; res0 <- run_tests()[, asd_adjustment := "none"]
ASD_ADJ <- TRUE;  res1 <- run_tests()[, asd_adjustment := fifelse(grepl("ASD-adjusted", comparison), "covariate", "stratified")][, comparison := sub(", ASD-adjusted", "", comparison)]
res <- rbind(res0[asd_adjustment == "none" & !grepl("ASD carriers|non-ASD carriers|within carriers", comparison)], res1, fill = TRUE)
# SSC alone: BMI and age x BMI in the control set (all, unaffected only, probands only), and ASD (probands) vs unaffected
ssc_rows <- rbindlist(fill = TRUE, lapply(setdiff(PHEN, CNB_PH), function(ph) {
  ctl <- D[pheno == ph & locus == "No CNV"]; un <- ctl[unaffected == TRUE]; pr <- ctl[role == "proband"]
  one <- function(x, g, comp) { b <- bmi_main(copy(x)); data.table(pheno = ph, group = g, comparison = comp, n_ctl = nrow(x), n_ctl_unaffected = sum(x$unaffected), carrier_bmi_slope = b[1], carrier_bmi_se = b[2], carrier_bmi_p = b[3], within_age_x_bmi_p = within_p(copy(x))) }
  out <- one(ctl, "SSC non-carriers", "SSC alone: BMI effect, ASD-adjusted")
  if (nrow(un) >= 30) out <- rbind(out, one(un, "SSC unaffected (siblings, parents)", "SSC alone: BMI effect within group"), fill = TRUE)
  if (nrow(pr) >= 30) out <- rbind(out, one(pr, "SSC probands (ASD)", "SSC alone: BMI effect within group"), fill = TRUE)
  if (nrow(un) >= 30 && nrow(pr) >= 30) { z <- rbind(un, pr)[, is_cnv := factor(ifelse(role == "proband", "ASD", "Unaffected"), levels = c("Unaffected", "ASD"))]; e <- cnvbmi(copy(z))
    out <- rbind(out, data.table(pheno = ph, group = "SSC probands vs unaffected", comparison = "SSC alone: ASD x BMI", n_cnv = nrow(pr), n_ctl = nrow(un), cnv_main_z = e[4], cnv_main_p = e[5], cnv_x_bmi = e[1], cnv_x_bmi_se = e[2], cnv_x_bmi_p = e[3], cnv_x_age_x_bmi_p = threeway_p(copy(z))), fill = TRUE) }
  out
}))
ssc_rows[, asd_adjustment := fifelse(grepl("ASD-adjusted", comparison), "covariate", "none (within-stratum)")]; ssc_rows[, comparison := sub(", ASD-adjusted", "", comparison)]
res <- rbind(res, ssc_rows, fill = TRUE)
setcolorder(res, intersect(c("pheno", "group", "comparison", "asd_adjustment", "n_cnv", "n_cnv_asd", "n_cnv_noasd", "n_cnv_hi", "n_cnv_lo", "n_ctl", "n_ctl_unaffected", "cnv_main_z", "cnv_main_p", "cnv_x_bmi", "cnv_x_bmi_se", "cnv_x_bmi_p", "cnv_x_age_x_bmi_p",
                             "carrier_bmi_slope", "carrier_bmi_se", "carrier_bmi_p", "within_age_x_bmi_p", "unaff_cnv_x_bmi", "unaff_cnv_x_bmi_p", "unaff_cnv_x_age_x_bmi_p", "note"), names(res)))
fwrite(res, file.path(HERE, "pheno_bmi_gam_tests.csv"))
cat("\n== tests\n")
print(res[, .(pheno, group = substr(group, 1, 18), comparison = substr(comparison, 1, 36), n_cnv, n_ctl, grp_bmi = round(carrier_bmi_slope, 3), p_grp_bmi = signif(carrier_bmi_p, 2), p_within = signif(within_age_x_bmi_p, 2),
              cnv_x_bmi = round(cnv_x_bmi, 3), p_cnv_x_bmi = signif(cnv_x_bmi_p, 2), p_3way = signif(cnv_x_age_x_bmi_p, 2), p_unaff_cxb = signif(unaff_cnv_x_bmi_p, 2))], nrows = 300)

cat("\n== tests (carrier rows: group BMI slope = within carriers; cnv_x_bmi = age-adjusted shift vs controls)\n")
print(res[, .(pheno, group = substr(group, 1, 26), n_cnv, n_ctl, grp_bmi = round(carrier_bmi_slope, 3), p_grp_bmi = signif(carrier_bmi_p, 2), p_within_age_bmi = signif(within_age_x_bmi_p, 2),
              cnv_x_bmi = round(cnv_x_bmi, 3), p_cnv_x_bmi = signif(cnv_x_bmi_p, 2), p_3way = signif(cnv_x_age_x_bmi_p, 2), p_unaff_cxb = signif(unaff_cnv_x_bmi_p, 2))], nrows = 200)

} else res <- fread(file.path(HERE, "pheno_bmi_gam_tests.csv"))   # figure-only run: reuse the saved table

# ---------------------------------------------------------------- figure: phenotype by age, BMI tails, per group (rows) x phenotype (columns)
if (toupper(Sys.getenv("PHENO_DRAW_FIG", "TRUE")) != "FALSE") {   # PHENO_DRAW_FIG=FALSE -> data + tests only (used by 22_16pdel_bmi90.R)
COLS <- c("BMI \u2265 70th centile" = "#F0645A", "BMI \u2264 30th centile" = "#00A5AD")
COLS3 <- c(COLS, "Non-carriers (no BMI)" = "#8C9A9A")
FIG_ADJ <- "none"   # which table rows the panels print: "none" (unadjusted) or "covariate" (ASD-adjusted)
YLAB <- c(SRS = "SRS total (Z vs SSC unaffected)", SCQ = "SCQ lifetime (Z vs SSC siblings)", `RBS-R` = "RBS-R 29 items (Z vs SSC probands)", FSIQ = "FSIQ (Z, norm 100/15)",
          `CNB accuracy` = "CNB accuracy (Z vs G2MH non-carriers)", `CNB speed` = "CNB speed (Z, higher = faster)")
GROUPS <- c("No CNV", LOCI)
panel <- function(ph, G, show_title) {
  ctl <- G == "No CNV"; iscnb <- ph %in% CNB_PH
  # panel statistics: the primary (unadjusted, cases + controls) carrier-vs-SSC contrast; the ASD-adjusted and stratified versions are in the table
  r <- res[pheno == ph & group == G & grepl("^carriers", comparison) & asd_adjustment == FIG_ADJ][1]
  x <- if (ctl && iscnb) D[pheno == ph & locus == G][, grp := factor("Non-carriers (no BMI)", names(COLS3))] else D[pheno == ph & locus == G & !is.na(tail)][, grp := factor(tail, names(COLS3))]
  lab <- if (ctl) "" else if (iscnb) gsub("= < ", "< ", sprintf("BMI p = %s\nage × BMI p = %s", pfx(r$carrier_bmi_p), pfx(r$within_age_x_bmi_p)), fixed = TRUE) else
    gsub("= < ", "< ", sprintf("CNV × BMI p = %s\nCNV × age × BMI p = %s", pfx(r$cnv_x_bmi_p), pfx(r$cnv_x_age_x_bmi_p)), fixed = TRUE)
  pv <- if (iscnb) c(r$carrier_bmi_p, r$within_age_x_bmi_p) else c(r$cnv_x_bmi_p, r$cnv_x_age_x_bmi_p)
  bold <- !ctl && isTRUE(suppressWarnings(min(pv, na.rm = TRUE)) < 0.05)
  yl <- range(quantile(D[pheno == ph, y], c(0.005, 0.995))); yl <- yl + c(-0.6, 0.9) * diff(yl) * 0.2
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = if (ctl) 0.5 else 1.2, alpha = if (ctl) 0.15 else 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 6), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = 65, y = yl[2], label = lab, hjust = 1, vjust = 1, size = 3, fontface = if (bold) "bold" else "plain") +
    scale_colour_manual(values = COLS3, breaks = names(COLS3)) + scale_fill_manual(values = COLS3, breaks = names(COLS3)) +
    coord_cartesian(xlim = c(4, 65), ylim = yl) +
    labs(x = "Age (years)", y = YLAB[[ph]], colour = NULL, fill = NULL, title = if (show_title) ph else NULL,
         subtitle = sprintf("%s (n = %d)", if (ctl) (if (iscnb) "G2MH non-carriers" else "SSC non-carriers") else G, nrow(x))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 7), legend.position = "inside", legend.position.inside = c(0.22, 0.1),
                                          legend.background = element_blank(), legend.key.height = unit(7, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
}
# JS 2026-10-24 (Jon): no control row for the CNB (controls have no BMI; CNB tests are within carriers)
pl <- list(); for (ph in PHEN) for (G in GROUPS) pl[[paste(ph, G)]] <- if (ph %in% CNB_PH && G == "No CNV") ggplot() + theme_void() + labs(title = ph) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) else panel(ph, G, G == "No CNV")
fig <- wrap_plots(pl, ncol = length(PHEN), byrow = FALSE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(HERE, "fig_supp_pheno_bmi_gam.png"), fig, width = 28, height = 20, dpi = 130)
tryCatch(ggsave(file.path(HERE, "fig_supp_pheno_bmi_gam.pdf"), fig, width = 28, height = 20, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_supp_pheno_bmi_gam.pdf"), fig, width = 28, height = 20))
cat("done\n")
}
