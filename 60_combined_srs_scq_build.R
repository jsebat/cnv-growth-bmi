# Combined SRS / SCQ carrier sample: G2MH + Simons VIP + Simons Searchlight + SSC carriers (JS / Claude 2026-11-14; Jon: combine
# SSC, Searchlight, SVIP and G2MH for the CNV x BMI analysis). Until 2026-11-14 the SRS/SCQ analyses used G2MH carriers only.
# Carriers of 16p11.2 (BP4-BP5) and 22q11.2 (A-D) deletions / duplications:
#   G2MH  : D rows (21_pheno_bmi_gam.R; G2MH release, BMI visit within 3 y of the scales age)
#   SVIP  : Simons_VIP__16p11.2_Dataset_v10.0 srs_parent / srs_adult / scq_life; status svip_subjects.genetic_status_16p
#   Searchlight: Simons_Searchlight_Dataset_v14.0/<16p11.2 deletion|duplication>/measure_data srs_sa / srs_adult / scq_life
#           (genetic_status), only individuals without an SVIP record of the same measure (the SVIP families appear in both)
#   SSC   : SSC carriers of the four CNVs in the growth genomatrix, srs_parent / srs_adult / scq_life (21: srs_ssc, scq_ssc)
#   BMI for SVIP / Searchlight / SSC from the growth genomatrix visits (gmc), nearest the evaluation within 3 y; z by 00_scoring.
# Non-carrier references: SSC unaffected relatives (main reference, as in Fig. 5); SVIP / Searchlight non-carrier relatives kept
# for a within-cohort check. Scale: Box-Cox with the SSC-relative lambdas (00_boxcox_srs_scq.R), z vs SSC relatives.
# Output (controlled, participant level): DATA/combined_srs_scq.rds; sample counts: combined_srs_scq_sample.csv (aggregate).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
num <- function(v) suppressWarnings(as.numeric(v)); LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
bcz <- function(raw, f) { ph <- if (f == "scq") "SCQ" else "SRS"; r <- D[pheno == ph & form2 == f & unaffected == TRUE & locus == "No CNV", raw]; L <- LAM[[f]]
  bc <- function(v) if (abs(L) < 1e-8) log(v + 1) else ((v + 1)^L - 1) / L; (bc(raw) - mean(bc(r))) / sd(bc(r)) }
near_bmi <- function(P, cohorts) {   # P: id, pheno, form2, raw, age_ph, locus, sex_ph; ids in genomatrix format (14700.x3)
  g <- gmc[cohort %in% cohorts & !is.na(BMI_raw) & !is.na(age), .(id = sampleID, age_v = age, BMI_raw, sex_g = sex)]
  m <- merge(P, g, by = "id", allow.cartesian = TRUE)[, dage := abs(age_v - age_ph)][order(id, pheno, dage)][, .SD[1], by = .(id, pheno)][dage <= 3]
  m[, sex := fifelse(!is.na(sex_g) & sex_g != "", sex_g, sex_ph)][, BMI := bmiZ(BMI_raw, age_ph, sex)][, age := age_ph][] }
# ---------------------------------------------------------------- SVIP
SV <- file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0")
sm <- fread(file.path(SV, "svip_subjects.csv"))[, .(id = sfari_id, status = genetic_status_16p, sex_ph = tolower(sex))]
rdsv <- function(f, tot, form, ph) { x <- fread(file.path(SV, f)); data.table(id = x$individual, raw = num(x[[tot]]), age_ph = num(x[[grep("eval_age_months", names(x), value = TRUE)[1]]]) / 12, form2 = form, pheno = ph) }
sv <- rbind(rdsv("srs_parent.csv", "srs_parent.total", "child", "SRS"), rdsv("srs_adult.csv", grep("^srs_adult[.]total$", names(fread(file.path(SV, "srs_adult.csv"), nrows = 1)), value = TRUE), "adult", "SRS"),
            rdsv("scq_life.csv", "scq_life.summary_score", "scq", "SCQ"))[!is.na(raw) & !is.na(age_ph)]
sv <- merge(sv, sm, by = "id")[, locus := fcase(status == "deletion", "16p11.2 del", status == "duplication", "16p11.2 dup", status == "negative", "Non-carrier relative", default = NA_character_)][!is.na(locus)]
sv <- sv[order(id, pheno, form2 == "adult")][, .SD[1], by = .(id, pheno)]
svb <- near_bmi(sv[, .(id, pheno, form2, raw, age_ph, locus, sex_ph)], c("SVIP", "Searchlight"))[, cohort := "SVIP"]
# ---------------------------------------------------------------- Searchlight v14
SL <- file.path(DB, "G2MH/Simons_Searchlight/Simons_Searchlight_Dataset_v14.0")
rdsl <- function(cn, f, tot, form, ph) { x <- fread(file.path(SL, cn, "measure_data", f))
  data.table(id = x$sfari_id, raw = num(x[[tot]]), age_raw = num(x$age_at_evaluation), status = x$genetic_status, sex_ph = tolower(x$sex), form2 = form, pheno = ph) }
sl <- rbindlist(lapply(c("16p11.2 deletion", "16p11.2 duplication"), function(cn) rbind(rdsl(cn, "srs_sa.csv", "total", "child", "SRS"), rdsl(cn, "srs_adult.csv", "total", "adult", "SRS"), rdsl(cn, "scq_life.csv", "summary_score", "scq", "SCQ"))))
cat("Searchlight age_at_evaluation range:", range(sl$age_raw, na.rm = TRUE), " status:\n"); print(sl[, .N, by = status])
sl[, age_ph := if (max(age_raw, na.rm = TRUE) > 120) age_raw / 12 else age_raw]   # months if values exceed 120
sl[, id := sub("-", ".", id, fixed = TRUE)]
sl[, locus := fcase(grepl("^16p11.2 deletion|deletion", status, ignore.case = TRUE) & !grepl("negative|not|false", status, ignore.case = TRUE), "16p11.2 del",
                    grepl("duplication", status, ignore.case = TRUE) & !grepl("negative|not|false", status, ignore.case = TRUE), "16p11.2 dup",
                    grepl("negative|familial.*negative|non.?carrier", status, ignore.case = TRUE), "Non-carrier relative", default = NA_character_)]
sl <- sl[!is.na(raw) & !is.na(age_ph) & !is.na(locus)][order(id, pheno, form2 == "adult")][, .SD[1], by = .(id, pheno)]
sl <- sl[!paste(id, pheno) %in% paste(sv$id, sv$pheno)]   # SVIP record takes precedence
slb <- near_bmi(sl[, .(id, pheno, form2, raw, age_ph, locus, sex_ph)], c("Searchlight", "SVIP"))[, cohort := "Searchlight"]
# ---------------------------------------------------------------- SSC carriers
scar <- gmc[cohort == "SSC" & CNV %in% names(lab4) & !is.na(BMI_raw) & !is.na(age), .(id = sampleID, locus = lab4[CNV], sex, age, BMI_raw)][order(id, -age)][, .SD[1], by = id]
sscp <- rbind(merge(scar, srs_ssc, by = "id")[, .(id, locus, sex, age, BMI_raw, raw = srs, form2 = form, pheno = "SRS")], merge(scar, scq_ssc, by = "id")[, .(id, locus, sex, age, BMI_raw, raw = scq, form2 = "scq", pheno = "SCQ")])
sscp[, BMI := bmiZ(BMI_raw, age, sex)][, cohort := "SSC"]
# ---------------------------------------------------------------- G2MH carriers and SSC relatives (as in Fig. 5)
g2c <- D[cohort == "G2MH" & locus %in% LOC & pheno %in% c("SRS", "SCQ"), .(id, cohort, locus, sex, age, BMI, raw, form2, pheno, asd)]
ref <- D[cohort == "SSC" & locus == "No CNV" & unaffected == TRUE & pheno %in% c("SRS", "SCQ"), .(id, cohort = "SSC", locus = "SSC relative", sex, age, BMI, raw, form2, pheno, asd)]
C <- rbind(g2c, svb[, .(id, cohort, locus, sex, age, BMI, raw, form2, pheno, asd = NA_integer_)], slb[, .(id, cohort, locus, sex, age, BMI, raw, form2, pheno, asd = NA_integer_)],
           sscp[, .(id, cohort, locus, sex, age, BMI, raw, form2, pheno, asd = NA_integer_)], ref)
C <- C[!is.na(BMI) & !is.na(raw) & age >= 4 & age <= 65]
C[, y := mapply(function(r, f) bcz(r, f), raw, form2)][, male := as.numeric(sex == "male")]
C[, is_cnv := as.numeric(locus %in% LOC)]
saveRDS(C, file.path(DATA, "combined_srs_scq.rds"))
smp <- dcast(C[, .N, by = .(pheno, cohort, locus)], pheno + locus ~ cohort, value.var = "N", fill = 0); fwrite(smp, file.path(HERE, "combined_srs_scq_sample.csv"))
cat("\ncombined sample:\n"); print(smp)
