# Non-carrier ASD vs no-ASD reference, within family (JS 2026-10-13): SSC proband vs designated unaffected sibling (.s1),
# measured at the same visit (Proband Data/ssc_hwhc.csv: height, head_circumference, weight and sibling1_*).
# Families restricted to non-carrier probands (genomatrix CNV == "No CNV"; carrier families excluded).
# Ages: proband = ssc_case_HC age at evaluation (months); sibling = proband age + (proband birth - sibling birth)
#   (birth year-month from SSC_Age_of_Subjects). Pairs kept when both are 2-18 y.
# Z-scores recomputed identically for both members with childsds: height WHO (< 5 y) / WHO 2007, head circumference
#   UK90, BMI WHO / WHO 2007 -- the same references as the carrier panels.
# Tests: paired difference (proband - sibling), paired t-test; and mixed model y ~ ASD + (1 | family) as a check.
# Output: ssc_sibpairs_results.csv (aggregates only)

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(childsds) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
SSC <- file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9")

hw <- fread(file.path(SSC, "Proband Data", "ssc_hwhc.csv"))
ag <- fread(file.path(SSC, "SSC_Age_of_Subjects_&_ParentalAgeAtBirth.csv"))
setnames(ag, "id()", "ind")
ev <- fread(file.path(SSC, "ssc_case_HC.csv"))[, .(proband = Proband_ID, age_eval = `Age at Evaluation (Non-PHI) :: Proband SSC Height, Weight, Head Circumference`)]
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
noncarrier_p <- unique(gm[cohort == "SSC" & trimws(CNV) == "No CNV" & grepl("[.]p1$", sampleID), sampleID])

ym <- function(s) { s <- as.character(s); as.numeric(substr(s, 1, 4)) * 12 + as.numeric(substr(s, 6, 7)) }
d <- hw[individual %in% noncarrier_p, .(proband = individual, family = sub("[.]p1$", "", individual),
        p_h = height, p_hc = head_circumference, p_w = weight, s_h = sibling1_height, s_hc = sibling1_head_circumference, s_w = sibling1_weight)]
d <- merge(d, ev, by = "proband")
d <- merge(d, ag[, .(proband = ind, p_sex = sex, p_birth = ym(birth))], by = "proband")
d <- merge(d, ag[grepl("[.]s1$", ind), .(family = sub("[.]s1$", "", ind), s_sex = sex, s_birth = ym(birth))], by = "family")
d[, `:=`(p_age = age_eval / 12, s_age = (age_eval + (p_birth - s_birth)) / 12)]
d <- d[p_age >= 2 & p_age < 18 & s_age >= 2 & s_age < 18]

source(file.path(HERE, "00_scoring.R"))   # shared scorer (same references as every other panel)
z <- function(val, age, sex, what) switch(what, height = score_height(val, age, sex), hc = score_hc(val, age, sex), bmi = score_bmi_child(val, age, sex))
d[, `:=`(p_hZ = z(p_h, p_age, p_sex, "height"), s_hZ = z(s_h, s_age, s_sex, "height"),
         p_hcZ = z(p_hc, p_age, p_sex, "hc"), s_hcZ = z(s_hc, s_age, s_sex, "hc"),
         p_bmiZ = z(p_w / (p_h / 100)^2, p_age, p_sex, "bmi"), s_bmiZ = z(s_w / (s_h / 100)^2, s_age, s_sex, "bmi"))]

res <- rbindlist(lapply(list(Height = c("p_hZ", "s_hZ"), `Head circumference` = c("p_hcZ", "s_hcZ"), BMI = c("p_bmiZ", "s_bmiZ")), function(v) {
  x <- d[!is.na(get(v[1])) & !is.na(get(v[2]))]
  pr <- x[[v[1]]]; sb <- x[[v[2]]]; tt <- t.test(pr, sb, paired = TRUE)
  data.table(pairs = nrow(x), mean_sibling = mean(sb), se_sibling = sd(sb) / sqrt(nrow(x)), mean_proband = mean(pr), se_proband = sd(pr) / sqrt(nrow(x)),
             paired_diff = unname(tt$estimate), paired_se = unname(tt$stderr), paired_p = tt$p.value,
             mean_age_proband = mean(x$p_age), mean_age_sibling = mean(x$s_age))
}), idcol = "outcome")
fwrite(res, file.path(HERE, "ssc_sibpairs_results.csv")); print(res, digits = 3)
