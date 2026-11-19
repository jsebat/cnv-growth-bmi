# Autism polygenic score by sex, and the three-way interaction with sex-standardised scores (JS / Claude 2026-11-14; Jon: in ASD
# family samples females carry higher autism scores in cases and controls). Data of 66 (combined genotyped sample, Box-Cox).
#  1. mean score (z within cohort) by sex in carriers, SSC relatives and SSC autistic probands (probands via ssc_pgs_scores + X of 24)
#  2. three-way CNV x BMI x score vs SSC relatives, (a) as in Supplementary Fig. 8, (b) score re-standardised within sex and cohort,
#     (c) (b) + sex x BMI, sex x score and sex x CNV x BMI terms; within-carrier BMI x score likewise. Output: pgs_sex_check.csv.
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s64c <- readLines(file.path(HERE, "64c_fig_supp_interaction_forest_age_asd.R")); eval(parse(text = s64c[1:grep("^SIBf", s64c)[1]]))
tm <- function(m, parts) { nm <- names(coef(m)); nm[sapply(strsplit(nm, ":"), function(z) setequal(z, parts))][1] }
out <- list(); add <- function(...) out[[length(out) + 1]] <<- data.table(...)
pg <- fread(file.path(DATA, "ssc_pgs_scores.tsv"))
sscid <- fread(file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9/SSC IDs Sanders.csv"))
pg[, role := fcase(grepl("[.]p[0-9]+$", individual), "SSC autistic probands", grepl("[.]s[0-9]+$", individual), "SSC siblings", grepl("[.](mo|fa)$", individual), "SSC parents", default = NA_character_)]
pg[, sex := fcase(grepl("[.]mo$", individual), "female", grepl("[.]fa$", individual), "male", default = NA_character_)]
sx <- unique(rbind(X[, .(id, sex)], C[, .(id, sex)]))[!duplicated(id)]
pg <- merge(pg, sx, by.x = "individual", by.y = "id", all.x = TRUE, suffixes = c("", ".x")); pg[is.na(sex), sex := sex.x]
for (r in c("SSC autistic probands", "SSC siblings", "SSC parents")) { d <- pg[role == r & !is.na(sex)]; if (nrow(d) < 50) next
  t <- t.test(pgs_asd ~ sex, d); add(set = r, n_female = sum(d$sex == "female"), n_male = sum(d$sex == "male"), mean_female = mean(d[sex == "female", pgs_asd]), mean_male = mean(d[sex == "male", pgs_asd]), p = t$p.value) }
car <- unique(X[is_cnv == 1, .(id, cohort, sex, pgs)])
for (k in c("all carriers", unique(car$cohort))) { d <- if (k == "all carriers") car else car[cohort == k]; if (min(table(d$sex)) < 10) next
  t <- t.test(pgs ~ sex, d); add(set = paste("carriers:", k), n_female = sum(d$sex == "female"), n_male = sum(d$sex == "male"), mean_female = mean(d[sex == "female", pgs]), mean_male = mean(d[sex == "male", pgs]), p = t$p.value) }
R1 <- rbindlist(out, fill = TRUE); out <- list()
X[, pgs_sx := as.numeric(scale(pgs)), by = .(cohort, sex, is_cnv)]
for (ph in c("SRS", "SCQ")) { X3 <- X[pheno == ph]; car3 <- X3[is_cnv == 1]; ref3 <- RELf(X3)
  z <- rbind(car3, ref3)[, .(y, BMI, pgs, pgs_sx, age, male, cg, is_cnv)]
  for (v in c("pgs", "pgs_sx")) {
    m1 <- lm(as.formula(sprintf("y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI * %s", v)), z); r1 <- rob(m1, tm(m1, c("is_cnv", "BMI", v)))
    m2 <- lm(as.formula(sprintf("y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI * %s + male * BMI * is_cnv + male * %s", v, v)), z); r2 <- rob(m2, tm(m2, c("is_cnv", "BMI", v)))
    w1 <- rob(lm(as.formula(sprintf("y ~ cg + poly(age, 2) + male + BMI * %s", v)), car3), tm(lm(as.formula(sprintf("y ~ cg + poly(age, 2) + male + BMI * %s", v)), car3), c("BMI", v)))
    w2m <- lm(as.formula(sprintf("y ~ cg + poly(age, 2) + male * BMI + male * %s + BMI * %s", v, v)), car3); w2 <- rob(w2m, tm(w2m, c("BMI", v)))
    lab <- if (v == "pgs") "score z within cohort (as Supp. Fig. 8)" else "score z within cohort and sex"
    add(trait = ph, score = lab, test = "three-way vs SSC relatives", b = r1[["b"]], se = r1[["se"]], p = r1[["p"]])
    add(trait = ph, score = lab, test = "three-way vs SSC relatives, + sex x BMI (x CNV) and sex x score", b = r2[["b"]], se = r2[["se"]], p = r2[["p"]])
    add(trait = ph, score = lab, test = "within carriers BMI x score", b = w1[["b"]], se = w1[["se"]], p = w1[["p"]])
    add(trait = ph, score = lab, test = "within carriers BMI x score, + sex x BMI and sex x score", b = w2[["b"]], se = w2[["se"]], p = w2[["p"]]) } }
R2 <- rbindlist(out); fwrite(rbind(R1, R2, fill = TRUE), file.path(HERE, "pgs_sex_check.csv")); options(width = 200)
print(R1[, .(set, n_female, n_male, mean_female = round(mean_female, 3), mean_male = round(mean_male, 3), p = signif(p, 2))]); print(R2[, .(trait, score, test, b = round(b, 3), se = round(se, 3), p = signif(p, 2))])
