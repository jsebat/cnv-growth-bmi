# Why does our analysis not reproduce the nominal ASD - head circumference association in 22q11.2 dup?
# Reference: Anjali's asd_stratified_trajectories.R (violin model) on genomatrix_clinical_expanded.csv
#   (asd_dx from add_asd_dx.R: G2MH/Cardiff = Feb-2026 freeze corepsych_asd_dx, 2 = subthreshold -> 0; SSC .p = 1,
#   SSC parents/siblings = 0), latest visit, age <= 45, Cardiff relabelled G2MH and NOT de-duplicated,
#   lm(HEADCIRCUM ~ ASD + age + age^2).
# Starting from her exact spec, change one thing at a time towards ours (supplementary panel C / main panel J):
#   ASD source (our G2MH subjects.csv psy_asd), Cardiff/G2MH duplicates dropped, age 2-40, + cohort + adult terms.
suppressPackageStartupMessages({ library(data.table) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
g <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/trajectory_data/genomatrix_clinical_expanded.csv"))
g <- g[last_timepoint == TRUE][, CNV := gsub("proximal", "BP4-BP5", trimws(CNV))]
g[, cohort_orig := cohort][, cohort := ifelse(cohort == "Cardiff", "G2MH", cohort)]
s <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd")); ours <- setNames(as.integer(s$asd), s$sfari_id)
g[, asd_ours := asd_dx][cohort_orig %in% c("G2MH", "Cardiff"), asd_ours := ours[GUID]]
g[cohort_orig == "SSC" & !grepl("[.]p\\d$", sampleID), asd_ours := NA_integer_]     # we left SSC parents/sibs out

fit <- function(x, asd_col, lab, adj = FALSE) {
  x <- x[!is.na(get(asd_col)) & !is.na(HEADCIRCUM)]; x[, a := get(asd_col)]
  f <- if (adj) HEADCIRCUM ~ a + cohort_orig + I(age >= 18) + age + I(age^2) else HEADCIRCUM ~ a + age + I(age^2)
  co <- summary(lm(f, data = x))$coefficients["a", ]
  data.table(variant = lab, n_ASD = sum(x$a == 1), n_noASD = sum(x$a == 0), beta = round(co[1], 3), p = signif(co[4], 3))
}
d <- g[CNV == "22q11.2 A-D DUP"]
dd <- d[!(cohort_orig == "G2MH" & Cardiff_mask == TRUE)]
res <- rbind(
  fit(d[age <= 45], "asd_dx", "1 Anjali exact (her asd_dx, dups kept, age <= 45)"),
  fit(dd[age <= 45], "asd_dx", "2 + Cardiff/G2MH duplicates dropped"),
  fit(d[age <= 45], "asd_ours", "3 her spec but our ASD source (psy_asd)"),
  fit(dd[age >= 2 & age <= 40], "asd_dx", "4 her ASD, dups dropped, age 2-40"),
  fit(dd[age >= 2 & age <= 40], "asd_dx", "5 = 4 + cohort + adult terms", adj = TRUE),
  fit(dd[age >= 2 & age <= 40], "asd_ours", "6 our ASD, dups dropped, age 2-40"),
  fit(dd[age >= 2 & age <= 40], "asd_ours", "7 = 6 + cohort + adult terms (ours)", adj = TRUE))
print(res); fwrite(res, file.path(HERE, "reconcile_22qdup_hc_asd.csv"))

# where do the two ASD sources disagree for 22q11.2 dup carriers with HC?
x <- dd[!is.na(HEADCIRCUM) & cohort_orig %in% c("G2MH", "Cardiff")]
cat("\n22q11.2 dup with HC (G2MH/Cardiff): her asd_dx x our psy_asd\n"); print(table(hers = x$asd_dx, ours = x$asd_ours, useNA = "ifany"))
cat("\nadults vs children, her ASD coding (dups dropped, age <= 45):\n")
print(dd[!is.na(HEADCIRCUM) & !is.na(asd_dx) & age <= 45, .(n = .N, n_asd = sum(asd_dx), mean_HC = round(mean(HEADCIRCUM), 2)), by = .(adult = age >= 18, asd_dx)][order(adult, asd_dx)])
