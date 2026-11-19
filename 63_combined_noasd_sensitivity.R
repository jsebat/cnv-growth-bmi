# Autism-ascertainment sensitivity for the combined CNV x BMI and three-way analyses (JS / Claude 2026-11-14; Jon: SVIP and
# Searchlight tend to recruit from autism cohorts, which could skew the effects; exclude carriers with an autism diagnosis).
# Autism diagnosis: G2MH subjects.csv (as in 21/23); SVIP diagnosis_summary.clinical_asd_dx; Searchlight v14
# summary_data/searchlight_core_descriptive_variables.csv (autism if clinical_asd = 1 or reported_asd = 1; 888 = not reported;
# unknown if neither is 0/1); SSC carriers: probands = autism, others not. SSC relatives (reference) are unaffected by definition.
# Re-runs (a) CNV x BMI vs SSC relatives (61 models) and (b) the within-carrier BMI x autism score interaction (62 models) in
# carriers without an autism diagnosis, and in carriers with one, combined and per cohort. Output: combined_noasd_sensitivity.csv.
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s61 <- readLines(file.path(HERE, "61_combined_cnv_x_bmi.R")); eval(parse(text = s61[1:(grep("^out <- list", s61) - 1)]))
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
sv <- fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"))[, .(id = individual, a_sv = as.integer(diagnosis_summary.clinical_asd_dx))]
sl <- rbindlist(lapply(c("16p11.2 deletion", "16p11.2 duplication"), function(cn) fread(file.path(DB, "G2MH/Simons_Searchlight/Simons_Searchlight_Dataset_v14.0", cn, "summary_data/searchlight_core_descriptive_variables.csv"))[, .(sfari_id, clinical_asd, reported_asd)]))
sl <- unique(sl)[, a_sl := fifelse(clinical_asd %in% 1L | reported_asd %in% 1L, 1L, fifelse(clinical_asd %in% 0L | reported_asd %in% 0L, 0L, NA_integer_))][, id := sub("-", ".", sfari_id, fixed = TRUE)][, .(a_sl = max(a_sl, na.rm = FALSE)), by = id]
C <- merge(merge(C, sv, by = "id", all.x = TRUE), sl, by = "id", all.x = TRUE)
C[, asd2 := fcase(cohort == "G2MH", as.integer(asd), cohort == "SVIP", fifelse(!is.na(a_sv), a_sv, a_sl), cohort == "Searchlight", fifelse(!is.na(a_sl), a_sl, a_sv),
                  cohort == "SSC" & is_cnv == 1, as.integer(grepl("[.]p[0-9]$", id)), locus == "SSC relative", 0L, default = NA_integer_)]
cat("autism diagnosis among carriers, by cohort (SRS rows):\n"); print(C[pheno == "SRS" & is_cnv == 1, .(n = .N, autism = sum(asd2 %in% 1L), no_autism = sum(asd2 %in% 0L), unknown = sum(is.na(asd2)), pct_autism_known = round(100 * mean(asd2[!is.na(asd2)]), 1)), by = cohort])
out <- list(); add <- function(ph, sub, coh, test, r) out[[length(out) + 1]] <<- data.table(trait = ph, subset = sub, cohort = coh, test = test, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
for (ph in c("SRS", "SCQ")) {
  ref <- C[pheno == ph & locus == "SSC relative"]
  for (sub in c("all carriers", "no autism diagnosis", "autism diagnosis")) for (coh in c("Combined", "G2MH", "SVIP", "Searchlight")) {
    xc <- C[pheno == ph & is_cnv == 1 & (coh == "Combined" | cohort == coh)]
    xc <- switch(sub, `all carriers` = xc, `no autism diagnosis` = xc[asd2 %in% 0L], `autism diagnosis` = xc[asd2 %in% 1L])
    if (nrow(xc) < 15) next
    add(ph, sub, coh, "BMI slope in carriers", slope(xc)); add(ph, sub, coh, "CNV x BMI vs SSC relatives", vs_ref(xc, ref))
  }
  for (L in c("22q11.2 del", "16p11.2 del", "16p11.2 dup")) { xc <- C[pheno == ph & is_cnv == 1 & locus == L & asd2 %in% 0L]; if (nrow(xc) >= 15) add(ph, "no autism diagnosis", paste("Combined,", L), "CNV x BMI vs SSC relatives", vs_ref(xc, ref)) }
}
# three-way (within-carrier BMI x autism score) in carriers without an autism diagnosis, 62 data
AN <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data"); s62 <- readLines(file.path(HERE, "62_combined_threeway.R")); eval(parse(text = s62[grep("^rob <- function", s62):grep("^X\\[, cg :=", s62)]))
for (ph in c("SRS", "SCQ")) for (sub in c("all carriers", "no autism diagnosis", "autism diagnosis")) {
  car <- X[pheno == ph & is_cnv == 1]; car <- switch(sub, `all carriers` = car, `no autism diagnosis` = car[asd2 %in% 0L], `autism diagnosis` = car[asd2 %in% 1L])
  r <- rob(lm(y ~ cg + poly(age, 2) + male + BMI * pgs + cohort:BMI + cohort:pgs, car), "BMI:pgs"); add(ph, sub, "Combined (G2MH + SVIP + SSC)", "BMI x autism score within carriers", c(n = nrow(car), r))
  for (coh in c("G2MH", "SVIP")) { d <- car[cohort == coh]; if (nrow(d) >= 15) add(ph, sub, coh, "BMI x autism score within carriers", c(n = nrow(d), rob(lm(as.formula(paste("y ~", if (uniqueN(d$cg) > 1) "cg +" else "", "poly(age, 2) + male + BMI * pgs")), d), "BMI:pgs"))) }
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "combined_noasd_sensitivity.csv"))
options(width = 220); print(R[, .(trait, subset, cohort, test, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
