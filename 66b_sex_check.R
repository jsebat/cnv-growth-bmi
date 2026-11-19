# Check of the sex difference in the CNV x BMI interaction on SRS (JS / Claude 2026-11-14): BMI slopes by sex in carriers and in
# SSC relatives (siblings, parents), same bam age-smooth model per group, and the sex difference with sex-specific age terms.
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s64c <- readLines(file.path(HERE, "64c_fig_supp_interaction_forest_age_asd.R")); eval(parse(text = s64c[1:grep("^SIBf", s64c)[1]]))
out <- list()
for (ph in c("SRS", "SCQ")) { C2 <- C[pheno == ph]
  for (g in c("carriers", "SSC siblings", "SSC parents")) for (s in c("female", "male")) { d <- switch(g, carriers = C2[is_cnv == 1], `SSC siblings` = SIBf(C2), `SSC parents` = PARf(C2))[sex == s]
    if (nrow(d) < 30) next; r <- slope(d); out[[length(out) + 1]] <- data.table(trait = ph, group = g, sex = s, n = r[["n"]], slope = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  for (rn in c("siblings", "parents")) { ref <- if (rn == "siblings") SIBf(C2) else PARf(C2); if (nrow(ref) < 50) next
    for (s in c("female", "male")) { r <- vs_ref(C2[is_cnv == 1 & sex == s & (if (rn == "parents") age >= 18 else age < 18)], ref[sex == s])
      out[[length(out) + 1]] <- data.table(trait = ph, group = paste("CNV x BMI vs SSC", rn, "(age-matched)"), sex = s, n = r[["n"]], slope = r[["b"]], se = r[["se"]], p = r[["p"]]) } }
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "sex_check_bmi_slopes.csv")); options(width = 200); print(R[, .(trait, group, sex, n, slope = round(slope, 3), se = round(se, 3), p = signif(p, 2))])
