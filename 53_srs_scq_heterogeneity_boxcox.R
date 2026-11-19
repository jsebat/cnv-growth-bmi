# Heterogeneity of the BMI slope on SRS / SCQ among the four CNVs, Box-Cox scale (JS / Claude 2026-11-14), same model as the
# Methods: within carriers y ~ s(age, by = locus) + locus + BMI + locus:BMI, Wald test of the three locus:BMI terms.
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
for (ph in c("SRS", "SCQ")) { x <- X[pheno == ph & grp %in% LOC][, locus := factor(grp)]
  m <- gam(y ~ s(age, by = locus, k = 6) + locus + BMI + locus:BMI, data = x, method = "REML"); a <- anova(m)$pTerms.table
  cat(ph, "heterogeneity locus:BMI P =", signif(a["locus:BMI", "p-value"], 2), "\n") }
