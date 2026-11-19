# Adult BMI scaling check (JS / Claude 2026-11-14): old rank-based vs new Box-Cox adult BMI z (as the growth panels), in the clinical
# non-carrier reference adults and in SSC parents. Output printed only.
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
sink(tempfile()); source(file.path(HERE, "21_pheno_bmi_gam.R")); sink()
r <- REFad[, .(sex = as.character(sex), age, BMI_raw)][, `:=`(z_rank = bmiZ_rank(BMI_raw, age, sex), z_bc = bmiZ_bc(BMI_raw, age, sex))]
print(r[, .(n = .N, mean_rank = round(mean(z_rank, na.rm = TRUE), 3), sd_rank = round(sd(z_rank, na.rm = TRUE), 3), mean_bc = round(mean(z_bc, na.rm = TRUE), 3), sd_bc = round(sd(z_bc, na.rm = TRUE), 3), cor = round(cor(z_rank, z_bc, use = "complete"), 4)), by = sex])
cat("children unchanged:", isTRUE(all.equal(bmiZ_rank(20, 10, "male"), bmiZ_bc(20, 10, "male"))), "\n")
