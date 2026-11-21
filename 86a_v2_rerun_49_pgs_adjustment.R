# 86a_v2_rerun_49_pgs_adjustment.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 3)
# What it does: re-runs 49_fig_supp_bmi_srs_scq_pgs.R (G2MH carriers with ten polygenic scores: BMI slope under score adjustment,
#   score -> trait reduced forms, Mendelian randomisation with the adiposity scores, BMI x score interactions) from its own code
#   (49 -> 48 -> 47 -> mr_reduced_form.R -> 23 / 24), unchanged, with its outputs redirected to autism/aux_rerun, and compares the
#   result with the saved table fig_supp_bmi_srs_scq_pgs_boxcox_table.csv. ADD_CHOP is FALSE here, as when 49 was first run: the
#   only growth input of 49 is the carriers' own measured height, and the CHOP overlap rule would remove the G2MH CHOP-site rows.
# Inputs : as 49. Outputs: aux_rerun/fig_supp_bmi_srs_scq_pgs_boxcox{.png,.pdf,_table.csv}; participant-level cache
#   DATA/v2_autism_cache/pgs49_2026-11-20.RData (W = G2MH carriers with ten scores; tables of 49) for 86 / 86b.
# null model: none; Wald tests with HC1 robust SE (ordinary and two-stage least squares), as in 49.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
Sys.setenv(ADD_CHOP = "FALSE", SRSSCQ_SCALE = "boxcox")
v2_source(file.path(HERE, "49_fig_supp_bmi_srs_scq_pgs.R"))
old <- fread(file.path(HERE, "fig_supp_bmi_srs_scq_pgs_boxcox_table.csv")); new <- fread(file.path(V2AUX, "fig_supp_bmi_srs_scq_pgs_boxcox_table.csv"))
cat("49 re-run vs saved table: rows", nrow(old), nrow(new), " same n:", all(old$n == new$n), " max |b diff|:", max(abs(old$b - new$b)), " max |p diff|:", max(abs(old$p - new$p)), "\n")
W49 <- copy(W); A49 <- copy(A); B49 <- copy(B); C49 <- copy(C); D49 <- copy(Dd); tab49 <- copy(tab); g2all <- Reduce(function(a, b) merge(a, b, by = "GUID", all = TRUE), list(g2, ex[, .(GUID, pgs_aam, pgs_hci)], nd[, .(GUID, pgs_adhd, pgs_edu, pgs_iq)]))
DATA0 <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
save(W49, A49, B49, C49, D49, tab49, g2all, LAB, ADI, GRO, ND, TEN, BASE, file = file.path(DATA0, "v2_autism_cache", "pgs49_2026-11-20.RData"))
print(names(W49)); print(names(g2all)); cat("86a done\n")
