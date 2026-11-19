#!/bin/bash
# 2026-11-14: rerun of the SRS/SCQ pipeline after adult BMI moved to Box-Cox (21_pheno_bmi_gam.R bmiZ_bc). Logs in rerun_logs/.
cd "$(dirname "$0")"; mkdir -p rerun_logs; R=$(ls /c/Program\ Files/R/*/bin/Rscript.exe | tail -1)
for s in 60_combined_srs_scq_build.R 61_combined_cnv_x_bmi.R 62_combined_threeway.R 63_combined_noasd_sensitivity.R 66_fig_supp_interactions_main.R 66b_sex_check.R 67_ssc_sex_x_bmi.R 68_pgs_sex_check.R 51b_cog_panels_boxcox_combined.R 49_fig_supp_bmi_srs_scq_pgs.R 30b_sex_pgs_sensitivity_boxcox.R 64_fig_supp_interaction_forest.R 64b_fig_supp_interaction_forest_sibpar.R 64c_fig_supp_interaction_forest_age_asd.R; do
  "$R" "$s" > "rerun_logs/${s%.R}.log" 2>&1; echo "$s rc $?"; done
