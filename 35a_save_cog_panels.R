# Save the SRS/SCQ panel machinery of 24_cognition_bmi_slopes.R (data X and R, panel functions and their constants) so that
# the six-column main figure with UK Biobank adults (35_main_figure_v5_ukbb.R) can draw the SRS/SCQ columns on Expanse,
# where the SSC and G2MH phenotype files are not available. JS 2026-10-27. Output: cog_panels.RData (participant-level
# rows without identifiers: id column dropped), copied to the Expanse working directory by hand.
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))   # 24 needs D (and helpers) from 23
source(file.path(HERE, "24_cognition_bmi_slopes.R"))
X <- data.table::copy(X); if ("id" %in% names(X)) X[, id := NULL]
objs <- c("X", "R", "panel_ctl", "panel_cnv", "fitline", "bins", "ylim_of", "ptxt", "REFC", "COLS3", "XL", "YL")
objs <- objs[sapply(objs, exists)]
cat("saving:", objs, "\n")
save(list = objs, file = file.path(HERE, "cog_panels.RData"))
cat("cog_panels.RData written; X rows:", nrow(X), " phenos:", paste(unique(X$pheno), collapse = ", "), "\n")
