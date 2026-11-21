# 83b_v2_verify_combined_build.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C)
# What it does: checks that the SRS / SCQ analysis table DATA/combined_srs_scq.rds (written 2026-11-14 by
#   60_combined_srs_scq_build.R) is what the current code and source files produce: 60 is re-run from its own code in memory (its
#   saveRDS / fwrite lines are dropped, nothing is written to DATA) and the result is compared row by row with the saved table.
#   ADD_CHOP is FALSE as when 60 was run (CHOP has no SRS / SCQ; the CHOP overlap rule only concerns the growth tables).
# Inputs : as 60. Output: autism/aux_rerun/verify_combined_build.txt (aggregate comparison only).
# null model: none; a data check, no test.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
Sys.setenv(ADD_CHOP = "FALSE")
L <- v2_lines(file.path(HERE, "60_combined_srs_scq_build.R")); L <- L[!grepl("^saveRDS\\(|fwrite\\(smp", L)]
L <- sub("; fwrite\\(smp.*$", "", L); v2_eval(L)
old <- readRDS(file.path(DATA, "combined_srs_scq.rds")); k <- c("id", "pheno", "cohort", "locus")
m <- merge(old, C, by = k, all = TRUE, suffixes = c(".saved", ".rerun"))
txt <- c(sprintf("rows saved %d, re-run %d, matched on id x trait x cohort x CNV %d", nrow(old), nrow(C), sum(!is.na(m$y.saved) & !is.na(m$y.rerun))),
         sprintf("max |difference| in the trait score y: %.3g; in BMI Z: %.3g; in age: %.3g", max(abs(m$y.saved - m$y.rerun), na.rm = TRUE), max(abs(m$BMI.saved - m$BMI.rerun), na.rm = TRUE), max(abs(m$age.saved - m$age.rerun), na.rm = TRUE)),
         "rows only in one of the two tables, by cohort and CNV:", capture.output(print(m[is.na(y.saved) | is.na(y.rerun), .N, by = .(pheno, cohort, locus, only_in = fifelse(is.na(y.saved), "re-run", "saved"))])))
writeLines(txt, file.path(V2AUX, "verify_combined_build.txt")); cat(txt, sep = "\n"); cat("83b done\n")
