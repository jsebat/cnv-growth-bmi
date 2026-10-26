# Build resources/ssc_pgs_scores.tsv (JS / Claude 2026-10-26): SSC polygenic scores for childhood obesity, childhood BMI and
# adult BMI, scored on Expanse (jsebat/pgs_obesity_ssc/scores/ssc.<trait>.sscore; plink2 --score on the SSC genotype set with
# Anjali's SBayesRC weights, the same weights used for the G2MH WGS cohort), residualised on the 10 ancestry PCs of the SSC
# samples in all_cohorts_adsi_pgs.tsv and z-scored within SSC, mapped from the SS genotype id to the SSC individual id
# (SSC IDs Sanders.csv, SS_RUID / RUID-1). Input .sscore files are copied to resources/ssc_scores_raw/ by the shell step.
suppressPackageStartupMessages(library(data.table))
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")   # participant-level tables stay in Jon/Claude
RAW <- file.path(DATA, "ssc_scores_raw")
rd <- function(nm) { x <- fread(file.path(RAW, sprintf("ssc.%s.sscore", nm))); setnames(x, "#IID", "IID"); x[, .(IID, score = SCORE1_AVG)] }
# height: GIANT 2022 SBayesRC weights (the file used for the G2MH WGS cohort); autism: April 2026 ASD GWAS leave-SSC-out weights
sc <- Reduce(function(a, b) merge(a, b, by = "IID"), list(setnames(rd("childhood_obesity"), "score", "cobes"), setnames(rd("childhood_bmi"), "score", "cbmi"), setnames(rd("bmi_2018"), "score", "bmi"),
                                                          setnames(rd("height_2022"), "score", "ht"), setnames(rd("asd2026_loSSC"), "score", "asd")))
pc <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/all_cohorts_adsi_pgs.tsv"), select = c("sample_id", "source", paste0("PC", 1:10)))[source == "ssc"]
sc <- merge(sc, pc[, c("sample_id", paste0("PC", 1:10)), with = FALSE], by.x = "IID", by.y = "sample_id")
cat("SSC scored:", nrow(sc), "with PCs\n")
resid_z <- function(v) { f <- as.formula(paste("v ~", paste0("PC", 1:10, collapse = " + "))); as.numeric(scale(residuals(lm(f, data = cbind(sc, v = v))))) }
sc[, `:=`(pgs_cobes = resid_z(cobes), pgs_cbmi = resid_z(cbmi), pgs_bmi = resid_z(bmi), pgs_ht = resid_z(ht), pgs_asd = resid_z(asd))]
san <- fread(file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9/SSC IDs Sanders.csv"))
map <- unique(rbind(san[, .(IID = SS_RUID, individual)], san[, .(IID = `RUID-1`, individual)])[IID != "." & IID != ""])
out <- merge(sc[, .(IID, pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht, pgs_asd)], map, by = "IID")[, .SD[1], by = individual]
cat("mapped to SSC individuals:", nrow(out), "\n")
cat("correlations among the scores:\n"); print(round(cor(out[, .(pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht, pgs_asd)]), 2))
fwrite(out, file.path(DATA, "ssc_pgs_scores.tsv"), sep = "\t")
