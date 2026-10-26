# Build the analysis table for CNV growth trajectories (cross-sectional + longitudinal).
#
# Source: Molly's clinical genomatrix (G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/
#   genomatrix.Clinical.02092026.csv, built by CreateGenomatrix2.R) -- the table behind Anjali's
#   trajectory figures. Cohorts G2MH, Cardiff, SVIP, SSC, Searchlight; SSC/SVIP non-carriers = controls.
#   Z-scores as supplied: < 18 y childsds (height WHO / WHO2007, BMI CDC <= 5 y / WHO2007, HC UK90);
#   >= 18 y within-sample age/age^2/sex residual -> Box-Cox -> scale.
# Fixes applied here:
#   * 71 G2MH rows with Cardiff_mask == TRUE are Cardiff participants entered twice -> dropped
#     (Anjali's compute_interaction_results.R keeps both copies).
#   * visit order taken from age, not `timepoint` (Cardiff timepoints are not chronological).
#   * Cardiff stimulant use per visit (column `ndri` = methylphenidate/amfetamines in Cardiff's coding)
#     merged from CardiffLongitudinalWeightandHeightDataFinalwithdrugs_2025-07-02.xlsx by GUID + age.
# Output: analysis_table.rds (one row per visit) -- written to DATA (Jon/Claude/CNV_growth_trajectories_data), not to the
#   manuscript folder, because it is participant-level phenotype data.

suppressPackageStartupMessages({ library(data.table); library(readxl) })
DB   <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")   # participant-level tables stay in Jon/Claude (controlled data)

g <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
g <- g[!(cohort == "G2MH" & Cardiff_mask == TRUE)]
g[, CNV := trimws(CNV)]
lab <- c("16p11.2 proximal DEL" = "16p11.2 del", "16p11.2 proximal DUP" = "16p11.2 dup",
         "22q11.2 A-D DEL" = "22q11.2 del", "22q11.2 A-D DUP" = "22q11.2 dup", "No CNV" = "No CNV")
g <- g[CNV %in% names(lab)]
g[, locus := factor(lab[CNV], levels = lab)]
g[, id := paste(cohort, MSassigned_ID, sep = ":")]
g[, sex := factor(ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA)))]
# JS 2026-10-24: (1) SSC parent ages in the genomatrix are a 2012-09 snapshot, ~3 y older than the visit; replace by visit age
# (proband age at evaluation + birth-month difference; see README_ages_and_line_endings.md in the SSC folder). (2) Plausibility
# QC of height / BMI / head circumference (00_qc_growth.R; manual corrections from implausible_height_weight_corrections.csv).
SSC <- file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9")
agef <- fread(file.path(SSC, "Age_individual_v14_parent_age.txt")); setnames(agef, c("id()", "individual_properties.role"), c("sid", "role0"))
agef[, bm := as.integer(substr(birth, 1, 4)) * 12L + as.integer(substr(birth, 6, 7))]
ssc_ev <- fread(file.path(SSC, "ssc_case_HC.csv"))[, .(proband = Proband_ID, age_eval = `Age at Evaluation (Non-PHI) :: Proband SSC Height, Weight, Head Circumference`)]
ssc_pro <- merge(agef[role0 == "proband", .(family, pbm = bm, proband = sid)], ssc_ev, by = "proband")
parage <- merge(agef[role0 %in% c("mother", "father"), .(sid, family, obm = bm)], ssc_pro, by = "family")[, .(sid, age = (age_eval + pbm - obm) / 12)][!is.na(age)]
n_par <- sum(g$cohort == "SSC" & g$sampleID %in% parage$sid)
g[cohort == "SSC" & sampleID %in% parage$sid, age := parage[match(sampleID, sid), age]]
cat("SSC parent ages replaced by visit age:", n_par, "rows\n")
source(file.path(HERE, "00_scoring.R")); source(file.path(HERE, "00_qc_growth.R"))
g <- qc_growth(g, label = "genomatrix (01_build)")
setorder(g, id, age)
g[, visit := seq_len(.N), by = id]; g[, n_visits := .N, by = id]
g[, latest := visit == n_visits]

# Cardiff stimulant flag
tmp <- file.path(tempdir(), "cardiff"); dir.create(tmp, showWarnings = FALSE)
unzip(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/Cardiff22qData.zip"), exdir = tmp)
cf <- as.data.table(read_excel(list.files(tmp, "CardiffLongitudinal.*xlsx", recursive = TRUE, full.names = TRUE)[1], sheet = 1))
cf <- unique(cf[, .(GUID = guid, age = as.numeric(Age), stimulant = ndri == "yes", antipsychotic = antipsychotic == "yes")])
cf <- cf[, .(stimulant = any(stimulant), antipsychotic = any(antipsychotic)), by = .(GUID, age)]
g <- merge(g, cf, by = c("GUID", "age"), all.x = TRUE)
g[cohort != "Cardiff", `:=`(stimulant = NA, antipsychotic = NA)]
setorder(g, id, age)

keep <- c("id", "cohort", "locus", "sex", "age", "visit", "n_visits", "latest", "HEIGHT", "BMI", "HEADCIRCUM",
          "HEIGHT_raw", "BMI_raw", "HEADCIRCUM_raw", "stimulant", "antipsychotic")
tab <- g[, ..keep]
saveRDS(tab, file.path(DATA, "analysis_table.rds"))

cat("rows", nrow(tab), " persons", uniqueN(tab$id), "\n")
print(dcast(tab[latest == TRUE], locus ~ cohort, fun.aggregate = length, value.var = "id"))
cat("\nCardiff visits with stimulant flag matched:", tab[cohort == "Cardiff", sum(!is.na(stimulant))], "of", tab[cohort == "Cardiff", .N],
    "; stimulant = yes:", tab[cohort == "Cardiff", sum(stimulant, na.rm = TRUE)], "\n")
lg <- tab[locus != "No CNV" & age < 18 & !is.na(HEIGHT) & !is.na(BMI), .N, by = .(id, cohort, locus)][N >= 2]
cat("\ncarriers < 18 y with >= 2 height+BMI visits:\n"); print(dcast(lg, locus ~ cohort, fun.aggregate = length, value.var = "id"))
