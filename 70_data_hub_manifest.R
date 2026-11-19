# Data hub manifest for the CNV x BMI / autistic-trait analyses of Fig. 5 and Supplementary Figs. 7-8 (JS / Claude 2026-11-14;
# Jon: all data in one centralised, easy-to-find place). Lists every Dropbox input the pipeline reads (exact files), the derived
# participant-level tables, the code and the aggregate results, with size and md5; writes hub_filelist_*.txt (paths relative to the
# Dropbox root, for rsync --files-from) and DATA_MANIFEST.tsv. The hub itself is on Expanse:
# /expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits (inputs/, derived/, code/, results/, expanse_links/).
suppressPackageStartupMessages(library(data.table)); library(tools)
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
SSCr <- "Autism/SSC_Version_15_Phenotype_Data_Set_9"; SVr <- "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0"
SLr <- "G2MH/Simons_Searchlight/Simons_Searchlight_Dataset_v14.0"; ANr <- "G2MH/NetworkPaper/AnalysisByAnjali/data"
inputs <- c(
  "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv", "G2MH/DataFreeze/Feb2026/g2mh_project2_medications_data_release_022026.csv",
  "G2MH/DataFreeze/Feb2026/g2mh_project2_data_dictionary_022026.xlsx",
  "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv",
  "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/Searchlight.WGS2.2023_05.sample_metadata.tsv",
  "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/CreateGenomatrix2.R",
  file.path(ANr, "pca_results/naive/cohort.sscore"), file.path(ANr, "all_cohorts_adsi_pgs.tsv"), file.path(ANr, "searchlight/pgs/analysis_dataset.tsv"),
  file.path(ANr, "searchlight/pgs/analysis_dataset_dictionary.tsv"), file.path(ANr, "pgs_scores/pgs_data_dictionary.xlsx"),
  list.files(file.path(DB, ANr, "pgs_scores/scores_wgs"), pattern = "[.]sscore$", full.names = FALSE) |> (\(f) file.path(ANr, "pgs_scores/scores_wgs", f))(),
  file.path(SVr, c("srs_parent.csv", "srs_adult.csv", "scq_life.csv", "svip_subjects.csv", "diagnosis_summary.csv")),
  file.path(SLr, "all_cohorts/measure_data", c("srs_sa.csv", "srs_adult.csv", "scq_life.csv", "medical_history_version_3.1.csv", "medical_history_version_3.2.csv")),
  file.path(SLr, "all_cohorts/summary_data/searchlight_core_descriptive_variables.csv"), file.path(SLr, "all_cohorts/carrier_count_by_measure.csv"),
  unlist(lapply(c("16p11.2 deletion", "16p11.2 duplication"), function(cn) c(file.path(SLr, cn, "measure_data", c("srs_sa.csv", "srs_adult.csv", "scq_life.csv")),
    file.path(SLr, cn, "summary_data/searchlight_core_descriptive_variables.csv")))),
  file.path(SSCr, c("Proband Data/srs_parent.csv", "Proband Data/scq_life.csv", "Proband Data/rbs_r_raw.csv", "Proband Data/ssc_core_descriptive.csv",
    "Designated Unaffected Sibling Data/srs_parent.csv", "Designated Unaffected Sibling Data/srs_adult.csv", "Designated Unaffected Sibling Data/scq_life.csv",
    "Designated Unaffected Sibling Data/ssc_commonly_used.csv", "Mother Data/srs_adult.csv", "Father Data/srs_adult.csv", "Other Sibling Data/scq_life.csv",
    "Age_individual_v14_parent_age.txt", "ssc_case_HC.csv", "SSC IDs Sanders.csv")),
  "Jon/Claude/G2MH_CNV_growth/subjects.csv")
derived <- file.path("Jon/Claude/CNV_growth_trajectories_data", list.files(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data"), recursive = TRUE))
code <- file.path("Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories",
  list.files(HERE, pattern = "[.](R|sh|py|md|sbatch)$"))
res <- file.path("Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories",
  list.files(HERE, pattern = "[.](csv|png|pdf|xlsx|RData)$"))
M <- rbind(data.table(set = "inputs", rel = inputs), data.table(set = "derived", rel = derived), data.table(set = "code", rel = code), data.table(set = "results", rel = res))
M[, path := file.path(DB, rel)][, exists := file.exists(path)][, bytes := file.size(path)][exists == TRUE, md5 := as.character(md5sum(path))]
cat("missing:\n"); print(M[exists == FALSE, rel])
print(M[, .(files = .N, missing = sum(!exists), MB = round(sum(bytes, na.rm = TRUE) / 1e6, 1)), by = set])
fwrite(M[, .(set, file = rel, bytes, md5)], file.path(HERE, "DATA_MANIFEST.tsv"), sep = "\t")
for (s in unique(M$set)) writeLines(M[set == s & exists == TRUE, rel], file.path(HERE, sprintf("hub_filelist_%s.txt", s)))
cat("written DATA_MANIFEST.tsv and hub_filelist_*.txt\n")
