# run_check_lab.R  (2026-11-20)
# What: validation of the federated package against the lab's own v2 results. Runs federated/site_export.R on the lab's
#   clinical table (read only, in memory; UK Biobank adults are not available locally, so height is clinical only) and
#   compares the aggregate output with the paper run:
#     growth_unaffctl/main/main_with_ukbb_bmi_slope_by_band.csv, pooled_band_slopes.csv, band_slopes_family_clustered.csv,
#     cnv_x_bmi_all_ages.csv, sex_cnv_bmi_growth.csv, main_with_ukbb_bmi_genotype_by_band.csv, and longitudinal/results_*.csv.
#   Three exports: (1) export_fixed_ukbb = the package default (adult Box-Cox exponents fixed at the UK Biobank values);
#   (2) export_estimate = adult exponents estimated in the table's own non-carrier adults, which is what
#   17_main_with_ukbb.R does for BMI in its head-circumference table (no UK Biobank rows there);
#   (3) export_hc_longitudinal_inputs = the head-circumference visit sources of 80_longitudinal_v2.R (hc_visits.rds for
#   Simons VIP and Cardiff + the CHOP workbook), longitudinal section only.
# Inputs: DATA/analysis_table_chop_unaffctl_2026-11-20.rds, DATA/hc_visits.rds, CHOP workbook (all controlled, read only).
# Outputs (this folder): export_*/ (aggregate CSVs), comparison_with_paper_run.csv, comparison_summary.csv.
# null model: none; this script compares estimates.
suppressPackageStartupMessages(library(data.table))
ROOT <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
OUT <- file.path(ROOT, "run_v2_2026-11-20", "federated_check_lab"); V2 <- file.path(ROOT, "run_v2_2026-11-20")
FEDERATED_SOURCE_ONLY <- TRUE; FEDERATED_DIR <- file.path(ROOT, "federated"); source(file.path(FEDERATED_DIR, "site_export.R"))
RERUN <- Sys.getenv("RERUN", "TRUE") == "TRUE"
wide <- list(height_cm = c(0, 1000), weight_kg = c(0, 1000), bmi = c(0, 1000), hc_cm = c(0, 1000))   # the lab table is already cleaned (00_qc_growth.R)
base <- list(site_label = "Sebat lab clinical cohorts", columns = list(id = "id", family = "family", sex = "sex", age = "age", height_cm = "HEIGHT_raw", weight_kg = NULL, bmi = "BMI_raw",
             hc_cm = "HEADCIRCUM_raw", genotype = "locus", asd = "asd", cohort = "cohort"), genotype_map = list(control = "No CNV"), plausible = wide)
if (RERUN) {
  tab <- readRDS(file.path(DATA, "analysis_table_chop_unaffctl_2026-11-20.rds"))
  cat("lab table groups (rows); only 'No CNV' and the four CNVs are mapped, as in 82_growth_blocks_v2.R:\n"); print(tab[, .N, by = locus])
  run_site_export(modifyList(base, list(out_dir = file.path(OUT, "export_fixed_ukbb"))), data = tab)
  run_site_export(modifyList(base, list(out_dir = file.path(OUT, "export_estimate"), adult_boxcox = "estimate", sections = c("sample", "cross_sectional", "sex"))), data = tab)
  # head-circumference visit sources of 80_longitudinal_v2.R
  source(file.path(ROOT, "00_add_chop.R"))
  hv <- as.data.table(readRDS(file.path(DATA, "hc_visits.rds")))[cohort %in% c("SVIP", "Cardiff") & locus %in% LOCI, .(id = as.character(id), cohort, locus, sex, age, HEADCIRCUM_raw = hc, BMI_raw = bmi)]
  cw <- as.data.table(chop_rows())[, .(id = paste("CHOP", sampleID, sep = ":"), cohort = "CHOP", locus = ifelse(CNV == "22q11.2 A-D DEL", "22q11.2 del", "22q11.2 dup"),
                                       sex = ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA_character_)), age, HEADCIRCUM_raw, BMI_raw)][!is.na(HEADCIRCUM_raw)]
  hin <- rbind(hv, cw)[, HEIGHT_raw := NA_real_]
  cfg <- modifyList(base, list(out_dir = file.path(OUT, "export_hc_longitudinal_inputs"), sections = "longitudinal")); cfg$columns$family <- NULL; cfg$columns$asd <- NULL
  run_site_export(cfg, data = hin)
}
rdE <- function(d) fread(file.path(OUT, d, "estimates.csv"), na.strings = c("", "NA"), colClasses = list(character = c("n_carriers", "n_controls", "n_pairs")))
EF <- rdE("export_fixed_ukbb"); EE <- rdE("export_estimate"); EH <- rdE("export_hc_longitudinal_inputs")
G <- file.path(V2, "growth_unaffctl", "main"); cmp <- list()
add <- function(block, key, paper, n_paper, fed, n_fed, fed_alt = NA_real_) cmp[[length(cmp) + 1]] <<- data.table(block = block, item = key, paper = paper, n_paper = as.character(n_paper), federated = fed, n_federated = as.character(n_fed), federated_estimate_mode = fed_alt)
pick <- function(E, an, o, g, st, tm, sx = "all", mod = NULL) { r <- E[analysis == an & outcome == o & group == g & stratum == st & term == tm & sex == sx]; if (!is.null(mod)) r <- r[grepl(mod, model, fixed = TRUE)]
  if (nrow(r) != 1) return(list(estimate = NA_real_, se = NA_real_, n_carriers = NA, n_controls = NA, se_family_clustered = NA_real_)); as.list(r) }
# 1. band slopes per genotype (band_slopes of 17)
b <- fread(file.path(G, "main_with_ukbb_bmi_slope_by_band.csv"))
for (i in seq_len(nrow(b))) { r <- b[i]
  if (r$group == "No CNV") { f <- pick(EF, "cross-sectional band", r$outcome, "No CNV", r$age_bin, "BMI slope (controls-only model)"); e <- pick(EE, "cross-sectional band", r$outcome, "No CNV", r$age_bin, "BMI slope (controls-only model)")
    add("band: control BMI slope", paste(r$outcome, r$group, r$age_bin, sep = " | "), r$slope_controls, r$n_controls, f$estimate, f$n_controls, e$estimate); next }
  for (tm in c("CNV x BMI", "BMI slope, carriers")) { f <- pick(EF, "cross-sectional band", r$outcome, r$group, r$age_bin, tm); e <- pick(EE, "cross-sectional band", r$outcome, r$group, r$age_bin, tm)
    add(paste("band:", tm), paste(r$outcome, r$group, r$age_bin, sep = " | "), if (tm == "CNV x BMI") r$diff else r$slope_carriers, paste(r$n_carriers, r$n_controls, sep = " / "), f$estimate, paste(f$n_carriers, f$n_controls, sep = " / "), e$estimate) }
  f <- pick(EF, "cross-sectional band", r$outcome, r$group, r$age_bin, "CNV x BMI"); e <- pick(EE, "cross-sectional band", r$outcome, r$group, r$age_bin, "CNV x BMI")
  add("band: CNV x BMI, SE", paste(r$outcome, r$group, r$age_bin, sep = " | "), r$se_diff, "", f$se, "", e$se) }
fc <- fread(file.path(G, "band_slopes_family_clustered.csv"))[group %in% LOCI]
for (i in seq_len(nrow(fc))) { r <- fc[i]; f <- pick(EF, "cross-sectional band", r$outcome, r$group, r$age_bin, "CNV x BMI"); e <- pick(EE, "cross-sectional band", r$outcome, r$group, r$age_bin, "CNV x BMI")
  add("band: CNV x BMI, family-clustered SE", paste(r$outcome, r$group, r$age_bin, sep = " | "), r$se_diff_family_clustered, r$n_families, f$se_family_clustered, f$n_families, e$se_family_clustered) }
# 2. pooled band slopes (82) and all-ages terms
pb <- fread(file.path(G, "pooled_band_slopes.csv"))
for (i in seq_len(nrow(pb))) { r <- pb[i]; all_ages <- r$age_bin == "all ages"
  mod <- if (!all_ages) { if (grepl("primary", r$model)) "(primary)" else "no locus term" } else if (grepl("no locus", r$model)) "no locus term" else if (grepl("primary", r$model)) "(primary)" else "s(age, by = locus)"
  an <- if (all_ages) "cross-sectional all ages" else "cross-sectional band"; f <- pick(EF, an, r$outcome, POOL, r$age_bin, "CNV x BMI", mod = mod); e <- pick(EE, an, r$outcome, POOL, r$age_bin, "CNV x BMI", mod = mod)
  add(paste("pooled four CNVs: CNV x BMI,", if (all_ages) "all ages" else "band"), paste(r$outcome, r$age_bin, r$model, sep = " | "), r$diff, paste(r$n_carriers, r$n_controls, sep = " / "), f$estimate, paste(f$n_carriers, f$n_controls, sep = " / "), e$estimate) }
aa <- fread(file.path(G, "cnv_x_bmi_all_ages.csv"))[grepl("est_main of 17", quantity)]
for (i in seq_len(nrow(aa))) { r <- aa[i]; f <- pick(EF, "cross-sectional all ages", r$outcome, r$group, "all ages", "CNV x BMI"); e <- pick(EE, "cross-sectional all ages", r$outcome, r$group, "all ages", "CNV x BMI")
  add("all ages: CNV x BMI per genotype (bam)", paste(r$outcome, r$group, sep = " | "), r$est, paste(r$n_carriers, r$n_controls, sep = " / "), f$estimate, paste(f$n_carriers, f$n_controls, sep = " / "), e$estimate) }
# 3. sex: three-way term
sx <- fread(file.path(G, "sex_cnv_bmi_growth.csv"))
for (i in seq_len(nrow(sx))) { r <- sx[i]; f <- pick(EF, "sex", r$outcome, r$group, r$scope, "CNV x BMI x sex (male minus female)", "male - female"); e <- pick(EE, "sex", r$outcome, r$group, r$scope, "CNV x BMI x sex (male minus female)", "male - female")
  add("sex: CNV x BMI x sex", paste(r$outcome, r$group, r$scope, sep = " | "), r$cnv_x_bmi_x_sex, r$n_carriers_male + r$n_carriers_female, f$estimate, f$n_carriers, e$estimate) }
# 4. genotype-by-band F test
TF <- fread(file.path(OUT, "export_fixed_ukbb", "tests.csv")); TE <- fread(file.path(OUT, "export_estimate", "tests.csv")); gt <- fread(file.path(G, "main_with_ukbb_bmi_genotype_by_band.csv"))
for (i in seq_len(nrow(gt))) { r <- gt[i]; g1 <- function(T) { v <- T[analysis == "cross-sectional band" & outcome == r$outcome & stratum == r$age_bin & test == "BMI effect differs by genotype (F)", statistic]; if (length(v) == 1) v else NA_real_ }
  add("band: BMI effect differs by genotype (F statistic)", paste(r$outcome, r$age_bin, sep = " | "), r$F, r$groups, g1(TF), "", g1(TE)) }
# 5. longitudinal
for (o in c("Height", "Head circumference")) { L <- fread(file.path(V2, "longitudinal", if (o == "Height") "results_height.csv" else "results_head_circumference.csv"))
  for (i in seq_len(nrow(L))) { r <- L[i]; f <- pick(EF, "longitudinal", o, r$genotype, r$scope, "BMI Z at the visit")
    add(paste0("longitudinal: BMI coefficient (", if (o == "Height") "same table" else "analysis table only", ")"), paste(o, r$genotype, r$scope, sep = " | "), r$beta, paste(r$children, r$pairs, sep = " / "), f$estimate, paste(f$n_carriers, f$n_pairs, sep = " / "))
    if (o == "Head circumference") { h <- pick(EH, "longitudinal", o, r$genotype, r$scope, "BMI Z at the visit")
      add("longitudinal: BMI coefficient (head-circumference sources of 80)", paste(o, r$genotype, r$scope, sep = " | "), r$beta, paste(r$children, r$pairs, sep = " / "), h$estimate, paste(h$n_carriers, h$n_pairs, sep = " / ")) } } }
C <- rbindlist(cmp); C[, `:=`(diff = federated - paper, diff_estimate_mode = federated_estimate_mode - paper)]
C[, outcome := sub(" [|].*", "", item)][, band := fifelse(grepl("2-9", item), "2-9", fifelse(grepl("10-18", item), "10-18", fifelse(grepl("19-39", item), "19-39", fifelse(grepl("40-73", item), "40-73", "all ages / other"))))]
fwrite(C, file.path(OUT, "comparison_with_paper_run.csv"))
S <- C[, .(rows = .N, rows_in_both = sum(!is.na(diff)), max_abs_diff_fixed_ukbb = if (any(!is.na(diff))) max(abs(diff), na.rm = TRUE) else NA_real_, exact_fixed_ukbb = sum(abs(diff) < 1e-8, na.rm = TRUE),
           max_abs_diff_estimate_mode = if (any(!is.na(diff_estimate_mode))) max(abs(diff_estimate_mode), na.rm = TRUE) else NA_real_, exact_estimate_mode = sum(abs(diff_estimate_mode) < 1e-8, na.rm = TRUE)), by = .(block, outcome, band)]
fwrite(S, file.path(OUT, "comparison_summary.csv")); options(width = 250); print(S, nrows = 300)
