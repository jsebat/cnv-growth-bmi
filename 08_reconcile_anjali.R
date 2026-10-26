# Why do our CNV x age x BMI p-values for 16p11.2 del differ from Anjali's (head circumference 0.0028, height 0.00092)?
# Start from her exact specification (compute_interaction_results.R: genomatrix last_timepoint rows, Cardiff relabelled
# G2MH and NOT de-duplicated, age < 45, default gam() smoothing selection) and change one thing at a time:
#   A  her spec exactly                     -> should reproduce 0.0028 / 0.00092
#   B  A + Cardiff/G2MH duplicates dropped
#   C  A + REML smoothing selection (ours)
#   D  our data selection (latest visit WITH the outcome and BMI, age 2-40, duplicates dropped), her default gam()
#   E  ours (= D + REML)                   -> 0.32 / 0.11
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")   # participant-level tables stay in Jon/Claude (controlled data)
g <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
g[, CNV := gsub("proximal", "BP4-BP5", trimws(CNV))]

fit3 <- function(d, outcome, reml) {
  d <- copy(d); d[, is_cnv := factor(ifelse(CNV == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  f <- as.formula(sprintf("%1$s ~ is_cnv + s(age, k = 15) + s(BMI, k = 10) + s(age, k = 15, by = is_cnv) +
                           ti(age, BMI, k = c(10, 10)) + ti(age, BMI, k = c(10, 10), by = is_cnv)", outcome))
  m <- if (reml) gam(f, data = d, method = "REML") else gam(f, data = d)
  summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]
}
sel_hers <- function(dedup) {
  x <- g[last_timepoint == TRUE]
  if (dedup) x <- x[!(cohort == "G2MH" & Cardiff_mask == TRUE)]
  x[, cohort := ifelse(cohort == "Cardiff", "G2MH", cohort)]; x
}
sel_ours <- function(outcome) {
  x <- g[!(cohort == "G2MH" & Cardiff_mask == TRUE) & !is.na(get(outcome)) & !is.na(BMI) & age >= 2 & age <= 40]
  x[order(MSassigned_ID, cohort, -age)][, .SD[1], by = .(cohort, MSassigned_ID)]
}
res <- list()
for (oc in c("HEADCIRCUM", "HEIGHT")) {
  hers <- function(x) x[CNV %in% c("16p11.2 BP4-BP5 DEL", "No CNV") & !is.na(get(oc)) & !is.na(BMI) & !is.na(age) & age < 45]
  for (v in list(list("A her spec", hers(sel_hers(FALSE)), FALSE), list("B + dedup", hers(sel_hers(TRUE)), FALSE),
                 list("C + REML", hers(sel_hers(FALSE)), TRUE),
                 list("D our selection, default gam", sel_ours(oc)[CNV %in% c("16p11.2 BP4-BP5 DEL", "No CNV")], FALSE),
                 list("E ours (REML)", sel_ours(oc)[CNV %in% c("16p11.2 BP4-BP5 DEL", "No CNV")], TRUE))) {
    t0 <- Sys.time(); p <- fit3(v[[2]], oc, v[[3]])
    r <- data.table(outcome = oc, variant = v[[1]], n_cnv = v[[2]][CNV != "No CNV", .N], n_ctrl = v[[2]][CNV == "No CNV", .N], p = p,
                    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
    print(r); res[[length(res) + 1]] <- r
    fwrite(rbindlist(res), file.path(HERE, "reconcile_anjali_16p_del.csv"))
  }
}
