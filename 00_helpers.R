# Shared data helpers for the longitudinal scripts (03, 05): cleaning, per-person visit table, intervals.
suppressPackageStartupMessages({ library(data.table) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")

# Growth-record cleaning (found 2026-10-12: SVIP carries the last height forward to later visits while weight is
# updated, which by itself makes BMI Z rise as height Z falls; some visits are entered twice a month apart; Cardiff
# and SVIP have impossible shrinkage).  Walking each person's visits in age order, a visit is dropped when, relative
# to the last kept visit, it is
#   (a) < 0.25 y later (duplicate entry),
#   (b) height: identical raw height (|diff| < 0.1 cm) >= 0.25 y later before age 16 (carried forward),
#   (c) height: > 2 cm shorter, or growing > 20 cm/y after age 3,
#   (d) outcome Z changed by > 2.5 SD.
CLEAN_LOG <- list()
clean_visits <- function(d, outcome) {
  raw <- if (outcome == "HEIGHT") "HEIGHT_raw" else "HEADCIRCUM_raw"
  keep <- logical(nrow(d)); why <- character(nrow(d))
  for (ix in split(seq_len(nrow(d)), d$id)) {
    last <- ix[1]; keep[last] <- TRUE
    for (i in ix[-1]) {
      dt <- d$age[i] - d$age[last]; dr <- d[[raw]][i] - d[[raw]][last]; dz <- d[[outcome]][i] - d[[outcome]][last]
      w <- if (dt < 0.25) "a_duplicate_visit"
           else if (outcome == "HEIGHT" && !is.na(dr) && abs(dr) < 0.1 && d$age[i] < 16) "b_carried_forward"
           else if (outcome == "HEIGHT" && !is.na(dr) && (dr < -2 || (d$age[i] > 3 && dr / dt > 20))) "c_impossible_change"
           else if (abs(dz) > 2.5) "d_jump_gt_2.5SD" else ""
      if (w == "") { keep[i] <- TRUE; last <- i } else why[i] <- w
    }
  }
  CLEAN_LOG[[length(CLEAN_LOG) + 1]] <<- data.table(outcome = outcome, cohort = d$cohort, why = why)[why != ""]
  d[keep]
}

long_data <- function(outcome, subset = "all") {
  d <- tab[locus %in% LOCI & age < 18 & !is.na(tab[[outcome]]) & !is.na(BMI)]
  if (subset == "exact_age") d <- d[cohort %in% c("SVIP", "Searchlight")]
  if (subset == "Cardiff") d <- d[cohort == "Cardiff"]
  if (subset == "no_meds") d <- d[!(stimulant %in% TRUE) & !(antipsychotic %in% TRUE)]
  d <- d[order(id, age)]
  if (subset != "uncleaned") d <- clean_visits(d, outcome)
  d <- unique(d, by = c("id", "age"))
  d <- d[, if (.N >= 2) .SD, by = id]
  if (nrow(d) == 0) return(d)
  d[, Z := d[[outcome]]][, locus := droplevels(locus)]
  d[, BMIb := mean(BMI), by = id][, BMIw := BMI - BMIb]
  d
}

intervals <- function(d) {
  iv <- d[, .(age0 = head(age, -1), dt = diff(age), Z0 = head(Z, -1), dZ = diff(Z), B0 = head(BMI, -1), dB = diff(BMI),
              sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id]
  iv[dt >= 0.5]
}

