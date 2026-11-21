# 00_control_set.R  (JS / Claude 2026-11-20): the non-carrier CONTROL SET of the growth tables, corrected.
#
# Problem found 2026-11-20 (coordinator, from the sex table of 82_growth_blocks_v2.R: control children 87% male): in Molly's
#   clinical growth genomatrix the SSC "No CNV" rows are fathers (.fa), mothers (.mo), PROBANDS (.p1, .p2) and siblings (.s1-.s5),
#   and only the parents and the probands carry height / weight / head circumference; the siblings have none. 01_build.R kept every
#   "No CNV" row, so the child "non-carrier controls" of 17_main_with_ukbb.R were the SSC autistic probands (2,074 of the 6,283
#   clinical height controls of the 2026-11-19 run), and no SSC sibling contributed. The manuscript describes the controls as "SSC
#   designated siblings and parents and Simons VIP non-carrier family members without autism".
# What apply_control_set() does (CONTROL_SET = "unaffected", the default):
#   * SSC probands without a CNV keep their rows but get CNV = "No CNV, autism (SSC proband)", so they are never in "No CNV".
#   * SSC designated siblings without a CNV get height, weight (-> BMI) and head circumference from the SSC sibling table
#     (Designated Unaffected Sibling Data/ssc_commonly_used.csv) and their age is replaced by the age at the visit
#     (proband age at evaluation + birth-month difference; the genomatrix age is the 2012-09 snapshot, ~3 y too old; same fix as
#     21_pheno_bmi_gam.R / 23_unaffected_controls.R). Siblings without a visit age keep no measurements.
#   * Simons VIP non-carriers with a clinical autism diagnosis (v10 diagnosis_summary, clinical_asd_dx = true) get
#     CNV = "No CNV, autism (SVIP)"; those without a recorded answer get "No CNV, autism status unknown (SVIP)".
#   * G2MH non-carriers (2 people with measurements; not part of the stated control set, autism status not established here) get
#     CNV = "No CNV, not in control set (G2MH)".
#   * SSC parents and Simons VIP non-carriers without autism stay "No CNV". Carriers are untouched.
#   CONTROL_SET = "paper_2026-11-19" leaves every "No CNV" row as it was (the control set of the runs up to v10).
#   Both settings add the columns family (SSC / SVIP family id from the sample id; own id otherwise) and control_role.
# Inputs: genomatrix-format table g (cohort, sampleID, CNV, sex_numeric, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw); SSC v15
#   phenotype folder; SVIP v10 diagnosis_summary.csv. Output: g, modified. Called from 01_build.R before the plausibility QC.
# null model: none; data preparation only.

CONTROL_FLAGS <- c("No CNV, autism (SSC proband)", "No CNV, autism (SVIP)", "No CNV, autism status unknown (SVIP)", "No CNV, not in control set (G2MH)")

apply_control_set <- function(g, DB, set = Sys.getenv("CONTROL_SET", "unaffected")) {
  stopifnot(set %in% c("unaffected", "paper_2026-11-19"))
  g <- data.table::copy(g)
  suf <- ifelse(g$cohort == "SSC", sub("^.*[.]", "", g$sampleID), NA_character_)
  g[, family := ifelse(cohort %in% c("SSC", "SVIP") & grepl("[.]", sampleID), paste0(cohort, ":", sub("[.].*$", "", sampleID)), NA_character_)]
  g[, control_role := data.table::fcase(cohort == "SSC" & suf == "fa", "SSC father", cohort == "SSC" & suf == "mo", "SSC mother",
                                        cohort == "SSC" & grepl("^p[0-9]+$", suf), "SSC proband", cohort == "SSC" & suf == "s1", "SSC designated sibling",
                                        cohort == "SSC" & grepl("^s[0-9]+$", suf), "SSC other sibling", cohort == "SVIP", "SVIP family member",
                                        default = NA_character_)]
  g[CNV != "No CNV", control_role := NA_character_]
  if (set == "paper_2026-11-19") { cat("CONTROL_SET = paper_2026-11-19: 'No CNV' rows left as supplied (SSC probands are among the controls)\n"); return(g) }
  nc <- g$CNV == "No CNV"
  # SSC probands
  pro <- nc & g$control_role %in% "SSC proband"; g[pro, CNV := CONTROL_FLAGS[1]]
  # Simons VIP non-carriers by clinical autism diagnosis
  dx <- data.table::fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"), colClasses = list(character = "diagnosis_summary.clinical_asd_dx"))
  sv <- setNames(tolower(dx$diagnosis_summary.clinical_asd_dx), dx$individual)
  sva <- nc & g$cohort == "SVIP"; a <- sv[g$sampleID]
  g[sva & a %in% "true", CNV := CONTROL_FLAGS[2]]; g[sva & !(a %in% c("true", "false")), CNV := CONTROL_FLAGS[3]]
  # G2MH non-carriers
  g2 <- nc & g$cohort == "G2MH"; g[g2, CNV := CONTROL_FLAGS[4]]
  # SSC designated siblings: measurements from the SSC sibling table, age at the visit
  SSCD <- file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9")
  rd <- function(sub, f) { txt <- paste(readLines(file.path(SSCD, sub, f), warn = FALSE), collapse = "\n"); txt <- gsub("\r\n|\r", "\n", txt); data.table::fread(text = txt) }   # CR-only line endings in some SSC v15 files
  cu <- rd("Designated Unaffected Sibling Data", "ssc_commonly_used.csv")[, .(sid = individual, h = suppressWarnings(as.numeric(height)), w = suppressWarnings(as.numeric(weight)), hc = suppressWarnings(as.numeric(head_circumference)))]
  agef <- data.table::fread(file.path(SSCD, "Age_individual_v14_parent_age.txt")); data.table::setnames(agef, c("id()", "individual_properties.role"), c("sid", "role0"))
  agef[, bm := as.integer(substr(birth, 1, 4)) * 12L + as.integer(substr(birth, 6, 7))]
  ev <- data.table::fread(file.path(SSCD, "ssc_case_HC.csv"))[, .(proband = Proband_ID, age_eval = `Age at Evaluation (Non-PHI) :: Proband SSC Height, Weight, Head Circumference`)]
  pr <- merge(agef[role0 == "proband", .(family, pbm = bm, proband = sid)], ev, by = "proband")
  sa <- merge(agef[role0 == "designated-sibling", .(sid, family, obm = bm)], pr, by = "family")[, .(sid, vage = (age_eval + pbm - obm) / 12)][!is.na(vage)]
  sib <- merge(cu, sa, by = "sid")[!duplicated(sid)]
  i <- which(g$CNV == "No CNV" & g$cohort == "SSC" & grepl("^s[0-9]+$", suf) & g$sampleID %in% sib$sid); j <- match(g$sampleID[i], sib$sid)
  had <- sum(!is.na(g$HEIGHT_raw[i]) | !is.na(g$HEADCIRCUM_raw[i]))
  data.table::set(g, i, "age", sib$vage[j]); data.table::set(g, i, "HEIGHT_raw", sib$h[j])
  data.table::set(g, i, "BMI_raw", ifelse(!is.na(sib$h[j]) & !is.na(sib$w[j]) & sib$h[j] > 0, sib$w[j] / (sib$h[j] / 100)^2, NA_real_))
  data.table::set(g, i, "HEADCIRCUM_raw", sib$hc[j])
  cat(sprintf("CONTROL_SET = unaffected: flagged %d SSC probands, %d SVIP non-carriers with autism, %d SVIP with unknown autism status, %d G2MH non-carriers (rows; kept, not controls)\n",
              sum(pro), sum(sva & a %in% "true"), sum(sva & !(a %in% c("true", "false"))), sum(g2)))
  cat(sprintf("  SSC siblings given measurements from the sibling table: %d rows (%d had a measurement before); with height %d, with BMI %d, with head circumference %d; age range %.1f-%.1f y\n",
              length(i), had, sum(!is.na(g$HEIGHT_raw[i])), sum(!is.na(g$BMI_raw[i])), sum(!is.na(g$HEADCIRCUM_raw[i])), min(g$age[i]), max(g$age[i])))
  g
}
