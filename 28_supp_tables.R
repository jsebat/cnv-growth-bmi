# Supplementary tables (JS / Claude 2026-10-24): all statistics behind the main figure v3 and the CNB figures, with carriers
# compared to (a) unaffected non-carrier relatives (SSC siblings/parents, + SVIP non-ASD for height/HC) and (b) SSC probands
# with ASD and no CNV. Table S1: height and head circumference, CNV x BMI (age-adjusted constant shift) and CNV x age x BMI
# vs each reference group, per locus and pooled, plus the within-group BMI effect and age x BMI. Table S2: SRS, SCQ, CNB
# accuracy and speed: age-adjusted BMI slope per group, and the difference of that slope from each reference group
# (CNV x BMI) with its age interaction (CNV x age x BMI); CNB has no reference group with BMI. Models as in 23/24
# (GAM/bam, mgcv). Data: QC'd build of 21 (00_qc_growth.R). Outputs: supp_tables_bmi_effects.xlsx (sheets S1, S2, notes),
# supp_table_S1_growth.csv, supp_table_S2_cognition.csv, supp_tables_bmi_effects.md
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))   # A (growth + SRS/SCQ with unaffected controls), gmc, car_g, mk_growth, LOCI, POOLS
source(file.path(HERE, "24_cognition_bmi_slopes.R"))  # X (SRS/SCQ/CNB with unaffected relatives and SSC ASD probands)
suppressPackageStartupMessages({ library(mgcv); library(openxlsx) })
t0 <- Sys.time()
REFS <- c("Unaffected controls", "SSC ASD probands")

# ---------------------------------------------------------------- growth rows for SSC ASD probands (non-carriers), same build as 23
pro_g <- gmc[cohort == "SSC" & CNV == "No CNV" & grepl("[.]p[0-9]$", sampleID), .(id, cohort, sampleID, GUID, CNV, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd = 1L)][, locus := "SSC ASD probands"]
GP <- rbind(mk_growth(pro_g, "HEIGHT_raw", "Height"), mk_growth(pro_g, "HEADCIRCUM_raw", "Head circumference"))
GP[, `:=`(unaffected = FALSE, role = "asd proband")]
G2 <- rbind(A[pheno %in% c("Height", "Head circumference"), .(id, locus, age, BMI, y, pheno, asd)], GP[, .(id, locus, age, BMI, y, pheno, asd)])
# cognition rows: X (grp = locus labels / "Unaffected relatives" / "SSC ASD probands")
C2 <- X[, .(id, locus = fifelse(grp == "Unaffected relatives", "Unaffected controls", grp), age, BMI, y, pheno, asd)]
ALL <- rbind(G2, C2)
cat("== rows by trait x group\n"); print(dcast(ALL[, .N, by = .(pheno, locus)], locus ~ pheno, value.var = "N", fill = 0))

# ---------------------------------------------------------------- model functions (as in 23 / 24)
within <- function(x) tryCatch({                       # BMI main effect (age-adjusted) and age x BMI within a group; locus covariate when pooled
  cov <- if (uniqueN(x$locus) > 1) "+ locus" else ""; x <- copy(x)[, locus := factor(locus)]
  m1 <- gam(as.formula(paste("y ~ BMI + s(age, k = 8)", cov)), data = x, method = "REML"); co <- summary(m1)$p.table["BMI", ]
  m2 <- gam(as.formula(paste("y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6))", cov)), data = x, method = "REML")
  data.table(n = nrow(x), bmi_slope = unname(co[1]), bmi_se = unname(co[2]), bmi_p = unname(co[4]), age_x_bmi_p = summary(m2)$s.table["ti(age,BMI)", "p-value"]) },
  error = function(e) data.table(n = nrow(x), bmi_slope = NA_real_, bmi_se = NA_real_, bmi_p = NA_real_, age_x_bmi_p = NA_real_))
contrast2 <- function(ref, xc) tryCatch({             # carriers vs a reference group: CNV main, CNV x BMI shift, CNV x age x BMI
  cov <- if (uniqueN(xc$locus) > 1) "locus + " else ""
  z <- rbind(ref[, .(id, locus = "ref", age, BMI, y)], xc[, .(id, locus, age, BMI, y)])[, is_cnv := factor(ifelse(locus == "ref", "ref", "CNV"), levels = c("ref", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  z[, locus := factor(ifelse(locus == "ref", levels(factor(xc$locus))[1], locus))]   # locus covariate: reference rows take the first carrier locus (absorbed by is_cnv)
  m1 <- bam(as.formula(paste("y ~", cov, "is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI")), data = z, method = "fREML"); co <- summary(m1)$p.table
  m2 <- bam(as.formula(paste("y ~", cov, "is_cnv + s(age, k = 10) + s(BMI, k = 6) + s(age, k = 10, by = is_cnv) + ti(age, BMI, k = c(6, 6)) + ti(age, BMI, k = c(6, 6), by = is_cnv)")), data = z, method = "fREML"); st <- summary(m2)$s.table
  data.table(n_ref = nrow(ref), cnv_main = co["is_cnvCNV", 1], cnv_main_p = co["is_cnvCNV", 4], cnv_x_bmi = co["cnvBMI", 1], cnv_x_bmi_se = co["cnvBMI", 2], cnv_x_bmi_p = co["cnvBMI", 4],
             cnv_x_age_x_bmi_p = st[grep("^ti\\(age,BMI\\):is_cnv", rownames(st))[1], "p-value"]) },
  error = function(e) data.table(n_ref = nrow(ref), cnv_main = NA_real_, cnv_main_p = NA_real_, cnv_x_bmi = NA_real_, cnv_x_bmi_se = NA_real_, cnv_x_bmi_p = NA_real_, cnv_x_age_x_bmi_p = NA_real_))

GROUPS <- c(setNames(as.list(LOCI), LOCI), POOLS)
PH <- c("Height", "Head circumference", "SRS", "SCQ", "CNB accuracy", "CNB speed")
rows <- list()
for (ph in PH) {
  d <- ALL[pheno == ph]
  for (r in REFS) { xr <- d[locus == r]; if (nrow(xr) >= 30) rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = r, members = r, reference = NA_character_), within(xr)) }
  for (g in names(GROUPS)) {
    xc <- d[locus %in% GROUPS[[g]]]; if (nrow(xc) < 15) next
    w <- within(xc)
    for (r in REFS) { xr <- d[locus == r]
      rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = g, members = paste(GROUPS[[g]], collapse = " + "), reference = r), w, if (nrow(xr) >= 30) contrast2(xr, xc) else data.table(n_ref = nrow(xr))) }
    if (ph %like% "CNB") rows[[length(rows) + 1]] <- cbind(data.table(pheno = ph, group = g, members = paste(GROUPS[[g]], collapse = " + "), reference = "none (no reference group with BMI)"), w)
  }
  cat(ph, "done", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "min\n")
}
T <- rbindlist(rows, fill = TRUE)
T[pheno %like% "CNB" & reference %in% REFS, (c("n_ref", "cnv_main", "cnv_main_p", "cnv_x_bmi", "cnv_x_bmi_se", "cnv_x_bmi_p", "cnv_x_age_x_bmi_p")) := NA]
T <- T[!(pheno %like% "CNB" & reference %in% REFS)]
T[, pheno := factor(pheno, PH)][, group := factor(group, c(REFS, names(GROUPS)))]; setorder(T, pheno, group, reference)
S1 <- T[pheno %in% c("Height", "Head circumference")]; S2 <- T[!pheno %in% c("Height", "Head circumference")]
fwrite(S1, file.path(HERE, "supp_table_S1_growth.csv")); fwrite(S2, file.path(HERE, "supp_table_S2_cognition.csv"))

# ---------------------------------------------------------------- formatted versions
fmt <- function(x) { x <- copy(x)
  x[, `:=`(`BMI slope (SE)` = ifelse(is.na(bmi_slope), "", sprintf("%+.3f (%.3f)", bmi_slope, bmi_se)), `BMI slope p` = pfx(bmi_p), `age x BMI p` = pfx(age_x_bmi_p),
           `CNV main effect (SD)` = ifelse(is.na(cnv_main), "", sprintf("%+.2f", cnv_main)), `CNV main p` = ifelse(is.na(cnv_main_p), "", pfx(cnv_main_p)),
           `CNV x BMI (SE)` = ifelse(is.na(cnv_x_bmi), "", sprintf("%+.3f (%.3f)", cnv_x_bmi, cnv_x_bmi_se)), `CNV x BMI p` = ifelse(is.na(cnv_x_bmi_p), "", pfx(cnv_x_bmi_p)), `CNV x age x BMI p` = ifelse(is.na(cnv_x_age_x_bmi_p), "", pfx(cnv_x_age_x_bmi_p)))]
  x[, .(Trait = as.character(pheno), Group = as.character(group), Members = members, `Reference group` = ifelse(is.na(reference), "", reference), n = n, `n reference` = ifelse(is.na(n_ref), "", as.character(n_ref)),
        `BMI slope (SE)`, `BMI slope p`, `age x BMI p`, `CNV main effect (SD)`, `CNV main p`, `CNV x BMI (SE)`, `CNV x BMI p`, `CNV x age x BMI p`)] }
F1 <- fmt(S1); F2 <- fmt(S2)
notes <- data.table(Note = c(
  "Table S1: height and head circumference (WHO / UK90 Z), clinical cohorts, ages 2-65, latest visit with the outcome and BMI (plausibility QC applied, 00_qc_growth.R). Groups: 16p11.2 proximal deletion / duplication, 22q11.2 A-D deletion / duplication; pooled overgrowth = 16p11.2 del + 22q11.2 dup, undergrowth = 16p11.2 dup + 22q11.2 del, all CNVs (locus as a covariate). Reference groups: unaffected controls = SSC designated siblings and parents (ages at the study visit) plus SVIP non-carrier family members without ASD; SSC ASD probands = SSC probands with no CNV.",
  "Table S2: SRS total (Z vs SSC unaffected relatives), SCQ lifetime (Z vs SSC siblings), Penn CNB accuracy and speed (Z vs G2MH non-carriers; speed = -RT z). Reference groups as in S1 (SCQ: siblings only; CNB: the G2MH non-carriers have no BMI, so no reference contrast).",
  "BMI: WHO child BMI Z; adults, centile vs pooled non-carrier adults mapped to Z. All tests use continuous BMI.",
  "BMI slope: within-group age-adjusted BMI effect, SD of trait per SD of BMI (GAM y ~ BMI + s(age) [+ locus]). age x BMI p: tensor-product age x BMI interaction within the group (GAM y ~ s(age) + s(BMI) + ti(age, BMI) [+ locus]).",
  "CNV main effect: age-adjusted mean difference carriers - reference (SD) from the CNV x BMI model. CNV x BMI: difference between carriers and the reference group in the BMI slope, adjusted for age (bam y ~ [locus +] CNV + s(age, by CNV) + s(age, by BMI) + CNV:BMI). CNV x age x BMI: whether that difference varies with age (bam y ~ [locus +] CNV + s(age) + s(BMI) + s(age, by CNV) + ti(age, BMI) + ti(age, BMI, by CNV); p of the carrier-specific tensor product). The three-way terms are sensitive to extreme observations (see unaffected_controls_report.md section 9).",
  "Reference-group rows give the within-group BMI slope and age x BMI for that reference group itself.",
  sprintf("Generated %s by 28_supp_tables.R.", format(Sys.Date()))))
wb <- createWorkbook(); addWorksheet(wb, "S1 growth"); addWorksheet(wb, "S2 cognition"); addWorksheet(wb, "S1 raw"); addWorksheet(wb, "S2 raw"); addWorksheet(wb, "notes")
writeData(wb, "S1 growth", F1); writeData(wb, "S2 cognition", F2); writeData(wb, "S1 raw", S1); writeData(wb, "S2 raw", S2); writeData(wb, "notes", notes)
for (s in c("S1 growth", "S2 cognition")) { setColWidths(wb, s, cols = 1:14, widths = c(18, 18, 34, 28, 7, 11, 16, 11, 11, 18, 10, 16, 12, 16)); addStyle(wb, s, createStyle(textDecoration = "bold"), rows = 1, cols = 1:14) }
setColWidths(wb, "notes", cols = 1, widths = 160)
saveWorkbook(wb, file.path(HERE, "supp_tables_bmi_effects.xlsx"), overwrite = TRUE)
mdt <- function(x, title) { h <- names(x); c(sprintf("## %s", title), "", paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", paste(rep("---", length(h)), collapse = "|"), "|"), apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "") }
writeLines(c("# Supplementary tables: BMI effects on growth and social-communication / cognitive measures in CNV carriers", "", paste0("- ", notes$Note), "", mdt(F1, "Table S1. Height and head circumference"), mdt(F2, "Table S2. SRS, SCQ and Penn CNB")), file.path(HERE, "supp_tables_bmi_effects.md"))
cat("\n== S1\n"); print(F1[, .(Trait, Group, `Reference group`, n, `n reference`, `BMI slope (SE)`, `BMI slope p`, `CNV x BMI (SE)`, `CNV x BMI p`, `CNV x age x BMI p`)], nrows = 100)
cat("\n== S2\n"); print(F2[, .(Trait, Group, `Reference group`, n, `n reference`, `BMI slope (SE)`, `BMI slope p`, `CNV x BMI (SE)`, `CNV x BMI p`, `CNV x age x BMI p`)], nrows = 100)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
