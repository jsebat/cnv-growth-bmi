# 90_supp_tables_v2_2026-11-20.R  (2026-11-20, BMI modifier paper v2: supplementary tables workbook)
#
# What it does: builds the supplementary tables workbook of the v2 draft. It copies "Sacks Supplementary Tables v3_chop.xlsx"
#   (never modified) into a NEW workbook "Sacks Supplementary Tables v4_v2draft.xlsx" and
#     replaces  T13 (same table; title and note rows added, note starts "v1 analysis; superseded by Supplementary Table 30"),
#               T17 (autism diagnosis by CNV, v2), T21 (BMI slope by age band, corrected controls, with the rerun 0-2 y block),
#               T22 (carriers versus controls with UK Biobank adults, corrected controls);
#     adds      T30 longitudinal, T31 control set, T32 growth sensitivity, T33 sex tests (growth), T34 autism diagnosis by BMI,
#               T35 polygenic score, T36 sex tests (autism), T37 puberty Mendelian randomization;
#     replaces  T14, T15, T16 with the 30b run of 2026-11-14 (supp_table_S4/S5/S6*_boxcox.csv; Jon 2026-11-20: v3 held the older
#               30_sex_pgs_sensitivity.R output of 2026-10-26; both runs use the unaffected control set of 23_unaffected_controls.R);
#     renames   the v3 sheet "T30 Extended CNV models" to "T38 Extended CNV models" (title row too; content untouched), so the
#               number 30 belongs to the new longitudinal table only;
#     leaves    every other sheet as in v3 (T18 as in v3).
#   Every sheet written here: row 1 title, row 2 plain-language note, table from row 4 (bold header); stacked blocks carry a
#   sub-title line each. Estimates rounded to 3 decimals, P to 2 significant figures. No participant-level rows. Counts of
#   people below 5 are written "<5" (and the proportions that would reveal them are blanked).
# Inputs (aggregate result tables only): run_v2_2026-11-20/{growth_unaffctl, longitudinal, autism, puberty_mr, infancy_0to2}/*.csv
#   (growth_unaffctl = corrected control set; run_v2_2026-11-20/growth/ = old controls, never read here), and the v3 workbook.
# Outputs: manuscript_materials_2026-10-24/Sacks Supplementary Tables v4_v2draft.xlsx; the same tables as CSV in
#   run_v2_2026-11-20/supplementary_tables/.
# null model: none; this script only formats result tables (the tests are those of the scripts that wrote the inputs).
suppressMessages({ library(data.table); library(openxlsx) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
RUN <- file.path(HERE, "run_v2_2026-11-20"); G <- file.path(RUN, "growth_unaffctl"); L <- file.path(RUN, "longitudinal")
A <- file.path(RUN, "autism"); PM <- file.path(RUN, "puberty_mr"); INF <- file.path(RUN, "infancy_0to2")
ST <- file.path(RUN, "supplementary_tables"); dir.create(ST, showWarnings = FALSE)
MM <- file.path(HERE, "manuscript_materials_2026-10-24")
wbf3 <- file.path(MM, "Sacks Supplementary Tables v3_chop.xlsx"); wbf4 <- file.path(MM, "Sacks Supplementary Tables v4_v2draft.xlsx")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
CTRL <- "Controls: SSC designated siblings and parents and Simons VIP non-carrier relatives without autism."
DATE <- " Built 2026-11-20."
CLIN <- "clinical cohorts (G2MH, Simons VIP, Searchlight, Cardiff, CHU Sainte-Justine, SSC carriers, and the CHOP 22q11.2 A-D records; G2MH carriers recruited at CHOP removed), one visit per person chosen to fill sparse age bands"
bold <- createStyle(textDecoration = "bold"); ital <- createStyle(textDecoration = "italic"); wrap <- createStyle(wrapText = FALSE)
rd <- function(...) fread(file.path(...), encoding = "UTF-8")

# ---------------------------------------------------------------- formatting helpers
NOCOUNT <- c("n_snps", "n_perm", "n_as_extreme", "n_abs_diff_gt_0.1", "n_outcome", "n_families")
is_p <- function(nm) grepl("(^|_)p($|_)", nm) | grepl("^P($|[ ,])", nm)
is_count <- function(nm) (grepl("^(n|N)($|_)|^children|_children$|^pairs|^with_2plus|^non_carriers$|^with_more|\\(n\\)$", nm)) & !(nm %in% NOCOUNT)
fmt_p <- function(v) { z <- !is.na(v) & v == 0; r <- signif(v, 2); if (any(z)) { r <- as.character(r); r[z] <- "<1e-300" }; r }
mask5 <- function(v) { if (!is.numeric(v) || !any(v < 5, na.rm = TRUE)) return(v); r <- as.character(v); r[!is.na(v) & v < 5] <- "<5"; r }
# round estimates to 3 decimals, P to 2 significant figures; mask person counts below 5 (mask = FALSE for tables that count
# visits or matched records rather than people)
fmt <- function(x, mask = TRUE) {
  x <- copy(as.data.table(x))
  for (nm in names(x)) { v <- x[[nm]]; if (!is.numeric(v)) next
    if (is_p(nm)) set(x, j = nm, value = fmt_p(v))
    else if (is_count(nm)) { if (mask) set(x, j = nm, value = mask5(v)) }
    else if (any(abs(v - round(v)) > 1e-9, na.rm = TRUE)) set(x, j = nm, value = round(v, 3)) }
  x
}
# plain-word column names: exact dictionary first, then a token-by-token translation
DICT <- c(outcome = "Outcome", genotype = "Genotype", scope = "Age window", cohorts = "Cohorts", cohort = "Cohort", children = "Children (n)", pairs = "Visit pairs (n)",
          median_age0 = "Median age at first visit of pair (years)", median_interval = "Median interval (years)", q25_age0 = "Age at first visit, 25th centile", q75_age0 = "Age at first visit, 75th centile",
          beta = "Change in Z to next visit per SD of BMI Z", se = "SE", lo = "Lower 95% CI", hi = "Upper 95% CI", p = "P", p_bmi_x_window = "P, BMI x age window",
          p_heterogeneity = "P, heterogeneity across genotypes", sensitivity = "Analysis", window = "Age window", group = "Group", age_bin = "Age band (years)", band = "Band",
          n_carriers = "Carriers (n)", n_controls = "Controls (n)", pheno = "Trait", trait = "Trait", analysis = "Analysis", model = "Model", note = "Note", n = "n",
          n_noasd = "No autism diagnosis (n)", n_asd = "Autism diagnosis (n)", mean_noasd = "Mean Z, no autism diagnosis", se_noasd = "SE (no autism diagnosis)", mean_asd = "Mean Z, autism diagnosis",
          se_asd = "SE (autism diagnosis)", age_med_noasd = "Median age, no autism diagnosis", age_med_asd = "Median age, autism diagnosis", delta = "Difference, autism diagnosis minus none (SD)",
          het_Q = "Cochran Q", het_df = "Heterogeneity df", het_p = "P, heterogeneity", het_I2 = "I2", asd_x_cnv_F_p = "P, autism x CNV (F test)", het_F = "Heterogeneity F",
          est = "Estimate", estimate = "Estimate", b = "Estimate", lower95 = "Lower 95% CI", upper95 = "Upper 95% CI", ci_lo = "Lower 95% CI", ci_hi = "Upper 95% CI", OR = "Odds ratio per SD",
          OR_lo = "Odds ratio, lower 95% CI", OR_hi = "Odds ratio, upper 95% CI", logOR = "Log odds ratio per SD", prop = "Share with autism diagnosis", mean_z = "Mean Z of the predictor",
          predictor = "Predictor", term = "Term", sample = "Sample", stratum = "Stratum", covariates = "Covariates", estimator = "Estimator", first_stage_F = "First-stage F",
          sargan_p = "P, Sargan over-identification", ar_lo = "Anderson-Rubin 95% set, lower", ar_hi = "Anderson-Rubin 95% set, upper", ar_note = "Anderson-Rubin set",
          mde80 = "Smallest effect detectable with 80% power", power_at_observational_slope = "Power at the observational slope", t_obs = "t", perm_p = "P, permutation",
          n_perm = "Permutations (n)", n_as_extreme = "Permutations as extreme (n)", robust_t = "Robust t (HC1)", subset = "Subset", kind = "Score bin type", bin = "Bin", 
          slope = "BMI slope", test = "Test", sex = "Sex", n_ref = "Reference group (n)", status = "Status", block = "Block", n_families = "Families (n)", se_model = "SE (model)", p_model = "P (model)",
          se_family_clustered = "SE (family-clustered)", p_family_clustered = "P (family-clustered)", scale = "Scale", se_ratio = "SE ratio, clustered / model", crosses_0.05 = "Crosses P = 0.05",
          conclusion_changes = "Conclusion changes", carriers_sharing_a_family = "Carriers sharing a family (n)", step = "Step", unit = "Unit", F_stat = "F", R2_first_stage = "R2, first stage",
          method = "Method", n_snps = "SNPs (n)", mean_F = "Mean F of instruments", n_outcome = "Outcome sample (n)", Q = "Cochran Q", Q_p = "P, Cochran Q", egger_intercept = "Egger intercept",
          egger_intercept_se = "SE (Egger intercept)", egger_intercept_p = "P, Egger intercept", check = "Check", what_was_done = "What was done", result = "Result", conclusion_holds = "Conclusion holds",
          source_file = "Source file", sexes = "Sexes", measure = "Measure", p_interaction_vs_noncarriers = "P, interaction versus non-carriers", quantity = "Quantity", n_women = "Women (n)",
          mean_age_at_menarche = "Mean age at menarche (years)", sd_age_at_menarche = "SD of age at menarche (years)", control_set = "Control set", role = "Role", n_male = "Males (n)", n_female = "Females (n)",
          median_age = "Median age (years)", pct_male = "Percent male", section = "Comparison", n_controls_old = "Controls, old set (n)", est_old = "Estimate, old controls", se_old = "SE, old controls",
          p_old = "P, old controls", n_controls_corrected = "Controls, corrected set (n)", est_corrected = "Estimate, corrected controls", se_corrected = "SE, corrected controls",
          p_corrected = "P, corrected controls", p_family_clustered_corrected = "P, corrected controls, family-clustered", change = "Change", variant = "Run", diff = "Carrier - control (CNV x BMI)",
          se_diff = "Difference SE", p_diff = "P", se_diff_family_clustered = "Difference SE (family-clustered)", p_diff_family_clustered = "P (family-clustered)", rule = "Matching rule",
          measure_ = "Measure", visits = "Visits (n)", correlation = "Correlation", children_f = "Girls (n)", children_m = "Boys (n)", pairs_f = "Visit pairs, girls (n)", pairs_m = "Visit pairs, boys (n)",
          beta_f = "Estimate, girls", se_f = "SE, girls", p_f = "P, girls", beta_m = "Estimate, boys", se_m = "SE, boys", p_m = "P, boys", diff_m_minus_f = "Difference, boys minus girls",
          p_bmi_x_sex = "P, BMI x sex", p_heterogeneity_bmi_x_sex = "P, heterogeneity of BMI x sex across genotypes", cohorts_used = "Cohorts in the model", score_mean = "Mean score in bin (SD)", score_name = "Polygenic score",
          scored_by_00_scoring_only = "Scored by 00_scoring.R only (n)", scored_by_childsds_only = "Scored by childsds only (n)", mean_diff_childsds_minus_00 = "Mean difference, childsds minus 00_scoring.R",
          sd_diff = "SD of the difference", max_abs_diff = "Largest absolute difference", "n_abs_diff_gt_0.1" = "Visits differing by more than 0.1 SD (n)",
          same_sex_child_pairs = "Same-sex child pairs matched (n)", same_sex_cardiff_children = "Cardiff children in a same-sex match (n)", same_sex_chop_children = "CHOP children in a same-sex match (n)",
          opposite_sex_child_pairs_negative_control = "Opposite-sex child pairs matched, negative control (n)", cardiff_22q_children = "Cardiff 22q11.2 children (n)", chop_children = "CHOP children (n)",
          n_controls_role = "Controls (n)")
TOK <- c(n = "n", se = "SE", p = "P", est = "estimate", lo = "lower", hi = "upper", asd = "autism", noasd = "no autism", bmi = "BMI", cnv = "CNV", hc = "head circumference", pgs = "polygenic score",
         het = "heterogeneity", med = "median", diff = "difference", ukbb = "UK Biobank", chop = "CHOP", ssc = "SSC", g2mh = "G2MH", svip = "Simons VIP", cardiff = "Cardiff", x = "x", f = "female", m = "male",
         sd = "SD", pct = "percent", "00" = "00_scoring.R", childsds = "childsds", abs = "absolute", gt = "over", max = "largest", "2plus" = "two or more")
pretty1 <- function(nm) { if (nm %in% names(DICT)) return(DICT[[nm]]); if (grepl("[ (]", nm)) return(nm)
  t <- strsplit(nm, "_", fixed = TRUE)[[1]]; t <- ifelse(t %in% names(TOK), TOK[t], t); s <- paste(t, collapse = " "); paste0(toupper(substr(s, 1, 1)), substring(s, 2)) }
pretty <- function(x) { nn <- vapply(names(x), pretty1, ""); if (anyDuplicated(nn)) stop("duplicate column names after renaming: ", paste(nn[duplicated(nn)], collapse = "; ")); setnames(x, nn); x }
# a column that says the same thing in every row moves to the sub-title of its block
const_out <- function(x, cols = c("note", "status")) { x <- copy(as.data.table(x)); extra <- character(0)
  for (cn in intersect(cols, names(x))) if (uniqueN(x[[cn]]) == 1) { extra <- c(extra, as.character(x[[cn]][1])); x[, (cn) := NULL] }
  list(x = x, extra = paste(extra, collapse = "; ")) }
blk <- function(sub, x, mask = TRUE, file = NULL) { co <- const_out(x); list(sub = paste0(sub, if (nzchar(co$extra)) paste0(" [", co$extra, "]") else ""), tab = pretty(fmt(co$x, mask = mask)), file = file) }

# ---------------------------------------------------------------- workbook writer
wb <- loadWorkbook(wbf3)
put <- function(sh, title, note, blocks, csv) {
  stopifnot(nchar(sh) <= 31)
  if (sh %in% names(wb)) removeWorksheet(wb, sh); addWorksheet(wb, sh)
  writeData(wb, sh, title, startRow = 1); writeData(wb, sh, note, startRow = 2); addStyle(wb, sh, bold, rows = 1, cols = 1)
  r <- 4; W <- numeric(0)
  for (i in seq_along(blocks)) { b <- blocks[[i]]; t <- as.data.table(b$tab)
    if (length(blocks) > 1) { writeData(wb, sh, paste0(LETTERS[i], ". ", b$sub), startRow = r - 1); addStyle(wb, sh, ital, rows = r - 1, cols = 1) }
    writeData(wb, sh, t, startRow = r, headerStyle = bold, keepNA = FALSE)
    w <- vapply(seq_along(t), function(j) max(nchar(names(t)[j]), nchar(as.character(t[[j]])), na.rm = TRUE), 0); length(W) <- max(length(W), length(w)); length(w) <- length(W)
    W <- pmax(W, w, na.rm = TRUE)
    fwrite(t, file.path(ST, sprintf("%s%s.csv", csv, if (length(blocks) > 1) paste0(letters[i], "_", b$file) else "")))
    r <- r + nrow(t) + 3 }
  setColWidths(wb, sh, seq_along(W), pmin(pmax(W + 2, 8), 60))   # widths from the table cells and headers (not the title and note lines)
  cat(sprintf("%-32s %d block(s), rows: %s\n", sh, length(blocks), paste(vapply(blocks, function(b) nrow(b$tab), 0L), collapse = " + ")))
}

# ================================================================ T13 (old table kept; v1 label)
s13 <- "T13 Prospective BMI effects"; t13 <- as.data.table(read.xlsx(wb, s13, sep.names = " "))
put(s13, "Supplementary Table 13. Prospective effect of BMI on growth to the next visit in CNV carriers under 18 years (v1 analysis)",
    paste0("v1 analysis; superseded by Supplementary Table 30. Change in height Z (WHO) or head-circumference Z (UK90) from one visit to the next per SD of BMI Z at the first visit, adjusted for Z at the visit, ",
           "interval, sex, CNV (pooled) and cohort, with a child random intercept; carriers only. Data: Simons VIP, Searchlight and Cardiff (before the CHOP records were added). Table unchanged since 2026-10-24; title and note rows added 2026-11-20."),
    list(list(sub = "", tab = t13)), "supp_table_S13_prospective_v1")

# ================================================================ T14, T15, T16 sex and polygenic background (Jon 2026-11-20: use the 30b run)
N30B <- "Run of 2026-11-14 (30b_sex_pgs_sensitivity_boxcox.R); SRS, SCQ and adult anthropometry on the Box-Cox scale; controls are SSC designated siblings and parents and Simons VIP non-carriers without autism; G2MH carriers and SSC controls with polygenic scores; CHOP records are not in this subset. "
x <- rd(HERE, "supp_table_S4_sex_pgs_strata_boxcox.csv")
setnames(x, c("Stratification", "Stratum", "Trait", "Group", "CNVs in the group", "Carriers (n)", "Carrier BMI slope", "SE (BMI slope)", "P (BMI slope)", "P, age x BMI", "Controls (n)",
              "CNV x BMI vs controls", "SE (CNV x BMI)", "P (CNV x BMI)", "P, CNV x age x BMI", "BMI x stratum, carriers", "SE (BMI x stratum)", "P (BMI x stratum)",
              "Polygenic score main effect, carriers", "P (score main effect)", "BMI x score (continuous), carriers", "SE (BMI x score)", "P (BMI x score)"))
put("T14 Sex and PGS strata", "Supplementary Table 14. BMI effects within strata of sex and of polygenic scores",
    paste0(N30B, "The BMI tests are repeated within males and females and within the lower and upper half of each polygenic score (childhood obesity, adult BMI, height, autism; median split within carriers and within controls). ",
           "Carrier BMI slope: y ~ BMI + smooth of age (+ CNV); CNV x BMI: difference from the controls of the same stratum; BMI x stratum: whether the carrier slope differs between the two strata. ",
           "Simons VIP, Searchlight and Cardiff carriers have no scores and are in the sex strata only. Penn CNB has no controls with BMI. Exploratory; P values are not corrected for multiple testing.", DATE),
    list(list(sub = "", tab = fmt(x))), "supp_table_S14_sex_pgs_strata_30b")
x <- rd(HERE, "supp_table_S5_pgs_adjusted_boxcox.csv")
setnames(x, c("Trait", "Group", "Model", "Covariates added (each with its interaction with BMI)", "Carriers (n)", "Carrier BMI slope", "SE (BMI slope)", "P (BMI slope)", "P, covariate x BMI (single-covariate models)",
              "Controls (n)", "CNV x BMI vs controls", "SE (CNV x BMI)", "P (CNV x BMI)"))
put("T15 PGS adjusted contrasts", "Supplementary Table 15. BMI effects with sex and polygenic scores as covariates",
    paste0(N30B, "For each trait and group, the carrier BMI slope and the CNV x BMI difference from controls are refitted unadjusted, with sex, with each polygenic score, with all four scores, and with all four scores plus sex ",
           "(each covariate with its interaction with BMI), in the people who have all four scores so that every row of a block uses the same individuals.", DATE),
    list(list(sub = "", tab = fmt(x))), "supp_table_S15_pgs_adjusted_30b")
x <- rd(HERE, "supp_table_S6_interaction_power_boxcox.csv")
setnames(x, c("Trait", "Group", "Modifier", "Carriers (n)", "Carrier BMI slope (reference effect)", "Modifier x BMI interaction", "SE (interaction)", "P (interaction)",
              "Smallest interaction detectable with 80% power", "Power for an interaction as large as the BMI slope (proportion)"))
put("T16 Modifier x BMI interactions", "Supplementary Table 16. Modifier x BMI interactions within carriers and the power to detect them",
    paste0(N30B, "Within carriers, y ~ smooth of age + BMI x modifier (+ CNV); modifiers are sex (male versus female) and each polygenic score per SD. Power (normal approximation from the standard error) is for an interaction ",
           "as large as the carriers' own BMI slope. 206 tests; P values are not corrected for multiple testing.", DATE),
    list(list(sub = "", tab = fmt(x))), "supp_table_S16_interaction_power_30b")

# ================================================================ T17 autism diagnosis by CNV (v2)
x <- rd(A, "asd_by_cnv_tests_v2.csv")
has_og <- any(grepl("Overgrowth|Undergrowth", x$group)); cat("T17: overgrowth / undergrowth rows present in the v2 output:", has_og, "\n")
if (!has_og) {   # appended from the v3 sheet and labelled, as instructed
  o <- as.data.table(read.xlsx(loadWorkbook(wbf3), "T17 ASD within CNV", startRow = 4))[grepl("Overgrowth|Undergrowth", group)]
  x <- rbind(x, o[, analysis := "earlier grouping (v3 sheet, 2026-11-19)"], fill = TRUE) }
b17 <- blk("", x); n17 <- b17$sub
put("T17 ASD within CNV", "Supplementary Table 17. Height, head circumference and BMI by autism diagnosis within CNV group",
    paste0("Difference = adjusted mean Z of carriers with an autism diagnosis minus carriers without (SD); model: y ~ autism diagnosis + smooth of age + sex + cohort; when CNVs are pooled, one intercept per cohort within CNV in place of separate cohort and CNV terms (cohort and CNV are partly confounded; the earlier additive model is a sensitivity row). ",
           "Rows: each CNV, all four CNVs pooled, the earlier overgrowth / undergrowth groupings, SSC non-carriers (probands versus designated siblings; the comparison group of the autism x CNV test, not growth controls), ",
           "heterogeneity across the four CNVs and across cohorts (Cochran Q), single-cohort and leave-one-cohort-out fits, under 18 years only, and with Searchlight diagnoses added. ",
           "Autism diagnosis is a secondary outcome. Data: ", CLIN, " (script 85_v2_asd_by_cnv.R).", DATE), list(list(sub = "", tab = b17$tab)), "supp_table_S17_asd_by_cnv_v2")

# ================================================================ T21 BMI slope by age band
a <- rd(G, "supp_table_bmi_slope_by_ageband_v2.csv"); fc <- rd(G, "band_slopes_family_clustered.csv")[variant == "main"]
a[, k := paste(Outcome, `CNV group`, `Age band (years)`)]; fc[, k := paste(outcome, group, age_bin)]; stopifnot(setequal(a$k, fc$k), !anyDuplicated(fc$k))
a <- merge(a, fc[, .(k, `Difference SE (family-clustered)` = round(se_diff_family_clustered, 3), `P (family-clustered)` = signif(p_diff_family_clustered, 2), `Families (n)` = n_families)], by = "k", sort = FALSE)[, k := NULL]
setcolorder(a, c(setdiff(names(a), "Model"), "Model")); t21a <- fmt(a)
nums <- rd(G, "fig6_numbers.csv")
t21b <- nums[quantity == "omnibus_5group_band", .(Outcome = outcome, `Age band (years)` = group, Test = "BMI effect differs among controls and the CNVs (genotype x BMI, F test)", F = round(F, 2), df, P = signif(p, 2),
                                                `CNVs in the model (at least 15 carriers in the band)` = groups, Model = "outcome Z ~ genotype x BMI Z + sex + age, controls and CNVs together")]
i2 <- rd(INF, "infancy_0to2_v2_bands.csv")
t21c <- fmt(i2[, .(Outcome = outcome, `CNV group` = group, `Age band (years)` = age_bin, `Carriers (n)` = n, `Carriers under 4 months (n)` = n_under_4_months, `Median age (years)` = round(median_age, 2),
                   `Carrier BMI slope` = slope, `Carrier SE` = se, P = p, `P, heterogeneity across CNVs (pooled rows)` = het_p, Model = model)])
put("T21 BMI slope by age band", "Supplementary Table 21. BMI effect on height and head circumference within age bands",
    paste0("Slopes are SD of the outcome per SD of BMI. A: one least-squares model per CNV and age band fitted to carriers and controls together (control slope = BMI term; carrier slope = BMI + CNV x BMI); ",
           "'All four CNVs' rows pool the carriers with a CNV term and give the heterogeneity test across CNVs; family-clustered standard errors are the check for relatives in the same family. ",
           "B: whether the BMI effect differs by genotype within each band. C: infancy (0-2 years), carriers only because there are no controls under 2 years; latest visit under 2 years per person; ",
           "head circumference adds the Simons VIP infant visits; rerun on 2026-11-20 with the CHOP infant lengths under 60 cm kept (the earlier version had lost them to a plausibility screen); the 22q11.2 duplication estimate rests on 30 infants. ",
           "Data: ", CLIN, "; children scored against WHO (height, BMI) and UK90 (head circumference) references. ", CTRL, DATE),
    list(list(sub = "Carriers versus controls by age band (2-73 years)", tab = t21a, file = "bmi_slope_by_ageband"),
         list(sub = "BMI effect differs by genotype within age band", tab = t21b, file = "genotype_by_ageband"),
         list(sub = "Infancy (0-2 years), carriers only (rerun 2026-11-20)", tab = t21c, file = "bmi_slope_0to2_carriers")), "supp_table_S21")

# ================================================================ T22 carriers versus controls with UK Biobank adults
bea <- rd(G, "main", "main_with_ukbb_bmi_effect_by_age.csv"); cxb <- unique(bea[, .(outcome, group, cnv_x_bmi = est_main, cnv_x_bmi_se = se_main, cnv_x_bmi_p = p_main)])
t22 <- rbindlist(lapply(c("Height", "Head circumference"), function(o) rbindlist(lapply(c("No CNV", LOCI), function(g) {
  n <- nums[quantity == "n" & outcome == o & group == (if (g == "No CNV") "Controls (tests)" else g)]
  sc <- nums[quantity == "shift_child" & outcome == o & group == g]; sa <- nums[quantity == "shift_adult" & outcome == o & group == g]
  cx <- cxb[outcome == o & group == g]; om <- nums[quantity == "omnibus_contrast" & outcome == o & group == g]; stopifnot(nrow(n) == 1, nrow(cx) == 1)
  data.table(Outcome = o, Group = if (g == "No CNV") "Non-carrier controls" else g, `Clinical cohorts (n)` = n$n_clinical, `UK Biobank (n)` = n$n_ukbb,
             `Shift vs controls, under 18 y (SD)` = if (nrow(sc)) round(sc$est, 3) else NA, `SE, under 18 y` = if (nrow(sc)) round(sc$se, 3) else NA, `P, under 18 y` = if (nrow(sc)) signif(sc$p, 2) else NA,
             `Shift vs controls, 18 y and over (SD)` = if (nrow(sa)) round(sa$est, 3) else NA, `SE, 18 y and over` = if (nrow(sa)) round(sa$se, 3) else NA, `P, 18 y and over` = if (nrow(sa)) signif(sa$p, 2) else NA,
             `BMI slope (controls) or CNV x BMI (carriers), all ages, SD per SD` = round(cx$cnv_x_bmi, 3), SE = round(cx$cnv_x_bmi_se, 3), P = fmt_p(cx$cnv_x_bmi_p),
             `CNV x BMI in the five-group model` = if (nrow(om)) round(om$est, 3) else NA, `SE (five-group)` = if (nrow(om)) round(om$se, 3) else NA, `P (five-group)` = if (nrow(om)) signif(om$p, 2) else NA)
}))))
t22[, P := as.character(P)]
aa <- rd(G, "cnv_x_bmi_all_ages.csv")[variant == "main" & group == "All four CNVs"]
t22b <- fmt(aa[, .(Outcome = outcome, Group = group, `Carriers (n)` = n_carriers, `Controls (n)` = n_controls, `CNV x BMI, all ages (SD per SD)` = est, SE = se, P = p, Model = model)])
tests <- nums[quantity %in% c("omnibus_5group", "omnibus_5group_band", "heterogeneity"),
              .(Outcome = outcome, Test = fcase(quantity == "omnibus_5group", "BMI slope differs among the five genotype groups, all ages (age-adjusted; Wald chi-square)",
                                                quantity == "omnibus_5group_band", paste0("Same, within age band ", group, " y (F; groups: ", groups, ")"),
                                                quantity == "heterogeneity", "BMI slope differs among the four CNVs, carriers only (Wald F)"),
                Statistic = round(fifelse(is.na(F), chisq, F), 2), df = df, P = signif(p, 2))]
ctl <- rd(G, "main", "main_with_ukbb_control_bmi_tests.csv")
tests <- rbind(ctl[, .(Outcome = outcome, Test = "Controls: BMI effect varies with age (age x BMI tensor interaction)", Statistic = NA_real_, df = NA_integer_, P = as.character(fmt_p(p_age_x_bmi)))], tests[, P := as.character(P)])
put("T22 Growth vs controls UKBB", "Supplementary Table 22. Height and head circumference in CNV carriers versus non-carrier controls, with UK Biobank adults",
    paste0("Children: external growth references (WHO height and BMI, UK90 head circumference); adults: Z relative to the non-carrier adults of the same cohort (Box-Cox). ",
           "A: shift = carrier minus control mean (y ~ CNV + smooth of age), separately under and from 18 years; CNV x BMI = age-adjusted difference from controls in the BMI slope over all ages, per CNV and in one five-group model. ",
           "B: all four CNVs pooled, all ages, under three model forms (the form with a CNV term is primary). C: tests of whether the BMI slope differs by genotype. ",
           "Data: ", CLIN, "; UK Biobank adults (16p11.2 and 22q11.2 carriers and 50,000 non-carriers) contribute to height only. ", CTRL, DATE),
    list(list(sub = "Per genotype", tab = t22, file = "growth_tests_ukbb"), list(sub = "All four CNVs pooled, all ages", tab = t22b, file = "pooled_all_ages"),
         list(sub = "Genotype tests", tab = tests, file = "genotype_tests")), "supp_table_S22")

# ================================================================ T30 longitudinal
put("T30 Longitudinal next-visit", "Supplementary Table 30. Longitudinal analysis: BMI at a visit and the change in head circumference and height by the next visit in CNV carriers",
    paste0("Estimate = change in Z to the next visit per SD of BMI Z at the visit. Mixed model per genotype, then all four CNVs with a genotype term, then heterogeneity across genotypes: ",
           "change in Z ~ Z at the visit + interval + sex + spline of age + BMI Z + child random intercept, with a cohort term in the models of one genotype and one intercept per cohort within genotype in the pooled models (cohort and genotype are partly confounded; other adjustments are sensitivity rows); age windows from one model with window x BMI. Carriers only: no non-carrier control has repeated visits. ",
           "Data: head circumference 0.1-18 years from Simons VIP, Cardiff and the CHOP 22q11.2 A-D records; height 2-18 years from Simons VIP, Searchlight, Cardiff and CHOP; all cleaned visits. ",
           "Scoring: WHO height and BMI, UK90 head circumference. Estimates from fewer than 15 children are not reported. Sex tests are exploratory; P values are not corrected for multiple testing. ",
           "H and I count visits and matched records, not people, and are not masked. Script 80_longitudinal_v2.R.", DATE),
    list(blk("Head circumference: main results", rd(L, "results_head_circumference.csv"), file = "results_head_circumference"),
         blk("Height: main results", rd(L, "results_height.csv"), file = "results_height"),
         blk("Separate model per age window", rd(L, "results_by_window_separate_models.csv"), file = "results_by_window_separate_models"),
         blk("BMI x sex (exploratory)", rd(L, "sex_interaction.csv"), file = "sex_interaction"),
         blk("Sensitivity analyses", rd(L, "sensitivity.csv"), file = "sensitivity"),
         blk("Sample, head circumference, by genotype and cohort (counts below 5 masked)", rd(L, "sample_head_circumference_by_genotype_cohort.csv"), file = "sample_head_circumference"),
         blk("Sample, height, by genotype and cohort (counts below 5 masked)", rd(L, "sample_height_by_genotype_cohort.csv"), file = "sample_height"),
         blk("Sample by genotype, cohort and age window (counts below 5 masked)", rd(L, "sample_by_genotype_cohort_window.csv"), file = "sample_by_window"),
         blk("Scoring agreement: childsds versus 00_scoring.R (visits)", rd(L, "scoring_agreement_childsds_vs_00_scoring.csv"), mask = FALSE, file = "scoring_agreement"),
         blk("Duplicate check, CHOP versus Cardiff 22q11.2 children (matched record pairs; opposite-sex pairs are the negative control)", rd(L, "duplicate_check_chop_cardiff.csv"), mask = FALSE, file = "duplicate_check")),
    "supp_table_S30")

# ================================================================ T31 control set
cn <- rd(G, "control_n_by_role_sex_band.csv")
struct0 <- cn$role %in% c("SSC father", "SSC mother")   # fathers have no females and mothers no males by definition: not a small cell
small <- !struct0 & ((!is.na(cn$n_male) & cn$n_male < 5) | (!is.na(cn$n_female) & cn$n_female < 5) | cn$n < 5)
cn2 <- copy(cn)[, c("n", "n_male", "n_female", "pct_male") := .(as.character(n), as.character(n_male), as.character(n_female), as.character(pct_male))]
cn2[small & cn$n < 5, n := "<5"]; cn2[small & cn$n_male < 5, n_male := "<5"]; cn2[small & cn$n_female < 5, n_female := "<5"]
cn2[small, pct_male := ""]; cn2[small & cn$n_male >= 5 & cn$n >= 5 & cn$n_female < 5, n_male := paste0(">", cn$n[small & cn$n_male >= 5 & cn$n >= 5 & cn$n_female < 5] - 5)]
cn2[small & cn$n_female >= 5 & cn$n >= 5 & cn$n_male < 5, n_female := paste0(">", cn$n[small & cn$n_female >= 5 & cn$n >= 5 & cn$n_male < 5] - 5)]
cat("T31: control count rows with a masked cell:", sum(small), "of", nrow(cn), "\n")
put("T31 Control set", "Supplementary Table 31. The non-carrier control set of the growth analyses, and results with the old and the corrected controls",
    paste0("Until 2026-11-20 the child controls of the growth analyses were, by error, the SSC probands (autistic children), because the growth table carried measurements for SSC probands and parents and not for siblings. ",
           "The corrected set takes sibling measurements from the SSC sibling table. A: number, sex and age of the controls by role and age band, old and corrected sets (one visit per person; counts below 5 masked). ",
           "B: every growth result with the old and the corrected controls side by side. Data: ", CLIN, ". ", CTRL, " Scripts 00_control_set.R, 84_control_set_check_v2.R, 85_old_vs_corrected_controls_v2.R.", DATE),
    list(list(sub = "Controls by role, sex and age band", tab = pretty(setnames(fmt(cn2, mask = FALSE), c("n", "band"), c("n_controls_role", "age_bin"))), file = "control_n_by_role_sex_band"),
         blk("Results with the old and the corrected control set", rd(G, "old_vs_corrected_controls.csv"), file = "old_vs_corrected_controls")), "supp_table_S31")

# ================================================================ T32 growth sensitivity
ss <- rd(G, "sensitivity_summary.csv")
VAR <- c(main = "main", visit_latest = "latest visit per person", no_chop = "CHOP records left out", chop_overlap_keep = "G2MH carriers recruited at CHOP kept",
         no_chop_g2mh203_kept = "CHOP left out, G2MH CHOP-site carriers kept", controls_ssc_probands = "old controls (SSC probands)")
STAT <- c(n_carriers = "Carriers (n)", est = "Estimate", se = "SE", p = "P", flag = "Flag")
ss <- fmt(ss); for (v in names(VAR)) for (s in names(STAT)) { o <- paste0(s, "_", v); if (o %in% names(ss)) setnames(ss, o, paste0(STAT[[s]], ", ", VAR[[v]])) }
left <- grep("^(n_carriers|est|se|p|flag)_", names(ss), value = TRUE); if (length(left)) stop("unrenamed sensitivity columns: ", paste(left, collapse = ", "))
put("T32 Growth sensitivity", "Supplementary Table 32. Sensitivity analyses of the cross-sectional growth results",
    paste0("A: each CNV x BMI estimate and genotype test of Supplementary Tables 21 and 22 repeated with the latest visit per person instead of the gap-filling visit rule, without the CHOP records, ",
           "with the G2MH carriers recruited at CHOP kept, and with the old control set (SSC probands); a flag marks an estimate that crosses P = 0.05 or changes sign. ",
           "B: CNV x BMI by age band with standard errors clustered by family. Data: ", CLIN, ". ", CTRL, " Script 83_growth_sensitivity_summary_v2.R.", DATE),
    list(list(sub = "Main result against each sensitivity run", tab = pretty(ss), file = "sensitivity_summary"),
         blk("CNV x BMI by age band with family-clustered standard errors", rd(G, "band_slopes_family_clustered.csv"), file = "band_slopes_family_clustered")), "supp_table_S32")

# ================================================================ T33 sex tests, growth
sx <- rd(G, "supp_table_sex_tests_v2.csv"); nn <- names(sx)
for (j in grep("^SE[0-9]+$", nn)) nn[j] <- paste0("SE (", nn[j - 1], ")"); setnames(sx, nn)
sg <- rd(G, "sex_cnv_bmi_growth.csv")[variant == "main"]
if (nrow(sg) == nrow(sx) && all(sg$outcome == sx$Outcome & sg$group == sx$`CNV group` & sg$scope == sx$`Age range`)) {
  sx[, `SE, three-way (family-clustered)` := round(sg$se_cnv_x_bmi_x_sex_family_clustered, 3)][, `P, three-way (family-clustered)` := signif(sg$p_cnv_x_bmi_x_sex_family_clustered, 2)]
} else cat("T33: family-clustered columns not added (row keys differ)\n")
put("T33 Sex tests growth", "Supplementary Table 33. Sex differences in the BMI effect on height and head circumference (exploratory)",
    paste0("BMI slopes (SD of the outcome per SD of BMI) in male and female carriers and controls, CNV x BMI within each sex, and the three-way term CNV x BMI x sex (male minus female), over all ages and within age bands. ",
           "Exploratory; P values are not corrected for multiple testing. Family-clustered columns are available for the within-band models. Data: ", CLIN, ", plus UK Biobank adults for height. ", CTRL, " Script 82_growth_blocks_v2.R.", DATE),
    list(list(sub = "", tab = fmt(sx))), "supp_table_S33_sex_tests_growth")

# ================================================================ T34 autism diagnosis by BMI
sb <- rd(A, "asd_by_bmi_stratum.csv"); m <- sb$n < 5 | sb$n_asd < 5 | (sb$n - sb$n_asd) < 5
sb2 <- copy(sb)[, c("n", "n_asd") := .(as.character(n), as.character(n_asd))]
sb2[sb$n < 5, n := "<5"]; sb2[sb$n_asd < 5, n_asd := "<5"]; sb2[sb$n >= 5 & sb$n_asd >= 5 & (sb$n - sb$n_asd) < 5, n_asd := paste0(">", sb$n[sb$n >= 5 & sb$n_asd >= 5 & (sb$n - sb$n_asd) < 5] - 5)]
sb2[m, c("prop", "ci_lo", "ci_hi") := NA_real_]; sb2[sb$n < 5, mean_z := NA_real_]
cat("T34: BMI-band rows with a masked cell:", sum(m), "of", nrow(sb), "\n")
b34a <- const_out(sb2)
put("T34 Autism diagnosis by BMI", "Supplementary Table 34. Autism diagnosis by BMI and head circumference in CNV carriers (secondary outcome)",
    paste0("A: share of carriers with an autism diagnosis in the low (30th centile or below), middle and high (70th centile or above) bands of BMI Z and of head-circumference Z, by CNV and cohort, with Wilson 95% intervals; ",
           "cells with fewer than 5 people, or with fewer than 5 with or without a diagnosis, are masked together with their proportions. ",
           "B: logistic regression of diagnosis on the Z score, odds ratio per SD, with age (spline), sex and cohort-within-CNV terms (the earlier additive cohort + CNV model is a sensitivity row); per CNV, pooled, heterogeneity, single-cohort and leave-one-cohort-out fits. Carriers only. ",
           "Autism diagnosis is secondary to the SRS and SCQ scores; how it was recorded differs by cohort. Data: ", CLIN, " (script 85_v2_asd_by_cnv.R).", DATE),
    list(list(sub = paste0("Share with an autism diagnosis by band [", b34a$extra, "]"), tab = pretty(fmt(b34a$x, mask = FALSE)), file = "asd_by_bmi_stratum"),
         blk("Logistic regression of autism diagnosis on BMI Z and head-circumference Z", rd(A, "asd_by_bmi_logistic.csv"), file = "asd_by_bmi_logistic")), "supp_table_S34")

# ================================================================ T35 polygenic score
put("T35 Polygenic score", "Supplementary Table 35. CNV x BMI x autism polygenic score on autistic traits (exploratory)",
    paste0("SRS and SCQ totals are Box-Cox transformed and scaled to SSC unaffected relatives; BMI is a Z score at the questionnaire visit. A: the three-way term CNV x BMI x autism polygenic score and its robustness ",
           "(permutation of the score within cohort and CNV, per CNV, per cohort, leave-one-cohort-out, alternative models). B: carrier and reference BMI slopes within fifths of the score. ",
           "C: logistic regression of autism diagnosis on BMI x score in carriers. D: Mendelian randomization of BMI on autistic traits with adiposity scores as instruments (inconclusive in carriers). ",
           "Data: carriers from G2MH (about four fifths of the genotyped carriers), Simons VIP and SSC with the April 2026 autism score; reference group: SSC unaffected relatives (siblings and parents). ",
           "Exploratory; P values are not corrected for multiple testing. Scripts 86_v2_cnv_bmi_pgs_stats.R and 86a.", DATE),
    list(blk("Three-way term and robustness", rd(A, "threeway_pgs_robustness.csv"), file = "threeway_pgs_robustness"),
         blk("BMI slope by fifth of the polygenic score", setnames(rd(A, "pgs_bmi_slope_by_score_fifth.csv"), "score", "score_mean"), file = "pgs_bmi_slope_by_score_fifth"),
         blk("Autism diagnosis: BMI x score, logistic", rd(A, "asd_dx_bmi_x_pgs_logistic.csv"), file = "asd_dx_bmi_x_pgs_logistic"),
         blk("Mendelian randomization: BMI on autistic traits", rd(A, "mr_bmi_autistic_traits_v2.csv"), file = "mr_bmi_autistic_traits")), "supp_table_S35")

# ================================================================ T36 sex, autism
put("T36 Sex tests autism", "Supplementary Table 36. Sex differences in the BMI effect on autistic traits and autism diagnosis (exploratory)",
    paste0("A: BMI slope per sex, CNV x BMI per sex against same-sex relatives, the three-way term CNV x BMI x sex with its robustness checks, and logistic regression of diagnosis on BMI x sex in carriers. ",
           "B: the main autism results refitted with standard errors clustered by family. SRS and SCQ on the Box-Cox scale. Data: carriers from G2MH, Simons VIP, Searchlight and SSC (diagnosis rows also CHOP, Cardiff and CHU Sainte-Justine); ",
           "reference group: SSC unaffected relatives (siblings and parents). Exploratory; P values are not corrected for multiple testing. Scripts 87_v2_sex_cnv_bmi_autism.R, 87b_v2_family_clustered_check.R.", DATE),
    list(blk("Sex-specific slopes and interaction tests", rd(A, "sex_cnv_bmi_autism.csv"), file = "sex_cnv_bmi_autism"),
         blk("Family-clustered check", rd(A, "family_clustered_check.csv"), file = "family_clustered_check")), "supp_table_S36")

# ================================================================ T37 puberty Mendelian randomization
E <- file.path(PM, "expanse_results")
SCN <- c(pgs_cobes = "childhood obesity score", pgs_cbmi = "childhood BMI score", pgs_bmi = "adult BMI score", pgs_aam = "age-at-menarche score (Guan 2025 weights)", score_perry = "age-at-menarche score (Perry 2014 SNPs)")
mrc <- rd(E, "mr_controls.csv"); mrc[, score := ifelse(score %in% names(SCN), SCN[score], score)]; setnames(mrc, "score", "score_name")
put("T37 Puberty MR", "Supplementary Table 37. UK Biobank: childhood body size, age at menarche and adult height, with Mendelian randomization",
    paste0("All numbers are from UK Biobank (European-ancestry cluster); non-carrier women for the Mendelian randomization (238,681 with age at menarche); 22q11.2 deletion has 5 carriers and is not analysed. ",
           "A: one-sample Mendelian randomization with genome-wide polygenic scores as instruments (childhood adiposity -> age at menarche -> adult height). B: two-sample SNP-level estimates. ",
           "C: negative-control outcomes and other samples. D: the 19 checks of the analysis. E: adult height per step of recalled body size at age 10 and of adult BMI group, by genotype. ",
           "F: age at menarche by genotype and adult height per year of later menarche. Data source: UK Biobank only; scripts 88a and 88b (run on Expanse; only aggregates leave it).", DATE),
    list(blk("One-sample Mendelian randomization (polygenic scores)", rd(E, "mr_onesample.csv"), file = "mr_onesample"),
         blk("Two-sample Mendelian randomization (SNP instruments)", rd(E, "mr_twosample.csv"), file = "mr_twosample"),
         blk("Negative controls and other samples", mrc, file = "mr_controls"),
         blk("Checks", rd(PM, "mr_checks.csv"), file = "mr_checks"),
         blk("Adult height by recalled body size at age 10 and adult BMI group", rd(E, "panel_j_trends.csv"), file = "body_size_trends"),
         blk("Age at menarche by genotype", rd(E, "panel_k_menarche.csv"), file = "menarche")), "supp_table_S37")

# ---------------------------------------------------------------- order and save
if ("T22 Fig6 growth tests UKBB" %in% names(wb)) removeWorksheet(wb, "T22 Fig6 growth tests UKBB")   # replaced by "T22 Growth vs controls UKBB"
ord <- names(wb); num <- suppressWarnings(as.numeric(sub("^T(\\d+).*", "\\1", ord))); num[is.na(num)] <- 999; num[grepl("^T19b", ord)] <- 19.5
# the v3 sheet "T30 Extended CNV models" becomes Supplementary Table 38 (coordinator 2026-11-20): sheet name and title row; content untouched
t38 <- read.xlsx(wb, "T30 Extended CNV models", colNames = FALSE, rows = 1, cols = 1)[[1]]; stopifnot(grepl("^Supplementary Table 30\\. ", t38))
writeData(wb, "T30 Extended CNV models", sub("^Supplementary Table 30\\. ", "Supplementary Table 38. ", t38), startRow = 1, startCol = 1)
renameWorksheet(wb, "T30 Extended CNV models", "T38 Extended CNV models"); ord[ord == "T30 Extended CNV models"] <- "T38 Extended CNV models"; num[ord == "T38 Extended CNV models"] <- 38
worksheetOrder(wb) <- order(num); saveWorkbook(wb, wbf4, overwrite = TRUE)
cat("written:", wbf4, "\n")

# ---------------------------------------------------------------- read back and check
chk <- loadWorkbook(wbf4); sn <- names(chk)[worksheetOrder(chk)]; cat("\nsheets, in order:\n", paste(sn, collapse = " | "), "\n\n")
NEW <- c(s13, "T14 Sex and PGS strata", "T15 PGS adjusted contrasts", "T16 Modifier x BMI interactions", "T38 Extended CNV models", "T17 ASD within CNV", "T21 BMI slope by age band", "T22 Growth vs controls UKBB", "T30 Longitudinal next-visit", "T31 Control set", "T32 Growth sensitivity", "T33 Sex tests growth",
         "T34 Autism diagnosis by BMI", "T35 Polygenic score", "T36 Sex tests autism", "T37 Puberty MR")
for (s in sn) { x <- read.xlsx(wbf4, sheet = s, colNames = FALSE, skipEmptyRows = FALSE); if (is.null(x) || !nrow(x)) stop("empty sheet: ", s)
  if (s %in% NEW) { cat("==", s, "| rows", nrow(x), "cols", ncol(x), "\n"); for (i in 1:5) cat("  ", i, ":", substr(paste(na.omit(unlist(x[i, ])), collapse = " | "), 1, 230), "\n")
    stopifnot(grepl("^Supplementary Table [0-9]+\\. ", x[1, 1]), nchar(x[2, 1]) > 50) } }
