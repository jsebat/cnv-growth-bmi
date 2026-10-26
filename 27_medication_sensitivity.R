# Medication sensitivity analysis (JS / Claude 2026-10-24): can psychotropic medication explain the BMI effects on SRS, SCQ
# and CNB in CNV carriers? Source: G2MH medications release (one row per medication per rarecnv_id; meds_current 1 = current,
# 0 = past, blank = unknown; ATC codes are empty, so classes are assigned from the drug name). Classes: antipsychotic,
# stimulant/ADHD (incl. atomoxetine), ADHD alpha-agonist, antidepressant, antiepileptic/mood stabiliser, anxiolytic/
# benzodiazepine; "weight-gain" = antipsychotic + valproate + lithium + mirtazapine + gabapentinoids; "weight-loss" =
# stimulants + atomoxetine + topiramate + metformin + GLP-1 agonists. Medication status is known for individuals with any
# row in the medications file (an empty row = asked, no medication); others are unknown and excluded from adjusted models.
# For each trait x group: medication prevalence; age-adjusted association of medication class with BMI and with the trait;
# BMI slope unadjusted (same subset), adjusted for classes, adjusted for classes + ASD, and in carriers on no psychotropic.
# Height and head circumference (G2MH carriers, same rows as the main figure) were added on 2026-10-26.
# Output: medication_sensitivity.csv, medication_prevalence.csv; table and write-up in medication_sensitivity_results.md (table rebuilt by medication_table.R).
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE", UNAFF_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
sink(tempfile()); source(file.path(HERE, "23_unaffected_controls.R")); sink()     # via 21: D (SRS/SCQ with asd), cnb (CNB with BMI), rel (release), DB; from 23: A (growth rows, height / HC Z with BMI Z), gmc (GUIDs)
suppressPackageStartupMessages(library(mgcv))
LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")

# ---------------------------------------------------------------- medication classes
m <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_medications_data_release_022026.csv"), colClasses = "character")
m[, nm := tolower(trimws(meds_name))][, current := meds_current == "1"]
RX <- list(
  antipsychotic = "risperid|risperdal|aripiprazol|abilify|quetiapin|seroquel|olanzapin|zyprexa|clozapin|clozaril|pipamperon|haloperidol|haldol|ziprasidon|geodon|lurasidon|latuda|paliperidon|invega|chlorpromazin|levomepromazin|perphenazin|fluphenazin|pimozid|tiaprid|sulpirid|amisulprid|brexpiprazol|rexulti|cariprazin|asenapin|iloperidon|loxapin|thioridazin",
  stimulant     = "methylphenid|concerta|ritalin|medikinet|equasym|vyvanse|lisdexamf|adderall|amphetamin|amfetamin|dexmethylphenid|focalin|dexedrin|elvanse|strattera|atomoxetin|modafinil|qelbree|viloxazin",
  adhd_alpha    = "guanfacin|intuniv|tenex|clonidin|kapvay|catapres",
  antidepressant = "sertralin|zoloft|fluoxetin|prozac|citalopram|celexa|escitalopram|lexapro|cipralex|paroxetin|paxil|fluvoxamin|luvox|duloxetin|cymbalta|venlafaxin|effexor|desvenlafaxin|bupropion|wellbutrin|mirtazapin|remeron|trazodon|amitriptylin|nortriptylin|imipramin|clomipramin|doxepin|vortioxetin|trintellix|buspiron|buspar",
  antiepileptic = "valpro|depakot|depakin|divalpro|lithium|lamotrigin|lamictal|carbamazepin|tegretol|oxcarbazepin|trileptal|levetiracetam|keppra|topiramat|topamax|lacosamid|vimpat|gabapentin|neurontin|pregabalin|lyrica|zonisamid|clobazam|onfi|phenytoin|dilantin|ethosuximid|phenobarb|rufinamid|cenobamat|brivaracetam|perampanel",
  anxiolytic    = "lorazepam|ativan|clonazepam|klonopin|diazepam|valium|oxazepam|alprazolam|xanax|hydroxyzin|atarax|temazepam|propranolol")
for (k in names(RX)) m[, (k) := grepl(RX[[k]], nm)]
m[, weight_gain := antipsychotic | grepl("valpro|depakot|depakin|divalpro|lithium|mirtazapin|remeron|gabapentin|neurontin|pregabalin|lyrica", nm)]
m[, weight_loss := grepl("methylphenid|concerta|ritalin|medikinet|equasym|vyvanse|lisdexamf|adderall|amphetamin|amfetamin|dexmethylphenid|focalin|dexedrin|elvanse|strattera|atomoxetin|topiramat|topamax|metformin|glucophage|semaglutid|ozempic|wegovy|liraglutid|saxenda|victoza|tirzepatid|mounjaro|zepbound|phentermin|naltrexone/bupropion|contrave", nm)]
m[, psychotropic := antipsychotic | stimulant | adhd_alpha | antidepressant | antiepileptic | anxiolytic]
CLS <- c("psychotropic", "antipsychotic", "stimulant", "adhd_alpha", "antidepressant", "antiepileptic", "anxiolytic", "weight_gain", "weight_loss")
cat("medication rows classified (current):\n"); print(m[current == TRUE, lapply(.SD, sum), .SDcols = CLS])
cat("unclassified current names (top 40, for a check that nothing psychotropic was missed):\n")
print(m[current == TRUE & nzchar(nm) & !psychotropic & !weight_loss & !weight_gain, .N, by = nm][order(-N)][1:40], nrows = 50)
# per-person flags: current use (primary) and ever use
Mc <- m[, lapply(.SD, function(v) any(v & current)), .SDcols = CLS, by = rarecnv_id]
Me <- m[, lapply(.SD, function(v) any(v)), .SDcols = CLS, by = rarecnv_id]; setnames(Me, CLS, paste0(CLS, "_ever"))
M <- merge(Mc, Me, by = "rarecnv_id")[, known := TRUE]
M <- merge(M, unique(rel[, .(rarecnv_id, GUID = guid)]), by = "rarecnv_id")
cat("individuals with medication information:", nrow(M), "\n")

# ---------------------------------------------------------------- analysis tables
X <- rbind(D[locus %in% LOC & pheno %in% c("SRS", "SCQ") & !is.na(BMI) & !is.na(y), .(id, locus, age, BMI, y, asd, pheno)],
           cnb[locus %in% LOC & !is.na(BMI) & !is.na(accuracy) & cnb_age >= 4 & cnb_age <= 65, .(id = GUID, locus, age = cnb_age, BMI, y = accuracy, asd = NA_integer_, pheno = "CNB accuracy")],
           cnb[locus %in% LOC & !is.na(BMI) & !is.na(speed) & cnb_age >= 4 & cnb_age <= 65, .(id = GUID, locus, age = cnb_age, BMI, y = speed, asd = NA_integer_, pheno = "CNB speed")])
gid <- unique(gmc[cohort == "G2MH" & !is.na(GUID) & GUID != "", .(id, GUID)])   # growth rows are keyed by genomatrix id; medications by GUID
GRW <- merge(A[cohort == "G2MH" & locus %in% LOC & pheno %in% c("Height", "Head circumference"), .(id, locus, age, BMI, y, asd, pheno)], gid, by = "id")[, .(id = GUID, locus, age, BMI, y, asd, pheno)]
X <- rbind(X, GRW)
X <- merge(X, M[, -"rarecnv_id"], by.x = "id", by.y = "GUID", all.x = TRUE)
X[is.na(known), known := FALSE]
X[, locus := factor(locus)]
# ASD for CNB rows from the G2MH subjects file (same source as 21)
X[pheno %like% "CNB", asd := asdmap[cohort == "G2MH"][match(X[pheno %like% "CNB", id], id), asd]]
cat("\n== medication information coverage by trait\n"); print(X[, .(n = .N, with_med_info = sum(known), pct = round(100 * mean(known))), by = pheno])

# ---------------------------------------------------------------- prevalence and association with BMI / trait
prev <- X[known == TRUE, c(list(n = .N), lapply(.SD, function(v) round(100 * mean(v), 1))), .SDcols = CLS, by = .(pheno, locus)][order(pheno, locus)]
prev_all <- X[known == TRUE, c(list(locus = "All CNVs", n = .N), lapply(.SD, function(v) round(100 * mean(v), 1))), .SDcols = CLS, by = pheno]
prev <- rbind(prev, prev_all, fill = TRUE); fwrite(prev, file.path(HERE, "medication_prevalence.csv"))
cat("\n== % on current medication, carriers with medication information\n"); print(prev, nrows = 100)
assoc <- function(x, cls, out) tryCatch({   # age- and locus-adjusted difference (out: BMI or trait) for medication class, SD units
  x <- copy(x)[, med := as.numeric(get(cls))]; if (sum(x$med) < 5 || sum(x$med == 0) < 5) return(c(NA, NA, NA))
  m <- gam(as.formula(paste(out, "~ med + s(age, k = 8)", if (uniqueN(x$locus) > 1) "+ locus" else "")), data = x, method = "REML"); co <- summary(m)$p.table["med", ]; c(co[1], co[2], co[4]) }, error = function(e) c(NA, NA, NA))
AS <- rbindlist(lapply(unique(X$pheno), function(ph) rbindlist(lapply(c("antipsychotic", "stimulant", "antidepressant", "antiepileptic", "weight_gain", "weight_loss", "psychotropic"), function(cls) {
  x <- X[pheno == ph & known == TRUE]; b <- assoc(x, cls, "BMI"); t <- assoc(x, cls, "y")
  data.table(pheno = ph, class = cls, n_on = sum(x[[cls]]), bmi_diff = b[1], bmi_se = b[2], bmi_p = b[3], trait_diff = t[1], trait_se = t[2], trait_p = t[3]) }))))
cat("\n== medication class vs BMI and vs trait (all carriers with medication info; age + locus adjusted; SD units)\n")
print(AS[, .(pheno, class, n_on, bmi_diff = round(bmi_diff, 2), bmi_p = signif(bmi_p, 2), trait_diff = round(trait_diff, 2), trait_p = signif(trait_p, 2))], nrows = 100)

# ---------------------------------------------------------------- BMI slopes with and without medication adjustment
slope <- function(x, extra = "") tryCatch({
  m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI", if (uniqueN(x$locus) > 1) "+ locus" else "", extra)), data = x, method = "REML"); co <- summary(m)$p.table["BMI", ]
  as.data.table(list(n = nrow(x), slope = unname(co[1]), se = unname(co[2]), p = unname(co[4]))) }, error = function(e) as.data.table(list(n = nrow(x), slope = NA_real_, se = NA_real_, p = NA_real_)))
MEDCOV <- "+ antipsychotic + stimulant + adhd_alpha + antidepressant + antiepileptic + anxiolytic"
rows <- list()
for (ph in unique(X$pheno)) for (g in c(LOC, "All CNVs")) {
  xa <- X[pheno == ph & (if (g == "All CNVs") TRUE else locus == g)]; xk <- xa[known == TRUE]
  if (nrow(xk) < 20) next
  add <- function(model, x, extra = "") rows[[length(rows) + 1]] <<- cbind(data.table(pheno = ph, group = g, model = model), slope(x, extra))
  add("all carriers (as in figures)", xa)
  add("with medication info, unadjusted", xk)
  add("+ any psychotropic", xk, "+ psychotropic")
  add("+ medication classes", xk, MEDCOV)
  add("+ weight-gain + weight-loss meds", xk, "+ weight_gain + weight_loss")
  xo <- xk[psychotropic == FALSE]; if (nrow(xo) >= 20) add("no current psychotropic", xo)
  xw <- xk[weight_gain == FALSE & weight_loss == FALSE]; if (nrow(xw) >= 20) add("no weight-affecting medication", xw)
  xs <- xk[!is.na(asd)]; if (nrow(xs) >= 20 && uniqueN(xs$asd) == 2) { xs[, asd_f := factor(asd)]; add("+ medication classes + ASD", xs, paste(MEDCOV, "+ asd_f")) }
}
S <- rbindlist(rows); fwrite(S, file.path(HERE, "medication_sensitivity.csv"))
cat("\n== BMI slope on trait (SD per SD) by model\n"); print(S[, .(pheno, group, model, n, slope = round(slope, 3), se = round(se, 3), p = signif(p, 2))], nrows = 300)
cat("\ndone\n")
