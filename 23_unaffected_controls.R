# Height, head circumference and behavioural / cognitive trajectories by BMI, with an UNAFFECTED control set (JS 2026-10-24, Jon).
# Controls (no ASD, no CNV): SSC designated siblings (height, weight, head circumference from the SSC sibling ssc_commonly_used
#   table; age/sex from Age_individual_v14_parent_age.txt), SSC parents (growth genomatrix), SVIP non-carrier relatives without a
#   clinical ASD diagnosis (v10 diagnosis_summary). SSC probands are excluded everywhere. No UK Biobank; runs locally.
# Carriers: 16p11.2 BP4-BP5 del/dup, 22q11.2 A-D del/dup, clinical cohorts (G2MH, Cardiff, SVIP, SSC; Searchlight for height only
#   as in the main figure), latest visit with the outcome and BMI, ages 2-65. ASD status as in panels F/L.
# Outcomes: height Z and head circumference Z (00_scoring.R), and the phenotypes of 21_pheno_bmi_gam.R (its data build is
#   sourced): SRS and SCQ have unaffected SSC controls; RBS-R, FSIQ and the CNB exist only in carriers (SSC has them for probands
#   only, G2MH non-carriers have no BMI), so they are carriers-only.
# Tests (REML; bam for the control contrasts): within carriers, BMI effect (y ~ BMI + s(age)) and age x BMI (ti); vs unaffected
#   controls, the age-adjusted CNV x BMI shift and CNV x age x BMI, for all carriers and for non-ASD carriers; per locus and pooled
#   (overgrowth, undergrowth, all) with a locus covariate. Control rows: BMI effect and age x BMI within the unaffected set.
# Figures: fig_unaffected_controls.{png,pdf} (Height, HC, SRS, SCQ; unaffected controls + 4 CNVs) and
#          fig_carriers_only_trajectories.{png,pdf} (SRS, SCQ, RBS-R, FSIQ, CNB accuracy, CNB speed; 4 CNVs, no controls).
# Table: unaffected_controls_tests.csv
Sys.setenv(PHENO_RUN_TESTS = "FALSE", PHENO_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "21_pheno_bmi_gam.R"))   # D (phenotypes), gm, rd(), bmiZ(), score_*(), pfx, LOCI, lab4, agef, SSC, DB, DATA
t0 <- Sys.time()

# ---------------------------------------------------------------- ASD status
s_asd <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd"))
g2asd <- setNames(as.integer(s_asd$asd), s_asd$sfari_id)
dx <- fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"), colClasses = list(character = "diagnosis_summary.clinical_asd_dx"))
sv <- setNames(tolower(dx$diagnosis_summary.clinical_asd_dx), dx$individual)
asd_of <- function(cohort, sampleID, GUID) fifelse(cohort %in% c("G2MH", "Cardiff"), suppressWarnings(as.integer(g2asd[GUID])),
                                                 fifelse(cohort == "SVIP", fifelse(sv[sampleID] == "true", 1L, fifelse(sv[sampleID] == "false", 0L, NA_integer_)),
                                                 fifelse(cohort == "SSC", as.integer(grepl("[.]p[0-9]$", sampleID)), NA_integer_)))

# ---------------------------------------------------------------- growth outcomes: carriers and unaffected controls
gmc <- gm[!(cohort == "G2MH" & Cardiff_mask %in% TRUE) & !is.na(BMI_raw) & !is.na(age) & age >= 2 & age <= 65]
gmc[, id := ifelse(cohort == "SVIP", paste0("SVIP:", sampleID), paste(cohort, ifelse(is.na(MSassigned_ID) | MSassigned_ID == "", sampleID, MSassigned_ID), sep = ":"))]
gmc[, asd := asd_of(cohort, sampleID, GUID)]
# SSC siblings (not in the genomatrix with measurements): height, weight, head circumference from the SSC sibling table
cu <- rd("Designated Unaffected Sibling Data", "ssc_commonly_used.csv")[, .(sampleID = individual, h = suppressWarnings(as.numeric(height)), w = suppressWarnings(as.numeric(weight)), hc = suppressWarnings(as.numeric(head_circumference)))]
sibage <- sibage[, .(sampleID = id, sex, age)]   # sibling age at the proband evaluation visit (see 21_pheno_bmi_gam.R)
sib <- merge(cu, sibage, by = "sampleID")[sampleID %in% gm[cohort == "SSC" & CNV == "No CNV", sampleID] & !is.na(age) & age >= 2]
sib <- sib[, .(id = paste0("SSC:", sampleID), cohort = "SSC", sampleID, GUID = NA_character_, CNV = "No CNV", sex, age, HEIGHT_raw = h, BMI_raw = w / (h / 100)^2, HEADCIRCUM_raw = hc, asd = 0L)]
cat("SSC siblings with BMI:", sum(!is.na(sib$BMI_raw)), " height:", sum(!is.na(sib$HEIGHT_raw)), " HC:", sum(!is.na(sib$HEADCIRCUM_raw)), "\n")
ctl_g <- rbind(gmc[CNV == "No CNV" & ((cohort == "SSC" & !grepl("[.]p[0-9]$", sampleID)) | (cohort == "SVIP" & asd %in% 0L)), .(id, cohort, sampleID, GUID, CNV, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd)], sib[!is.na(BMI_raw)])
ctl_g[, locus := "Unaffected controls"]
car_g <- gmc[CNV %in% names(lab4), .(id, cohort, sampleID, GUID, CNV, sex, age, HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd)][, locus := lab4[CNV]]
mk_growth <- function(x, raw, ph) {
  x <- x[!is.na(get(raw))][if (ph == "Height") cohort != "Searchlight" | TRUE else cohort != "Searchlight"]   # Searchlight has no HC
  x <- x[order(id, -age)][, .SD[1], by = id]
  x[, y := if (ph == "Height") score_height(get(raw), age, sex) else score_hc(get(raw), age, sex)][, BMI := bmiZ(BMI_raw, age, sex)]
  x[!is.na(y) & !is.na(BMI), .(id, cohort, locus, sex, age, BMI, raw = get(raw), y, pheno = ph, asd, unaffected = locus == "Unaffected controls", role = ifelse(locus == "Unaffected controls", "unaffected", "carrier"))]
}
G <- rbind(mk_growth(rbind(ctl_g, car_g), "HEIGHT_raw", "Height"), mk_growth(rbind(ctl_g, car_g), "HEADCIRCUM_raw", "Head circumference"))

# ---------------------------------------------------------------- phenotypes: carriers + unaffected SSC controls (no probands)
P <- D[locus != "No CNV" | unaffected == TRUE]
P[locus == "No CNV", locus := "Unaffected controls"]
P <- P[, .(id, cohort, locus, sex, age, BMI, raw, y, pheno, asd, unaffected, role)]   # asd already attached by 21_pheno_bmi_gam.R (G2MH subjects.csv; SSC probands = 1)
A <- rbind(G, P)
A[, tail := fifelse(BMI >= qnorm(.7), "BMI \u2265 70th centile", fifelse(BMI <= qnorm(.3), "BMI \u2264 30th centile", NA_character_))]
PHEN_C <- c("Height", "Head circumference", "SRS", "SCQ")                           # have unaffected controls
PHEN_O <- c("SRS", "SCQ", "RBS-R", "FSIQ", "CNB accuracy", "CNB speed")             # carriers-only figure
cat("\n== sample (unaffected controls and carriers)\n")
print(A[, .(n = .N, n_asd = sum(asd %in% 1L), n_noasd = sum(asd %in% 0L), hi = sum(tail == "BMI \u2265 70th centile", na.rm = TRUE), lo = sum(tail == "BMI \u2264 30th centile", na.rm = TRUE), age_med = as.numeric(median(age))), by = .(pheno, locus)][order(pheno, locus)], nrows = 100)

# ---------------------------------------------------------------- tests
POOLS <- list("Overgrowth CNVs" = c("16p11.2 del", "22q11.2 dup"), "Undergrowth CNVs" = c("16p11.2 dup", "22q11.2 del"), "All CNVs" = LOCI)
cov_locus <- function(x) if (uniqueN(x$locus[x$locus != "Unaffected controls"]) > 1) "locus + " else ""
within_p <- function(x) tryCatch({ f <- as.formula(paste("y ~", cov_locus(x), "s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6))")); summary(gam(f, data = x, method = "REML"))$s.table["ti(age,BMI)", "p-value"] }, error = function(e) NA_real_)
bmi_main <- function(x) tryCatch({ f <- as.formula(paste("y ~", cov_locus(x), "BMI + s(age, k = 8)")); co <- summary(gam(f, data = x, method = "REML"))$p.table["BMI", ]; c(co[1], co[2], co[4]) }, error = function(e) rep(NA_real_, 3))
contrast <- function(ctl, xc) tryCatch({   # all with bam; returns cnv main (z, p), CNV x BMI shift (est, se, p), CNV x age x BMI p
  if (nrow(xc) < 30 || nrow(ctl) < 30) stop()
  z <- rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == "Unaffected controls", "ctl", "CNV"), levels = c("ctl", "CNV"))][, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  cl <- cov_locus(z)
  m1 <- bam(as.formula(paste("y ~", cl, "is_cnv + s(age, k = 10, by = is_cnv) + s(age, k = 8, by = BMI) + cnvBMI")), data = z, method = "fREML"); co <- summary(m1)$p.table
  m2 <- bam(as.formula(paste("y ~", cl, "is_cnv + s(age, k = 10) + s(BMI, k = 6) + s(age, k = 10, by = is_cnv) + ti(age, BMI, k = c(6, 6)) + ti(age, BMI, k = c(6, 6), by = is_cnv)")), data = z, method = "fREML")
  st <- summary(m2)$s.table
  c(co["is_cnvCNV", 1], co["is_cnvCNV", 4], co["cnvBMI", 1], co["cnvBMI", 2], co["cnvBMI", 4], st[grep("^ti\\(age,BMI\\):is_cnv", rownames(st))[1], "p-value"]) }, error = function(e) rep(NA_real_, 6))
test_group <- function(ph, gname, loci_in) {
  xc <- A[pheno == ph & locus %in% loci_in]; ctl <- A[pheno == ph & locus == "Unaffected controls"]
  base <- data.table(pheno = ph, group = gname, n_cnv = nrow(xc), n_cnv_asd = sum(xc$asd %in% 1L), n_cnv_noasd = sum(xc$asd %in% 0L),
                     n_cnv_hi = sum(xc$tail == "BMI \u2265 70th centile", na.rm = TRUE), n_cnv_lo = sum(xc$tail == "BMI \u2264 30th centile", na.rm = TRUE), n_ctl = nrow(ctl))
  if (nrow(xc) < 30) return(base[, comparison := "within carriers"])
  bm <- bmi_main(copy(xc)); out <- copy(base)[, `:=`(comparison = "within carriers: BMI effect and age x BMI", carrier_bmi_slope = bm[1], carrier_bmi_se = bm[2], carrier_bmi_p = bm[3], within_age_x_bmi_p = within_p(copy(xc)))]
  xn <- xc[asd %in% 0L]; if (nrow(xn) >= 30) { b <- bmi_main(copy(xn)); out <- rbind(out, copy(base)[, `:=`(comparison = "within non-ASD carriers: BMI effect and age x BMI", n_cnv = nrow(xn), carrier_bmi_slope = b[1], carrier_bmi_se = b[2], carrier_bmi_p = b[3], within_age_x_bmi_p = within_p(copy(xn)))], fill = TRUE) }
  if (nrow(ctl) >= 30) {
    e <- contrast(ctl, xc); out <- rbind(out, copy(base)[, `:=`(comparison = "all carriers vs unaffected controls", cnv_main_z = e[1], cnv_main_p = e[2], cnv_x_bmi = e[3], cnv_x_bmi_se = e[4], cnv_x_bmi_p = e[5], cnv_x_age_x_bmi_p = e[6])], fill = TRUE)
    if (nrow(xn) >= 30) { en <- contrast(ctl, xn); out <- rbind(out, copy(base)[, `:=`(comparison = "non-ASD carriers vs unaffected controls", n_cnv = nrow(xn), cnv_main_z = en[1], cnv_main_p = en[2], cnv_x_bmi = en[3], cnv_x_bmi_se = en[4], cnv_x_bmi_p = en[5], cnv_x_age_x_bmi_p = en[6])], fill = TRUE) }
  }
  cat(ph, gname, "done", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "min\n"); out
}
ALLPH <- unique(A$pheno)
res <- rbindlist(fill = TRUE, lapply(ALLPH, function(ph) rbind(rbindlist(fill = TRUE, lapply(LOCI, function(L) test_group(ph, L, L))), rbindlist(fill = TRUE, lapply(names(POOLS), function(Pn) test_group(ph, Pn, POOLS[[Pn]]))))))
ctl_rows <- rbindlist(lapply(PHEN_C, function(ph) { x <- A[pheno == ph & locus == "Unaffected controls"]; b <- bmi_main(copy(x))
  data.table(pheno = ph, group = "Unaffected controls", comparison = "within unaffected controls: BMI effect and age x BMI", n_ctl = nrow(x), carrier_bmi_slope = b[1], carrier_bmi_se = b[2], carrier_bmi_p = b[3], within_age_x_bmi_p = within_p(copy(x))) }))
res <- rbind(res, ctl_rows, fill = TRUE)
setcolorder(res, intersect(c("pheno", "group", "comparison", "n_cnv", "n_cnv_asd", "n_cnv_noasd", "n_cnv_hi", "n_cnv_lo", "n_ctl", "cnv_main_z", "cnv_main_p", "cnv_x_bmi", "cnv_x_bmi_se", "cnv_x_bmi_p", "cnv_x_age_x_bmi_p", "carrier_bmi_slope", "carrier_bmi_se", "carrier_bmi_p", "within_age_x_bmi_p"), names(res)))
fwrite(res, file.path(HERE, "unaffected_controls_tests.csv"))
cat("\n== tests\n"); print(res[, .(pheno, group = substr(group, 1, 16), comparison = substr(comparison, 1, 40), n_cnv, n_ctl, bmi = round(carrier_bmi_slope, 3), p_bmi = signif(carrier_bmi_p, 2), p_age_bmi = signif(within_age_x_bmi_p, 2), cxb = round(cnv_x_bmi, 3), p_cxb = signif(cnv_x_bmi_p, 2), p_3way = signif(cnv_x_age_x_bmi_p, 2))], nrows = 300)

# ---------------------------------------------------------------- figures
COLS <- c("BMI \u2265 70th centile" = "#F0645A", "BMI \u2264 30th centile" = "#00A5AD")
YLAB <- c(Height = "Height (SD)", `Head circumference` = "Head circumference (SD)", SRS = "SRS total (Z vs SSC unaffected)", SCQ = "SCQ lifetime (Z vs SSC siblings)", `RBS-R` = "RBS-R 29 items (Z vs SSC probands)", FSIQ = "FSIQ (Z, norm 100/15)", `CNB accuracy` = "CNB accuracy (Z vs G2MH non-carriers)", `CNB speed` = "CNB speed (Z, higher = faster)")
pick <- function(ph, Gn, comp) res[pheno == ph & group == Gn & startsWith(comparison, comp)][1]
panel <- function(ph, Gn, show_title, mode) {
  x <- A[pheno == ph & locus == Gn & !is.na(tail)][, grp := factor(tail, names(COLS))]; ctl <- Gn == "Unaffected controls"
  if (mode == "controls") { r <- pick(ph, Gn, "all carriers vs unaffected"); lab <- if (ctl) "" else gsub("= < ", "< ", sprintf("CNV × BMI p = %s\nCNV × age × BMI p = %s", pfx(r$cnv_x_bmi_p), pfx(r$cnv_x_age_x_bmi_p)), fixed = TRUE); pv <- c(r$cnv_x_bmi_p, r$cnv_x_age_x_bmi_p) }
  else { r <- pick(ph, Gn, "within carriers"); lab <- gsub("= < ", "< ", sprintf("BMI p = %s\nage × BMI p = %s", pfx(r$carrier_bmi_p), pfx(r$within_age_x_bmi_p)), fixed = TRUE); pv <- c(r$carrier_bmi_p, r$within_age_x_bmi_p) }
  ll <- if (ctl) c("", "") else strsplit(lab, "\n", fixed = TRUE)[[1]]; b1 <- !ctl && isTRUE(pv[1] < 0.05); b2 <- !ctl && isTRUE(pv[2] < 0.05)   # each p-value line bold only if it is itself < 0.05
  yl <- range(quantile(A[pheno == ph, y], c(0.005, 0.995))); yl <- yl + c(-0.6, 0.9) * diff(yl) * 0.2
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = if (ctl) 0.5 else 1.2, alpha = if (ctl) 0.15 else 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 6), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) +
    annotate("text", x = 65, y = yl[2], label = ll[1], hjust = 1, vjust = 1, size = 3, fontface = if (b1) "bold" else "plain") +
    annotate("text", x = 65, y = yl[2] - 0.075 * diff(yl), label = ll[2], hjust = 1, vjust = 1, size = 3, fontface = if (b2) "bold" else "plain") +
    scale_colour_manual(values = COLS, breaks = names(COLS)) + scale_fill_manual(values = COLS, breaks = names(COLS)) + coord_cartesian(xlim = c(2, 65), ylim = yl) +
    labs(x = "Age (years)", y = YLAB[[ph]], colour = NULL, fill = NULL, title = if (show_title) ph else NULL, subtitle = sprintf("%s (n = %d)", if (ctl) "Unaffected (SSC sibs/parents, SVIP)" else Gn, nrow(x))) +
    theme_classic(base_size = 13) + theme(legend.text = element_text(size = 7), legend.position = "inside", legend.position.inside = c(0.22, 0.1), legend.background = element_blank(), legend.key.height = unit(7, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
}
if (toupper(Sys.getenv("UNAFF_DRAW_FIG", "TRUE")) != "FALSE") {   # 25_main_figure_v3.R sources this file for A, res and panel() and skips the figures
GR1 <- c("Unaffected controls", LOCI)
p1 <- list(); for (ph in PHEN_C) for (Gn in GR1) p1[[paste(ph, Gn)]] <- panel(ph, Gn, Gn == GR1[1], "controls")
f1 <- wrap_plots(p1, ncol = length(PHEN_C), byrow = FALSE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(HERE, "fig_unaffected_controls.png"), f1, width = 20, height = 20, dpi = 140)
tryCatch(ggsave(file.path(HERE, "fig_unaffected_controls.pdf"), f1, width = 20, height = 20, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_unaffected_controls.pdf"), f1, width = 20, height = 20))
p2 <- list(); for (ph in PHEN_O) for (Gn in LOCI) p2[[paste(ph, Gn)]] <- panel(ph, Gn, Gn == LOCI[1], "carriers")
f2 <- wrap_plots(p2, ncol = length(PHEN_O), byrow = FALSE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(HERE, "fig_carriers_only_trajectories.png"), f2, width = 28, height = 16, dpi = 130)
tryCatch(ggsave(file.path(HERE, "fig_carriers_only_trajectories.pdf"), f2, width = 28, height = 16, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_carriers_only_trajectories.pdf"), f2, width = 28, height = 16))
}
cat("done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "min\n")
