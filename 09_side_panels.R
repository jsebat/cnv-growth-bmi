# Candidate right-hand panels for the cross-sectional figure (fig_cross_sectional_split):
#  (1) ICV (BrainChart Z) by BMI stratum, carriers < 18 y, same split rule (16p11.2 del > 90th centile vs rest;
#      other loci >= 70th vs <= 30th).  Source: Kuldeep's MRI table (Jon/Claude/MRI_CNV_growth/mri_subjects.csv),
#      first scan per subject (repeat scans dropped).  NB the file holds only the BMI tails (<= 30th and >= 70th), so the
#      16p11.2 del "rest" group is the <= 30th tail plus the 70-90th band, not everyone below the 90th.
#  (2) Head circumference by ASD vs no ASD, combined sample, ages 2-40 (JS: adults added back), one record per person
#      (latest visit). < 18 y: UK90 HC Z from hc_visits.rds (07: full SVIP v10 series + Cardiff/G2MH/SSC);
#      >= 18 y: genomatrix HEADCIRCUM (Molly's within-sample Box-Cox Z), the same scales as the trajectory panels.
#      The adjusted model adds an adult indicator because the two Z scales differ and adults are mostly no-ASD.
#      ASD: G2MH and Cardiff = G2MH psy_asd by GUID (Cardiff participants are enrolled in G2MH); SVIP = v10
#      clinical_asd_dx; SSC carrier probands = ASD (SSC carrier mothers excluded, never assessed).
#      Tests: Welch t (as in a bar plot) and a regression with a cohort term, because ASD rates differ by cohort
#      (SSC 100%, G2MH ~50%, SVIP ~20%).
# Output: side_panel_icv.csv, side_panel_hc_asd.csv, fig_side_panels_preview.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
SPLIT <- c("16p11.2 del" = "30_70", "16p11.2 dup" = "30_70", "22q11.2 del" = "30_70", "22q11.2 dup" = "30_70")   # JS 2026-10-12: 30/70 tails for every locus (was > 90th vs rest for 16p11.2 del)
welch <- function(lo, hi) if (length(lo) >= 2 && length(hi) >= 2) t.test(hi, lo)$p.value else NA_real_
pfmt <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.2g", p)))

# ---------------------------------------------------------------- (1) ICV
mri <- fread(file.path(DB, "Jon/Claude/MRI_CNV_growth/mri_subjects.csv"))
mri <- mri[age < 18][order(source, subject_id, age)][, .SD[1], by = .(source, subject_id)]
mri[, locus := sub("_", " ", sub("DEL", "del", sub("DUP", "dup", group)))]
mri <- mri[locus %in% LOCI]
mri[, grp := ifelse(SPLIT[locus] == "p90", ifelse(bmi_z > qnorm(.9), "high", "low"),
                    ifelse(bmi_z >= qnorm(.7), "high", ifelse(bmi_z <= qnorm(.3), "low", NA)))]
mri <- mri[!is.na(grp)]
icv <- mri[, .(n_low = sum(grp == "low"), n_high = sum(grp == "high"), mean_low = mean(icv_z[grp == "low"]), mean_high = mean(icv_z[grp == "high"]),
               se_low = sd(icv_z[grp == "low"]) / sqrt(sum(grp == "low")), se_high = sd(icv_z[grp == "high"]) / sqrt(sum(grp == "high")),
               p = welch(icv_z[grp == "low"], icv_z[grp == "high"]), sources = paste(unique(source), collapse = "+")), by = locus]
icv[, diff := mean_high - mean_low]; setorder(icv, locus)
fwrite(icv, file.path(HERE, "side_panel_icv.csv")); cat("(1) ICV < 18 y by BMI stratum (first scan per subject)\n"); print(icv, digits = 3)

# ---------------------------------------------------------------- (2) HC by ASD
# head circumference: hc_visits.rds from 07 (SVIP v10 full series + Cardiff/G2MH/SSC, UK90 Z recomputed), latest visit 2-18 y
kid <- readRDS(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/hc_visits.rds"))[age >= 2 & age < 18]
kid <- kid[, .(id, cohort, locus, age, HEADCIRCUM = hcz)]
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
gm <- gm[!(cohort == "G2MH" & Cardiff_mask == TRUE)][, CNV := trimws(CNV)]
gm[, id := ifelse(cohort == "SVIP", paste0("SVIP:", sampleID), paste(cohort, MSassigned_ID, sep = ":"))]
lab4 <- c("16p11.2 proximal DEL" = "16p11.2 del", "16p11.2 proximal DUP" = "16p11.2 dup", "22q11.2 A-D DEL" = "22q11.2 del", "22q11.2 A-D DUP" = "22q11.2 dup")
adult <- gm[CNV %in% names(lab4) & !is.na(HEADCIRCUM) & age >= 18 & age <= 40, .(id, cohort, locus = lab4[CNV], age, HEADCIRCUM)]
g <- rbind(kid, adult)[order(id, -age)][, .SD[1], by = id]
g <- merge(g, unique(gm[, .(id, GUID, sampleID)]), by = "id", all.x = TRUE)
g[cohort == "SVIP", sampleID := sub("^SVIP:", "", id)]
g[, adult := age >= 18]
s <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd"))
g2 <- setNames(s$asd, s$sfari_id)
dx <- fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"),
            colClasses = list(character = "diagnosis_summary.clinical_asd_dx"))
sv <- setNames(tolower(dx$diagnosis_summary.clinical_asd_dx), dx$individual)
g[, asd := NA_integer_]
g[cohort %in% c("G2MH", "Cardiff"), asd := suppressWarnings(as.integer(g2[GUID]))]
g[cohort == "SVIP", asd := ifelse(sv[sampleID] == "true", 1L, ifelse(sv[sampleID] == "false", 0L, NA_integer_))]
g[cohort == "SSC" & grepl("[.]p1$", sampleID), asd := 1L]
cat("
ASD status available, by cohort:
"); print(g[, .(persons = .N, with_asd_status = sum(!is.na(asd)), asd = sum(asd == 1, na.rm = TRUE)), by = cohort])
h <- g[!is.na(asd)]
hc <- h[, {
  lo <- HEADCIRCUM[asd == 0]; hi <- HEADCIRCUM[asd == 1]
  fit <- lm(as.formula(paste("HEADCIRCUM ~ asd + cohort + age + I(age^2)", if (uniqueN(adult) > 1) "+ adult" else "")), data = .SD)
  co <- summary(fit)$coefficients["asd", ]
  .(n_noASD = length(lo), n_ASD = length(hi), mean_noASD = mean(lo), mean_ASD = mean(hi), se_noASD = sd(lo) / sqrt(length(lo)),
    se_ASD = sd(hi) / sqrt(length(hi)), diff = mean(hi) - mean(lo), p_welch = welch(lo, hi), beta_cohort_adj = co[1], p_cohort_adj = co[4],
    cohorts = paste(sort(unique(cohort)), collapse = "+"))
}, by = locus]
setorder(hc, locus)
for (grp in list(c("16p11.2 del", "22q11.2 dup"), c("16p11.2 dup", "22q11.2 del"))) {
  x <- h[locus %in% grp]; fit <- lm(HEADCIRCUM ~ asd + locus + cohort + adult + age + I(age^2), data = x); co <- summary(fit)$coefficients["asd", ]
  lo <- x[asd == 0, HEADCIRCUM]; hi <- x[asd == 1, HEADCIRCUM]
  hc <- rbind(hc, data.table(locus = if (grp[1] == "16p11.2 del") "Overgrowth (16p del + 22q dup)" else "Undergrowth (16p dup + 22q del)",
    n_noASD = length(lo), n_ASD = length(hi), mean_noASD = mean(lo), mean_ASD = mean(hi), se_noASD = sd(lo) / sqrt(length(lo)), se_ASD = sd(hi) / sqrt(length(hi)),
    diff = mean(hi) - mean(lo), p_welch = welch(lo, hi), beta_cohort_adj = co[1], p_cohort_adj = co[4], cohorts = paste(sort(unique(x$cohort)), collapse = "+")))
}
fwrite(hc, file.path(HERE, "side_panel_hc_asd.csv")); cat("\n(2) Head circumference by ASD (one record per person, latest visit, ages 2-40)\n"); print(hc, digits = 3)
cat('
children vs adults with ASD status:
'); print(g[!is.na(asd), .(n = .N, asd = sum(asd)), by = .(locus, adult)][order(locus, adult)])
cat("\nHC by ASD within cohort (16p11.2 del):\n"); print(h[locus == "16p11.2 del", .(n = .N, n_asd = sum(asd), diff = mean(HEADCIRCUM[asd == 1]) - mean(HEADCIRCUM[asd == 0]),
                                                                                   p = welch(HEADCIRCUM[asd == 0], HEADCIRCUM[asd == 1])), by = cohort], digits = 3)

# ---------------------------------------------------------------- preview figure (bars, one column per panel type)
bars <- function(dt, lo_col, hi_col, lo_se, hi_se, n_lo, n_hi, labs, cols, ylab, title, pcol, plab) {
  d <- rbind(dt[, .(locus, grp = labs[1], m = get(lo_col), se = get(lo_se), n = get(n_lo))], dt[, .(locus, grp = labs[2], m = get(hi_col), se = get(hi_se), n = get(n_hi))])
  d[, grp := factor(grp, labs)]; d[, locus := factor(locus, unique(dt$locus))]
  pl <- dt[, .(locus = factor(locus, unique(dt$locus)), lab = plab(.SD))]
  ggplot(d, aes(locus, m, fill = grp)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(y = -Inf, label = paste0("n=", n)), position = position_dodge(0.8), vjust = -0.4, size = 2.4) +
    geom_text(data = pl, aes(locus, Inf, label = lab), inherit.aes = FALSE, vjust = 1.2, size = 2.6) +
    scale_fill_manual(values = cols) + labs(x = NULL, y = ylab, fill = NULL, title = title) +
    theme_classic(base_size = 10) + theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
}
p1 <- bars(icv, "mean_low", "mean_high", "se_low", "se_high", "n_low", "n_high", c("Lower BMI stratum", "Higher BMI stratum"), c("#00A5AD", "#F0645A"),
           "ICV Z (BrainChart)", "(1) ICV by BMI stratum, < 18 y", "p", function(x) sprintf("\u0394 %+.2f\n%s", x$diff, pfmt(x$p)))
p2 <- bars(hc, "mean_noASD", "mean_ASD", "se_noASD", "se_ASD", "n_noASD", "n_ASD", c("No ASD", "ASD"), c("#8C9A9A", "#F28E2B"),
           "Head circumference Z", "(2) Head circumference by ASD, ages 2-40", "p",
           function(x) sprintf("\u0394 %+.2f, %s\nadj. %+.2f, %s", x$diff, pfmt(x$p_welch), x$beta_cohort_adj, pfmt(x$p_cohort_adj)))
ggsave(file.path(HERE, "fig_side_panels_preview.png"), p1 + p2 + plot_layout(widths = c(1, 1.5)), width = 12, height = 5, dpi = 200)
