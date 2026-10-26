# Penn CNB by CNV group and by BMI (JS 2026-10-23). Companion to panels F and L of the main growth figure.
#
# Data: Anjali's processed G2MH release, NetworkPaper/AnalysisByAnjali/data/cnb_processing/cnb_normalized.csv
#   (process_cnb.R): every CNB score set to NA unless the test is valid ("V"), then z-scored against the G2MH
#   non-carriers (lm on age + age^2 + sex + site fitted in no_cnv; everyone scaled by the control residual SD).
#   One CNB session per person. Groups = final_group_class (16p11.2 BP4-BP5 del/dup, 22q11.2 A-D del/dup, no_cnv);
#   people with two CNVs or a triplication are excluded.
# Composites (per person, >= 5 contributing tests):
#   accuracy = mean z of the accuracy measures (PCET acc, CPT true-positive rate, LNB correct, ER40, PMAT, VOLT, CPF,
#              MEDF, ADT, PLOT); higher = better
#   speed    = -(mean z of the median correct-response RTs of MPRAXIS, PCET, CPT, LNB, ER40, PMAT, VOLT, CPF, MEDF,
#              ADT, PLOT); higher = faster
# CNV groups as in F/L: non-carriers; "overgrowth" = 16p11.2 del + 22q11.2 dup; "undergrowth" = 16p11.2 dup + 22q11.2
#   del; plus each locus. Group effect = mean difference vs non-carriers (Welch) and lm adjusted for CNB age, age^2, sex.
# BMI: raw height/weight from (1) the growth genomatrix visit nearest the CNB age (G2MH rows, carriers), else (2) the
#   release medical form (med_height, med_weight at dem_age; the only source for non-carriers). BMI Z by 00_scoring.R
#   (WHO < 19 y; adults: centile among clinical non-carrier adults). Tails: >= 70th vs <= 30th centile, as in the main
#   figure. BMI effect = high - low tail difference within each CNV group (Welch, and lm adjusted for age, age^2, sex),
#   plus continuous BMI Z slope, and the CNV x BMI interaction vs non-carriers.
# Output: cnb_results.csv, cnb_trajectory_tests.csv, fig_cnb_trajectories.{png,pdf}; cnb_hc_results.csv, fig_cnb_by_hc.{png,pdf} (HC section at the end)
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
DB   <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")
source(file.path(HERE, "00_scoring.R"))
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))
welch <- function(a, b) if (length(a) > 1 && length(b) > 1) t.test(a, b)$p.value else NA_real_

# ---------------------------------------------------------------- CNB
cnb <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/cnb_processing/cnb_normalized.csv"))
LAB <- c(no_cnv = "Non-carriers", `16p11.2_Deletion` = "16p11.2 del", `16p11.2_Duplication` = "16p11.2 dup",
         `22q11.2_Deletion` = "22q11.2 del", `22q11.2_Duplication` = "22q11.2 dup")
cnb <- cnb[final_group_class %in% names(LAB) & !is.na(cnb_age)][, locus := LAB[final_group_class]]
ACC <- c(PCET = "cnb_pcet_acc2_z", CPT = "cnb_cpt_ptp_z", LNB = "cnb_lnb_mcr_z", ER40 = "cnb_er40_cr_z", PMAT = "cnb_pmat_pc_z",
         VOLT = "cnb_volt_cr_z", CPF = "cnb_cpf_cr_z", MEDF = "cnb_medf_pc_z", ADT = "cnb_adt_pc_z", PLOT = "cnb_plot_pc_z")
SPD <- c(MPRAXIS = "cnb_mpraxis_rtcr_z", PCET = "cnb_pcet_rtcr_z", CPT = "cnb_cpt_tprt_z", LNB = "cnb_lnb_mrtc_z", ER40 = "cnb_er40_rtcr_z",
         PMAT = "cnb_pmat_rtcr_z", VOLT = "cnb_volt_rtcr_z", CPF = "cnb_cpf_rtcr_z", MEDF = "cnb_medf_rtcr_z", ADT = "cnb_adt_rtcr_z", PLOT = "cnb_plot_rtcr_z")
cnb[, n_acc := rowSums(!is.na(.SD)), .SDcols = ACC][, n_spd := rowSums(!is.na(.SD)), .SDcols = SPD]
cnb[, accuracy := ifelse(n_acc >= 5, rowMeans(.SD, na.rm = TRUE), NA_real_), .SDcols = ACC]
cnb[, speed    := ifelse(n_spd >= 5, -rowMeans(.SD, na.rm = TRUE), NA_real_), .SDcols = SPD]
cnb[, sex := ifelse(sex_final == "Male", "male", ifelse(sex_final == "Female", "female", NA))]
cnb[, pooled := fcase(locus == "Non-carriers", "Non-carriers", locus %in% c("16p11.2 del", "22q11.2 dup"), "Overgrowth CNVs",
                      locus %in% c("16p11.2 dup", "22q11.2 del"), "Undergrowth CNVs")]
cat("CNB subjects by locus:\n"); print(cnb[, .(n = .N, acc = sum(!is.na(accuracy)), spd = sum(!is.na(speed)), age_med = as.numeric(median(cnb_age))), by = locus])

# ---------------------------------------------------------------- BMI
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
g2 <- gm[cohort == "G2MH" & !(Cardiff_mask %in% TRUE) & !is.na(GUID) & GUID != "" & !is.na(BMI_raw), .(GUID, age, BMI_raw)]
v <- merge(cnb[, .(GUID, cnb_age)], g2, by = "GUID", allow.cartesian = TRUE)[, dage := abs(age - cnb_age)][order(GUID, dage)][, .SD[1], by = GUID]
v <- v[dage <= 3]                                                   # visit within 3 y of the CNB session
cnb <- merge(cnb, v[, .(GUID, bmi_age = age, bmi_raw = BMI_raw, bmi_src = "genomatrix visit")], by = "GUID", all.x = TRUE)
rel <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"),
             select = c("rarecnv_id", "dem_age", "med_height_m", "med_height_cm", "med_weight"))
rel[, h := suppressWarnings(as.numeric(med_height_m)) * 100 + suppressWarnings(as.numeric(med_height_cm))][, w := suppressWarnings(as.numeric(med_weight))]
rel <- rel[!is.na(h) & !is.na(w) & h > 50 & h < 230 & w > 5 & w < 300 & !is.na(dem_age), .(rarecnv_id, rel_age = as.numeric(dem_age), rel_bmi = w / (h / 100)^2)]
cnb <- merge(cnb, rel, by = "rarecnv_id", all.x = TRUE)
cnb[is.na(bmi_raw) & !is.na(rel_bmi) & abs(rel_age - cnb_age) <= 3, `:=`(bmi_age = rel_age, bmi_raw = rel_bmi, bmi_src = "release medical form")]
cat("BMI source:\n"); print(cnb[, .N, by = .(locus, bmi_src)][order(locus, bmi_src)])
REFad <- readRDS(file.path(DATA, "analysis_table.rds"))[locus == "No CNV" & age >= 19 & !is.na(BMI_raw)][order(id, -age)][, .SD[1], by = id]
cnb[, bmiZ := score_bmi(bmi_raw, bmi_age, sex, REFad$BMI_raw, REFad$age, as.character(REFad$sex))]
cnb[, tail := fifelse(bmiZ >= qnorm(.7), "High BMI", fifelse(bmiZ <= qnorm(.3), "Low BMI", NA_character_))]
cat("BMI tails by locus:\n"); print(dcast(cnb[!is.na(tail) & !is.na(accuracy), .N, by = .(locus, tail)], locus ~ tail, value.var = "N", fill = 0))

# ---------------------------------------------------------------- group contrasts vs non-carriers
GROUPS <- c("Non-carriers", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup", "Overgrowth CNVs", "Undergrowth CNVs")
grp_rows <- function(y) rbindlist(fill = TRUE, lapply(GROUPS, function(G) {
  d <- cnb[!is.na(get(y))]; d[, yy := get(y)]
  x <- if (G %in% c("Overgrowth CNVs", "Undergrowth CNVs")) d[pooled == G] else d[locus == G]
  ctl <- d[locus == "Non-carriers"]
  r <- data.table(measure = y, group = G, n = nrow(x), mean = mean(x$yy), se = sd(x$yy) / sqrt(nrow(x)), age_med = median(x$cnb_age))
  if (G != "Non-carriers") {
    z <- rbind(x[, .(yy, cnb_age, sex, cnv = 1)], ctl[, .(yy, cnb_age, sex, cnv = 0)])
    co <- summary(lm(yy ~ cnv + cnb_age + I(cnb_age^2) + sex, data = z))$coefficients["cnv", ]
    r[, `:=`(diff = mean(x$yy) - mean(ctl$yy), p_welch = welch(x$yy, ctl$yy), beta_adj = co[1], se_adj = co[2], p_adj = co[4])]
  }
  r
}))
G <- rbind(grp_rows("accuracy"), grp_rows("speed"))

# ---------------------------------------------------------------- BMI contrasts within groups
bmi_rows <- function(y) rbindlist(fill = TRUE, lapply(GROUPS, function(Gn) {
  d <- cnb[!is.na(get(y)) & !is.na(bmiZ)]; d[, yy := get(y)]
  x <- if (Gn %in% c("Overgrowth CNVs", "Undergrowth CNVs")) d[pooled == Gn] else d[locus == Gn]
  lo <- x[tail == "Low BMI"]; hi <- x[tail == "High BMI"]
  r <- data.table(measure = y, group = Gn, n_bmi = nrow(x), n_lo = nrow(lo), n_hi = nrow(hi),
                  mean_lo = mean(lo$yy), mean_hi = mean(hi$yy), se_lo = sd(lo$yy) / sqrt(nrow(lo)), se_hi = sd(hi$yy) / sqrt(nrow(hi)),
                  diff_hi_lo = mean(hi$yy) - mean(lo$yy), p_welch = welch(lo$yy, hi$yy))
  if (nrow(x) >= 10) {
    xt <- x[!is.na(tail)][, hi := as.numeric(tail == "High BMI")]
    if (nrow(lo) >= 3 && nrow(hi) >= 3) { co <- summary(lm(yy ~ hi + cnb_age + I(cnb_age^2) + sex, data = xt))$coefficients["hi", ]; r[, `:=`(beta_adj = co[1], se_adj = co[2], p_adj = co[4])] }
    cs <- summary(lm(yy ~ bmiZ + cnb_age + I(cnb_age^2) + sex, data = x))$coefficients["bmiZ", ]; r[, `:=`(slope_bmiZ = cs[1], se_slope = cs[2], p_slope = cs[4])]
  }
  if (Gn != "Non-carriers") {   # CNV x BMI (continuous) vs non-carriers
    ctl <- d[locus == "Non-carriers"]
    if (nrow(ctl) >= 10) {
      z <- rbind(x[, .(yy, bmiZ, cnb_age, sex, cnv = 1)], ctl[, .(yy, bmiZ, cnb_age, sex, cnv = 0)])
      ci <- summary(lm(yy ~ cnv * bmiZ + cnb_age + I(cnb_age^2) + sex, data = z))$coefficients["cnv:bmiZ", ]
      r[, `:=`(cnv_x_bmi = ci[1], se_cnv_x_bmi = ci[2], p_cnv_x_bmi = ci[4])]
    }
  }
  r
}))
Bm <- rbind(bmi_rows("accuracy"), bmi_rows("speed"))
res <- rbind(G[, type := "CNV group vs non-carriers"], Bm[, type := "High vs low BMI within group"], fill = TRUE)
fwrite(res, file.path(HERE, "cnb_results.csv"))
cat("\n== group contrasts\n"); print(G[, .(measure, group, n, mean = round(mean, 2), diff = round(diff, 2), p_welch = signif(p_welch, 2), beta_adj = round(beta_adj, 2), p_adj = signif(p_adj, 2))])
cat("\n== BMI contrasts\n"); print(Bm[, .(measure, group, n_lo, n_hi, diff_hi_lo = round(diff_hi_lo, 2), p_welch = signif(p_welch, 2), beta_adj = round(beta_adj, 2), p_adj = signif(p_adj, 2), slope = round(slope_bmiZ, 2), p_slope = signif(p_slope, 2))])   # CNV x BMI vs controls not estimable (controls lack BMI)

# ---------------------------------------------------------------- per-test means by locus (supplement)
long <- melt(cnb[, c("GUID", "locus", ACC, SPD), with = FALSE], id.vars = c("GUID", "locus"), variable.name = "var", value.name = "z")[!is.na(z)]
long[, domain := fifelse(var %in% ACC, "Accuracy", "Speed")][, test := fifelse(domain == "Accuracy", names(ACC)[match(var, ACC)], names(SPD)[match(var, SPD)])]
long[domain == "Speed", z := -z]
bt <- long[, .(n = .N, mean_z = mean(z), se = sd(z) / sqrt(.N)), by = .(domain, test, locus)]
bt <- merge(bt, bt[locus == "Non-carriers", .(domain, test, ctl_mean = mean_z)], by = c("domain", "test"))[, diff_vs_ctl := mean_z - ctl_mean]
bt[, p_welch := mapply(function(dm, te, lo) { a <- long[domain == dm & test == te & locus == lo, z]; b <- long[domain == dm & test == te & locus == "Non-carriers", z]; welch(a, b) }, domain, test, locus)]
fwrite(bt, file.path(HERE, "cnb_by_test.csv"))

# ---------------------------------------------------------------- figure: CNB trajectories by BMI tail, carriers only
# JS 2026-10-23 (Jon): CNB replaces height in the trajectory view. One panel per locus and measure: composite z by CNB age,
# BMI tails as in the main figure (>= 70th centile red, <= 30th blue; mid-BMI not shown), GAM curve per tail with 95% band.
# No control comparison (G2MH non-carriers have almost no BMI). Printed: "BMI p" = age-adjusted continuous-BMI slope within
# the group (lm on bmiZ + age + age^2 + sex) and "age x BMI p" = within-group tensor interaction (gam: s(age) + s(BMI) +
# ti(age, BMI), REML). NB the z-scores are already normalised to age-matched non-carriers, so the curves describe how the
# deficit relative to norms changes with age. Separate figure, not yet part of the main figure.
COLS <- c("BMI ≥ 70th centile" = "#F0645A", "BMI ≤ 30th centile" = "#00A5AD")
within_p <- function(x) tryCatch(summary(gam(yy ~ s(cnb_age, k = 8) + s(bmiZ, k = 5) + ti(cnb_age, bmiZ, k = c(5, 5)), data = x, method = "REML"))$s.table["ti(cnb_age,bmiZ)", "p-value"], error = function(e) NA_real_)
traj_stats <- rbindlist(lapply(c("accuracy", "speed"), function(y) rbindlist(lapply(LAB[-1], function(L) {
  x <- cnb[locus == L & !is.na(get(y)) & !is.na(bmiZ)]; x[, yy := get(y)]
  cs <- summary(lm(yy ~ bmiZ + cnb_age + I(cnb_age^2) + sex, data = x))$coefficients["bmiZ", ]
  data.table(measure = y, locus = L, n = nrow(x), n_hi = sum(x$tail == "High BMI", na.rm = TRUE), n_lo = sum(x$tail == "Low BMI", na.rm = TRUE),
             slope_bmiZ = cs[1], se_slope = cs[2], p_bmi = cs[4], p_age_x_bmi = within_p(x))
}))))
fwrite(traj_stats, file.path(HERE, "cnb_trajectory_tests.csv")); cat("\n== trajectory tests\n"); print(traj_stats[, .(measure, locus, n, n_hi, n_lo, slope = round(slope_bmiZ, 3), p_bmi = signif(p_bmi, 2), p_age_x_bmi = signif(p_age_x_bmi, 2))])
traj <- function(y, L, ylab, title, xmax = 65) {
  x <- cnb[locus == L & !is.na(get(y)) & !is.na(tail)]; x[, yy := get(y)]
  x[, grp := factor(ifelse(tail == "High BMI", names(COLS)[1], names(COLS)[2]), levels = names(COLS))]
  st <- traj_stats[measure == y & locus == L]
  lab <- gsub("= < ", "< ", sprintf("BMI p = %s\nage × BMI p = %s", pfx(st$p_bmi), pfx(st$p_age_x_bmi)), fixed = TRUE)
  yl <- c(-4, 2.5)
  ggplot(x, aes(cnb_age, yy, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 1.4, alpha = 0.7) +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 5), method.args = list(method = "REML"), linewidth = 1, alpha = 0.2) +
    annotate("text", x = xmax, y = yl[1], label = lab, hjust = 1, vjust = 0, size = 3.3, fontface = ifelse(any(c(st$p_bmi, st$p_age_x_bmi) < 0.05, na.rm = TRUE), "bold", "plain")) +
    scale_colour_manual(values = COLS, breaks = names(COLS)) + scale_fill_manual(values = COLS, breaks = names(COLS)) +
    coord_cartesian(xlim = c(6, xmax), ylim = yl) +
    labs(x = "Age at CNB (years)", y = ylab, colour = NULL, fill = NULL, title = title, subtitle = sprintf("%s (n = %d)", L, nrow(x))) +
    theme_classic(base_size = 14) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.position.inside = c(0.75, 0.93),
                                          legend.background = element_blank(), legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
}
LOCI4 <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
pl <- c(lapply(seq_along(LOCI4), function(i) traj("accuracy", LOCI4[i], "CNB accuracy (z vs non-carriers)", if (i == 1) "Accuracy" else NULL)),
        lapply(seq_along(LOCI4), function(i) traj("speed", LOCI4[i], "CNB speed (z, higher = faster)", if (i == 1) "Speed" else NULL)))
fig <- wrap_plots(pl, ncol = 2, byrow = FALSE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20, face = "bold"))
ggsave(file.path(HERE, "fig_cnb_trajectories.png"), fig, width = 12, height = 16, dpi = 170)
tryCatch(ggsave(file.path(HERE, "fig_cnb_trajectories.pdf"), fig, width = 12, height = 16, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_cnb_trajectories.pdf"), fig, width = 12, height = 16))
cat("done\n")

# ---------------------------------------------------------------- JS 2026-10-23: CNB by head circumference (F/L-style bars)
# HC Z (UK90, adults at the 18-y reference; 00_scoring.R) from the growth genomatrix visit nearest the CNB age (<= 3 y), else
# the release medical form (med_headcircum at dem_age). Tails: HC Z <= 30th vs >= 70th centile. Bars = mean CNB composite
# +/- SE per tail within each CNV group (loci and the overgrowth/undergrowth pools); Delta and Welch p = large - small HC,
# adj. = lm with CNB age, age^2, sex; also the continuous HC Z slope. Non-carriers shown only if both tails have >= 10.
gh <- gm[cohort == "G2MH" & !(Cardiff_mask %in% TRUE) & !is.na(GUID) & GUID != "" & !is.na(HEADCIRCUM_raw), .(GUID, age, HEADCIRCUM_raw)]
vh <- merge(cnb[, .(GUID, cnb_age)], gh, by = "GUID", allow.cartesian = TRUE)[, dage := abs(age - cnb_age)][order(GUID, dage)][, .SD[1], by = GUID][dage <= 3]
cnb <- merge(cnb, vh[, .(GUID, hc_age = age, hc_raw = HEADCIRCUM_raw, hc_src = "genomatrix visit")], by = "GUID", all.x = TRUE)
relh <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"), select = c("rarecnv_id", "dem_age", "med_headcircum"))
relh <- relh[, .(rarecnv_id, rel_hage = as.numeric(dem_age), rel_hc = suppressWarnings(as.numeric(med_headcircum)))][!is.na(rel_hc) & rel_hc >= 40 & rel_hc <= 70 & !is.na(rel_hage)]
cnb <- merge(cnb, relh, by = "rarecnv_id", all.x = TRUE)
cnb[is.na(hc_raw) & !is.na(rel_hc) & abs(rel_hage - cnb_age) <= 3, `:=`(hc_age = rel_hage, hc_raw = rel_hc, hc_src = "release medical form")]
cnb[, hcZ := score_hc(hc_raw, hc_age, sex)]
cnb[, hc_tail := fifelse(hcZ >= qnorm(.7), "Large HC", fifelse(hcZ <= qnorm(.3), "Small HC", NA_character_))]
cat("HC source:\n"); print(cnb[, .N, by = .(locus, hc_src)][order(locus, hc_src)])
cat("HC tails by locus (with accuracy):\n"); print(dcast(cnb[!is.na(hc_tail) & !is.na(accuracy), .N, by = .(locus, hc_tail)], locus ~ hc_tail, value.var = "N", fill = 0))
hc_rows <- function(y) rbindlist(fill = TRUE, lapply(GROUPS, function(Gn) {
  d <- cnb[!is.na(get(y)) & !is.na(hcZ)]; d[, yy := get(y)]
  x <- if (Gn %in% c("Overgrowth CNVs", "Undergrowth CNVs")) d[pooled == Gn] else d[locus == Gn]
  lo <- x[hc_tail == "Small HC"]; hi <- x[hc_tail == "Large HC"]
  r <- data.table(measure = y, group = Gn, n_hc = nrow(x), n_small = nrow(lo), n_large = nrow(hi),
                  mean_small = mean(lo$yy), mean_large = mean(hi$yy), se_small = sd(lo$yy) / sqrt(nrow(lo)), se_large = sd(hi$yy) / sqrt(nrow(hi)),
                  diff_large_small = mean(hi$yy) - mean(lo$yy), p_welch = welch(lo$yy, hi$yy))
  if (nrow(lo) >= 3 && nrow(hi) >= 3) { xt <- x[!is.na(hc_tail)][, big := as.numeric(hc_tail == "Large HC")]
    co <- summary(lm(yy ~ big + cnb_age + I(cnb_age^2) + sex, data = xt))$coefficients["big", ]; r[, `:=`(beta_adj = co[1], se_adj = co[2], p_adj = co[4])] }
  if (nrow(x) >= 10) { cs <- summary(lm(yy ~ hcZ + cnb_age + I(cnb_age^2) + sex, data = x))$coefficients["hcZ", ]; r[, `:=`(slope_hcZ = cs[1], se_slope = cs[2], p_slope = cs[4])] }
  r
}))
Hc <- rbind(hc_rows("accuracy"), hc_rows("speed"))
fwrite(Hc, file.path(HERE, "cnb_hc_results.csv"))
cat("\n== HC contrasts\n"); print(Hc[, .(measure, group, n_small, n_large, diff = round(diff_large_small, 2), p_welch = signif(p_welch, 2), beta_adj = round(beta_adj, 2), p_adj = signif(p_adj, 2), slope = round(slope_hcZ, 2), p_slope = signif(p_slope, 2))])
HCCOL <- c("Small HC" = "#8C9A9A", "Large HC" = "#4B6EA9")
ORD2 <- c("Non-carriers", "Overgrowth CNVs", "Undergrowth CNVs", "16p11.2 del", "22q11.2 dup", "16p11.2 dup", "22q11.2 del")
hc_bar <- function(meas, ylab) {
  d <- Hc[measure == meas & n_small >= 10 & n_large >= 10][, group := factor(group, ORD2)]
  d[, group := factor(as.character(group), ORD2[ORD2 %in% as.character(group)])]
  long2 <- rbind(d[, .(group, tail = "Small HC", m = mean_small, se = se_small, n = n_small)], d[, .(group, tail = "Large HC", m = mean_large, se = se_large, n = n_large)])
  long2[, tail := factor(tail, names(HCCOL))]
  d[, lab := sprintf("Δ %+.2f, p = %s\nadj. %+.2f\np = %s", diff_large_small, pfx(p_welch), beta_adj, pfx(p_adj))][, bold := !is.na(p_adj) & p_adj < 0.05]
  yr <- range(c(long2$m - long2$se, long2$m + long2$se, 0)); top <- yr[2] + 0.6 * diff(yr); bot <- yr[1] - 0.14 * diff(yr)
  ggplot(long2, aes(group, m, fill = tail)) + geom_hline(yintercept = 0, colour = "grey50") + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(y = bot, label = n), position = position_dodge(0.8), size = 3.2) +
    geom_text(data = d, aes(group, top, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 3.3, vjust = 1, lineheight = 0.9) +
    scale_fill_manual(values = HCCOL, name = NULL, labels = c("Head circumference ≤ 30th centile", "Head circumference ≥ 70th centile")) +
    coord_cartesian(ylim = c(bot, top)) + labs(x = NULL, y = ylab, subtitle = "Small vs large head circumference") +
    theme_classic(base_size = 15) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
                                           legend.text = element_text(size = 9), legend.key.size = unit(10, "pt"), axis.text.x = element_text(size = 10.5))
}
fig_hc <- (hc_bar("accuracy", "CNB accuracy (z vs non-carriers)") / hc_bar("speed", "CNB speed (z, higher = faster)")) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))
ggsave(file.path(HERE, "fig_cnb_by_hc.png"), fig_hc, width = 11, height = 10, dpi = 170)
tryCatch(ggsave(file.path(HERE, "fig_cnb_by_hc.pdf"), fig_hc, width = 11, height = 10, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_cnb_by_hc.pdf"), fig_hc, width = 11, height = 10))
cat("HC figure done\n")
