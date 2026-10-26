# Head-circumference trajectories 0-18 y, keeping every SVIP visit (the clinical genomatrix keeps one SVIP HC record
# per child; SVIP v10 growth_measurements.csv has the full longitudinal series, mostly birth-5 y).
#
# Data: SVIP v10 (all carrier visits with HC, replacing the genomatrix SVIP rows) + Cardiff, G2MH, SSC from the
#   genomatrix (analysis_table.rds; Cardiff duplicates in G2MH already removed).
# One scale for everyone, recomputed from raw values with childsds (as in Molly's CreateGenomatrix2.R):
#   HC Z  = UK90 (uk1990.ref) from birth to 18 y
#   BMI Z = WHO (who.ref) < 5 y, WHO 2007 (who2007.ref) 5-18 y
# Cleaning (per child, age order, vs the last kept visit): duplicate entry (< 0.04 y apart under 2 y, < 0.25 y after);
#   identical HC >= 0.1 y later before age 3 (carried forward); HC smaller by > 1.5 cm or growing > 15 cm/y after 1 y;
#   HC Z jump > 2.5 SD.  Visits before 0.1 y (birth records) are excluded.
# Analyses:
#   A. BMI-stratified trajectories (split rule: 16p11.2 del BMI > 90th vs rest; other loci <= 30th vs >= 70th),
#      stratum = child's mean BMI Z over visits from age 2 (all visits if none), mixed GAM with child random intercept.
#   B. Age-varying prospective effect: change in HC Z to the next visit per 1 SD of BMI Z at a visit, by age at that visit.
#   C. Fixed windows (< 2, 2-5, 5-18 y) for B.
# Output: hc_visits.rds, hc_cleaning_log.csv, hc_stratified_tests.csv, hc_prospective_windows.csv,
#         fig_hc_trajectories.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(lmerTest); library(splines); library(ggplot2); library(patchwork) })
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressPackageStartupMessages(library(childsds))
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")   # participant-level tables stay in Jon/Claude (controlled data)
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
SPLIT <- c("16p11.2 del" = "30_70", "16p11.2 dup" = "30_70", "22q11.2 del" = "30_70", "22q11.2 dup" = "30_70")   # JS 2026-10-12: 30/70 tails for every locus (was > 90th vs rest for 16p11.2 del)
COL <- c("#00A5AD", "#F0645A")

# ---------------------------------------------------------------- data
gm <- readRDS(file.path(DATA, "analysis_table.rds"))[cohort %in% c("Cardiff", "G2MH", "SSC") & locus %in% LOCI & age < 18 & !is.na(HEADCIRCUM_raw)]
gm <- gm[, .(id, cohort, locus = as.character(locus), sex = as.character(sex), age, hc = HEADCIRCUM_raw, bmi = BMI_raw)]
SV <- file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0")   # JS 2026-10-24: v10 folder moved here (was Autism/searchlight_multimodal/...)
summ <- fread(file.path(SV, "svip_summary_variables.csv"))
summ <- summ[, .(individual, status = svip_summary_variables.genetic_status_16p, sex = svip_summary_variables.sex)]
sv <- fread(file.path(SV, "growth_measurements.csv"))
sv <- merge(sv[, .(individual, age = eval_age_months / 12, hc = growth_measurements_raw.head_circumference, bmi = growth_measurements_raw.bmi)],
            summ[status %in% c("deletion", "duplication")], by = "individual")
sv <- sv[!is.na(hc) & age < 18, .(id = paste0("SVIP:", individual), cohort = "SVIP",
                                  locus = ifelse(status == "deletion", "16p11.2 del", "16p11.2 dup"), sex, age, hc, bmi)]
d <- rbind(gm, sv)[age >= 0.1 & sex %in% c("male", "female")]
# JS 2026-10-24: plausibility bounds as in 00_qc_growth.R (genomatrix rows are already QC'd in 01_build; SVIP v10 rows screened here)
d[!is.na(hc) & (hc < 30 | hc > 70), hc := NA_real_]; d[!is.na(bmi) & (bmi < 10 | bmi > 70), bmi := NA_real_]; d <- d[!is.na(hc)]
d[, hcz := sds(hc, age = age, sex = sex, male = "male", female = "female", item = "hc", ref = uk1990.ref)]
d[, bmiz := NA_real_]
d[age < 5 & !is.na(bmi), bmiz := sds(bmi, age = age, sex = sex, male = "male", female = "female", item = "bmi", ref = who.ref)]
d[age >= 5 & !is.na(bmi), bmiz := sds(bmi, age = age, sex = sex, male = "male", female = "female", item = "bmi", ref = who2007.ref)]
d <- d[!is.na(hcz) & abs(hcz) <= 5][order(id, age)]; d[!is.na(bmiz) & abs(bmiz) > 5, bmiz := NA_real_]   # JS 2026-10-24: |Z| <= 5 (was < 10), as in 00_qc_growth.R

# ---------------------------------------------------------------- cleaning
keep <- logical(nrow(d)); why <- character(nrow(d))
for (ix in split(seq_len(nrow(d)), d$id)) {
  last <- ix[1]; keep[last] <- TRUE
  for (i in ix[-1]) {
    dt <- d$age[i] - d$age[last]; dh <- d$hc[i] - d$hc[last]; dz <- d$hcz[i] - d$hcz[last]
    w <- if (dt < ifelse(d$age[i] < 2, 0.04, 0.25)) "duplicate_visit"
         else if (abs(dh) < 0.05 && dt >= 0.1 && d$age[i] < 3) "carried_forward"
         else if (dh < -1.5 || (d$age[i] > 1 && dh / dt > 15)) "impossible_change"
         else if (abs(dz) > 2.5) "jump_gt_2.5SD" else ""
    if (w == "") { keep[i] <- TRUE; last <- i } else why[i] <- w
  }
}
cl <- data.table(cohort = d$cohort, why = why)[why != "", .N, by = .(cohort, why)]
fwrite(cl, file.path(HERE, "hc_cleaning_log.csv")); cat("HC visits dropped by cleaning:\n"); print(cl)
d <- d[keep]
saveRDS(d, file.path(DATA, "hc_visits.rds"))
d[, n_vis := .N, by = id]
cat("\nchildren / visits / children with repeat HC, by locus and cohort:\n")
print(d[, .(children = uniqueN(id), visits = .N, repeat_children = uniqueN(id[n_vis > 1])), by = .(locus, cohort)][order(locus, cohort)])

# ---------------------------------------------------------------- A. stratified trajectories
strat_of <- function(x, L) {
  if (SPLIT[[L]] == "p90") { lv <- c("BMI \u2264 90th centile", "BMI > 90th centile"); factor(ifelse(x > qnorm(.9), lv[2], lv[1]), levels = lv) }
  else { lv <- c("BMI \u2264 30th centile", "BMI \u2265 70th centile")
         factor(ifelse(x >= qnorm(.7), lv[2], ifelse(x <= qnorm(.3), lv[1], NA)), levels = lv) }
}
resA <- list(); plotsA <- list()
for (L in LOCI) {
  x <- d[locus == L]
  mb <- x[!is.na(bmiz), .(b2 = mean(bmiz[age >= 2]), ball = mean(bmiz)), by = id][, bm := fifelse(is.nan(b2), ball, b2)]
  x <- merge(x, mb[, .(id, bm)], by = "id")[, strat := strat_of(bm, L)][!is.na(strat)]
  x[, `:=`(fid = factor(id), cohort = factor(cohort), strat_ord = as.ordered(strat))]
  ct <- if (nlevels(x$cohort) > 1) "+ cohort" else ""
  ka <- min(12, x[, uniqueN(round(age, 2)), by = strat][, min(V1)] - 1)   # Cardiff/G2MH ages are whole years
  m <- gam(as.formula(sprintf("hcz ~ strat_ord + sex %s + s(age, k = %d) + s(age, by = strat_ord, k = %d) + s(fid, bs = 're')", ct, ka, ka)), data = x, method = "REML")
  s <- summary(m)
  r <- data.table(locus = L, split = SPLIT[[L]], children = uniqueN(x$id), visits = nrow(x), repeat_children = x[, .N, by = id][N > 1, .N],
                  children_high = x[as.integer(strat) == 2, uniqueN(id)], offset = s$p.coeff[2], p_offset = s$p.pv[2],
                  p_curve_shape = s$s.table[grep("strat_ord", rownames(s$s.table)), "p-value"])
  resA[[L]] <- r
  nd <- CJ(age = seq(0.1, 18, 0.1), strat = levels(x$strat), sex = c("male", "female"), cohort = levels(x$cohort))
  rng <- x[, .(lo = min(age), hi = max(age)), by = .(strat = as.character(strat))]
  nd <- merge(nd, rng, by = "strat")[age >= lo & age <= hi]
  nd[, `:=`(strat_ord = factor(strat, levels = levels(x$strat), ordered = TRUE), fid = x$fid[1], cohort = factor(cohort, levels(x$cohort)))]
  pr <- predict(m, nd, exclude = "s(fid)", se.fit = TRUE); nd[, `:=`(fit = pr$fit, se = pr$se.fit)]
  cur <- nd[, .(fit = mean(fit), se = mean(se)), by = .(strat, age)][, strat := factor(strat, levels(x$strat))]
  lab <- sprintf("high \u2212 low offset %+.2f, p = %s\ncurve-shape difference p = %s\n%d children (%d with repeat visits)",
                 r$offset, format.pval(r$p_offset, 2, eps = 1e-4), format.pval(r$p_curve_shape, 2, eps = 1e-4), r$children, r$repeat_children)
  plotsA[[L]] <- ggplot() +
    annotate("rect", xmin = 8, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey93") + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(data = x, aes(age, hcz, colour = strat), size = 0.7, alpha = 0.3, shape = 16) +
    geom_ribbon(data = cur, aes(age, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se, fill = strat), alpha = 0.25) +
    geom_line(data = cur, aes(age, fit, colour = strat), linewidth = 1) +
    annotate("text", x = 18, y = -4.5, label = lab, hjust = 1, vjust = 0, size = 2.5) +
    scale_colour_manual(values = COL) + scale_fill_manual(values = COL) + coord_cartesian(xlim = c(0, 18), ylim = c(-4.5, 4.5)) +
    labs(x = "Age (years)", y = "Head circumference Z (UK90)", colour = NULL, fill = NULL, subtitle = L) +
    theme_classic(base_size = 10) + theme(legend.position = "inside", legend.position.inside = c(0.3, 0.93), legend.background = element_blank(),
                                           legend.key.height = unit(8, "pt"))
}
resA <- rbindlist(resA); fwrite(resA, file.path(HERE, "hc_stratified_tests.csv")); print(resA, digits = 3)

# ---------------------------------------------------------------- B/C. prospective, age-varying and windows
iv <- d[!is.na(bmiz)][order(id, age)][, .(age0 = head(age, -1), dt = diff(age), Z0 = head(hcz, -1), dZ = diff(hcz), B0 = head(bmiz, -1),
                                          sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id][dt >= 0.1]
iv[, `:=`(fid = factor(id), cohort = factor(cohort), sex = factor(sex))]
curves <- list(); win <- list(); GRID <- seq(0.1, 17, 0.1)
for (L in c("Pooled", LOCI)) {
  ii <- if (L == "Pooled") iv else iv[locus == L]; ii[, cohort := droplevels(cohort)]
  if (uniqueN(ii$id) < 20) next
  adj <- paste0(if (L == "Pooled") "+ locus " else "", if (nlevels(ii$cohort) > 1) "+ cohort" else "")
  kk <- min(if (uniqueN(ii$id) < 60) 4 else 6, uniqueN(round(ii$age0, 2)) - 1)
  m <- gam(as.formula(sprintf("dZ ~ Z0 + dt + sex %s + s(age0, k = %d) + s(age0, by = B0, k = %d) + s(fid, bs = 're')", adj, kk, kk)), data = ii, method = "REML")
  nd1 <- data.table(age0 = GRID, Z0 = 0, dt = 1, B0 = 1, sex = factor("male", levels(ii$sex)), locus = ii$locus[1], cohort = ii$cohort[1], fid = ii$fid[1])
  nd0 <- copy(nd1)[, B0 := 0]
  X <- predict(m, nd1, type = "lpmatrix", exclude = "s(fid)") - predict(m, nd0, type = "lpmatrix", exclude = "s(fid)")
  X[, -grep("s(age0):B0", colnames(X), fixed = TRUE)] <- 0
  cb <- data.table(age = GRID, beta = drop(X %*% coef(m)), se = sqrt(rowSums((X %*% vcov(m)) * X)))
  cb <- cb[age >= quantile(ii$age0, .02) & age <= quantile(ii$age0, .98)]
  curves[[L]] <- cb[, `:=`(locus = L, p = summary(m)$s.table[grep("B0", rownames(summary(m)$s.table)), "p-value"],
                           children = uniqueN(ii$id), intervals = nrow(ii))]
  ii[, w := cut(age0, c(0, 2, 5, 18), right = FALSE, labels = c("<2 y", "2-5 y", "5-18 y"))]
  adjl <- paste0(if (L == "Pooled") "+ locus " else "", if (nlevels(ii$cohort) > 1) "+ cohort" else "")
  ml <- tryCatch(suppressMessages(lmer(as.formula(sprintf("dZ ~ B0:w + w + Z0 + dt + sex * ns(age0, 3) %s + (1 | fid)", adjl)), data = ii)), error = function(e) NULL)
  if (!is.null(ml)) {
    co <- summary(ml)$coefficients; r <- grep(":B0|B0:", rownames(co))
    win[[L]] <- merge(data.table(locus = L, w = sub("w", "", sub(":B0|B0:", "", rownames(co)[r])), est = co[r, 1], se = co[r, 2], p = co[r, 5]),
                      ii[, .(children = uniqueN(id), intervals = .N), by = .(w = as.character(w))], by = "w")
  }
}
curves <- rbindlist(curves); win <- rbindlist(win)
fwrite(curves, file.path(HERE, "hc_prospective_age_varying.csv")); fwrite(win, file.path(HERE, "hc_prospective_windows.csv"))
print(unique(curves[, .(locus, p, children, intervals)])); print(win[order(locus, w)], digits = 3)

curves[, locus := factor(locus, levels = c("Pooled", LOCI))]
labB <- unique(curves[, .(locus, txt = sprintf("age-varying BMI effect p = %s\n%d children, %d intervals", format.pval(p, 2, eps = 1e-4), children, intervals))])
pB <- ggplot(curves, aes(age, beta)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), fill = "#F0645A", alpha = 0.25) + geom_line(colour = "#C0392B", linewidth = 0.9) +
  geom_text(data = labB, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.3, size = 2.6) +
  facet_wrap(~locus, nrow = 1, scales = "free_x") +
  labs(x = "Age at visit (years)", y = "Change in HC Z to next visit\nper 1 SD of BMI Z (\u00b1 95% CI)",
       title = "B. Prospective: does high BMI at a visit precede faster head growth?") + theme_bw(base_size = 10)
fig <- (wrap_plots(plotsA, nrow = 1) + plot_annotation()) / pB +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(title = "Head-circumference trajectories by BMI in 16p11.2 and 22q11.2 carriers, 0-18 y (all visits, incl. full SVIP longitudinal series)",
                  subtitle = paste("A (top). Mixed GAM \u00b1 95% CI, child random intercept; 16p11.2 del BMI > 90th vs rest, other loci \u2265 70th vs \u2264 30th centile (child's mean BMI Z from age 2).",
                                   "\nUK90 head-circumference Z and WHO BMI Z recomputed from raw values for all cohorts (SVIP v10, Cardiff, G2MH, SSC). Grey = 8-13 y."))
ggsave(file.path(HERE, "fig_hc_trajectories.png"), fig, width = 16, height = 9, dpi = 200)
ggsave(file.path(HERE, "fig_hc_trajectories.pdf"), fig, width = 16, height = 9, device = cairo_pdf)
