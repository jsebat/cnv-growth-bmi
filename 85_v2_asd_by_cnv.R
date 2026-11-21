# 85_v2_asd_by_cnv.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 2)
# What it does: autism diagnosis within CNV group, revised per Jon ("replace undergrowth/overgrowth with pooled (all CNVs)";
#   "does stratifying the pooled analysis by BMI show an increase in autism?").
#   a-c  height Z, head-circumference Z and BMI Z in carriers with vs without an autism diagnosis: each CNV, all four CNVs pooled
#        (adjusted for cohort within CNV), SSC probands vs designated siblings. Same data and model as 32_supp_asd_by_cnv.R
#        (ADD_CHOP=TRUE, VISIT_RULE=fill_gaps; gam y ~ ASD + s(age, k = 8) + sex [+ cohort] [+ CNV], REML). Heterogeneity of the
#        ASD difference across the four CNVs: Cochran Q on the four per-CNV estimates (and the ASD x CNV F test of the pooled model).
#        Overgrowth / undergrowth groupings are kept in the table only.
#   d    NEW: proportion with an autism diagnosis by BMI band (<= 30th, 30th-70th, >= 70th centile of the BMI Z reference), pooled
#        and per CNV; test = logistic regression diagnosis ~ BMI Z + ns(age, 3) + sex + cohort within CNV; per cohort, leave one
#        cohort out, without cohort terms, children only, linear age, and with Searchlight diagnoses added.
#   e    NEW: the same for head-circumference Z bands.
# Inputs : growth table of 23_unaffected_controls.R / 32 (G2MH, Cardiff, Simons VIP, SSC, CHU Sainte-Justine, CHOP), via v2_base().
#   Autism diagnosis: G2MH / Cardiff subjects.csv; Simons VIP diagnosis_summary; SSC proband = diagnosis; CHUSJ staging sheet;
#   CHOP index sheet ("ASD diagnosis" Yes / No; the sheet does not say whether "No" means evaluated).
# Outputs: run_v2_2026-11-20/autism/fig_asd_by_cnv_v2.{png,pdf}, asd_by_cnv_tests_v2.csv, asd_by_bmi_stratum.csv,
#          asd_by_bmi_logistic.csv, asd_by_cnv_n_by_cohort_sex.csv.
# Cohort adjustment (JS 2026-11-21): cohort and CNV are partly confounded (CHOP has only 22q11.2 carriers), so a model that pools
#   CNVs carries one intercept per cohort-within-CNV cell (cg = cohort x CNV) instead of additive cohort + CNV terms. In a model of
#   one CNV this is the cohort term as before. ADDITIVE = TRUE restores the earlier pooled adjustment (kept as sensitivity rows).
# null model: none; Wald tests (gam coefficient of ASD; logistic coefficient of BMI Z or head-circumference Z), Cochran Q for
#   heterogeneity across the four CNVs and across cohorts.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
suppressPackageStartupMessages(library(splines))
v2_base()
NCG <- "Non-carriers (SSC)"; PH <- c("Height", "Head circumference", "BMI"); POOL <- "All four CNVs"
X <- XA[pheno %in% PH & locus %in% c(LOC4, NCG)]
GROUPS <- c(setNames(as.list(LOC4), LOC4), setNames(list(LOC4), POOL), list(`Overgrowth CNVs (16p11.2 del, 22q11.2 dup)` = c("16p11.2 del", "22q11.2 dup"), `Undergrowth CNVs (16p11.2 dup, 22q11.2 del)` = c("16p11.2 dup", "22q11.2 del")), setNames(list(NCG), NCG))
# Searchlight / Simons VIP diagnoses missing from 32 (sensitivity only): Searchlight v14 clinical or reported ASD (as in 63)
XS <- merge(A[pheno %in% c("Height", "Head circumference") & locus %in% LOC4], idmap, by = c("id", "cohort"), all.x = TRUE)
XS <- merge(XS, sl, by.x = "sampleID", by.y = "id", all.x = TRUE)[, asd := fifelse(!is.na(asd), as.integer(asd), a_sl)][!is.na(asd)]
XS <- rbind(XS[, .(id, cohort, locus, sex, age, BMI, y, pheno, asd)], XS[pheno == "Height", .(id, cohort, locus, sex, age, BMI, y = BMI, pheno = "BMI", asd)])
# ---------------------------------------------------------------- a-c: outcome by diagnosis (32 model)
ADDITIVE <- FALSE
additive <- function(expr) { ADDITIVE <<- TRUE; on.exit(ADDITIVE <<- FALSE); force(expr) }
fit <- function(x) tryCatch({
  x <- copy(x)[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))][, cg := factor(paste(cohort, locus))]
  cov <- if (ADDITIVE) paste(c(if (nlevels(x$cohort) > 1) "+ cohort", if (nlevels(x$locus) > 1) "+ locus"), collapse = " ") else if (nlevels(x$cg) > 1) "+ cg" else ""
  m <- gam(as.formula(paste("y ~ asd + s(age, k = 8) + sex", cov)), data = x, method = "REML"); co <- summary(m)$p.table["asd", ]
  data.table(delta = unname(co[1]), se = unname(co[2]), p = unname(co[4])) }, error = function(e) data.table(delta = NA_real_, se = NA_real_, p = NA_real_))
desc <- function(x) data.table(n_noasd = sum(x$asd == 0L), n_asd = sum(x$asd == 1L), mean_noasd = x[asd == 0L, mean(y)], se_noasd = x[asd == 0L, sd(y) / sqrt(.N)], mean_asd = x[asd == 1L, mean(y)], se_asd = x[asd == 1L, sd(y) / sqrt(.N)],
                               age_med_noasd = x[asd == 0L, median(age)], age_med_asd = x[asd == 1L, median(age)], cohorts = paste(sort(unique(x$cohort)), collapse = "+"))
ok <- function(x) sum(x$asd == 1L) >= 5 && sum(x$asd == 0L) >= 5
res <- rbindlist(lapply(PH, function(ph) rbindlist(lapply(names(GROUPS), function(g) { x <- X[pheno == ph & locus %in% GROUPS[[g]]]; if (!ok(x)) return(NULL)
  cbind(data.table(analysis = "main (as 32)", pheno = ph, group = g), desc(x), fit(x)) }))))
het <- rbindlist(lapply(PH, function(ph) { r <- res[pheno == ph & group %in% LOC4]; q <- cochranQ(r$delta, r$se); x <- copy(X[pheno == ph & locus %in% LOC4])[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))][, cg := factor(paste(cohort, locus))]
  m0 <- gam(y ~ asd + s(age, k = 8) + sex + cg, data = x, method = "ML"); m1 <- gam(y ~ asd:locus + s(age, k = 8) + sex + cg, data = x, method = "ML")
  data.table(analysis = "heterogeneity across the four CNVs", pheno = ph, group = POOL, n_noasd = sum(x$asd == 0L), n_asd = sum(x$asd == 1L), het_Q = q[["Q"]], het_df = q[["df"]], het_p = q[["p"]], het_I2 = q[["I2"]], asd_x_cnv_F_p = anova(m0, m1, test = "F")$`Pr(>F)`[2]) }))
# pooled: per cohort, leave one cohort out, with Searchlight diagnoses, children only
sens <- rbindlist(lapply(PH, function(ph) { x <- X[pheno == ph & locus %in% LOC4]; o <- list()
  for (k in sort(unique(x$cohort))) { d <- x[cohort == k]; if (ok(d)) o[[length(o) + 1]] <- cbind(data.table(analysis = paste("pooled, cohort:", k), pheno = ph, group = POOL), desc(d), fit(d))
    d <- x[cohort != k]; o[[length(o) + 1]] <- cbind(data.table(analysis = paste("pooled, leave out", k), pheno = ph, group = POOL), desc(d), fit(d)) }
  d <- x[age < 18]; o[[length(o) + 1]] <- cbind(data.table(analysis = "pooled, age < 18 y", pheno = ph, group = POOL), desc(d), fit(d))
  o[[length(o) + 1]] <- cbind(data.table(analysis = "pooled, additive cohort + CNV terms (main until 2026-11-21)", pheno = ph, group = POOL), desc(x), additive(fit(x)))
  for (g in c(LOC4[1:2], POOL)) { d <- XS[pheno == ph & locus %in% GROUPS[[g]]]; if (ok(d)) o[[length(o) + 1]] <- cbind(data.table(analysis = "with Searchlight v14 diagnoses added (Searchlight has no head circumference)", pheno = ph, group = g), desc(d), fit(d)) }
  rbindlist(o) }))
cq <- rbindlist(lapply(PH, function(ph) { r <- sens[pheno == ph & grepl("^pooled, cohort:", analysis)]; q <- cochranQ(r$delta, r$se)
  data.table(analysis = "heterogeneity across cohorts (pooled CNVs)", pheno = ph, group = POOL, het_Q = q[["Q"]], het_df = q[["df"]], het_p = q[["p"]], het_I2 = q[["I2"]]) }))
TAB <- rbind(res, het, sens, cq, fill = TRUE)
TAB[, note := "delta = adjusted difference (SD) autism diagnosis minus none; gam y ~ ASD + s(age) + sex + cohort within CNV (one intercept per cohort x CNV cell); cohorts G2MH, Cardiff, Simons VIP, SSC, CHUSJ, CHOP"]
fwrite(TAB, file.path(V2OUT, "asd_by_cnv_tests_v2.csv"))
options(width = 230); print(TAB[, .(analysis = substr(analysis, 1, 34), pheno = substr(pheno, 1, 6), group = substr(group, 1, 18), n_noasd, n_asd, delta = round(delta, 3), se = round(se, 3), p = signif(p, 3), het_p = signif(het_p, 3), Fp = signif(asd_x_cnv_F_p, 3))], nrows = 300)
ncs <- dcast(X[locus %in% LOC4, .N, by = .(outcome = pheno, locus, cohort, sex, asd = fifelse(asd == 1L, "autism", "no_autism"))], outcome + locus + cohort ~ sex + asd, value.var = "N", fill = 0)
fwrite(ncs, file.path(V2OUT, "asd_by_cnv_n_by_cohort_sex.csv"))
# ---------------------------------------------------------------- d, e: diagnosis by BMI band / head-circumference band
wil <- function(k, n) { z <- 1.96; p <- k / n; den <- 1 + z^2 / n; ce <- (p + z^2 / (2 * n)) / den; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den; list(lo = ce - h, hi = ce + h) }
BANDS <- c("\u2264 30th", "30th-70th", "\u2265 70th"); bandz <- function(z) factor(fifelse(z >= qnorm(.7), BANDS[3], fifelse(z <= qnorm(.3), BANDS[1], BANDS[2])), BANDS)
DB <- list(BMI = X[pheno == "BMI" & locus %in% LOC4][, z := y], `Head circumference` = X[pheno == "Head circumference" & locus %in% LOC4][, z := y])
DBS <- list(BMI = XS[pheno == "BMI"][, z := y], `Head circumference` = XS[pheno == "Head circumference"][, z := y])
strat <- rbindlist(lapply(names(DB), function(pr) { d <- copy(DB[[pr]])[, band := bandz(z)]
  s <- function(x, g, k) x[, .(n = .N, n_asd = sum(asd), mean_z = mean(z)), by = band][, `:=`(predictor = pr, group = g, cohort = k)]
  o <- rbind(s(d, POOL, "all cohorts"), rbindlist(lapply(LOC4, function(L) s(d[locus == L], L, "all cohorts"))), rbindlist(lapply(sort(unique(d$cohort)), function(k) s(d[cohort == k], POOL, k))),
             rbindlist(lapply(LOC4, function(L) rbindlist(lapply(sort(unique(d[locus == L, cohort])), function(k) s(d[locus == L & cohort == k], L, k))))))
  o[, prop := n_asd / n][, c("ci_lo", "ci_hi") := wil(n_asd, n)][order(group, cohort, band)] }))
setcolorder(strat, c("predictor", "group", "cohort", "band", "n", "n_asd", "prop", "ci_lo", "ci_hi", "mean_z"))
strat[, note := "band = centile band of the Z-score reference (Z <= -0.524, between, >= 0.524); prop = share with an autism diagnosis; Wilson 95% CI"]
fwrite(strat, file.path(V2OUT, "asd_by_bmi_stratum.csv"))
lg <- function(d, age = "ns", fe = TRUE) tryCatch({
  d <- copy(d)[, cg := paste(cohort, locus)]; kv <- if (ADDITIVE) "cohort" else "cg"; if (fe) { keep <- d[, .(v = sum(asd) > 0 & sum(asd) < .N), by = kv][v == TRUE][[kv]]; d <- d[get(kv) %in% keep] }   # a cohort without both outcomes carries no information once it has its own intercept
  if (sum(d$asd) < 8 || sum(1 - d$asd) < 8) stop("too few")
  tm <- c("z", if (age == "ns") "ns(age, 3)" else "age", if (uniqueN(d$sex) > 1) "sex", if (fe && ADDITIVE && uniqueN(d$cohort) > 1) "cohort", if ((ADDITIVE || !fe) && uniqueN(d$locus) > 1) "locus", if (fe && !ADDITIVE && uniqueN(d$cg) > 1) "cg")
  m <- suppressWarnings(glm(as.formula(paste("asd ~", paste(tm, collapse = " + "))), data = d, family = binomial)); co <- coef(summary(m))["z", ]
  data.table(n = nrow(d), n_asd = sum(d$asd), cohorts_used = paste(sort(unique(d$cohort)), collapse = "+"), logOR = co[[1]], se = co[[2]], OR = exp(co[[1]]), OR_lo = exp(co[[1]] - 1.96 * co[[2]]), OR_hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) },
  error = function(e) data.table(n = nrow(d), n_asd = sum(d$asd), cohorts_used = NA_character_, logOR = NA_real_, se = NA_real_, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
LG <- rbindlist(lapply(names(DB), function(pr) { d <- DB[[pr]]; o <- list(); add <- function(an, g, r) o[[length(o) + 1]] <<- cbind(data.table(predictor = pr, analysis = an, group = g), r)
  add("main: + ns(age,3) + sex + cohort within CNV", POOL, lg(d)); for (L in LOC4) add("main: + ns(age,3) + sex + cohort within CNV", L, lg(d[locus == L]))
  add("additive cohort + CNV terms (main until 2026-11-21)", POOL, additive(lg(d))); add("no cohort terms", POOL, lg(d, fe = FALSE)); add("linear age", POOL, lg(d, age = "lin")); add("age < 18 y", POOL, lg(d[age < 18])); add("age >= 18 y", POOL, lg(d[age >= 18]))
  for (L in LOC4) add("age < 18 y", L, lg(d[locus == L & age < 18]))
  for (k in sort(unique(d$cohort))) { add(paste("cohort:", k), POOL, lg(d[cohort == k], fe = FALSE)); add(paste("leave out", k), POOL, lg(d[cohort != k])) }
  for (L in LOC4) for (k in sort(unique(d[locus == L, cohort]))) add(paste("cohort:", k), L, lg(d[locus == L & cohort == k], fe = FALSE))
  add("without CHOP and SSC (diagnosis recorded by a research protocol)", POOL, lg(d[!cohort %in% c("CHOP", "SSC")]))
  ds <- DBS[[pr]]; add("with Searchlight v14 diagnoses added", POOL, lg(ds)); for (L in LOC4[1:2]) add("with Searchlight v14 diagnoses added", L, lg(ds[locus == L]))
  r <- rbindlist(o); q <- cochranQ(r[analysis == "main: + ns(age,3) + sex + cohort within CNV" & group %in% LOC4, logOR], r[analysis == "main: + ns(age,3) + sex + cohort within CNV" & group %in% LOC4, se])
  qc <- cochranQ(r[grepl("^cohort:", analysis) & group == POOL, logOR], r[grepl("^cohort:", analysis) & group == POOL, se])
  dd <- copy(d)[, cg := paste(cohort, locus)]; keep <- dd[, .(v = sum(asd) > 0 & sum(asd) < .N), by = cg][v == TRUE, cg]; dd <- dd[cg %in% keep]
  m0 <- glm(asd ~ z + ns(age, 3) + sex + cg, data = dd, family = binomial); m1 <- glm(asd ~ z:locus + ns(age, 3) + sex + cg, data = dd, family = binomial)
  rbind(r, data.table(predictor = pr, analysis = c("heterogeneity across the four CNVs (Cochran Q)", "heterogeneity across the four CNVs (LRT of predictor x CNV)", "heterogeneity across cohorts, pooled CNVs (Cochran Q)"), group = POOL,
                      n = nrow(dd), p = c(q[["p"]], anova(m0, m1, test = "Chisq")$`Pr(>Chi)`[2], qc[["p"]])), fill = TRUE) }))
LG[, note := "logistic regression of autism diagnosis on the Z score (per SD); OR per SD; cohort x CNV cells without both outcomes are dropped when cohort-within-CNV terms are fitted"]
fwrite(LG, file.path(V2OUT, "asd_by_bmi_logistic.csv"))
print(LG[, .(predictor = substr(predictor, 1, 4), analysis = substr(analysis, 1, 46), group = substr(group, 1, 13), n, n_asd, OR = round(OR, 3), lo = round(OR_lo, 2), hi = round(OR_hi, 2), p = signif(p, 3), cohorts = substr(cohorts_used, 1, 40))], nrows = 300)
print(strat[cohort == "all cohorts", .(predictor, group, band, n, n_asd, prop = round(prop, 3))], nrows = 100)
print(dcast(strat[group == POOL & cohort != "all cohorts"], predictor + cohort ~ band, value.var = c("n", "n_asd")), nrows = 100)
# ---------------------------------------------------------------- figure
GL <- c(LOC4, POOL, NCG); GLAB <- c(sub(" ", "\n", LOC4), "All four\nCNVs", "No CNV\n(SSC)")
COLS <- c("No autism diagnosis" = C_NON, "Autism diagnosis" = C_CAR)
YL <- c(Height = "Height (Z)", `Head circumference` = "Head circumference (Z)", BMI = "BMI (Z)")
r0 <- res[group %in% GL]
bars <- rbind(r0[, .(pheno, group, grp = names(COLS)[1], m = mean_noasd, se = se_noasd, n = n_noasd)], r0[, .(pheno, group, grp = names(COLS)[2], m = mean_asd, se = se_asd, n = n_asd)])
bars[, grp := factor(grp, names(COLS))][, group := factor(group, GL, GLAB)]
pan <- function(ph, tag) {
  b <- bars[pheno == ph]; r <- copy(r0[pheno == ph])[, group := factor(group, GL, GLAB)]; h <- het[pheno == ph]
  top <- b[, max(m + se, 0)]; bot <- b[, min(m - se, 0)]; rng <- top - bot
  r[, `:=`(l1 = sprintf("\u0394 %+.2f", delta), l2 = plab(p))]
  ggplot(b, aes(group, m, fill = grp)) + geom_vline(xintercept = 4.5, colour = "grey80") + geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(position = position_dodge(0.8), width = 0.75) + geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = n, y = bot - 0.07 * rng, group = grp), position = position_dodge(0.8), size = 3.6, colour = "grey25") +
    geom_text(data = r, aes(group, top + 0.24 * rng, label = l1), inherit.aes = FALSE, size = 4.4) +
    geom_text(data = r, aes(group, top + 0.12 * rng, label = l2, fontface = face(p)), inherit.aes = FALSE, size = 4.1) +
    scale_fill_manual(values = COLS, name = NULL) + coord_cartesian(ylim = c(bot - 0.1 * rng, top + 0.3 * rng)) +
    labs(x = NULL, y = YL[[ph]], tag = tag, title = sprintf("%s by autism diagnosis", ph), subtitle = sprintf("\u0394 = adjusted difference (SD); heterogeneity across the four CNVs %s", plab(h$het_p))) +
    th_v2() + theme(axis.text.x = element_text(size = 12.5))
}
pstrat <- function(pr, tag) {
  s <- strat[predictor == pr & cohort == "all cohorts"][, group := factor(group, c(POOL, LOC4))]
  l <- LG[predictor == pr & analysis == "main: + ns(age,3) + sex + cohort within CNV"][, group := factor(group, c(POOL, LOC4))]
  hq <- LG[predictor == pr & analysis == "heterogeneity across the four CNVs (Cochran Q)", p]; ymax <- max(s$ci_hi) * 100
  l[, `:=`(l1 = sprintf("OR per SD = %.2f (%.2f-%.2f)", OR, OR_lo, OR_hi), l2 = plab(p))]
  nm <- if (pr == "BMI") "BMI" else "head-circumference"
  ggplot(s, aes(band, 100 * prop, fill = band)) + geom_col(width = 0.72) + geom_errorbar(aes(ymin = 100 * ci_lo, ymax = 100 * ci_hi), width = 0.18) +
    geom_text(aes(y = -0.055 * ymax, label = sprintf("%d / %d", n_asd, n)), size = 3.9, colour = "grey20") +
    geom_text(data = l, aes(x = 2, y = ymax * 1.26, label = l1), inherit.aes = FALSE, size = 4.2) + geom_text(data = l, aes(x = 2, y = ymax * 1.14, label = l2, fontface = face(p)), inherit.aes = FALSE, size = 4.2) +
    facet_wrap(~group, nrow = 1) + scale_fill_manual(values = setNames(c(C_LO, C_MID, C_HI), BANDS), guide = "none") + coord_cartesian(ylim = c(-0.07 * ymax, ymax * 1.3)) +
    labs(x = sprintf("%s centile band", if (pr == "BMI") "BMI" else "Head-circumference"), y = "With an autism diagnosis (%)", tag = tag, title = sprintf("Autism diagnosis by %s band in CNV carriers", nm),
         subtitle = sprintf("Bars: share with a diagnosis (95%% CI), diagnosed / all below. OR per SD of %s Z from logistic regression adjusted for age, sex and cohort within CNV; heterogeneity across the four CNVs %s", nm, plab(hq))) +
    th_v2() + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 15))
}
top <- (pan("Height", "a") | pan("Head circumference", "b") | pan("BMI", "c")) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 14))
fig <- wrap_elements(full = top) / pstrat("BMI", "d") / pstrat("Head circumference", "e") + plot_layout(heights = c(1.25, 1, 1))
save_v2(fig, "fig_asd_by_cnv_v2", 20, 17.5)
fwrite(rbind(bars[, .(panel = "a-c", pheno, group = gsub("\n", " ", group), series = grp, value = m, se, n)], strat[cohort == "all cohorts", .(panel = fifelse(predictor == "BMI", "d", "e"), pheno = predictor, group, series = band, value = prop, se = NA_real_, n)]), file.path(V2OUT, "fig_asd_by_cnv_v2_plotted.csv"))
cat("85 done\n")
