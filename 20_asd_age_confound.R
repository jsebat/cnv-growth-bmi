# Is the ASD vs no-ASD difference in head circumference (panel L) and height (panel F) confounded with age? (JS 2026-10-24)
# ASD status in G2MH/Cardiff/SVIP comes mostly from probands, who are children, while non-ASD carriers include adult relatives,
# so an age trend in the Z scores could masquerade as an ASD effect. Panels F/L adjust for cohort, adult, age and age^2 within
# each locus and pool the loci (10_supplement.R). Here, same records as F/L:
#   1. age distribution by ASD within each locus and pool
#   2. ASD effect within children (< 18 y) and within adults (>= 18 y) separately, per locus and pooled (fixed-effect over loci)
#   3. GAM adjustment: y ~ asd + s(age) + cohort within each locus (REML), pooled by fixed effect
#   4. age-matched check: ASD effect in the age window where both groups have support (5th-95th percentile overlap)
# Output: asd_age_confound.csv, fig_asd_age_confound.{png,pdf}. Runs locally (clinical data only).
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")
source(file.path(HERE, "00_scoring.R"))
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
POOL <- list("Overgrowth CNVs" = c("16p11.2 del", "22q11.2 dup"), "Undergrowth CNVs" = c("16p11.2 dup", "22q11.2 del"))
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))

# ---- same records as 10_supplement.R panels B/C ----
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
gm <- gm[!(cohort == "G2MH" & Cardiff_mask == TRUE)][, CNV := trimws(CNV)]
gm[, id := ifelse(cohort == "SVIP", paste0("SVIP:", sampleID), paste(cohort, MSassigned_ID, sep = ":"))]
gm[, sexc := ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA))]
idmap <- unique(gm[, .(id, GUID, sampleID)])
s <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd")); g2 <- setNames(s$asd, s$sfari_id)
dx <- fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"), colClasses = list(character = "diagnosis_summary.clinical_asd_dx"))
sv <- setNames(tolower(dx$diagnosis_summary.clinical_asd_dx), dx$individual)
add_asd <- function(g) {
  g <- merge(g, idmap, by = "id", all.x = TRUE); g[cohort == "SVIP", sampleID := sub("^SVIP:", "", id)]
  g[, asd := NA_integer_]
  g[cohort %in% c("G2MH", "Cardiff"), asd := suppressWarnings(as.integer(g2[GUID]))]
  g[cohort == "SVIP", asd := ifelse(sv[sampleID] == "true", 1L, ifelse(sv[sampleID] == "false", 0L, NA_integer_))]
  g[cohort == "SSC" & grepl("[.]p1$", sampleID), asd := 1L]
  g[!is.na(asd)]
}
lab4 <- c("16p11.2 proximal DEL" = "16p11.2 del", "16p11.2 proximal DUP" = "16p11.2 dup", "22q11.2 A-D DEL" = "22q11.2 del", "22q11.2 A-D DUP" = "22q11.2 dup")
hgt <- gm[CNV %in% names(lab4) & !is.na(HEIGHT_raw) & age >= 2 & age <= 40 & cohort != "Searchlight", .(id, cohort, locus = lab4[CNV], age, y = score_height(HEIGHT_raw, age, sexc))][!is.na(y)]
hgt <- add_asd(hgt[order(id, -age)][, .SD[1], by = id])[, outcome := "Height"]
kid <- readRDS(file.path(DATA, "hc_visits.rds"))[age >= 2 & age < 18, .(id, cohort, locus, age, y = hcz)]
adult <- gm[CNV %in% names(lab4) & !is.na(HEADCIRCUM_raw) & age >= 18 & age <= 40, .(id, cohort, locus = lab4[CNV], age, y = score_hc(HEADCIRCUM_raw, age, sexc))][!is.na(y)]
hc <- add_asd(rbind(kid, adult)[order(id, -age)][, .SD[1], by = id])[, outcome := "Head circumference"]
d <- rbind(hgt, hc)[, `:=`(adult = age >= 18, asd_lab = ifelse(asd == 1, "ASD", "No ASD"))]

# ---- 1. age by ASD ----
agetab <- d[, .(n = .N, age_median = median(age), age_q25 = quantile(age, .25), age_q75 = quantile(age, .75), pct_adult = 100 * mean(adult)), by = .(outcome, locus, asd_lab)][order(outcome, locus, asd_lab)]
cat("== age by ASD status\n"); print(agetab[, .(outcome, locus, asd_lab, n, age_median, pct_adult = round(pct_adult))])

# ---- 2-4. ASD effect under different age handling ----
fe <- function(b, se) { w <- 1 / se^2; data.table(beta = sum(w * b) / sum(w), se = sqrt(1 / sum(w))) }
est <- function(x, how) {
  if (nrow(x) < 8 || uniqueN(x$asd) < 2) return(c(NA, NA, x[asd == 0, .N], x[asd == 1, .N]))
  cov_cohort <- if (uniqueN(x$cohort) > 1) "+ cohort" else ""
  co <- tryCatch(switch(how,
    raw      = { m <- lm(y ~ asd, data = x); summary(m)$coefficients["asd", ] },
    adj_FL   = { f <- as.formula(paste("y ~ asd", cov_cohort, if (uniqueN(x$adult) > 1) "+ adult" else "", "+ age + I(age^2)")); summary(lm(f, data = x))$coefficients["asd", ] },
    gam      = { f <- as.formula(paste("y ~ asd", cov_cohort, "+ s(age, k = 6)")); summary(gam(f, data = x, method = "REML"))$p.table["asd", ] },
    children = { xx <- x[adult == FALSE]; if (nrow(xx) < 8 || uniqueN(xx$asd) < 2) stop(); f <- as.formula(paste("y ~ asd", if (uniqueN(xx$cohort) > 1) "+ cohort" else "", "+ age + I(age^2)")); summary(lm(f, data = xx))$coefficients["asd", ] },
    adults   = { xx <- x[adult == TRUE];  if (nrow(xx) < 8 || uniqueN(xx$asd) < 2) stop(); f <- as.formula(paste("y ~ asd", if (uniqueN(xx$cohort) > 1) "+ cohort" else "", "+ age")); summary(lm(f, data = xx))$coefficients["asd", ] },
    overlap  = { lo <- max(quantile(x[asd == 0, age], .05), quantile(x[asd == 1, age], .05)); hi <- min(quantile(x[asd == 0, age], .95), quantile(x[asd == 1, age], .95))
                 xx <- x[age >= lo & age <= hi]; if (nrow(xx) < 8 || uniqueN(xx$asd) < 2) stop(); f <- as.formula(paste("y ~ asd", if (uniqueN(xx$cohort) > 1) "+ cohort" else "", "+ age + I(age^2)")); summary(lm(f, data = xx))$coefficients["asd", ] }),
    error = function(e) rep(NA_real_, 4))
  n <- switch(how, children = x[adult == FALSE, .(sum(asd == 0), sum(asd == 1))], adults = x[adult == TRUE, .(sum(asd == 0), sum(asd == 1))], x[, .(sum(asd == 0), sum(asd == 1))])
  c(co[1], co[2], unlist(n))
}
HOWS <- c(raw = "Unadjusted", adj_FL = "F/L adjustment (cohort, adult, age, age²)", gam = "GAM: cohort + s(age)", children = "Children only (< 18 y)", adults = "Adults only (18-40 y)", overlap = "Age-overlap window only")
res <- rbindlist(lapply(unique(d$outcome), function(o) rbindlist(lapply(names(HOWS), function(h) {
  loc <- rbindlist(lapply(LOCI, function(L) { e <- est(d[outcome == o & locus == L], h); data.table(outcome = o, method = HOWS[[h]], group = L, beta = e[1], se = e[2], n_noasd = e[3], n_asd = e[4]) }))
  pool <- rbindlist(lapply(names(POOL), function(P) { r <- loc[group %in% POOL[[P]] & is.finite(se)]; if (nrow(r) == 0) return(NULL)
    f <- fe(r$beta, r$se); data.table(outcome = o, method = HOWS[[h]], group = P, beta = f$beta, se = f$se, n_noasd = sum(r$n_noasd), n_asd = sum(r$n_asd)) }))
  rbind(loc, pool)
}))))
res[, p := 2 * pnorm(-abs(beta / se))][, method := factor(method, HOWS)]
fwrite(res, file.path(HERE, "asd_age_confound.csv"))
cat("\n== ASD effect (ASD minus no ASD, SD) by age handling, pooled groups\n")
print(res[group %in% names(POOL), .(outcome, group, method, n_noasd, n_asd, beta = round(beta, 2), se = round(se, 2), p = signif(p, 2))][order(outcome, group, method)], nrows = 100)
cat("\n== per locus, head circumference\n"); print(res[outcome == "Head circumference" & group %in% LOCI, .(group, method, n_noasd, n_asd, beta = round(beta, 2), p = signif(p, 2))][order(group, method)], nrows = 100)

# ---- figure ----
d[, pool := fcase(locus %in% POOL[[1]], names(POOL)[1], locus %in% POOL[[2]], names(POOL)[2])]
ACOL <- c("No ASD" = "#8C9A9A", "ASD" = "#F28E2B")
pA <- ggplot(d, aes(age, fill = asd_lab)) + geom_histogram(binwidth = 2, position = "identity", alpha = 0.55, colour = NA) + geom_vline(xintercept = 18, linetype = 2, colour = "grey40") +
  facet_grid(outcome ~ pool) + scale_fill_manual(values = ACOL, name = NULL) + labs(x = "Age (years)", y = "Carriers", subtitle = "Age distribution by ASD status") +
  theme_classic(base_size = 13) + theme(legend.position = "top", legend.justification = "left", strip.background = element_blank())
pB <- ggplot(d, aes(age, y, colour = asd_lab, fill = asd_lab)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") + geom_point(size = 1, alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) + geom_vline(xintercept = 18, linetype = 3, colour = "grey60") +
  facet_grid(outcome ~ pool) + scale_colour_manual(values = ACOL, name = NULL) + scale_fill_manual(values = ACOL, name = NULL) + coord_cartesian(ylim = c(-3.5, 3.5)) +
  labs(x = "Age (years)", y = "Z (SD)", subtitle = "Z by age and ASD status; GAM per group") + theme_classic(base_size = 13) + theme(legend.position = "top", legend.justification = "left", strip.background = element_blank())
rp <- res[group %in% names(POOL)][, group := factor(group, names(POOL))][, method := factor(method, rev(HOWS))]
pC <- ggplot(rp, aes(beta, method, colour = outcome)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se), position = position_dodge(width = 0.5), size = 0.5) +
  geom_text(aes(label = sprintf("n = %d / %d", n_noasd, n_asd), x = 1.55), position = position_dodge(width = 0.5), size = 2.6, hjust = 1, show.legend = FALSE) +
  facet_wrap(~ group, nrow = 1) + scale_colour_manual(values = c(Height = "#4E79A7", `Head circumference` = "#E15759"), name = NULL) + coord_cartesian(xlim = c(-0.8, 1.6)) +
  labs(x = "ASD minus no ASD (SD), 95% CI", y = NULL, subtitle = "Pooled (fixed effect over the two loci) ASD effect under different handling of age; n = no ASD / ASD") +
  theme_classic(base_size = 13) + theme(legend.position = "top", legend.justification = "left", strip.background = element_blank())
fig <- (pA | pB) / pC + plot_layout(heights = c(1.1, 0.8)) + plot_annotation(tag_levels = "A", title = "Is the ASD difference in height / head circumference (panels F, L) confounded with age?") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
ggsave(file.path(HERE, "fig_asd_age_confound.png"), fig, width = 15, height = 11, dpi = 170)
tryCatch(ggsave(file.path(HERE, "fig_asd_age_confound.pdf"), fig, width = 15, height = 11, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_asd_age_confound.pdf"), fig, width = 15, height = 11))
cat("done\n")
