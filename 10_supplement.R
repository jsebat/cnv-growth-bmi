# Supplementary bar figure (JS 2026-10-12): per-locus bars plus pooled overgrowth (16p11.2 del + 22q11.2 dup) and
# undergrowth (16p11.2 dup + 22q11.2 del) groups for
#   A. ICV (BrainChart Z) by BMI stratum, BMI <= 30th vs >= 70th centile, ages 2-40 (as B, C), 16p11.2 (BP4-5) and 22q11.2 (A-D) only.
#      Kuldeep's MRI table (Jon/Claude/MRI_CNV_growth/mri_subjects.csv), first scan per subject. The file carries no
#      breakpoint field: 16p11.2 = Clinical_16p (SVIP/CIBM, BP4-5 by ascertainment) + UKBB 16p11.2 (no distal group in the
#      file); 22q11.2 = Bearden cohort (A-D assumed). The file holds only the BMI tails, so 30/70 is the natural split.
#   B. Height by ASD vs no ASD, ages 2-40, one record per person (latest visit), genomatrix HEIGHT Z
#      (childsds WHO < 18 y; within-sample Box-Cox >= 18 y). Searchlight carriers have no linked ASD status -> excluded.
#   C. Head circumference by ASD vs no ASD, ages 2-40 (as in the main figure: < 18 y UK90 Z incl. full SVIP series;
#      >= 18 y genomatrix Box-Cox Z).
# ASD: G2MH and Cardiff = G2MH psy_asd by GUID; SVIP = v10 clinical_asd_dx; SSC carrier probands = ASD.
# Statistics: Welch t (bars) and an adjusted regression -- ICV: Z ~ group + source (+ locus for pooled);
#   ASD panels: Z ~ ASD + cohort + adult + age + age^2 (+ locus for pooled).
# Output: supp_bars_results.csv, fig_supp_bars.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data")
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
POOL <- list("Overgrowth\n(16p del + 22q dup)" = c("16p11.2 del", "22q11.2 dup"), "Undergrowth\n(16p dup + 22q del)" = c("16p11.2 dup", "22q11.2 del"))
welch <- function(lo, hi) if (length(lo) >= 2 && length(hi) >= 2) t.test(hi, lo)$p.value else NA_real_

summarise <- function(d, adj_terms, panel) {
  # d: locus, grp (0/1), y, plus covariates named in adj_terms
  one <- function(x, lab, pooled) {
    lo <- x[grp == 0, y]; hi <- x[grp == 1, y]
    terms <- c("grp", adj_terms[sapply(adj_terms, function(t) uniqueN(x[[t]]) > 1)], if (pooled) "locus")
    f <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    if ("age" %in% names(x) && "age" %in% adj_terms) f <- update(f, . ~ . + I(age^2))
    co <- tryCatch(summary(lm(f, data = x))$coefficients["grp", ], error = function(e) rep(NA, 4))
    data.table(panel = panel, group = lab, n_0 = length(lo), n_1 = length(hi), mean_0 = mean(lo), mean_1 = mean(hi),
               se_0 = sd(lo) / sqrt(length(lo)), se_1 = sd(hi) / sqrt(length(hi)), diff = mean(hi) - mean(lo), p_welch = welch(lo, hi),
               beta_adj = co[1], se_adj = co[2], p_adj = co[4])
  }
  loc <- rbindlist(lapply(LOCI, function(L) one(d[locus == L], L, FALSE)))
  # Pooled groups, locus-standardised (JS 2026-10-12): weights w_i = 1 / se_adj_i^2 (fixed-effect, k = 2 loci per group).
  #   bars     = w-weighted means of each locus's raw group means, so the gap = w-weighted mean of the within-locus raw
  #              differences (no between-locus confounding from different BMI / ASD mixes);
  #   p_welch  = z-test of that pooled raw difference (SE from the Welch SEs);
  #   beta_adj = fixed-effect meta-analysis of the per-locus adjusted betas (the covariate-adjusted pooled effect).
  #   The single pooled regression (locus as a covariate) is kept as beta_reg / p_reg for reference.
  pool <- rbindlist(lapply(names(POOL), function(P) {
    r <- loc[group %in% POOL[[P]]]; w <- (1 / r$se_adj^2) / sum(1 / r$se_adj^2)
    se_d <- sqrt(sum(w^2 * (r$se_0^2 + r$se_1^2))); dd <- sum(w * r$diff)
    b <- sum(w * r$beta_adj); se_b <- sqrt(1 / sum(1 / r$se_adj^2))
    reg <- one(d[locus %in% POOL[[P]]], P, TRUE)
    data.table(panel = panel, group = P, n_0 = sum(r$n_0), n_1 = sum(r$n_1), mean_0 = sum(w * r$mean_0), mean_1 = sum(w * r$mean_1),
               se_0 = sqrt(sum(w^2 * r$se_0^2)), se_1 = sqrt(sum(w^2 * r$se_1^2)), diff = dd, p_welch = 2 * pnorm(-abs(dd / se_d)),
               beta_adj = b, se_adj = se_b, p_adj = 2 * pnorm(-abs(b / se_b)), beta_reg = reg$beta_adj, p_reg = reg$p_adj)
  }))
  rbind(loc, pool, fill = TRUE)
}

# ---------------------------------------------------------------- A. ICV, all ages, 30/70
mri <- fread(file.path(DB, "Jon/Claude/MRI_CNV_growth/mri_subjects.csv"))[group %in% c("16p11.2_DEL", "16p11.2_DUP", "22q11.2_DEL", "22q11.2_DUP")]
mri <- mri[order(source, subject_id, age)][, .SD[1], by = .(source, subject_id)][age >= 2 & age <= 40]   # same ages as B, C (JS)
mri[, locus := sub("_", " ", sub("DEL", "del", sub("DUP", "dup", group)))]
mri[, grp := ifelse(bmi_z >= qnorm(.7), 1L, ifelse(bmi_z <= qnorm(.3), 0L, NA_integer_))]
A <- summarise(mri[!is.na(grp), .(locus, grp, y = icv_z, source)], c("source"), "A. ICV by BMI (\u2264 30th vs \u2265 70th centile), ages 2-40")

# ---------------------------------------------------------------- ASD lookup
gm <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv"))
gm <- gm[!(cohort == "G2MH" & Cardiff_mask == TRUE)][, CNV := trimws(CNV)]
gm[, id := ifelse(cohort == "SVIP", paste0("SVIP:", sampleID), paste(cohort, MSassigned_ID, sep = ":"))]
idmap <- unique(gm[, .(id, GUID, sampleID)])
s <- fread(file.path(DB, "Jon/Claude/G2MH_CNV_growth/subjects.csv"), colClasses = list(character = "asd")); g2 <- setNames(s$asd, s$sfari_id)
dx <- fread(file.path(DB, "G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/diagnosis_summary.csv"),   # JS 2026-10-21: the Autism/searchlight_multimodal copy no longer exists; same SVIP v10.0 release
            colClasses = list(character = "diagnosis_summary.clinical_asd_dx"))
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

# ---------------------------------------------------------------- B. height by ASD
source(file.path(HERE, "00_scoring.R"))   # JS 2026-10-13: everyone scored the same way (external references, from raw values)
gm[, sexc := ifelse(sex_numeric == 1, "male", ifelse(sex_numeric == 2, "female", NA))]
hgt <- gm[CNV %in% names(lab4) & !is.na(HEIGHT_raw) & age >= 2 & age <= 40 & cohort != "Searchlight", .(id, cohort, locus = lab4[CNV], age, y = score_height(HEIGHT_raw, age, sexc))][!is.na(y)]
hgt <- add_asd(hgt[order(id, -age)][, .SD[1], by = id])[, `:=`(grp = asd, adult = age >= 18)]
B <- summarise(hgt, c("cohort", "adult", "age"), "B. Height by ASD, ages 2-40")

# ---------------------------------------------------------------- C. head circumference by ASD
kid <- readRDS(file.path(DATA, "hc_visits.rds"))[age >= 2 & age < 18, .(id, cohort, locus, age, y = hcz)]
adult <- gm[CNV %in% names(lab4) & !is.na(HEADCIRCUM_raw) & age >= 18 & age <= 40, .(id, cohort, locus = lab4[CNV], age, y = score_hc(HEADCIRCUM_raw, age, sexc))][!is.na(y)]
hc <- add_asd(rbind(kid, adult)[order(id, -age)][, .SD[1], by = id])[, `:=`(grp = asd, adult = age >= 18)]
C <- summarise(hc, c("cohort", "adult", "age"), "C. Head circumference by ASD, ages 2-40")

# ---------------------------------------------------------------- non-carrier reference group for B and C (JS 2026-10-13)
# SSC non-carriers: probands (.p1/.p2) = ASD vs their siblings (.s#) = no ASD; parents excluded (other generation,
# ~no ASD, so ASD would be confounded with adulthood).
# SVIP non-carriers (family members): v10 clinical_asd_dx. Same ages, record choice, Z scales and adjusted model
# (cohort + adult + age + age^2) as the carrier groups; HC Z from the genomatrix (childsds UK90 < 18 y, Box-Cox >= 18 y).
nc_rows <- function(outcome, panel) {
  x <- gm[CNV == "No CNV" & cohort %in% c("SSC", "SVIP") & !is.na(get(outcome)) & age >= 2 & age <= 40,
          .(id, cohort, sampleID, age, y = get(outcome))][order(id, -age)][, .SD[1], by = id]
  x[, asd := NA_integer_]
  x[cohort == "SSC", asd := fifelse(grepl("[.]p[0-9]$", sampleID), 1L, fifelse(grepl("[.]s[0-9]$", sampleID), 0L, NA_integer_))]
  x[cohort == "SVIP", asd := ifelse(sv[sampleID] == "true", 1L, ifelse(sv[sampleID] == "false", 0L, NA_integer_))]
  x <- x[!is.na(asd)][, `:=`(grp = asd, adult = age >= 18, locus = "Non-carriers")]
  lo <- x[grp == 0, y]; hi <- x[grp == 1, y]
  co <- summary(lm(y ~ grp + cohort + adult + age + I(age^2), data = x))$coefficients["grp", ]
  data.table(panel = panel, group = "Non-carriers", n_0 = length(lo), n_1 = length(hi), mean_0 = mean(lo), mean_1 = mean(hi),
             se_0 = sd(lo) / sqrt(length(lo)), se_1 = sd(hi) / sqrt(length(hi)), diff = mean(hi) - mean(lo), p_welch = welch(lo, hi),
             beta_adj = co[1], se_adj = co[2], p_adj = co[4])
}
# JS 2026-10-13: non-carrier group = within-family SSC proband vs designated-sibling pairs (18_ssc_sibpairs.R).
# nc_rows() above is kept for reference but not used (SSC siblings have no height/HC in the genomatrix).
sp <- fread(file.path(HERE, "ssc_sibpairs_results.csv"))
sib_row <- function(o, panel) sp[outcome == o, .(panel = panel, group = "Non-carriers", n_0 = pairs, n_1 = pairs, mean_0 = mean_sibling, mean_1 = mean_proband,
  se_0 = se_sibling, se_1 = se_proband, diff = paired_diff, p_welch = paired_p, beta_adj = paired_diff, se_adj = paired_se, p_adj = paired_p, within_family = TRUE)]
B <- rbind(sib_row("Height", B$panel[1]), B, fill = TRUE)
C <- rbind(sib_row("Head circumference", C$panel[1]), C, fill = TRUE)

res <- rbind(A, B, C, fill = TRUE)
fwrite(res, file.path(HERE, "supp_bars_results.csv"))
print(res[, .(panel = substr(panel, 1, 2), group = gsub("\n", " ", group), n_0, n_1, diff = round(diff, 2), p_welch = signif(p_welch, 2),
              beta_adj = round(beta_adj, 2), p_adj = signif(p_adj, 2))])

# ---------------------------------------------------------------- figure
pf <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))
bar_panel <- function(r, labs, cols, ylab, rev_legend = FALSE) {
  r <- copy(r)[, group := factor(group, intersect(c("Non-carriers", LOCI, names(POOL)), group))]
  d <- rbind(r[, .(group, grp = labs[1], m = mean_0, se = se_0, n = n_0)], r[, .(group, grp = labs[2], m = mean_1, se = se_1, n = n_1)])
  d[, grp := factor(grp, labs)]
  r[, lab := sprintf("\u0394 %+.2f, p = %s\nadj. %+.2f, p = %s", diff, pf(p_welch), beta_adj, pf(p_adj))]
  if ("within_family" %in% names(r)) r[within_family %in% TRUE, lab := sprintf("sib-pair \u0394 %+.2f\np = %s", diff, pf(p_welch))]
  r[, bold := !is.na(p_adj) & p_adj < 0.05]
  yr <- range(c(d$m - d$se, d$m + d$se), na.rm = TRUE); top <- yr[2] + 0.25 * diff(yr); bot <- yr[1] - 0.12 * diff(yr)
  ggplot(d, aes(group, m, fill = grp)) + geom_hline(yintercept = 0, colour = "grey50") +
    annotate("rect", xmin = which(levels(r$group) == names(POOL)[1]) - 0.5, xmax = nlevels(r$group) + 0.6, ymin = -Inf, ymax = Inf, fill = "grey95") +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2, na.rm = TRUE) +
    geom_text(aes(y = bot, label = paste0("n=", n)), position = position_dodge(0.8), size = 3.6) +
    geom_text(data = r, aes(group, top, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 3.8, vjust = 1) +
    scale_fill_manual(values = cols) + guides(fill = guide_legend(reverse = rev_legend)) + coord_cartesian(ylim = c(bot, top)) +
    labs(x = NULL, y = ylab, fill = NULL, title = r$panel[1]) +
    theme_classic(base_size = 15) + theme(legend.position = "top", legend.text = element_text(size = 8), legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 15.75))
}
pA <- bar_panel(A, c("BMI \u2264 30th centile", "BMI \u2265 70th centile"), c("#00A5AD", "#F0645A"), "ICV Z (BrainChart)", rev_legend = TRUE)
pB <- bar_panel(B, c("No ASD", "ASD"), c("#8C9A9A", "#F28E2B"), "Height Z")
pC <- bar_panel(C, c("No ASD", "ASD"), c("#8C9A9A", "#F28E2B"), "Head circumference Z")
fig <- pA / pB / pC + plot_annotation(theme = theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 11), plot.caption = element_text(size = 9)),
  title = "Supplementary: ICV by BMI and growth by ASD in 16p11.2 and 22q11.2 carriers, per locus and pooled by growth direction",
  caption = paste0("Bars = mean \u00b1 SE; n under bars; one record per person. Pooled (shaded): overgrowth = 16p11.2 del + 22q11.2 dup; undergrowth = 16p11.2 dup + 22q11.2 del.\n",
                   "Pooled bars are locus-standardised (inverse-variance-weighted locus means), so their gap is the pooled within-locus difference; pooled adj. = fixed-effect meta-analysis of per-locus adjusted effects.\n",
                   "\u0394 and Welch p = higher-BMI or ASD minus comparison group (pooled: z-test). adj. = regression adjusting for MRI source (A) or cohort, adult, age, age\u00b2 (B, C). Bold = adj. p < 0.05.\n",
                   "A: Kuldeep MRI table, first scan per subject, BMI tails only; 16p11.2 = SVIP/CIBM (BP4-5) + UKBB, 22q11.2 = Bearden cohort (A-D assumed; no breakpoint field).\n",
                   "B, C: G2MH, Cardiff, SVIP (full v10 series for HC < 18 y), SSC; ASD from G2MH psy_asd, SVIP clinical_asd_dx, SSC probands. Z < 18 y growth references; \u2265 18 y within-sample Box-Cox."))
ggsave(file.path(HERE, "fig_supp_bars.png"), fig, width = 14, height = 15, dpi = 200)
ggsave(file.path(HERE, "fig_supp_bars.pdf"), fig, width = 14, height = 15, device = cairo_pdf)
