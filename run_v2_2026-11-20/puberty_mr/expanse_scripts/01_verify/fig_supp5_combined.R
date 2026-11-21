# Combined Supplementary Fig. 5 (JS 2026-11-02, Jon: "combine most of this into 1 figure"): body mass in childhood vs adulthood,
# pubertal timing, and genetic instruments, all in UK Biobank CNV carriers (EUR/blue cluster, QC'd).
#  a  adult height per ordered step of recalled body size at 10 vs adult BMI group (ukbb_bodysize/ukbb_bodysize10_vs_adult_trends.csv)
#  b  age at menarche by genotype (women; field 2714, first valid answer 8-18 y): mean difference from non-carriers, adjusted for
#     age and PC1-10, 95% CI, P
#  c  adult height (Z) per year of later menarche within genotype, 95% CI, P; interaction vs non-carriers
#  d  adult height (Z) per s.d. of polygenic score within genotype (childhood obesity PGS, adult BMI PGS; + age-at-menarche PGS
#     when jsebat/pgs_menarche/ukbb/PGS_guan_2025_ext_agemenarche_genomewide.tsv exists), height PGS as covariate; 95% CI, P;
#     first-stage association of each score with the exposure printed in the table (fig_supp5_combined_tables.csv)
# Groups: 16p11.2 DEL, No CNV, 16p11.2 DUP, 22q11.2 DUP (22q11.2 DEL: 5 carriers, omitted). Output: fig_supp5_combined.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
WD <- getwd(); options(width = 220)
GM <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"; MS <- "/expanse/projects/sebat1/msacks/UKBB"
BS <- "/expanse/projects/sebat1/jsebat/growth_tables_20260920/ukbb_bodysize/ukbb_bodysize10_vs_adult_trends.csv"
MR <- "/expanse/projects/sebat1/jsebat/growth_tables_20260920/ukbb_mr/ukbb_mr_childhood_obesity.csv"
AAM <- "/expanse/projects/sebat1/jsebat/pgs_menarche/ukbb/PGS_guan_2025_ext_agemenarche_cojo_genomewide.tsv"
GROUPS <- c("16p11.2 DEL", "No CNV", "16p11.2 DUP", "22q11.2 DUP")
GCOL <- c("16p11.2 DEL" = "#C0392B", "No CNV" = "grey40", "16p11.2 DUP" = "#2E86C1", "22q11.2 DUP" = "#8E44AD")
pfx <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-300, "p < 1e-300", ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.2g", p))))
th <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5, colour = "grey30"), legend.position = "top", legend.key.size = unit(0.35, "cm"))
# ---------------------------------------------------------------- data (as in ukbb_bodysize10_vs_adult.R / ukbb_mr_childhood_obesity.R)
gm <- fread(GM, select = c("project_40980_IID", "f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", paste0("PC", 1:10), "Height", "BMI_raw", "age", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"))
setnames(gm, c("project_40980_IID", "f.eid", "sex_combined", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"), c("IID", "eid", "sex", "l16", "l22"))
gm <- gm[duplicate_drop == "no" & withdrawn_consent == "no" & ancestry == "blue" & sex %in% 1:2 & !is.na(Height) & !is.na(age) & l16 %in% 1:3 & l22 %in% 1:3]
gm <- gm[complete.cases(gm[, paste0("PC", 1:10), with = FALSE])]
gm[, locus := fcase(l16 == 1, "16p11.2 DEL", l16 == 3, "16p11.2 DUP", l22 == 1, "22q11.2 DEL", l22 == 3, "22q11.2 DUP", l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
gm <- gm[!is.na(locus) & locus %in% GROUPS][, locus := factor(locus, GROUPS)]
gm[, age_sq := age^2]; gm[, sexf := factor(sex)]; pcs <- paste0("PC", 1:10, collapse = " + ")
gm[, hZ := { r <- resid(lm(as.formula(paste("Height ~ sexf + age + age_sq +", pcs)), gm)); r / sd(r) }]
men <- fread(file.path(MS, "comparitiveBodySize/field2714_age_at_menarche.tsv")); setnames(men, "f.eid", "eid")
inst <- grep("^f\\.2714\\.", names(men), value = TRUE); men[, menarche := NA_real_]; for (v in inst) men[is.na(menarche) & get(v) >= 8 & get(v) <= 18, menarche := as.numeric(get(v))]
w <- merge(gm[sex == 2], men[!is.na(menarche), .(eid, menarche)], by = "eid")
bs1 <- fread(file.path(MS, "comparitiveBodySize/field1687_comparative_body_size_age10.tsv")); setnames(bs1, "f.eid", "eid")
i1 <- grep("^f[.]1687[.]", names(bs1), value = TRUE); bs1[, code := NA_integer_]; for (v in i1) bs1[is.na(code) & get(v) %in% 1:3, code := as.integer(get(v))]
a1 <- merge(gm[!is.na(BMI_raw)], bs1[!is.na(code), .(eid, body_size_10 = factor(code, levels = c(1, 3, 2), labels = c("Thinner", "About average", "Plumper")))], by = "eid")
a1[, bZ := { r <- resid(lm(as.formula(paste("log(BMI_raw) ~ sexf + age + age_sq +", pcs)), a1)); r / sd(r) }]
a1[, bmi_pct := frank(bZ, ties.method = "average") / .N][, bmi_group := factor(fifelse(bmi_pct <= 0.30, "Low", fifelse(bmi_pct >= 0.70, "High", "Middle")), levels = c("Low", "Middle", "High"))]
al <- rbind(a1[, .(group = locus, hZ, measure = "Body size at age 10", level = as.character(body_size_10))], a1[, .(group = locus, hZ, measure = "Adult BMI group", level = as.character(bmi_group))])
al[, level := factor(level, c("Thinner", "About average", "Plumper", "Low", "Middle", "High"))][, measure := factor(measure, c("Body size at age 10", "Adult BMI group"))]
set.seed(1); alp <- rbind(al[group != "No CNV"], al[group == "No CNV"][sample(.N, 600)])   # points: all carriers, 600 random non-carriers
LCOL <- c(Thinner = "#1AA6B0", `About average` = "grey55", Plumper = "#F0645A", Low = "#1AA6B0", Middle = "grey55", High = "#F0645A")
cat("women with menarche:", nrow(w), "\n"); print(w[, .(n = .N, mean_menarche = round(mean(menarche), 2)), by = locus])
# ---------------------------------------------------------------- b: menarche by genotype
w[, locus0 := relevel(locus, "No CNV")]; mb <- lm(as.formula(paste("menarche ~ locus0 + age +", pcs)), w); s <- summary(mb)$coefficients; ci <- confint(mb)
tb <- rbindlist(lapply(GROUPS, function(G) { if (G == "No CNV") data.table(group = G, n = w[locus == G, .N], mean = w[locus == G, mean(menarche)], diff = 0, lo = NA_real_, hi = NA_real_, p = NA_real_) else { k <- paste0("locus0", G); data.table(group = G, n = w[locus == G, .N], mean = w[locus == G, mean(menarche)], diff = s[k, 1], lo = ci[k, 1], hi = ci[k, 2], p = s[k, 4]) } }))
tb[, group := factor(group, GROUPS)]
# ---------------------------------------------------------------- c: height per year of menarche within genotype
tc <- rbindlist(lapply(GROUPS, function(G) { x <- w[locus == G]; m <- lm(hZ ~ menarche, x); s <- summary(m)$coefficients; ci <- confint(m)
  ip <- if (G == "No CNV") NA_real_ else summary(lm(hZ ~ menarche * I(locus == G), w[locus %in% c(G, "No CNV")]))$coefficients["menarche:I(locus == G)TRUE", 4]
  data.table(group = G, n = nrow(x), slope = s["menarche", 1], lo = ci["menarche", 1], hi = ci["menarche", 2], p = s["menarche", 4], p_interaction = ip) }))
tc[, group := factor(group, GROUPS)]
# ---------------------------------------------------------------- d: height per s.d. of PGS within genotype (from the MR table; + menarche PGS if scored)
mr <- fread(MR)
td <- rbind(mr[group %in% GROUPS, .(group, pgs, n = n_carriers, est = height_per_pgs_carriers, se = se_carriers, p = p_carriers, fs = fs_step10_carriers, p_int = interaction_p)],
            unique(mr[, .(group = "No CNV", pgs, n = n_noncarriers, est = height_per_pgs_noncarriers, se = NA_real_, p = p_noncarriers, fs = fs_step10_noncarriers, p_int = NA_real_)]))
USE_AAM <- TRUE   # JS 2026-11-03: menarche PGS included, plus panel e (MR of age at menarche on adult height). Note: the Guan 2025 ext
                  # weights are the within-family estimate (N ~11k), so the instrument is weak; its first stage is reported.
if (USE_AAM && file.exists(AAM)) {
  aam <- fread(AAM)[n_chr == 22, .(IID, pgs_aam = PGS)]; dd <- merge(gm, aam, by = "IID")
  ht <- fread(file.path(MS, "comparitiveBodySize/PGS_Height_genomewide.tsv"))[n_chr == 22, .(IID, pgs_ht = PGS)]; dd <- merge(dd, ht, by = "IID")
  for (v in c("pgs_aam", "pgs_ht")) { r <- resid(lm(as.formula(paste(v, "~", pcs)), dd)); set(dd, j = v, value = r / sd(r)) }
  dd <- merge(dd, men[!is.na(menarche), .(eid, menarche)], by = "eid", all.x = TRUE)
  td <- rbind(td, rbindlist(lapply(GROUPS, function(G) { x <- dd[locus == G]; s <- summary(lm(hZ ~ pgs_aam + pgs_ht, x))$coefficients
    fs <- if (x[!is.na(menarche), .N] > 10) coef(lm(menarche ~ pgs_aam, x[!is.na(menarche)]))["pgs_aam"] else NA_real_
    ip <- if (G == "No CNV") NA_real_ else summary(lm(hZ ~ pgs_aam * I(locus == G) + pgs_ht * I(locus == G), dd[locus %in% c(G, "No CNV")]))$coefficients["pgs_aam:I(locus == G)TRUE", 4]
    data.table(group = G, pgs = "age at menarche PGS", n = nrow(x), est = s["pgs_aam", 1], se = s["pgs_aam", 2], p = s["pgs_aam", 4], fs = fs, p_int = ip) })))
  # ---- e: Mendelian randomization of age at menarche on adult height, women, within genotype.
  # Wald ratio = reduced form (height ~ menarche PGS + height PGS, within genotype) / first stage (menarche ~ menarche PGS
  # + height PGS + age + PCs) estimated in all women (the carrier groups are too small for their own first stage; theirs is
  # reported as fs_group). Delta-method SE. p_diff: reduced-form difference from non-carriers (PGS x genotype interaction).
  dw <- dd[sex == 2 & !is.na(menarche)]
  fsa <- summary(lm(as.formula(paste("menarche ~ pgs_aam + pgs_ht + age +", pcs)), dw))$coefficients["pgs_aam", ]
  te <- rbindlist(lapply(GROUPS, function(G) { x <- dw[locus == G]
    rf <- summary(lm(hZ ~ pgs_aam + pgs_ht, x))$coefficients["pgs_aam", ]
    fg <- summary(lm(menarche ~ pgs_aam, x))$coefficients["pgs_aam", ]
    est <- rf[1] / fsa[1]; se <- sqrt(rf[2]^2 / fsa[1]^2 + rf[1]^2 * fsa[2]^2 / fsa[1]^4)
    pdf_ <- if (G == "No CNV") NA_real_ else summary(lm(hZ ~ pgs_aam * I(locus == G) + pgs_ht * I(locus == G), dw[locus %in% c(G, "No CNV")]))$coefficients["pgs_aam:I(locus == G)TRUE", 4]
    data.table(group = G, n = nrow(x), mr = est, mr_se = se, mr_p = 2 * pnorm(-abs(est / se)), rf = rf[1], rf_se = rf[2],
               fs_group = fg[1], fs_group_se = fg[2], p_diff = pdf_) }))
  te[, `:=`(fs_all = fsa[1], fs_all_se = fsa[2], fs_all_F = (fsa[1] / fsa[2])^2, n_women_all = nrow(dw))]
  fwrite(te, file.path(WD, "mr_menarche_height_ukbb.csv"))
  cat("menarche PGS included; first stage in all women:", round(fsa[1], 3), "years per s.d., F =", round((fsa[1] / fsa[2])^2), "\n"); print(te)
} else cat("menarche PGS not yet scored (", AAM, "): panel d has childhood obesity and adult BMI PGS only\n")
PGSL <- c("childhood obesity PGS", "adult BMI PGS", "age at menarche PGS"); td <- td[pgs %in% PGSL][, pgs := factor(pgs, PGSL)][, group := factor(group, GROUPS)]
td[is.na(se) & group == "No CNV", se := abs(est / qnorm(p / 2))]
# ---------------------------------------------------------------- panel a: box plots by level (Molly's format), trend bracket from the existing table
ta <- fread(BS)[group %in% GROUPS][, group := factor(group, GROUPS)][, measure := factor(measure, c("Body size at age 10", "Adult BMI group"))]
ta[, lab := sprintf("%+.2f SD/step
%s", trend, pfx(p))][, bold := p < 0.05]
nl <- al[, .N, by = .(group, measure, level)]
ytop <- 3.6
pa <- ggplot(al, aes(level, hZ)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(aes(fill = level), outlier.shape = NA, width = 0.7, colour = "grey25", linewidth = 0.3) +
  geom_jitter(data = alp, aes(colour = level), width = 0.18, size = 0.6, shape = 1, alpha = 0.6, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2, fill = "white", colour = "grey20") +
  geom_text(data = nl, aes(level, -3.55, label = paste0("n=", ifelse(N >= 10000, paste0(round(N / 1000), "k"), N))), size = 2.3, colour = "grey30") +
  geom_segment(data = ta, aes(x = 1, xend = 3, y = ytop - 0.25, yend = ytop - 0.25), inherit.aes = FALSE, colour = "grey30", linewidth = 0.3) +
  geom_text(data = ta, aes(x = 2, y = ytop + 0.05, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 2.5, vjust = 0, lineheight = 0.9) +
  facet_grid(measure ~ group, scales = "free_x") +
  scale_fill_manual(values = LCOL, guide = "none") + scale_colour_manual(values = LCOL, guide = "none") +
  coord_cartesian(ylim = c(-3.8, 4.5)) +
  labs(title = "a  Body size at age 10 vs adult BMI group, and adult height", subtitle = "same individuals; box = all, points = all carriers and 600 random non-carriers; diamond = mean; bracket = linear trend per step", x = NULL, y = "Adult height (SD)") +
  th + theme(axis.text.x = element_text(size = 7.5, angle = 30, hjust = 1), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9))
# ---------------------------------------------------------------- panel b: menarche by genotype, box plots
set.seed(2); wp <- rbind(w[locus != "No CNV"], w[locus == "No CNV"][sample(.N, 600)])
tb[, lab := ifelse(is.na(p), "", sprintf("%+.2f y
%s", diff, pfx(p)))][, bold := !is.na(p) & p < 0.05]
pb <- ggplot(w, aes(x = locus, y = menarche)) +
  geom_boxplot(aes(fill = locus), outlier.shape = NA, width = 0.7, colour = "grey25", linewidth = 0.3) +
  geom_jitter(data = wp, aes(colour = locus), width = 0.18, height = 0.15, size = 0.6, shape = 1, alpha = 0.6, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2, fill = "white", colour = "grey20") +
  geom_text(data = tb, aes(group, 18.9, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 2.5, vjust = 0, lineheight = 0.9) +
  geom_text(data = tb, aes(group, 7.3, label = paste0("n=", ifelse(n >= 10000, paste0(round(n / 1000), "k"), n), "
mean ", sprintf("%.1f", mean))), inherit.aes = FALSE, size = 2.3, colour = "grey30", lineheight = 0.9) +
  scale_fill_manual(values = GCOL, guide = "none") + scale_colour_manual(values = GCOL, guide = "none") +
  coord_cartesian(ylim = c(6.8, 20.2)) +
  labs(title = "b  Age at menarche by genotype", subtitle = "women; difference from non-carriers adjusted for
age and PC1-10 (years), diamond = mean", x = NULL, y = "Age at menarche (years)") + th
# ---------------------------------------------------------------- panel c
pc <- ggplot(tc, aes(group, slope, colour = group)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.7) + geom_point(size = 3) +
  geom_text(aes(y = hi + 0.01, label = ifelse(is.na(p_interaction), pfx(p), paste0(pfx(p), "\nvs non-carriers ", pfx(p_interaction))), fontface = ifelse(p < 0.05, "bold", "plain")), size = 2.4, vjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = GCOL, guide = "none") + scale_x_discrete(labels = setNames(paste0(levels(tc$group), "\nn = ", tc[match(levels(tc$group), group), ifelse(n >= 10000, paste0(round(n / 1000), "k"), n)]), levels(tc$group))) +
  coord_cartesian(ylim = c(min(tc$lo) - 0.02, max(tc$hi) + 0.16)) +
  labs(title = "c  Later menarche, taller adult", subtitle = "adult height Z per year of later menarche
within genotype, 95% CI", x = NULL, y = "Adult height (SD) per year") + th
# ---------------------------------------------------------------- panel d
pdd <- position_dodge(width = 0.65)
pd4 <- ggplot(td, aes(group, est, colour = pgs)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), position = pdd, width = 0.2, linewidth = 0.7) + geom_point(position = pdd, size = 3) +
  geom_text(aes(y = est + 1.96 * se + 0.02, label = pfx(p), fontface = ifelse(p < 0.05, "bold", "plain"), group = pgs), position = pdd, size = 2.3, vjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = c("childhood obesity PGS" = "#F28E2B", "adult BMI PGS" = "#4B6EA9", "age at menarche PGS" = "#27AE60"), name = NULL) +
  scale_x_discrete(labels = setNames(paste0(levels(td$group), "\nn = ", td[pgs == "childhood obesity PGS"][match(levels(td$group), group), ifelse(n >= 10000, paste0(round(n / 1000), "k"), n)]), levels(td$group))) +
  labs(title = "d  Genetic instruments", subtitle = "adult height per s.d. of polygenic score within genotype,
adjusted for height PGS and PC1-10, 95% CI", x = NULL, y = "Adult height (SD) per s.d. of PGS") + th
# ---------------------------------------------------------------- panel e: observational vs MR estimate of later menarche on adult height
ee <- rbind(tc[, .(group, type = "observational", est = slope, lo, hi, p)],
            te[, .(group = factor(group, GROUPS), type = "MR (menarche PGS)", est = mr, lo = mr - 1.96 * mr_se, hi = mr + 1.96 * mr_se, p = mr_p)])
ee[, type := factor(type, c("observational", "MR (menarche PGS)"))]
pe <- ggplot(ee, aes(group, est, shape = type, colour = group)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = pdd, width = 0.2, linewidth = 0.7) + geom_point(position = pdd, size = 3, fill = "white") +
  geom_text(aes(y = hi + 0.02, label = pfx(p), fontface = ifelse(p < 0.05, "bold", "plain"), group = type), position = pdd, size = 2.3, vjust = 0, show.legend = FALSE) +
  scale_shape_manual(values = c("observational" = 16, "MR (menarche PGS)" = 24), name = NULL) + scale_colour_manual(values = GCOL, guide = "none") +
  labs(title = "e  Menarche and adult height: observational vs MR", subtitle = sprintf("women; adult height Z per year of later menarche, 95%% CI;
MR first stage in all women %.2f y per s.d. of PGS (F = %.0f)", te$fs_all[1], te$fs_all_F[1]),
       x = NULL, y = "Adult height (SD) per year") + th
fig <- pa / (pb | pc) / (pd4 | pe) + plot_layout(heights = c(1.25, 1, 1))
ggsave(file.path(WD, "fig_supp5_combined.png"), fig, width = 14, height = 15.5, dpi = 250, bg = "white")
ggsave(file.path(WD, "fig_supp5_combined.pdf"), fig, width = 14, height = 15.5)
fwrite(rbindlist(list(ta[, .(panel = "a", group, measure, n, est = trend, lo = lower95, hi = upper95, p, p_interaction = p_interaction_vs_noCNV)],
                      tb[, .(panel = "b", group, measure = "age at menarche", n, est = diff, lo, hi, p, p_interaction = NA_real_, mean_years = mean)],
                      tc[, .(panel = "c", group, measure = "height per year of menarche", n, est = slope, lo, hi, p, p_interaction)],
                      td[, .(panel = "d", group, measure = as.character(pgs), n, est, lo = est - 1.96 * se, hi = est + 1.96 * se, p, p_interaction = p_int, first_stage = fs)],
                      te[, .(panel = "e", group, measure = "MR: height per year of menarche", n, est = mr, lo = mr - 1.96 * mr_se, hi = mr + 1.96 * mr_se, p = mr_p, p_interaction = p_diff, first_stage = fs_all)]), fill = TRUE),
       file.path(WD, "fig_supp5_combined_tables.csv"))
print(tb); print(tc); print(td)
cat("fig_supp5_combined written\n")
