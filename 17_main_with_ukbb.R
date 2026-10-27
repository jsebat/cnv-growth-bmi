# Main growth figure with individual-level UK Biobank adults added to the height panels (JS 2026-10-13).
# RUNS ON EXPANSE: UKBB participant data stay in /expanse/projects/sebat1 (only the figure and result tables are copied
# back to Dropbox). Inputs in WD (/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb):
#   analysis_table.rds (clinical, from 01_build.R), supp_bars_results.csv (10_supplement.R),
#   within_vs_controls_tests.csv (14_within_vs_controls.R, HC statistics), and Molly's UKBB genomatrix.
#
# UKBB: EUR ("blue"), duplicate_drop == "no", consent retained; carriers of 16p11.2 proximal (BP4-BP5) and 22q11.2 A-D
#   (1 = del, 3 = dup); non-carriers = reference at both loci. Height and BMI Z = residual on age + age^2 within sex,
#   scaled by the non-carrier residual SD (adult reference). One record per person (ages 40-73).
# Clinical: as in 02b (latest visit with height and BMI, ages 2-40; < 18 y growth-reference Z, >= 18 y within-sample
#   Box-Cox Z). NB the two adult Z scales are both non-growth-chart, within-sample adult references.
# Height panels: clinical + UKBB points (UKBB open circles; controls show a random 3,000 UKBB non-carriers), GAM curves
#   per BMI tail (>= 70th vs <= 30th centile) over ages 2-73.
# Height statistics (clinical + UKBB):
#   within-carrier age x BMI  : y ~ s(age) + s(BMI) + ti(age, BMI), carriers only, GCV and REML
#   vs controls (CNV x age x BMI): bam(), controls = clinical SSC/SVIP + 50,000 random UKBB non-carriers, GCV and fREML
#   offset: high - low BMI tail shift of the carrier age curve (gam, REML)
# HC panels and ASD bars: unchanged (clinical only; UKBB has no head circumference).
# Output (WD): fig_main_growth_ukbb.{png,pdf}, main_with_ukbb_height_tests.csv; v2 with BMI-effect-by-age panels (JS 2026-10-23, see end)

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork) })
set.seed(1)
WD <- getwd()
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup"); GROUPS <- c("No CNV", LOCI)
COLS <- c("BMI ≥ 70th centile" = "#F0645A", "BMI ≤ 30th centile" = "#00A5AD")
K_AGE <- 10; K_BMI <- 6; K_TI <- c(6, 5); GAMMA <- 1.5   # JS 2026-10-27: smoother age x BMI surfaces (Jon: the UKBB prime panels were too wiggly at k = 15/10/c(10,10)); used by threeway_p and bmi_effect_by_age
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))
tail_of <- function(z) factor(ifelse(z >= qnorm(.7), names(COLS)[1], ifelse(z <= qnorm(.3), names(COLS)[2], NA)), levels = names(COLS))

# ---------------------------------------------------------------- data
# JS 2026-10-13: every subject scored the same way from RAW measurements (00_scoring.R): external references for height,
# HC and child BMI. Ages 2-73 for every cohort.
# JS 2026-10-27 (Jon's decision): children keep the external growth references; ADULTS are scored relative to the
# non-carrier adults of their own cohort (clinical cohorts pooled = SSC parents + SVIP non-carriers; UK Biobank = its own
# non-carriers), adjusted for sex, age and age^2 with Molly's Box-Cox scaling (score_*_cohort in 00_scoring.R). The
# pooled cross-cohort adult BMI centile scale of 2026-10-13 is replaced. Non-carrier adults therefore have mean 0 in
# every cohort by construction; the panels read carrier minus same-cohort non-carrier at each age.
source(file.path(WD, "00_scoring.R"))
cl <- readRDS(file.path(WD, "analysis_table.rds"))[!is.na(sex) & age >= 2 & age <= 73]
latest <- function(raw) cl[!is.na(cl[[raw]]) & !is.na(BMI_raw)][order(id, -age)][, .SD[1], by = id]

UKF <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"   # JS 2026-10-27: optional, so the script also runs locally (clinical only) for checks
if (file.exists(UKF)) {
  u <- fread(UKF, select = c("f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D", "Height", "BMI_raw", "age"))
  setnames(u, c("eid", "sex", "dup", "wd", "anc", "l16", "l22", "h", "bmi", "age"))
  u <- u[dup == "no" & wd == "no" & anc == "blue" & !is.na(h) & !is.na(bmi) & !is.na(age) & sex %in% c(1, 2)]
  u[, locus := fcase(l16 == 1, "16p11.2 del", l16 == 3, "16p11.2 dup", l22 == 1, "22q11.2 del", l22 == 3, "22q11.2 dup", l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
  u <- u[!is.na(locus)]
  u[, sex := ifelse(sex == 1, "male", "female")]
} else { cat("NOTE: UK Biobank genomatrix not found; running with the clinical cohorts only\n"); u <- data.table(eid = integer(), sex = character(), locus = character(), h = numeric(), bmi = numeric(), age = numeric()) }
# one scoring table: clinical (latest visit with height and BMI) + UK Biobank, scored together so each cohort's non-carrier adults are its reference
chr <- latest("HEIGHT_raw")[, .(id, source = "Clinical", locus = as.character(locus), sex = as.character(sex), age, h = HEIGHT_raw, bmi = BMI_raw)]
ukr <- u[, .(id = paste0("UKBB:", eid), source = "UK Biobank", locus, sex, age, h, bmi)]
A <- rbind(chr, ukr)[, is_ref := locus == "No CNV"]
A[, y := score_height_cohort(h, age, sex, source, is_ref)]; cat("height adult scaling:\n"); print(attr(score_height_cohort(A$h, A$age, A$sex, A$source, A$is_ref), "info"))
A[, BMI := score_bmi_cohort(bmi, age, sex, source, is_ref)];  cat("BMI adult scaling:\n");    print(attr(score_bmi_cohort(A$bmi, A$age, A$sex, A$source, A$is_ref), "info"))
A <- A[!is.na(y) & !is.na(BMI)]
ch <- A[source == "Clinical", .(id, source, locus, sex, age, y, BMI)]
uk <- A[source == "UK Biobank", .(id, source, locus, sex, age, y, BMI)]
cat("mean height Z and BMI Z (children: external reference; adults: own-cohort Box-Cox) by source x locus x age group:\n")
print(A[, .(n = .N, height_Z = round(mean(y), 2), BMI_Z = round(mean(BMI), 2)), by = .(source, locus, adult = age >= 19)][order(source, locus, adult)])
cat("UKBB:\n"); print(uk[, .N, by = locus])
ukc <- uk[locus != "No CNV"]
ukn <- uk[locus == "No CNV"]
H <- rbind(ch, ukc)                                         # carriers + clinical controls
Hctl_test <- rbind(ch[locus == "No CNV"], if (nrow(ukn) > 50000) ukn[sample(.N, 50000)] else ukn)
Hctl_plot <- rbind(ch[locus == "No CNV"], if (nrow(ukn) > 3000) ukn[sample(.N, 3000)] else ukn)

# ---------------------------------------------------------------- height statistics
within_p <- function(x, meth) tryCatch(summary(gam(y ~ s(age, k = 10) + s(BMI, k = 6) + ti(age, BMI, k = c(6, 6)), data = x, method = meth, gamma = GAMMA))$s.table["ti(age,BMI)", "p-value"], error = function(e) NA_real_)
threeway_p <- function(x, meth) {
  x <- copy(x)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]
  m <- bam(y ~ is_cnv + s(age, k = K_AGE) + s(BMI, k = K_BMI) + s(age, k = K_AGE, by = is_cnv) + ti(age, BMI, k = K_TI) +
             ti(age, BMI, k = K_TI, by = is_cnv), data = x, method = meth, gamma = GAMMA)
  p <- summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]
  if (!is.finite(p) || p <= 0 || p > 1) NA_real_ else p   # guard against numerical failure
}
offset <- function(x) {
  x <- copy(x)[, grp := tail_of(BMI)][!is.na(grp)][, grp := ordered(grp, levels = rev(names(COLS)))]
  s <- summary(gam(y ~ grp + s(age, k = 10) + s(age, k = 10, by = grp), data = x, method = "REML"))
  c(offset = unname(s$p.coeff[2]), p = unname(s$p.pv[2]))
}
# JS 2026-10-13: does the carriers' high-low BMI offset differ from the controls' offset? (1 df: CNV x BMI tail)
offset_vs_ctl <- function(xc, ctl) {
  x <- rbind(xc, ctl)[, grp := tail_of(BMI)][!is.na(grp)]
  x[, `:=`(hi = as.numeric(grp == names(COLS)[1]), is_cnv = factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV")))]
  m <- bam(y ~ is_cnv * hi + s(age, k = 10, by = is_cnv), data = x, method = "REML")
  co <- summary(m)$p.table["is_cnvCNV:hi", ]
  c(diff = unname(co[1]), p = unname(co[4]))
}
ht <- rbindlist(lapply(GROUPS, function(G) {
  x <- if (G == "No CNV") Hctl_test else H[locus == G]
  r <- data.table(group = G, n_clinical = x[source == "Clinical", .N], n_ukbb = x[source == "UK Biobank", .N],
                  within_GCV = within_p(x, "GCV.Cp"), within_REML = within_p(x, "REML"),
                  threeway_GCV = if (G == "No CNV") NA_real_ else threeway_p(rbind(Hctl_test, x), "GCV.Cp"),
                  threeway_fREML = if (G == "No CNV") NA_real_ else threeway_p(rbind(Hctl_test, x), "REML"))
  o <- offset(x); r[, `:=`(offset = o[["offset"]], offset_p = o[["p"]])]
  if (G != "No CNV") { v <- offset_vs_ctl(x, Hctl_test); r[, `:=`(offset_vs_ctl = v[["diff"]], offset_vs_ctl_p = v[["p"]])] } else r[, `:=`(offset_vs_ctl = NA_real_, offset_vs_ctl_p = NA_real_)]
  cat(G, "done\n"); r
}))
fwrite(ht, file.path(WD, "main_with_ukbb_height_tests.csv")); print(ht)

# ---------------------------------------------------------------- panels
base_theme <- theme_classic(base_size = 15) + theme(legend.text = element_text(size = 8), legend.position = "inside", legend.background = element_blank(),
                                                    legend.key.height = unit(8, "pt"), plot.title = element_text(face = "bold", hjust = 0.5))
traj <- function(x, G, ylab, title, lab, xmax, lpos) {
  x <- copy(x)[, grp := tail_of(BMI)][!is.na(grp)]
  ctl <- G == "No CNV"
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(aes(shape = source), size = if (ctl) 0.5 else 1.1, alpha = if (ctl) 0.2 else 0.55, stroke = 0.5) +
    scale_shape_manual(values = c(Clinical = 16, `UK Biobank` = 1), guide = "none") +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 8), method.args = list(method = "REML"), linewidth = 0.9, alpha = 0.2) +   # log(age): childhood resolution over 2-73 y
    annotate("text", x = xmax, y = -3.6, label = lab, hjust = 1, vjust = 0, size = 3.5) +
    scale_colour_manual(values = COLS, breaks = names(COLS)) + scale_fill_manual(values = COLS, breaks = names(COLS)) +
    coord_cartesian(xlim = c(2, xmax), ylim = c(-3.6, 3.6)) +
    labs(x = "Age (years)", y = ylab, colour = NULL, fill = NULL, title = title,
         subtitle = sprintf("%s (n = %d)", if (ctl) "Non-carrier controls" else G, nrow(x))) +
    base_theme + theme(legend.position.inside = lpos)
}
plots <- list()
for (G in GROUPS) {
  r <- ht[group == G]
  # JS 2026-10-13: REML only (GCV dropped)
  lab <- paste0(sprintf("age × BMI p = %s", pfx(r$within_REML)),
                if (G != "No CNV") sprintf("\nCNV × age × BMI vs controls p = %s", pfx(r$threeway_fREML)) else "",
                sprintf("\nhigh − low BMI offset %+.2f, p = %s", r$offset, pfx(r$offset_p)),
                if (G != "No CNV") sprintf("\noffset vs controls %+.2f, p = %s", r$offset_vs_ctl, pfx(r$offset_vs_ctl_p)) else "")
  lab <- gsub("= < ", "< ", lab, fixed = TRUE)
  plots[[paste("HEIGHT", G)]] <- traj(if (G == "No CNV") Hctl_plot else H[locus == G], G, "Height (SD)", if (G == "No CNV") "Height" else NULL, lab, 60, c(0.75, 0.93))
}
# JS 2026-10-27: HC adults (>= 18 y) scored relative to the clinical non-carrier adults (SSC parents, SVIP) with the Box-Cox
# scaling, BMI as above (clinical only: UK Biobank has no head circumference)
dh <- latest("HEADCIRCUM_raw")[, .(id, source = "Clinical", locus = as.character(locus), sex = as.character(sex), age, hc = HEADCIRCUM_raw, bmi = BMI_raw)][, is_ref := locus == "No CNV"]
dh[, y := score_hc_cohort(hc, age, sex, source, is_ref)]; cat("HC adult scaling:\n"); print(attr(score_hc_cohort(dh$hc, dh$age, dh$sex, dh$source, dh$is_ref), "info"))
dh[, BMI := score_bmi_cohort(bmi, age, sex, source, is_ref)]
dh <- dh[!is.na(y) & !is.na(BMI), .(id, source, locus, sex, age, y, BMI)]
hcw <- rbindlist(lapply(GROUPS, function(G) {
  x <- dh[locus == G]
  data.table(group = G, n = nrow(x), within_GCV = within_p(x, "GCV.Cp"), within_REML = within_p(x, "REML"),
             threeway_GCV = if (G == "No CNV") NA_real_ else threeway_p(dh[locus %in% c(G, "No CNV")], "GCV.Cp"),
             threeway_REML = if (G == "No CNV") NA_real_ else threeway_p(dh[locus %in% c(G, "No CNV")], "REML"),
             offset_vs_ctl = if (G == "No CNV") NA_real_ else offset_vs_ctl(dh[locus == G], dh[locus == "No CNV"])[["diff"]],
             offset_vs_ctl_p = if (G == "No CNV") NA_real_ else offset_vs_ctl(dh[locus == G], dh[locus == "No CNV"])[["p"]])
}))
fwrite(hcw, file.path(WD, "main_with_ukbb_hc_tests.csv")); print(hcw)
for (G in GROUPS) {
  x <- dh[locus == G]; w <- hcw[group == G]; o <- offset(x)
  # JS 2026-10-13: HC shows a constant BMI offset; age-interaction statistics dropped from the panel (kept in the CSV)
  lab <- paste0(sprintf("high − low BMI offset %+.2f, p = %s", o[["offset"]], pfx(o[["p"]])),
                if (G != "No CNV") sprintf("\noffset vs controls %+.2f, p = %s", w$offset_vs_ctl, pfx(w$offset_vs_ctl_p)) else "")
  lab <- gsub("= < ", "< ", lab, fixed = TRUE)
  plots[[paste("HC", G)]] <- traj(x, G, "Head circumference (SD)", if (G == "No CNV") "Head circumference" else NULL, lab, 60, c(0.75, 0.93))
}
sb <- fread(file.path(WD, "supp_bars_results.csv"))
asd_bar <- function(prefix, ylab) {
  r <- sb[startsWith(panel, prefix) & grepl("growth|Non-carriers", group)][order(group != "Non-carriers")]   # non-carriers first (JS)
  r[group == "Non-carriers", group := "Non-carriers\n(SSC sib-pairs)"][grepl("growth", group), group := sub("\n.*", "\nCNVs", group)][, group := factor(group, unique(group))]
  d <- rbind(r[, .(group, grp = "No ASD", m = mean_0, se = se_0, n = n_0)], r[, .(group, grp = "ASD", m = mean_1, se = se_1, n = n_1)])[, grp := factor(grp, c("No ASD", "ASD"))]
  r[, lab := sprintf("Δ %+.2f, p = %s\nadj. %+.2f\np = %s", diff, pfx(p_welch), beta_adj, pfx(p_adj))][, bold := p_adj < 0.05]
  r[within_family %in% TRUE, `:=`(lab = sprintf("sib-pair\nΔ %+.2f\np = %s", diff, pfx(p_welch)), bold = p_welch < 0.05)]
  yr <- range(c(d$m - d$se, d$m + d$se)); top <- yr[2] + 0.75 * diff(yr); bot <- yr[1] - 0.18 * diff(yr)
  ggplot(d, aes(group, m, fill = grp)) + geom_hline(yintercept = 0, colour = "grey50") + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(y = bot, label = n), position = position_dodge(0.8), size = 3.3) +
    geom_text(data = r, aes(group, top, label = lab, fontface = ifelse(bold, "bold", "plain")), inherit.aes = FALSE, size = 3.4, vjust = 1, lineheight = 0.9) +
    scale_fill_manual(values = c("#8C9A9A", "#F28E2B")) + coord_cartesian(ylim = c(bot, top)) +
    labs(x = NULL, y = ylab, fill = NULL, subtitle = "ASD vs no ASD") +
    theme_classic(base_size = 15) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
                                           legend.text = element_text(size = 9), legend.key.size = unit(10, "pt"))
}
plots[["ASD HEIGHT"]] <- asd_bar("B.", "Height Z"); plots[["ASD HC"]] <- asd_bar("C.", "Head circumference Z")
ord <- c(paste("HEIGHT", GROUPS), "ASD HEIGHT", paste("HC", GROUPS), "ASD HC")
fig <- wrap_plots(plots[ord], ncol = 2, byrow = FALSE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(WD, "fig_main_growth_ukbb.png"), fig, width = 12, height = 20.5, dpi = 200)
ggsave(file.path(WD, "fig_main_growth_ukbb.pdf"), fig, width = 12, height = 20.5)

# ---------------------------------------------------------------- life-course height effect (carriers vs controls by period)
P <- rbind(H, Hctl_test)[, period := factor(fcase(age < 13, "child", age < 18, "adolescent", default = "adult"), c("child", "adolescent", "adult"))]
lc <- rbindlist(lapply(LOCI, function(L) {
  x <- P[locus %in% c(L, "No CNV")][, CNV := as.numeric(locus != "No CNV")]
  rbindlist(lapply(c("all", "high", "low"), function(s) {
    xs <- if (s == "all") x else x[tail_of(BMI) == names(COLS)[ifelse(s == "high", 1, 2)]]
    co <- summary(lm(y ~ CNV * period + sex, data = xs))$coefficients
    data.table(locus = L, bmi = s, child = co["CNV", 1], child_p = co["CNV", 4],
               adult = co["CNV", 1] + co["CNV:periodadult", 1], change = co["CNV:periodadult", 1], change_p = co["CNV:periodadult", 4],
               n_child = xs[CNV == 1 & period == "child", .N], n_adult = xs[CNV == 1 & period == "adult", .N])
  }))
}))
fwrite(lc, file.path(WD, "main_with_ukbb_lifecourse.csv")); print(lc, digits = 3)


# ---------------------------------------------------------------- JS 2026-10-23: figure version 2 — BMI effect by age
# The prime panels show, for each row, the BMI effect on the outcome as a continuous function of age:
#   slope(age) = [E(y | BMI = +1 SD, age) − E(y | BMI = −1 SD, age)] / 2
# evaluated on the fitted age × BMI surface of the very model behind the "CNV × age × BMI vs controls" test (threeway_p:
#   y ~ is_cnv + s(age) + s(BMI) + s(age, by = is_cnv) + ti(age, BMI) + ti(age, BMI, by = is_cnv), REML), for carriers and
# for controls, with 95% CI from the coefficient covariance. That tensor term is what the three-way test evaluates, so the
# panel shows where in the age range the carriers' BMI effect departs from the controls'. The control row uses the
# within-group model (y ~ s(age) + s(BMI) + ti(age, BMI)) and its age × BMI p. A rug marks carrier ages (data support).
# Replaces the 2026-10-21 bar panels (high − low BMI-tail offsets by age group), kept in
#   17_main_with_ukbb.R.bak_bars_20261023 (Jon's scratch) and main_with_ukbb_offset_by_age.csv.
# Output: main_with_ukbb_bmi_effect_by_age.csv, fig_main_growth_ukbb_v2.{png,pdf}
GRID <- seq(2.5, 60, by = 0.5)
grid_of <- function(a) GRID[GRID >= quantile(a, 0.025) & GRID <= quantile(a, 0.975)]   # JS 2026-10-27: draw each curve only over the central 95% of that group's ages (the tails beyond the data support bent downward: smoother edge effect, Jon)
support_of <- function(a_grid, ages, w = 3, nmin = 5) sapply(a_grid, function(a) sum(abs(ages - a) <= w) >= nmin)   # JS 2026-10-27: carrier curve drawn only where >= 5 carriers lie within 3 y (removes the dip across the 18-30 y gap of the clinical carriers)
slope_rows <- function(m, grid, groups) {
  V <- vcov(m); B <- coef(m)
  rbindlist(lapply(groups, function(g) {
    nd <- function(b) { d <- data.frame(age = grid, BMI = b); if (!is.null(g)) d$is_cnv <- factor(g, levels = c("No CNV", "CNV")); d }
    X <- (predict(m, nd(1), type = "lpmatrix") - predict(m, nd(-1), type = "lpmatrix")) / 2
    data.table(who = if (is.null(g) || g == "No CNV") "Non-carriers" else "Carriers", age = grid,
               slope = drop(X %*% B), se = sqrt(rowSums((X %*% V) * X)))
  }))
}
bmi_effect_by_age <- function(xc, ctl) {
  if (is.null(xc)) {   # controls only: within-group age × BMI model
    m <- bam(y ~ s(age, k = K_AGE) + s(BMI, k = K_BMI) + ti(age, BMI, k = K_TI), data = ctl, method = "REML", gamma = GAMMA)
    r <- slope_rows(m, grid_of(ctl$age), list(NULL)); r[, p := summary(m)$s.table["ti(age,BMI)", "p-value"]]
    # JS 2026-10-23: age-constant BMI effect (continuous BMI, 1 df): y ~ BMI + s(age)
    m0 <- bam(y ~ BMI + s(age, k = 15), data = ctl, method = "REML"); c0 <- summary(m0)$p.table["BMI", ]
    r[, `:=`(est_main = unname(c0[1]), se_main = unname(c0[2]), p_main = unname(c0[4]))]; return(r)
  }
  x <- rbind(ctl, xc)[, is_cnv := factor(ifelse(locus == "No CNV", "No CNV", "CNV"), levels = c("No CNV", "CNV"))]   # controls first: same row order as threeway_p, so the p matches the trajectory panel
  m <- bam(y ~ is_cnv + s(age, k = K_AGE) + s(BMI, k = K_BMI) + s(age, k = K_AGE, by = is_cnv) + ti(age, BMI, k = K_TI) +
             ti(age, BMI, k = K_TI, by = is_cnv), data = x, method = "REML", gamma = GAMMA)
  pv <- summary(m)$s.table["ti(age,BMI):is_cnvCNV", "p-value"]; if (!is.finite(pv) || pv <= 0 || pv > 1) pv <- NA_real_   # JS 2026-10-27: same guard as threeway_p (a negative p appeared for 22q11.2 del height)
  r <- slope_rows(m, grid_of(xc$age), list("No CNV", "CNV")); r[, p := pv]; r[who == "Carriers" & !support_of(age, xc$age), `:=`(slope = NA_real_, se = NA_real_)]
  # JS 2026-10-23: age-constant CNV x BMI (continuous BMI, 1 df), AGE-ADJUSTED: the BMI slope is allowed to vary with age
  # in a way shared by both groups (varying-coefficient term s(age, by = BMI); a numeric by-smooth is not centred, so it
  # carries the BMI main effect) and cnvBMI = BMI x carrier tests a constant shift of the carriers' slope relative to the
  # controls' at every age, i.e. the mean vertical gap between the two curves. A model with one slope per group across all
  # ages (y ~ is_cnv * BMI + s(age, by = is_cnv)) is confounded by age composition (controls mostly UKBB adults, where the
  # BMI slope is negative; carriers mostly children, where it is positive) and gave p = 3e-5 for 16p11.2 del height.
  x[, cnvBMI := as.numeric(is_cnv == "CNV") * BMI]
  m0 <- bam(y ~ is_cnv + s(age, k = 15, by = is_cnv) + s(age, k = 10, by = BMI) + cnvBMI, data = x, method = "REML")
  c0 <- summary(m0)$p.table["cnvBMI", ]
  r[, `:=`(est_main = unname(c0[1]), se_main = unname(c0[2]), p_main = unname(c0[4]))]; r
}
bea <- rbindlist(lapply(GROUPS, function(G) {
  cat("BMI effect by age:", G, "\n")
  rbind(cbind(outcome = "Height", group = G, bmi_effect_by_age(if (G == "No CNV") NULL else H[locus == G], Hctl_test)),
        cbind(outcome = "Head circumference", group = G, bmi_effect_by_age(if (G == "No CNV") NULL else dh[locus == G], dh[locus == "No CNV"])))
}))
fwrite(bea, file.path(WD, "main_with_ukbb_bmi_effect_by_age.csv"))
print(bea[age %in% c(5, 10, 15, 20, 30, 50), .(outcome, group, who, age, slope = round(slope, 2), se = round(se, 2), p = signif(p, 2))], nrows = 200)
cat("age-constant CNV x BMI (continuous BMI):\n"); print(unique(bea[, .(outcome, group, est_main = round(est_main, 3), se_main = round(se_main, 3), p_main = signif(p_main, 3), p_threeway = signif(p, 3))]))
# JS 2026-10-23: control-row BMI tests, reported in a table rather than in panels A, A', G, G'
ctl_tab <- unique(bea[group == "No CNV", .(outcome, n_controls = ifelse(outcome == "Height", nrow(Hctl_test), dh[locus == "No CNV", .N]),
                                            bmi_effect_sd_per_sd = est_main, se = se_main, p_bmi = p_main, p_age_x_bmi = p)])
fwrite(ctl_tab, file.path(WD, "main_with_ukbb_control_bmi_tests.csv")); print(ctl_tab)

BEACOLS <- c("Non-carriers" = "#8C9A9A", "Carriers" = "#4B6EA9")
eff_panel <- function(r, carriers, ylab, ctl_only, xmax = 60, ylim = c(-1, 1.2)) {
  r <- copy(r)[, who := factor(who, names(BEACOLS))]
  pv <- unique(r$p)[1]; pm <- unique(r$p_main)[1]
  # two lines, top right: line 1 = age-constant (CNV x) BMI effect, line 2 = its age modification; each bold when P < 0.05
  lab1 <- gsub("= < ", "< ", if (ctl_only) sprintf("BMI p = %s", pfx(pm)) else sprintf("CNV × BMI p = %s", pfx(pm)), fixed = TRUE)
  lab2 <- if (ctl_only) "" else gsub("= < ", "< ", sprintf("CNV × age × BMI p = %s", pfx(pv)), fixed = TRUE)   # JS 2026-10-23: control row carries no p-values (its tests go to main_with_ukbb_control_bmi_tests.csv)
  bold <- function(p) ifelse(is.finite(p) && p < 0.05, "bold", "plain")
  brk <- if (ctl_only) "Non-carriers" else names(BEACOLS)
  g <- ggplot(r, aes(age, slope, colour = who, fill = who)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_ribbon(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.9)
  if (!ctl_only) g <- g + geom_rug(data = carriers, aes(x = age), inherit.aes = FALSE, sides = "b", alpha = 0.35,
                                   length = unit(0.05, "npc"), colour = BEACOLS[["Carriers"]])
  g + annotate("text", x = xmax, y = ylim[1] + 0.1 * diff(ylim), label = lab2, hjust = 1, vjust = 0, size = 3.3, fontface = bold(pv)) +   # JS 2026-10-23: the age-constant test moved to the trajectory panel
    scale_colour_manual(values = BEACOLS, breaks = brk) + scale_fill_manual(values = BEACOLS, breaks = brk) +
    coord_cartesian(xlim = c(2, xmax), ylim = ylim) +
    labs(x = "Age (years)", y = ylab, colour = NULL, fill = NULL, subtitle = "BMI effect by age (95% CI)") +
    base_theme + theme(legend.position.inside = c(0.25, 0.9), legend.text = element_text(size = 7.5), plot.subtitle = element_text(size = 11))   # key bottom-left, above the rug; p-value top-right
}
# trajectory panel with slightly more prominent points (JS 2026-10-21, for v2)
traj2 <- function(x, G, ylab, title, lab, xmax, lpos, bold = FALSE) {
  x <- copy(x)[, grp := tail_of(BMI)][!is.na(grp)]
  ctl <- G == "No CNV"
  ggplot(x, aes(age, y, colour = grp, fill = grp)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(aes(shape = source), size = if (ctl) 0.6 else 1.4, alpha = if (ctl) 0.25 else 0.7, stroke = 0.6) +
    scale_shape_manual(values = c(Clinical = 16, `UK Biobank` = 1), guide = "none") +
    geom_smooth(method = "gam", formula = y ~ s(log(x), k = 8), method.args = list(method = "REML"), linewidth = 1, alpha = 0.2) +
    annotate("text", x = xmax, y = -3.6, label = lab, hjust = 1, vjust = 0, size = 3.5, fontface = if (bold) "bold" else "plain") +
    scale_colour_manual(values = COLS, breaks = names(COLS)) + scale_fill_manual(values = COLS, breaks = names(COLS)) +
    coord_cartesian(xlim = c(2, xmax), ylim = c(-3.6, 3.6)) +
    labs(x = "Age (years)", y = ylab, colour = NULL, fill = NULL, title = title,
         subtitle = sprintf("%s (n = %d)", if (ctl) "Non-carrier controls" else G, nrow(x))) +
    base_theme + theme(legend.position.inside = lpos)
}
HCMAX <- ceiling(quantile(dh$age, 0.975))   # JS 2026-10-27: HC axes end where the HC data end (Jon)
yl_of <- function(b) { v <- b[!is.na(slope), c(slope - 1.96 * se, slope + 1.96 * se)]; r <- range(v); r + c(-0.12, 0.25) * diff(r) }   # JS 2026-10-27: prime-panel y axis from the curves of that column (Jon: the fixed -1..1.2 window wasted space); headroom for the p label
YLH <- yl_of(bea[outcome == "Height"]); YLC <- yl_of(bea[outcome == "Head circumference"])
v2 <- list(); L1 <- LETTERS[1:6]; L2 <- LETTERS[7:12]
for (i in seq_along(GROUPS)) {
  G <- GROUPS[i]; r <- ht[group == G]; ctl <- G == "No CNV"
  # JS 2026-10-23: trajectory panel carries the age-constant test (continuous BMI): CNV × BMI p for carriers, BMI p for controls;
  # the age modification (CNV × age × BMI, or age × BMI for controls) is printed in the prime panel
  pmh <- unique(bea[outcome == "Height" & group == G, p_main])[1]; pmc <- unique(bea[outcome == "Head circumference" & group == G, p_main])[1]
  labh <- if (ctl) "" else gsub("= < ", "< ", sprintf("CNV × BMI p = %s", pfx(pmh)), fixed = TRUE)   # JS 2026-10-23: no p in the control trajectory panels (A, G): a within-control BMI effect is not the carrier-vs-control test shown in B-E
  labc <- if (ctl) "" else gsub("= < ", "< ", sprintf("CNV × BMI p = %s", pfx(pmc)), fixed = TRUE)
  v2[[paste("HEIGHT", G)]] <- traj2(if (ctl) Hctl_plot else H[locus == G], G, "Height (SD)", NULL, labh, 60, c(0.75, 0.93), bold = isTRUE(pmh < 0.05)) + labs(tag = L1[i])
  v2[[paste("EFF HEIGHT", G)]] <- eff_panel(bea[outcome == "Height" & group == G], if (ctl) NULL else H[locus == G], "Height per +1 SD BMI (SD)", ctl, ylim = YLH) + labs(tag = paste0(L1[i], "'"))
  v2[[paste("HC", G)]] <- traj2(dh[locus == G], G, "Head circumference (SD)", NULL, labc, HCMAX, c(0.75, 0.93), bold = isTRUE(pmc < 0.05)) + labs(tag = L2[i])
  v2[[paste("EFF HC", G)]] <- eff_panel(bea[outcome == "Head circumference" & group == G], if (ctl) NULL else dh[locus == G], "Head circ. per +1 SD BMI (SD)", ctl, xmax = HCMAX, ylim = YLC) + labs(tag = paste0(L2[i], "'"))
}
v2[["ASD HEIGHT"]] <- plots[["ASD HEIGHT"]] + labs(tag = L1[6]); v2[["ASD HC"]] <- plots[["ASD HC"]] + labs(tag = L2[6])
# JS 2026-10-23: column headers above the panel letters, each spanning the trajectory + effect columns
hdr <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 11, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
hdr_row <- hdr("Height") | hdr("Head circumference")   # stacked above the body as its own patchwork (nested, so no design-letter mapping)
des <- "ABCD\nEFGH\nIJKL\nMNOP\nQRST\nUUVV"   # rows 1-5: traj H | effect H | traj HC | effect HC; row 6: ASD bars span two columns
body <- wrap_plots(v2, design = des, widths = c(2, 1.35, 2, 1.35)) & theme(plot.tag = element_text(size = 22, face = "bold"))
fig2 <- hdr_row / body + plot_layout(heights = c(0.05, 1))
ggsave(file.path(WD, "fig_main_growth_ukbb_v2.png"), fig2, width = 20, height = 20.5, dpi = 200)
tryCatch(ggsave(file.path(WD, "fig_main_growth_ukbb_v2.pdf"), fig2, width = 20, height = 20.5, device = cairo_pdf),
         error = function(e) ggsave(file.path(WD, "fig_main_growth_ukbb_v2.pdf"), fig2, width = 20, height = 20.5))
cat("figure v2 written\n")
