# UK Biobank: effect of body mass on adult height, comparing body mass at age 10 (self-reported comparative body size,
# field 1687: Thinner / About average / Plumper) with adult BMI group (Low / Middle / High) IN THE SAME INDIVIDUALS,
# by CNV genotype (16p11.2 BP4-BP5 del/dup, 22q11.2 A-D dup; 22q11.2 del has too few EUR carriers).
# JS 2026-10-21, after Molly's comparitiveBodySize/03 and /05 scripts (msacks/UKBB/comparitiveBodySize).
#
# Scaling (Molly's Box-Cox scaling, Methods of the Nature submission): for height and BMI separately, within sex,
#   fit lm(y ~ age + age^2 + PC1-10), choose the Box-Cox lambda maximising the profile likelihood, transform y,
#   take the residuals of the same model on the transformed scale and divide by their SD -> Z (SD units), then
#   recombine the sexes. (Molly's UKBB plots residualise on sex + age + age^2 + PC1-10; the PCs are kept here.)
# Groups, each coded as an ORDERED factor -1 / 0 / +1:
#   age 10 : Thinner (-1) < About average (0) < Plumper (+1)        [first valid answer across instances 0-2]
#   adult  : BMI Z rank percentile <= 30 (-1), 30-70 (0), >= 70 (+1)  [as in Molly's 05 script, on the Box-Cox Z]
# Effect = linear trend of height Z per step (lm(hZ ~ step) within each CNV group), 95% CI and p; the same people
#   enter both estimates (sample = blue/EUR cluster, QC'd, with height, BMI, age, PCs and a valid field-1687 answer).
#   Also: the two measures mutually adjusted (lm(hZ ~ step10 + stepAdult)), and trend x CNV interaction vs No CNV.
# Sample: blue ancestry cluster (EUR), duplicate_drop == no, consent retained.
# Output (WD): ukbb_bodysize10_vs_adult_trends.csv, ukbb_bodysize10_vs_adult_counts.csv, fig_ukbb_bodysize10_vs_adult.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(MASS) })
WD <- getwd()
GM <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"
BS <- "/expanse/projects/sebat1/msacks/UKBB/comparitiveBodySize/field1687_comparative_body_size_age10.tsv"

gm <- fread(GM, select = c("f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", paste0("PC", 1:10),
                           "Height", "BMI_raw", "age", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"))
setnames(gm, c("f.eid", "sex_combined", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"), c("eid", "sex", "l16", "l22"))
gm <- gm[duplicate_drop == "no" & withdrawn_consent == "no" & ancestry == "blue" & sex %in% 1:2 &
         !is.na(Height) & !is.na(BMI_raw) & !is.na(age) & l16 %in% 1:3 & l22 %in% 1:3]
gm <- gm[complete.cases(gm[, paste0("PC", 1:10), with = FALSE])]
gm[, locus := fcase(l16 == 1, "16p11.2 DEL", l16 == 3, "16p11.2 DUP", l22 == 1, "22q11.2 DEL", l22 == 3, "22q11.2 DUP",
                    l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
gm <- gm[!is.na(locus)]
cat("UKBB blue, QC, height+BMI+age+PCs:", nrow(gm), "\n"); print(gm[, .N, by = locus])

# field 1687 (coding 100426: 1 Thinner, 2 Plumper, 3 About average; -1/-3 missing), first valid answer across instances
bs <- fread(BS); setnames(bs, "f.eid", "eid")
inst <- grep("^f\\.1687\\.", names(bs), value = TRUE)
bs[, code := NA_integer_]
for (v in inst) bs[is.na(code) & get(v) %in% 1:3, code := as.integer(get(v))]
bs <- bs[!is.na(code), .(eid, body_size_10 = factor(code, levels = c(1, 3, 2), labels = c("Thinner", "About average", "Plumper")))]
d <- merge(gm, bs, by = "eid")                       # SAME individuals for both measures from here on
cat("with a valid body size at age 10:", nrow(d), "\n")

# ---- Molly's Box-Cox scaling, within sex, covariates age + age^2 + PC1-10 ----
d[, age_sq := age^2]
covs <- paste("age + age_sq +", paste0("PC", 1:10, collapse = " + "))
bcz <- function(y, dd) {
  z <- rep(NA_real_, nrow(dd)); lam <- c()
  for (s in 1:2) {
    i <- which(dd$sex == s); x <- copy(dd[i])[, y := y[i]]
    f0 <- as.formula(paste("y ~", covs))
    b <- boxcox(f0, data = x, lambda = seq(-2, 2, 0.05), plotit = FALSE); l <- b$x[which.max(b$y)]; lam[s] <- l   # formula + data: boxcox.lm cannot see a local data frame
    x[, yb := if (abs(l) < 1e-8) log(y) else (y^l - 1) / l]
    r <- resid(lm(as.formula(paste("yb ~", covs)), data = x))
    z[i] <- r / sd(r)
  }
  cat("  lambda female/male:", round(lam[2], 2), "/", round(lam[1], 2), "\n"); z
}
cat("height:"); d[, hZ := bcz(Height, d)]
cat("BMI:");    d[, bZ := bcz(BMI_raw, d)]
cat("Shapiro-Wilk W of 5000 random hZ / bZ:", round(shapiro.test(sample(d$hZ, 5000))$statistic, 3), round(shapiro.test(sample(d$bZ, 5000))$statistic, 3), "\n")

# ---- ordered factors coded -1 / 0 / +1 ----
d[, bmi_pct := frank(bZ, ties.method = "average") / .N]
d[, bmi_group := factor(fifelse(bmi_pct <= 0.30, "Low", fifelse(bmi_pct >= 0.70, "High", "Middle")), levels = c("Low", "Middle", "High"))]
d[, step10 := c(Thinner = -1, "About average" = 0, Plumper = 1)[as.character(body_size_10)]]
d[, stepAd := c(Low = -1, Middle = 0, High = 1)[as.character(bmi_group)]]
cat("age-10 body size x adult BMI group (all):\n"); print(table(d$body_size_10, d$bmi_group))
cat("Spearman(step10, stepAd) =", round(cor(d$step10, d$stepAd, method = "spearman"), 3), "\n")
GROUPS <- c("16p11.2 DEL", "No CNV", "16p11.2 DUP", "22q11.2 DUP")
cnt <- rbind(d[, .(measure = "Body size at age 10", level = as.character(body_size_10)), by = .(group = locus)][, .N, by = .(group, measure, level)],
             d[, .(measure = "Adult BMI group", level = as.character(bmi_group)), by = .(group = locus)][, .N, by = .(group, measure, level)])
fwrite(cnt[order(group, measure, level)], file.path(WD, "ukbb_bodysize10_vs_adult_counts.csv")); print(cnt[group != "No CNV"])
cat("22q11.2 DEL carriers (not plotted):", d[locus == "22q11.2 DEL", .N], "\n")

# ---- trends: height Z per step, within each CNV group; marginal and mutually adjusted; interaction vs No CNV ----
MEAS <- c(step10 = "Body size at age 10", stepAd = "Adult BMI group")
tr <- rbindlist(lapply(GROUPS, function(G) {
  x <- d[locus == G]
  rbindlist(lapply(names(MEAS), function(v) {
    m <- lm(as.formula(paste("hZ ~", v)), data = x); co <- summary(m)$coefficients; ci <- confint(m)
    mj <- lm(hZ ~ step10 + stepAd, data = x); cj <- summary(mj)$coefficients; cij <- confint(mj)
    if (G == "No CNV") ip <- NA_real_ else {
      xi <- d[locus %in% c(G, "No CNV")][, CNV := as.numeric(locus == G)]
      ip <- summary(lm(as.formula(paste("hZ ~", v, "* CNV")), data = xi))$coefficients[paste0(v, ":CNV"), 4]
    }
    data.table(group = G, measure = MEAS[[v]], n = nrow(x),
               trend = co[v, 1], se = co[v, 2], lower95 = ci[v, 1], upper95 = ci[v, 2], p = co[v, 4],
               trend_adj = cj[v, 1], se_adj = cj[v, 2], lower95_adj = cij[v, 1], upper95_adj = cij[v, 2], p_adj = cj[v, 4],
               p_interaction_vs_noCNV = ip)
  }))
}))
fwrite(tr, file.path(WD, "ukbb_bodysize10_vs_adult_trends.csv")); print(tr, digits = 3)

# ---- plot: the ordered-factor effect (SD of height per step) with 95% CI, by CNV group, both measures ----
pfx <- function(p) ifelse(p < 1e-300, "p < 1e-300", ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.2g", p)))
PLOT_GROUPS <- c("16p11.2 DEL", "No CNV", "16p11.2 DUP")   # JS 2026-10-21: 22q11.2 dup kept in the CSV, not plotted
tr_all <- copy(tr); tr <- tr[group %in% PLOT_GROUPS][, group := factor(group, PLOT_GROUPS)][, measure := factor(measure, MEAS)]
tr[, lab := pfx(p)][, bold := p < 0.05]
COLS <- c("Body size at age 10" = "#F28E2B", "Adult BMI group" = "#4B6EA9")
pd <- position_dodge(width = 0.6)
yr <- range(c(tr$lower95, tr$upper95, 0)); top <- yr[2] + 0.12 * diff(yr); bot <- yr[1] - 0.10 * diff(yr)
nlab <- tr[measure == MEAS[[1]], .(group, n)]
p <- ggplot(tr, aes(group, trend, colour = measure)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), position = pd, width = 0.2, linewidth = 0.8) +
  geom_point(position = pd, size = 4) +
  geom_text(aes(y = ifelse(measure == MEAS[[1]], top, top - 0.07 * diff(yr)), label = lab, fontface = ifelse(bold, "bold", "plain"), group = measure), position = pd, size = 3.8, show.legend = FALSE) +   # group = measure: fontface must not define the dodge groups
  geom_text(data = nlab, aes(group, bot, label = paste0("n = ", ifelse(n >= 10000, paste0(round(n / 1000), "k"), n))), inherit.aes = FALSE, size = 4, colour = "grey30") +
  scale_colour_manual(values = COLS, name = NULL) +
  coord_cartesian(ylim = c(bot, top + 0.04 * diff(yr))) +
  labs(title = "Body mass and adult height: body size at age 10 vs adult BMI
in the same individuals (UKBB EUR)",
       subtitle = "Linear trend per ordered step (Thinner < Average < Plumper; Low < Middle < High BMI)
height Box-Cox Z per step, 95% CI",
       x = NULL, y = "Adult height (SD) per step") +
  theme_classic(base_size = 16) + theme(legend.position = "top", legend.justification = "left", plot.title = element_text(face = "bold", size = 14),
                                        plot.subtitle = element_text(size = 11, colour = "grey30"))
ggsave(file.path(WD, "fig_ukbb_bodysize10_vs_adult.png"), p, width = 10, height = 6.5, dpi = 200)
tryCatch(ggsave(file.path(WD, "fig_ukbb_bodysize10_vs_adult.pdf"), p, width = 10, height = 6.5, device = cairo_pdf),
         error = function(e) ggsave(file.path(WD, "fig_ukbb_bodysize10_vs_adult.pdf"), p, width = 10, height = 6.5))
cat("done\n")
