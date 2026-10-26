# Age-varying BMI effects on growth (carriers < 18 y, cleaned visits).
# The linear models in 03 force one BMI coefficient across childhood. If high BMI speeds growth before puberty and
# the advantage is lost afterwards (earlier puberty, earlier epiphyseal closure), a single coefficient averages a
# positive and a negative phase to ~0.  Here the BMI coefficient is a smooth function of age (mgcv varying-coefficient
# terms, s(age, by = BMI)), with a person random intercept:
#   VW within-person : Z ~ sex + locus + cohort + s(age, by = sex) + BMIb + s(age, by = BMIw) + s(id, re)
#                      beta(age) = height Z difference per 1 SD of BMI above the child's own mean, at that age
#   VP prospective   : dZ ~ Z0 + dt + sex + locus + cohort + s(age0, by = sex) + s(age0, by = B0) + s(id, re)
#                      beta(age0) = change in height Z to the next visit per 1 SD of BMI at a visit at age0
# Fitted pooled and per locus (one by-smooth per locus).  Also a fixed-window check (< 8, 8-13, >= 13 y) with lmer.
# Output: age_varying_beta.csv, age_window_coefficients.csv, fig_age_varying.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(lmerTest); library(splines); library(ggplot2) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
source(file.path(HERE, "00_helpers.R"))
GRID <- seq(3, 17, by = 0.25)

beta_curve <- function(m, newbase, var, pattern, agevar) {
  # beta(age) = d(prediction)/d(var), restricted to the by-smooth columns matching `pattern`
  nd1 <- copy(newbase); nd0 <- copy(newbase); nd1[[var]] <- 1; nd0[[var]] <- 0
  X <- predict(m, nd1, type = "lpmatrix", exclude = "s(id)") - predict(m, nd0, type = "lpmatrix", exclude = "s(id)")
  cols <- grep(pattern, colnames(X), fixed = TRUE); X[, -cols] <- 0
  b <- drop(X %*% coef(m)); se <- sqrt(rowSums((X %*% vcov(m)) * X))
  data.table(age = newbase[[agevar]], beta = b, se = se)
}

base_nd <- function(d, agevar, extra) {
  nd <- data.table(GRID); setnames(nd, agevar)
  nd[, `:=`(sex = factor("male", levels = levels(d$sex)), locus = factor(levels(d$locus)[1], levels = levels(d$locus)),
            cohort = d$cohort[1], id = d$id[1])]
  for (k in names(extra)) nd[[k]] <- extra[[k]]
  nd
}

curves <- list(); smry <- list()
for (oc in c("HEIGHT", "HEADCIRCUM")) {
  d <- long_data(oc); d[, id := factor(id)]; d[, cohort := factor(cohort)]
  iv <- intervals(d); iv[, id := factor(id)]; iv[, cohort := factor(cohort)]; iv[, locus := factor(locus)]
  mc <- if (nlevels(d$cohort) > 1) "+ cohort" else ""
  ml <- if (nlevels(d$locus) > 1) "+ locus" else ""
  for (L in c("Pooled", LOCI)) {
    dd <- if (L == "Pooled") d else d[locus == L]; ii <- if (L == "Pooled") iv else iv[locus == L]
    if (uniqueN(dd$id) < 20) next
    lterm <- if (L == "Pooled") ml else ""
    dd[, cohort := droplevels(cohort)]; ii[, cohort := droplevels(cohort)]
    cterm <- if (nlevels(dd$cohort) > 1) "+ cohort" else ""
    kk <- if (uniqueN(dd$id) < 60) 4 else 6
    fW <- as.formula(sprintf("Z ~ sex %s %s + s(age, by = sex, k = 8) + BMIb + s(age, by = BMIw, k = %d) + s(id, bs = 're')", lterm, cterm, kk))
    mW <- gam(fW, data = dd, method = "REML")
    cW <- beta_curve(mW, base_nd(dd, "age", list(BMIb = 0, BMIw = 0)), "BMIw", "s(age):BMIw", "age")
    cW <- cW[age >= dd[, quantile(age, .02)] & age <= dd[, quantile(age, .98)]]
    pW <- summary(mW)$s.table[grep("BMIw", rownames(summary(mW)$s.table)), "p-value"]
    curves[[length(curves) + 1]] <- cW[, `:=`(outcome = oc, locus = L, model = "Within-person: BMI above own mean \u2192 Z", p_smooth = pW)]
    if (uniqueN(ii$id) >= 20) {
      fP <- as.formula(sprintf("dZ ~ Z0 + dt + sex %s %s + s(age0, by = sex, k = 6) + s(age0, by = B0, k = %d) + s(id, bs = 're')", lterm, cterm, kk))
      mP <- gam(fP, data = ii, method = "REML")
      cP <- beta_curve(mP, base_nd(ii, "age0", list(Z0 = 0, dt = 1, B0 = 0)), "B0", "s(age0):B0", "age0")
      cP <- cP[age >= ii[, quantile(age0, .02)] & age <= ii[, quantile(age0, .98)]]
      pP <- summary(mP)$s.table[grep("B0", rownames(summary(mP)$s.table)), "p-value"]
      curves[[length(curves) + 1]] <- cP[, `:=`(outcome = oc, locus = L, model = "Prospective: BMI at visit \u2192 change in Z by next visit", p_smooth = pP)]
    }
  }
  # fixed age windows (lmer), pooled and per locus
  d[, win := cut(age, c(0, 8, 13, 18), right = FALSE, labels = c("<8 y", "8-13 y", "13-18 y"))]
  iv[, win := cut(age0, c(0, 8, 13, 18), right = FALSE, labels = c("<8 y", "8-13 y", "13-18 y"))]
  for (L in c("Pooled", LOCI)) {
    dd <- if (L == "Pooled") d else d[locus == L]; ii <- if (L == "Pooled") iv else iv[locus == L]
    if (uniqueN(dd$id) < 20) next
    adj <- paste(c(if (L == "Pooled" && uniqueN(dd$locus) > 1) "locus", if (uniqueN(dd$cohort) > 1) "cohort"), collapse = " + ")
    adj <- if (adj == "") "" else paste("+", adj)
    for (spec in list(list(m = "W", f = sprintf("Z ~ BMIw:win + BMIb + win + sex * ns(age, 3) %s + (1 | id)", adj), dat = dd, t = "BMIw:win"),
                      list(m = "P", f = sprintf("dZ ~ B0:win + Z0 + dt + win + sex * ns(age0, 3) %s + (1 | id)", adj), dat = ii, t = "win.*:B0|B0:win"))) {
      m <- tryCatch(suppressMessages(lmer(as.formula(spec$f), data = spec$dat)), error = function(e) NULL); if (is.null(m)) next
      co <- summary(m)$coefficients; r <- grep("BMIw:win|win.*:B0|B0:win", rownames(co))
      n <- spec$dat[, .(persons = uniqueN(id), rows = .N), by = win]
      out <- data.table(outcome = oc, locus = L, model = spec$m, win = sub(".*win", "", sub(":B0", "", sub("BMIw:", "", rownames(co)[r]))),
                        est = co[r, 1], se = co[r, 2], p = co[r, 5])
      smry[[length(smry) + 1]] <- merge(out, n[, win := as.character(win)], by = "win", all.x = TRUE)
    }
  }
}
curves <- rbindlist(curves); smry <- rbindlist(smry)
fwrite(curves, file.path(HERE, "age_varying_beta.csv")); fwrite(smry, file.path(HERE, "age_window_coefficients.csv"))
print(unique(curves[, .(outcome, locus, model, p_smooth = signif(p_smooth, 2))]))
print(smry[, .(outcome, locus, model, win, est = round(est, 3), se = round(se, 3), p = signif(p, 2), persons, rows)][order(outcome, model, locus, win)])

cv <- curves[outcome == "HEIGHT"][, locus := factor(locus, levels = c("Pooled", LOCI))]
MODEL_LEVELS <- c("Prospective: BMI at visit → change in Z by next visit", "Within-person: BMI above own mean → Z")
cv[, model := factor(model, levels = MODEL_LEVELS)]
lab <- unique(cv[, .(locus, model, p_smooth)])[, txt := sprintf("age-varying BMI effect p = %s", format.pval(p_smooth, digits = 2, eps = 1e-4))]
fig <- ggplot(cv, aes(age, beta)) +
  annotate("rect", xmin = 8, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), fill = "#F0645A", alpha = 0.25) +
  geom_line(colour = "#C0392B", linewidth = 0.9) +
  geom_text(data = lab, aes(x = 3, y = Inf, label = txt), hjust = 0, vjust = 1.4, size = 2.8) +
  facet_grid(model ~ locus, scales = "free_y", labeller = label_wrap_gen(28)) +
  labs(x = "Age (years)", y = "Height Z per 1 SD of BMI Z (\u00b1 95% CI)",
       title = "Age-varying effect of BMI on height in CNV carriers < 18 y (cleaned longitudinal visits; SVIP, Searchlight, Cardiff)",
       subtitle = "Top row (prospective): change in height Z to the next visit per 1 SD of BMI Z at a visit, by age at that visit.
Bottom row (within-person): difference in a child's height Z at visits where their BMI Z is 1 SD above their own average, by age.\n
Grey band = 8-13 y. Person random intercept; adjusted for sex-specific age curves, locus and cohort.") +
  theme_bw(base_size = 10)
ggsave(file.path(HERE, "fig_age_varying.png"), fig, width = 14, height = 7, dpi = 200)
ggsave(file.path(HERE, "fig_age_varying.pdf"), fig, width = 14, height = 7, device = cairo_pdf)
