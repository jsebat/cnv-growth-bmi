# 91_cohort_adjustment_check_v2.R (JS / Claude 2026-11-21; BMI modifier paper v2): pooled longitudinal BMI slope under three adjustments: genotype only (current primary),
# additive cohort + genotype (primary until 2026-11-20), cohort within genotype (one intercept per cohort x genotype cell).
# Input: DATA/longitudinal_v2_panels.rds (80_longitudinal_v2.R). Output: CSV given as argument 1 (run_v2_2026-11-20/longitudinal/cohort_adjustment_check.csv). Model: dZ ~ Z0 + dt + sex + ns(age0, 3) [+ adjustment] + B0 (or win + win:B0) + (1 | id).
# null model: none; Wald (Satterthwaite) test of the BMI coefficients.
suppressMessages({library(data.table); library(lmerTest); library(splines)})
P <- readRDS("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/longitudinal_v2_panels.rds")
ADJ <- c("genotype only" = "+ locus", "cohort + genotype (additive)" = "+ cohort + locus", "cohort within genotype" = "+ cg")
R <- rbindlist(lapply(c("Height", "Head circumference"), function(o) { I <- copy(if (o == "Height") P$IH else P$IC)[, cg := paste(cohort, locus)]
  rbindlist(lapply(names(ADJ), function(a) { f0 <- paste("dZ ~ Z0 + dt + sex + ns(age0, 3)", ADJ[[a]])
    m <- suppressMessages(suppressWarnings(lmer(as.formula(paste(f0, "+ B0 + (1 | id)")), I))); co <- coef(summary(m))
    mw <- suppressMessages(suppressWarnings(lmer(as.formula(paste(f0, "+ win + win:B0 + (1 | id)")), I))); cw <- coef(summary(mw)); k <- grep(":B0$", rownames(cw))
    rbind(data.table(outcome = o, adjustment = a, scope = "All ages", beta = co["B0", 1], se = co["B0", 2], p = co["B0", ncol(co)]),
          data.table(outcome = o, adjustment = a, scope = sub("^win", "", sub(":B0$", "", rownames(cw)[k])), beta = cw[k, 1], se = cw[k, 2], p = cw[k, ncol(cw)])) })) }))
R[, `:=`(beta = round(beta, 4), se = round(se, 4), p = signif(p, 3))]; options(width = 200); print(dcast(R, outcome + scope ~ adjustment, value.var = c("beta", "p")), nrows = 50)
cat("cells (cohort x genotype) with pairs, height:\n"); print(P$IH[, .N, by = .(locus, cohort)][order(locus, cohort)])
a <- commandArgs(TRUE); if (length(a)) fwrite(R, a[1])
