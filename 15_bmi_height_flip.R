# Targeted test of the paper's claim (JS 2026-10-13): in 16p11.2 deletion carriers, higher BMI goes with TALLER
# stature in childhood but SHORTER stature in adulthood (the biobank finding), i.e. the BMI-height slope reverses sign.
# The tensor-product CNV x age x BMI test is an omnibus test spread over many df; this is the 1-df contrast the
# hypothesis actually specifies.
#
# Data: one record per person (latest visit with height and BMI), ages 2-40; continuous BMI Z.
# Periods: child 2-12 y, adolescent 13-17 y, adult 18-40 y.
#  (1) within each group:  HEIGHT ~ BMI * period + sex + cohort      -> BMI slope per period; flip = adult - child slope
#  (2) carriers vs controls: HEIGHT ~ BMI * period * CNV + sex + cohort -> BMI:adult:CNV = (adult - child slope) in
#      carriers minus the same in controls (1 df)
#  (3) same with the BMI tails (>= 70th vs <= 30th centile) as a binary predictor, matching the figure strata.
# HC3 heteroskedasticity-robust SEs (sandwich) alongside OLS SEs.
# Output: bmi_height_flip.csv, fig_bmi_height_flip.{png,pdf}

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(sandwich) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
t <- readRDS(file.path(DATA, "analysis_table.rds"))
d <- t[!is.na(HEIGHT) & !is.na(BMI) & age >= 2 & age <= 40 & !is.na(sex)][order(id, -age)][, .SD[1], by = id]
d[, period := factor(fcase(age < 13, "child", age < 18, "adolescent", default = "adult"), levels = c("child", "adolescent", "adult"))]
d[, hiBMI := fifelse(BMI >= qnorm(.7), 1, fifelse(BMI <= qnorm(.3), 0, NA_real_))]
LOCI <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")

coefs <- function(m, terms) {
  V <- vcovHC(m, type = "HC3"); b <- coef(m)
  rbindlist(lapply(names(terms), function(nm) {
    L <- setNames(rep(0, length(b)), names(b)); L[names(terms[[nm]])] <- terms[[nm]]
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% vcov(m) %*% L)); se_r <- sqrt(drop(t(L) %*% V %*% L)); dfr <- m$df.residual
    data.table(term = nm, est = est, se = se, p = 2 * pt(-abs(est / se), dfr), se_robust = se_r, p_robust = 2 * pt(-abs(est / se_r), dfr))
  }))
}
res <- list()
for (pred in c("BMI", "hiBMI")) {
  for (G in c("No CNV", LOCI)) {                                   # (1) within group
    x <- d[locus == G & !is.na(get(pred))]; x[, P := get(pred)]
    m <- lm(HEIGHT ~ P * period + sex + cohort, data = droplevels(x))
    r <- coefs(m, list(`slope child` = c(P = 1), `slope adolescent` = c(P = 1, `P:periodadolescent` = 1),
                       `slope adult` = c(P = 1, `P:periodadult` = 1), `flip (adult - child)` = c(`P:periodadult` = 1)))
    res[[length(res) + 1]] <- r[, `:=`(predictor = pred, group = G, comparison = "within group",
                                       n_child = x[period == "child", .N], n_adult = x[period == "adult", .N])]
  }
  for (G in LOCI) {                                                # (2)/(3) carriers vs controls
    x <- d[locus %in% c(G, "No CNV") & !is.na(get(pred))]; x[, P := get(pred)][, CNV := as.numeric(locus != "No CNV")]
    m <- lm(HEIGHT ~ P * period * CNV + sex + cohort, data = droplevels(x))
    r <- coefs(m, list(`flip, carriers minus controls` = c(`P:periodadult:CNV` = 1),
                       `child slope, carriers minus controls` = c(`P:CNV` = 1),
                       `adult slope, carriers minus controls` = c(`P:CNV` = 1, `P:periodadult:CNV` = 1)))
    res[[length(res) + 1]] <- r[, `:=`(predictor = pred, group = G, comparison = "vs controls",
                                       n_child = x[CNV == 1 & period == "child", .N], n_adult = x[CNV == 1 & period == "adult", .N])]
  }
}
res <- rbindlist(res); setcolorder(res, c("predictor", "comparison", "group", "term"))
fwrite(res, file.path(HERE, "bmi_height_flip.csv"))
print(res[, .(predictor, comparison, group, term, est = round(est, 3), se = round(se, 3), p = signif(p, 2), p_robust = signif(p_robust, 2), n_child, n_adult)], nrows = 200)

# figure: BMI-height slope by period, per group (continuous BMI, within group)
f <- res[predictor == "BMI" & comparison == "within group" & grepl("^slope", term)]
f[, period := factor(sub("slope ", "", term), c("child", "adolescent", "adult"))][, group := factor(group, c("No CNV", LOCI), c("Controls", LOCI))]
pl <- ggplot(f, aes(period, est, colour = group, group = group)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), position = position_dodge(0.5), size = 0.6) +
  geom_line(position = position_dodge(0.5)) +
  scale_colour_manual(values = c(Controls = "grey40", `16p11.2 del` = "#C0392B", `16p11.2 dup` = "#2471A3", `22q11.2 del` = "#7D3C98", `22q11.2 dup` = "#D68910")) +
  labs(x = NULL, y = "Height Z per 1 SD BMI Z (\u00b1 95% CI)", colour = NULL,
       title = "BMI-height slope by age period (child 2-12, adolescent 13-17, adult 18-40 y)") +
  theme_classic(base_size = 15) + theme(plot.title = element_text(size = 13))
ggsave(file.path(HERE, "fig_bmi_height_flip.png"), pl, width = 9, height = 6, dpi = 200)
ggsave(file.path(HERE, "fig_bmi_height_flip.pdf"), pl, width = 9, height = 6, device = cairo_pdf)
