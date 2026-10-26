# Life-course effect of the 16p11.2 BP4-BP5 deletion on height: clinical cohorts by age period (carriers vs non-carrier
# controls), overall and by BMI stratum, alongside Molly's adult biobank estimates (Sacks et al. Supplementary Tables:
# T2 per cohort, T3 meta-analysis; "Sacks Supplementary Tables.xlsx", NatureSubmission_04_13_2026).
#
# Clinical: one record per person (latest visit with height and BMI), ages 2-40; HEIGHT ~ CNV * period + sex, fitted
#   (a) all BMI, (b) BMI >= 70th centile only, (c) BMI <= 30th centile only; carriers and controls in the same stratum.
#   CNV effect per period = carriers minus controls (height Z); change = adolescent/adult minus child.
# Biobank: betaHeight (SD units, adults), per cohort/ancestry with >= 10 carriers, and the meta-analysis.
# Output: lifecourse_16pdel_height.csv, fig_lifecourse_16pdel_height.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(readxl) })
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
HERE <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
t <- readRDS(file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/analysis_table.rds"))
d <- t[!is.na(HEIGHT) & !is.na(BMI) & age >= 2 & age <= 40 & !is.na(sex) & locus %in% c("No CNV", "16p11.2 del")][order(id, -age)][, .SD[1], by = id]
d[, period := factor(fcase(age < 13, "Child 2-12 y", age < 18, "Adolescent 13-17 y", default = "Adult 18-40 y"),
                     levels = c("Child 2-12 y", "Adolescent 13-17 y", "Adult 18-40 y"))]
d[, CNV := as.numeric(locus != "No CNV")]
d[, stratum := fifelse(BMI >= qnorm(.7), "BMI \u2265 70th centile", fifelse(BMI <= qnorm(.3), "BMI \u2264 30th centile", "middle"))]

clin <- rbindlist(lapply(list("All BMI" = d, "BMI \u2265 70th centile" = d[stratum == "BMI \u2265 70th centile"],
                              "BMI \u2264 30th centile" = d[stratum == "BMI \u2264 30th centile"]), function(x) {
  m <- lm(HEIGHT ~ CNV * period + sex, data = x); b <- coef(m); V <- vcov(m); lv <- levels(x$period)
  rbindlist(lapply(seq_along(lv), function(i) {
    L <- setNames(rep(0, length(b)), names(b)); L["CNV"] <- 1; if (i > 1) L[paste0("CNV:period", lv[i])] <- 1
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L))
    chg <- if (i > 1) c(b[paste0("CNV:period", lv[i])], summary(m)$coefficients[paste0("CNV:period", lv[i]), 4]) else c(NA, NA)
    data.table(period = lv[i], est = est, se = se, p = 2 * pt(-abs(est / se), m$df.residual),
               n_carriers = x[CNV == 1 & period == lv[i], .N], change_vs_child = chg[1], p_change = chg[2])
  }))
}), idcol = "stratum")
clin[, source := "Clinical cohorts (carriers vs controls)"]

xl <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/NatureSubmission_04_13_2026/Sacks Supplementary Tables.xlsx")
sh <- excel_sheets(xl)
num <- c("count", "betaHeight", "seHeight", "pHeight")
rd <- function(s) { x <- as.data.table(suppressMessages(read_excel(xl, sheet = sh[startsWith(sh, s)]))); for (cc in num) x[, (cc) := suppressWarnings(as.numeric(get(cc)))]; x }
t2 <- rd("T2")[locus == "16p11.2_BP4_BP5_DEL" & count >= 10 & !is.na(betaHeight)]
t3 <- rd("T3")[locus == "16p11.2_BP4_BP5_DEL"]
bio <- rbind(t2[, .(stratum = cohort, est = betaHeight, se = seHeight, p = pHeight, n_carriers = count)],
             t3[, .(stratum = "Meta-analysis", est = betaHeight, se = seHeight, p = pHeight, n_carriers = count)])
bio[, `:=`(period = "Adult biobanks", source = "Biobanks, adults (Sacks et al.)")]
out <- rbind(clin, bio, fill = TRUE)
fwrite(out, file.path(HERE, "lifecourse_16pdel_height.csv"))
print(out[, .(source, stratum, period, est = round(est, 2), se = round(se, 2), p = signif(p, 2), n_carriers, change_vs_child = round(change_vs_child, 2), p_change = signif(p_change, 2))])

# figure
cl <- clin[, .(x = period, grp = stratum, est, se, n = n_carriers, p)]
cl[, x := factor(x, c(levels(d$period), "Adult biobanks"))]
bb <- bio[, .(x = factor("Adult biobanks", levels(cl$x)), grp = stratum, est, se, n = n_carriers, p)]
cols <- c("All BMI" = "black", "BMI \u2265 70th centile" = "#F0645A", "BMI \u2264 30th centile" = "#00A5AD")
pd <- position_dodge(0.45)
chg <- clin[!is.na(p_change)]
lab_chg <- paste(chg[, sprintf("%s, %s vs child: %+.2f, p = %s", stratum, sub(" .*", "", tolower(period)), change_vs_child, signif(p_change, 2))], collapse = "
")
bb[, `:=`(meta = grp == "Meta-analysis", ord = rank(-se, ties.method = "first"))]
fig <- ggplot() + geom_hline(yintercept = 0, linetype = 2) +
  geom_line(data = cl, aes(x, est, colour = grp, group = grp), position = pd, linewidth = 0.6) +
  geom_pointrange(data = cl, aes(x, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se, colour = grp), position = pd, size = 0.6) +
  geom_pointrange(data = bb, aes(as.numeric(x) - 0.3 + 0.6 * (ord - 1) / (nrow(bb) - 1), est, ymin = est - 1.96 * se, ymax = est + 1.96 * se,
                                 size = meta, colour = I(ifelse(meta, "black", "grey55")))) +
  scale_size_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
  geom_text(data = bb[meta == TRUE], aes(as.numeric(x) + 0.32, est, label = sprintf("meta-analysis
%.2f, n = %d
p = %.0e", est, n, p)), hjust = 0, size = 3.5) +
  scale_colour_manual(values = cols) +
  annotate("text", x = 0.55, y = 0.98, label = lab_chg, hjust = 0, vjust = 1, size = 3.3) +
  coord_cartesian(ylim = c(-1.6, 1.0), xlim = c(0.6, 4.9)) +
  labs(x = NULL, y = "16p11.2 deletion effect on height\n(carriers \u2212 non-carriers, SD)", colour = NULL,
       title = "16p11.2 BP4-BP5 deletion: from average height in childhood to short stature in adulthood",
       subtitle = "Clinical cohorts (G2MH, Cardiff, SVIP, SSC, Searchlight) by age period, overall and within BMI strata; adult biobanks from Sacks et al. (grey = cohort/ancestry with ≥ 10 carriers; black = meta-analysis)") +
  theme_classic(base_size = 14) + theme(legend.position = "bottom", plot.title = element_text(size = 13), plot.subtitle = element_text(size = 9))
ggsave(file.path(HERE, "fig_lifecourse_16pdel_height.png"), fig, width = 11, height = 6.5, dpi = 200)
ggsave(file.path(HERE, "fig_lifecourse_16pdel_height.pdf"), fig, width = 11, height = 6.5, device = cairo_pdf)
