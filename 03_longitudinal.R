# Longitudinal: does a rise in BMI go with (or precede) faster growth in CNV carriers < 18 y?
#
# Carriers with >= 2 visits < 18 y that have both the outcome Z and BMI Z (SVIP, Searchlight, Cardiff; G2MH and
# SSC are single-visit). Duplicate ages within a person are collapsed to the first record.
#
# Models (lmerTest, person random intercept; covariates sex x natural spline of age, locus, cohort):
#  W  within-person : Z_it ~ BMIw_it + BMIb_i + ...        BMIw = BMI Z minus the person's own mean BMI Z
#                     "when a child's BMI is above their own usual level, is their height Z higher?"
#  P  prospective   : dZ ~ BMI_t0 + Z_t0 + dt + age_t0 ...  consecutive visits, dt >= 0.5 y
#                     "does high BMI at one visit precede a gain in height Z by the next?"
#  R  reverse       : dBMI ~ Z_t0 + BMI_t0 + dt + age_t0 ...
#                     "does being tall at one visit precede a gain in BMI Z?"
#  C  concurrent    : dZ ~ dBMI + Z_t0 + BMI_t0 + dt + age_t0 ...
#                     "do intervals with a larger BMI gain show a larger height gain?"
# Each fitted pooled (one coefficient) and by locus (coefficient x locus).
# Datasets: all; exact-age cohorts only (SVIP + Searchlight; Cardiff ages are whole years); Cardiff only;
#           all excluding Cardiff visits on stimulants or antipsychotics.
# Output: longitudinal_coefficients.csv, fig_longitudinal.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(lmerTest); library(splines); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
source(file.path(HERE, "00_helpers.R"))
OUT <- c(HEIGHT = "Height", HEADCIRCUM = "Head circumference")

fit <- function(f, data, term, by_locus) {
  mult <- length(unique(data$cohort)) > 1
  if (!mult) f <- update(f, . ~ . - cohort)
  if (length(unique(data$locus)) < 2) f <- update(f, . ~ . - locus)
  if (by_locus) f <- update(f, as.formula(sprintf(". ~ . - %s + %s:locus", term, term)))
  m <- tryCatch(suppressMessages(lmer(f, data = data)), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients
  rows <- grep(sprintf("^%s$|^%s:locus|^locus.*:%s$", term, term, term), rownames(co))
  data.table(term_row = rownames(co)[rows], est = co[rows, 1], se = co[rows, 2], p = co[rows, ncol(co)])
}

specs <- list(
  W = list(f = Z ~ BMIw + BMIb + sex * ns(age, 3) + locus + cohort + (1 | id), term = "BMIw", unit = "visits",
           label = "Within-person: BMI above own mean \u2192 outcome Z"),
  P = list(f = dZ ~ B0 + Z0 + dt + sex * ns(age0, 3) + locus + cohort + (1 | id), term = "B0", unit = "intervals",
           label = "Prospective: BMI Z at visit \u2192 change in outcome Z"),
  R = list(f = dB ~ Z0 + B0 + dt + sex * ns(age0, 3) + locus + cohort + (1 | id), term = "Z0", unit = "intervals",
           label = "Reverse: outcome Z at visit \u2192 change in BMI Z"),
  C = list(f = dZ ~ dB + Z0 + B0 + dt + sex * ns(age0, 3) + locus + cohort + (1 | id), term = "dB", unit = "intervals",
           label = "Concurrent: change in BMI Z \u2192 change in outcome Z"))

res <- list()
for (oc in names(OUT)) for (ss in c("all", "uncleaned", "exact_age", "Cardiff", "no_meds")) {
  d <- long_data(oc, ss); if (nrow(d) == 0 || uniqueN(d$id) < 15) next
  iv <- intervals(d)
  for (k in names(specs)) {
    sp <- specs[[k]]; dat <- if (sp$unit == "visits") d else iv
    if (sp$unit == "intervals" && uniqueN(dat$id) < 15) next
    for (bl in c(FALSE, TRUE)) {
      r <- fit(sp$f, dat, sp$term, bl); if (is.null(r)) next
      r[, locus := if (bl) sub(".*locus", "", sub(paste0(":", sp$term, "|", sp$term, ":"), "", term_row)) else "Pooled"]
      cnt <- dat[, .(n_persons = uniqueN(id), n_rows = .N), by = locus]
      r <- merge(r, rbind(cnt, data.table(locus = "Pooled", n_persons = uniqueN(dat$id), n_rows = nrow(dat))), by = "locus", all.x = TRUE)
      res[[length(res) + 1]] <- r[, `:=`(outcome = OUT[oc], subset = ss, model = k, model_label = sp$label, unit = sp$unit)]
    }
  }
}
res <- rbindlist(res, fill = TRUE)[, term_row := NULL]
setcolorder(res, c("outcome", "subset", "model", "locus", "est", "se", "p", "n_persons", "n_rows", "unit", "model_label"))
fwrite(res, file.path(HERE, "longitudinal_coefficients.csv"))
cl <- rbindlist(CLEAN_LOG)[, .N, by = .(outcome, cohort, why)][order(outcome, cohort, why)]
fwrite(cl, file.path(HERE, "cleaning_log.csv")); cat("
Visits dropped by cleaning (summed over subsets):
"); print(cl)
print(res[, .(outcome, subset, model, locus, est = round(est, 3), se = round(se, 3), p = signif(p, 2), n_persons, n_rows)], nrows = 300)

# ---------------------------------------------------------------- figure
d <- long_data("HEIGHT")
pA <- ggplot(d, aes(age, Z, group = id)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(colour = "grey70", linewidth = 0.3) + geom_point(aes(colour = BMI), size = 1.1) +
  scale_colour_gradient2(low = "#2C7BB6", mid = "grey85", high = "#D7191C", midpoint = 0, limits = c(-3, 4), oob = scales::squish, name = "BMI Z") +
  facet_wrap(~locus, nrow = 1) + coord_cartesian(ylim = c(-4, 4)) +
  labs(x = "Age (years)", y = "Height Z", title = "A. Height trajectories of carriers with repeat visits (< 18 y), points coloured by BMI Z at that visit") +
  theme_classic(base_size = 10)
iv <- intervals(d)
pB <- ggplot(iv, aes(dB, dZ)) + geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(colour = cohort), size = 1.2, alpha = 0.7) + geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8) +
  facet_wrap(~locus, nrow = 1) + coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
  scale_colour_manual(values = c(Cardiff = "#E69F00", SVIP = "#0072B2", Searchlight = "#009E73")) +
  labs(x = "Change in BMI Z between consecutive visits", y = "Change in height Z", colour = NULL,
       title = "B. Change in height Z vs change in BMI Z over the same interval (\u2265 0.5 y apart)") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
fd <- res[subset == "all"][, model := factor(model, levels = c("W", "P", "C", "R"), labels = c("Within-person\nBMI \u2192 Z", "Prospective\nBMI(t0) \u2192 \u0394Z", "Concurrent\n\u0394BMI \u2192 \u0394Z", "Reverse\nZ(t0) \u2192 \u0394BMI"))]
fd[, locus := factor(locus, levels = rev(c("Pooled", LOCI)))]
pC <- ggplot(fd, aes(est, locus, colour = locus == "Pooled")) + geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), width = 0.2, orientation = "y") + geom_point(size = 2) +
  geom_text(aes(label = sprintf("%+.2f%s", est, ifelse(p < 0.05, "*", ""))), nudge_y = 0.35, size = 2.6, show.legend = FALSE) +
  facet_grid(outcome ~ model, scales = "free_x") + scale_colour_manual(values = c("grey30", "#C0392B"), guide = "none") +
  labs(x = "Coefficient (Z per 1 SD of predictor), 95% CI", y = NULL,
       title = "C. Mixed-model coefficients, carriers < 18 y (all cohorts, cleaned visits); * p < 0.05") + theme_bw(base_size = 10)
fig <- pA / pB / pC + plot_layout(heights = c(1, 1, 1.5))
ggsave(file.path(HERE, "fig_longitudinal.png"), fig, width = 13, height = 13, dpi = 200)
ggsave(file.path(HERE, "fig_longitudinal.pdf"), fig, width = 13, height = 13, device = cairo_pdf)
