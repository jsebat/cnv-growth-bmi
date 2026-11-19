# Sex difference in the BMI slope on SRS / SCQ in SSC without CNVs (JS / Claude 2026-11-14): autistic probands, designated siblings,
# parents (SRS only). Box-Cox scale (00_boxcox_srs_scq.R). Per group: BMI slope by sex (gam y ~ s(age) + BMI within each sex) and
# the sex difference, from the separately fitted sexes (as for carriers in 66) and from one model (y ~ s(age) + male * BMI).
# Output: ssc_sex_x_bmi.csv (aggregate statistics only).
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE", SRSSCQ_SCALE = "boxcox")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R")); source(file.path(HERE, "24_cognition_bmi_slopes.R"))
S <- X[cohort == "SSC" & pheno %in% c("SRS", "SCQ")]
S[, g := fcase(grp == "SSC ASD probands", "Autistic probands", grepl("[.]s[0-9]+$", id), "Siblings", grepl("[.](mo|fa)$", id), "Parents", default = NA_character_)][, male := as.numeric(sex == "male")]
fit <- function(d, f, v) { m <- gam(f, data = d, method = "REML"); co <- summary(m)$p.table[v, ]; c(b = co[[1]], se = co[[2]], p = co[[4]]) }
out <- list()
for (ph in c("SRS", "SCQ")) for (gg in c("Autistic probands", "Siblings", "Parents")) {
  d <- S[pheno == ph & g == gg]; if (nrow(d) < 100) next
  f <- fit(d[male == 0], y ~ s(age, k = 6) + BMI, "BMI"); m <- fit(d[male == 1], y ~ s(age, k = 6) + BMI, "BMI")
  dif <- m[["b"]] - f[["b"]]; se <- sqrt(f[["se"]]^2 + m[["se"]]^2); j <- fit(d, y ~ s(age, k = 6) + male * BMI, "male:BMI")
  out[[length(out) + 1]] <- data.table(trait = ph, group = gg, n_female = sum(d$male == 0), n_male = sum(d$male == 1), slope_female = f[["b"]], p_female = f[["p"]], slope_male = m[["b"]], p_male = m[["p"]],
    diff_male_minus_female = dif, diff_p_separate = 2 * pnorm(-abs(dif / se)), diff_joint = j[["b"]], diff_p_joint = j[["p"]]) }
R <- rbindlist(out); fwrite(R, file.path(HERE, "ssc_sex_x_bmi.csv")); options(width = 220)
print(R[, .(trait, group, nF = n_female, nM = n_male, slopeF = round(slope_female, 3), pF = signif(p_female, 2), slopeM = round(slope_male, 3), pM = signif(p_male, 2), diff = round(diff_male_minus_female, 3), p_sep = signif(diff_p_separate, 2), p_joint = signif(diff_p_joint, 2))])
