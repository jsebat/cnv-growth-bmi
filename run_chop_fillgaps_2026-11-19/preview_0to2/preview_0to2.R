# Head circumference 0-2 y also uses the Simons VIP infant visits of hc_visits.rds (16p11.2 carriers).
# PREVIEW (JS 2026-11-19, Jon: "show me first, then we'll decide"): a 0-2 y band for the BMI effect on head circumference and height.
# Not part of the paper pipeline. Data: analysis_table_chop_2026-11-19.rds (clinical cohorts incl. CHOP A-D). Infants are scored as in
# 00_scoring.R (UK90 head circumference, WHO 2006 length/height and BMI-for-age from birth); |Z| > 5 dropped. One record per person
# within 0-2 y (latest visit under 2 y with the outcome and BMI); these infants may also appear in the older bands of the figure
# (a different visit). Model per genotype: Z ~ BMI Z + sex + age. Older bands are copied from the paper run (expanse_v9e
# main_with_ukbb_bmi_slope_by_band.csv) for comparison. Output: preview_0to2_bands.csv, preview_0to2.png
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
OUT <- file.path(HERE, "run_chop_fillgaps_2026-11-19", "preview_0to2"); setwd(HERE); source("00_scoring.R")
cl <- readRDS("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/analysis_table_chop_2026-11-19.rds")[!is.na(sex) & age >= 0 & age < 2]
cl[, sex := as.character(sex)][, locus := as.character(locus)]
res <- list(); cnt <- list()
for (o in c("Head circumference", "Height")) {
  raw <- if (o == "Height") "HEIGHT_raw" else "HEADCIRCUM_raw"
  x <- cl[!is.na(get(raw)) & !is.na(BMI_raw)][order(id, -age)][, .SD[1], by = id]
  if (o == "Head circumference") {   # Simons VIP infant head-circumference visits (hc_visits.rds, 07_hc_trajectories.R); not in the growth table
    hv <- as.data.table(readRDS("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/hc_visits.rds"))[age >= 0 & age < 2 & !is.na(hc) & !is.na(bmi)]
    hv <- hv[order(id, -age)][, .SD[1], by = id][, .(id = paste0("hcv:", id), cohort = paste0(cohort, " (hc_visits)"), locus, sex = as.character(sex), age, HEADCIRCUM_raw = hc, BMI_raw = bmi)]
    x <- rbind(x, hv, fill = TRUE)
  }
  x[, y := if (o == "Height") score_height(get(raw), age, sex) else score_hc(get(raw), age, sex)][, BMI := score_bmi_child(BMI_raw, age, sex)]
  x <- x[!is.na(y) & !is.na(BMI) & abs(y) <= 5 & abs(BMI) <= 5]
  cnt[[o]] <- x[, .N, by = .(locus, cohort)][, outcome := o]
  for (G in unique(x$locus)) { d <- x[locus == G]; if (nrow(d) < 15) next
    cf <- summary(lm(y ~ BMI + sex + age, data = d))$coefficients["BMI", ]
    res[[length(res) + 1]] <- data.table(outcome = o, group = G, age_bin = "0-2", n = nrow(d), slope = cf[1], se = cf[2], p = cf[4], median_age = median(d$age)) }
}
res <- rbindlist(res); cnt <- rbindlist(cnt); fwrite(res, file.path(OUT, "preview_0to2_bands.csv"))
cat("== 0-2 y sample by cohort\n"); print(dcast(cnt, outcome + locus ~ cohort, value.var = "N", fill = 0))
cat("\n== 0-2 y BMI slope within group\n"); print(res[, .(outcome, group, n, slope = round(slope, 3), se = round(se, 3), p = signif(p, 2), median_age = round(median_age, 2))])
# older bands from the paper run
b <- fread(file.path(HERE, "run_chop_fillgaps_2026-11-19", "expanse_v9e", "main_with_ukbb_bmi_slope_by_band.csv"))
old <- rbind(b[group == "No CNV", .(outcome, group, age_bin, n = n_controls, slope = slope_controls, se = se_controls, who = "Non-carriers")],
             b[group != "No CNV", .(outcome, group, age_bin, n = n_carriers, slope = slope_carriers, se = se_carriers, who = "Carriers")])
new <- res[, .(outcome, group, age_bin, n, slope, se, who = ifelse(group == "No CNV", "Non-carriers", "Carriers"))]
P <- rbind(new, old)[, age_bin := factor(age_bin, c("0-2", "2-9", "10-18", "19-39", "40-73"))][, who := factor(who, c("Non-carriers", "Carriers"))]
P[, band_new := age_bin == "0-2"]
COL <- c(`Non-carriers` = "grey60", Carriers = "#3B6FB6")
pan <- function(o, G) {
  d <- P[outcome == o & group == G]; if (!nrow(d)) return(plot_spacer())
  ggplot(d, aes(age_bin, slope, fill = who, alpha = band_new)) + geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(position = position_dodge(0.8, preserve = "single"), width = 0.75) +
    geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), position = position_dodge(0.8, preserve = "single"), width = 0.25, colour = "grey20", alpha = 1) +
    geom_text(data = d[band_new == TRUE], aes(age_bin, slope + sign(slope) * (1.96 * se + 0.05), label = paste0("n=", n)), size = 2.6, alpha = 1, vjust = ifelse(d[band_new == TRUE]$slope > 0, 0, 1)) +
    scale_fill_manual(values = COL, drop = FALSE, name = NULL) + scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.55), guide = "none") +
    labs(title = if (G == "No CNV") "Non-carrier controls" else G, x = NULL, y = paste(o, "per +1 SD BMI")) + theme_classic(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(size = 9.5))
}
GR <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
fig <- wrap_plots(c(lapply(GR, function(G) pan("Head circumference", G)), lapply(GR, function(G) pan("Height", G))), ncol = 5, byrow = TRUE) +
  plot_annotation(title = "PREVIEW: BMI effect by age band with a new 0-2 y band (lighter bars; within-group slope, no control comparison)",
                  caption = "0-2 y: latest visit under 2 y per person, Z ~ BMI + sex + age within group (groups with >= 15 people). Other bands: paper run v9e (carrier slopes from the pooled carrier-control model).")
ggsave(file.path(OUT, "preview_0to2.png"), fig, width = 15, height = 6.5, dpi = 170, bg = "white"); cat("preview written\n")
