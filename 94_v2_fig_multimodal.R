# 94_v2_fig_multimodal.R  (JS / Claude 2026-11-21; BMI modifier paper v2, finale figure)
# What it does: draws main Fig. 6 (autistic traits in CNV carriers from growth, sex and genetic background: the core "modifiers"
#   model) and Supplementary Fig. 7 (separation within genotype x cohort, and the autism-diagnosis version of the model), from the
#   aggregate tables written by 92_v2_srs_multimodal.R and 93_v2_asd_multimodal.R. Drawing only: no statistics are computed here.
#   Fig. 6  a  out-of-sample R2 of the benchmark (genotype x cohort + age + SRS form), the core model and the full model
#           b  importance: loss in out-of-sample R2 when each factor is dropped (core model), and for the measures that are in the
#              full model only
#           c  BMI slope on SRS in girls and boys, and at low and high autism polygenic score (core model; estimate +/- s.e.)
#           d  the two modifiers per CNV and pooled (estimate +/- s.e.), with heterogeneity across CNVs
#   Supp. Fig. 7  a  SRS beyond the benchmark by fifth of the core added-factor prediction (within genotype x cohort, age, form)
#                 b  autism diagnosis: out-of-sample AUC of benchmark, core and full models
#                 c  autism diagnosis: loss in AUC when each factor is dropped (core model)
#                 d  share with an autism diagnosis by fifth of the core diagnosis prediction
# Inputs : run_v2_2026-11-20/multimodal/{srs_core_*, srs_multimodal_*, asd_multimodal_*}.csv
# Outputs: run_v2_2026-11-20/multimodal/fig_multimodal_main.{png,pdf}, fig_supp_multimodal.{png,pdf}, and *_plotted.csv for each.
# null model: none here (the permutation null of the R2 gain is in 92_v2_srs_multimodal.R).
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
MM <- file.path(HERE, "run_v2_2026-11-20", "multimodal"); V2OUT <- MM
rd <- function(f) fread(file.path(MM, f))
BS <- 16; TS <- 4.6; C_LOW <- "#4B6EA9"; C_HIGH <- "#D9822B"; C_EXTRA <- "grey60"
LOC4 <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
# ---------------------------------------------------------------- a: out-of-sample R2
r2c <- rd("srs_core_r2.csv"); pr <- rd("srs_multimodal_prediction.csv")[quantity == "out-of-sample R2"]; acc <- rd("srs_core_accuracy.csv")
A <- rbind(data.table(m = "Genotype \u00d7 cohort, age,\nSRS form (benchmark)", est = r2c$r2[1], lo = r2c$lo[1], hi = r2c$hi[1], k = "bench"),
           data.table(m = "+ sex, BMI, autism PGS,\nBMI \u00d7 sex, BMI \u00d7 PGS (core)", est = r2c$r2[2], lo = r2c$lo[2], hi = r2c$hi[2], k = "core"),
           data.table(m = "+ also height, head size,\nadiposity PGS, BMI \u00d7 CNV (full)", est = pr[grepl("^full", model), est], lo = pr[grepl("^full", model), lo], hi = pr[grepl("^full", model), hi], k = "full"))
stopifnot(abs(pr[grepl("^benchmark: genotype", model), est] - r2c$r2[1]) < 0.002)   # the core and full runs use different random fold draws (benchmark 0.2160 vs 0.2162); each gain is against its own benchmark
A[, m := factor(m, rev(m))]
gc <- r2c[3]; gf <- pr[grepl("^gain", model)]
pa <- ggplot(A, aes(est, m, colour = k)) + geom_vline(xintercept = A[k == "bench", est], colour = "grey70", linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 1.3) + geom_point(size = 5) +
  annotate("text", x = 0.2445, y = 2, hjust = 1, vjust = -0.9, size = TS, fontface = "bold", label = sprintf("gain +%.3f (%.3f\u2013%.3f)\npermutation P < 0.001", gc$r2, gc$lo, gc$hi)) +
  annotate("text", x = 0.2445, y = 1, hjust = 1, vjust = -0.9, size = TS, colour = "grey35", label = sprintf("gain +%.3f (%.3f to %.3f)", gf$est, gf$lo, gf$hi)) +
  scale_colour_manual(values = c(bench = "grey25", core = C_CAR, full = C_EXTRA), guide = "none") + coord_cartesian(xlim = c(0.19, 0.245)) +
  labs(x = "Out-of-sample R\u00b2 for SRS", y = NULL, tag = "a", title = "Prediction in carriers left out of the fit",
       subtitle = sprintf("%d carriers; 10-fold cross-validation by family, 50 repeats\n(range = 2.5\u201397.5%% of repeats)", acc[1, n])) + th_v2(BS)
# ---------------------------------------------------------------- b: importance (drop-one)
dc <- rd("srs_core_drop_analysis.csv"); df <- rd("srs_multimodal_drop_analysis.csv")
LBL <- c(BMI = "BMI (with its interactions)", `ASD PGS` = "Autism PGS (with BMI \u00d7 PGS)", `BMI x ASD PGS` = "BMI \u00d7 autism PGS", `BMI x sex` = "BMI \u00d7 sex",
         Sex = "Sex (with BMI \u00d7 sex)", Height = "Height", `Head circumference` = "Head circumference", `Adiposity and height PGS` = "Adiposity and height PGS", `BMI x CNV` = "BMI \u00d7 CNV")
B <- rbind(dc[, .(f = factor, est = loss_r2, lo, hi, set = "Core model")],
           df[factor %in% c("Height", "Head circumference", "Adiposity and height PGS", "BMI x CNV"), .(f = factor, est = loss_r2, lo, hi, set = "Added in the full model")])
stopifnot(all(B$f %in% names(LBL))); B[, lab := factor(LBL[f], rev(LBL[f]))][, set := factor(set, c("Core model", "Added in the full model"))]
pb <- ggplot(B, aes(1000 * est, lab, fill = set)) + geom_col(width = 0.68) + geom_vline(xintercept = 0, colour = "grey30") +
  geom_errorbar(aes(xmin = 1000 * lo, xmax = 1000 * hi), width = 0.22, orientation = "y", linewidth = 0.7) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") + scale_fill_manual(values = c(`Core model` = C_CAR, `Added in the full model` = C_EXTRA), guide = "none") +
  labs(x = "Loss in out-of-sample R\u00b2 when dropped (\u00d710\u207b\u00b3)", y = NULL, tag = "b", title = "Which factors carry the prediction",
       subtitle = "Below zero: the model predicts better without it") + th_v2(BS) +
  theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, size = 13))
# ---------------------------------------------------------------- c: BMI slope by sex and by autism PGS
sl <- rd("srs_core_bmi_slopes_by_sex_and_pgs.csv")
S <- sl[group != "interaction"][, g := c(girls = "Girls", boys = "Boys", `low ASD PGS (-1 SD)` = "Low PGS\n(\u22121 SD)", `high ASD PGS (+1 SD)` = "High PGS\n(+1 SD)")[group]]
stopifnot(!anyNA(S$g)); S[, g := factor(g, g)][, panel := factor(fifelse(panel == "sex", "By sex", "By autism polygenic score"), c("By sex", "By autism polygenic score"))]
I <- sl[group == "interaction"][, panel := factor(fifelse(panel == "sex", "By sex", "By autism polygenic score"), levels(S$panel))]
I[, lab := sprintf("%s = %+.3f\n%s", fifelse(panel == "By sex", "BMI \u00d7 sex", "BMI \u00d7 PGS"), est, plab(p))]
CS <- c(Girls = C_F, Boys = C_M); CS[levels(S$g)[3]] <- C_LOW; CS[levels(S$g)[4]] <- C_HIGH
pc <- ggplot(S, aes(g, est, fill = g)) + geom_col(width = 0.66) + geom_hline(yintercept = 0, colour = "grey30") +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.16, linewidth = 0.8) +
  geom_text(aes(y = pmax(est + se, 0) + 0.012, label = plab(p), fontface = face(p)), size = TS - 0.6, vjust = 0) +
  geom_text(data = I, aes(x = 1.5, y = 0.255, label = lab, fontface = face(p)), inherit.aes = FALSE, size = TS - 0.2, lineheight = 0.95) +
  facet_wrap(~panel, scales = "free_x") + scale_fill_manual(values = CS, guide = "none") + coord_cartesian(ylim = c(-0.06, 0.29)) +
  labs(x = NULL, y = "BMI slope on SRS (SD per SD, \u00b1 s.e.)", tag = "c", title = "BMI tracks autistic traits in girls\nand at a low autism PGS",
       subtitle = "Core model; slopes at the other factor's sample mean") + th_v2(BS) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14))
# ---------------------------------------------------------------- d: modifiers per CNV and pooled
pc_ <- rd("srs_core_per_cnv.csv")[factor %in% c("BMI x sex (male)", "BMI x ASD PGS")]
D <- rbind(pc_[, .(term = factor, g = locus, n, est, se, p, p_het)],
           data.table(term = c("BMI x sex (male)", "BMI x ASD PGS"), g = "All four CNVs", n = c(sl[panel == "sex" & group == "interaction", n], sl[panel == "ASD PGS" & group == "interaction", n]),
                      est = c(sl[panel == "sex" & group == "interaction", est], sl[panel == "ASD PGS" & group == "interaction", est]),
                      se = c(sl[panel == "sex" & group == "interaction", se], sl[panel == "ASD PGS" & group == "interaction", se]),
                      p = c(sl[panel == "sex" & group == "interaction", p], sl[panel == "ASD PGS" & group == "interaction", p]), p_het = NA_real_))
D[, g := factor(g, c(LOC4, "All four CNVs"))][, term := factor(c(`BMI x sex (male)` = "BMI \u00d7 sex (boys minus girls)", `BMI x ASD PGS` = "BMI \u00d7 autism PGS")[term], c("BMI \u00d7 sex (boys minus girls)", "BMI \u00d7 autism PGS"))]
H <- D[!is.na(p_het), .(p_het = p_het[1]), by = term][, lab := sprintf("heterogeneity across CNVs %s", plab(p_het))]
D[, xl := factor(ifelse(g == "All four CNVs", "All four\nCNVs", sub(" ", "\n", as.character(g))), c(sub(" ", "\n", LOC4), "All four\nCNVs"))]
pd <- ggplot(D, aes(xl, est, fill = g == "All four CNVs")) + geom_col(width = 0.66) + geom_hline(yintercept = 0, colour = "grey30") + geom_vline(xintercept = 4.5, colour = "grey80") +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.16, linewidth = 0.8) +
  geom_text(aes(y = 0.33, label = plab(p), fontface = face(p)), size = TS - 0.9) +
  geom_text(data = H, aes(x = 2.5, y = -0.42, label = lab), inherit.aes = FALSE, size = TS - 0.6, colour = "grey25") +
  facet_wrap(~term) + scale_fill_manual(values = c(`FALSE` = "#9DB1D3", `TRUE` = C_CAR), guide = "none") + coord_cartesian(ylim = c(-0.45, 0.36)) +
  labs(x = NULL, y = "Interaction on SRS (SD per SD, \u00b1 s.e.)", tag = "d", title = "The two modifiers in each CNV", subtitle = "Core model fitted within each CNV and pooled") +
  th_v2(BS) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14), axis.text.x = element_text(size = 11.5))
fig <- ((pa | pb) + plot_layout(widths = c(1, 1.1))) / ((pc | pd) + plot_layout(widths = c(0.85, 1.25))) + plot_layout(heights = c(1, 1.1))
save_v2(fig, "fig_multimodal_main", 20, 13)
fwrite(rbind(A[, .(panel = "a", series = gsub("\n", " ", m), est, lo, hi, se = NA_real_, p = NA_real_)], B[, .(panel = "b", series = paste(set, f), est, lo, hi, se = NA_real_, p = NA_real_)],
             S[, .(panel = "c", series = gsub("\n", " ", g), est, lo = NA_real_, hi = NA_real_, se, p)], D[, .(panel = "d", series = paste(term, g), est, lo = NA_real_, hi = NA_real_, se, p)]),
       file.path(MM, "fig_multimodal_main_plotted.csv"))
# ---------------------------------------------------------------- supplementary figure
wc <- rd("srs_core_within_cell_by_fifth.csv"); sep <- rd("srs_core_separation.csv")
sa <- ggplot(wc, aes(factor(fifth), srs_vs_bench)) + geom_col(width = 0.66, fill = C_CAR) + geom_hline(yintercept = 0, colour = "grey30") +
  geom_errorbar(aes(ymin = srs_vs_bench - se, ymax = srs_vs_bench + se), width = 0.16) +
  annotate("text", x = 1, y = 0.42, hjust = 0, size = TS - 0.3, label = sprintf("top minus bottom fifth +%.2f SD (95%% CI %.2f\u2013%.2f; family bootstrap)", sep$est[1], sep$lo[1], sep$hi[1])) +
  coord_cartesian(ylim = c(-0.35, 0.45)) + scale_x_discrete(labels = c("1\nlowest", "2", "3", "4", "5\nhighest")) +
  labs(x = "Fifth of the prediction from sex, BMI, autism PGS and the two modifiers", y = "SRS beyond genotype \u00d7 cohort,\nage and form (SD, \u00b1 s.e.)", tag = "a",
       title = "Separation among carriers of the same CNV and cohort") + th_v2(BS)
ap <- rd("asd_multimodal_performance.csv")[metric == "AUC"]; ap[, m := factor(c(benchmark = "Genotype \u00d7 cohort, age\n(benchmark)", core = "+ core factors", full = "+ all measures (full)")[model], rev(c("Genotype \u00d7 cohort, age\n(benchmark)", "+ core factors", "+ all measures (full)")))]
ag <- rd("asd_multimodal_performance.csv")[model == "core - benchmark" & metric == "AUC gain"]; smp <- rd("asd_multimodal_sample.csv")
sb <- ggplot(ap, aes(mean, m, colour = model)) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 1.3) + geom_point(size = 5) +
  scale_colour_manual(values = c(benchmark = "grey25", core = C_CAR, full = C_EXTRA), guide = "none") +
  annotate("text", x = min(ap$lo), y = 3.42, hjust = 0, size = TS - 0.3, label = sprintf("core gain +%.3f (%.3f\u2013%.3f)", ag$mean, ag$lo, ag$hi)) +
  labs(x = "Out-of-sample AUC for autism diagnosis", y = NULL, tag = "b", title = "Autism diagnosis: prediction",
       subtitle = sprintf("%s carriers, %s with a diagnosis; SSC carriers excluded", format(smp[cohort == "All", carriers], big.mark = ","), smp[cohort == "All", autism])) + th_v2(BS)
ad <- rd("asd_multimodal_drop_analysis.csv")[model == "core"][, lab := factor(LBL[factor], rev(LBL[factor]))]
sc <- ggplot(ad, aes(1000 * auc_loss, lab)) + geom_col(width = 0.66, fill = C_CAR) + geom_vline(xintercept = 0, colour = "grey30") +
  geom_errorbar(aes(xmin = 1000 * auc_loss_lo, xmax = 1000 * auc_loss_hi), width = 0.22, orientation = "y") +
  labs(x = "Loss in out-of-sample AUC when dropped (\u00d710\u207b\u00b3)", y = NULL, tag = "c", title = "Autism diagnosis: the gain is almost all sex") + th_v2(BS)
af <- rd("asd_multimodal_by_fifth.csv"); asep <- rd("asd_multimodal_separation.csv")
sd_ <- ggplot(af, aes(factor(fifth), 100 * autism)) + geom_col(width = 0.66, fill = C_CAR) +
  geom_errorbar(aes(ymin = 100 * (autism - sqrt(autism * (1 - autism) / n)), ymax = 100 * (autism + sqrt(autism * (1 - autism) / n))), width = 0.16) +
  geom_text(aes(y = 1.5, label = paste0("n = ", n)), colour = "white", size = TS - 1) + scale_x_discrete(labels = c("1\nlowest", "2", "3", "4", "5\nhighest")) +
  labs(x = "Fifth of the core added-factor prediction", y = "With an autism diagnosis (%, \u00b1 s.e.)", tag = "d", title = "Autism diagnosis by fifth of the prediction",
       subtitle = sprintf("Top minus bottom fifth %+.1f percentage points (95%% CI %.1f\u2013%.1f)", 100 * asep$est[1], 100 * asep$lo[1], 100 * asep$hi[1])) + th_v2(BS)
sfig <- (sa | sb) / (sc | sd_)
save_v2(sfig, "fig_supp_multimodal", 20, 12)
fwrite(rbind(wc[, .(panel = "a", series = paste("fifth", fifth), est = srs_vs_bench, se, n)], ap[, .(panel = "b", series = model, est = mean, se = NA_real_, n = NA_integer_)],
             ad[, .(panel = "c", series = factor, est = auc_loss, se = NA_real_, n = NA_integer_)], af[, .(panel = "d", series = paste("fifth", fifth), est = autism, se = NA_real_, n)]),
       file.path(MM, "fig_supp_multimodal_plotted.csv"))
cat("94 done\n")
