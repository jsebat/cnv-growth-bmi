# Grant figure (JS / Claude 2026-11-15): Fig. 5 v8 (35d_main_figure_v8_adultbc.R) without the height and SCQ columns:
#   Head circumference | HC BMI effect by age band | SRS. Row labels (Non-carrier controls, 16p11.2 del, ...) bold and larger.
#   Top HC panel (controls) carries the four CNVs pooled high - low BMI offset vs controls for head circumference (offset_vs_ctl of
#   17_main_with_ukbb.R: CNV x BMI tail, age-adjusted; analogue of the pooled SRS p in panel k).
# RUNS ON EXPANSE in main_with_ukbb/v8_grant (sources 17, which needs the UKBB genomatrix). Output: fig_grant_hc_srs.{png,pdf},
#   grant_hc_pooled_cnv_x_bmi.csv
source("17_main_with_ukbb.R")
load(file.path(WD, "cog_panels_boxcox_combined.RData"))
LOCI_C <- LOCI; LOC <- LOCI; BR <- c(-Inf, -1, 0, 1, 2, Inf)
# pooled BMI offset vs controls for HC (Jon 2026-11-15: "that is on the BMI offset vs controls"): offset_vs_ctl() of 17, i.e. the
# CNV x BMI-tail (>= 70th vs <= 30th centile) interaction, carriers of all four CNVs pooled vs non-carrier controls
po <- offset_vs_ctl(dh[locus %in% LOCI], dh[locus == "No CNV"]); cp <- c(po[["diff"]], NA, NA, po[["p"]])
# heterogeneity of the offset across the four CNVs (carriers only, locus x BMI tail, age-adjusted): CSV only
xc <- dh[locus %in% LOCI][, grp := tail_of(BMI)][!is.na(grp)][, `:=`(hi = as.numeric(grp == names(COLS)[1]), loc = factor(locus))]
mh <- bam(y ~ loc * hi + s(age, k = 10, by = loc), data = xc, method = "REML")
ph <- anova(mh)$pTerms.table["loc:hi", "p-value"]
oc <- offset(dh[locus %in% LOCI]); o0 <- offset(dh[locus == "No CNV"])
pool <- data.table(outcome = "Head circumference", n_controls = dh[locus == "No CNV", .N], n_carriers = dh[locus %in% LOCI, .N],
                   offset_carriers = oc[["offset"]], offset_controls = o0[["offset"]], offset_vs_ctl = po[["diff"]], p = po[["p"]], p_heterogeneity_4cnv = ph)
fwrite(pool, file.path(WD, "grant_hc_pooled_offset_vs_controls.csv")); print(pool)
# per-CNV values already in the figure source: hcw (main_with_ukbb_hc_tests.csv)
print(hcw[, .(group, offset_vs_ctl, offset_vs_ctl_p)])
SUBT <- theme(plot.subtitle = element_text(face = "bold", size = 17))
tr_c <- lapply(GROUPS, function(G) v2[[paste("HC", G)]] + labs(tag = NULL) + SUBT)
plab <- gsub("= < ", "< ", sprintf("4 CNVs pooled: high − low BMI offset\nvs controls %+.2f SD, p = %s", cp[[1]], pfx(cp[[4]])), fixed = TRUE)
tr_c[[1]] <- tr_c[[1]] + annotate("rect", xmin = 2.3, xmax = 36, ymin = 2.15, ymax = 3.65, fill = "white", alpha = 0.85, colour = NA) +
  annotate("text", x = 3, y = 3.55, label = plab, hjust = 0, vjust = 1, size = 4, lineheight = 0.95, fontface = if (cp[[4]] < 0.05) "bold" else "plain") +
  theme(legend.position.inside = c(0.8, 0.1))   # legend moved to bottom right so the pooled test can sit at the top
ef_c <- lapply(GROUPS, function(G) v2[[paste("EFF HC", G)]] + labs(tag = NULL))
srs <- c(list(panel_ctl("SRS", probands = TRUE, show_title = FALSE)), lapply(LOCI_C, function(g) panel_cnv("SRS", g)))
srs <- lapply(srs, function(p) p + theme(plot.subtitle = element_text(face = "bold", size = 14)))
cols <- list(tr_c, ef_c, srs)
tagsets <- list(letters[1:5], paste0(letters[1:5], "'"), letters[6:10])
plots <- unlist(lapply(seq_along(cols), function(j) Map(function(p, tg) p + labs(tag = tg), cols[[j]], tagsets[[j]])), recursive = FALSE)
hdr5 <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 9, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
W <- c(2, 1.3, 1.3)
hdr_row <- (hdr5("Head circumference") | hdr5("") | hdr5("SRS")) + plot_layout(widths = W)
body <- wrap_plots(plots, ncol = 3, byrow = FALSE, widths = W) & theme(plot.tag = element_text(size = 18, face = "bold"))
fig <- hdr_row / body + plot_layout(heights = c(0.04, 1))
ggsave(file.path(WD, "fig_grant_hc_srs.png"), fig, width = 13.5, height = 20.5, dpi = 150)
tryCatch(ggsave(file.path(WD, "fig_grant_hc_srs.pdf"), fig, width = 13.5, height = 20.5, device = cairo_pdf), error = function(e) ggsave(file.path(WD, "fig_grant_hc_srs.pdf"), fig, width = 13.5, height = 20.5))
cat("fig_grant_hc_srs written\n")
