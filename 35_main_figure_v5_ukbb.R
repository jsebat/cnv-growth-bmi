# Main figure v5 (JS / Claude 2026-10-27): the v4 six-column layout (Height | Height' | Head circumference | HC' | SRS | SCQ)
# built from the UK Biobank version of the growth panels (17_main_with_ukbb.R v2: UKBB adults in the height panels, adults
# scored relative to their own cohort's non-carriers with the UKBB Box-Cox lambdas, prime panels = BMI effect by age) plus
# the SRS/SCQ columns of 24_cognition_bmi_slopes.R, whose data and panel functions are loaded from cog_panels.RData
# (written locally by 35a_save_cog_panels.R and copied here). RUNS ON EXPANSE in the main_with_ukbb working directory:
# sourcing 17 rebuilds every object it needs (H, dh, bea, ht, traj2, eff_panel, YLH, YLC, HCMAX, ...).
# Tags: A-E / A'-E' height, F-J / F'-J' head circumference, K-O SRS, P-T SCQ (as in 34_main_figure_v4.R).
# Output (WD): fig_main_v5_ukbb.{png,pdf}
source("17_main_with_ukbb.R")
load(file.path(WD, "cog_panels.RData"))   # X, R, panel_ctl, panel_cnv, fitline, bins, ylim_of, ptxt, REFC, COLS3, XL, YL
LOCI_C <- LOCI
cog <- function(ph) c(list(panel_ctl(ph, probands = TRUE, show_title = FALSE)), lapply(LOCI_C, function(g) panel_cnv(ph, g)))
# growth columns: trajectory + prime panels from 17 (v2 objects), in the v4 order; the UKBB trajectory panels keep their
# own labels (CNV x BMI p from the age-adjusted model of 17; CNV x age x BMI p in the prime panel)
tr_h <- lapply(GROUPS, function(G) v2[[paste("HEIGHT", G)]] + labs(tag = NULL)); ef_h <- lapply(GROUPS, function(G) v2[[paste("EFF HEIGHT", G)]] + labs(tag = NULL))
tr_c <- lapply(GROUPS, function(G) v2[[paste("HC", G)]] + labs(tag = NULL));     ef_c <- lapply(GROUPS, function(G) v2[[paste("EFF HC", G)]] + labs(tag = NULL))
cols <- list(tr_h, ef_h, tr_c, ef_c, cog("SRS"), cog("SCQ"))
tagsets <- list(LETTERS[1:5], paste0(LETTERS[1:5], "'"), LETTERS[6:10], paste0(LETTERS[6:10], "'"), LETTERS[11:15], LETTERS[16:20])
plots <- unlist(lapply(seq_along(cols), function(j) Map(function(p, tg) p + labs(tag = tg), cols[[j]], tagsets[[j]])), recursive = FALSE)
hdr5 <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 9, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
W <- c(2, 1.3, 2, 1.3, 1, 1)
hdr_row5 <- (hdr5("Height") | hdr5("") | hdr5("Head circumference") | hdr5("") | hdr5("SRS") | hdr5("SCQ")) + plot_layout(widths = W)
body5 <- wrap_plots(plots, ncol = 6, byrow = FALSE, widths = W) & theme(plot.tag = element_text(size = 18, face = "bold"))
fig5 <- hdr_row5 / body5 + plot_layout(heights = c(0.04, 1))
ggsave(file.path(WD, "fig_main_v5_ukbb.png"), fig5, width = 22, height = 20.5, dpi = 150)
tryCatch(ggsave(file.path(WD, "fig_main_v5_ukbb.pdf"), fig5, width = 22, height = 20.5, device = cairo_pdf), error = function(e) ggsave(file.path(WD, "fig_main_v5_ukbb.pdf"), fig5, width = 22, height = 20.5))
cat("fig_main_v5_ukbb written\n")
