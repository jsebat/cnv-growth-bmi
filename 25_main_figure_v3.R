# Main figure v3 (JS / Claude 2026-10-24): Height | Head circumference | SRS | SCQ, rows = unaffected controls, 16p11.2 del,
# 16p11.2 dup, 22q11.2 del, 22q11.2 dup. Clinical cohorts only (no UKBB); controls = unaffected (non-ASD) SSC siblings and
# parents (+ SVIP non-ASD family members for height/HC). Height/HC: trajectories by BMI tail vs age, p-values for CNV x BMI
# (age-adjusted, continuous BMI) and CNV x age x BMI vs unaffected (23_unaffected_controls.R). SRS/SCQ: age-adjusted trait vs
# continuous BMI Z with slope, binned means and the unaffected slope as a dashed reference; p-values for the within-carrier
# BMI slope and for the slope difference vs unaffected (24_cognition_bmi_slopes.R). Each p-value line is bold only if < 0.05.
# Output: fig_main_v3.{png,pdf}, fig_main_v3_legend.txt
Sys.setenv(UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "23_unaffected_controls.R"))      # A, res, panel(), LOCI, COLS, pfx, D
source(file.path(HERE, "24_cognition_bmi_slopes.R"))     # X, R, panel_ctl(), panel_cnv()
GR1 <- c("Unaffected controls", LOCI)
growth <- function(ph) lapply(GR1, function(Gn) panel(ph, Gn, FALSE, "controls"))
cog <- function(ph) c(list(panel_ctl(ph, probands = TRUE, show_title = FALSE)), lapply(LOCI, function(g) panel_cnv(ph, g)))
blank <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 3.6, colour = "grey40") + theme_void()
cog_nc <- function(ph) c(list(blank("No non-carrier
controls with BMI")), lapply(LOCI, function(g) panel_cnv(ph, g)))   # CNB: carriers only
plots <- c(growth("Height"), growth("Head circumference"), cog("SRS"), cog("SCQ"))   # JS 2026-10-24: CNB columns tried and dropped (26_supp_cnb_bmi.R)
tags <- c(LETTERS, paste0("A", LETTERS))[seq_along(plots)]   # 30 panels: A-Z, AA-AD
plots <- Map(function(p, tg) p + labs(tag = tg), plots, tags)
hdr <- function(txt) ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, size = 9, fontface = "bold") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
hdr_row <- (hdr("Height") | hdr("Head circumference") | hdr("SRS") | hdr("SCQ")) + plot_layout(widths = c(2, 2, 1, 1))
body <- wrap_plots(plots, ncol = 4, byrow = FALSE, widths = c(2, 2, 1, 1)) & theme(plot.tag = element_text(size = 18, face = "bold"))
fig <- hdr_row / body + plot_layout(heights = c(0.04, 1))
ggsave(file.path(HERE, "fig_main_v3.png"), fig, width = 16, height = 20.5, dpi = 150)
tryCatch(ggsave(file.path(HERE, "fig_main_v3.pdf"), fig, width = 16, height = 20.5, device = cairo_pdf), error = function(e) ggsave(file.path(HERE, "fig_main_v3.pdf"), fig, width = 16, height = 20.5))


# ---------------------------------------------------------------- heterogeneity of the BMI slope across the four CNV groups
# JS 2026-10-24: within carriers, y ~ s(age, by = locus) + locus + BMI + locus:BMI; Wald test of the locus:BMI parametric
# terms (mgcv anova.gam). Height/HC from A (23), SRS/SCQ from X (24; age-adjusted trait, same model on the raw trait y).
het_test <- function(x) tryCatch({
  x <- copy(x)[, locus := factor(locus)]
  m <- gam(y ~ s(age, k = 8, by = locus) + locus + BMI + locus:BMI, data = x, method = "REML")
  a <- anova(m)$pTerms.table; c(F = unname(a["locus:BMI", "F"]), df = unname(a["locus:BMI", "df"]), p = unname(a["locus:BMI", "p-value"])) }, error = function(e) c(F = NA, df = NA, p = NA))
HET <- rbind(data.table(pheno = "Height", t(het_test(A[pheno == "Height" & locus %in% LOCI]))),
             data.table(pheno = "Head circumference", t(het_test(A[pheno == "Head circumference" & locus %in% LOCI]))),
             data.table(pheno = "SRS", t(het_test(X[pheno == "SRS" & grp %in% LOCI, .(y, age, BMI, locus = grp)]))),
             data.table(pheno = "SCQ", t(het_test(X[pheno == "SCQ" & grp %in% LOCI, .(y, age, BMI, locus = grp)]))),
             data.table(pheno = "CNB accuracy", t(het_test(X[pheno == "CNB accuracy" & grp %in% LOCI, .(y, age, BMI, locus = grp)]))),
             data.table(pheno = "CNB speed", t(het_test(X[pheno == "CNB speed" & grp %in% LOCI, .(y, age, BMI, locus = grp)]))))
fwrite(HET, file.path(HERE, "main_v3_bmi_slope_heterogeneity.csv")); cat("\n== heterogeneity of BMI slope across CNV groups\n"); print(HET)
het_txt <- paste(sprintf("%s p = %s", HET$pheno, sapply(HET$p, pfx)), collapse = ", ")
# ---------------------------------------------------------------- legend text
n_of <- function(ph, Gn) nrow(A[pheno == ph & locus == Gn & !is.na(tail)])
nX <- function(ph, g) nrow(X[pheno == ph & grp == g])
leg <- c(
"Figure. Body mass and growth and social-communication phenotypes in carriers of 16p11.2 and 22q11.2 CNVs versus unaffected relatives.",
"Columns: height (A-E), head circumference (F-J), SRS (K-O) and SCQ (P-T). Rows: unaffected non-carrier controls (top; SSC designated siblings and parents, plus SVIP non-carrier family members without ASD for height and head circumference), then 16p11.2 deletion, 16p11.2 duplication, 22q11.2 deletion and 22q11.2 duplication carriers (G2MH, SVIP, Searchlight, Cardiff). No UK Biobank data are included.",
sprintf("A-J: height and head circumference Z (WHO/UK90 references) against age, for individuals with BMI at or above the 70th centile (red) or at or below the 30th centile (teal); GAM smooths with 95%% confidence bands. Sibling and parent ages are ages at the study visit. n: controls %d (height) and %d (head circumference); carriers %s (height) and %s (head circumference). P-values: CNV x BMI is the age-adjusted difference between carriers and controls in the slope of the outcome on continuous BMI Z (GAM y ~ CNV + s(age, by CNV) + s(age, by BMI) + CNV:BMI); CNV x age x BMI tests whether that difference varies with age (tensor-product interaction). Tests use all individuals with BMI, not only the plotted tails.",
        n_of("Height", "Unaffected controls"), n_of("Head circumference", "Unaffected controls"), paste(sapply(LOCI, function(L) n_of("Height", L)), collapse = "/"), paste(sapply(LOCI, function(L) n_of("Head circumference", L)), collapse = "/")),
sprintf("K-T: SRS total (Z relative to SSC unaffected relatives) and SCQ lifetime total (Z relative to SSC siblings), adjusted for age within group, against BMI Z. Because the BMI effects on these traits do not vary with age, the trait is plotted against BMI rather than age. Points are coloured by BMI tail as in A-J (grey: 30th-70th centile); black line and band: linear fit with 95%% CI; diamonds: mean +/- SE in BMI bins (< -1, -1 to 0, 0 to 1, 1 to 2, > 2 SD); dashed grey line: the fit in unaffected controls. n: controls %d (SRS) and %d (SCQ); carriers %s (SRS) and %s (SCQ). Printed values: the BMI slope (SD of trait per SD of BMI; GAM y ~ s(age) + BMI) and the p-value for the difference of that slope from unaffected controls (model as in A-J). K and P also show SSC probands with ASD and no CNV (black; n = %d SRS, %d SCQ), whose slope does not differ from controls, so the association is specific to CNV carriers rather than to ASD.",
        nX("SRS", "Unaffected relatives"), nX("SCQ", "Unaffected relatives"), paste(sapply(LOCI, function(L) nX("SRS", L)), collapse = "/"), paste(sapply(LOCI, function(L) nX("SCQ", L)), collapse = "/"), nX("SRS", "SSC ASD probands"), nX("SCQ", "SSC ASD probands")),

sprintf("Heterogeneity of the age-adjusted BMI slope across the four CNV groups (within carriers, Wald test of the locus x BMI interaction): %s.", paste(sprintf("%s p = %s", HET$pheno, sapply(HET$p, pfx))[HET$pheno %in% c("Height", "Head circumference", "SRS", "SCQ")], collapse = ", ")),
"Bold: nominal p < 0.05 for that line. Carriers with and without ASD are combined (ASD status in the carrier samples: see unaffected_controls_tests.csv); ASD-stratified and ASD-adjusted estimates are in unaffected_controls_tests.csv and cognition_bmi_slopes.csv.")
writeLines(leg, file.path(HERE, "fig_main_v3_legend.txt"))
cat("fig_main_v3 written\n")
