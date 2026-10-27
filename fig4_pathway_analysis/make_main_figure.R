#!/usr/bin/env Rscript

# =============================================================================
# MAIN FIGURE (Figure 4), five panels  -- JS 2026-10-21 revision of the 6-panel version
#   (the 6-panel script is kept as make_main_figure_v1_6panel.R)
#
#   A  asymmetric vs mirror CNVs: pathway enrichment contrast (two-sided)
#   B  pathway enrichment in height genes       (Kosmicki et al., 207 exome genes)
#   C  height-decreasing bias of the pathway's Kosmicki genes vs all mapped Kosmicki genes
#   D  the asymmetric-CNV genes contributing the GO:0008284 overlap, by locus
#   E  height effect of the 17 Kosmicki genes in GO:0008284, ascending
#
# WHY ONE CNV PANEL. The CNV gene universe contains only genes of the mirror and
# asymmetric CNV classes, so the former panels A (mirror enrichment) and B
# (asymmetric enrichment) were near-complements: enrichment in one class is
# depletion in the other (log OR correlation 0.82 over the analysis set; not exact
# because 13 genes in nested 22q11.2 intervals belong to both classes). Panel A
# now shows the single contrast. x = log2 OR of the pathway's CNV-universe genes
# being asymmetric-class rather than not (the pipeline's asymmetric_cnv_odds_ratio;
# OR > 1 = enriched in asymmetric CNVs, OR < 1 = enriched in mirror CNVs). y = the
# TWO-SIDED Fisher exact P of that same 2x2 table, computed here from the counts
# stored in the master table (the only statistic this script computes; the
# pipeline stored one-sided hypergeometric P values for each class separately).
# Labels sit on both tails: red = asymmetric-enriched (OR > 1, P < 0.05),
# blue = mirror-enriched (OR < 1, P < 0.05).
#
# ANALYSIS SET (panels A-C): the union of the pathways for which the one-sided
# asymmetric test can reach P < 0.05 (>= 6 genes in the CNV universe) and those
# for which the mirror test can (>= 4); the asymmetric set is a strict subset of
# the mirror set, so the union is the 102-pathway mirror set. That single set
# fixes the multiple-testing burden: Bonferroni = 0.05/102, BH FDR within it.
# The black ring marks the pathway significant in A, B and C.
#
# Panels D and E are gene-level and use their own palettes (D by recurrent locus;
# E grey = height-decreasing, orange = height-increasing).
#
# Usage:  Rscript make_main_figure.R [master_csv] [gene207_csv] [cnvgenes_csv] [outdir]
# Output: <outdir>/Figure_4_convergence.{png,pdf}, <outdir>/panelA_asym_vs_mirror_two_sided.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "kosmicki_pathway_cnv_master_results.csv"
g207in <- if (length(args) >= 2) args[2] else "kosmicki_207_genes_annotated.csv"
cnvin  <- if (length(args) >= 3) args[3] else "cnv_genes_annotated.csv"
outdir <- if (length(args) >= 4) args[4] else "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x    <- fread(input)
k207 <- fread(g207in)
cnvg <- fread(cnvin)

alpha    <- 0.05
go_label <- "Positive Regulation Of Cell Population Proliferation (GO:0008284)"

SHORT <- c(
  "Positive Regulation Of Cell Population Proliferation (GO:0008284)" = "Pos Reg of Cell Proliferation",
  "Regulation Of Nucleic Acid-Templated Transcription (GO:1903506)"   = "Transcription Reg",
  "Negative Regulation Of Nucleic Acid-Templated Transcription (GO:1903507)" = "Neg Transcription Reg",
  "Regulation of actin cytoskeleton"                                  = "Actin Cytoskeleton",
  "Protein Ubiquitination (GO:0016567)"                               = "Ubiquitination",
  "Phosphate-Containing Compound Metabolic Process (GO:0006796)"      = "Phosphate Metabolism")
short_name <- function(p) { s <- unname(SHORT[p]); ifelse(is.na(s), sub(" \\(GO:[0-9]+\\)$", "", p), s) }

# ---- analysis set -----------------------------------------------------------
att <- function(n, k, u) phyper(pmin(n, k) - 1, k, u - k, n, lower.tail = FALSE)
x[, minP_asym   := att(asymmetric_cnv_pathway_n_in_universe,
                       asymmetric_cnv_target_n, cnv_library_mapped_universe_n)]
x[, minP_mirror := att(mirror_cnv_pathway_n_in_universe,
                       mirror_cnv_target_n,     cnv_library_mapped_universe_n)]
U  <- x[minP_asym < alpha | minP_mirror < alpha]
nU <- nrow(U)
bonf <- alpha / nU

# ---- panel A statistic: two-sided Fisher exact P on the asymmetric test's 2x2 ----
#   rows: pathway genes in the CNV universe / other CNV-universe genes
#   cols: asymmetric class / not asymmetric (= mirror-only)
f2 <- function(a, b, c, d) fisher.test(matrix(c(a, b, c, d), 2))$p.value
U[, a_ := asymmetric_cnv_overlap_n]
U[, b_ := asymmetric_cnv_pathway_n_in_universe - asymmetric_cnv_overlap_n]
U[, c_ := asymmetric_cnv_target_n - asymmetric_cnv_overlap_n]
U[, d_ := cnv_library_mapped_universe_n - asymmetric_cnv_target_n - b_]
U[, asym_vs_mirror_p_two_sided := mapply(f2, a_, b_, c_, d_)]

nom <- function(or, p) is.finite(or) & or > 1 & is.finite(p) & p < alpha
U[, sig_asym   := nom(asymmetric_cnv_odds_ratio,       asym_vs_mirror_p_two_sided)]   # red
U[, sig_mirror := is.finite(asymmetric_cnv_odds_ratio) & asymmetric_cnv_odds_ratio < 1 &
                  asym_vs_mirror_p_two_sided < alpha]                                  # blue
U[, sig_kos    := nom(kosmicki_odds_ratio,  kosmicki_p)]                                # panel B
U[, sig_bias   := nom(direction_odds_ratio, direction_fisher_p_one_sided)]             # panel C
U[, conv_asym  := sig_asym & sig_kos & sig_bias]

fwrite(U[order(asym_vs_mirror_p_two_sided),
         .(library, pathway, n_cnv_universe_genes = asymmetric_cnv_pathway_n_in_universe,
           n_asymmetric = asymmetric_cnv_overlap_n, n_mirror = mirror_cnv_overlap_n,
           asym_vs_mirror_odds_ratio = asymmetric_cnv_odds_ratio,
           asym_one_sided_p = asymmetric_cnv_p, mirror_one_sided_p = mirror_cnv_p,
           asym_vs_mirror_p_two_sided, sig_asym, sig_mirror, sig_kos, sig_bias, conv_asym)],
       file.path(outdir, "panelA_asym_vs_mirror_two_sided.csv"))

# ---- colour groups and labels (panels A-C) ----------------------------------
LV <- c("Enriched in asymmetric CNVs (P < 0.05)",
        "Enriched in mirror CNVs (P < 0.05)",
        "Other pathways")
COL_ASYM  <- "#D7191C"
COL_MIR   <- "#2C7FB8"
COL_OTHER <- "grey72"

U[, grp := fifelse(sig_asym, LV[1], fifelse(sig_mirror, LV[2], LV[3]))]
U[, grp := factor(grp, levels = LV)]
U[, lbl := fifelse(grp != LV[3], short_name(pathway), NA_character_)]

grp_scale <- scale_colour_manual(
  values = c(COL_ASYM, COL_MIR, COL_OTHER),
  breaks = LV[1:2], limits = LV, drop = FALSE)

th <- theme_bw(base_size = 18) +   # JS 2026-10-22: larger type throughout
  theme(
    panel.grid    = element_blank(),
    panel.border  = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13.5, colour = "grey30"),
    axis.title    = element_text(size = 18),
    axis.text     = element_text(size = 15, colour = "grey15"),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 13.5))

RING <- "Significant in A, B and C"

# ---- volcano panel ----------------------------------------------------------
panel <- function(orcol, pcol, title, subtitle, xlab, pad_bottom = 0.08, pin = NULL) {   # JS 2026-10-22: tight y padding
  # pin: optional named list of c(dx, dy) nudges (data units) that place a given label away from the points,
  #      e.g. pin = list("Actin Cytoskeleton" = c(1.6, -0.4)); other labels are placed by ggrepel as before
  d <- copy(U)
  setnames(d, c(orcol, pcol), c("OR", "P"))
  d <- d[is.finite(OR) & OR > 0 & is.finite(P) & P > 0]
  lab_free <- d[!is.na(lbl) & !(lbl %in% names(pin))]
  lab_pin  <- d[!is.na(lbl) &   lbl %in% names(pin)]

  xr <- range(log2(d$OR)); yr <- range(-log10(d$P))
  xe <- diff(xr) * 0.12;   ye <- diff(yr) * 0.10;   yb <- diff(yr) * pad_bottom

  g <- ggplot(d, aes(log2(OR), -log10(P), colour = grp)) +
    geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0,             linetype = 2, colour = "grey45", linewidth = 0.4)
  if (-log10(bonf) <= yr[2] + ye)
    g <- g + geom_hline(yintercept = -log10(bonf), linetype = 3, colour = "grey25", linewidth = 0.5)
  g +
    geom_point(data = d[grp == LV[3]], alpha = 0.75, size = 2.6) +
    geom_point(data = d[grp != LV[3]], alpha = 0.95, size = 4.2) +
    geom_point(data = d[pathway == go_label], aes(shape = RING),
               fill = NA, colour = "black", stroke = 1.5, size = 8) +
    scale_shape_manual(name = NULL, values = setNames(21, RING)) +
    geom_text_repel(data = lab_free, aes(label = lbl),
                    size = 5.2, fontface = "bold",
                    segment.colour = "grey40", segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.9,
                    point.padding = 0.6, max.overlaps = Inf,
                    force = 8, force_pull = 0.6, seed = 7, show.legend = FALSE) +
    (if (nrow(lab_pin)) geom_text_repel(data = lab_pin, aes(label = lbl),
                    nudge_x = vapply(lab_pin$lbl, function(l) pin[[l]][1], numeric(1)),
                    nudge_y = vapply(lab_pin$lbl, function(l) pin[[l]][2], numeric(1)),
                    size = 5.2, fontface = "bold", segment.colour = "grey40", segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.3, point.padding = 0.4, max.overlaps = Inf,
                    force = 0, force_pull = 0, direction = "both", seed = 7, show.legend = FALSE) else NULL) +
    grp_scale +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1, label = "")),
           shape  = guide_legend(override.aes = list(colour = "black", fill = NA,
                                                     stroke = 1.3, size = 5))) +
    coord_cartesian(xlim = c(xr[1] - xe, xr[2] + xe), ylim = c(yr[1] - yb, yr[2] + ye)) +
    labs(title = title, subtitle = subtitle,
         x = xlab, y = expression(paste(-log[10], "(P)"))) +
    th
}

nA_asym <- sum(U$sig_asym); nA_mir <- sum(U$sig_mirror)
nB <- sum(U$sig_kos); bB <- sum(U$sig_kos & U$kosmicki_p < bonf)
nC <- sum(U$sig_bias); bC <- sum(U$sig_bias & U$direction_fisher_p_one_sided < bonf)

pA <- panel("asymmetric_cnv_odds_ratio", "asym_vs_mirror_p_two_sided",
  "A   Asymmetric vs mirror CNVs",
  sprintf("%d pathways; two-sided P < 0.05: %d asymmetric-, %d mirror-enriched", nU, nA_asym, nA_mir),
  expression(paste("Pathway genes in asymmetric vs mirror CNVs, ", log[2], "(OR)")), pad_bottom = 0.22) +
  # JS 2026-10-22: direction of the contrast spelled out inside the panel; the colour key is collected into its own
  # cell (guide_area) above panel D rather than drawn inside A
  annotate("text", x = -Inf, y = -Inf, label = "←  enriched in mirror CNVs", hjust = -0.05, vjust = -0.7,
           colour = COL_MIR, fontface = "bold", size = 5.4) +
  annotate("text", x = Inf, y = -Inf, label = "enriched in asymmetric CNVs  →", hjust = 1.05, vjust = -0.7,
           colour = COL_ASYM, fontface = "bold", size = 5.4) +
  theme(legend.direction = "vertical", legend.box = "horizontal", legend.box.just = "top",
        legend.text = element_text(size = 15), legend.key.size = unit(18, "pt"),
        legend.spacing.x = unit(18, "pt"), legend.margin = margin(0, 0, 0, 0))

pB <- panel("kosmicki_odds_ratio", "kosmicki_p",
  "B   Enrichment in height genes from Kosmicki et al.",
  sprintf("%d pathways; %d at P < 0.05, %d Bonferroni", nU, nB, bB),
  expression(paste("Kosmicki height-gene enrichment, ", log[2], "(OR)")),
  pin = list("Actin Cytoskeleton" = c(1.7, -0.35))) +   # label to the lower right, clear of the point cloud
  guides(colour = "none", shape = "none")

pC <- panel("direction_odds_ratio", "direction_fisher_p_one_sided",
  "C   Bias of pathway genes toward short stature in Kosmicki et al.",
  sprintf("%d pathways; %d at P < 0.05, %d Bonferroni", nU, nC, bC),
  expression(paste("Height-decreasing bias, ", log[2], "(OR)")),
  pin = list("Actin Cytoskeleton" = c(0.85, -0.5))) +   # label to the lower right, clear of the point cloud
  guides(colour = "none", shape = "none")

# =============================================================================
# Gene-level panels for GO:0008284
# =============================================================================
go <- U[pathway == go_label]
ov   <- strsplit(go$overlap_genes, ";")[[1]]
asym <- strsplit(go$asymmetric_cnv_overlap_genes, ";")[[1]]
stopifnot(length(ov)   == go$kosmicki_overlap_n,
          length(asym) == go$asymmetric_cnv_overlap_n)

# ---- D: contributing asymmetric-CNV genes, by locus -------------------------
dl <- unique(cnvg[Gene %in% asym & Dose_response_class == "asymmetric",
                  .(Gene, locus = CNV_cytoband_hg38)])
dl <- dl[, .SD[1], by = Gene]
dl[, locus := factor(locus, levels = sort(unique(locus)))]
dl <- dl[order(locus, Gene)]
dl[, Gene := factor(Gene, levels = rev(Gene))]

pD <- ggplot(dl, aes(x = locus, y = Gene, fill = locus)) +
  geom_point(shape = 21, size = 4.2, stroke = 0.5, colour = "grey20") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_x_discrete(labels = function(l) sub("-", "-\n", l)) +   # two-line locus labels
  labs(
    title = "D   Asymmetric-CNV genes in GO:0008284",
    subtitle = sprintf("%d of %d asymmetric-CNV genes in %d loci; locus-block permutation P = %.3f",
                       go$asymmetric_cnv_overlap_n, go$asymmetric_cnv_target_n,
                       uniqueN(dl$locus), go$cnv_locus_perm_p_asymmetric_greater),
    x = "Recurrent CNV locus", y = NULL) +
  th + theme(axis.text.y = element_text(size = 11, face = "italic"),
             axis.text.x = element_text(size = 12),
             plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 11, colour = "grey30"),
             axis.title = element_text(size = 14), legend.position = "none")

# ---- E: height effect of the 17 Kosmicki genes, ascending -------------------
COL_DEC <- "#4D4D4D"
COL_INC <- "#E8850C"

gk <- k207[Gene %in% ov][order(Effect_SD)]
gk[, Gene := factor(Gene, levels = rev(Gene))]

pE <- ggplot(gk, aes(x = Effect_SD, y = Gene, colour = Effect_direction)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_segment(aes(x = 0, xend = Effect_SD, yend = Gene), linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(Negative = COL_DEC, Positive = COL_INC), guide = "none") +
  labs(
    title = sprintf("E   Height effect of the %d Kosmicki genes in GO:0008284", nrow(gk)),
    subtitle = sprintf("%d decreasing / %d increasing (%.0f%% vs %.0f%% background); OR %.2f, P = %.3f",
                       go$direction_negative_n, go$direction_positive_n,
                       100 * go$direction_negative_fraction,
                       100 * go$direction_background_negative_fraction,
                       go$direction_odds_ratio, go$direction_fisher_p_one_sided),
    x = "Height effect (SD per allele)", y = NULL) +
  th + theme(axis.text.y = element_text(size = 9.5, face = "italic"),
             plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 11, colour = "grey30"),
             axis.title = element_text(size = 14), legend.position = "none")

nXa <- sum(U$conv_asym)

# ---- layout (JS 2026-10-22) ---------------------------------------------------
# A | B on top; C below A at the same height, with the colour key (collected from A) under C; D and E stacked on the right.
# The figure legend text is NOT drawn in the image; it lives in results/Figure_4_legend.txt for the manuscript.
# D and E are nested into one element (wrap_elements) so their title/axis gutters do not stretch C's row:
# the outer grid is then two columns with rows 1 and 2 of equal height, making A, B and C identical in size.
right <- wrap_elements(full = (pD / pE) + plot_layout(heights = c(0.36, 0.64)))
design <- "AB\nCR\nKR"   # row 1: A | B; row 2: C | D+E; row 3: key under C (R spans rows 2-3)
fig <- pA + pB + pC + guide_area() + right +   # plots are matched to design letters alphabetically (A B C K R)
  plot_layout(design = design, heights = c(1, 1, 0.14), guides = "collect") &
  theme(legend.direction = "vertical", legend.box = "horizontal", legend.box.just = "top", legend.spacing.x = unit(24, "pt"),
        legend.text = element_text(size = 15), legend.key.size = unit(18, "pt"), legend.margin = margin(4, 0, 0, 0))   # applies to the collected key

ggsave(file.path(outdir, "Figure_4_convergence.png"), fig, width = 16, height = 12.5, dpi = 400)
ggsave(file.path(outdir, "Figure_4_convergence.pdf"), fig, width = 16, height = 12.5, device = cairo_pdf)
ggsave(file.path(outdir, "Figure_4_convergence_preview.png"), fig, width = 16, height = 12.5, dpi = 200)   # screen copy

fig_old_annotation <- function() {   # kept for reference; not called
  plot_annotation(
    title = paste0(
      "Convergence of pathway enrichment in asymmetric CNVs (A) with enrichment in height genes from Kosmicki et al. (B)\n",
      "and a height-decreasing effect bias among the pathway's Kosmicki genes (C); gene-level evidence for GO:0008284 (D, E)"),
    subtitle = sprintf(paste0(
      "Panels A-C plot the same %d pathways (those with enough genes in the CNV gene universe for a CNV test to reach P < 0.05); Bonferroni = 0.05/%d = %.1e.\n",
      "A contrasts the two CNV classes directly: OR > 1 = the pathway's CNV genes are over-represented in asymmetric CNVs, OR < 1 = in mirror CNVs (two-sided P).\n",
      "Red = asymmetric-enriched (n = %d), blue = mirror-enriched (n = %d); these pathways carry labels in the same colours in A-C. The black ring marks the %d pathway\n",
      "significant in A, B and C: GO:0008284. Panels D and E are gene-level and use their own palettes."),
      nU, nU, bonf, nA_asym, nA_mir, nXa),
    caption = paste0(
      "Dashed lines mark OR = 1 and nominal P = 0.05; dotted line marks the Bonferroni threshold where it falls within the plotted range. The CNV gene universe holds only mirror- and\n",
      "asymmetric-class genes, so enrichment in one class is depletion in the other; panel A therefore replaces separate mirror and asymmetric panels. In D, genes within one recurrent\n",
      "locus travel on the same CNV and are not independent, which is why the locus-block permutation, not the gene-level test, is the statistic of record for the CNV layer."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 17, lineheight = 1.2),
      plot.subtitle = element_text(size = 12.5, colour = "grey25", lineheight = 1.3),
      plot.caption  = element_text(size = 11, colour = "grey45", lineheight = 1.25))
  )
}

cat("Figure 4 written to ", outdir, "\n", sep = "")
cat("Analysis set n = ", nU, "; Bonferroni = ", signif(bonf, 3), "\n", sep = "")
cat("A two-sided P<0.05: asymmetric-enriched ", nA_asym, " (", paste(U[sig_asym == TRUE, lbl], collapse = "; "),
    "); mirror-enriched ", nA_mir, " (", paste(U[sig_mirror == TRUE, lbl], collapse = "; "), ")\n", sep = "")
cat("B Kosmicki P<0.05: ", nB, " (", bB, " Bonferroni);  C bias P<0.05: ", nC, ";  convergent (A,B,C): ", nXa, "\n", sep = "")
print(U[pathway == go_label, .(n_cnv_universe = asymmetric_cnv_pathway_n_in_universe, n_asym = asymmetric_cnv_overlap_n, n_mirror = mirror_cnv_overlap_n,
                               OR_A = round(asymmetric_cnv_odds_ratio, 2), P_A_two_sided = signif(asym_vs_mirror_p_two_sided, 3), P_asym_one_sided = signif(asymmetric_cnv_p, 3),
                               perm_P = signif(cnv_locus_perm_p_asymmetric_greater, 3), kos = paste0(kosmicki_overlap_n, "/", kosmicki_mapped_n),
                               OR_B = round(kosmicki_odds_ratio, 2), P_B = signif(kosmicki_p, 3), q_B = NA,
                               neg_pos = paste0(direction_negative_n, "/", direction_positive_n), OR_C = round(direction_odds_ratio, 2), P_C = signif(direction_fisher_p_one_sided, 3))])
cat("universe: ", U$cnv_library_mapped_universe_n[1], " CNV genes mapped to the GO library (", U$asymmetric_cnv_target_n[1], " asymmetric, ",
    U$mirror_cnv_target_n[1], " mirror); Kosmicki genes mapped: ", U$kosmicki_mapped_n[1], "\n", sep = "")
