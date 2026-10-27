#!/usr/bin/env Rscript

# =============================================================================
# MAIN FIGURE (3 rows x 2 columns)
#
#   A  pathway enrichment in MIRROR CNVs        (Sacks et al.)
#   B  pathway enrichment in ASYMMETRIC CNVs    (Sacks et al.)
#   C  pathway enrichment in height genes       (Kosmicki et al.)
#   D  negative-effect bias vs all 207 Kosmicki genes
#   E  the asymmetric-CNV genes contributing the GO:0008284 overlap, by locus
#   F  height effect of the 17 Kosmicki genes in GO:0008284, ascending
#
# ANALYSIS SET (panels A-D): the union of the pathways for which the asymmetric
# CNV test can reach P < 0.05 (>= 6 genes in the CNV universe) and those for
# which the mirror CNV test can (>= 4). The asymmetric set is a strict subset of
# the mirror set, so the union is the 102-pathway mirror set. Every volcano panel
# plots the same 102 pathways and that single set fixes the multiple-testing
# burden throughout: Bonferroni = 0.05/102.
#
# COLOUR AND LABELS (panels A-D)
#   Colour marks membership of the two CNV-enriched groups, not significance in
#   the panel being viewed:
#     blue  - the pathways enriched in mirror CNVs      (panel A, n = 3)
#     red   - the pathways enriched in asymmetric CNVs  (panel B, n = 2)
#     grey  - every other pathway
#   The two groups are disjoint. Those five pathways, and only those five, carry
#   text labels, in the same colours, in all four volcano panels. Significance
#   within a panel is read off the threshold lines. The black ring marks the
#   pathway significant in all three tests of the asymmetric analysis (B, C, D).
#
#   Panels E and F are gene-level and use their own palettes: E colours by
#   recurrent locus, F uses grey for height-decreasing and orange for
#   height-increasing. Those deliberately avoid the blue/red of A-D, which carry
#   a different meaning.
#
# No statistic is recomputed; all values are read from the master results table
# and the two annotated gene tables.
#
# Usage:  Rscript make_main_figure.R [master_csv] [gene207_csv] [cnvgenes_csv] [outdir]
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
  "Regulation of actin cytoskeleton"                                  = "Actin Cytoskeleton",
  "Protein Ubiquitination (GO:0016567)"                               = "Ubiquitination",
  "Phosphate-Containing Compound Metabolic Process (GO:0006796)"      = "Phosphate Metabolism")

# ---- analysis set -----------------------------------------------------------
att <- function(n, k, u) phyper(pmin(n, k) - 1, k, u - k, n, lower.tail = FALSE)
x[, minP_asym   := att(asymmetric_cnv_pathway_n_in_universe,
                       asymmetric_cnv_target_n, cnv_library_mapped_universe_n)]
x[, minP_mirror := att(mirror_cnv_pathway_n_in_universe,
                       mirror_cnv_target_n,     cnv_library_mapped_universe_n)]
U  <- x[minP_asym < alpha | minP_mirror < alpha]
nU <- nrow(U)
bonf <- alpha / nU

nom <- function(or, p) is.finite(or) & or > 1 & is.finite(p) & p < alpha
U[, sig_mirror := nom(mirror_cnv_odds_ratio,     mirror_cnv_p)]      # panel A
U[, sig_asym   := nom(asymmetric_cnv_odds_ratio, asymmetric_cnv_p)]  # panel B
U[, sig_kos    := nom(kosmicki_odds_ratio,       kosmicki_p)]        # panel C
U[, sig_bias   := nom(direction_odds_ratio, direction_fisher_p_one_sided)] # panel D
U[, conv_asym  := sig_asym & sig_kos & sig_bias]

stopifnot(sum(U$sig_mirror & U$sig_asym) == 0)   # the two groups must be disjoint

# ---- colour groups and labels (panels A-D) ----------------------------------
LV <- c("Enriched in mirror CNVs (A)",
        "Enriched in asymmetric CNVs (B)",
        "Other pathways")
COL_MIR   <- "#2C7FB8"
COL_ASYM  <- "#D7191C"
COL_OTHER <- "grey72"

U[, grp := fifelse(sig_mirror, LV[1], fifelse(sig_asym, LV[2], LV[3]))]
U[, grp := factor(grp, levels = LV)]
U[, lbl := fifelse(grp != LV[3], unname(SHORT[pathway]), NA_character_)]
if (any(U$grp != LV[3] & is.na(U$lbl)))
  stop("unnamed pathway selected for labelling: ",
       paste(U[grp != LV[3] & is.na(lbl)]$pathway, collapse = "; "))

grp_scale <- scale_colour_manual(
  values = c(COL_MIR, COL_ASYM, COL_OTHER),
  breaks = LV[1:2], limits = LV, drop = FALSE)

th <- theme_bw(base_size = 15) +
  theme(
    panel.grid    = element_blank(),
    panel.border  = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
    plot.title    = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 12, colour = "grey30"),
    axis.title    = element_text(size = 13),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 13.5))

# ---- volcano panel ----------------------------------------------------------
panel <- function(orcol, pcol, title, xlab) {
  d <- copy(U)
  setnames(d, c(orcol, pcol), c("OR", "P"))
  d <- d[is.finite(OR) & OR > 0 & is.finite(P) & P > 0]

  nsig  <- sum(d$OR > 1 & d$P < alpha)
  nbonf <- sum(d$OR > 1 & d$P < bonf)

  xr <- range(log2(d$OR)); yr <- range(-log10(d$P))
  xe <- diff(xr) * 0.12;   ye <- diff(yr) * 0.24

  g <- ggplot(d, aes(log2(OR), -log10(P), colour = grp)) +
    geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0,             linetype = 2, colour = "grey45", linewidth = 0.4)
  if (-log10(bonf) <= yr[2] + ye)
    g <- g + geom_hline(yintercept = -log10(bonf), linetype = 3, colour = "grey25", linewidth = 0.5)
  g +
    geom_point(data = d[grp == LV[3]], alpha = 0.75, size = 2.2) +
    geom_point(data = d[grp != LV[3]], alpha = 0.95, size = 3.6) +
    geom_point(data = d[pathway == go_label],
               aes(shape = "Significant in all three tests (B, C, D)"),
               fill = NA, colour = "black", stroke = 1.4, size = 7) +
    scale_shape_manual(name = NULL,
                       values = c("Significant in all three tests (B, C, D)" = 21)) +
    geom_text_repel(data = d[!is.na(lbl)], aes(label = lbl),
                    size = 4.0, fontface = "bold",
                    segment.colour = "grey40", segment.size = 0.4,
                    min.segment.length = 0, box.padding = 0.9,
                    point.padding = 0.6, max.overlaps = Inf,
                    force = 8, force_pull = 0.6, seed = 7, show.legend = FALSE) +
    grp_scale +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1, label = "")),
           shape  = guide_legend(override.aes = list(colour = "black", fill = NA,
                                                     stroke = 1.3, size = 5))) +
    coord_cartesian(xlim = c(xr[1] - xe, xr[2] + xe), ylim = c(yr[1] - ye, yr[2] + ye)) +
    labs(title = title,
         subtitle = sprintf("%d pathways; %d at P < 0.05, %d Bonferroni", nrow(d), nsig, nbonf),
         x = xlab, y = expression(paste(-log[10], "(P)"))) +
    th
}

pA <- panel("mirror_cnv_odds_ratio", "mirror_cnv_p",
  "A   Enrichment in mirror CNVs",
  expression(paste("Mirror-CNV enrichment, ", log[2], "(OR)"))) +
  guides(colour = "none", shape = "none")

pB <- panel("asymmetric_cnv_odds_ratio", "asymmetric_cnv_p",
  "B   Enrichment in asymmetric CNVs",
  expression(paste("Asymmetric-CNV enrichment, ", log[2], "(OR)"))) +
  guides(colour = "none", shape = "none")

pC <- panel("kosmicki_odds_ratio", "kosmicki_p",
  "C   Enrichment in Kosmicki genes",
  expression(paste("Kosmicki height-gene enrichment, ", log[2], "(OR)")))   # carries legends

pD <- panel("direction_odds_ratio", "direction_fisher_p_one_sided",
  "D   Negative-effect bias vs all 207",
  expression(paste("Negative-effect bias, ", log[2], "(OR)"))) +
  guides(colour = "none", shape = "none")

# =============================================================================
# Gene-level panels for GO:0008284
# =============================================================================
go <- U[pathway == go_label]
ov   <- strsplit(go$overlap_genes, ";")[[1]]
asym <- strsplit(go$asymmetric_cnv_overlap_genes, ";")[[1]]
# the master's gene-name strings are unreliable in general; check this row
stopifnot(length(ov)   == go$kosmicki_overlap_n,
          length(asym) == go$asymmetric_cnv_overlap_n)

# ---- E: contributing asymmetric-CNV genes, by locus -------------------------
dl <- unique(cnvg[Gene %in% asym & Dose_response_class == "asymmetric",
                  .(Gene, locus = CNV_cytoband_hg38)])
dl <- dl[, .SD[1], by = Gene]
dl[, locus := factor(locus, levels = sort(unique(locus)))]
dl <- dl[order(locus, Gene)]
dl[, Gene := factor(Gene, levels = rev(Gene))]

pE <- ggplot(dl, aes(x = locus, y = Gene, fill = locus)) +
  geom_point(shape = 21, size = 5, stroke = 0.5, colour = "grey20") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "E   Asymmetric-CNV genes contributing the GO:0008284 overlap",
    subtitle = sprintf("%d of %d asymmetric-CNV genes, %d recurrent intervals; locus-block permutation P = %.3f",
                       go$asymmetric_cnv_overlap_n, go$asymmetric_cnv_target_n,
                       uniqueN(dl$locus), go$cnv_locus_perm_p_asymmetric_greater),
    x = "Recurrent CNV locus", y = NULL) +
  th + theme(axis.text.y = element_text(size = 11, face = "italic"),
             axis.text.x = element_text(size = 10.5),
             legend.position = "none")

# ---- F: height effect of the 17 Kosmicki genes, ascending -------------------
COL_DEC <- "#4D4D4D"   # height-decreasing; deliberately not the blue of A-D
COL_INC <- "#E8850C"   # height-increasing; deliberately not the red of A-D

gk <- k207[Gene %in% ov][order(Effect_SD)]
gk[, Gene := factor(Gene, levels = rev(Gene))]   # ascending down the axis

pF <- ggplot(gk, aes(x = Effect_SD, y = Gene, colour = Effect_direction)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_segment(aes(x = 0, xend = Effect_SD, yend = Gene), linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c(Negative = COL_DEC, Positive = COL_INC), guide = "none") +
  labs(
    title = "F   Height effect of the 17 Kosmicki genes in GO:0008284",
    subtitle = sprintf("%d decreasing / %d increasing (%.0f%% vs %.0f%% background); OR %.2f, P = %.3f",
                       go$direction_negative_n, go$direction_positive_n,
                       100 * go$direction_negative_fraction,
                       100 * go$direction_background_negative_fraction,
                       go$direction_odds_ratio, go$direction_fisher_p_one_sided),
    x = "Height effect (SD per allele)", y = NULL) +
  th + theme(axis.text.y = element_text(size = 11, face = "italic"),
             legend.position = "none")

nXa <- sum(U$conv_asym)

fig <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0(
      "Three-way convergence on a single analysis set: pathway enrichment in mirror (A) and asymmetric (B) CNVs from\n",
      "Sacks et al., pathway enrichment in height genes from Kosmicki et al. (C), and enrichment of height-decreasing\n",
      "effect bias in the pathway's Kosmicki genes relative to all 207 Kosmicki genes (D), with the gene-level evidence\n",
      "behind the convergent pathway GO:0008284 (E, F)"),
    subtitle = sprintf(paste0(
      "Panels A-D plot the same %d pathways: the union of those for which the asymmetric CNV test can reach P < 0.05 (>= 6 genes in the CNV gene universe) and those for\n",
      "which the mirror CNV test can (>= 4). That single set fixes the multiple-testing burden throughout, Bonferroni = 0.05/%d = %.1e. In A-D colour marks the three\n",
      "pathways enriched in mirror CNVs (blue) and the two enriched in asymmetric CNVs (red); these groups are disjoint, and only these five carry labels, in the same\n",
      "colours in every panel so each can be followed across the tests. Significance within a panel is read off the threshold lines. The black ring marks the %d pathway\n",
      "significant in all three tests of the asymmetric analysis (B, C and D): GO:0008284. Panels E and F are gene-level and use their own palettes."),
      nU, nU, bonf, nXa),
    caption = paste0(
      "Dashed lines mark OR = 1 and nominal P = 0.05; dotted line marks the Bonferroni threshold where it falls within the plotted range. Neither CNV test shows a screen-level\n",
      "excess (asymmetric 2 hits, mirror 3; chance expects 2.4 for the asymmetric test, which only 48 of the 102 pathways can reach). The contrast between the two CNV classes is a\n",
      "specificity check rather than a formal test. In E, genes within one recurrent locus travel on the same CNV and are not independent, which is why the locus-block permutation,\n",
      "not the gene-level Fisher test, is the statistic of record for the CNV layer."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 17, lineheight = 1.2),
      plot.subtitle = element_text(size = 12.5, colour = "grey25", lineheight = 1.3),
      plot.caption  = element_text(size = 11, colour = "grey45", lineheight = 1.25))
  ) & theme(legend.position = "bottom")

ggsave(file.path(outdir, "Figure_1_three_way_convergence.png"), fig,
       width = 15.5, height = 19.5, dpi = 400)
ggsave(file.path(outdir, "Figure_1_three_way_convergence.pdf"), fig,
       width = 15.5, height = 19.5, device = cairo_pdf)

cat("Figure 1 (main) written to ", outdir, "\n", sep = "")
cat("Analysis set n = ", nU, "; Bonferroni = ", signif(bonf, 3), "\n", sep = "")
cat("A mirror P<0.05: ", sum(U$sig_mirror), "   B asymmetric P<0.05: ", sum(U$sig_asym),
    "   convergent (B,C,D): ", nXa, "\n", sep = "")
cat("E: ", nrow(dl), " asymmetric-CNV genes;  F: ", nrow(gk),
    " Kosmicki genes (ascending)\n", sep = "")
