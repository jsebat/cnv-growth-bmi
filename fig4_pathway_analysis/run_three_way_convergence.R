#!/usr/bin/env Rscript

# =============================================================================
# THREE-WAY CONVERGENCE ANALYSIS  (primary analysis script)
#
# Four tests, all on one analysis set:
#   A  pathway enrichment in mirror CNVs        (Sacks et al.)
#   B  pathway enrichment in asymmetric CNVs    (Sacks et al.)
#   C  pathway enrichment in height genes       (Kosmicki et al., 207 exome genes)
#   D  negative-effect bias among the pathway's Kosmicki genes, relative to all 207
#
# ANALYSIS SET
#   The union of the pathways for which the asymmetric CNV test can reach
#   P < 0.05 and those for which the mirror CNV test can. For each pathway the
#   best-case attainable P (all of its genes in the CNV universe belonging to the
#   target class) is computed from the hypergeometric distribution; a pathway is
#   admitted if that is below 0.05 for either class. This works out to >= 6 genes
#   in the CNV universe for the asymmetric test and >= 4 for the mirror test.
#   The asymmetric set is a strict subset of the mirror set, so the union is the
#   102-pathway mirror set.
#
#   That single set is used for all four tests and fixes the multiple-testing
#   burden throughout: Bonferroni = 0.05/102, with BH FDR computed within it.
#   Tests C and D do not involve CNVs, so they are computed once and shared by
#   the two convergence analyses; only the CNV test differs between them.
#
# CONVERGENCE
#   asymmetric : significant in B, C and D
#   mirror     : significant in A, C and D
#
# No pathway statistic is recomputed. All enrichment, direction and permutation
# values are read unchanged from kosmicki_pathway_cnv_master_results.csv.
#
# Usage:  Rscript run_three_way_convergence.R [master_csv] [outdir]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "kosmicki_pathway_cnv_master_results.csv"
outdir <- if (length(args) >= 2) args[2] else "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- fread(input)
alpha    <- 0.05
go_label <- "Positive Regulation Of Cell Population Proliferation (GO:0008284)"

log <- character(0)
say <- function(...) { s <- paste0(...); log <<- c(log, s); cat(s, "\n", sep = "") }

# -----------------------------------------------------------------------------
# Analysis set
# -----------------------------------------------------------------------------
att <- function(n, k, u) phyper(pmin(n, k) - 1, k, u - k, n, lower.tail = FALSE)
x[, minP_asym   := att(asymmetric_cnv_pathway_n_in_universe,
                       asymmetric_cnv_target_n, cnv_library_mapped_universe_n)]
x[, minP_mirror := att(mirror_cnv_pathway_n_in_universe,
                       mirror_cnv_target_n,     cnv_library_mapped_universe_n)]
x[, in_set := (minP_asym < alpha) | (minP_mirror < alpha)]

S  <- x[in_set == TRUE]
nS <- nrow(S)
bonf <- alpha / nS

min_asym   <- min(S[minP_asym   < alpha]$asymmetric_cnv_pathway_n_in_universe)
min_mirror <- min(S[minP_mirror < alpha]$mirror_cnv_pathway_n_in_universe)

# -----------------------------------------------------------------------------
# The four tests
# -----------------------------------------------------------------------------
nom <- function(or, p) is.finite(or) & or > 1 & is.finite(p) & p < alpha
bft <- function(or, p) is.finite(or) & or > 1 & is.finite(p) & p < bonf

S[, sig_A_mirror     := nom(mirror_cnv_odds_ratio,     mirror_cnv_p)]
S[, sig_B_asymmetric := nom(asymmetric_cnv_odds_ratio, asymmetric_cnv_p)]
S[, sig_C_kosmicki   := nom(kosmicki_odds_ratio,       kosmicki_p)]
S[, sig_D_bias       := nom(direction_odds_ratio, direction_fisher_p_one_sided)]

S[, bonf_A := bft(mirror_cnv_odds_ratio,     mirror_cnv_p)]
S[, bonf_B := bft(asymmetric_cnv_odds_ratio, asymmetric_cnv_p)]
S[, bonf_C := bft(kosmicki_odds_ratio,       kosmicki_p)]
S[, bonf_D := bft(direction_odds_ratio, direction_fisher_p_one_sided)]

S[, q_A := p.adjust(mirror_cnv_p,                 method = "BH")]
S[, q_B := p.adjust(asymmetric_cnv_p,             method = "BH")]
S[, q_C := p.adjust(kosmicki_p,                   method = "BH")]
S[, q_D := p.adjust(direction_fisher_p_one_sided, method = "BH")]

S[, convergent_asymmetric := sig_B_asymmetric & sig_C_kosmicki & sig_D_bias]
S[, convergent_mirror     := sig_A_mirror     & sig_C_kosmicki & sig_D_bias]

cnt <- function(v) sum(v)
nA <- cnt(S$sig_A_mirror); nB <- cnt(S$sig_B_asymmetric)
nC <- cnt(S$sig_C_kosmicki); nD <- cnt(S$sig_D_bias)
bA <- cnt(S$bonf_A); bB <- cnt(S$bonf_B); bC <- cnt(S$bonf_C); bD <- cnt(S$bonf_D)
nXa <- cnt(S$convergent_asymmetric); nXm <- cnt(S$convergent_mirror)

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
say("=============================================================")
say("THREE-WAY CONVERGENCE ANALYSIS")
say("No statistic recomputed; all values read from the master results table.")
say("=============================================================")
say("")
say("-- PATHWAY SET --")
say("  GMT terms retained at 5-500 unique gene symbols")
say("  (min_set_size=5, max_set_size=500 in kosmicki_pathway_cnv_pipeline.R):")
say("    GO Biological Process 2023 : 5407 terms -> ",
    sum(x$library == "GO_Biological_Process_2023"), " retained")
say("    KEGG 2021 Human            :  320 terms -> ",
    sum(x$library == "KEGG_2021_Human"), " retained")
say("    Pathways with statistics   : ", nrow(x))
say("")
say("  ANALYSIS SET = union of the two CNV attainability sets.")
say("    asymmetric test attainable (>= ", min_asym, " genes in the CNV universe) : ",
    sum(x$minP_asym < alpha))
say("    mirror test attainable     (>= ", min_mirror, " genes in the CNV universe) : ",
    sum(x$minP_mirror < alpha))
say("    asymmetric set is a subset of the mirror set : ",
    all(x$minP_asym >= alpha | x$minP_mirror < alpha))
say("    UNION = ", nS, " pathways   (GO:BP ", sum(S$library == "GO_Biological_Process_2023"),
    ", KEGG ", sum(S$library == "KEGG_2021_Human"), ")")
say("")
say("  One set for all four tests. Bonferroni = 0.05/", nS, " = ", signif(bonf, 3),
    "; BH FDR within the set.")
say("")
say("-- THE FOUR TESTS --")
say(sprintf("  %-38s %8s %12s", "", "P<0.05", "Bonferroni"))
say(sprintf("  %-38s %8d %12d", "A  mirror-CNV enrichment", nA, bA))
say(sprintf("  %-38s %8d %12d", "B  asymmetric-CNV enrichment", nB, bB))
say(sprintf("  %-38s %8d %12d", "C  Kosmicki height-gene enrichment", nC, bC))
say(sprintf("  %-38s %8d %12d", "D  negative-effect bias", nD, bD))
say("")
say("-- CONVERGENCE --")
say("  asymmetric (B & C & D) : ", nXa)
say("  mirror     (A & C & D) : ", nXm)
g <- S[convergent_asymmetric == TRUE]
for (i in seq_len(nrow(g))) {
  say("")
  say("  ", g$pathway[i], "   [", g$library[i], "]")
  say("    A  mirror CNV     : ", g$mirror_cnv_overlap_n[i], " genes, OR ",
      round(g$mirror_cnv_odds_ratio[i], 2), ", P ", signif(g$mirror_cnv_p[i], 3),
      "   (not significant)")
  say("    B  asymmetric CNV : ", g$asymmetric_cnv_overlap_n[i], " of ",
      g$asymmetric_cnv_pathway_n_in_universe[i], " universe genes, OR ",
      round(g$asymmetric_cnv_odds_ratio[i], 2), ", P ", signif(g$asymmetric_cnv_p[i], 3),
      ", q ", signif(g$q_B[i], 3))
  say("       locus-block permutation P ", signif(g$cnv_locus_perm_p_asymmetric_greater[i], 3),
      "   <-- statistic of record for the CNV layer")
  say("    C  Kosmicki       : ", g$kosmicki_overlap_n[i], "/", g$kosmicki_mapped_n[i],
      " genes, OR ", round(g$kosmicki_odds_ratio[i], 2), ", P ", signif(g$kosmicki_p[i], 3),
      ", q ", signif(g$q_C[i], 3), ifelse(g$bonf_C[i], "   *** BONFERRONI SIGNIFICANT ***", ""))
  say("    D  negative bias  : ", g$direction_negative_n[i], " neg / ",
      g$direction_positive_n[i], " pos, OR ", round(g$direction_odds_ratio[i], 2),
      ", P ", signif(g$direction_fisher_p_one_sided[i], 3), ", q ", signif(g$q_D[i], 3))
}

say("")
say("-- RANK OF THE CONVERGENT PATHWAY WITHIN THE ANALYSIS SET --")
rlo <- function(v, val) sum(v < val, na.rm = TRUE) + 1
rhi <- function(v, val) sum(v > val, na.rm = TRUE) + 1
rk <- data.table(
  test = c("A mirror-CNV", "B asymmetric-CNV", "C Kosmicki", "D negative-effect bias"),
  P = c(g$mirror_cnv_p[1], g$asymmetric_cnv_p[1], g$kosmicki_p[1],
        g$direction_fisher_p_one_sided[1]),
  rank_by_P = c(rlo(S$mirror_cnv_p, g$mirror_cnv_p[1]),
                rlo(S$asymmetric_cnv_p, g$asymmetric_cnv_p[1]),
                rlo(S$kosmicki_p, g$kosmicki_p[1]),
                rlo(S$direction_fisher_p_one_sided, g$direction_fisher_p_one_sided[1])),
  odds_ratio = c(g$mirror_cnv_odds_ratio[1], g$asymmetric_cnv_odds_ratio[1],
                 g$kosmicki_odds_ratio[1], g$direction_odds_ratio[1]),
  rank_by_OR = c(rhi(S$mirror_cnv_odds_ratio, g$mirror_cnv_odds_ratio[1]),
                 rhi(S$asymmetric_cnv_odds_ratio, g$asymmetric_cnv_odds_ratio[1]),
                 rhi(S$kosmicki_odds_ratio, g$kosmicki_odds_ratio[1]),
                 rhi(S$direction_odds_ratio, g$direction_odds_ratio[1])),
  n_analysis_set = nS)
for (i in seq_len(nrow(rk)))
  say(sprintf("  %-24s P %-9.3g rank %3d/%d    OR %6.2f  rank %3d/%d",
              rk$test[i], rk$P[i], rk$rank_by_P[i], nS,
              rk$odds_ratio[i], rk$rank_by_OR[i], nS))
say("=============================================================")

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
keycols <- c(
  "library", "pathway", "pathway_size",
  "mirror_cnv_pathway_n_in_universe", "mirror_cnv_overlap_n",
  "mirror_cnv_odds_ratio", "mirror_cnv_p", "q_A", "sig_A_mirror", "bonf_A",
  "asymmetric_cnv_pathway_n_in_universe", "asymmetric_cnv_overlap_n",
  "asymmetric_cnv_odds_ratio", "asymmetric_cnv_p", "q_B", "sig_B_asymmetric", "bonf_B",
  "cnv_locus_perm_p_asymmetric_greater",
  "kosmicki_overlap_n", "kosmicki_mapped_n", "kosmicki_odds_ratio",
  "kosmicki_p", "q_C", "sig_C_kosmicki", "bonf_C",
  "direction_overlap_n", "direction_negative_n", "direction_positive_n",
  "direction_odds_ratio", "direction_fisher_p_one_sided", "q_D", "sig_D_bias", "bonf_D",
  "convergent_asymmetric", "convergent_mirror")

fwrite(S[order(asymmetric_cnv_p), ..keycols],
       file.path(outdir, "analysis_set_102_pathways.csv"))
fwrite(S[convergent_asymmetric == TRUE | convergent_mirror == TRUE, ..keycols],
       file.path(outdir, "convergent_pathways.csv"))
fwrite(S[sig_A_mirror     == TRUE, ..keycols][order(mirror_cnv_p)],
       file.path(outdir, "testA_mirror_significant.csv"))
fwrite(S[sig_B_asymmetric == TRUE, ..keycols][order(asymmetric_cnv_p)],
       file.path(outdir, "testB_asymmetric_significant.csv"))
fwrite(S[sig_C_kosmicki   == TRUE, ..keycols][order(kosmicki_p)],
       file.path(outdir, "testC_kosmicki_significant.csv"))
fwrite(S[sig_D_bias       == TRUE, ..keycols][order(direction_fisher_p_one_sided)],
       file.path(outdir, "testD_bias_significant.csv"))
fwrite(rk, file.path(outdir, "convergent_pathway_ranks.csv"))

univ <- x[library == "GO_Biological_Process_2023", cnv_library_mapped_universe_n][1]
attain <- data.table(genes_in_cnv_universe = 1:14)
attain[, best_case_P_asymmetric := att(genes_in_cnv_universe,
         x[library == "GO_Biological_Process_2023", asymmetric_cnv_target_n][1], univ)]
attain[, best_case_P_mirror := att(genes_in_cnv_universe,
         x[library == "GO_Biological_Process_2023", mirror_cnv_target_n][1], univ)]
fwrite(attain, file.path(outdir, "analysis_set_attainability.csv"))

fwrite(data.table(
  quantity = c("pathways with statistics", "asymmetric-attainable", "mirror-attainable",
               "ANALYSIS SET (union)", "Bonferroni threshold",
               "A mirror P<0.05", "A Bonferroni",
               "B asymmetric P<0.05", "B Bonferroni",
               "C Kosmicki P<0.05", "C Bonferroni",
               "D bias P<0.05", "D Bonferroni",
               "convergent asymmetric (B&C&D)", "convergent mirror (A&C&D)"),
  value = c(nrow(x), sum(x$minP_asym < alpha), sum(x$minP_mirror < alpha),
            nS, signif(bonf, 4), nA, bA, nB, bB, nC, bC, nD, bD, nXa, nXm)),
  file.path(outdir, "analysis_summary.csv"))

writeLines(log, file.path(outdir, "00_analysis_log.txt"))
cat("\nWrote outputs to ", outdir, "\n", sep = "")
