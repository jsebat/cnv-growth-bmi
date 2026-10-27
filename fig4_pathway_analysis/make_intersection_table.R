#!/usr/bin/env Rscript

# =============================================================================
# Intersection table for Figure 1: the pathways at each exclusive intersection
#
# of the three tests, per analysis (asymmetric = B|C|D, mirror = A|C|D), on the
# single analysis set, and checks whether any gene in those pathways is BOTH a
# Kosmicki height gene and a Sacks CNV gene of the relevant class.
#
# NOTE ON PROVENANCE: the *_overlap_genes string columns in
# kosmicki_pathway_cnv_master_results.csv are unreliable - 1,286 asymmetric and
# 967 mirror rows carry a gene string whose length disagrees with that row's own
# overlap_n (many are blank). The gene sets are therefore RECONSTRUCTED here from
# the source GMT files plus kosmicki_207_genes_annotated.csv and
# cnv_genes_annotated.csv. The reconstruction is validated against the master
# overlap_n columns: 0 mismatches in 5,712 pathways for all three gene sets.
#
# Usage:  Rscript make_intersection_table.R [master_csv] [gmt_dir] [outdir]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "kosmicki_pathway_cnv_master_results.csv"
gmtdir <- if (length(args) >= 2) args[2] else "../Old/original_pipeline_review/original_bundle"
outdir <- if (length(args) >= 3) args[3] else "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05

read_gmt <- function(f) {
  l <- readLines(f, warn = FALSE); out <- list()
  for (li in l) {
    fl <- strsplit(li, "\t", fixed = TRUE)[[1]]; if (length(fl) < 3) next
    g <- unique(toupper(trimws(fl[3:length(fl)]))); g <- g[nzchar(g)]
    if (length(g) >= 5 && length(g) <= 500) out[[trimws(fl[1])]] <- g
  }
  out
}
gmt <- c(read_gmt(file.path(gmtdir, "GO_Biological_Process_2023.gmt")),
         read_gmt(file.path(gmtdir, "KEGG_2021_Human.gmt")))

k207 <- toupper(fread("kosmicki_207_genes_annotated.csv")$Gene)
cnvg <- fread("cnv_genes_annotated.csv")
gset <- list(asymmetric = unique(toupper(cnvg[Dose_response_class == "asymmetric"]$Gene)),
             mirror     = unique(toupper(cnvg[Dose_response_class == "mirror"]$Gene)))

x <- fread(input)

# ---- validate the reconstruction before relying on it -----------------------
x[, rec_kos  := vapply(pathway, function(p) length(intersect(gmt[[p]], k207)), 1L)]
x[, rec_asym := vapply(pathway, function(p) length(intersect(gmt[[p]], gset$asymmetric)), 1L)]
x[, rec_mirr := vapply(pathway, function(p) length(intersect(gmt[[p]], gset$mirror)), 1L)]
stopifnot(all(x$rec_kos  == x$kosmicki_overlap_n),
          all(x$rec_asym == x$asymmetric_cnv_overlap_n),
          all(x$rec_mirr == x$mirror_cnv_overlap_n))
cat("Reconstruction validated against master overlap counts (0 mismatches).\n")

# ---- the single analysis set: union of the two attainability sets -----------
att <- function(n, k, u) phyper(pmin(n, k) - 1, k, u - k, n, lower.tail = FALSE)
x[, mpA := att(asymmetric_cnv_pathway_n_in_universe, asymmetric_cnv_target_n, cnv_library_mapped_universe_n)]
x[, mpM := att(mirror_cnv_pathway_n_in_universe,     mirror_cnv_target_n,     cnv_library_mapped_universe_n)]
U <- x[mpA < alpha | mpM < alpha]
cat("Analysis set (union of the asymmetric and mirror attainability sets): ", nrow(U), "
", sep = "")

nom <- function(or, p) is.finite(or) & or > 1 & is.finite(p) & p < alpha
A <- nom(U$mirror_cnv_odds_ratio,     U$mirror_cnv_p)       # panel A, mirror CNV
B <- nom(U$asymmetric_cnv_odds_ratio, U$asymmetric_cnv_p)   # panel B, asymmetric CNV
C <- nom(U$kosmicki_odds_ratio,       U$kosmicki_p)         # panel C, Kosmicki
D <- nom(U$direction_odds_ratio, U$direction_fisher_p_one_sided) # panel D, bias

# Asymmetric analysis uses B, C, D; the mirror analysis uses A, C, D.
analyses <- list(
  list(cls = "asymmetric", S1 = B, L = c("B", "C", "D")),
  list(cls = "mirror",     S1 = A, L = c("A", "C", "D")))

out <- list()
for (an in analyses) {
  S1 <- an$S1; L <- an$L
  regions <- list(
    list(lab = paste(L, collapse = "|"),     sel = S1 & C & D),
    list(lab = paste(L[1], L[2], sep = "|"), sel = S1 & C & !D),
    list(lab = paste(L[1], L[3], sep = "|"), sel = S1 & D & !C),
    list(lab = paste(L[2], L[3], sep = "|"), sel = C & D & !S1))
  for (rg in regions) {
    if (!sum(rg$sel)) {
      out[[length(out) + 1]] <- data.table(
        analysis = an$cls, intersection = rg$lab, pathway_number = NA_integer_,
        pathway = "(none)", database = "", kosmicki_genes = NA_integer_,
        cnv_genes = NA_integer_, shared_genes = "")
      next
    }
    sub <- U[rg$sel]
    for (i in seq_len(nrow(sub))) {
      pg  <- gmt[[sub$pathway[i]]]
      kg  <- intersect(pg, k207)
      cg  <- intersect(pg, gset[[an$cls]])
      shr <- intersect(kg, cg)
      out[[length(out) + 1]] <- data.table(
        analysis = an$cls, intersection = rg$lab, pathway_number = i,
        pathway = sub$pathway[i],
        database = c(GO_Biological_Process_2023 = "GO:BP",
                     KEGG_2021_Human            = "KEGG")[[sub$library[i]]],
        kosmicki_genes = length(kg), cnv_genes = length(cg),
        shared_genes = if (length(shr)) paste(sort(shr), collapse = ";") else "none")
    }
  }
}
tab <- rbindlist(out)

fwrite(tab, file.path(outdir, "intersection_pathways.csv"))

cat("\n")
print(tab, nrows = 100)
cat("\nAny pathway with a gene that is both a Kosmicki height gene and a Sacks CNV gene? ",
    ifelse(any(tab$shared_genes != "none" & tab$shared_genes != ""), "YES", "NO"), "\n", sep = "")
