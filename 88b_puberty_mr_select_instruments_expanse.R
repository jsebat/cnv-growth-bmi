# 88b_puberty_mr_select_instruments_expanse.R  (work package E, BMI modifier paper v2; 2026-11-20). Runs on Expanse.
# What it does: picks SNP instruments for the two-sample Mendelian randomisation of childhood body mass -> age at menarche ->
#   adult height, from three GWAS that do NOT include UK Biobank:
#     Perry 2014 (age at menarche, ReproGen, up to 182,416 women; PMID 25231870; downloaded to Expanse 2026-11-20 from reprogen.org),
#     Vogelezang 2020 (childhood BMI, EGG, up to 61,111 children; GCST90002409),
#     Bradfield 2019 (childhood obesity, EGG trans-ancestral; European-ancestry columns used here; PMID 31504550).
#   Genome-wide significant SNPs (P < 5e-8) present in the UK Biobank imputed bed files are thinned greedily by P, one SNP per
#   +/- 1 Mb (LD pruning on UK Biobank genotypes happens in 88a). For every selected SNP the effect in all three GWAS is looked up
#   (needed for the multivariable model).
# Inputs: the three summary-statistic files, UK Biobank .bim files (j3guevar/archive/ukb_imp_chrN_bgen1).
# Outputs (working directory): instrument_candidates.csv, extract_chr<N>.txt (SNP lists for plink2 --extract).
# null model: none; instrument selection only (no test).
suppressPackageStartupMessages(library(data.table))
WD <- getwd()
BIM <- "/expanse/projects/sebat1/j3guevar/archive/ukb_imp_chrN_bgen1/ukb_imp_chr%d_v3_qc_uniq.bim"
PERRY <- "/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120/00_external_sumstats/perry2014/Menarche_Nature2014_GWASMetaResults_17122014.txt"
RAW <- "/expanse/projects/sebat1/a1sriniv/pgs_sumstats_organized/raw_sumstats"
bim <- rbindlist(lapply(1:22, function(i) fread(sprintf(BIM, i), header = FALSE, col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))))[, cm := NULL]
bim <- bim[nchar(a1) == 1 & nchar(a2) == 1]
cat("UKB bim SNVs:", nrow(bim), "\n")
num <- function(x) suppressWarnings(as.numeric(x))

pe <- fread(PERRY); setnames(pe, c("snp", "ea", "oa", "frq", "beta_gwas", "p_gwas", "icogs", "beta", "p"))
pe[, `:=`(beta = num(beta), p = num(p), frq = num(frq))]; pe <- pe[!is.na(beta) & !is.na(p) & p > 0]
pe[, se := abs(beta) / qnorm(p / 2, lower.tail = FALSE)]; pe <- pe[is.finite(se) & se > 0]      # no SE in the file: derived from beta and P
pe <- pe[, .(snp, ea = toupper(ea), oa = toupper(oa), frq, beta, se, p)]
vo <- fread(file.path(RAW, "GCST90002409_childhood_BMI.tsv.gz"))
vo <- vo[, .(snp = variant_id, ea = toupper(effect_allele), oa = toupper(other_allele), frq = NA_real_, beta = num(beta), se = num(standard_error), p = num(p_value), n = TotalSampleSize)]
vo <- vo[!is.na(beta) & !is.na(p)]
br <- fread(file.path(RAW, "CHILDHOOD_OBESITY.TRANS_ANCESTRAL_2019.RESULTS.txt.gz"), select = c("CHR", "POS", "EA", "OA", "EUR_N", "EUR_FRQ", "EUR_BETA", "EUR_SE", "EUR_P"))
br <- br[, .(chr = num(CHR), pos = num(POS), ea = toupper(EA), oa = toupper(OA), n = num(EUR_N), frq = num(EUR_FRQ), beta = num(EUR_BETA), se = num(EUR_SE), p = num(EUR_P))]
br <- br[!is.na(beta) & !is.na(p) & !is.na(chr)]
br <- merge(br, bim[, .(chr, pos, snp)], by = c("chr", "pos"))[, .(snp, ea, oa, frq, beta, se, p, n)]
cat("rows: Perry", nrow(pe), " Vogelezang", nrow(vo), " Bradfield EUR (in UKB bim)", nrow(br), "\n")
cat("Vogelezang max N", max(vo$n, na.rm = TRUE), " Bradfield EUR max N", max(br$n, na.rm = TRUE), "\n")

thin <- function(x, nm) {                       # greedy: best P first, drop everything within +/- 1 Mb
  x <- merge(x[p < 5e-8], bim[, .(snp, chr, pos)], by = "snp")[order(p)]
  cat(nm, ": genome-wide significant SNPs in UKB bim:", nrow(x), "\n")
  keep <- logical(nrow(x))
  for (i in seq_len(nrow(x))) if (!any(keep & x$chr == x$chr[i] & abs(x$pos - x$pos[i]) < 1e6)) keep[i] <- TRUE
  cat(nm, ": after 1 Mb thinning:", sum(keep), "\n"); x[keep, snp]
}
s_pe <- thin(pe, "Perry menarche"); s_vo <- thin(vo, "Vogelezang childhood BMI"); s_br <- thin(br, "Bradfield childhood obesity")
u <- unique(c(s_pe, s_vo, s_br))
C <- bim[snp %in% u]; C <- C[!duplicated(snp)]
C[, `:=`(sel_perry = snp %in% s_pe, sel_vogelezang = snp %in% s_vo, sel_bradfield = snp %in% s_br)]
add <- function(C, x, tag) { y <- x[snp %in% C$snp][!duplicated(snp), .(snp, ea, oa, frq, beta, se, p)]; setnames(y, c("ea", "oa", "frq", "beta", "se", "p"), paste0(c("ea", "oa", "frq", "beta", "se", "p"), "_", tag)); merge(C, y, by = "snp", all.x = TRUE) }
C <- add(C, pe, "perry"); C <- add(C, vo, "vogelezang"); C <- add(C, br, "bradfield")
fwrite(C[order(chr, pos)], file.path(WD, "instrument_candidates.csv"))
for (i in 1:22) writeLines(C[chr == i, snp], file.path(WD, sprintf("extract_chr%d.txt", i)))
print(C[, .(n = .N, with_perry_beta = sum(!is.na(beta_perry)), with_vog_beta = sum(!is.na(beta_vogelezang)), with_brad_beta = sum(!is.na(beta_bradfield))), by = .(sel_perry, sel_vogelezang, sel_bradfield)])
cat("done\n")
