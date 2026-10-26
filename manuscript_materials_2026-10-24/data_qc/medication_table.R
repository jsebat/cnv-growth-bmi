# Rebuilds medication_sensitivity_table.md from medication_sensitivity.csv + medication_prevalence.csv (all traits incl. height / HC)
suppressMessages(library(data.table))
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
S <- fread(file.path(HERE, "medication_sensitivity.csv")); P <- fread(file.path(HERE, "medication_prevalence.csv"))
pf <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), as.character(signif(p, 2)))
cell <- function(s, e, p, n = NULL) ifelse(is.na(s), "n/a", paste0(sprintf("%+.2f (%.2f), p = %s", s, e, pf(p)), if (!is.null(n)) sprintf(" (%d)", n) else ""))
get <- function(ph, g, mod) S[pheno == ph & group == g & model == mod]
PH <- c("Height", "Head circumference", "SRS", "SCQ", "CNB accuracy", "CNB speed"); GR <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup", "All CNVs")
rows <- list()
for (ph in PH) for (g in GR) {
  u <- get(ph, g, "with medication info, unadjusted"); if (!nrow(u)) next
  a <- get(ph, g, "+ medication classes"); o <- get(ph, g, "no current psychotropic"); s <- get(ph, g, "+ medication classes + ASD"); pr <- P[pheno == ph & locus == g]
  rows[[length(rows) + 1]] <- data.table(Trait = ph, Group = g, n = u$n, `% psychotropic` = round(pr$psychotropic), `% antipsychotic` = round(pr$antipsychotic),
    Unadjusted = cell(u$slope, u$se, u$p), `+ medication classes` = if (nrow(a)) cell(a$slope, a$se, a$p) else "n/a",
    `No psychotropic (n)` = if (nrow(o)) cell(o$slope, o$se, o$p, o$n) else "n/a", `+ classes + ASD` = if (nrow(s)) cell(s$slope, s$se, s$p) else "n/a") }
T <- rbindlist(rows)
md <- c(paste0("| ", paste(names(T), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(T)), collapse = "|"), "|"), apply(T, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(HERE, "medication_sensitivity_table.md"))
# splice the table into medication_sensitivity_results.md (between the "| Trait |" header and the first blank line after it)
R <- file.path(HERE, "medication_sensitivity_results.md"); L <- readLines(R)
i <- grep("^\\| Trait \\| Group \\|", L)[1]; j <- i; while (j < length(L) && nzchar(L[j + 1])) j <- j + 1
L <- c(L[seq_len(i - 1)], md, L[seq(j + 1, length(L))]); writeLines(L, R)
cat("rows:", nrow(T), "\n"); print(T[Trait %in% c("Height", "Head circumference"), .(Trait, Group, n, Unadjusted, `+ medication classes`, `No psychotropic (n)`)])
