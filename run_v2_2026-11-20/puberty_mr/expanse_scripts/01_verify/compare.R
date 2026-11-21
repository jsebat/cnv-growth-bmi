# Compare re-run CSVs with the originals: max absolute and relative difference over numeric columns. null model: none; reproduction check.
suppressPackageStartupMessages(library(data.table))
O <- commandArgs(TRUE)[1]
pairs <- list(c("ukbb_bodysize10_vs_adult_trends.csv", "ukbb_bodysize"), c("ukbb_bodysize10_vs_adult_counts.csv", "ukbb_bodysize"),
              c("ukbb_mr_childhood_obesity.csv", "ukbb_mr"), c("ukbb_mr_puberty_pathway_A.csv", "ukbb_mr"), c("ukbb_mr_puberty_pathway_B.csv", "ukbb_mr"),
              c("ukbb_mr_puberty_pathway_C.csv", "ukbb_mr"), c("mr_menarche_height_ukbb.csv", "ukbb_mr"), c("fig_supp5_combined_tables.csv", "ukbb_mr"))
res <- rbindlist(lapply(pairs, function(p) {
  if (!file.exists(p[1])) return(data.table(file = p[1], status = "missing in re-run"))
  a <- fread(p[1]); b <- fread(file.path(O, p[2], p[1]))
  if (!identical(dim(a), dim(b))) return(data.table(file = p[1], status = paste("dim differs", paste(dim(a), collapse = "x"), "vs", paste(dim(b), collapse = "x"))))
  key <- names(a)[!sapply(a, is.numeric)]; if (length(key)) { setorderv(a, key); setorderv(b, key) }
  num <- names(a)[sapply(a, is.numeric)]
  d <- unlist(lapply(num, function(v) abs(a[[v]] - b[[v]]))); r <- unlist(lapply(num, function(v) abs(a[[v]] - b[[v]]) / pmax(abs(b[[v]]), 1e-300)))
  data.table(file = p[1], status = "compared", n_numbers = sum(!is.na(d)), max_abs_diff = max(d, na.rm = TRUE), max_rel_diff = max(r, na.rm = TRUE),
             nonnumeric_identical = all(sapply(key, function(v) identical(as.character(a[[v]]), as.character(b[[v]])))))
}), fill = TRUE)
fwrite(res, "reproduction_check.csv"); print(res)
