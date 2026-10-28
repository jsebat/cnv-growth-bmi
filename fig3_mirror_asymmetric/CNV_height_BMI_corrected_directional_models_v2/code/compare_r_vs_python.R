#!/usr/bin/env Rscript
# Compare the R pipeline outputs (R_outputs/) against the shipped Python outputs
# (height/outputs/, bmi/outputs/). Writes qa/r_vs_python_comparison.csv.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", flag)
  args[hit + 1]
}
project_dir <- normalizePath(get_arg("--project-dir", "."), mustWork = TRUE)
tol <- as.numeric(get_arg("--tol", "1e-4"))

read_csv_base <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

# Shipped file -> R counterpart. The R driver names the sparse model "single";
# the Python driver names it "single_hit".
pairs <- list()
for (tr in c("height", "bmi")) {
  shipped_dir <- file.path(project_dir, tr, "outputs")
  r_dir <- file.path(project_dir, "R_outputs", tr)
  for (f in list.files(shipped_dir, pattern = "\\.csv$")) {
    cand <- file.path(r_dir, f)
    if (!file.exists(cand)) cand <- file.path(r_dir, sub("_single_hit", "_single", f))
    pairs[[length(pairs) + 1]] <- list(trait = tr, file = f, shipped = file.path(shipped_dir, f), r = cand)
  }
}

rows <- list()
for (p in pairs) {
  if (!file.exists(p$r)) {
    rows[[length(rows) + 1]] <- data.frame(
      trait = p$trait, file = p$file, status = "MISSING_IN_R",
      n_numeric_cells = NA_integer_, max_abs_diff = NA_real_, max_rel_diff = NA_real_,
      cols_only_in_python = NA_character_, cols_only_in_R = NA_character_,
      stringsAsFactors = FALSE)
    next
  }
  a <- read_csv_base(p$shipped); b <- read_csv_base(p$r)
  only_py <- setdiff(names(a), names(b)); only_r <- setdiff(names(b), names(a))
  shared <- intersect(names(a), names(b))
  status <- "OK"
  if (nrow(a) != nrow(b)) status <- "ROW_COUNT_DIFFERS"

  n_cells <- 0L; max_abs <- 0; max_rel <- 0
  if (nrow(a) == nrow(b)) {
    for (cn in shared) {
      x <- a[[cn]]; y <- b[[cn]]
      if (is.numeric(x) && is.numeric(y)) {
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {
          d <- abs(x[ok] - y[ok])
          denom <- pmax(abs(x[ok]), abs(y[ok]), 1e-12)
          n_cells <- n_cells + sum(ok)
          max_abs <- max(max_abs, max(d))
          max_rel <- max(max_rel, max(d / denom))
        }
      } else if (!identical(as.character(x), as.character(y))) {
        status <- if (status == "OK") "TEXT_DIFFERS" else status
      }
    }
    # A file with no shared numeric columns was never actually compared;
    # reporting it as OK is a false pass.
    if (n_cells == 0) status <- "NO_NUMERIC_OVERLAP"
    else if (status == "OK" && max_abs > tol) status <- "NUMERIC_DIFFERS"
  }
  rows[[length(rows) + 1]] <- data.frame(
    trait = p$trait, file = p$file, status = status,
    n_numeric_cells = n_cells, max_abs_diff = max_abs, max_rel_diff = max_rel,
    cols_only_in_python = paste(only_py, collapse = ";"),
    cols_only_in_R = paste(only_r, collapse = ";"),
    stringsAsFactors = FALSE)
}

out <- do.call(rbind, rows)
dir.create(file.path(project_dir, "qa"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(project_dir, "qa", "r_vs_python_comparison.csv"), row.names = FALSE)
print(out[, c("trait", "file", "status", "n_numeric_cells", "max_abs_diff", "max_rel_diff")], row.names = FALSE)
cat("\nColumns present in Python outputs but absent from R outputs:\n")
for (i in seq_len(nrow(out))) if (nzchar(out$cols_only_in_python[i]))
  cat("  ", out$file[i], ": ", out$cols_only_in_python[i], "\n", sep = "")
cat("\nWrote ", file.path(project_dir, "qa", "r_vs_python_comparison.csv"), "\n", sep = "")
