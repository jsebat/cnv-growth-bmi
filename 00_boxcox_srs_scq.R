# Box-Cox scale for SRS and SCQ (JS / Claude 2026-11-14; Jon: Box-Cox is the scaling recommended by Peter Visscher and used for the
# biobank traits). Sourced by 24_cognition_bmi_slopes.R when Sys.getenv("SRSSCQ_SCALE") == "boxcox", after D exists (21).
# Raw totals + 1 are Box-Cox transformed with lambda chosen by profile likelihood in the unaffected SSC relatives of the same form
# (raw + 1 ~ age + age^2 + sex; SRS child form, SRS adult form, SCQ) and z-scored against those relatives. The linear z is kept as
# y_lin. Same transform as 51_cog_panels_boxcox.R.
if (!"y_lin" %in% names(D)) {
  .rs <- list(child = ssc[form == "child" & role == "sibling", srs], adult = ssc[form == "adult" & role == "parent", srs])
  D[pheno == "SRS", form2 := fifelse(abs(y - zref(raw, .rs$child)) < 1e-8, "child", fifelse(abs(y - zref(raw, .rs$adult)) < 1e-8, "adult", NA_character_))]
  D[pheno == "SCQ", form2 := "scq"]
  stopifnot(!anyNA(D[pheno %in% c("SRS", "SCQ"), form2]))
  D[, y_lin := y]; LAM <- list()
  for (.f in c("child", "adult", "scq")) {
    .ph <- if (.f == "scq") "SCQ" else "SRS"
    .r <- D[pheno == .ph & form2 == .f & unaffected == TRUE & locus == "No CNV"]; .r[, `:=`(yy = raw + 1, a2 = age^2)]
    .b <- MASS::boxcox(yy ~ age + a2 + sex, data = .r, lambda = seq(-2, 2, 0.01), plotit = FALSE); .L <- .b$x[which.max(.b$y)]; LAM[[.f]] <- .L
    .bc <- function(x) if (abs(.L) < 1e-8) log(x + 1) else ((x + 1)^.L - 1) / .L
    D[pheno == .ph & form2 == .f, y := zref(.bc(raw), .bc(.r$raw))]
  }
  cat("SRS/SCQ on the Box-Cox scale; lambda:", paste(names(LAM), round(unlist(LAM), 2), collapse = "; "), "\n")
}
