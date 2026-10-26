# Plausibility QC for growth measures (JS / Claude 2026-10-24). Applied to the growth genomatrix (and to any other
# height/weight source) before scoring. Rules, in order:
#   1. height < 3            -> entered in metres: multiply by 100 (G2MH site 213 block; BMI there is already valid)
#   2. height outside 60-230 -> height NA, and BMI NA (genomatrix BMI is derived from that height)
#   3. BMI outside 10-70     -> BMI NA
#   4. |height Z| > 5 (WHO)  -> height NA
#   5. child (< 19 y) |BMI Z| > 5 (WHO) -> BMI NA
#   6. head circumference outside 30-70 cm or |HC Z| > 5 (UK90, adults scored at 18) -> HC NA
# Values are set to NA rather than rows dropped, so an individual with a bad height keeps a valid HC, etc.
# Flagged rows are listed in implausible_height_weight.{csv,md}. Requires 00_scoring.R (score_height, score_bmi_child, score_hc).
# Manual corrections: implausible_height_weight_corrections.csv (same folder) lists every flagged visit with blank columns
# action (correct / exclude / keep), corrected_height_cm, corrected_weight_kg, corrected_hc_cm. Rows with action "correct"
# replace the raw values (BMI recomputed from corrected or original height and weight) before the rules run; "exclude"
# sets height, BMI and HC to NA; "keep" exempts the row from the rules; blank = untouched. Matching: cohort + sampleID +
# age rounded to 1 decimal (G2MH release rows are matched by GUID in 21_pheno_bmi_gam.R's relb and are not handled here).
QC_CORRECTIONS <- file.path(if (exists("HERE")) HERE else getwd(), "implausible_height_weight_corrections.csv")
apply_corrections <- function(x, path = QC_CORRECTIONS, verbose = TRUE) {
  if (!file.exists(path)) return(x)
  k <- data.table::fread(path, colClasses = "character")[action %in% c("correct", "exclude", "keep")]
  if (!nrow(k) || !all(c("cohort", "sampleID", "age") %in% names(x))) return(x)
  k[, age1 := round(as.numeric(age), 1)]; x[, age1 := round(age, 1)]
  n_c <- 0L; n_e <- 0L; n_k <- 0L
  for (i in seq_len(nrow(k))) {
    j <- which(x$cohort == k$cohort[i] & x$sampleID == k$sampleID[i] & x$age1 == k$age1[i]); if (!length(j)) next
    if (k$action[i] == "exclude") { x[j, c("HEIGHT_raw", "BMI_raw", "HEADCIRCUM_raw") := .(NA_real_, NA_real_, NA_real_)]; n_e <- n_e + length(j) }
    else if (k$action[i] == "correct") {
      h <- suppressWarnings(as.numeric(k$corrected_height_cm[i])); w <- suppressWarnings(as.numeric(k$corrected_weight_kg[i])); c <- suppressWarnings(as.numeric(k$corrected_hc_cm[i]))
      h0 <- x$HEIGHT_raw[j]; w0 <- x$BMI_raw[j] * (h0 / 100)^2
      hn <- if (!is.na(h)) h else h0; wn <- if (!is.na(w)) w else w0
      x[j, HEIGHT_raw := hn]; x[j, BMI_raw := wn / (hn / 100)^2]; if (!is.na(c)) x[j, HEADCIRCUM_raw := c]
      x[j, qc_keep := TRUE]; n_c <- n_c + length(j) }
    else { x[j, qc_keep := TRUE]; n_k <- n_k + length(j) }
  }
  x[, age1 := NULL]
  if (verbose) cat(sprintf("manual corrections applied: %d corrected, %d excluded, %d kept as is\n", n_c, n_e, n_k))
  x
}
qc_growth <- function(x, height = "HEIGHT_raw", bmi = "BMI_raw", hc = "HEADCIRCUM_raw", age = "age", sex = "sex", label = "growth table", verbose = TRUE) {
  x <- data.table::copy(x)
  if (!"qc_keep" %in% names(x)) x[, qc_keep := FALSE]
  x <- apply_corrections(x, verbose = verbose)
  keep <- x$qc_keep; x[, qc_keep := NULL]
  if (any(keep)) { kept <- x[keep]; x <- x[!keep] }   # rows marked keep/correct bypass the rules and are re-attached below
  h <- x[[height]]; b <- if (bmi %in% names(x)) x[[bmi]] else rep(NA_real_, nrow(x)); c <- if (!is.null(hc) && hc %in% names(x)) x[[hc]] else NULL
  a <- x[[age]]; s <- x[[sex]]
  n_m <- sum(!is.na(h) & h < 3); h[!is.na(h) & h < 3] <- h[!is.na(h) & h < 3] * 100
  bad_h <- !is.na(h) & (h < 60 | h > 230); n_hr <- sum(bad_h); b[bad_h] <- NA; h[bad_h] <- NA
  bad_b <- !is.na(b) & (b < 10 | b > 70); n_br <- sum(bad_b); b[bad_b] <- NA
  hz <- score_height(h, a, s); bad_hz <- !is.na(hz) & abs(hz) > 5; n_hz <- sum(bad_hz); h[bad_hz] <- NA
  bz <- ifelse(!is.na(a) & a < 19, score_bmi_child(b, a, s), NA_real_); bad_bz <- !is.na(bz) & abs(bz) > 5; n_bz <- sum(bad_bz); b[bad_bz] <- NA
  n_hc <- 0L
  if (!is.null(c)) { cz <- score_hc(c, a, s); bad_c <- !is.na(c) & (c < 30 | c > 70 | (!is.na(cz) & abs(cz) > 5)); n_hc <- sum(bad_c); c[bad_c] <- NA; data.table::set(x, j = hc, value = c) }
  data.table::set(x, j = height, value = h); if (bmi %in% names(x)) data.table::set(x, j = bmi, value = b)
  if (verbose) cat(sprintf("QC %s: %d heights in metres x100; set NA: %d heights outside 60-230 cm (+BMI), %d BMI outside 10-70, %d |height Z| > 5, %d child |BMI Z| > 5, %d head circumference\n", label, n_m, n_hr, n_br, n_hz, n_bz, n_hc))
  if (any(keep)) x <- data.table::rbindlist(list(x, kept), use.names = TRUE)
  x
}
