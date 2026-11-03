# New subjects: CHU Sainte-Justine (Montreal) archival carriers, ARPA-H project (JS 2026-11-02).
# Source workbook (participant-level, kept in the controlled data folder): Jon/Claude/CNV_growth_trajectories_data/
#   chusj/16p_22q-BMI_ASD_data.xlsx (37 carriers, one visit each; sex, CNV, age "15y9m", height, weight, HC, BMI, ASD).
# Jon (2026-11-02): enrolled in another study, no overlap with G2MH (CHUSJ is also a G2MH site, but these are archival).
# Provided by Marilyse Gallant (Jacquemont lab, CHU Sainte-Justine), e-mail 2026-11-02: 37 patients with a 16p/22q variant
#   (22q11.2 A-D only) not in SFARI/G2MH; weight and height noted even where no BMI was calculated; HC and ASD where available.
#   ASD coding in the sheet: "Not tested" = extensively evaluated at Sainte-Justine for other neurodevelopmental disorders but not
#   referred for ASD assessment (usually no signs of ASD); "Unknown" = little information, evaluation elsewhere not excluded.
# Mapping decisions (Jon's defaults, 2026-11-02): every 16p11.2 label = BP4-BP5 proximal; every 22q label
#   (22q11.1 / 22q11.2 / 22q11.21) = A-D; ASD "Not tested"/"Unknown" = NA, "No ASD" = 0, "ASD level ..." = 1;
#   one row (F, 16y11m, Dup 16p11.2) had height 76 / weight 174.3, consistent with its BMI 25.0 only when swapped -> swapped.
# BMI is recomputed from height and weight (the sheet's BMI column is checked against it).
# add_chusj(g) appends the rows to a genomatrix-format table (cohort "CHUSJ", sampleID "CHUSJ_01".., GUID NA).

chusj_rows <- function() {
  f <- file.path("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data/chusj", "16p_22q-BMI_ASD_data.xlsx")
  x <- as.data.table(readxl::read_excel(f))
  setnames(x, c("sex", "group", "age_txt", "height", "weight", "hc", "bmi_sheet", "asd_txt"))
  x[, sampleID := sprintf("CHUSJ_%02d", seq_len(.N))]
  for (v in c("height", "weight", "hc", "bmi_sheet")) x[, (v) := suppressWarnings(as.numeric(get(v)))]
  # age "15y9m", "6m", "1y0m" -> years
  a <- trimws(x$age_txt)
  yr <- suppressWarnings(ifelse(grepl("y", a), as.numeric(sub("y.*", "", a)), 0))
  mo <- suppressWarnings(ifelse(grepl("m", a), as.numeric(sub("m.*", "", sub(".*y", "", a))), 0))
  x[, age := yr + mo / 12]
  # the swapped row
  sw <- x[, which(!is.na(height) & !is.na(weight) & height < weight & age > 2)]
  if (length(sw)) { h <- x$height[sw]; x[sw, height := weight]; x[sw, weight := h] }
  x[, HEIGHT_raw := height][, HEADCIRCUM_raw := hc]
  x[, BMI_raw := ifelse(!is.na(height) & !is.na(weight), weight / (height / 100)^2, NA_real_)]
  x[, bmi_check := round(BMI_raw - bmi_sheet, 1)]
  g <- trimws(x$group)
  x[, CNV := fcase(grepl("^Del 16p", g), "16p11.2 proximal DEL", grepl("^Dup 16p", g), "16p11.2 proximal DUP",
                   grepl("^Del 22q", g), "22q11.2 A-D DEL", grepl("^Dup 22q", g), "22q11.2 A-D DUP")]
  x[, locus_label_source := g]
  x[, sex_numeric := fifelse(sex == "M", 1L, fifelse(sex == "F", 2L, NA_integer_))]
  # Jon (2026-11-02): "Not tested" (evaluated for other NDDs, not referred for ASD assessment) = no ASD, flagged in asd_source so the
  # ASD-by-CNV analysis can be rerun without them; "Unknown" = NA.
  x[, asd := fcase(asd_txt == "No ASD", 0L, asd_txt == "Not tested", 0L, grepl("^ASD", asd_txt), 1L, default = NA_integer_)]
  x[, asd_source := fcase(asd_txt == "No ASD", "evaluated", asd_txt == "Not tested", "not tested", grepl("^ASD", asd_txt), "evaluated", default = NA_character_)]
  x[, .(sampleID, MSassigned_ID = sampleID, GUID = NA_character_, cohort = "CHUSJ", CNV, locus_label_source, sex_numeric, age,
        HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, bmi_check, asd, asd_txt, timepoint = 1L, last_timepoint = TRUE, first_timepoint = TRUE,
        Cardiff_mask = FALSE, swapped_hw = seq_len(.N) %in% sw)]
}
add_chusj <- function(g) {
  n <- chusj_rows()
  cat("CHUSJ rows appended:", nrow(n), "| by CNV:", paste(names(table(n$CNV)), table(n$CNV), collapse = ", "),
      "| ASD known:", sum(!is.na(n$asd)), "| swapped height/weight:", sum(n$swapped_hw), "\n")
  if (any(abs(n$bmi_check) > 0.3, na.rm = TRUE)) cat("  NOTE: sheet BMI differs from height/weight BMI by > 0.3 in", sum(abs(n$bmi_check) > 0.3, na.rm = TRUE), "rows\n")
  rbind(g, n[, setdiff(names(n), c("locus_label_source", "bmi_check", "asd_txt", "swapped_hw", "asd")), with = FALSE], fill = TRUE)
}
