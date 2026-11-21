# New subjects: CHOP 22q11.2 growth records (Children's Hospital of Philadelphia 22q and You Center; T. Blaine Crowley), received
# 2026-11-19 (JS). Workbook (participant-level, controlled): G2MH/NetworkPaper/AnalysisByMolly/Clinical/Growth Dataset for Sebat.xlsx
#   (Expanse copy: /expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits/inputs/G2MH/NetworkPaper/AnalysisByMolly/Clinical/
#   CHOP_growth_2026-11-19/Growth_Dataset_for_Sebat.xlsx, md5 5f1e7cb44e55ef8a91fdfc212af00fc8).
#   Sheets: Index (1,042 carriers: study ID, deletion/duplication, breakpoints, sex, ASD diagnosis), Height (cm), Weight (kg), OFC (cm);
#   long format, one row per measurement with age in years.
# Mapping decisions (Jon 2026-11-19):
#   * Only the standard A-D CNV is used for the BMI-modifier paper ("22q11.2 A-D DEL/DUP"); nested, distal, FISH-only and
#     non-LCR classes are left out (A-B, A-C, B-D, C-D carriers are analysed separately in the dosage-architecture paper).
#   * One row per encounter (study ID x age). Height and weight are paired when recorded at the same age (same encounter; repeated
#     entries at one encounter are averaged); BMI = weight / height^2. Head circumference is joined at the same encounter; OFC-only
#     encounters are kept as rows without height/BMI.
#   * ASD diagnosis: Index column "ASD diagnosis" (Yes = 1, No = 0). The sheet does not say whether "No" means evaluated.
# Overlap: CHOP is also a G2MH recruitment site (G2MH site_id "penn_chop"; MSassigned_ID prefix 203), and no CHOP <-> G2MH identifier
#   map exists yet, so the same child could be in both. add_chop(g, overlap = "drop_g2mh_chop_site") removes the G2MH 22q11.2 A-D
#   carriers recruited at that site from the growth table (their growth is expected to be in the CHOP records); overlap = "keep"
#   keeps both (sensitivity). Non-carrier controls are untouched. See Jon/Claude/22q_subregions_clinical/SAMPLE_OVERLAP_AND_FIG1_PLAN.md.
# add_chop(g) appends the rows to a genomatrix-format table (cohort "CHOP", sampleID / MSassigned_ID "CHOP_<study id>", GUID NA).

CHOP_XLSX <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/G2MH/NetworkPaper/AnalysisByMolly/Clinical/Growth Dataset for Sebat.xlsx"

chop_rows <- function(f = CHOP_XLSX) {
  rd <- function(s) { v <- data.table::as.data.table(readxl::read_excel(f, sheet = s)); v[, .(sid = as.character(study_id), age = as.numeric(vitals_age_years), val = as.numeric(vitals_value))] }
  I <- data.table::as.data.table(readxl::read_excel(f, sheet = "Index"))
  data.table::setnames(I, c("Study ID", "22q Type", "22q Brkpts", "Sex", "ASD diagnosis"), c("sid", "type", "brk", "sex0", "asd0"))
  I <- I[brk == "A-D", .(sid = as.character(sid), CNV = ifelse(type == "Deletion", "22q11.2 A-D DEL", "22q11.2 A-D DUP"),
                         sex_numeric = ifelse(sex0 == "Male", 1L, ifelse(sex0 == "Female", 2L, NA_integer_)),
                         asd = ifelse(asd0 == "Yes", 1L, ifelse(asd0 == "No", 0L, NA_integer_)))]
  H <- rd("Height")[, .(h = mean(val)), by = .(sid, age)]; W <- rd("Weight")[, .(w = mean(val)), by = .(sid, age)]
  O <- rd("OFC")[, .(hc = mean(val)), by = .(sid, age)]
  e <- merge(merge(H, W, by = c("sid", "age"), all = TRUE), O, by = c("sid", "age"), all = TRUE)
  e <- merge(e, I, by = "sid")
  e[, `:=`(HEIGHT_raw = h, BMI_raw = ifelse(!is.na(h) & !is.na(w) & h > 0, w / (h / 100)^2, NA_real_), HEADCIRCUM_raw = hc)]
  e <- e[!is.na(HEIGHT_raw) | !is.na(HEADCIRCUM_raw)]
  data.table::setorder(e, sid, age)
  e[, .(sampleID = paste0("CHOP_", sid), MSassigned_ID = paste0("CHOP_", sid), GUID = NA_character_, cohort = "CHOP", CNV, sex_numeric, age,
        HEIGHT_raw, BMI_raw, HEADCIRCUM_raw, asd, timepoint = NA_integer_, Cardiff_mask = FALSE)]
}

add_chop <- function(g, overlap = Sys.getenv("CHOP_OVERLAP", "drop_g2mh_chop_site")) {
  n <- chop_rows()
  cat("CHOP rows appended:", nrow(n), "encounters,", data.table::uniqueN(n$sampleID), "carriers | by CNV:",
      paste(names(table(unique(n[, .(sampleID, CNV)])$CNV)), table(unique(n[, .(sampleID, CNV)])$CNV), collapse = ", "),
      "| ASD yes/no:", sum(unique(n[, .(sampleID, asd)])$asd == 1, na.rm = TRUE), "/", sum(unique(n[, .(sampleID, asd)])$asd == 0, na.rm = TRUE), "\n")
  if (overlap == "drop_g2mh_chop_site") {
    drop <- g$cohort == "G2MH" & grepl("^203-", g$MSassigned_ID) & grepl("^22q11.2", trimws(g$CNV))
    cat("CHOP overlap rule 'drop_g2mh_chop_site': removed", sum(drop), "G2MH rows (", data.table::uniqueN(g$MSassigned_ID[drop]),
        "22q11.2 carriers recruited at the G2MH CHOP site, MSassigned_ID 203-*)\n")
    g <- g[!drop]
  } else cat("CHOP overlap rule 'keep': G2MH CHOP-site carriers kept (possible double counting)\n")
  rbind(g, n, fill = TRUE)
}
