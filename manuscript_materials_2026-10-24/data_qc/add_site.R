# Add a G2MH site column (release site_id, matched by GUID) to the flagged-record tables (JS 2026-10-27).
suppressMessages(library(data.table))
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
H <- file.path(DB, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")
rel <- unique(fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"), select = c("guid", "site_id")))
SITE <- c(penn_chop = "Penn/CHOP", toro = "Toronto", leuv = "Leuven", mass = "Massachusetts (MGH)", card = "Cardiff", mont = "Montreal", ucla = "UCLA")
PREFIX <- c(`192` = "penn_chop", `207` = "toro", `208` = "leuv", `209` = "mass", `210` = "card", `211` = "mont", `213` = "ucla")
add_site <- function(f) {
  x <- fread(f, colClasses = "character")
  if ("site" %in% names(x)) x[, site := NULL]
  s <- rel$site_id[match(x$GUID, rel$guid)]
  pre <- sub("-.*", "", x$sampleID); s[is.na(s) & x$cohort %like% "G2MH"] <- PREFIX[pre[is.na(s) & x$cohort %like% "G2MH"]]
  x[, site := ifelse(cohort %like% "G2MH", ifelse(is.na(s), "", s), "")]
  setcolorder(x, c("cohort", "site", setdiff(names(x), c("cohort", "site"))))
  fwrite(x, f); x
}
x <- add_site(file.path(H, "implausible_height_weight.csv"))
xc <- add_site(file.path(H, "implausible_height_weight_corrections.csv"))
cat("G2MH records by site:\n"); print(x[cohort %like% "G2MH", .N, by = site][order(-N)])
# rebuild the markdown table (keep the text above it)
md <- readLines(file.path(H, "implausible_height_weight.md"), encoding = "UTF-8")
hdr <- md[seq_len(grep("^\\| cohort", md)[1] - 1)]
cols <- c("cohort", "site", "sampleID", "GUID", "CNV", "age", "sex", "height_cm", "weight_kg_backcalc", "BMI", "hc_cm", "height_Z", "child_BMI_Z", "hc_Z", "flag")
lab <- c("cohort", "site", "sampleID", "GUID", "CNV", "age", "sex", "height_cm", "weight_kg*", "BMI", "hc_cm", "height_Z", "child_BMI_Z", "hc_Z", "flag")
rows <- apply(x[, ..cols], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
tail_txt <- md[grep("^\\| ", md)]; after <- md[seq(max(grep("^\\| ", md)) + 1, length(md))]
writeLines(c(hdr, paste0("| ", paste(lab, collapse = " | "), " |"), paste0("|", paste(rep("---", length(lab)), collapse = "|"), "|"), rows, after), file.path(H, "implausible_height_weight.md"), useBytes = TRUE)
for (f in c("implausible_height_weight.csv", "implausible_height_weight.md", "implausible_height_weight_corrections.csv")) file.copy(file.path(H, f), file.path(H, "manuscript_materials_2026-10-24/data_qc", f), overwrite = TRUE)
cat("done; md rows:", length(rows), "\n")
