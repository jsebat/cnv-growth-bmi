# federated/00_select_visit.R  (copied 2026-11-20)
# Copied UNCHANGED (everything below this header) from 00_select_visit.R in the root of this repository
# (jsebat/cnv-growth-bmi), for the federated replication package (site_export.R): one visit per person for the
# cross-sectional models, rule "fill_gaps" (default) or "latest".
# Inputs: a data.table of visits with columns id, age, BMI_raw, the raw outcome column and locus. Outputs: none (functions only).
# null model: none; this file only selects visits.

# One visit per person for the cross-sectional growth models (JS 2026-11-19, Jon: "22q11.2 deletion is missing the youngest age
# group (< 7 y) and we also have a gap between 18 and 40 for most groups. We want to fill these gaps so that we can get a full
# trajectory ... choose the youngest age instead of the oldest for some people ... record the selection criteria").
#
# Rule "fill_gaps" (default; VISIT_RULE env, "latest" restores the previous rule):
#   1. Eligible visits: the outcome (height or head circumference) and BMI both measured, age 2-73.
#   2. Age bands: 2-7, 7-13, 13-19, 19-40, 40-73 years (BAND_BREAKS).
#   3. A person whose eligible visits all fall in one band keeps the latest visit in that band (= the previous rule).
#   4. A person with eligible visits in several bands is assigned to ONE of those bands, chosen to fill the emptiest band of their
#      genotype: people are processed in order of fewest candidate bands first, then by ID (deterministic); each is placed in the
#      candidate band with the fewest people of the same genotype so far (counting step 3 first); ties go to the oldest band, so
#      the previous "latest visit" choice is kept unless a younger band is emptier. Within the chosen band one visit is drawn at
#      random (seed VISIT_SEED, default 20261119), so ages spread across the band instead of piling at its upper edge.
#   5. Still one record per person, so the cross-sectional models are unchanged; only which visit represents a person changes.
# select_visit() returns the chosen rows with columns band (chosen band) and rule_changed (TRUE if the visit differs from the
# latest eligible visit). visit_selection_table() tabulates people by genotype x band under both rules.
BAND_BREAKS <- c(2, 7, 13, 19, 40, 73.0001); BAND_LABELS <- c("2-7", "7-13", "13-19", "19-40", "40-73")

select_visit <- function(cl, raw, rule = Sys.getenv("VISIT_RULE", "fill_gaps"), group = "locus") {
  x <- data.table::copy(cl[!is.na(cl[[raw]]) & !is.na(BMI_raw) & age >= 2 & age <= 73])
  x[, band := as.character(cut(age, BAND_BREAKS, labels = BAND_LABELS, right = FALSE))]
  lat <- x[order(id, -age)][, .SD[1], by = id]
  if (rule == "latest") return(lat[, rule_changed := FALSE][])
  bands <- unique(x[, .(id, grp = as.character(get(group)), band)])
  nb <- bands[, .(nb = .N), by = .(id, grp)]
  cnt <- list(); key <- function(g, b) paste(g, b, sep = "|"); inc <- function(g, b) cnt[[key(g, b)]] <<- (if (is.null(cnt[[key(g, b)]])) 0 else cnt[[key(g, b)]]) + 1
  getc <- function(g, b) if (is.null(cnt[[key(g, b)]])) 0 else cnt[[key(g, b)]]
  choice <- list()
  for (i in nb[nb == 1, which = TRUE]) { r <- nb[i]; b <- bands[id == r$id, band]; inc(r$grp, b); choice[[r$id]] <- b }
  multi <- nb[nb > 1][order(nb, id)]
  for (i in seq_len(nrow(multi))) {
    r <- multi[i]; cb <- bands[id == r$id, band]
    cb <- cb[order(match(cb, BAND_LABELS), decreasing = TRUE)]               # oldest band first, so ties go to the oldest
    b <- cb[which.min(sapply(cb, function(bb) getc(r$grp, bb)))]; inc(r$grp, b); choice[[r$id]] <- b
  }
  ch <- data.table::data.table(id = names(choice), band_sel = unlist(choice))
  # JS 2026-11-19 (Jon: points clustered): for people placed by step 4, taking the latest visit of the band piled their ages at the
  # upper band edges (7, 13, 19 y); they now get a random visit within the band (fixed seed). Step-3 people keep their latest visit.
  set.seed(as.integer(Sys.getenv("VISIT_SEED", "20261119")))
  sel <- merge(x, ch, by = "id")[band == band_sel]
  multi_ids <- multi$id
  pick <- sel[, if (.BY$id %in% multi_ids) .SD[sample.int(.N, 1)] else .SD[which.max(age)], by = id]
  out <- pick[, band_sel := NULL]
  out[, rule_changed := !(paste(id, age) %in% paste(lat$id, lat$age))]
  out[]
}

visit_selection_table <- function(cl, raw, group = "locus") {
  a <- select_visit(cl, raw, "latest", group)[, rule := "latest visit"]; b <- select_visit(cl, raw, "fill_gaps", group)[, rule := "fill gaps"]
  z <- rbind(a, b, fill = TRUE)[, band := factor(band, BAND_LABELS)]
  data.table::dcast(z[, .N, by = .(grp = get(group), rule, band)], grp + rule ~ band, value.var = "N", fill = 0)
}
