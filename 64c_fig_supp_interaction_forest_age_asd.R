# Supplementary figure: two-way (CNV x BMI) and three-way (CNV x BMI x autism score) interactions on SRS / SCQ stratified by age
# and by autism diagnosis (JS / Claude 2026-11-14). Same data and models as 64 (combined G2MH + Simons VIP + Searchlight + SSC
# carriers, Box-Cox scale; three-way in genotyped carriers). Reference = SSC unaffected relatives of matching age: designated
# siblings for carriers under 18, parents for adults (SRS only; parents have no SCQ); all relatives for the all-age rows.
# Age bands: 4-11, 12-17, 18+ y. Autism diagnosis as in 63 (G2MH subjects.csv; SVIP diagnosis_summary; Searchlight clinical or
# reported ASD; SSC probands). Strata with < 15 carriers omitted. Output: fig_supp_interaction_forest_age_asd.{png,pdf} + _table.csv.
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s64 <- readLines(file.path(HERE, "64_fig_supp_interaction_forest.R")); s63 <- readLines(file.path(HERE, "63_combined_noasd_sensitivity.R"))
eval(parse(text = s64[1:(grep("^s62 <- ", s64) - 1)]))
eval(parse(text = s63[grep("^sv <- fread", s63):(grep("^cat\\(\"autism diagnosis among", s63) - 1)]))
eval(parse(text = s64[grep("^s62 <- ", s64):grep("^COH <- ", s64)]))
SIBf <- function(d) d[locus == "SSC relative" & grepl("[.]s[0-9]+$", id)]; PARf <- function(d) d[locus == "SSC relative" & grepl("[.](mo|fa)$", id)]; RELf <- function(d) d[locus == "SSC relative"]
rows <- list(); addr <- function(ph, sec, lab, typ, r, nref) if (!is.null(r)) rows[[length(rows) + 1]] <<- data.table(trait = ph, section = sec, label = lab, type = typ, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]], n_ref = nref)
BANDS <- list(`4-11 y` = c(4, 12), `12-17 y` = c(12, 18), `18+ y` = c(18, 99))
for (ph in c("SRS", "SCQ")) {
  C2 <- C[pheno == ph]; X3 <- X[pheno == ph]
  car2 <- C2[is_cnv == 1]; car3 <- X3[is_cnv == 1]
  addr(ph, "All ages", "All carriers", "Two-way: CNV x BMI", vs_ref(car2, RELf(C2)), nrow(RELf(C2))); addr(ph, "All ages", "All carriers", "Three-way: CNV x BMI x autism score", three(car3, RELf(X3)), nrow(RELf(X3)))
  for (bn in names(BANDS)) { a <- BANDS[[bn]]; RF <- if (a[1] >= 18) PARf else SIBf
    r2 <- RF(C2)[age >= a[1] & age < a[2]]; r3 <- RF(X3)[age >= a[1] & age < a[2]]; if (nrow(r2) < 50) next
    c2 <- car2[age >= a[1] & age < a[2]]; c3 <- car3[age >= a[1] & age < a[2]]; lab <- sprintf("%s (vs SSC %s)", bn, if (a[1] >= 18) "parents" else "siblings")
    if (nrow(c2) >= 15) addr(ph, "By age", lab, "Two-way: CNV x BMI", vs_ref(c2, r2), nrow(r2)); addr(ph, "By age", lab, "Three-way: CNV x BMI x autism score", three(c3, r3), nrow(r3)) }
  for (k in c(0L, 1L)) { lab <- if (k == 1L) "Autism diagnosis" else "No autism diagnosis"
    c2 <- car2[asd2 %in% k]; c3 <- car3[asd2 %in% k]
    if (nrow(c2) >= 15) addr(ph, "By autism diagnosis", lab, "Two-way: CNV x BMI", vs_ref(c2, RELf(C2)), nrow(RELf(C2))); addr(ph, "By autism diagnosis", lab, "Three-way: CNV x BMI x autism score", three(c3, RELf(X3)), nrow(RELf(X3))) }
  for (k in c(0L, 1L)) for (bn in c("< 18 y", "18+ y")) { lab <- sprintf("%s, %s", if (k == 1L) "Autism" else "No autism", bn); ad <- bn == "18+ y"; RF <- if (ad) PARf else SIBf
    r2 <- RF(C2)[if (ad) age >= 18 else age < 18]; r3 <- RF(X3)[if (ad) age >= 18 else age < 18]; if (nrow(r2) < 50) next
    c2 <- car2[asd2 %in% k & (if (ad) age >= 18 else age < 18)]; c3 <- car3[asd2 %in% k & (if (ad) age >= 18 else age < 18)]
    if (nrow(c2) >= 15) addr(ph, "Autism x age", lab, "Two-way: CNV x BMI", vs_ref(c2, r2), nrow(r2)); addr(ph, "Autism x age", lab, "Three-way: CNV x BMI x autism score", three(c3, r3), nrow(r3)) }
  cat(ph, "done\n")
}
R <- rbindlist(rows); fwrite(R, file.path(HERE, "fig_supp_interaction_forest_age_asd_table.csv"))
options(width = 200); print(R[, .(trait, section, label, type = substr(type, 1, 9), n, n_ref, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
R <- R[!is.na(b)]
ordl <- c("All carriers", sprintf("%s (vs SSC %s)", names(BANDS), c("siblings", "siblings", "parents")), "No autism diagnosis", "Autism diagnosis", "No autism, < 18 y", "No autism, 18+ y", "Autism, < 18 y", "Autism, 18+ y")
R[, label := factor(label, levels = rev(ordl))][, section := factor(section, c("All ages", "By age", "By autism diagnosis", "Autism x age"))][, trait := factor(trait, c("SRS", "SCQ"))]
R[, type := factor(type, c("Two-way: CNV x BMI", "Three-way: CNV x BMI x autism score"))]
OC <- c(SRS = "#C0392B", SCQ = "#2E86C1"); pd <- position_dodge(width = 0.6)
fig <- ggplot(R, aes(b, label, colour = trait, group = trait)) + geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd, linewidth = 0.55) +
  geom_point(aes(shape = p < 0.05), size = 2.1, position = pd) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("P >= 0.05", "P < 0.05"), name = NULL) +
  geom_text(aes(x = Inf, label = paste0("n=", n), group = trait), position = pd, hjust = 1.05, size = 2.4, show.legend = FALSE) +
  facet_grid(section ~ type, scales = "free", space = "free_y") + scale_colour_manual(values = OC, name = NULL) +
  labs(x = "Interaction estimate (95% CI); two-way: s.d. of trait per s.d. of BMI, carriers minus SSC relatives; three-way: per s.d. of autism score", y = NULL,
       title = "CNV x BMI and CNV x BMI x autism polygenic score interactions, by age and autism diagnosis",
       subtitle = "Combined sample; reference SSC siblings (< 18 y) or parents (18+ y; SRS only), all relatives for all-age rows; Box-Cox scale; n = carriers") +
  theme_classic(base_size = 10.5) + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0), strip.text.x = element_text(face = "bold"),
    panel.spacing.y = unit(0.6, "lines"), legend.position = "top", plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.5, colour = "grey30"))
ggsave(file.path(HERE, "fig_supp_interaction_forest_age_asd.png"), fig, width = 12, height = 7.5, dpi = 250, bg = "white")
ggsave(file.path(HERE, "fig_supp_interaction_forest_age_asd.pdf"), fig, width = 12, height = 7.5)
