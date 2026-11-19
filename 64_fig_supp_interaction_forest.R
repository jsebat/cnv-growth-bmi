# Supplementary figure: forest plot of the two-way (CNV x BMI) and three-way (CNV x BMI x autism polygenic score) interactions on
# SRS and SCQ, pooled and by CNV, by cohort and by site (JS / Claude 2026-11-14). Combined sample of 60 (G2MH + Simons VIP +
# Searchlight + SSC carriers; Box-Cox scale; reference = SSC unaffected relatives).
#  two-way : y ~ is_cnv + s(age, by = is_cnv) + s(age, by = BMI) + cnvBMI [+ carrier cohort:CNV] (as Fig. 5k-t; 61 vs_ref)
#  three-way: y ~ is_cnv * (poly(age, 2) + male) + carrier cohort:CNV + is_cnv * BMI * pgs, robust SE (62); genotyped carriers only
#             (G2MH, Simons VIP genomes, SSC); April 2026 autism score (lo-SSC weights for SSC), residualised on 10 PCs within cohort.
# Sites: G2MH site_id (7 sites); the Simons VIP release carries no site field and Searchlight is a remote registry, so each is one
# group. Subgroups with fewer than 15 carriers are omitted. Output: fig_supp_interaction_forest.{png,pdf} + _table.csv (aggregate).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s61 <- readLines(file.path(HERE, "61_combined_cnv_x_bmi.R")); eval(parse(text = s61[1:(grep("^out <- list", s61) - 1)]))
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; AN <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data")
suppressPackageStartupMessages({ library(ggplot2) })
st <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"), select = c("guid", "site_id"))[site_id != ""][, .SD[1], by = guid]
SITE <- c(penn_chop = "Penn/CHOP", card = "Cardiff", ucla = "UCLA", mont = "Montreal", toro = "Toronto", mass = "Maastricht", leuv = "Leuven")
C <- merge(C, st[, .(id = guid, site_id)], by = "id", all.x = TRUE)
C[, site := fcase(cohort == "G2MH" & !is.na(site_id), paste("G2MH", SITE[site_id]), cohort == "G2MH", "G2MH site unknown", default = cohort)]
s62 <- readLines(file.path(HERE, "62_combined_threeway.R")); eval(parse(text = s62[grep("^rob <- function", s62):grep("^X\\[, cg :=", s62)]))
X[site == "Searchlight", site := "SVIP"]   # genotyped Searchlight records are SVIP-family genomes: shown with Simons VIP
three <- function(car, ref) { if (nrow(car) < 15) return(NULL); z <- rbind(car[, .(y, BMI, pgs, age, male, cg, is_cnv)], ref[, .(y, BMI, pgs, age, male, cg, is_cnv)])
  f <- if (uniqueN(car$cg) > 1) y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI * pgs else y ~ is_cnv * (poly(age, 2) + male) + is_cnv * BMI * pgs
  r <- tryCatch(rob(lm(f, z), "is_cnv:BMI:pgs"), error = function(e) c(b = NA, se = NA, p = NA)); c(n = nrow(car), r) }
LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
COH <- c("G2MH", "SVIP", "Searchlight", "SSC"); COHL <- c(G2MH = "G2MH", SVIP = "Simons VIP", Searchlight = "Simons Searchlight", SSC = "SSC carriers")
SITES <- c(paste("G2MH", SITE), "SVIP", "Searchlight")
rows <- list(); addr <- function(ph, sec, lab, typ, r) if (!is.null(r)) rows[[length(rows) + 1]] <<- data.table(trait = ph, section = sec, label = lab, type = typ, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
for (ph in c("SRS", "SCQ")) {
  ref2 <- C[pheno == ph & locus == "SSC relative"]; ref3 <- X[pheno == ph & locus == "SSC relative"]
  car2 <- C[pheno == ph & is_cnv == 1]; car3 <- X[pheno == ph & is_cnv == 1]
  grp <- list(list("Pooled", "All carriers", car2, car3))
  for (L in LOC) grp[[length(grp) + 1]] <- list("By CNV", L, car2[locus == L], car3[locus == L])
  for (k in COH) grp[[length(grp) + 1]] <- list("By cohort", COHL[[k]], car2[cohort == k], car3[cohort == k])
  for (s in SITES) grp[[length(grp) + 1]] <- list("By site", if (s %in% c("SVIP", "Searchlight")) paste(COHL[[s]], "(no site field)") else s, car2[site == s], car3[site == s])
  for (g in grp) { if (nrow(g[[3]]) >= 15) addr(ph, g[[1]], g[[2]], "Two-way: CNV x BMI", vs_ref(g[[3]], ref2)); addr(ph, g[[1]], g[[2]], "Three-way: CNV x BMI x autism score", three(g[[4]], ref3)) }
  cat(ph, "done\n")
}
R <- rbindlist(rows); fwrite(R, file.path(HERE, "fig_supp_interaction_forest_table.csv"))
options(width = 200); print(R[, .(trait, section, label, type = substr(type, 1, 9), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
# ---------------------------------------------------------------- figure
R <- R[!is.na(b)]
ordl <- c("All carriers", LOC, COHL, paste("G2MH", SITE), "Simons VIP (no site field)", "Simons Searchlight (no site field)")
R[, label := factor(label, levels = rev(ordl))][, section := factor(section, c("Pooled", "By CNV", "By cohort", "By site"))][, trait := factor(trait, c("SRS", "SCQ"))]
R[, type := factor(type, c("Two-way: CNV x BMI", "Three-way: CNV x BMI x autism score"))]
R[, nl := sprintf("%d", n)]
OC <- c(SRS = "#C0392B", SCQ = "#2E86C1"); pd <- position_dodge(width = 0.6)
fig <- ggplot(R, aes(b, label, colour = trait, group = trait)) + geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd, linewidth = 0.55) +
  geom_point(aes(shape = p < 0.05), size = 2.1, position = pd) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("P >= 0.05", "P < 0.05"), name = NULL) +
  facet_grid(section ~ type, scales = "free", space = "free_y") + scale_colour_manual(values = OC, name = NULL) +
  coord_cartesian(clip = "off") +
  labs(x = "Interaction estimate (95% CI); two-way: s.d. of trait per s.d. of BMI, carriers minus SSC relatives; three-way: per s.d. of autism score", y = NULL,
       title = "CNV x BMI and CNV x BMI x autism polygenic score interactions on autistic traits",
       subtitle = "Combined sample (G2MH, Simons VIP, Simons Searchlight, SSC carriers) versus SSC unaffected relatives; Box-Cox scale; three-way: genotyped carriers only") +
  theme_classic(base_size = 10.5) + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0), strip.text.x = element_text(face = "bold"),
    panel.spacing.y = unit(0.6, "lines"), legend.position = "top", plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.5, colour = "grey30"))
ggsave(file.path(HERE, "fig_supp_interaction_forest.png"), fig, width = 11, height = 9, dpi = 250, bg = "white")
ggsave(file.path(HERE, "fig_supp_interaction_forest.pdf"), fig, width = 11, height = 9)
