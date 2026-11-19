# CNV x BMI x autism polygenic score on SRS / SCQ in the combined genotyped carrier sample: G2MH (WGS) + Simons VIP (WGS of the
# SVIP families, Anjali's searchlight_393 set) + SSC carriers (arrays) (JS / Claude 2026-11-14; Jon: G2MH + SVIP + SSC genotypes
# until the Searchlight WGS arrives). Autism score: April 2026 ASD GWAS, full SBayesRC weights for G2MH and SVIP (not GWAS cohorts;
# Expanse pgs_obesity_ssc/scores/{g2mh_wgs,svip_wgs}.asd2026_full.sscore), leave-SSC-out weights for SSC; each residualised on 10
# PCs and z-scored within cohort. Phenotypes from the combined table (60; Box-Cox scale). Reported for the combined sample and
# per cohort; three-way test against SSC unaffected relatives and SSC autistic probands (lo-SSC score).
# Models: within carriers y ~ cohort:CNV + poly(age, 2) + male + BMI * pgs + cohort x BMI + cohort x pgs (robust SE);
# three-way y ~ is_cnv * (poly(age,2) + male) + carrier cohort:CNV + is_cnv * BMI * pgs. Output: combined_threeway.csv.
suppressPackageStartupMessages({ library(data.table); library(sandwich) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DB <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"; DATA <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data"); AN <- file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data")
C <- readRDS(file.path(DATA, "combined_srs_scq.rds")); LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
rob <- function(m, v) { V <- vcovHC(m, type = "HC1"); b <- unname(coef(m)[v]); s <- sqrt(V[v, v]); c(b = b, se = s, p = 2 * pnorm(-abs(b / s))) }
rzp <- function(v, P) as.numeric(scale(residuals(lm(as.formula(paste("v ~", paste(names(P), collapse = " + "))), data = cbind(P, v = v)))))
# G2MH
g <- fread(file.path(DATA, "g2mh_scores_raw/g2mh_wgs.asd2026_full.sscore")); setnames(g, "#IID", "IID"); g[, id := sub("^[0-9]+-", "", IID)]
pcw <- fread(file.path(AN, "pca_results/naive/cohort.sscore")); setnames(pcw, "#IID", "IID"); pcw[, id := sub("^[0-9]+-", "", IID)]
g <- merge(g[, .(id, s = SCORE1_AVG)], pcw[, c("id", paste0("PC", 1:10, "_AVG")), with = FALSE], by = "id"); g[, pgs := rzp(s, g[, paste0("PC", 1:10, "_AVG"), with = FALSE])]
# SVIP
sv <- fread(file.path(DATA, "g2mh_scores_raw/svip_wgs.asd2026_full.sscore")); setnames(sv, "#IID", "IID")
pcs <- fread(file.path(AN, "searchlight/pgs/analysis_dataset.tsv"), select = c("IID", paste0("GLOBAL_PC", 1:10)))
meta <- unique(fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/Searchlight.WGS2.2023_05.sample_metadata.tsv"))[, .(searchlight_id, sfari_iid)])
sv <- merge(merge(sv[, .(IID, s = SCORE1_AVG)], pcs, by = "IID"), meta, by.x = "IID", by.y = "searchlight_id"); sv[, pgs := rzp(s, sv[, paste0("GLOBAL_PC", 1:10), with = FALSE])][, id := sub("-", ".", sfari_iid, fixed = TRUE)]
# SSC (lo-SSC, already residualised and z-scored within SSC by 30a)
ss <- fread(file.path(DATA, "ssc_pgs_scores.tsv"))[, .(id = individual, pgs = pgs_asd)]
S <- rbind(g[, .(id, pgs, src = "G2MH")], sv[, .(id, pgs, src = "SVIP")], ss[, .(id, pgs, src = "SSC")])
X <- merge(C, S, by = "id")[(cohort == src) | (cohort == "Searchlight" & src == "SVIP")]
X[cohort == "Searchlight", cohort := "SVIP"]   # Searchlight records of SVIP-family genomes
X[, cg := fifelse(is_cnv == 1, paste(cohort, locus, sep = ":"), "ref")]
cat("genotyped sample:\n"); print(dcast(X[, .N, by = .(pheno, cohort, locus)], pheno + locus ~ cohort, value.var = "N", fill = 0))
# SSC autistic probands with the lo-SSC score (reference 2), from the Fig. 5 table built by 24 (not in C): read from 52's inputs
out <- list(); add <- function(ph, block, spec, r, n) out[[length(out) + 1]] <<- data.table(trait = ph, block = block, spec = spec, n = n, b = r[["b"]], se = r[["se"]], p = r[["p"]])
for (ph in c("SRS", "SCQ")) {
  car <- X[pheno == ph & is_cnv == 1]
  add(ph, "within carriers", "combined: G2MH + SVIP + SSC carriers", rob(lm(y ~ cg + poly(age, 2) + male + BMI * pgs + cohort:BMI + cohort:pgs, car), "BMI:pgs"), nrow(car))
  add(ph, "within carriers", "combined, without cohort x BMI and cohort x score terms", rob(lm(y ~ cg + poly(age, 2) + male + BMI * pgs, car), "BMI:pgs"), nrow(car))
  for (coh in c("G2MH", "SVIP", "SSC")) { d <- car[cohort == coh]; if (nrow(d) < 15) { cat(ph, coh, "n =", nrow(d), "(too few)\n"); next }
    add(ph, "within carriers", paste("cohort:", coh), rob(lm(as.formula(paste("y ~", if (uniqueN(d$cg) > 1) "cg +" else "", "poly(age, 2) + male + BMI * pgs")), d), "BMI:pgs"), nrow(d)) }
  m1 <- lm(y ~ cg + poly(age, 2) + male + cohort * BMI * pgs, car); m0 <- lm(y ~ cg + poly(age, 2) + male + cohort * BMI + cohort * pgs + BMI:pgs, car)
  add(ph, "within carriers", "cohort heterogeneity of BMI x score (F test)", c(b = NA, se = NA, p = anova(m0, m1)$`Pr(>F)`[2]), nrow(car))
  for (L in LOC) { d <- car[locus == L]; if (nrow(d) >= 15) add(ph, "within carriers, combined cohorts", paste("CNV:", L), rob(lm(as.formula(paste("y ~", if (uniqueN(d$cg) > 1) "cg +" else "", "poly(age, 2) + male + BMI * pgs")), d), "BMI:pgs"), nrow(d)) }
  med <- median(car$pgs); for (h in c("below", "above")) { d <- car[if (h == "below") pgs < med else pgs >= med]
    add(ph, "within carriers, combined cohorts", paste("BMI slope, autism score", h, "median"), rob(lm(y ~ cg + poly(age, 2) + male + BMI, d), "BMI"), nrow(d)) }
  ref <- X[pheno == ph & locus == "SSC relative"]
  z <- rbind(car[, .(y, BMI, pgs, age, male, cg, is_cnv)], ref[, .(y, BMI, pgs, age, male, cg, is_cnv)])
  add(ph, "three-way", "CNV x BMI x score vs SSC relatives", rob(lm(y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI * pgs, z), "is_cnv:BMI:pgs"), nrow(z))
  add(ph, "three-way", "BMI x score in SSC relatives", rob(lm(y ~ poly(age, 2) + male + BMI * pgs, ref), "BMI:pgs"), nrow(ref))
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "combined_threeway.csv"))
options(width = 220); print(R[, .(trait, block, spec, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 100)
