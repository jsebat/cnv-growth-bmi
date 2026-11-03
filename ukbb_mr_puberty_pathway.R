# UK Biobank: childhood adiposity -> puberty timing -> adult height (JS 2026-11-03). Runs on Expanse in growth_tables_20260920/ukbb_mr.
#  A  step 1 (women): age at menarche ~ PGS (childhood obesity, childhood BMI, adult BMI) + height PGS + age, by genotype.
#     Does genetic liability to childhood adiposity advance puberty?
#  B  two-step mediation (women): indirect effect of each adiposity PGS on adult height through menarche = A x (MR effect of one
#     year of later menarche on height, menarche PGS / first stage in all women), compared with the total effect of the PGS on height.
#  C  multivariable reduced form (both sexes, and women): adult height ~ childhood obesity PGS + menarche PGS + adult BMI PGS +
#     height PGS, by genotype, with the interaction vs non-carriers for each score.
# CAVEAT: the menarche weights (Guan et al. 2025, Nat Genet, PMID 40065166) were estimated in UK Biobank itself, so the menarche
# instrument is not independent of this sample (first stage inflated, MR biased toward the observational estimate). Adult BMI
# (Yengo 2018) and height (Yengo 2022) GWAS also include UK Biobank; the childhood obesity (Bradfield 2019) and childhood BMI
# (Vogelezang 2020) GWAS do not. Output: ukbb_mr_puberty_pathway_{A,B,C}.csv
suppressPackageStartupMessages({ library(data.table) })
WD <- getwd(); options(width = 220)
GM <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"; MS <- "/expanse/projects/sebat1/msacks/UKBB"
AAM <- "/expanse/projects/sebat1/jsebat/pgs_menarche/ukbb/PGS_guan_2025_ext_agemenarche_cojo_genomewide.tsv"
gm <- fread(GM, select = c("project_40980_IID", "f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", paste0("PC", 1:10), "Height", "age", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"))
setnames(gm, c("project_40980_IID", "f.eid", "sex_combined", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"), c("IID", "eid", "sex", "l16", "l22"))
gm <- gm[duplicate_drop == "no" & withdrawn_consent == "no" & ancestry == "blue" & sex %in% 1:2 & !is.na(Height) & !is.na(age) & l16 %in% 1:3 & l22 %in% 1:3]
gm <- gm[complete.cases(gm[, paste0("PC", 1:10), with = FALSE])]
gm[, locus := fcase(l16 == 1, "16p11.2 DEL", l16 == 3, "16p11.2 DUP", l22 == 1, "22q11.2 DEL", l22 == 3, "22q11.2 DUP", l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
GROUPS <- c("16p11.2 DEL", "No CNV", "16p11.2 DUP", "22q11.2 DUP")
gm <- gm[locus %in% GROUPS]
rp <- function(f, nm) { p <- fread(f)[n_chr == 22]; setnames(p, "PGS", nm); p[, c("IID", nm), with = FALSE] }
d <- merge(gm, rp(file.path(MS, "PRS_childhood_obesity/PGS_CHILDHOOD_OBESITY.TRANS_ANCESTRAL_2019.RESULTS_cojo_genomewide.tsv"), "pgs_cobes"), by = "IID")
d <- merge(d, rp(file.path(MS, "PRS_childhood_obesity/PGS_GCST90002409_childhood_BMI_cojo_genomewide.tsv"), "pgs_cbmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_BMI_Yengo_genomewide.tsv"), "pgs_bmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_Height_genomewide.tsv"), "pgs_ht"), by = "IID")
d <- merge(d, rp(AAM, "pgs_aam"), by = "IID")
d[, age_sq := age^2]; d[, sexf := factor(sex)]; pcs <- paste0("PC", 1:10, collapse = " + ")
for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi", "pgs_ht", "pgs_aam")) { r <- resid(lm(as.formula(paste(v, "~", pcs)), d)); set(d, j = v, value = r / sd(r)) }
d[, hZ := { r <- resid(lm(as.formula(paste("Height ~ sexf + age + age_sq +", pcs)), d)); r / sd(r) }]
men <- fread(file.path(MS, "comparitiveBodySize/field2714_age_at_menarche.tsv")); setnames(men, "f.eid", "eid")
inst <- grep("^f[.]2714[.]", names(men), value = TRUE); men[, menarche := NA_real_]; for (v in inst) men[is.na(menarche) & get(v) >= 8 & get(v) <= 18, menarche := as.numeric(get(v))]
w <- merge(d[sex == 2], men[!is.na(menarche), .(eid, menarche)], by = "eid")
cat("sample: all", nrow(d), " women with menarche", nrow(w), "\n"); print(w[, .N, by = locus])
print(round(cor(d[, .(pgs_cobes, pgs_cbmi, pgs_bmi, pgs_aam, pgs_ht)]), 3))
co <- function(m, term) { s <- summary(m)$coefficients; s[term, c(1, 2, 4)] }
PG <- c(pgs_cobes = "childhood obesity PGS", pgs_cbmi = "childhood BMI PGS", pgs_bmi = "adult BMI PGS")
# ---- A: step 1, menarche ~ adiposity PGS
A <- rbindlist(lapply(GROUPS, function(G) rbindlist(lapply(names(PG), function(v) { x <- w[locus == G]
  b <- co(lm(as.formula(paste("menarche ~", v, "+ pgs_ht + age")), x), v)
  pi <- if (G == "No CNV") NA_real_ else co(lm(as.formula(paste("menarche ~", v, "* I(locus == G) + pgs_ht + age")), w[locus %in% c(G, "No CNV")]), paste0(v, ":I(locus == G)TRUE"))[3]
  data.table(group = G, pgs = PG[[v]], n_women = nrow(x), years_per_sd = b[1], se = b[2], p = b[3], p_vs_noncarriers = pi) }))))
fwrite(A, file.path(WD, "ukbb_mr_puberty_pathway_A.csv")); cat("\n== A: years of menarche per s.d. of adiposity PGS (women) ==\n"); print(A)
# ---- B: two-step mediation in women
fsa <- co(lm(as.formula(paste("menarche ~ pgs_aam + pgs_ht + age +", pcs)), w), "pgs_aam")
B <- rbindlist(lapply(GROUPS, function(G) { x <- w[locus == G]
  rf <- co(lm(hZ ~ pgs_aam + pgs_ht, x), "pgs_aam"); mrb <- rf[1] / fsa[1]; mrse <- sqrt(rf[2]^2 / fsa[1]^2 + rf[1]^2 * fsa[2]^2 / fsa[1]^4)
  rbindlist(lapply(names(PG), function(v) { a <- A[group == G & pgs == PG[[v]]]
    tot <- co(lm(as.formula(paste("hZ ~", v, "+ pgs_ht")), x), v)
    ind <- a$years_per_sd * mrb; ind_se <- sqrt(a$years_per_sd^2 * mrse^2 + mrb^2 * a$se^2)
    data.table(group = G, pgs = PG[[v]], n_women = nrow(x), step1_years_per_sd = a$years_per_sd, mr_height_per_year = mrb, mr_se = mrse,
               indirect_via_menarche = ind, indirect_se = ind_se, indirect_p = 2 * pnorm(-abs(ind / ind_se)), total_height_per_sd = tot[1], total_se = tot[2], total_p = tot[3]) })) }))
fwrite(B, file.path(WD, "ukbb_mr_puberty_pathway_B.csv")); cat("\n== B: two-step mediation (women); first stage menarche PGS", round(fsa[1], 3), "y/sd ==\n"); print(B)
# ---- C: multivariable reduced form
C <- rbindlist(lapply(list(both = d, women = w), function(D) rbindlist(lapply(GROUPS, function(G) { x <- D[locus == G]
  m <- lm(hZ ~ pgs_cobes + pgs_aam + pgs_bmi + pgs_ht, x)
  rbindlist(lapply(c("pgs_cobes", "pgs_aam", "pgs_bmi"), function(v) { b <- co(m, v)
    pi <- if (G == "No CNV") NA_real_ else { z <- D[locus %in% c(G, "No CNV")]; z[, g := as.numeric(locus == G)]
      co(lm(hZ ~ (pgs_cobes + pgs_aam + pgs_bmi + pgs_ht) * g, z), paste0(v, ":g"))[3] }
    data.table(group = G, n = nrow(x), score = v, height_per_sd = b[1], se = b[2], p = b[3], p_vs_noncarriers = pi) })) }))), idcol = "sample")
C[, score := c(pgs_cobes = "childhood obesity PGS", pgs_aam = "age at menarche PGS", pgs_bmi = "adult BMI PGS")[score]]
fwrite(C, file.path(WD, "ukbb_mr_puberty_pathway_C.csv")); cat("\n== C: multivariable, adult height per s.d. of each score ==\n"); print(C)
cat("done\n")
