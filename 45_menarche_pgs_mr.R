# Age-at-menarche polygenic score in the clinical cohorts (JS / Claude 2026-11-02; reviewer MR request, puberty).
# Score: Guan 2025 'ext' age at menarche, SBayesRC/COJO weights from James (pgs_library_2026-11-02/age_menarche), plink2 --score
# on Expanse (jobs 55513232: SSC genotypes and G2MH WGS -> jsebat/pgs_menarche/scores; raw .sscore copied to
# Jon/Claude/CNV_growth_trajectories_data/menarche_pgs). Residualised on 10 ancestry PCs and z-scored within cohort, as in 30a.
# Higher score = later menarche (checked below against height at 10-18 in controls: later puberty -> shorter in adolescence).
#  A. outcome ~ pgs_aam x carrier + sex + age + pgs_ht, by outcome band and carrier group (both sexes; then females only)
#  B. childhood obesity PGS effect in carriers with and without the menarche score as covariate (does puberty timing carry it?)
# Setup (scores, phenotype table U) reused from the reduced-form script (mr_reduced_form.R in this folder, up to its analysis block).
# Outputs: mr_menarche_pgs_A.csv, mr_menarche_pgs_B.csv (aggregate statistics only).
SCR <- file.path("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories", "mr_reduced_form.R")
src <- readLines(SCR); stop_at <- grep("res <- list()", src, fixed = TRUE)[1] - 1
eval(parse(text = src[1:stop_at]))
MP <- file.path(DB, "Jon/Claude/CNV_growth_trajectories_data/menarche_pgs")
rdm <- function(f) { x <- fread(file.path(MP, f)); setnames(x, "#IID", "IID"); x[, .(IID, aam = SCORE1_AVG)] }
# G2MH WGS
gm <- rdm("g2mh_wgs.age_menarche_guan2025.sscore")[, GUID := sub("^[0-9]+-", "", IID)]
gm <- merge(gm, pcw, by = "GUID"); gm[, pgs_aam := rz(aam, gm)]
# SSC
sm <- rdm("ssc.age_menarche_guan2025.sscore")
pcs <- fread(file.path(DB, "G2MH/NetworkPaper/AnalysisByAnjali/data/all_cohorts_adsi_pgs.tsv"), select = c("sample_id", "source", paste0("PC", 1:10)))[source == "ssc"]
sm <- merge(sm, pcs[, c("sample_id", paste0("PC", 1:10)), with = FALSE], by.x = "IID", by.y = "sample_id")
sm[, pgs_aam := as.numeric(scale(residuals(lm(as.formula(paste("aam ~", paste0("PC", 1:10, collapse = " + "))), data = sm))))]
san <- fread(file.path(DB, "Autism/SSC_Version_15_Phenotype_Data_Set_9/SSC IDs Sanders.csv"))
map <- unique(rbind(san[, .(IID = SS_RUID, individual)], san[, .(IID = `RUID-1`, individual)])[IID != "." & IID != ""])
sm <- merge(sm[, .(IID, pgs_aam)], map, by = "IID")[, .SD[1], by = individual]
PA <- rbind(gm[, .(key = GUID, pgs_aam)], sm[, .(key = paste0("SSC:", individual), pgs_aam)])
cat("menarche PGS: G2MH", nrow(gm), " SSC", nrow(sm), "\n")
U <- merge(U, PA, by = "key", all.x = FALSE)
cat("phenotype rows with menarche PGS:", nrow(U), "; individuals: carriers", U[carrier == TRUE, uniqueN(key)], " controls", U[carrier == FALSE, uniqueN(key)], "\n")
print(U[, .(individuals = uniqueN(key)), by = .(locus, sex)][order(locus, sex)])
cat("correlation of menarche PGS with other scores (one row per person):\n")
print(round(cor(unique(U[, .(key, pgs_aam, pgs_cobes, pgs_cbmi, pgs_bmi, pgs_ht)])[, -1]), 2))

bands <- list(`Height 2-9` = quote(pheno == "Height" & age >= 2 & age < 10), `Height 10-14` = quote(pheno == "Height" & age >= 10 & age < 15),
              `Height 15-18` = quote(pheno == "Height" & age >= 15 & age < 19), `Height 10-18` = quote(pheno == "Height" & age >= 10 & age < 19),
              `Height adult` = quote(pheno == "Height" & age >= 19), `BMI 10-18` = quote(pheno == "Height" & age >= 10 & age < 19),
              `Head circumference child` = quote(pheno == "Head circumference" & age < 18))
GROUPS <- c("all carriers", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
fitA <- function(d, outcome_bmi, sexes) {
  d <- d[sex %in% sexes]; if (outcome_bmi) d[, y := BMI]
  if (d[carrier == TRUE, .N] < 15 || d[carrier == FALSE, .N] < 15) return(NULL)
  cov <- if (length(sexes) > 1) "+ sex" else ""
  m <- lm(as.formula(paste("y ~ pgs_aam * carrier", cov, "+ age + pgs_ht")), data = d); s <- summary(m)$coefficients
  mc <- summary(lm(as.formula(paste("y ~ pgs_aam", cov, "+ age + pgs_ht")), data = d[carrier == TRUE]))$coefficients["pgs_aam", ]
  mk <- summary(lm(as.formula(paste("y ~ pgs_aam", cov, "+ age + pgs_ht")), data = d[carrier == FALSE]))$coefficients["pgs_aam", ]
  data.table(n_carriers = d[carrier == TRUE, .N], n_controls = d[carrier == FALSE, .N],
             slope_controls = mk[1], p_controls = mk[4], slope_carriers = mc[1], se_carriers = mc[2], p_carriers = mc[4],
             interaction = s["pgs_aam:carrierTRUE", 1], int_se = s["pgs_aam:carrierTRUE", 2], int_p = s["pgs_aam:carrierTRUE", 4])
}
resA <- list()
for (sx in list(c("male", "female"), "female")) for (nm in names(bands)) for (G in GROUPS) {
  d <- U[eval(bands[[nm]])]; d <- if (G == "all carriers") d else d[locus == G | !carrier]
  r <- fitA(copy(d), grepl("^BMI", nm), sx)
  if (!is.null(r)) resA[[length(resA) + 1]] <- cbind(data.table(sexes = paste(sx, collapse = "+"), outcome = nm, group = G), r)
}
RA <- rbindlist(resA); fwrite(RA, file.path(HERE, "mr_menarche_pgs_A.csv"))
cat("\n== A. menarche PGS (per s.d.; + = later menarche) ==\n")
print(RA[, .(sexes, outcome, group, n_car = n_carriers, n_ctl = n_controls, ctl = round(slope_controls, 3), p_ctl = signif(p_controls, 2),
             car = round(slope_carriers, 3), p_car = signif(p_carriers, 2), int = round(interaction, 3), int_p = signif(int_p, 2))], nrows = 300)

# B. childhood obesity PGS in carriers, with and without the menarche score
resB <- list()
for (nm in c("Height 2-9", "Height 10-14", "Height 15-18", "Height 10-18", "Height adult", "Head circumference child")) for (G in GROUPS) for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi")) {
  d <- U[eval(bands[[nm]]) & carrier == TRUE]; if (G != "all carriers") d <- d[locus == G]
  if (nrow(d) < 15) next
  b0 <- summary(lm(as.formula(paste("y ~", v, "+ sex + age + pgs_ht")), data = d))$coefficients[v, ]
  b1 <- summary(lm(as.formula(paste("y ~", v, "+ pgs_aam + sex + age + pgs_ht")), data = d))$coefficients
  resB[[length(resB) + 1]] <- data.table(outcome = nm, group = G, pgs = v, n = nrow(d), beta = b0[1], p = b0[4],
                                         beta_adj_aam = b1[v, 1], p_adj_aam = b1[v, 4], beta_aam = b1["pgs_aam", 1], p_aam = b1["pgs_aam", 4])
}
RB <- rbindlist(resB); fwrite(RB, file.path(HERE, "mr_menarche_pgs_B.csv"))
cat("\n== B. obesity / BMI PGS in carriers, without and with the menarche PGS ==\n")
print(RB[, .(outcome, group, pgs, n, beta = round(beta, 3), p = signif(p, 2), beta_adj = round(beta_adj_aam, 3), p_adj = signif(p_adj_aam, 2), aam = round(beta_aam, 3), p_aam = signif(p_aam, 2))], nrows = 300)
