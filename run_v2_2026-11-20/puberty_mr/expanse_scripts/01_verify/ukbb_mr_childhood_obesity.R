# UK Biobank one-sample MR / genetic-exposure test (JS 2026-11-02, reviewer request): does genetic liability to CHILDHOOD
# obesity (vs ADULT BMI) predict adult height differently in 16p11.2 / 22q11.2 CNV carriers than in non-carriers?
# Exposure instruments: childhood obesity PGS (Bradfield 2019 trans-ancestral, SBayesRC), childhood BMI PGS (GCST90002409),
# adult BMI PGS (Yengo 2018); height PGS as pleiotropy covariate. Scores from Molly's folders (PRS_childhood_obesity,
# comparitiveBodySize), residualised on PC1-10 within the EUR (blue) sample and z-scored.
# Outcome: adult height residual (Height ~ sex + age + age^2 + PC1-10, divided by its SD) as in Molly's 06 script.
# Tests: (1) first stage, step10 (body size at 10: -1/0/+1) ~ PGS and BMI z ~ PGS, by CNV group; (2) reduced form,
# hZ ~ PGS x group (+ height PGS); (3) within-carrier slopes; (4) Wald ratio = reduced form / first stage (delta SE).
# Sample: blue ancestry, QC, height+age+PCs; the 1687 answer is needed only for the step10 first stage.
suppressPackageStartupMessages({ library(data.table) })
WD <- getwd(); options(width = 220)
GM <- "/expanse/projects/sebat1/msacks/UKBB/2024/genomatrix_UKBB_2025.csv"; MS <- "/expanse/projects/sebat1/msacks/UKBB"
gm <- fread(GM, select = c("project_40980_IID", "f.eid", "sex_combined", "duplicate_drop", "withdrawn_consent", "ancestry", paste0("PC", 1:10), "Height", "BMI_raw", "age", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"))
setnames(gm, c("project_40980_IID", "f.eid", "sex_combined", "locus_16p11.2_proximal", "locus_22q11.2_BP1_BP4_A_D"), c("IID", "eid", "sex", "l16", "l22"))
gm <- gm[duplicate_drop == "no" & withdrawn_consent == "no" & ancestry == "blue" & sex %in% 1:2 & !is.na(Height) & !is.na(age) & l16 %in% 1:3 & l22 %in% 1:3]
gm <- gm[complete.cases(gm[, paste0("PC", 1:10), with = FALSE])]
gm[, locus := fcase(l16 == 1, "16p11.2 DEL", l16 == 3, "16p11.2 DUP", l22 == 1, "22q11.2 DEL", l22 == 3, "22q11.2 DUP", l16 == 2 & l22 == 2, "No CNV", default = NA_character_)]
gm <- gm[!is.na(locus)]
rp <- function(f, nm) { p <- fread(f)[n_chr == 22]; setnames(p, "PGS", nm); p[, c("IID", nm), with = FALSE] }
d <- merge(gm, rp(file.path(MS, "PRS_childhood_obesity/PGS_CHILDHOOD_OBESITY.TRANS_ANCESTRAL_2019.RESULTS_cojo_genomewide.tsv"), "pgs_cobes"), by = "IID")
d <- merge(d, rp(file.path(MS, "PRS_childhood_obesity/PGS_GCST90002409_childhood_BMI_cojo_genomewide.tsv"), "pgs_cbmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_BMI_Yengo_genomewide.tsv"), "pgs_bmi"), by = "IID")
d <- merge(d, rp(file.path(MS, "comparitiveBodySize/PGS_Height_genomewide.tsv"), "pgs_ht"), by = "IID")
bs <- fread(file.path(MS, "comparitiveBodySize/field1687_comparative_body_size_age10.tsv")); setnames(bs, "f.eid", "eid")
inst <- grep("^f\\.1687\\.", names(bs), value = TRUE); bs[, code := NA_integer_]; for (v in inst) bs[is.na(code) & get(v) %in% 1:3, code := as.integer(get(v))]
bs <- bs[!is.na(code)]; bs[, step10 := fcase(code == 1L, -1, code == 3L, 0, code == 2L, 1)]
d <- merge(d, bs[, .(eid, step10)], by = "eid", all.x = TRUE)
cat("UKBB EUR with height, PCs and all four PGS:", nrow(d), "\n"); print(d[, .(n = .N, n_step10 = sum(!is.na(step10)), n_bmi = sum(!is.na(BMI_raw))), by = locus])
d[, age_sq := age^2]; d[, sexf := factor(sex)]; pcs <- paste0("PC", 1:10, collapse = " + ")
for (v in c("pgs_cobes", "pgs_cbmi", "pgs_bmi", "pgs_ht")) { r <- resid(lm(as.formula(paste(v, "~", pcs)), d)); set(d, j = v, value = r / sd(r)) }
d[, hZ := { r <- resid(lm(as.formula(paste("Height ~ sexf + age + age_sq +", pcs)), d)); r / sd(r) }]
d[, bZ := NA_real_]; ib <- which(!is.na(d$BMI_raw)); rb <- resid(lm(as.formula(paste("log(BMI_raw) ~ sexf + age + age_sq +", pcs)), d[ib])); set(d, i = ib, j = "bZ", value = rb / sd(rb))
d[, carrier := locus != "No CNV"]
LOCI <- c("16p11.2 DEL", "16p11.2 DUP", "22q11.2 DEL", "22q11.2 DUP"); PG <- c(pgs_cobes = "childhood obesity PGS", pgs_cbmi = "childhood BMI PGS", pgs_bmi = "adult BMI PGS")
co <- function(m, term) { s <- summary(m)$coefficients; if (term %in% rownames(s)) s[term, c(1, 2, 4)] else c(NA, NA, NA) }
wr <- function(rf, fs) { b <- rf[1] / fs[1]; se <- sqrt(rf[2]^2 / fs[1]^2 + rf[1]^2 * fs[2]^2 / fs[1]^4); c(b, se, 2 * pnorm(-abs(b / se))) }
out <- list()
for (G in c(LOCI, "all carriers")) for (v in names(PG)) {
  x <- if (G == "all carriers") copy(d) else d[locus %in% c(G, "No CNV")]; x[, cnv := as.numeric(carrier)]
  nc <- x[cnv == 1, .N]; if (nc < 5) next
  f10c <- co(lm(as.formula(paste("step10 ~", v)), x[cnv == 1]), v); f10n <- co(lm(as.formula(paste("step10 ~", v)), x[cnv == 0]), v)
  fbc <- co(lm(as.formula(paste("bZ ~", v)), x[cnv == 1]), v);      fbn <- co(lm(as.formula(paste("bZ ~", v)), x[cnv == 0]), v)
  r0 <- lm(as.formula(paste("hZ ~", v, "* cnv")), x); r1 <- lm(as.formula(paste("hZ ~", v, "* cnv + pgs_ht * cnv")), x)
  wc <- co(lm(as.formula(paste("hZ ~", v, "+ pgs_ht")), x[cnv == 1]), v); wn <- co(lm(as.formula(paste("hZ ~", v, "+ pgs_ht")), x[cnv == 0]), v)
  i0 <- co(r0, paste0(v, ":cnv")); i1 <- co(r1, paste0(v, ":cnv"))
  w10 <- wr(wc, f10c); wb <- wr(wc, fbc)
  out[[length(out) + 1]] <- data.table(group = G, pgs = PG[[v]], n_carriers = nc, n_carriers_step10 = x[cnv == 1 & !is.na(step10), .N], n_noncarriers = x[cnv == 0, .N],
    fs_step10_carriers = f10c[1], fs_step10_carriers_p = f10c[3], fs_step10_noncarriers = f10n[1], fs_bmiZ_carriers = fbc[1], fs_bmiZ_carriers_p = fbc[3], fs_bmiZ_noncarriers = fbn[1],
    height_per_pgs_carriers = wc[1], se_carriers = wc[2], p_carriers = wc[3], height_per_pgs_noncarriers = wn[1], p_noncarriers = wn[3],
    interaction = i0[1], interaction_p = i0[3], interaction_htPGSadj = i1[1], interaction_htPGSadj_p = i1[3],
    wald_height_per_step10 = w10[1], wald_step10_se = w10[2], wald_step10_p = w10[3], wald_height_per_sd_bmi = wb[1], wald_bmi_se = wb[2], wald_bmi_p = wb[3])
}
R <- rbindlist(out); fwrite(R, file.path(WD, "ukbb_mr_childhood_obesity.csv"))
print(R[, .(group, pgs, n_carriers, fs10_car = round(fs_step10_carriers, 3), fs10_non = round(fs_step10_noncarriers, 3), fsB_car = round(fs_bmiZ_carriers, 3), fsB_non = round(fs_bmiZ_noncarriers, 3),
            h_car = round(height_per_pgs_carriers, 3), p_car = signif(p_carriers, 2), h_non = round(height_per_pgs_noncarriers, 3), int = round(interaction, 3), int_p = signif(interaction_p, 2), int_adj_p = signif(interaction_htPGSadj_p, 2),
            wald10 = round(wald_height_per_step10, 2), wald10_p = signif(wald_step10_p, 2))], nrows = 50)
cat("ukbb MR table written\n")
