# 87b_v2_family_clustered_check.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C: follow-up check)
# What it does: for each headline test of the package, refits the model exactly as reported and adds a family-clustered sandwich
#   SE and P (cluster = family), to check whether relatedness (SSC parents and siblings; carriers of one family) changes a conclusion.
#   Family: SSC family number (relatives, probands, SSC carriers); Simons VIP / Searchlight family number (SFARI id prefix); G2MH and
#   Cardiff: rarecnv_family_id of the G2MH release, else site + local_family_id, else the individual; CHOP and CHU Sainte-Justine
#   records carry no family identifier, so each person is their own cluster.
#   lm / glm: sandwich::vcovCL (HC1 type). bam / gam (penalised): A M A with A = (X'X + S)^-1 = Vp / scale and
#   M = sum over families of (X_g' e_g)(X_g' e_g)', times G / (G - 1); z test.
#   Tests: (1) CNV x BMI on SRS and SCQ pooled, SRS 22q11.2 del; (2) CNV x BMI x autism score on SRS (vs relatives; within carriers)
#   and SCQ; (3) CNV x BMI x sex on SRS, plain and with sex-specific cohort x CNV intercepts; (4) pooled diagnosis contrasts for head
#   circumference and BMI, and the logistic odds ratios of diagnosis on BMI Z and head-circumference Z.
# Inputs : v2_base() objects; G2MH release (family columns). Output: autism/family_clustered_check.csv.
# null model: none; Wald tests with model-based and family-clustered sandwich variance.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
suppressPackageStartupMessages(library(splines)); v2_base(); options(width = 230)
DBX <- "C:/Users/Jonathan Sebat/sebatlab Dropbox"
rl <- fread(file.path(DBX, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"), select = c("guid", "rarecnv_family_id", "site_id", "local_family_id"), colClasses = "character")
rl <- rl[guid != ""][, f := fifelse(!is.na(rarecnv_family_id) & rarecnv_family_id != "", paste0("G2MH:r", rarecnv_family_id), fifelse(!is.na(local_family_id) & local_family_id != "", paste0("G2MH:", site_id, ":", local_family_id), NA_character_))]
rl <- rl[!is.na(f)][, .SD[1], by = guid]; g2f <- setNames(rl$f, rl$guid)
cat("G2MH release: GUIDs with a family id", length(g2f), "; families", uniqueN(g2f), "\n")
pre <- function(x) sub("[.].*$", "", x)
fam_srs <- function(d) d[, fifelse(cohort == "SSC", paste0("SSC:", pre(id)), fifelse(cohort %in% c("SVIP", "Searchlight"), paste0("SF:", pre(id)), fifelse(!is.na(g2f[id]), g2f[id], paste0("ind:", id))))]
fam_growth <- function(d) { m <- idmap[match(paste(d$id, d$cohort), paste(idmap$id, idmap$cohort))]
  fifelse(d$cohort == "SSC", paste0("SSC:", pre(m$sampleID)), fifelse(d$cohort %in% c("SVIP", "Searchlight"), paste0("SF:", pre(m$sampleID)), fifelse(d$cohort %in% c("G2MH", "Cardiff") & !is.na(g2f[m$GUID]), g2f[m$GUID], paste0("ind:", d$id)))) }
clgam <- function(m, v, fam) { X <- predict(m, type = "lpmatrix"); e <- residuals(m, type = "response"); A <- m$Vp / m$sig2; S <- rowsum(X * e, fam); G <- nrow(S); V <- A %*% crossprod(S) %*% A * G / (G - 1); j <- match(v, names(coef(m))); sqrt(V[j, j]) }
share <- function(fam, car) { n <- table(fam); sum(n[fam[car]] > 1) }
out <- list(); add <- function(block, test, n, ncar, nshare, nfam, b, se, p, sec, scale = "SD per SD") out[[length(out) + 1]] <<- data.table(block = block, test = test, n = n, n_carriers = ncar, carriers_sharing_a_family = nshare, n_families = nfam, estimate = b, se_model = se, p_model = p, se_family_clustered = sec, p_family_clustered = 2 * pnorm(-abs(b / sec)), scale = scale)
C[, fam := fam_srs(C)]; XG[, fam := fam_srs(XG)]
# ---------------------------------------------------------------- (1) CNV x BMI (model of 61 / 51b)
f1 <- function(xc, ref, lab) { z <- rbind(xc, ref)[, is_c := factor(ifelse(is_cnv == 1, "CNV", "ref"), c("ref", "CNV"))][, cnvBMI := is_cnv * BMI]
  m <- bam(as.formula(paste("y ~ is_c + s(age, k = 10, by = is_c) + s(age, k = 8, by = BMI) + cnvBMI", if (uniqueN(xc$cg) > 1) "+ cg" else "")), data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]
  add("1 CNV x BMI", lab, nrow(z), nrow(xc), share(z$fam, z$is_cnv == 1), uniqueN(z$fam), co[[1]], co[[2]], co[[4]], clgam(m, "cnvBMI", z$fam)) }
f1(C[pheno == "SRS" & is_cnv == 1], RELf(C[pheno == "SRS"]), "SRS, all four CNVs vs SSC relatives"); f1(C[pheno == "SCQ" & is_cnv == 1], RELf(C[pheno == "SCQ"]), "SCQ, all four CNVs vs SSC relatives")
f1(C[pheno == "SRS" & is_cnv == 1 & locus == "22q11.2 del"], RELf(C[pheno == "SRS"]), "SRS, 22q11.2 del vs SSC relatives")
# ---------------------------------------------------------------- (2) CNV x BMI x autism score (62 / 64)
f2 <- function(z, f, v, lab) { m <- lm(f, z); r <- robHC(m, v); sc <- sqrt(vcovCL(m, cluster = z$fam, type = "HC1")[v, v]); add("2 CNV x BMI x autism score", lab, nrow(z), sum(z$is_cnv == 1), share(z$fam, z$is_cnv == 1), uniqueN(z$fam), r[["b"]], r[["se"]], r[["p"]], sc) }
for (ph in c("SRS", "SCQ")) { car <- XG[pheno == ph & is_cnv == 1]; z <- rbind(car, RELf(XG[pheno == ph]))
  f2(z, y ~ is_cnv * (poly(age, 2) + male) + cg + is_cnv * BMI * pgs, "is_cnv:BMI:pgs", paste(ph, "three-way vs SSC relatives (model SE = HC1 robust)"))
  if (ph == "SRS") { f2(car, y ~ cg + poly(age, 2) + male + BMI * pgs, "BMI:pgs", "SRS within carriers, BMI x score (model SE = HC1 robust)"); f2(car, y ~ cg + poly(age, 2) + male + BMI * pgs + cohort:BMI + cohort:pgs, "BMI:pgs", "SRS within carriers, BMI x score, with cohort x BMI and cohort x score (model SE = HC1 robust)") } }
# ---------------------------------------------------------------- (3) CNV x BMI x sex (87)
f3 <- function(car, ref, extra, lab) { z <- rbind(car, ref, fill = TRUE)[, sg := factor(paste(is_cnv, sex))][, `:=`(cnvBMI = is_cnv * BMI, maleBMI = male * BMI, cnvmaleBMI = is_cnv * male * BMI)]
  m <- bam(as.formula(paste("y ~ sg + s(age, k = 10, by = sg) + s(age, k = 8, by = BMI) + maleBMI + cnvBMI + cnvmaleBMI + cg", extra)), data = z, method = "fREML"); co <- summary(m)$p.table["cnvmaleBMI", ]
  add("3 CNV x BMI x sex", lab, nrow(z), nrow(car), share(z$fam, z$is_cnv == 1), uniqueN(z$fam), co[[1]], co[[2]], co[[4]], clgam(m, "cnvmaleBMI", z$fam)) }
car <- C[pheno == "SRS" & is_cnv == 1]; ref <- RELf(C[pheno == "SRS"])
f3(car, ref, "", "SRS, one model"); f3(copy(car)[, cgm := factor(ifelse(male == 1, cg, "f"))], copy(ref)[, cgm := factor("f")], "+ cgm", "SRS, one model with sex-specific cohort x CNV intercepts")
# ---------------------------------------------------------------- (4) autism diagnosis (85)
D <- XA[pheno %in% c("Height", "Head circumference", "BMI") & locus %in% LOC4]; D[, fam := fam_growth(D)]
for (ph in c("Head circumference", "BMI")) { x <- copy(D[pheno == ph])[, `:=`(sex = factor(sex), cohort = factor(cohort), locus = factor(locus))]
  m <- gam(y ~ asd + s(age, k = 8) + sex + cohort + locus, data = x, method = "REML"); co <- summary(m)$p.table["asd", ]
  add("4 autism diagnosis", paste(ph, "Z, diagnosis vs none, all four CNVs"), nrow(x), nrow(x), share(x$fam, rep(TRUE, nrow(x))), uniqueN(x$fam), co[[1]], co[[2]], co[[4]], clgam(m, "asd", x$fam), "SD")
  d <- copy(D[pheno == ph])[, z := y]; keep <- d[, .(v = sum(asd) > 0 & sum(asd) < .N), by = cohort][v == TRUE, cohort]; d <- d[cohort %in% keep]
  g <- glm(asd ~ z + ns(age, 3) + sex + cohort + locus, data = d, family = binomial); cg <- coef(summary(g))["z", ]
  add("4 autism diagnosis", paste("logistic: diagnosis on", ph, "Z, all four CNVs"), nrow(d), nrow(d), share(d$fam, rep(TRUE, nrow(d))), uniqueN(d$fam), cg[[1]], cg[[2]], cg[[4]], sqrt(vcovCL(g, cluster = d$fam, type = "HC1")["z", "z"]), "log odds ratio per SD") }
R <- rbindlist(out)
R[, se_ratio := se_family_clustered / se_model][, crosses_0.05 := (p_model < 0.05) != (p_family_clustered < 0.05)][, conclusion_changes := fifelse(crosses_0.05, "yes: P crosses 0.05", "no")]
R[, note := "carriers_sharing_a_family = carriers whose family has another person (carrier or relative) in the same analysis; CHOP and CHUSJ have no family identifier (own cluster)"]
fwrite(R, file.path(V2OUT, "family_clustered_check.csv"))
print(R[, .(test = substr(test, 1, 70), n, ncar = n_carriers, share = carriers_sharing_a_family, nfam = n_families, est = round(estimate, 4), se = round(se_model, 4), p = signif(p_model, 3), se_cl = round(se_family_clustered, 4), p_cl = signif(p_family_clustered, 3), ratio = round(se_ratio, 2), changes = conclusion_changes)])
cat("carriers sharing a family, by cohort (SRS sample):\n"); z <- rbind(C[pheno == "SRS" & is_cnv == 1], RELf(C[pheno == "SRS"])); nn <- table(z$fam); print(z[is_cnv == 1, .(carriers = .N, sharing = sum(nn[fam] > 1), own_cluster_no_family_id = sum(grepl("^ind:", fam))), by = cohort])
cat("growth table (BMI rows):\n"); x <- D[pheno == "BMI"]; nn <- table(x$fam); print(x[, .(carriers = .N, sharing = sum(nn[fam] > 1), own_cluster_no_family_id = sum(grepl("^ind:", fam))), by = cohort])
cat("87b done\n")
