# CNV x BMI on SRS / SCQ in the combined carrier sample (G2MH + SVIP + Searchlight + SSC carriers; 60_combined_srs_scq_build.R)
# versus SSC unaffected relatives (JS / Claude 2026-11-14). Same model as Fig. 5k-t (24 vs_unaff): y ~ is_cnv + s(age, by = is_cnv)
# + s(age, by = BMI) + cnvBMI, plus carrier cohort (and CNV, when pooled) indicators; Box-Cox scale; ordinal check (scale-free;
# per form, inverse-variance combined). Reported for the combined sample, per cohort, per CNV, and the cohort heterogeneity of the
# carrier BMI slope; within-cohort check in SVIP / Searchlight against their own non-carrier relatives.
# Output: combined_cnv_x_bmi.csv (aggregate statistics only).
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ordinal); library(sandwich) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
C <- readRDS(file.path(DATA, "combined_srs_scq.rds")); LOC <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
C[, cg := fifelse(is_cnv == 1, paste(cohort, locus, sep = ":"), "ref")]
vs_ref <- function(xc, ref) tryCatch({ z <- rbind(xc, ref)[, is_c := factor(ifelse(is_cnv == 1, "CNV", "ref"), c("ref", "CNV"))][, cnvBMI := is_cnv * BMI]
  ex <- if (uniqueN(xc$cg) > 1) "+ cg" else ""
  m <- bam(as.formula(paste("y ~ is_c + s(age, k = 10, by = is_c) + s(age, k = 8, by = BMI) + cnvBMI", ex)), data = z, method = "fREML"); co <- summary(m)$p.table["cnvBMI", ]
  c(n = nrow(xc), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) c(n = nrow(xc), b = NA, se = NA, p = NA))
slope <- function(xc) tryCatch({ ex <- if (uniqueN(xc$cg) > 1) "+ cg" else ""; m <- gam(as.formula(paste("y ~ s(age, k = 8) + male + BMI", ex)), data = xc, method = "REML")
  co <- summary(m)$p.table["BMI", ]; c(n = nrow(xc), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) c(n = nrow(xc), b = NA, se = NA, p = NA))
ord <- function(xc, ref) { z <- rbind(xc, ref); r <- t(sapply(split(z, z$form2), function(d) tryCatch({ if (sum(d$is_cnv) < 15) stop()
    d <- copy(d); br <- unique(quantile(d$raw, seq(0, 1, length.out = 41), type = 1)); d[, ry := factor(findInterval(raw, br, rightmost.closed = TRUE), ordered = TRUE)][, `:=`(a1 = as.numeric(scale(age)), cnvBMI = is_cnv * BMI)]
    m <- clm(as.formula(paste("ry ~ is_cnv * (a1 + I(a1^2)) + BMI + cnvBMI", if (uniqueN(d$cg) > 2) "+ cg" else "")), data = d); co <- coef(summary(m)); c(co["cnvBMI", 1], co["cnvBMI", 2]) }, error = function(e) c(NA, NA))))
  r <- r[!is.na(r[, 1]), , drop = FALSE]; if (!nrow(r)) return(c(n = nrow(xc), b = NA, se = NA, p = NA)); w <- 1 / r[, 2]^2; b <- sum(w * r[, 1]) / sum(w); se <- sqrt(1 / sum(w)); c(n = nrow(xc), b = b, se = se, p = 2 * pnorm(-abs(b / se))) }
out <- list(); add <- function(ph, grp, coh, test, r) out[[length(out) + 1]] <<- data.table(trait = ph, cnv = grp, cohort = coh, test = test, n = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]])
for (ph in c("SRS", "SCQ")) {
  ref <- C[pheno == ph & locus == "SSC relative"]
  add(ph, "-", "SSC relatives", "BMI slope", slope(ref))
  for (g in c("All four CNVs", LOC)) for (coh in c("Combined", "G2MH", "SVIP", "Searchlight", "SSC")) {
    xc <- C[pheno == ph & is_cnv == 1 & (g == "All four CNVs" | locus == g) & (coh == "Combined" | cohort == coh)]
    if (nrow(xc) < 15) next
    add(ph, g, coh, "BMI slope in carriers", slope(xc)); add(ph, g, coh, "CNV x BMI vs SSC relatives (Box-Cox)", vs_ref(xc, ref))
    if (coh == "Combined") add(ph, g, coh, "CNV x BMI vs SSC relatives (ordinal)", ord(xc, ref))
  }
  # cohort heterogeneity of the carrier BMI slope (16p11.2 carriers, the CNVs shared by cohorts)
  x16 <- C[pheno == ph & is_cnv == 1 & locus %in% LOC[1:2]][, coh := factor(cohort)]
  m1 <- lm(y ~ cg + poly(age, 2) + male + BMI * coh, x16); m0 <- lm(y ~ cg + poly(age, 2) + male + BMI + coh, x16)
  add(ph, "16p11.2 del + dup", "Combined", "cohort heterogeneity of the carrier BMI slope (F test)", c(n = nrow(x16), b = NA, se = NA, p = anova(m0, m1)$`Pr(>F)`[2]))
  xa <- C[pheno == ph & is_cnv == 1][, coh := factor(cohort)]
  m1 <- lm(y ~ cg + poly(age, 2) + male + BMI * locus + BMI:coh, xa); m0 <- lm(y ~ cg + poly(age, 2) + male + BMI + locus, xa)
  add(ph, "All four CNVs", "Combined", "CNV and cohort heterogeneity of the carrier BMI slope (F test)", c(n = nrow(xa), b = NA, se = NA, p = anova(m0, m1)$`Pr(>F)`[2]))
  # within-cohort references (SVIP / Searchlight non-carrier relatives)
  own <- C[pheno == ph & locus == "Non-carrier relative" & cohort %in% c("SVIP", "Searchlight")]
  if (nrow(own) >= 15) { add(ph, "-", "SVIP/Searchlight relatives", "BMI slope", slope(own[, cg := "ref"]))
    xc <- C[pheno == ph & is_cnv == 1 & cohort %in% c("SVIP", "Searchlight")]; add(ph, "16p11.2 del + dup", "SVIP + Searchlight", "CNV x BMI vs own non-carrier relatives", vs_ref(xc, own[, is_cnv := 0])) }
  cat(ph, "done\n")
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "combined_cnv_x_bmi.csv"))
options(width = 220); print(R[, .(trait, cnv, cohort, test = substr(test, 1, 48), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
