# 87_v2_sex_cnv_bmi_autism.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C, deliverable 4; EXPLORATORY)
# What it does: one coherent set of sex tests for SRS, SCQ and autism diagnosis.
#   (i)   BMI slope in female and male carriers and in female and male SSC unaffected relatives (gam y ~ s(age) + BMI + cohort x CNV)
#   (ii)  CNV x BMI within each sex (model of 61 / Fig. 5: bam y ~ carrier + s(age, by = carrier) + s(age, by = BMI) + carrier:BMI
#         + cohort x CNV; carriers of one sex against relatives of the same sex)
#   (iii) three-way CNV x BMI x sex from ONE model with all lower-order terms (bam y ~ carrier x sex + s(age, by = carrier x sex)
#         + s(age, by = BMI) + sex:BMI + carrier:BMI + carrier:sex:BMI + cohort x CNV); sensitivity: linear model with robust SE
#         (y ~ carrier * male * BMI + carrier * male * poly(age, 2) + cohort x CNV); per CNV and pooled, Cochran Q across the four CNVs;
#         within carriers only: BMI x sex (one model) and the 66 version (slopes fitted separately by sex); leave one cohort out
#   (iv)  autism diagnosis within carriers (growth table incl. CHOP): logistic diagnosis ~ BMI x sex + ns(age, 3) + cohort + CNV;
#         per sex, per CNV, pooled, per cohort, leave one cohort out; proportions by BMI band and sex
#   (v)   autism polygenic score (genotyped carriers G2MH + Simons VIP + SSC): the sex terms with the score, BMI x score and their
#         carrier interactions added; BMI x score within female and within male carriers and against same-sex relatives (two strata;
#         no four-way model), and their difference
# Inputs : v2_base() objects (C = combined SRS / SCQ sample, XG = genotyped subset, XA = diagnosis table), 86 cache (P1).
# Outputs: autism/sex_cnv_bmi_autism.csv, autism/block_sex_autism.{png,pdf} (tags e, f, g), autism/block_sex_autism_plotted.csv,
#          autism/sample_n_by_cnv_cohort_sex.csv (n per outcome x CNV x cohort x sex for every analysis of the package).
# null model: none; Wald tests of the interaction coefficients (bam / gam / HC1-robust lm / logistic), Cochran Q for heterogeneity.
base::source("C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/83_v2_autism_common.R")
suppressPackageStartupMessages(library(splines))
v2_base(); load(file.path("C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data", "v2_autism_cache", "pgs_stats_2026-11-20.RData"))
POOL <- "All four CNVs"; GR <- c(POOL, LOC4); options(width = 230)
out <- list(); add <- function(outcome, test, group, sex, r, n_ref = NA, note = "") out[[length(out) + 1]] <<- data.table(outcome = outcome, test = test, group = group, sex = sex, n = r[["n"]], n_ref = n_ref, estimate = r[["b"]], se = r[["se"]], p = r[["p"]], note = note)
NA4 <- function(n) c(n = n, b = NA_real_, se = NA_real_, p = NA_real_)
sslope <- function(d) tryCatch({ if (nrow(d) < 15) stop(); m <- gam(as.formula(paste("y ~ s(age, k = 8) + BMI", if (uniqueN(d$cg) > 1) "+ cg" else "")), data = d, method = "REML"); co <- summary(m)$p.table["BMI", ]; c(n = nrow(d), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) NA4(nrow(d)))
three_sex <- function(car, ref, extra = "") tryCatch({ if (min(table(car$sex)) < 15) stop()
  z <- rbind(car, ref, fill = TRUE)[, sg := factor(paste(is_cnv, sex))][, `:=`(cnvBMI = is_cnv * BMI, maleBMI = male * BMI, cnvmaleBMI = is_cnv * male * BMI)]
  m <- bam(as.formula(paste("y ~ sg + s(age, k = 10, by = sg) + s(age, k = 8, by = BMI) + maleBMI + cnvBMI + cnvmaleBMI", if (uniqueN(car$cg) > 1) "+ cg" else "", extra)), data = z, method = "fREML"); co <- summary(m)$p.table["cnvmaleBMI", ]
  c(n = nrow(car), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) NA4(nrow(car)))
three_sex_lm <- function(car, ref) tryCatch({ if (min(table(car$sex)) < 15) stop(); z <- rbind(car, ref, fill = TRUE)
  m <- lm(as.formula(paste("y ~ is_cnv * male * BMI + is_cnv * male * poly(age, 2)", if (uniqueN(car$cg) > 1) "+ cg" else "")), z); c(n = nrow(car), robHC(m, term_of(m, c("is_cnv", "male", "BMI")))) }, error = function(e) NA4(nrow(car)))
within_sex <- function(car, extra = "") tryCatch({ if (min(table(car$sex)) < 15) stop(); d <- copy(car)[, sexf := factor(sex)][, maleBMI := male * BMI]
  m <- gam(as.formula(paste("y ~ sexf + s(age, k = 8, by = sexf) + BMI + maleBMI", if (uniqueN(d$cg) > 1) "+ cg" else "", extra)), data = d, method = "REML"); co <- summary(m)$p.table["maleBMI", ]; c(n = nrow(d), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) NA4(nrow(car)))
dif <- function(f, m, n) { b <- m[["b"]] - f[["b"]]; se <- sqrt(f[["se"]]^2 + m[["se"]]^2); c(n = n, b = b, se = se, p = 2 * pnorm(-abs(b / se))) }
SX <- c("female", "male")
# ================================================================ (i)-(iii) SRS, SCQ
for (ph in c("SRS", "SCQ")) { C2 <- C[pheno == ph]; ref <- RELf(C2); carA <- C2[is_cnv == 1]
  for (s in SX) add(ph, "(i) BMI slope", "SSC unaffected relatives", s, sslope(ref[sex == s]))
  rs <- lapply(SX, function(s) sslope(ref[sex == s])); add(ph, "(i) BMI slope, male - female (separately fitted)", "SSC unaffected relatives", "male - female", dif(rs[[1]], rs[[2]], nrow(ref)))
  for (g in GR) { car <- if (g == POOL) carA else carA[locus == g]; if (nrow(car) < 30) next
    sl <- lapply(SX, function(s) sslope(car[sex == s])); for (i in 1:2) add(ph, "(i) BMI slope", g, SX[i], sl[[i]])
    add(ph, "within carriers: BMI x sex, slopes fitted separately by sex (66 method)", g, "male - female", dif(sl[[1]], sl[[2]], nrow(car)))
    add(ph, "within carriers: BMI x sex, one model", g, "male - female", within_sex(car))
    for (s in SX) { xc <- car[sex == s]; add(ph, "(ii) CNV x BMI within sex (vs same-sex SSC relatives)", g, s, if (nrow(xc) >= 15) vs_ref(xc, ref[sex == s]) else NA4(nrow(xc)), n_ref = nrow(ref[sex == s])) }
    add(ph, "(iii) CNV x BMI x sex, one model (bam)", g, "male - female", three_sex(car, ref), n_ref = nrow(ref))
    if (uniqueN(car$cg) > 1) add(ph, "(iii) CNV x BMI x sex, one model (bam) with sex-specific cohort x CNV intercepts (sensitivity)", g, "male - female", three_sex(copy(car)[, cgm := factor(ifelse(male == 1, cg, "f"))], copy(ref)[, cgm := factor("f")], "+ cgm"), n_ref = nrow(ref))
    add(ph, "(iii) CNV x BMI x sex, linear model with robust SE (sensitivity)", g, "male - female", three_sex_lm(car, ref), n_ref = nrow(ref)) }
  for (k in sort(unique(carA$cohort))) { add(ph, paste("(iii) CNV x BMI x sex, one model (bam); leave out carriers of", k), POOL, "male - female", three_sex(carA[cohort != k], ref), n_ref = nrow(ref))
    add(ph, paste("within carriers: BMI x sex, one model; leave out", k), POOL, "male - female", within_sex(carA[cohort != k])); add(ph, paste("(iii) CNV x BMI x sex, one model (bam); carriers of", k, "only"), POOL, "male - female", three_sex(carA[cohort == k], ref), n_ref = nrow(ref)) }
  for (k in c(0L, 1L)) add(ph, paste("(iii) CNV x BMI x sex, one model (bam); carriers", if (k == 1L) "with" else "without", "an autism diagnosis"), POOL, "male - female", three_sex(carA[asd2 %in% k], ref), n_ref = nrow(ref))
  for (a in list(c(4, 18), c(18, 99))) { rr <- if (a[1] >= 18) PARf(C2) else SIBf(C2); add(ph, sprintf("(iii) CNV x BMI x sex, one model (bam); age %s (vs SSC %s)", if (a[1] >= 18) "18+" else "< 18", if (a[1] >= 18) "parents" else "siblings"), POOL, "male - female", if (nrow(rr) > 50) three_sex(carA[age >= a[1] & age < a[2]], rr[age >= a[1] & age < a[2]]) else NA4(0), n_ref = nrow(rr)) }
  cat(ph, "done\n") }
R <- rbindlist(out); out <- list()
for (ph in c("SRS", "SCQ")) for (tt in c("(iii) CNV x BMI x sex, one model (bam)", "(ii) CNV x BMI within sex (vs same-sex SSC relatives)", "within carriers: BMI x sex, one model")) for (s in unique(R[outcome == ph & test == tt, sex])) {
  r <- R[outcome == ph & test == tt & sex == s & group %in% LOC4]; q <- cochranQ(r$estimate, r$se); out[[length(out) + 1]] <- data.table(outcome = ph, test = paste(tt, "- heterogeneity across the four CNVs (Cochran Q)"), group = POOL, sex = s, p = q[["p"]]) }
for (ph in c("SRS", "SCQ")) { r <- R[outcome == ph & test == "(iii) CNV x BMI x sex, one model (bam)" & group %in% LOC4 & !is.na(estimate)]; w <- 1 / r$se^2; b <- sum(w * r$estimate) / sum(w); se <- sqrt(1 / sum(w))   # fixed-effect combination of the four per-CNV estimates
  out[[length(out) + 1]] <- data.table(outcome = ph, test = "(iii) CNV x BMI x sex: inverse-variance mean of the four per-CNV estimates", group = POOL, sex = "male - female", n = sum(r$n), estimate = b, se = se, p = 2 * pnorm(-abs(b / se))) }
R <- rbind(R, rbindlist(out, fill = TRUE), fill = TRUE); out <- list()
# ================================================================ (iv) autism diagnosis within carriers
D <- XA[pheno == "BMI" & locus %in% LOC4][, male := as.numeric(sex == "male")]
dl <- function(d, f, v) tryCatch({ d <- copy(d); keep <- d[, .(v = sum(asd) > 0 & sum(asd) < .N), by = cohort][v == TRUE, cohort]; d <- d[cohort %in% keep]; if (sum(d$asd) < 10 || sum(1 - d$asd) < 10) stop()
  tm <- c(f, "ns(age, 3)", if (uniqueN(d$cohort) > 1) "cohort", if (uniqueN(d$locus) > 1) "locus"); m <- glm(as.formula(paste("asd ~", paste(tm, collapse = " + "))), data = d, family = binomial); co <- coef(summary(m))[v, ]; c(n = nrow(d), b = co[[1]], se = co[[2]], p = co[[4]]) }, error = function(e) NA4(nrow(d)))
DXN <- "Autism diagnosis"; LN <- "log odds ratio per SD of BMI"
for (g in GR) { d <- if (g == POOL) D else D[locus == g]
  for (s in SX) add(DXN, "(iv) diagnosis ~ BMI within sex (logistic)", g, s, dl(d[sex == s], "BMI", "BMI"), note = LN)
  add(DXN, "(iv) diagnosis ~ BMI x sex (logistic, one model)", g, "male - female", dl(d, "BMI * male", "BMI:male"), note = LN) }
for (k in sort(unique(D$cohort))) { add(DXN, paste("(iv) diagnosis ~ BMI x sex; cohort", k, "only"), POOL, "male - female", dl(D[cohort == k], "BMI * male", "BMI:male"), note = LN); add(DXN, paste("(iv) diagnosis ~ BMI x sex; leave out", k), POOL, "male - female", dl(D[cohort != k], "BMI * male", "BMI:male"), note = LN)
  for (s in SX) add(DXN, paste("(iv) diagnosis ~ BMI within sex; cohort", k, "only"), POOL, s, dl(D[cohort == k & sex == s], "BMI", "BMI"), note = LN) }
add(DXN, "(iv) diagnosis ~ BMI x sex; age < 18 y", POOL, "male - female", dl(D[age < 18], "BMI * male", "BMI:male"), note = LN)
R4 <- rbindlist(out); out <- list()
for (tt in c("(iv) diagnosis ~ BMI x sex (logistic, one model)", "(iv) diagnosis ~ BMI within sex (logistic)")) for (s in unique(R4[test == tt, sex])) { r <- R4[test == tt & sex == s & group %in% LOC4]; q <- cochranQ(r$estimate, r$se)
  out[[length(out) + 1]] <- data.table(outcome = DXN, test = paste(tt, "- heterogeneity across the four CNVs (Cochran Q)"), group = POOL, sex = s, p = q[["p"]]) }
R <- rbind(R, R4, rbindlist(out, fill = TRUE), fill = TRUE); out <- list()
wil <- function(k, n) { z <- 1.96; p <- k / n; den <- 1 + z^2 / n; ce <- (p + z^2 / (2 * n)) / den; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den; list(lo = ce - h, hi = ce + h) }
BANDS <- c("\u2264 30th", "30th-70th", "\u2265 70th"); D[, band := factor(fifelse(BMI >= qnorm(.7), BANDS[3], fifelse(BMI <= qnorm(.3), BANDS[1], BANDS[2])), BANDS)]
PS <- rbind(D[, .(n = .N, n_asd = sum(asd)), by = .(sex, band)][, group := POOL], D[, .(n = .N, n_asd = sum(asd)), by = .(group = locus, sex, band)], use.names = TRUE)[, prop := n_asd / n][, c("ci_lo", "ci_hi") := wil(n_asd, n)][order(group, sex, band)]
# ================================================================ (v) autism polygenic score
PT <- "+ pgs + I(pgs * BMI) + I(pgs * is_cnv) + I(pgs * is_cnv * BMI)"
wsex_lm <- function(car, pg = FALSE) { m <- lm(as.formula(paste("y ~", if (uniqueN(car$cg) > 1) "cg +" else "", "male * poly(age, 2) + male * BMI", if (pg) "+ BMI * pgs" else "")), car); c(n = nrow(car), robHC(m, term_of(m, c("male", "BMI")))) }
w3 <- function(d) c(n = nrow(d), rob(lm(as.formula(paste("y ~", if (uniqueN(d$cg) > 1) "cg +" else "", "poly(age, 2) + BMI * pgs")), d), "BMI:pgs"))
for (ph in c("SRS", "SCQ")) { X3 <- XG[pheno == ph]; car <- X3[is_cnv == 1]; ref <- RELf(X3)
  add(ph, "(v) CNV x BMI x sex in the genotyped sample, without the autism score", POOL, "male - female", three_sex(car, ref), n_ref = nrow(ref))
  add(ph, "(v) CNV x BMI x sex in the genotyped sample, adjusted for score, BMI x score and their carrier interactions", POOL, "male - female", three_sex(car, ref, PT), n_ref = nrow(ref))
  add(ph, "(v) within carriers BMI x sex (linear, robust SE), genotyped sample, without the score", POOL, "male - female", wsex_lm(car)); add(ph, "(v) within carriers BMI x sex (linear, robust SE), genotyped sample, adjusted for BMI x score", POOL, "male - female", wsex_lm(car, TRUE))
  wf <- w3(car[sex == "female"]); wm <- w3(car[sex == "male"]); add(ph, "(v) BMI x autism score within carriers of one sex", POOL, "female", wf); add(ph, "(v) BMI x autism score within carriers of one sex", POOL, "male", wm)
  add(ph, "(v) BMI x autism score within carriers, male - female (separately fitted)", POOL, "male - female", dif(wf, wm, nrow(car)))
  tf <- three(car[sex == "female"], ref[sex == "female"]); tm <- three(car[sex == "male"], ref[sex == "male"]); add(ph, "(v) CNV x BMI x autism score within one sex (vs same-sex relatives)", POOL, "female", tf, n_ref = nrow(ref[sex == "female"])); add(ph, "(v) CNV x BMI x autism score within one sex (vs same-sex relatives)", POOL, "male", tm, n_ref = nrow(ref[sex == "male"]))
  add(ph, "(v) CNV x BMI x autism score, male - female (separately fitted)", POOL, "male - female", dif(tf, tm, nrow(car))) }
P1[, sex := fifelse(male == 1, "male", "female")]
add(DXN, "(v) diagnosis ~ BMI x sex, genotyped carriers (SRS / SCQ sample), without the score", POOL, "male - female", dl(P1, "BMI * male", "BMI:male"), note = LN); add(DXN, "(v) diagnosis ~ BMI x sex, genotyped carriers, adjusted for BMI x score", POOL, "male - female", dl(P1, "BMI * male + BMI * pgs", "BMI:male"), note = LN)
for (s in SX) add(DXN, "(v) diagnosis ~ BMI x autism score within one sex (genotyped carriers)", POOL, s, dl(P1[sex == s], "BMI * pgs", "BMI:pgs"), note = "log odds ratio per SD x SD")
R <- rbind(R, rbindlist(out, fill = TRUE), fill = TRUE)
R[note == "" | is.na(note), note := "SD of trait (Box-Cox) per SD of BMI"][, status := "exploratory; no correction for multiple testing"]
fwrite(R, file.path(V2OUT, "sex_cnv_bmi_autism.csv"))
print(R[, .(outcome = substr(outcome, 1, 6), test = substr(test, 1, 95), group = substr(group, 1, 14), sex = substr(sex, 1, 6), n, est = round(estimate, 3), se = round(se, 3), p = signif(p, 3))], nrows = 400)
print(PS[group == POOL])
# ---------------------------------------------------------------- sample counts for the README
cnt <- function(d, lab) d[, .N, by = .(locus, cohort, sex)][, outcome := lab]
SN <- rbind(cnt(C[pheno == "SRS" & is_cnv == 1], "SRS (CNV x BMI; sex)"), cnt(C[pheno == "SCQ" & is_cnv == 1], "SCQ (CNV x BMI; sex)"), cnt(XG[pheno == "SRS" & is_cnv == 1], "SRS with autism polygenic score"), cnt(XG[pheno == "SCQ" & is_cnv == 1], "SCQ with autism polygenic score"),
            cnt(D, "Autism diagnosis with BMI (growth table)"), cnt(XA[pheno == "Head circumference" & locus %in% LOC4], "Autism diagnosis with head circumference"), cnt(XA[pheno == "Height" & locus %in% LOC4], "Autism diagnosis with height"), cnt(P1, "Autism diagnosis with autism polygenic score (SRS / SCQ sample)"),
            C[pheno %in% c("SRS", "SCQ") & locus == "SSC relative", .N, by = .(outcome = paste(pheno, "reference"), cohort, sex)][, locus := "SSC unaffected relatives"], XG[pheno %in% c("SRS", "SCQ") & locus == "SSC relative", .N, by = .(outcome = paste(pheno, "reference with autism polygenic score"), cohort, sex)][, locus := "SSC unaffected relatives"], use.names = TRUE)
SNw <- dcast(SN, outcome + locus + cohort ~ sex, value.var = "N", fill = 0)[, total := female + male]; fwrite(SNw, file.path(V2OUT, "sample_n_by_cnv_cohort_sex.csv")); print(SNw, nrows = 200)
# ================================================================ figure block (tags e, f, g)
SC <- c(female = C_F, male = C_M); BS <- 15; TS <- 4.3
# e: SRS against BMI by sex, carriers and relatives (partial residuals: cohort x CNV and sex-specific quadratic age removed)
C2 <- C[pheno == "SRS"]; pr <- function(d, grp) { d <- copy(d); m <- lm(as.formula(paste("y ~", if (uniqueN(d$cg) > 1) "cg +" else "", "male * poly(age, 2) + male * BMI")), d); Xm <- model.matrix(m); nu <- grep("poly\\(age|^cg", colnames(Xm)); cf <- coef(m); cf[is.na(cf)] <- 0
  ct <- as.vector(Xm[, nu, drop = FALSE] %*% cf[nu]); d[, ya := y - ct + mean(ct)][, grp := grp]; d }
Z <- rbind(pr(C2[is_cnv == 1], "CNV carriers"), pr(RELf(C2), "SSC unaffected relatives"))[, grp := factor(grp, c("CNV carriers", "SSC unaffected relatives"))]
fl <- Z[, { f <- lm(ya ~ BMI); g <- seq(quantile(BMI, 0.02), quantile(BMI, 0.98), length.out = 40); p <- predict(f, data.frame(BMI = g), se.fit = TRUE); .(BMI = g, fit = p$fit, lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit) }, by = .(grp, sex)]
s1 <- R[outcome == "SRS" & test == "(i) BMI slope" & group %in% c(POOL, "SSC unaffected relatives")][, grp := factor(fifelse(group == POOL, "CNV carriers", "SSC unaffected relatives"), levels(Z$grp))]
s1[, lab := sprintf("%s: slope %+.2f (%s, n = %d)", sex, estimate, plab(p), n)][, yy := 6.3 - 0.55 * (sex == "male")]
i3 <- R[outcome == "SRS" & test == "(iii) CNV x BMI x sex, one model (bam)" & group == POOL]; i3s <- R[outcome == "SRS" & test == "(iii) CNV x BMI x sex, one model (bam) with sex-specific cohort x CNV intercepts (sensitivity)" & group == POOL]; w1 <- R[outcome == "SRS" & test == "within carriers: BMI x sex, one model" & group == POOL]
pe <- ggplot(Z, aes(BMI, ya, colour = sex)) + geom_point(data = Z[is_cnv == 0], size = 0.4, alpha = 0.06) + geom_point(data = Z[is_cnv == 1], size = 1, alpha = 0.4) +
  geom_ribbon(data = fl, aes(BMI, ymin = lo, ymax = hi, fill = sex), inherit.aes = FALSE, alpha = 0.22) + geom_line(data = fl, aes(BMI, fit), linewidth = 1.3) +
  geom_label(data = s1, aes(x = -3, y = yy, label = lab, colour = sex), hjust = 0, vjust = 1, size = TS - 0.5, fontface = "bold", fill = "white", alpha = 0.8, linewidth = 0, show.legend = FALSE) +
  facet_wrap(~grp, nrow = 1) + scale_colour_manual(values = SC, name = NULL) + scale_fill_manual(values = SC, guide = "none") + coord_cartesian(xlim = c(-3, 4.2), ylim = c(-2.5, 6.3)) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(x = "BMI (SD)", y = "Social Responsiveness Scale (SD)", tag = "e", title = "BMI and autistic traits by sex",
       subtitle = sprintf("CNV \u00d7 BMI \u00d7 sex %s;\n%s with sex-specific cohort \u00d7 CNV terms", plab(i3$p), plab(i3s$p))) + th_v2(BS) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.position = "bottom", legend.text = element_text(size = 14))
# f: CNV x BMI by sex, SRS and SCQ, per CNV and pooled
Ff <- R[test == "(ii) CNV x BMI within sex (vs same-sex SSC relatives)" & !is.na(estimate)][, group := factor(group, rev(GR))][, outcome := factor(outcome, c("SRS", "SCQ"))]
Fd <- R[test == "(iii) CNV x BMI x sex, one model (bam)" & !is.na(p)][, group := factor(group, rev(GR))][, outcome := factor(outcome, c("SRS", "SCQ"))][, lab := paste0("sex difference", "
", plab(p))]
hq3 <- R[grepl("^\\(iii\\) CNV x BMI x sex, one model \\(bam\\) - heterogeneity", test)]
xr <- c(min(Ff$estimate - 1.96 * Ff$se), max(Ff$estimate + 1.96 * Ff$se)); pdg <- position_dodge(width = 0.55)
pf <- ggplot(Ff, aes(estimate, group, colour = sex, group = sex)) + geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + geom_errorbar(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), width = 0, orientation = "y", position = pdg, linewidth = 0.9) +
  geom_point(aes(shape = p < 0.05), size = 3.3, position = pdg) + geom_text(data = Fd, aes(x = xr[2] + 0.04 * diff(xr), y = group, label = lab, fontface = face(p)), inherit.aes = FALSE, hjust = 0, size = TS - 0.6, lineheight = 0.9, colour = "grey15") +
  facet_wrap(~outcome, nrow = 1) + scale_colour_manual(values = SC, name = NULL) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") + coord_cartesian(xlim = c(xr[1], xr[2] + 0.62 * diff(xr))) + scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  labs(x = "CNV \u00d7 BMI (carriers minus same-sex relatives; SD per SD of BMI, 95% CI)", y = NULL, tag = "f", title = "CNV \u00d7 BMI by sex",
       subtitle = sprintf("Filled: p < 0.05. Sex difference from one CNV \u00d7 BMI \u00d7 sex model.\nHeterogeneity of the sex difference across CNVs: SRS %s, SCQ %s", plab(hq3[outcome == "SRS", p]), plab(hq3[outcome == "SCQ", p]))) +
  th_v2(BS) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.position = "bottom", legend.text = element_text(size = 14))
# g: diagnosis by BMI band and sex
g1 <- PS[group == POOL]; r4 <- R[outcome == DXN & group == POOL]; rf <- r4[test == "(iv) diagnosis ~ BMI within sex (logistic)"]; ri <- r4[test == "(iv) diagnosis ~ BMI x sex (logistic, one model)"]; ym <- max(g1$ci_hi) * 100
lg <- sprintf("OR per SD of BMI:\nfemale %.2f (%s)\nmale %.2f (%s)\nBMI \u00d7 sex %s", exp(rf[sex == "female", estimate]), plab(rf[sex == "female", p]), exp(rf[sex == "male", estimate]), plab(rf[sex == "male", p]), plab(ri$p))
pg <- ggplot(g1, aes(band, 100 * prop, fill = sex)) + geom_col(position = position_dodge(0.8), width = 0.72) + geom_errorbar(aes(ymin = 100 * ci_lo, ymax = 100 * ci_hi), position = position_dodge(0.8), width = 0.18) +
  geom_text(aes(y = -0.06 * ym, label = sprintf("%d/%d", n_asd, n)), position = position_dodge(0.8), size = TS - 1.2, colour = "grey20") + annotate("text", x = 0.5, y = ym * 1.45, label = lg, hjust = 0, vjust = 1, size = TS - 0.3, lineheight = 1, fontface = face(ri$p)) +
  scale_fill_manual(values = SC, name = NULL) + coord_cartesian(ylim = c(-0.07 * ym, ym * 1.45)) + labs(x = "BMI centile band", y = "With an autism diagnosis (%)", tag = "g", title = "Autism diagnosis by BMI and sex", subtitle = sprintf("All four CNVs (n = %d carriers);\ndiagnosed / all below the bars", nrow(D))) +
  th_v2(BS) + theme(legend.position = "bottom", legend.text = element_text(size = 14))
fig <- (pe | pf | pg) + plot_layout(widths = c(1.15, 1.1, 0.75))
save_v2(fig, "block_sex_autism", 20, 6.8)
fwrite(rbind(s1[, .(panel = "e", outcome = "SRS", group = as.character(grp), sex, n, value = estimate, se, p)], Ff[, .(panel = "f", outcome, group, sex, n, value = estimate, se, p)], Fd[, .(panel = "f (sex difference)", outcome, group, sex, n, value = estimate, se, p)],
             g1[, .(panel = "g", outcome = DXN, group = paste(group, band), sex, n, value = prop, se = NA_real_, p = NA_real_)]), file.path(V2OUT, "block_sex_autism_plotted.csv"))
cat("87 done\n")
