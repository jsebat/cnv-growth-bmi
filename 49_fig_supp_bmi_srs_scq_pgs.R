# Supplementary Fig. 7 (JS / Claude 2026-11-14): is the BMI -> SRS/SCQ association in CNV carriers explained by genetics shared
# between body mass, growth and neurodevelopment? G2MH carriers of the four CNVs pooled (the only carriers with polygenic scores);
# every model has CNV, CNV-specific quadratic age and sex (as in 47/48); scores residualised on 10 PCs, z within G2MH; HC1 robust SE.
#  a  BMI slope with no adjustment and after adding growth, adiposity, neurodevelopmental, all ten scores; measured height (subset)
#  b  each score's association with the trait (reduced form, 95% CI) next to the value expected if the score acted only through
#     BMI (its first-stage association with BMI x the unadjusted BMI slope)
#  c  observational (OLS) vs instrumented (2SLS; childhood obesity, childhood BMI and adult BMI scores) BMI effect, both adjusted for
#     the growth and neurodevelopmental scores; per-instrument Wald ratios; Sargan over-identification P
# Output: fig_supp_bmi_srs_scq_pgs.{png,pdf} + fig_supp_bmi_srs_scq_pgs_table.csv (aggregate statistics only).
# JS 2026-11-14: SRS/SCQ on the Box-Cox scale by default (SRSSCQ_SCALE=linear for the original linear z); outputs carry the scale suffix.
if (Sys.getenv("SRSSCQ_SCALE") == "") Sys.setenv(SRSSCQ_SCALE = "boxcox"); SUF <- if (Sys.getenv("SRSSCQ_SCALE") == "boxcox") "_boxcox" else ""
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s48 <- readLines(file.path(HERE, "48_bmi_srs_scq_asd_pgs.R")); eval(parse(text = s48[1:grep("^ND <- ", s48)]))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
LAB <- c(pgs_cobes = "Childhood obesity", pgs_cbmi = "Childhood BMI", pgs_bmi = "Adult BMI", pgs_ht = "Height", pgs_hci = "Infant head circumference",
         pgs_aam = "Age at menarche", pgs_asd = "Autism", pgs_adhd = "ADHD", pgs_edu = "Educational attainment", pgs_iq = "Cognitive ability (IQ)")
CLS <- c(ADI = "Adiposity", GRO = "Growth", ND = "Neurodevelopmental"); cls <- c(setNames(rep("Adiposity", 3), ADI), setNames(rep("Growth", 3), GRO), setNames(rep("Neurodevelopmental", 4), ND))
TEN <- c(ADI, GRO, ND); OC <- c(SRS = "#C0392B", SCQ = "#2E86C1")
tsls <- function(d, cov, ins) { Xc <- model.matrix(fm("y", cov), d); Z <- cbind(Xc, as.matrix(d[, ..ins])); X <- cbind(Xc, BMI = d$BMI); y <- d$y
  Xh <- Z %*% solve(crossprod(Z), crossprod(Z, X)); A <- solve(crossprod(Xh, X)); b <- A %*% crossprod(Xh, y); e <- as.vector(y - X %*% b)
  V <- A %*% crossprod(Xh * e) %*% t(A) * nrow(X) / (nrow(X) - ncol(X)); k <- ncol(X); s <- nrow(d) * summary(lm(e ~ Z - 1))$r.squared
  c(b = b[k], se = sqrt(V[k, k]), sargan_p = pchisq(s, length(ins) - 1, lower.tail = FALSE)) }
A <- B <- C <- list()
for (ph in c("SRS", "SCQ")) {
  d <- W[pheno == ph]
  sets <- list(`None` = NULL, `+ growth scores (3)` = GRO, `+ adiposity scores (3)` = ADI, `+ neurodevelopmental scores (4)` = ND, `+ all ten scores` = TEN)
  for (s in names(sets)) { r <- rob(lm(fm("y", c("BMI", sets[[s]])), d), "BMI"); A[[length(A) + 1]] <- data.table(outcome = ph, set = s, n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  dh <- d[!is.na(ht_z)]
  for (s in c("None (height subset)", "+ measured height z")) { r <- rob(lm(fm("y", c("BMI", if (s != "None (height subset)") "ht_z")), dh), "BMI")
    A[[length(A) + 1]] <- data.table(outcome = ph, set = s, n = nrow(dh), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  obs <- rob(lm(fm("y", "BMI"), d), "BMI")[["b"]]
  for (v in TEN) { f1 <- rob(lm(fm("BMI", v), d), v); rf <- rob(lm(fm("y", v), d), v)
    B[[length(B) + 1]] <- data.table(outcome = ph, score = v, class = cls[[v]], n = nrow(d), fs_b = f1[["b"]], fs_F = (f1[["b"]] / f1[["se"]])^2, b = rf[["b"]], se = rf[["se"]], p = rf[["p"]], expected = f1[["b"]] * obs) }
  COV <- c(GRO, ND); o <- rob(lm(fm("y", c("BMI", COV)), d), "BMI"); iv <- tsls(d, COV, ADI)
  fs <- anova(lm(fm("BMI", COV), d), lm(fm("BMI", c(COV, ADI)), d))
  C[[length(C) + 1]] <- data.table(outcome = ph, est = "Observational (OLS)", n = nrow(d), b = o[["b"]], se = o[["se"]], p = o[["p"]], F = NA, sargan_p = NA)
  C[[length(C) + 1]] <- data.table(outcome = ph, est = "2SLS, three adiposity scores", n = nrow(d), b = iv[["b"]], se = iv[["se"]], p = 2 * pnorm(-abs(iv[["b"]] / iv[["se"]])), F = fs$F[2], sargan_p = iv[["sargan_p"]])
  for (v in ADI) { f1 <- rob(lm(fm("BMI", c(v, COV)), d), v); rf <- rob(lm(fm("y", c(v, COV)), d), v)
    C[[length(C) + 1]] <- data.table(outcome = ph, est = paste0("Wald ratio, ", LAB[[v]], " score"), n = nrow(d), b = rf[["b"]] / f1[["b"]], se = abs(rf[["se"]] / f1[["b"]]), p = rf[["p"]], F = (f1[["b"]] / f1[["se"]])^2, sargan_p = NA) }
}
A <- rbindlist(A); B <- rbindlist(B); C <- rbindlist(C)
# d (JS 2026-11-14, Jon: liability-threshold reading): does the BMI slope depend on polygenic background? BMI x score interaction
# within carriers for each of the ten scores (continuous, per SD of score; robust SE), the autism-score interaction also adjusted for
# BMI x each adiposity score, and the BMI slope in carriers below / above the median autism score.
Dd <- list()
for (ph in c("SRS", "SCQ")) { d <- copy(W[pheno == ph])
  for (v in TEN) { d[, gB := get(v) * BMI]; r <- rob(lm(fm("y", c("BMI", v, "gB")), d), "gB")
    Dd[[length(Dd) + 1]] <- data.table(outcome = ph, term = paste(LAB[[v]], "score x BMI"), score = v, n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
  d[, gB := pgs_asd * BMI]; for (v in ADI) d[, paste0("bx_", v) := get(v) * BMI]
  r <- rob(lm(fm("y", c("BMI", TEN, "gB", paste0("bx_", ADI))), d), "gB")
  Dd[[length(Dd) + 1]] <- data.table(outcome = ph, term = "Autism score x BMI, adjusted for ten scores and BMI x adiposity scores", score = "pgs_asd_adj", n = nrow(d), b = r[["b"]], se = r[["se"]], p = r[["p"]])
  med <- median(d$pgs_asd)
  for (h in c("below", "above")) { dh <- d[if (h == "below") pgs_asd < med else pgs_asd >= med]; r <- rob(lm(fm("y", "BMI"), dh), "BMI")
    Dd[[length(Dd) + 1]] <- data.table(outcome = ph, term = sprintf("BMI slope, autism score %s median", h), score = paste0("half_", h), n = nrow(dh), b = r[["b"]], se = r[["se"]], p = r[["p"]]) }
}
Dd <- rbindlist(Dd)
tab <- rbind(A[, .(panel = "a BMI slope", outcome, term = set, n, b, se, p, F = NA_real_, expected = NA_real_, sargan_p = NA_real_)],
             B[, .(panel = "b score -> trait", outcome, term = LAB[score], n, b, se, p, F = fs_F, expected, sargan_p = NA_real_)],
             C[, .(panel = "c causal estimate", outcome, term = est, n, b, se, p, F, expected = NA_real_, sargan_p)],
             Dd[, .(panel = "d BMI x score within carriers", outcome, term, n, b, se, p, F = NA_real_, expected = NA_real_, sargan_p = NA_real_)])
fwrite(tab, file.path(HERE, paste0("fig_supp_bmi_srs_scq_pgs", SUF, "_table.csv")))
options(width = 220); print(tab[, .(panel, outcome, term, n, b = round(b, 3), se = round(se, 3), p = signif(p, 2), F = round(F, 1), expected = round(expected, 3), sargan_p = signif(sargan_p, 2))], nrows = 100)
# ---------------------------------------------------------------- figure
th <- theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5, colour = "grey30"),
  legend.position = "top", legend.key.size = unit(0.35, "cm"), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10), axis.title.y = element_blank())
pd <- position_dodge(width = 0.55); pf <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))
A[, set := factor(set, levels = rev(unique(set)))]; A[, outcome := factor(outcome, levels = c("SCQ", "SRS"))]
B[, outcome := factor(outcome, levels = c("SRS", "SCQ"))]
pa <- ggplot(A, aes(b, set, colour = outcome)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), width = 0, orientation = "y", position = pd) + geom_point(size = 2.2, position = pd) +
  scale_colour_manual(values = OC, name = NULL, breaks = c("SRS", "SCQ")) + labs(x = "s.d. of trait per s.d. of BMI (95% CI)", title = "a  BMI slope under adjustment",
  subtitle = sprintf("Carriers pooled across the four CNVs (SRS n = %d, SCQ n = %d); last two rows: carriers with measured height", A[outcome == "SRS"]$n[1], A[outcome == "SCQ"]$n[1])) + th
B[, lab := factor(LAB[score], levels = rev(LAB))]
pb <- ggplot(B, aes(b, lab)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se, colour = outcome), width = 0, orientation = "y") + geom_point(aes(colour = outcome), size = 2.2) +
  geom_point(aes(expected, lab), shape = 4, size = 2.6, stroke = 0.9) + facet_wrap(~outcome, nrow = 1) +
  scale_colour_manual(values = OC, guide = "none") + labs(x = "s.d. of trait per s.d. of polygenic score (95% CI)", title = "b  Polygenic scores and autistic traits",
  subtitle = "Dot, observed; cross, expected if the score acted on the trait only through BMI (first stage x observed BMI slope)") + th
C[, est := factor(est, levels = rev(unique(est)))]
XLC <- { m <- max(abs(c(C$b - 1.96 * C$se, C$b + 1.96 * C$se))); r <- min(3, ceiling(m * 2) / 2); c(-r, r) }
cl <- C[est == "2SLS, three adiposity scores", .(outcome, lbl = sprintf("%s (Sargan P = %s; first-stage F = %.1f)", outcome, pf(sargan_p), F))]
C[, panel := factor(cl$lbl[match(outcome, cl$outcome)], levels = cl$lbl)]
pc <- ggplot(C, aes(b, est)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se, colour = outcome), width = 0, orientation = "y") + geom_point(aes(colour = outcome), size = 2.2) +
  facet_wrap(~panel, nrow = 1) + scale_colour_manual(values = OC, guide = "none") + coord_cartesian(xlim = XLC) +
  labs(x = sprintf("s.d. of trait per s.d. of BMI (95%% CI; axis clipped at +/-%s)", XLC[2]), title = "c  Mendelian randomization",
  subtitle = "Adiposity scores as instruments for BMI; growth and neurodevelopmental scores as covariates") + th
DP <- Dd[!score %in% c("pgs_asd_adj")]
DP[, lab := fifelse(grepl("^half_", score), sub("BMI slope, autism score ", "BMI slope, autism score ", term), LAB[score])]
DP[, lab := factor(lab, levels = rev(c(LAB, "BMI slope, autism score below median", "BMI slope, autism score above median")))][, outcome := factor(outcome, levels = c("SRS", "SCQ"))]
DP[, hl := score %in% c("pgs_asd", "half_below", "half_above")]
dl <- Dd[score == "pgs_asd", .(outcome = factor(outcome, levels = c("SRS", "SCQ")), lbl = sprintf("autism score x BMI P = %s", pf(p)))]
pdd <- ggplot(DP, aes(b, lab)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 2.5, colour = "grey80", linewidth = 0.3) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se, colour = outcome, alpha = hl), width = 0, orientation = "y") + geom_point(aes(colour = outcome, alpha = hl), size = 2.2) +
  geom_text(data = dl, aes(x = Inf, y = Inf, label = lbl), hjust = 1.05, vjust = 1.5, size = 3, inherit.aes = FALSE) +
  facet_wrap(~outcome, nrow = 1) + scale_colour_manual(values = OC, guide = "none") + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  labs(x = "BMI x score interaction (s.d. of trait per s.d. of BMI per s.d. of score); bottom rows, BMI slope (95% CI)", title = "d  Does the BMI slope depend on polygenic background?",
  subtitle = "Within carriers: BMI x score interaction for each score; bottom, BMI slope in carriers below and above the median autism score") + th
fig <- pa / pb / pc / pdd + plot_layout(heights = c(1, 1.5, 0.8, 1.15))
ggsave(file.path(HERE, paste0("fig_supp_bmi_srs_scq_pgs", SUF, ".png")), fig, width = 10, height = 15.5, dpi = 250, bg = "white")
ggsave(file.path(HERE, paste0("fig_supp_bmi_srs_scq_pgs", SUF, ".pdf")), fig, width = 10, height = 15.5)
