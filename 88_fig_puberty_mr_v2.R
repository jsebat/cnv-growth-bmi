# 88_fig_puberty_mr_v2.R  (work package E, BMI modifier paper v2; 2026-11-20). Runs locally on aggregate tables only.
# What it does: draws the figure block "childhood body size, puberty and adult height" (lower part of main Figure 1):
#   i  adult height by recalled body size at age 10 in 16p11.2 del, 16p11.2 dup, 22q11.2 dup carriers and non-carriers
#      (cell means, 95% CI, linear trend per category with P; no individual UK Biobank points)
#   j  age at menarche by genotype (difference from non-carriers, years) and adult height per year of later menarche within genotype
#   k  path diagram of the Mendelian-randomisation result in non-carrier women
#   l  forest of the Mendelian-randomisation estimates behind k, by step and method, with instrument strength (F)
#   (mr_checks.csv in the output folder is the written verdict per reviewer check; its numbers come from the same input tables.)
# Inputs (aggregates made on Expanse by 88a_puberty_mr_ukbb_expanse.R, copied to run_v2_2026-11-20/puberty_mr/expanse_results):
#   panel_j_means.csv, panel_j_trends.csv, panel_k_menarche.csv, mr_onesample.csv, mr_twosample.csv, mr_controls.csv,
#   mr_sample_counts.csv, reproduction_check.csv
# Outputs (run_v2_2026-11-20/puberty_mr): block_puberty_mr.{png,pdf,csv}, panel_{i,j,k,l}.png (no tags)
# Panel tags i-l (coordinator, 2026-11-20: the block sits under panels a-h of main Figure 1); TAGS below re-letters everything., mr_checks.csv
# Data source of every number: UK Biobank (European-ancestry cluster); SNP-exposure effects in m from the named external GWAS.
# null model: none; this script only plots estimates and tests computed in 88a (regression Wald tests, delta-method / 2SLS SE).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
OUT <- file.path(HERE, "run_v2_2026-11-20/puberty_mr"); RES <- file.path(OUT, "expanse_results")
RED <- "#F0645A"; TEAL <- "#00A5AD"; CAR <- "#4B6EA9"; NON <- "#8C9A9A"; INK <- "grey15"
th <- theme_classic(base_size = 15) + theme(plot.title = element_text(size = 15, face = "bold"), plot.subtitle = element_text(size = 12.5, colour = "grey25"),
                                             strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold"), axis.text = element_text(colour = "grey15"))
GL <- c("16p11.2 DEL" = "16p11.2 del", "No CNV" = "Non-carriers", "16p11.2 DUP" = "16p11.2 dup", "22q11.2 DUP" = "22q11.2 dup")
sup <- function(k) { d <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079", "-" = "\u207b")
  paste(d[strsplit(as.character(k), "")[[1]]], collapse = "") }
pf <- function(p) sapply(p, function(x) { if (is.na(x)) return("")
  if (x < 1e-300) return(paste0("p < 10", sup(-300)))
  if (x >= 0.001) return(paste0("p = ", formatC(x, digits = if (x < 0.01) 1 else 2, format = "g")))
  e <- floor(log10(x)); m <- round(x / 10^e); if (m == 10) { m <- 1; e <- e + 1 }; paste0("p = ", m, " \u00d7 10", sup(e)) })
nf <- function(n) format(n, big.mark = ",", trim = TRUE)
mn <- function(x, d = 2) sub("-", "\u2212", formatC(x, digits = d, format = "f"), fixed = TRUE)

# ----------------------------------------------------------------------------------------------------- j
Jm <- fread(file.path(RES, "panel_j_means.csv")); Jt <- fread(file.path(RES, "panel_j_trends.csv"))[sexes == "all" & measure == "Body size at age 10"]
LV <- c("Thinner", "About average", "Plumper")
Jm[, x := match(body_size_at_10, LV)][, g := factor(GL[group], GL)]; Jt[, g := factor(GL[group], GL)]
Jt[, lab := paste0(mn(estimate), " per category\n", pf(p))][, face := ifelse(p < 0.05, "bold", "plain")][, nlab := paste0("n = ", nf(n))]
Jl <- Jm[, .(a = weighted.mean(mean_height_sd, n)), by = g]; Jl <- merge(Jl, Jt[, .(g, slope = estimate)], by = "g")
Jl <- merge(Jl, Jm[, .(xm = weighted.mean(x, n)), by = g], by = "g")                         # trend line through the group mean
yr <- range(c(Jm$lower95, Jm$upper95)); ytop <- yr[2] + 0.30 * diff(yr)
pj <- ggplot(Jm, aes(x, mean_height_sd)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_segment(data = Jl, aes(x = 0.8, xend = 3.2, y = a + slope * (0.8 - xm), yend = a + slope * (3.2 - xm)), colour = "grey35", linewidth = 0.7, linetype = "22", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95, colour = body_size_at_10), width = 0, linewidth = 1.1) + geom_point(aes(colour = body_size_at_10), size = 4) +
  geom_text(data = Jt, aes(x = 2, y = ytop, label = lab, fontface = face), size = 4.3, vjust = 1, lineheight = 0.95, inherit.aes = FALSE) +
  geom_text(data = Jt, aes(x = 2, y = yr[1] - 0.10 * diff(yr), label = nlab), size = 4.1, colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~g, nrow = 2) + scale_colour_manual(values = c(Thinner = TEAL, "About average" = NON, Plumper = RED), guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("Thinner", "Average", "Plumper"), limits = c(0.55, 3.45)) + coord_cartesian(ylim = c(yr[1] - 0.14 * diff(yr), ytop)) +
  labs(x = "Recalled body size at age 10", y = "Adult height (SD)", title = "Adult height by childhood body size", subtitle = "Mean and 95% CI; trend in SD of height") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.spacing = unit(0.9, "lines"))

# ----------------------------------------------------------------------------------------------------- k
K <- fread(file.path(RES, "panel_k_menarche.csv")); K[, g := factor(GL[group], rev(GL))]
K1 <- K[grepl("^menarche difference", quantity)]; K2 <- K[grepl("^adult height", quantity)]
K1[, lab := ifelse(group == "No CNV", "reference", paste0(mn(estimate, 1), " y, ", pf(p)))][, face := ifelse(!is.na(p) & p < 0.05, "bold", "plain")]
K2[, lab := paste0(mn(estimate), ", ", pf(p))][, face := ifelse(p < 0.05, "bold", "plain")]
gcol <- c("Non-carriers" = NON, "16p11.2 del" = CAR, "16p11.2 dup" = CAR, "22q11.2 dup" = CAR)
ylabs <- setNames(paste0(levels(K$g), "\n(n = ", nf(K1[match(levels(K$g), g), n_women]), ")"), levels(K$g))
fk <- function(D, xl, ttl, xlim) ggplot(D, aes(estimate, g, colour = g)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lower95, xmax = upper95), linewidth = 1.1, na.rm = TRUE) + geom_point(size = 4) +
  geom_text(aes(x = xlim[2], label = lab, fontface = face), colour = INK, size = 4.1, hjust = 1, vjust = -1.1) +
  scale_colour_manual(values = gcol, guide = "none") + scale_y_discrete(labels = ylabs, expand = expansion(add = c(0.45, 0.85))) + coord_cartesian(xlim = xlim) + labs(x = xl, y = NULL, title = ttl) + th +
  theme(plot.title = element_text(size = 14, face = "bold"), plot.title.position = "plot", axis.text.y = element_text(size = 12.5, lineheight = 0.9), axis.title.x = element_text(size = 13.5))
pk1 <- fk(K1, "Difference from non-carriers
(years)", "Age at menarche", c(-2.1, 2.6))
pk2 <- fk(K2, "Adult height (SD) per year
of later menarche", "Adult height by age at menarche", c(-0.25, 0.75))
pk <- pk1 / pk2

# ----------------------------------------------------------------------------------------------------- MR tables
O <- fread(file.path(RES, "mr_onesample.csv")); T2 <- fread(file.path(RES, "mr_twosample.csv")); CT <- fread(file.path(RES, "mr_controls.csv")); CN <- fread(file.path(RES, "mr_sample_counts.csv"))
A_PERRY <- "menarche instrument without UK Biobank: Perry 2014 allele score"; A_GUAN <- "main: childhood-obesity score, menarche score (Guan 2025), height score adjusted"
A_PERRY_NOHT <- "Perry 2014 allele score, without height-score adjustment"; A_PERRY_BMI <- "Perry 2014 allele score, plus adult-BMI score adjustment"
A_CBMI_PERRY <- "childhood-BMI score and Perry 2014 allele score (no UK Biobank in either)"
S_A <- "exposure score -> age at menarche"; S_FS <- "menarche instrument -> age at menarche (first stage)"; S_W <- "age at menarche -> adult height (Wald ratio)"
S_B <- "age at menarche -> adult height (2SLS, exposure score as covariate)"; S_OBS <- "age at menarche -> adult height (observational)"
S_TOT <- "exposure score -> adult height, total"; S_DIR <- "exposure score -> adult height, direct (menarche instrumented)"; S_IND <- "exposure score -> adult height, through menarche (a x b)"
g1 <- function(a, s) { r <- O[analysis == a & step == s]; stopifnot(nrow(r) == 1); r }
t2 <- function(a, m) { r <- T2[grepl(a, analysis, fixed = TRUE) & method == m]; stopifnot(nrow(r) >= 1); r }
ci <- function(r, d = 2) paste0("(", mn(r$lower95, d), " to ", mn(r$upper95, d), ")")
sg <- function(x, d) paste0(ifelse(x > 0, "+", ""), mn(x, d))
PRIMARY <- Sys.getenv("PRIMARY_INSTRUMENT", "perry")       # "perry" (allele score from a GWAS without UK Biobank) or "guan" (score trained in UK Biobank, as in the previous draft)
AP <- if (PRIMARY == "perry") A_PERRY else A_GUAN
nS <- g1(AP, S_A)$n

# ----------------------------------------------------------------------------------------------------- l: path diagram
ra <- g1(AP, S_A); rb <- g1(AP, S_B); rd <- g1(AP, S_DIR); ri <- g1(AP, S_IND); rt <- g1(AP, S_TOT)
bx <- data.table(x0 = c(0.1, 3.7, 7.5), x1 = c(3.4, 6.3, 9.9), y0 = c(1.35, 7.75, 1.35), y1 = c(3.25, 9.25, 3.25),
                 lab = c("Childhood-obesity\npolygenic score", "Age at\nmenarche", "Adult\nheight"), fill = c("#FBE3E0", "#E4E8EF", "#E4E8EF"))
ar <- arrow(length = unit(0.16, "inches"), type = "closed")
pl <- ggplot() + geom_rect(data = bx, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill), colour = "grey25", linewidth = 0.5) + scale_fill_identity() +
  geom_text(data = bx, aes((x0 + x1) / 2, (y0 + y1) / 2, label = lab), size = 4.5, fontface = "bold", lineheight = 0.95) +
  annotate("segment", x = 2.3, y = 3.4, xend = 4.25, yend = 7.6, arrow = ar, linewidth = 1.3, colour = RED) +
  annotate("segment", x = 5.75, y = 7.6, xend = 7.8, yend = 3.4, arrow = ar, linewidth = 1.3, colour = CAR) +
  annotate("segment", x = 3.55, y = 2.3, xend = 7.35, yend = 2.3, arrow = ar, linewidth = 1.3, colour = "grey35") +
  annotate("text", x = 2.75, y = 6.1, hjust = 1, size = 4.1, lineheight = 1.0, colour = INK, label = paste0("Earlier menarche\n", mn(ra$estimate), " years\nper SD of score\n", ci(ra), "\n", pf(ra$p))) +
  annotate("text", x = 7.25, y = 6.1, hjust = 0, size = 4.1, lineheight = 1.0, colour = INK, label = paste0("Later menarche,\ntaller adult\n", sg(rb$estimate, 3), " SD\nper year later\n", ci(rb, 3), "\n", pf(rb$p))) +
  annotate("text", x = 5.45, y = 3.5, size = 4.1, lineheight = 1.0, colour = INK, label = paste0("Direct effect\n", sg(rd$estimate, 3), " SD per SD\n", ci(rd, 3), "\n", pf(rd$p))) +
  annotate("text", x = 5.0, y = 6.05, size = 4.0, lineheight = 1.0, colour = "grey25", label = paste0("Through\nmenarche\n", sg(ri$estimate, 3), " SD\n", ci(ri, 3))) +
  annotate("text", x = 5.0, y = 0.4, size = 4.1, lineheight = 1.0, colour = INK, label = paste0("Total effect on adult height: ", sg(rt$estimate, 3), " SD per SD of score\n", ci(rt, 3), ", ", pf(rt$p))) +
  coord_cartesian(xlim = c(0, 10), ylim = c(-0.25, 9.4), expand = FALSE) +
  labs(title = "Mendelian randomisation in non-carriers", subtitle = paste0(nf(nS), " women, UK Biobank")) + theme_void(base_size = 15) +
  theme(plot.title = element_text(size = 15, face = "bold"), plot.subtitle = element_text(size = 12.5, colour = "grey25"), plot.margin = margin(5.5, 8, 5.5, 8))

# ----------------------------------------------------------------------------------------------------- m: forest
Fq <- function(a, s) nf(round(O[analysis == a & step == s, F_stat]))
mk <- function(block, label, r, src, kind = "score") data.table(block = block, label = label, estimate = r$estimate, lower95 = r$lower95, upper95 = r$upper95, p = r$p, kind = kind, source = src)
V <- "childhood BMI (Vogelezang 2020)"; P <- "age at menarche (Perry 2014) -> adult height"
kV <- t2(paste0(V, " -> age at menarche"), "Inverse-variance weighted"); kP <- t2(P, "Inverse-variance weighted"); kM <- T2[grepl("^multivariable: childhood BMI", analysis)]
B1 <- "Childhood body size → age at menarche  (years per SD of score)"
B2 <- paste0("Childhood BMI → age at menarche  (years per SD; ", kV$n_snps, " SNPs, mean F = ", round(kV$mean_F), ")")
B3 <- paste0("Age at menarche → adult height  (SD per year later; SNPs: mean F = ", round(kP$mean_F), ")")
B4 <- "Childhood body size → adult height, direct  (SD per SD of score)"
B5 <- paste0("Childhood BMI → adult height, direct  (SD per SD; ", kM$n_snps[1], " SNPs)")
Fm <- rbind(
  mk(B1, paste0("Childhood-obesity score (F = ", Fq(A_PERRY, S_A), ")"), g1(A_PERRY, S_A), "mr_onesample.csv"),
  mk(B1, paste0("Childhood-BMI score (F = ", Fq(A_CBMI_PERRY, S_A), ")"), g1(A_CBMI_PERRY, S_A), "mr_onesample.csv"),
  mk(B2, "Inverse-variance weighted", kV, "mr_twosample.csv", "snp"),
  mk(B2, "Weighted median", t2(paste0(V, " -> age at menarche"), "Weighted median"), "mr_twosample.csv", "snp"),
  mk(B2, "MR-Egger", t2(paste0(V, " -> age at menarche"), "MR-Egger"), "mr_twosample.csv", "snp"),
  mk(B3, paste0("Allele score, ", kP$n_snps, " SNPs (F = ", Fq(A_PERRY, S_FS), ")"), g1(A_PERRY, S_B), "mr_onesample.csv"),
  mk(B3, "  same, height score not adjusted", g1(A_PERRY_NOHT, S_B), "mr_onesample.csv"),
  mk(B3, paste0("Score trained in UK Biobank (F = ", Fq(A_GUAN, S_FS), ")"), g1(A_GUAN, S_B), "mr_onesample.csv", "overlap"),
  mk(B3, paste0(kP$n_snps, " SNPs, inverse-variance weighted"), kP, "mr_twosample.csv", "snp"),
  mk(B3, paste0(kP$n_snps, " SNPs, weighted median"), t2(P, "Weighted median"), "mr_twosample.csv", "snp"),
  mk(B3, paste0(kP$n_snps, " SNPs, MR-Egger"), t2(P, "MR-Egger"), "mr_twosample.csv", "snp"),
  mk(B3, "Multivariable, with childhood BMI", kM[grepl("per year later", unit)], "mr_twosample.csv", "snp"),
  mk(B4, "Childhood-obesity score", g1(A_PERRY, S_DIR), "mr_onesample.csv"),
  mk(B4, "  same, height score not adjusted", g1(A_PERRY_NOHT, S_DIR), "mr_onesample.csv"),
  mk(B4, "Childhood-BMI score", g1(A_CBMI_PERRY, S_DIR), "mr_onesample.csv"),
  mk(B5, "Multivariable, with age at menarche", kM[grepl("direct", unit)], "mr_twosample.csv", "snp"))
Fm[, block := factor(block, c(B1, B2, B3, B4, B5))][, ord := .I]
fcol <- c(score = INK, snp = CAR, overlap = NON)
fm1 <- function(b) { D <- Fm[block == b]; D[, lab := factor(label, rev(label))]; d <- if (max(abs(c(D$lower95, D$upper95))) < 0.1) 3 else 2
  D[, txt := paste0(mn(estimate, d), "  (", mn(lower95, d), " to ", mn(upper95, d), ")")]
  rg <- range(c(D$lower95, D$upper95, 0)); w <- diff(rg); xt <- rg[2] + 0.06 * w
  ggplot(D, aes(estimate, lab, colour = kind)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_linerange(aes(xmin = lower95, xmax = upper95), linewidth = 1) + geom_point(size = 3.3) +
    geom_text(aes(x = xt, label = txt), hjust = 0, size = 4.0, colour = INK) + scale_colour_manual(values = fcol, guide = "none") +
    scale_x_continuous(breaks = function(l) { b <- scales::breaks_pretty(2)(rg); b[b >= rg[1] - 0.03 * w & b <= rg[2] + 0.02 * w] }, labels = function(x) sub("-", "−", format(x, drop0trailing = TRUE, trim = TRUE), fixed = TRUE)) +
    coord_cartesian(xlim = c(rg[1] - 0.03 * w, rg[2] + 1.7 * w), clip = "off") + labs(x = NULL, y = NULL, title = b) + theme_classic(base_size = 15) +
    theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0), plot.title.position = "plot", axis.text.y = element_text(size = 12.5, colour = "grey10", hjust = 0), axis.text.x = element_text(size = 11.5),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(3, 6, 2, 4)) }
nb <- Fm[, .N, by = block]$N
pm <- wrap_plots(lapply(levels(Fm$block), fm1), ncol = 1, heights = nb + 0.9)

# ----------------------------------------------------------------------------------------------------- assemble and save
TAGS <- c("i", "j", "k", "l")                      # body size, menarche, path diagram, forest
tg <- theme(plot.tag = element_text(face = "bold", size = 22))
W <- c(4.3, 3.7, 5.3, 6.7); H <- as.numeric(Sys.getenv("BLOCK_H", "7"))
pkt <- wrap_elements(full = (pk1 + labs(tag = TAGS[2]) + tg) / pk2)
pmt <- wrap_elements(full = wrap_plots(c(list(fm1(levels(Fm$block)[1]) + labs(tag = TAGS[4]) + tg), lapply(levels(Fm$block)[-1], fm1)), ncol = 1, heights = nb + 0.9))
blk <- wrap_elements(full = pj + labs(tag = TAGS[1]) + tg) + pkt + wrap_elements(full = pl + labs(tag = TAGS[3]) + tg) + pmt + plot_layout(nrow = 1, widths = W)
sfx <- if (PRIMARY == "perry") "" else "_ukb_trained_instrument"
ggsave(file.path(OUT, paste0("block_puberty_mr", sfx, ".png")), blk, width = 20, height = H, dpi = 200, bg = "white")
ggsave(file.path(OUT, paste0("block_puberty_mr", sfx, ".pdf")), blk, width = 20, height = H, device = cairo_pdf, bg = "white")
pn <- function(k, s = "") file.path(OUT, paste0("panel_", TAGS[k], s, ".png"))
if (PRIMARY == "perry") {
  ggsave(pn(1), pj, width = W[1], height = H, dpi = 200, bg = "white")
  ggsave(pn(2), pk, width = W[2], height = H, dpi = 200, bg = "white")
  ggsave(pn(4), pm, width = W[4], height = H, dpi = 200, bg = "white")
}
ggsave(pn(3, sfx), pl, width = W[3], height = H, dpi = 200, bg = "white")

# every plotted number
if (PRIMARY == "perry") {
  blkcsv <- rbind(
    Jm[, .(panel = TAGS[1], group = GL[group], item = paste("mean adult height (SD), body size at 10:", body_size_at_10), n, estimate = mean_height_sd, lower95, upper95, p = NA_real_, source = "panel_j_means.csv; UK Biobank, both sexes")],
    Jt[, .(panel = TAGS[1], group = GL[group], item = "trend: adult height (SD) per category of body size at 10 (thinner < average < plumper)", n, estimate, lower95, upper95, p, source = "panel_j_trends.csv; UK Biobank, both sexes")],
    K[, .(panel = TAGS[2], group = GL[group], item = quantity, n = n_women, estimate, lower95, upper95, p, source = "panel_k_menarche.csv; UK Biobank women")],
    rbindlist(lapply(list(ra, rb, rd, ri, rt), function(r) data.table(panel = TAGS[3], group = "Non-carriers", item = paste0(r$step, " [", r$unit, "]"), n = r$n, estimate = r$estimate, lower95 = r$lower95, upper95 = r$upper95, p = r$p, source = paste0("mr_onesample.csv; ", r$analysis)))),
    Fm[, .(panel = TAGS[4], group = "Non-carriers", item = paste0(block, " | ", trimws(label)), n = nS, estimate, lower95, upper95, p, source = paste0(source, "; UK Biobank non-carrier women"))])
  fwrite(blkcsv, file.path(OUT, "block_puberty_mr.csv"))
}
cat("figure written\n")
