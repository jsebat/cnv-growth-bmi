# Fig. 1F-G rebuilt with faint grey 95% CI whiskers on both axes (JS 2026-11-03; reviewer request).
# Molly's original plotting code could not be found (Dropbox, Drive, GitHub), so the panels are redrawn from
# Supplementary Table 3 (meta-analysis DEL and DUP effects, SEs, P) for the 34 reciprocal loci.
#  F (height): two parallel dose-response lines from a 2-component mixture regression (common slope, EM),
#    which reproduces the lines of the original panel. Its printed statistics (Rsq 0.712, P 8.64e-18) could not
#    be reproduced from these data and are kept verbatim from the original panel (flagged to JS).
#  G (BMI): one least-squares line; Rsq 0.532, P 1.48e-06 reproduce the original panel exactly.
# Colours: P < 0.001 for both alleles (purple), deletion only (red), duplication only (blue), neither (black);
# only coloured points are labelled; 16p11.2 and 22q11.2 loci labelled in bold. Bands: 95% CI of the fitted lines.
# Axis ranges as in the original panels; whiskers that extend past them are clipped at the panel edge.
# Output: fig1FG_whiskers.png (panels only, letters added when composited into the Figure 1 image).
suppressMessages({ library(openxlsx); library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
t3 <- as.data.table(read.xlsx(file.path(HERE, "manuscript_materials_2026-10-24/Sacks Supplementary Tables v2.xlsx"), sheet = "T3 CNV main effects and frequen"))
t3 <- t3[!is.na(locus)]
num <- c("count", "betaHeight", "seHeight", "pHeight", "betaBMI", "seBMI", "pBMI")
t3[, (num) := lapply(.SD, as.numeric), .SDcols = num]
t3[, type := substr(locus, nchar(locus) - 2, nchar(locus))][, region := substr(locus, 1, nchar(locus) - 4)]
w <- merge(t3[type == "DEL", c("region", num), with = FALSE], t3[type == "DUP", c("region", num), with = FALSE],
           by = "region", suffixes = c("_DEL", "_DUP"))
stopifnot(nrow(w) == 34)
pretty_name <- function(r) {
  s <- gsub("_", " ", r)
  s <- gsub("(BP[0-9]) (BP[0-9])", "\\1-\\2", s)
  s <- gsub("^(22q11\\.21) ([A-H]) ([A-H])$", "\\1 \\2-\\3", s)
  s
}
w[, label := pretty_name(region)]
w[, bold := grepl("^(16p11\\.2|22q11\\.21)", region)]
COLS <- c(both = "#7b2d8e", `deletion only` = "#d62728", `duplication only` = "#1f4e9c", neither = "black")

# two-component mixture regression for height (common slope, separate intercepts, equal variance; EM, 50 starts)
em2 <- function(x, y, iters = 2000) {
  best <- NULL
  for (s in 1:50) {
    set.seed(s); z <- runif(length(x)); pr <- cbind(z, 1 - z)
    for (it in 1:iters) {
      X <- rbind(cbind(1, 0, x), cbind(0, 1, x)); b <- lm.wfit(X, c(y, y), c(pr[, 1], pr[, 2]))$coefficients
      mu <- cbind(b[1] + b[3] * x, b[2] + b[3] * x); s2 <- sum(pr * (y - mu)^2) / length(y); p1 <- mean(pr[, 1])
      d <- cbind(p1 * dnorm(y, mu[, 1], sqrt(s2)), (1 - p1) * dnorm(y, mu[, 2], sqrt(s2))); ll <- sum(log(rowSums(d))); pr <- d / rowSums(d)
    }
    if (is.null(best) || ll > best$ll) best <- list(ll = ll, b = b, pr = pr)
  }
  best
}

panel <- function(tr, title, stats_txt, two_lines, xl, yl) {
  d <- data.table(region = w$region, label = w$label, bold = w$bold,
                  x = w[[paste0("beta", tr, "_DEL")]], y = w[[paste0("beta", tr, "_DUP")]],
                  sx = w[[paste0("se", tr, "_DEL")]], sy = w[[paste0("se", tr, "_DUP")]],
                  px = w[[paste0("p", tr, "_DEL")]], py = w[[paste0("p", tr, "_DUP")]])
  d <- d[is.finite(x) & is.finite(y)]   # 7q11.23 distal has no BMI deletion estimate
  d[, sig := fifelse(px < 1e-3 & py < 1e-3, "both", fifelse(px < 1e-3, "deletion only", fifelse(py < 1e-3, "duplication only", "neither")))]
  d[, sig := factor(sig, levels = names(COLS))]
  grid <- data.table(x = seq(min(d$x), max(d$x), length.out = 200))
  if (two_lines) {
    f <- em2(d$x, d$y); d[, cl := factor(max.col(f$pr))]
    m <- lm(y ~ x + cl, d)
    band <- rbindlist(lapply(levels(d$cl), function(k) { g <- copy(grid)[, cl := factor(k, levels = levels(d$cl))]
      pp <- predict(m, g, interval = "confidence"); g[, `:=`(fit = pp[, 1], lo = pp[, 2], hi = pp[, 3])] }))
  } else {
    m <- lm(y ~ x, d); pp <- predict(m, grid, interval = "confidence"); band <- copy(grid)[, `:=`(fit = pp[, 1], lo = pp[, 2], hi = pp[, 3], cl = factor(1))]
  }
  lab <- d[sig != "neither"]
  ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi, group = cl), inherit.aes = FALSE, fill = "grey75", alpha = 0.55) +
    geom_line(data = band, aes(x = x, y = fit, group = cl), inherit.aes = FALSE, colour = "grey45", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = x - 1.96 * sx, xmax = x + 1.96 * sx), height = 0, colour = "grey80", linewidth = 0.3) +
    geom_errorbar(aes(ymin = y - 1.96 * sy, ymax = y + 1.96 * sy), width = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = sig), size = 2.2) +
    geom_text_repel(data = lab, aes(label = label, fontface = ifelse(bold, "bold", "plain")), size = 2.3, colour = "black",
                    min.segment.length = 0.3, segment.colour = "grey60", segment.size = 0.2, box.padding = 0.25, max.overlaps = Inf, seed = 1) +
    coord_cartesian(xlim = xl, ylim = yl) +
    scale_colour_manual(values = COLS, drop = FALSE, name = "p < 0.001") +
    annotate("text", x = Inf, y = Inf, label = stats_txt, hjust = 1.05, vjust = 1.3, size = 3.2, lineheight = 1.1) +
    labs(title = title, x = sprintf("Deletion Effect Size (%s s.d.)", tr), y = sprintf("Duplication Effect Size (%s s.d.)", tr)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 10.5), legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
          legend.background = element_blank(), legend.key.size = unit(8, "pt"), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4))
}
pF <- panel("Height", "Height DUP vs. DEL effect size", "Rsq =  0.712\nP-value =  8.64e-18", TRUE, c(-1.12, 0.42), c(-0.47, 0.42))
mB <- lm(betaBMI_DUP ~ betaBMI_DEL, w); sB <- summary(mB)
pG <- panel("BMI", "BMI DUP vs. DEL effect size",
            sprintf("Rsq =  %.3f\nP-value =  %s", sB$r.squared, formatC(anova(mB)[1, 5], format = "e", digits = 2)), FALSE, c(-0.56, 1.02), c(-0.62, 0.66))
png(file.path(HERE, "fig1FG_whiskers.png"), width = 1446, height = 603, res = 150)
print(pF + pG)
dev.off()
cat("height mixture classes and BMI stats written\n")
