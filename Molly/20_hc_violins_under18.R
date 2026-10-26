# Head circumference Z by BMI group, violins per CNV and non-carrier controls, ages 2 to < 18 (MS 2026-10-14).
# One record per person: visits at age >= 18 dropped, then each person's latest remaining visit (as in 19_hc_16p_under20.R).
# BMI split (centiles of the BMI Z reference): 16p11.2 del <= 90th vs > 90th; all other groups <= 30th vs >= 70th (middle dropped).
# Test: high - low BMI difference in mean HC Z (Welch t-test) and Wilcoxon rank-sum.
# Output: hc_violins_under18_tests.csv, fig_hc_violins_under18.png (age-binned version at the bottom overwrites both when BY_AGE is TRUE)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
AGE_MAX <- 18
GROUPS <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
SPLIT <- list("16p11.2 del" = c(lo = 90, hi = 90))   # others default to 30/70
split_of <- function(L) if (is.null(SPLIT[[L]])) c(lo = 30, hi = 70) else SPLIT[[L]]
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)))

d <- tab[!is.na(HEADCIRCUM) & !is.na(BMI) & !is.na(age) & age >= 2 & age < AGE_MAX][order(id, -age)][, .SD[1], by = id]

res <- list(); plots <- list()
for (L in GROUPS) {
  s <- split_of(L); lo <- s[["lo"]]; hi <- s[["hi"]]; gt <- if (hi == lo) ">" else "\u2265"
  lv <- c(sprintf("\u2264 %dth", lo), sprintf("%s %dth", gt, hi))
  x <- d[locus == L]
  x[, grp := factor(ifelse(if (hi == lo) BMI > qnorm(hi / 100) else BMI >= qnorm(hi / 100), lv[2], ifelse(BMI <= qnorm(lo / 100), lv[1], NA)), levels = lv)]
  x <- x[!is.na(grp)]
  tt <- t.test(HEADCIRCUM ~ grp, data = x); wt <- wilcox.test(HEADCIRCUM ~ grp, data = x, exact = FALSE)
  r <- data.table(locus = L, split = sprintf("<=%d / %s%d", lo, gt, hi), n_low = x[grp == lv[1], .N], n_high = x[grp == lv[2], .N],
                  mean_low = x[grp == lv[1], mean(HEADCIRCUM)], mean_high = x[grp == lv[2], mean(HEADCIRCUM)],
                  diff = unname(diff(tt$estimate)), p_welch = tt$p.value, p_wilcox = wt$p.value)
  res[[L]] <- r
  plots[[L]] <- ggplot(x, aes(grp, HEADCIRCUM, fill = grp)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_violin(colour = NA, alpha = 0.45, scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = grp), width = 0.12, height = 0, size = if (L == "No CNV") 0.4 else 0.9, alpha = if (L == "No CNV") 0.25 else 0.6, shape = 16) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.8) +
    annotate("text", x = 1.5, y = 4.9, vjust = 1, size = 3.4, label = sprintf("\u0394 %+.2f SD\nt p = %s, Wilcoxon p = %s", r$diff, pfx(r$p_welch), pfx(r$p_wilcox))) +
    scale_x_discrete(labels = c(sprintf("%s\n(n = %d)", lv[1], r$n_low), sprintf("%s\n(n = %d)", lv[2], r$n_high))) +
    scale_fill_manual(values = c("#00A5AD", "#F0645A"), guide = "none") + scale_colour_manual(values = c("#00A5AD", "#F0645A"), guide = "none") +
    coord_cartesian(ylim = c(-5, 5)) +
    labs(x = "BMI centile", y = "Head circumference (SD)", title = if (L == "No CNV") "Non-carrier controls" else L) +
    theme_classic(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
res <- rbindlist(res)
fwrite(res, file.path(OUT, "hc_violins_under18_tests.csv"))
print(res, digits = 3)

fig <- wrap_plots(plots, nrow = 1) +
  plot_annotation(title = sprintf("Head circumference by BMI group, ages 2-%d (one record per person, latest visit)", AGE_MAX - 1),
                  theme = theme(plot.title = element_text(size = 14)))
ggsave(file.path(OUT, "fig_hc_violins_under18.png"), fig, width = 17, height = 5, dpi = 200)

# ---------------------------------------------------------------- age-binned version (MS 2026-10-15)
# Each panel split into age bins on the x-axis, BMI groups side by side within each bin; test per bin.
# Set BY_AGE <- FALSE to get the single-violin-pair version above back.
BY_AGE <- TRUE
AGE_BREAKS <- c(2, 6, 12, 18); AGE_LABS <- c("2-5", "6-11", "12-17")
if (BY_AGE) {
  res <- list(); plots <- list()
  for (L in GROUPS) {
    s <- split_of(L); lo <- s[["lo"]]; hi <- s[["hi"]]; gt <- if (hi == lo) ">" else "≥"
    lv <- c(sprintf("BMI ≤ %dth", lo), sprintf("BMI %s %dth", gt, hi))
    x <- d[locus == L]
    x[, grp := factor(ifelse(if (hi == lo) BMI > qnorm(hi / 100) else BMI >= qnorm(hi / 100), lv[2], ifelse(BMI <= qnorm(lo / 100), lv[1], NA)), levels = lv)]
    x <- x[!is.na(grp)][, age_bin := cut(age, AGE_BREAKS, labels = AGE_LABS, right = FALSE)]
    st <- x[, {
      ok <- uniqueN(grp) == 2 && min(table(grp)) >= 2
      .(n_low = sum(grp == lv[1]), n_high = sum(grp == lv[2]),
        diff = if (ok) mean(HEADCIRCUM[grp == lv[2]]) - mean(HEADCIRCUM[grp == lv[1]]) else NA_real_,
        p_welch = if (ok) t.test(HEADCIRCUM ~ grp)$p.value else NA_real_,
        p_wilcox = if (ok) wilcox.test(HEADCIRCUM ~ grp, exact = FALSE)$p.value else NA_real_)
    }, by = age_bin][order(age_bin)]
    res[[L]] <- cbind(locus = L, split = sprintf("<=%d / %s%d", lo, gt, hi), st)
    st[, lab := ifelse(is.na(diff), "n/a", sprintf("Δ %+.2f\np = %s", diff, pfx(p_wilcox)))]
    st[, nlab := sprintf("%d | %d", n_low, n_high)]
    plots[[L]] <- ggplot(x, aes(age_bin, HEADCIRCUM, fill = grp)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_violin(colour = NA, alpha = 0.45, scale = "width", trim = FALSE, position = position_dodge(0.85), width = 0.8) +
      geom_point(aes(colour = grp), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85, seed = 1),
                 size = if (L == "No CNV") 0.3 else 0.8, alpha = if (L == "No CNV") 0.2 else 0.6, shape = 16) +
      geom_boxplot(aes(group = interaction(grp, age_bin)), width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8, position = position_dodge(0.85)) +
      geom_text(data = st, aes(age_bin, 5.6, label = lab), inherit.aes = FALSE, size = 3, vjust = 1, lineheight = 0.9) +
      geom_text(data = st, aes(age_bin, -5.4, label = nlab), inherit.aes = FALSE, size = 2.8, colour = "grey30") +
      scale_fill_manual(values = c("#00A5AD", "#F0645A")) + scale_colour_manual(values = c("#00A5AD", "#F0645A"), guide = "none") +
      coord_cartesian(ylim = c(-5.5, 5.6)) +
      labs(x = "Age (years)", y = "Head circumference (SD)", fill = NULL, title = if (L == "No CNV") "Non-carrier controls" else L) +
      theme_classic(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom",
                                            legend.text = element_text(size = 9), legend.margin = margin(0, 0, 0, 0))
  }
  res <- rbindlist(res)
  fwrite(res, file.path(OUT, "hc_violins_under18_tests.csv"))
  print(res, digits = 3)
  fig <- wrap_plots(plots, nrow = 1) +
    plot_annotation(title = sprintf("Head circumference by BMI group and age, ages 2-%d (one record per person, latest visit)", AGE_MAX - 1),
                    caption = "Numbers under each bin: n low BMI | n high BMI. Δ = high − low mean (SD); p = Wilcoxon rank-sum.",
                    theme = theme(plot.title = element_text(size = 14)))
  ggsave(file.path(OUT, "fig_hc_violins_under18.png"), fig, width = 20, height = 5.5, dpi = 200)
}
