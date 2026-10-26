# Age-varying effects of BMI and height on head circumference, per CNV (MS 2026-10-15).
# Varying-coefficient GAM, fitted per group (non-carriers and each CNV), one record per person (latest visit before 18 with
# height, BMI and HC; ages 2-17), REML:
#   HC ~ s(age) + s(age, by = BMI) + s(age, by = HEIGHT)
# s(age, by = BMI) = slope of HC Z on BMI Z at each age, holding height fixed (and vice versa for height).
# Stats per panel:
#   varies with age p -- ML likelihood-ratio test of the varying-coefficient model vs constant slopes (HC ~ s(age) + BMI + HEIGHT),
#                        dropping the age variation for that trait only
#   vs controls p     -- carriers + non-carriers together: ... + s(age, by = BMI) + s(age, by = BMI x carrier) (same for height);
#                        p of the carrier term = does the carriers' age-varying slope differ from non-carriers' (level or shape)
# Output: fig_hc_varying_coefficients.png, hc_varying_coefficients_tests.csv, hc_varying_coefficients_curves.csv

suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2) })
DBX <- if (dir.exists("C:/Users/Jonathan Sebat/sebatlab Dropbox")) "C:/Users/Jonathan Sebat/sebatlab Dropbox" else path.expand("~/sebatlab Dropbox")   # Jon (Windows) or Molly (Mac)
HERE <- file.path(DBX, "Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories")   # shared inputs (Jon's results)
OUT <- file.path(HERE, "Molly")   # Molly's exploratory outputs
DATA <- Sys.getenv("GROWTH_DATA", file.path(DBX, "Jon/Claude/CNV_growth_trajectories_data"))   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
AGE_MAX <- 18
K <- 6
GROUPS <- c("No CNV", "16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
GLAB <- c("No CNV" = "Non-carrier controls", "16p11.2 del" = "16p11.2 BP4-BP5 DEL", "16p11.2 dup" = "16p11.2 BP4-BP5 DUP",
          "22q11.2 del" = "22q11.2 A-D DEL", "22q11.2 dup" = "22q11.2 A-D DUP")
TERMS <- c(BMI = "BMI \u2192 HC slope\n(SD HC per SD BMI, height fixed)", HEIGHT = "Height \u2192 HC slope\n(SD HC per SD height, BMI fixed)")
pfx <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))))

d <- tab[!is.na(HEIGHT) & !is.na(BMI) & !is.na(HEADCIRCUM) & !is.na(age) & age >= 2 & age < AGE_MAX][order(id, -age)][, .SD[1], by = id]
d[, HC := HEADCIRCUM]

# slope curve f(age) for the by-variable term: lpmatrix columns of "s(age):<v>" evaluated with v = 1
slope_curve <- function(m, x, v) {
  ag <- seq(min(x$age), max(x$age), length.out = 100)
  nd <- data.frame(age = ag, BMI = 0, HEIGHT = 0); nd[[v]] <- 1
  X <- predict(m, nd, type = "lpmatrix"); cols <- grep(sprintf("^s\\(age\\):%s\\.", v), colnames(X))
  f <- drop(X[, cols] %*% coef(m)[cols]); se <- sqrt(rowSums((X[, cols] %*% m$Vp[cols, cols]) * X[, cols]))
  data.table(age = ag, est = f, lo = f - 1.96 * se, hi = f + 1.96 * se)
}

curves <- list(); res <- list()
for (L in GROUPS) {
  x <- d[locus == L]
  m <- gam(HC ~ s(age, k = K) + s(age, by = BMI, k = K) + s(age, by = HEIGHT, k = K), data = x, method = "REML")
  mML <- gam(HC ~ s(age, k = K) + s(age, by = BMI, k = K) + s(age, by = HEIGHT, k = K), data = x, method = "ML")
  for (v in names(TERMS)) {
    curves[[paste(L, v)]] <- slope_curve(m, x, v)[, `:=`(locus = L, term = v)]
    other <- setdiff(names(TERMS), v)
    m0 <- gam(as.formula(sprintf("HC ~ s(age, k = %d) + %s + s(age, by = %s, k = %d)", K, v, other, K)), data = x, method = "ML")   # constant slope for v
    p_vary <- anova(m0, mML, test = "Chisq")$`Pr(>Chi)`[2]
    p_ctrl <- NA_real_
    if (L != "No CNV") {
      xx <- copy(d[locus %in% c("No CNV", L)])[, `:=`(cnv = as.numeric(locus == L), grp = factor(locus == L, levels = c(FALSE, TRUE), ordered = TRUE))]
      xx[, `:=`(BMI_c = BMI * cnv, HEIGHT_c = HEIGHT * cnv)]
      mc <- gam(HC ~ grp + s(age, k = K) + s(age, by = grp, k = K) + s(age, by = BMI, k = K) + s(age, by = HEIGHT, k = K) +
                  s(age, by = BMI_c, k = K) + s(age, by = HEIGHT_c, k = K), data = xx, method = "REML")
      st <- summary(mc)$s.table; p_ctrl <- st[sprintf("s(age):%s_c", v), "p-value"]
    }
    res[[paste(L, v)]] <- data.table(locus = L, term = v, n = nrow(x), mean_slope = mean(curves[[paste(L, v)]]$est), p_varies_with_age = p_vary, p_vs_controls = p_ctrl)
  }
}
curves <- rbindlist(curves); res <- rbindlist(res)
fwrite(res, file.path(OUT, "hc_varying_coefficients_tests.csv")); fwrite(curves, file.path(OUT, "hc_varying_coefficients_curves.csv"))
print(res, digits = 3)

for (z in list(curves, res)) z[, `:=`(locus_lab = factor(GLAB[locus], levels = GLAB), term_lab = factor(TERMS[term], levels = TERMS))]
res[, lab := gsub("p = < ", "p < ", paste0(sprintf("varies with age p = %s", pfx(p_varies_with_age)),
                                           ifelse(locus == "No CNV", "", sprintf("\nvs controls p = %s", pfx(p_vs_controls)))), fixed = TRUE)]
res[, nlab := sprintf("n = %d", n)]
ref <- rbindlist(lapply(GLAB[-1], function(g) copy(curves[locus == "No CNV"])[, locus_lab := factor(g, levels = GLAB)]))   # non-carrier curve in every column

fig <- ggplot(curves, aes(age, est)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_line(data = ref, colour = "grey45", linetype = "22", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1F3A5F", alpha = 0.2) +
  geom_line(colour = "#1F3A5F", linewidth = 0.9) +
  geom_text(data = res, aes(x = AGE_MAX, y = -Inf, label = lab), hjust = 1, vjust = -0.4, size = 3.2, lineheight = 0.9, inherit.aes = FALSE) +
  geom_text(data = res, aes(x = 2.5, y = Inf, label = nlab), hjust = 0, vjust = 1.4, size = 3.2, inherit.aes = FALSE) +
  facet_grid(term_lab ~ locus_lab, switch = "y", scales = "free_y") +
  coord_cartesian(xlim = c(2, AGE_MAX)) +
  labs(x = "Age (years)", y = NULL,
       title = "Age-varying effects of BMI and height on head circumference",
       subtitle = sprintf("HC ~ s(age) + s(age, by = BMI) + s(age, by = height), per group; one record per person, ages 2-%d. Solid: slope \u00b1 95%% CI; dashed grey: non-carriers.", AGE_MAX - 1)) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(10, "pt"), plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_hc_varying_coefficients.png"), fig, width = 18, height = 7.5, dpi = 200)
