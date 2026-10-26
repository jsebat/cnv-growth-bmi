# Supplementary figure: prospective, age-varying effect of BMI on subsequent growth (the "top row" analysis) for
# height and head circumference, carriers < 18 y.
#   Height: from 05_age_varying.R (age_varying_beta.csv, model "Prospective"); cleaned longitudinal visits of SVIP,
#           Searchlight and Cardiff (genomatrix Z).
#   Head circumference: from 07_hc_trajectories.R (hc_prospective_age_varying.csv); full SVIP v10 series + Cardiff
#           (UK90 Z), so it reaches down to infancy for 16p11.2.
# Model (both): change in Z to the next visit ~ Z at the visit + interval + sex (+ locus, cohort) + s(age) +
#   s(age, by = BMI Z at the visit) + child random intercept; the curve is the BMI coefficient as a function of age.
# 22q11.2 dup is omitted (12 children with repeat height, 6 with repeat HC).
# Output: fig_supp_prospective.{png,pdf}

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
source(file.path(HERE, "00_helpers.R"))
COLS <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "Pooled")   # pooled on the right (JS)

h <- fread(file.path(HERE, "age_varying_beta.csv"), encoding = "UTF-8")[outcome == "HEIGHT" & grepl("^Prospective", model)]
# interval counts for height (same data as 05: cleaned visits, intervals >= 0.5 y)
ivh <- intervals(long_data("HEIGHT"))
nh <- rbind(ivh[, .(children = uniqueN(id), intervals = .N), by = .(locus = as.character(locus))],
            ivh[, .(locus = "Pooled", children = uniqueN(id), intervals = .N)])
h <- merge(h[, .(age, beta, se, locus, p = p_smooth)], nh, by = "locus")[, outcome := "Height"]
hc <- fread(file.path(HERE, "hc_prospective_age_varying.csv"))[, outcome := "Head circumference"]
cv <- rbind(h, hc[, .(locus, age, beta, se, p, children, intervals, outcome)])[locus %in% COLS]
cv[, `:=`(locus = factor(locus, COLS), outcome = factor(outcome, c("Height", "Head circumference")))]
# All-ages (constant) BMI coefficient, same data and covariates as the age-varying models but s(age, by = BMI)
# replaced by a single linear BMI term.  Shown on the head-circumference row (JS: 22q11.2 del HC looks flat but
# consistently positive, which the "differs from 0 at any age" test is not built to detect).
suppressPackageStartupMessages(library(mgcv))
const_fit <- function(ii) {
  ii <- copy(ii)[, `:=`(fid = factor(id), cohort = droplevels(factor(cohort)), locus = droplevels(factor(locus)), sex = factor(sex))]
  adj <- paste0(if (nlevels(ii$locus) > 1) "+ locus " else "", if (nlevels(ii$cohort) > 1) "+ cohort" else "")
  kk <- min(6, uniqueN(round(ii$age0, 2)) - 1)
  m <- gam(as.formula(sprintf("dZ ~ B0 + Z0 + dt + sex %s + s(age0, k = %d) + s(fid, bs = 're')", adj, kk)), data = ii, method = "REML")
  co <- summary(m)$p.table["B0", ]; data.table(beta_all = co[1], se_all = co[2], p_all = co[4])
}
hcv <- readRDS(file.path(DATA, "hc_visits.rds"))
ivc <- hcv[!is.na(bmiz)][order(id, age)][, .(age0 = head(age, -1), dt = diff(age), Z0 = head(hcz, -1), dZ = diff(hcz), B0 = head(bmiz, -1),
                                            sex = sex[1], locus = locus[1], cohort = cohort[1]), by = id][dt >= 0.1]   # as in 07
cst <- rbindlist(lapply(c(setdiff(COLS, "Pooled"), "Pooled"), function(L) rbind(
  cbind(data.table(locus = L, outcome = "Height"), const_fit(if (L == "Pooled") ivh else ivh[locus == L])),
  cbind(data.table(locus = L, outcome = "Head circumference"), const_fit(if (L == "Pooled") ivc else ivc[locus == L])))))
fwrite(cst, file.path(HERE, "supp_prospective_all_ages.csv")); print(cst, digits = 3)
lab <- unique(cv[, .(locus, outcome, p, children, intervals)])
lab <- merge(lab[, `:=`(locus = as.character(locus), outcome = as.character(outcome))], cst, by = c("locus", "outcome"))
lab[, txt := sprintf("age-varying p = %s\n%d children, %d intervals", format.pval(p, digits = 2, eps = 1e-4), children, intervals)]
lab[outcome == "Head circumference", txt := sprintf("%s\nall ages: %+.2f ± %.2f, p = %s", txt, beta_all, se_all, format.pval(p_all, digits = 2, eps = 1e-4))]
lab[, `:=`(locus = factor(locus, COLS), outcome = factor(outcome, c("Height", "Head circumference")))]

fig <- ggplot(cv, aes(age, beta)) +
  annotate("rect", xmin = 8, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), fill = "#F0645A", alpha = 0.25) +
  geom_line(colour = "#C0392B", linewidth = 0.9) +
  geom_text(data = lab, aes(x = 0, y = Inf, label = txt), hjust = 0, vjust = 1.3, size = 4.1) +
  facet_grid(outcome ~ locus, scales = "free_y") + coord_cartesian(xlim = c(0, 17)) +
  labs(x = "Age at visit (years)", y = "Change in Z to the next visit per 1 SD of BMI Z (\u00b1 95% CI)",
       title = "Supplementary: does higher BMI at a visit precede faster growth? Age-varying prospective effect in CNV carriers < 18 y",
       subtitle = paste("Change in height Z (top) or head-circumference Z (bottom) to the next visit per 1 SD of BMI Z at a visit, as a smooth function of age at that visit.",
                        "\nMixed GAM with child random intercept, adjusted for Z at the visit, interval length, sex, locus and cohort; age-varying p = BMI effect differs from 0 at any age; all ages = single constant BMI coefficient (± SE). Grey = 8-13 y.",
                        "\nHeight: SVIP, Searchlight, Cardiff (cleaned visits, intervals \u2265 0.5 y). Head circumference: full SVIP v10 series + Cardiff (UK90 Z). 22q11.2 dup omitted (\u2264 12 children).")) +
  theme_bw(base_size = 15) + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5))
ggsave(file.path(HERE, "fig_supp_prospective.png"), fig, width = 15, height = 8.5, dpi = 200)
ggsave(file.path(HERE, "fig_supp_prospective.pdf"), fig, width = 15, height = 8.5, device = cairo_pdf)
print(lab)
