# Is the negative within-person height-BMI coupling a measurement artifact?
# BMI = weight / height^2, so an error in measured height moves height Z and BMI Z in opposite directions.
# Weight does not use height, so a real "more adiposity -> faster growth" effect should show up with weight too,
# while the artifact should not.  Weight Z is sample-internal: log(weight) standardised against non-carrier
# controls < 18 y (SSC/SVIP/G2MH) by sex and age (GAM for the mean, GAM on |residual| for the spread).
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(lmerTest); library(splines) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
DATA <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"   # participant-level tables stay in Jon/Claude (controlled data)
tab <- readRDS(file.path(DATA, "analysis_table.rds"))
tab[, W := BMI_raw * (HEIGHT_raw / 100)^2]
tab <- tab[age < 18 & !is.na(sex) & !is.na(W) & W > 0]
ref <- tab[locus == "No CNV"]
m_mu <- gam(log(W) ~ sex + s(age, by = sex, k = 12), data = ref)
ref[, ar := abs(log(W) - predict(m_mu, newdata = ref))]
m_sd <- gam(ar ~ sex + s(age, by = sex, k = 8), data = ref)
tab[, WZ := (log(W) - predict(m_mu, newdata = tab)) / (predict(m_sd, newdata = tab) * sqrt(pi / 2))]
cat("controls weight Z mean/sd:", round(tab[locus == "No CNV", c(mean(WZ), sd(WZ))], 2), "\n")

d <- tab[locus != "No CNV" & !is.na(HEIGHT) & !is.na(BMI)]
d <- unique(d[order(id, age)], by = c("id", "age")); d <- d[, if (.N >= 2) .SD, by = id]
for (v in c("HEIGHT", "BMI", "WZ")) d[, paste0(v, "_w") := get(v) - mean(get(v)), by = id]
out <- list()
for (s in list(c("SVIP", "Searchlight"), "Cardiff")) for (L in c("all", "16p11.2 del")) {
  x <- d[cohort %in% s & (L == "all" | locus == L)]
  out[[length(out) + 1]] <- data.table(cohorts = paste(s, collapse = "+"), locus = L, persons = uniqueN(x$id), visits = nrow(x),
    sd_w_height = sd(x$HEIGHT_w), sd_w_BMI = sd(x$BMI_w), sd_w_weight = sd(x$WZ_w),
    r_height_BMI = cor(x$HEIGHT_w, x$BMI_w), r_height_weight = cor(x$HEIGHT_w, x$WZ_w), r_BMI_weight = cor(x$BMI_w, x$WZ_w))
}
out <- rbindlist(out); print(out, digits = 3); fwrite(out, file.path(HERE, "artifact_check_correlations.csv"))

# within-person model with weight in place of BMI, and with both (BMI given weight = the height-denominator part)
d[, WZb := mean(WZ), by = id]
for (f in list(HEIGHT ~ BMI_w + sex * ns(age, 3) + locus + cohort + (1 | id),
               HEIGHT ~ WZ_w + WZb + sex * ns(age, 3) + locus + cohort + (1 | id),
               HEIGHT ~ BMI_w + WZ_w + sex * ns(age, 3) + locus + cohort + (1 | id))) {
  m <- lmer(f, data = d); co <- summary(m)$coefficients
  cat("\n", deparse(f[[3]][[2]]), "\n"); print(round(co[grep("_w$", rownames(co)), c(1, 2, 5), drop = FALSE], 4))
}
