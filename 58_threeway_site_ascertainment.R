# Site and ascertainment sensitivity of the within-carrier BMI x autism score interaction (and of the BMI slope) on SRS / SCQ
# (JS / Claude 2026-11-14; Jon: could the effect be related to cohort?). The scored carriers are all G2MH (no Simons VIP); G2MH has
# seven sites (site_id), records whether the carrier is the family proband (first member ascertained) and why genetic testing was
# done (dem_ascertainment_reason: 1 medical, 2 developmental/learning, 3 mental health, 4 family member with the variant,
# 5 prenatal, 6 research). Setup and models as 57 (Box-Cox scale, robust SE).
#  - composition: carriers by site x CNV, proband status, ascertainment reason
#  - interaction with site fixed effects, + site x BMI and site x score; leave-one-site-out; site heterogeneity (site x BMI x score)
#  - by proband status and by ascertainment (developmental / mental-health reason vs not; cascade = family member reason);
#    interaction adjusted for ascertainment x BMI and x score; difference test (group x BMI x score)
# Output: threeway_site_ascertainment.csv (aggregate statistics only).
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
s57 <- readLines(file.path(HERE, "57_cnv_bmi_pgs_threeway.R")); eval(parse(text = s57[1:grep("^out <- list\\(\\); add", s57)]))
ra <- fread(file.path(DB, "G2MH/DataFreeze/Feb2026/g2mh_project2_data_release_022026.csv"), select = c("guid", "site_id", "proband", paste0("dem_ascertainment_reason___", 1:6)))
setnames(ra, paste0("dem_ascertainment_reason___", 1:6), c("asc_medical", "asc_developmental", "asc_mental", "asc_family", "asc_prenatal", "asc_research"))
ra[, proband := as.numeric(proband)]; ra <- ra[, lapply(.SD, function(v) if (is.character(v)) v[v != ""][1] else suppressWarnings(as.numeric(max(v, na.rm = TRUE)))), by = guid]
W <- merge(W, ra, by.x = "key", by.y = "guid", all.x = TRUE)
for (v in c("asc_medical", "asc_developmental", "asc_mental", "asc_family", "asc_prenatal", "asc_research")) W[!is.finite(get(v)), (v) := NA]
W[, asc_neurodev := as.numeric(asc_developmental %in% 1 | asc_mental %in% 1)][is.na(asc_developmental) & is.na(asc_mental), asc_neurodev := NA]
cat("composition (SRS rows):\n"); print(dcast(W[pheno == "SRS", .N, by = .(site_id, locus)], site_id ~ locus, value.var = "N", fill = 0))
print(W[pheno == "SRS", .(n = .N, proband = sum(proband %in% 1), nonproband = sum(proband %in% 0), proband_NA = sum(is.na(proband)),
  developmental = sum(asc_developmental %in% 1), mental = sum(asc_mental %in% 1), family = sum(asc_family %in% 1), medical = sum(asc_medical %in% 1), prenatal = sum(asc_prenatal %in% 1), reason_NA = sum(is.na(asc_neurodev)))])
for (ph in c("SRS", "SCQ")) {
  d <- W[pheno == ph & !is.na(site_id)]; d[, site := factor(site_id)]
  B <- "site"
  add(ph, B, "main (same rows)", ix(d))
  add(ph, B, "+ site", ix(d, "site"))
  add(ph, B, "+ site, site x BMI, site x score", ix(d, c("site", "site:BMI", "site:pgs_asd")))
  for (s in levels(d$site)) add(ph, B, paste("leave out", s), ix(d[site != s]))
  big <- d[, .N, by = site][N >= 40, site]; dh <- d[site %in% big]; dh[, site := droplevels(site)]
  m1 <- lm(y ~ locus + locus:poly(age, 2) + male + site * BMI * pgs_asd, dh); m0 <- lm(y ~ locus + locus:poly(age, 2) + male + site * BMI + site * pgs_asd + BMI:pgs_asd, dh)
  add(ph, B, sprintf("site heterogeneity of the interaction (F test, %d sites with n >= 40)", length(big)), c(n = nrow(dh), b = NA, se = NA, p = anova(m0, m1)$`Pr(>F)`[2]))
  for (s in big) { ds <- d[site == s]; add(ph, B, paste("only", s), ix(ds, base = if (uniqueN(ds$locus) > 1) BASEW else "poly(age, 2) + male")) }
  B <- "ascertainment"
  da <- W[pheno == ph]
  add(ph, B, "family proband (first member ascertained)", ix(da[proband %in% 1]))
  add(ph, B, "not the family proband (ascertained through a relative)", ix(da[proband %in% 0]))
  dp <- da[!is.na(proband)]; m1 <- lm(y ~ locus + locus:poly(age, 2) + male + factor(proband) * BMI * pgs_asd, dp); r <- rob(m1, "factor(proband)1:BMI:pgs_asd")
  add(ph, B, "difference proband vs not (proband x BMI x score)", c(n = nrow(dp), r))
  add(ph, B, "tested for a developmental or mental-health problem", ix(da[asc_neurodev %in% 1]))
  add(ph, B, "tested for other reasons (medical, family member, prenatal, research)", ix(da[asc_neurodev %in% 0]))
  dn <- da[!is.na(asc_neurodev)]; m1 <- lm(y ~ locus + locus:poly(age, 2) + male + asc_neurodev * BMI * pgs_asd, dn); r <- rob(m1, "asc_neurodev:BMI:pgs_asd")
  add(ph, B, "difference developmental/mental-health vs other (x BMI x score)", c(n = nrow(dn), r))
  add(ph, B, "+ proband and reason, each x BMI and x score", ix(da[!is.na(proband) & !is.na(asc_neurodev)], c("factor(proband)", "asc_neurodev", "factor(proband):BMI", "asc_neurodev:BMI", "factor(proband):pgs_asd", "asc_neurodev:pgs_asd")))
  add(ph, B, "tested because a family member carried the variant", ix(da[asc_family %in% 1]))
  # BMI slope itself (the CNV x BMI result) by site and ascertainment, for completeness
  B <- "BMI slope by site / ascertainment"
  bs <- function(x, base = BASEW) { m <- lm(as.formula(paste("y ~", base, "+ BMI")), x); c(n = nrow(x), rob(m, "BMI")) }
  add(ph, B, "all", bs(da)); add(ph, B, "+ site", bs(d, paste(BASEW, "+ site")))
  for (s in big) { ds <- d[site == s]; add(ph, B, paste("only", s), bs(ds, if (uniqueN(ds$locus) > 1) BASEW else "poly(age, 2) + male")) }
  add(ph, B, "family proband", bs(da[proband %in% 1])); add(ph, B, "not the family proband", bs(da[proband %in% 0]))
  add(ph, B, "developmental/mental-health reason", bs(da[asc_neurodev %in% 1])); add(ph, B, "other reasons", bs(da[asc_neurodev %in% 0]))
}
R <- rbindlist(out); fwrite(R, file.path(HERE, "threeway_site_ascertainment.csv"))
options(width = 230); print(R[, .(trait, block, spec = substr(spec, 1, 75), n, b = round(b, 3), se = round(se, 3), p = signif(p, 2))], nrows = 200)
