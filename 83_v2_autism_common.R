# 83_v2_autism_common.R  (JS / Claude 2026-11-20; BMI modifier paper v2, work package C: autistic traits and autism diagnosis)
# What it does: shared helpers for 84-87. (1) Re-runs the existing scripts (23, 24, 51b, 61-66, 32) from their own code, unchanged,
#   with every file they write redirected to run_v2_2026-11-20/autism/aux_rerun (nothing in HERE or in older run folders is
#   overwritten). (2) Figure style of the v2 brief (theme_classic 15, bold lowercase tags size 22, colours, "p = " labels).
#   (3) Small statistics helpers (robust SE, Cochran Q heterogeneity).
# Inputs : the numbered scripts in HERE; participant-level tables in DATA (read only).
# Outputs: none by itself (sourced by 84-87).
# null model: none; helpers only (each analysis script states its own test).
suppressPackageStartupMessages({ library(data.table); library(mgcv); library(ggplot2); library(patchwork); library(sandwich) })
HERE <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
V2OUT <- file.path(HERE, "run_v2_2026-11-20", "autism"); V2AUX <- file.path(V2OUT, "aux_rerun")
dir.create(V2AUX, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(OUT_DIR = V2AUX, ADD_CHOP = "TRUE", VISIT_RULE = "fill_gaps", UNAFF_DRAW_FIG = "FALSE", COG_DRAW_FIG = "FALSE")

# ---- re-running older scripts without touching their outputs
v2_lines <- function(f) {
  L <- readLines(f, warn = FALSE, encoding = "UTF-8")
  w <- grepl("fwrite\\(|ggsave\\(|save\\(list|saveRDS\\(|writeLines\\(|saveWorkbook\\(", L)
  L[w] <- gsub("file.path(HERE,", "file.path(V2AUX,", L[w], fixed = TRUE); L
}
v2_eval <- function(L) invisible(eval(parse(text = L, encoding = "UTF-8"), envir = globalenv()))
v2_source <- function(f, from = 1L, to = NULL) { L <- v2_lines(f); if (is.null(to)) to <- length(L); v2_eval(L[from:to]) }
# nested source() calls inside the older scripts go through the same redirection
source <- function(file, ...) v2_source(file)
v2_line <- function(f, pattern, which = 1L) grep(pattern, readLines(f, warn = FALSE, encoding = "UTF-8"))[which]

# ---- style (V2_BRIEF.md)
C_HI <- "#F0645A"; C_LO <- "#00A5AD"; C_MID <- "#B8BEC4"; C_CAR <- "#4B6EA9"; C_NON <- "#8C9A9A"; C_F <- "#C2528B"; C_M <- "#3C7FB1"
LOC4 <- c("16p11.2 del", "16p11.2 dup", "22q11.2 del", "22q11.2 dup")
th_v2 <- function(bs = 15) theme_classic(base_size = bs) + theme(plot.tag = element_text(size = 22, face = "bold"),
  plot.title = element_text(face = "bold", size = bs), plot.subtitle = element_text(size = bs - 3, colour = "grey25"), legend.background = element_blank())
pfmt <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), ifelse(p < 0.01, formatC(p, format = "f", digits = 4), formatC(p, format = "f", digits = 3))))
plab <- function(p, pre = "p") ifelse(is.na(p), paste(pre, "= n/a"), paste(pre, "=", pfmt(p)))
face <- function(p) ifelse(!is.na(p) & p < 0.05, "bold", "plain")
save_v2 <- function(fig, name, w, h) { ggsave(file.path(V2OUT, paste0(name, ".png")), fig, width = w, height = h, dpi = 200, bg = "white")
  tryCatch(ggsave(file.path(V2OUT, paste0(name, ".pdf")), fig, width = w, height = h, device = cairo_pdf), error = function(e) ggsave(file.path(V2OUT, paste0(name, ".pdf")), fig, width = w, height = h)) }

# ---- statistics helpers
robHC <- function(m, v) { V <- vcovHC(m, type = "HC1"); b <- unname(coef(m)[v]); s <- sqrt(V[v, v]); c(b = b, se = s, p = 2 * pnorm(-abs(b / s))) }
term_of <- function(m, parts) { nm <- names(coef(m)); nm[sapply(strsplit(nm, ":"), function(z) setequal(z, parts))][1] }
cochranQ <- function(b, se) { ok <- !is.na(b) & !is.na(se); b <- b[ok]; se <- se[ok]; if (length(b) < 2) return(c(Q = NA, df = NA, p = NA, I2 = NA))
  w <- 1 / se^2; bb <- sum(w * b) / sum(w); Q <- sum(w * (b - bb)^2); df <- length(b) - 1; c(Q = Q, df = df, p = pchisq(Q, df, lower.tail = FALSE), I2 = max(0, (Q - df) / Q)) }
BAND3 <- c("BMI \u2264 30th centile", "BMI 30th-70th centile", "BMI \u2265 70th centile")
band3 <- function(z) factor(fifelse(z >= qnorm(.7), BAND3[3], fifelse(z <= qnorm(.3), BAND3[1], BAND3[2])), BAND3)

# ---- base build (about 4 min): re-runs 51b (SRS / SCQ panels), the data block of 32 (autism diagnosis by CNV, with CHOP, fill_gaps)
# and 66 (CNV x BMI and CNV x BMI x autism score; via 64c / 64 / 63 / 62 / 61) and keeps their objects under new names:
#   XP, RP  : 51b panel data (relatives, SSC autistic probands, combined carriers; Box-Cox scale, y_adj) and its result table
#   A, XA   : growth / phenotype table of 23 and the diagnosis table of 32 (Height, Head circumference, BMI rows; asd)
#   C       : combined SRS / SCQ sample (combined_srs_scq.rds) with asd2 (63) and site (64)
#   XG, R66 : genotyped subset with the autism polygenic score (62) and the result table of 66
#   S, idmap, sl : score table (id, pgs, src), growth-table id map, Searchlight autism status (63)
# The participant-level cache stays in the controlled data folder (DATA/v2_autism_cache); V2_REBUILD=TRUE forces a rebuild.
v2_base <- function(rebuild = Sys.getenv("V2_REBUILD", "FALSE") == "TRUE") {
  DATA0 <- "C:/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/CNV_growth_trajectories_data"
  cf <- file.path(DATA0, "v2_autism_cache", "base_2026-11-20.RData")
  if (file.exists(cf) && !rebuild) { load(cf, envir = globalenv()); cat("base objects loaded from cache\n"); return(invisible(TRUE)) }
  f51 <- file.path(HERE, "51b_cog_panels_boxcox_combined.R"); v2_source(f51, to = v2_line(f51, "^# local check") - 1)
  assign("XP", copy(get("X", globalenv())), globalenv()); assign("RP", copy(get("R", globalenv())), globalenv())
  f32 <- file.path(HERE, "32_supp_asd_by_cnv.R"); v2_source(f32, from = v2_line(f32, "^pro_g <- "), to = v2_line(f32, "^X <- rbind[(]X, XB[)]"))
  assign("XA", copy(get("X", globalenv())), globalenv())
  v2_eval("idmap <- unique(rbind(gmc[, .(id, cohort, sampleID, GUID)], ctl_g[, .(id, cohort, sampleID, GUID)]))")
  v2_source(file.path(HERE, "66_fig_supp_interactions_main.R"))
  assign("XG", copy(get("X", globalenv())), globalenv()); assign("R66", copy(get("R", globalenv())), globalenv())
  dir.create(dirname(cf), showWarnings = FALSE, recursive = TRUE)
  save(list = c("XP", "RP", "LAM", "A", "XA", "C", "XG", "R66", "S", "idmap", "sl", "vs_ref", "three", "rob", "RELf", "SIBf", "PARf", "COHL"), file = cf, envir = globalenv())
  cat("base objects built and cached\n"); invisible(TRUE)
}
