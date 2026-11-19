# CNV x BMI and autistic traits: data hub

**Central location (Expanse, readable by jsebat-group):** `/expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits/`
Built 2026-11-14 (JS / Claude). Covers Fig. 5 (height, head circumference, SRS, SCQ by BMI in 16p11.2 / 22q11.2 carriers) and
Supplementary Figs. 7-8 of Sacks et al. (Google Doc 1G2AqOS6vkc03H64O0WCx1hu_RObLDHE9rl83t8P-l3o). A copy of this file is the hub's
`README.md`; the working copy of code and results stays in Dropbox
`Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/`.

Contains controlled participant-level data (G2MH, SSC, Simons VIP, Simons Searchlight). Do not copy outside approved storage or into
a public repository. Directories are setgid jsebat-group, mode 2750; files are not world-readable.

## Layout

| Folder | What | Files |
|---|---|---|
| `inputs/` | every raw input file the pipeline reads, copied with its Dropbox path below the Dropbox root (e.g. `inputs/G2MH/DataFreeze/Feb2026/...`) | 129 (107 MB) |
| `derived/` | participant-level tables built by the pipeline (`Jon/Claude/CNV_growth_trajectories_data/`): growth `analysis_table.rds`, `combined_srs_scq.rds` (combined SRS/SCQ sample), SSC polygenic scores, G2MH and Simons VIP autism scores (`g2mh_scores_raw/`), menarche scores, CHUSJ | 82 (17 MB) |
| `code/` | all scripts (`*.R`, `*.sh`, `*.md` handoffs) | 108 |
| `results/` | aggregate tables (`*.csv`), figures (`*.png`, `*.pdf`), panel objects (`*.RData`) | 217 (41 MB) |
| `expanse_links/` | symlinks to Expanse-native resources (not duplicated) | see below |
| `DATA_MANIFEST.tsv` | every file above with size and md5 (all 536 verified on Expanse 2026-11-14) | |

## Inputs and their canonical sources

| Data | Canonical source (Dropbox unless stated) |
|---|---|
| G2MH phenotypes, medications, dictionary (Feb 2026 freeze; site_id, ascertainment, SRS/SCQ items) | `G2MH/DataFreeze/Feb2026/` |
| Growth genomatrix (G2MH, SVIP, Searchlight, Cardiff, SSC: height, weight, HC, age, CNV) | `G2MH/NetworkPaper/AnalysisByMolly/Clinical/genomatrix/genomatrix.Clinical.02092026.csv` (built by `CreateGenomatrix2.R`) |
| Searchlight WGS ID map (searchlight_id SV##### -> sfari_iid) | same folder, `Searchlight.WGS2.2023_05.sample_metadata.tsv` |
| G2MH PGS panel (Anjali), WGS ancestry PCs | `G2MH/NetworkPaper/AnalysisByAnjali/data/pgs_scores/scores_wgs/`, `pca_results/naive/cohort.sscore` |
| SVIP-family genome PGS panel (Anjali; 393 genomes) | `AnalysisByAnjali/data/searchlight/pgs/analysis_dataset.tsv` |
| SSC PCs | `AnalysisByAnjali/data/all_cohorts_adsi_pgs.tsv` |
| Simons VIP SRS, SCQ, subjects, diagnoses | `G2MH/Simons_Searchlight/SVIP/SVIP/Simons_VIP__16p11.2_Dataset_v10.0/` |
| Simons Searchlight SRS, SCQ, height/weight (medical history), autism diagnosis | `G2MH/Simons_Searchlight/Simons_Searchlight_Dataset_v14.0/` (`all_cohorts/` and the 16p11.2 folders) |
| SSC phenotypes (SRS, SCQ, RBS-R, IQ, sibling height/weight, ages, ID map) | `Autism/SSC_Version_15_Phenotype_Data_Set_9/` (use `SSC IDs Sanders.csv` for roles) |
| G2MH autism diagnosis | `Jon/Claude/G2MH_CNV_growth/subjects.csv` |

## Expanse-native resources (`expanse_links/`)

| Link | Target | Contents |
|---|---|---|
| `pgs_scores_g2mh_svip_ssc` | `jsebat/pgs_obesity_ssc/scores` | April 2026 autism scores (G2MH WGS full weights; SVIP genomes full weights, job 55738525; SSC leave-SSC-out), SSC adiposity / height scores |
| `asd2026_sbayesrc_weights` | `jsebat/pgs_asd2026_loo/output/weights` | SBayesRC weights (full and leave-one-cohort-out) |
| `sumstats` | `jsebat/sumstats` | ASD 2026 GWAS summary statistics (full and LOO) |
| `svip_family_genomes_searchlight_393` | `a1sriniv/searchlight/pgs/results/03_missingness` | genotypes (pgen) of the 393 SVIP-family genomes |
| `fig5_build_main_with_ukbb` | `jsebat/growth_tables_20260920/main_with_ukbb` | Fig. 5 builds (`v6_boxcox`, `v7_combined`, `v8_adultbc` = current); UK Biobank adults enter here and stay on Expanse |

UK Biobank participant data are not copied into the hub (they stay in the Fig. 5 build directory under the UK Biobank access terms).

## Pipeline (current, 2026-11-14)

SRS/SCQ: `21_pheno_bmi_gam.R` (phenotypes, BMI z; adults Box-Cox) -> `23`, `24` (+ `00_boxcox_srs_scq.R` with `SRSSCQ_SCALE=boxcox`)
-> `60_combined_srs_scq_build.R` (combined G2MH + SVIP + Searchlight + SSC carriers -> `derived/.../combined_srs_scq.rds`) ->
`61` two-way, `62` three-way, `63` autism exclusion, `66` Supplementary Fig. 8, `49` Supplementary Fig. 7, `51b` Fig. 5k-t panels ->
Expanse `35d_main_figure_v8_adultbc.R` (Fig. 5). Full history and open items: `code/.../HANDOFF_2026-11-14.md`.

## Refreshing the hub

1. Rerun `70_data_hub_manifest.R` (Dropbox folder): rebuilds `DATA_MANIFEST.tsv` and `hub_filelist_*.txt`. Add any new input there.
2. From WSL: `bash hub_sync.sh` (rsync of the file lists into the hub), then on Expanse `bash hub_perm.sh` (group and mode) and
   `bash hub_verify.sh` (md5 check against the manifest). The three helper scripts are in `code/` and in
   `/expanse/projects/sebat1/jsebat/`.
