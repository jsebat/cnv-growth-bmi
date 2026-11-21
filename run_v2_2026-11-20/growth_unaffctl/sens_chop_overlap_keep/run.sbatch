#!/bin/bash
#SBATCH -J v11_sens_chop_overlap_keep
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nice=100
#SBATCH -N 1 -n 1 --cpus-per-task=4 --mem=32G -t 02:00:00
#SBATCH -o /expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v11_unaffected_controls_2026-11-20/sens_chop_overlap_keep/v11_sens_chop_overlap_keep.%j.log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=$HOME/R/dnmt3a-r
export ANALYSIS_TABLE=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v11_unaffected_controls_2026-11-20/sens_chop_overlap_keep/analysis_table_chop_unaffctl_overlapkeep_2026-11-20.rds VISIT_RULE=fill_gaps VARIANT=sens_chop_overlap_keep DROP_COHORT= MAKE_BLOCKS=0
cd /expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v11_unaffected_controls_2026-11-20/sens_chop_overlap_keep
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
md5sum *.R *.rds *.csv > inputs.md5
$RS 82_growth_blocks_v2.R && echo STEP82_OK
if [ "0" = "1" ]; then   # Supplementary Fig. 6 on the job's own five-group table written by 82
  export ANALYSIS_TABLE=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v11_unaffected_controls_2026-11-20/sens_chop_overlap_keep/analysis_table_job.rds DROP_COHORT=
  $RS 43_supp_main_effects_by_band.R && echo STEP43_OK
fi
