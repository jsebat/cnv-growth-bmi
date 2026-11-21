#!/bin/bash
#SBATCH -J v10_main_r3
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nice=100
#SBATCH -N 1 -n 1 --cpus-per-task=4 --mem=32G -t 02:00:00
#SBATCH -o /expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v10_v2draft_2026-11-20/main_r3/v10_main_r3.%j.log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=$HOME/R/dnmt3a-r
export ANALYSIS_TABLE=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v10_v2draft_2026-11-20/main_r3/analysis_table_chop_2026-11-19.rds VISIT_RULE=fill_gaps VARIANT=main_r3 DROP_COHORT= MAKE_BLOCKS=1
cd /expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v10_v2draft_2026-11-20/main_r3
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
md5sum *.R *.rds *.csv > inputs.md5
$RS 82_growth_blocks_v2.R && echo STEP82_OK
