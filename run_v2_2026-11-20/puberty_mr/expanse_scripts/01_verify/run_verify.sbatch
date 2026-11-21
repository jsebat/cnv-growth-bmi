#!/bin/bash
#SBATCH -J pubmr_verify
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1 -n 1 --cpus-per-task=2 --mem=32G -t 01:00:00
#SBATCH -o verify.%j.log
# 2026-11-20 work package E: re-run frozen copies of the existing UK Biobank body-size / MR scripts in a new folder and
# compare every output CSV with the 2026-11-02 originals (growth_tables_20260920/ukbb_mr, ukbb_bodysize).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=$HOME/R/dnmt3a-r
W=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120/01_verify
cd $W
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
for s in ukbb_bodysize10_vs_adult ukbb_mr_childhood_obesity ukbb_mr_puberty_pathway fig_supp5_combined; do
  $RS $s.R > $s.log 2>&1 && echo "$s OK" || echo "$s FAILED"
done
O=/expanse/projects/sebat1/jsebat/growth_tables_20260920
$RS compare.R $O > compare.log 2>&1; cat compare.log
