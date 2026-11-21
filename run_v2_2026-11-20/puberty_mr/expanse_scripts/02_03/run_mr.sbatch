#!/bin/bash
#SBATCH -J pubmr_v2
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1 -n 1 --cpus-per-task=4 --mem=64G -t 03:00:00
#SBATCH -o pubmr_v2.%j.log
# 2026-11-20 work package E: (02) pick SNP instruments from GWAS without UK Biobank, extract their UK Biobank genotypes,
# (03) run 88a (figure aggregates, one-sample score MR, two-sample SNP MR, controls). STAGES env: "02 03" (default) or "03".
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=$HOME/R/dnmt3a-r
B=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
STAGES=${STAGES:-"02 03"}; OUT03=${OUT03:-03_mr}
if [[ "$STAGES" == *02* ]]; then
  cd $B/02_instruments
  $RS 88b_puberty_mr_select_instruments_expanse.R > select.log 2>&1 && echo STEP02_OK || { echo STEP02_FAILED; exit 1; }
  for i in {1..22}; do
    [ -s extract_chr$i.txt ] || continue
    /home/j3guevar/bin/plink2 --bfile /expanse/projects/sebat1/j3guevar/archive/ukb_imp_chrN_bgen1/ukb_imp_chr${i}_v3_qc_uniq \
      --extract extract_chr$i.txt --export A --threads 4 --memory 40000 --out dosage_chr$i > plink_chr$i.out 2>&1 || echo "plink chr$i FAILED"
  done
  ls dosage_chr*.raw | wc -l
fi
mkdir -p $B/$OUT03; cd $B/$OUT03
export INSTR_DIR=$B/02_instruments
$RS 88a_puberty_mr_ukbb_expanse.R > run.log 2>&1 && echo STEP03_OK || echo STEP03_FAILED
tail -5 run.log
