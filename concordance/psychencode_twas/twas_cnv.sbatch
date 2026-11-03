#!/bin/bash
#SBATCH -J cnv_twas
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o /expanse/projects/sebat1/jsebat/cnv_twas_bmi/logs/twas_%A_%a.out
# FUSION TWAS of BMI (Yengo 2018) and height (Yengo 2022) on PsychENCODE cortex weights (PECfull), for the CNV-TWAS concordance
# (Sacks et al.); same settings as ConGA_GWAS/code/twas.sbatch, output kept outside the ConGA tree.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
B=/expanse/projects/sebat1/jsebat/ConGA_GWAS; O=/expanse/projects/sebat1/jsebat/cnv_twas_bmi; F=$B/tools/fusion/fusion_twas-master
T=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'BMI\nHEIGHT'); W=PECfull; mkdir -p $O/raw/$T; cd $F
for c in $(seq 1 22); do
  out=$O/raw/$T/$W.$c.dat; [ -s $out ] && continue
  /home/jsebat/micromamba/envs/conga/bin/Rscript FUSION.assoc_test.R --sumstats $O/sumstats/$T.sumstats \
    --weights $B/tools/fusion/weights/PECfull.pos --weights_dir $B/tools/fusion/weights --ref_ld_chr $B/tools/fusion/LDREF_PEC/1000G.EUR. \
    --max_impute 0.5 --chr $c --out $out > /dev/null 2> $O/raw/$T/$W.$c.log || echo "FAILED $T chr$c"
  echo "$T chr$c $(grep -vc '^PANEL' $out 2>/dev/null) genes"
done
echo "done $T"
