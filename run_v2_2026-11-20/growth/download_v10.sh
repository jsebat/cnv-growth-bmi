#!/bin/bash
# Copy the aggregate outputs of the v10 runs back from Expanse (result tables, figure blocks, logs; never the participant
# tables or any UK Biobank individual data). Run inside WSL from the laptop. Usage: bash download_v10.sh [variant ...]
set -euo pipefail
G="/mnt/c/Users/Jonathan Sebat/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories/run_v2_2026-11-20/growth"
R=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v10_v2draft_2026-11-20
VARS="${*:-main sens_visit_latest sens_no_chop sens_chop_overlap_keep sens_no_chop_g2mh203_kept}"
for v in $VARS; do
  mkdir -p "$G/$v"
  ssh -o BatchMode=yes expanse "cd $R/$v && tar cf - --exclude=growth_refs_lms.csv --exclude=supp_bars_results.csv *.csv block_*.png block_*.pdf v10_*.log run.sbatch inputs.md5 2>/dev/null || true" | tar xf - -C "$G/$v"
  echo "$v: $(ls "$G/$v" | wc -l) files"
done
