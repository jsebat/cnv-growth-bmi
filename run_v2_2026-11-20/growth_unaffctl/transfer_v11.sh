#!/bin/bash
# Laptop <-> Expanse transfer for the v11 runs, through ssh (scp is closed by the login banner). Run inside WSL.
#   bash transfer_v11.sh up                 upload 82, the setup script and the two corrected participant tables
#   bash transfer_v11.sh down [variant ...] copy AGGREGATE outputs back (tables, figure blocks, logs; never .rds, never UK Biobank rows)
set -euo pipefail
D="/mnt/c/Users/Jonathan Sebat/sebatlab Dropbox"
H="$D/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
G="$H/run_v2_2026-11-20/growth_unaffctl"
R=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v11_unaffected_controls_2026-11-20
if [ "$1" = "up" ]; then
  ssh -o BatchMode=yes expanse "mkdir -p $R/upload"
  up() { ssh -o BatchMode=yes expanse "cat > $R/upload/$(basename "$1")" < "$1"; }
  up "$H/82_growth_blocks_v2.R"; up "$G/expanse_setup_v11.sh"
  for f in analysis_table_chop_unaffctl_2026-11-20.rds analysis_table_chop_unaffctl_overlapkeep_2026-11-20.rds; do
    ssh -o BatchMode=yes expanse "test -f $R/upload/$f" || up "$D/Jon/Claude/CNV_growth_trajectories_data/$f"
  done
  ssh -o BatchMode=yes expanse "cd $R/upload && md5sum *"
else
  shift; VARS="${*:-main sens_visit_latest sens_no_chop sens_chop_overlap_keep sens_no_chop_g2mh203_kept sens_controls_ssc_probands}"
  for v in $VARS; do
    mkdir -p "$G/$v"
    ssh -o BatchMode=yes expanse "cd $R/$v && tar cf - --exclude=growth_refs_lms.csv --exclude=supp_bars_results.csv *.csv block_*.png block_*.pdf fig_supp_main_effects_by_band.* v11_*.log run.sbatch inputs.md5 2>/dev/null || true" | tar xf - -C "$G/$v"
    echo "$v: $(ls "$G/$v" | wc -l) files"
  done
fi
