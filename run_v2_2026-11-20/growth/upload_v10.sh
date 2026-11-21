#!/bin/bash
# Upload the v10 inputs to Expanse through ssh (scp is closed by the login banner). Run inside WSL from the laptop.
set -euo pipefail
D="/mnt/c/Users/Jonathan Sebat/sebatlab Dropbox"
H="$D/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
R=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb/v10_v2draft_2026-11-20/upload
up() { ssh -o BatchMode=yes expanse "cat > $R/$(basename "$1")" < "$1"; md5sum "$1" | cut -d' ' -f1; }
up "$H/82_growth_blocks_v2.R"
up "$H/run_v2_2026-11-20/growth/expanse_setup_v10.sh"
up "$D/Jon/Claude/CNV_growth_trajectories_data/analysis_table_chop_overlapkeep_2026-11-20.rds"
ssh -o BatchMode=yes expanse "cd $R && md5sum *"
