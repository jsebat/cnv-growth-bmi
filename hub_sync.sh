#!/bin/bash
# run from WSL: copy the Dropbox files listed in hub_filelist_*.txt into the Expanse hub, preserving relative paths
SRC="/mnt/c/Users/Jonathan Sebat/sebatlab Dropbox"
L="$SRC/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/CNV_growth_trajectories"
H=/expanse/projects/sebat1/jsebat/cnv_bmi_autistic_traits
cd "$SRC" || exit 1
for s in inputs derived code results; do rsync -a --files-from="$L/hub_filelist_$s.txt" ./ "expanse:$H/$s/" && echo "$s synced"; done
rsync -a "$L/DATA_MANIFEST.tsv" "expanse:$H/"
