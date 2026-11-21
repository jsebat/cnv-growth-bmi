B=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120
cp /tmp/pubmr_up/88b_puberty_mr_select_instruments_expanse.R $B/02_instruments/
cp /tmp/pubmr_up/88a_puberty_mr_ukbb_expanse.R $B/03_mr/
cp /tmp/pubmr_up/run_mr.sbatch $B/
cd $B && sbatch run_mr.sbatch
squeue -u jsebat; cat $B/01_verify/verify.*.log 2>/dev/null | tail -20
