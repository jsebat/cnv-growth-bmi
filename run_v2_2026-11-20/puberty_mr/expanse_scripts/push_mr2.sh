B=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120
mkdir -p $B/03_mr_r2 && cp /tmp/pubmr_up/88a_puberty_mr_ukbb_expanse.R $B/03_mr_r2/
cd $B && sbatch --export=ALL,STAGES=03,OUT03=03_mr_r2 -o pubmr_v2_r2.%j.log run_mr.sbatch
