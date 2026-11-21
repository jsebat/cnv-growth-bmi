B=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120
squeue -u jsebat; cat $B/pubmr_v2.*.log; cat $B/02_instruments/select.log; ls $B/02_instruments | wc -l; tail -3 $B/02_instruments/plink_chr1.out 2>/dev/null; ls $B/03_mr*; tail -${1:-30} $B/03_mr/run.log 2>/dev/null
