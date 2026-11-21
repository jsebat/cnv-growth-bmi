#!/bin/bash
# run inside WSL from the expanse_results folder: copies aggregate tables and logs back (no individual-level data)
B=/expanse/projects/sebat1/jsebat/growth_tables_20260920/puberty_mr_v2_20261120
D=${1:-03_mr}
scp -q "expanse:$B/$D/*.csv" "expanse:$B/$D/run.log" "expanse:$B/01_verify/reproduction_check.csv" "expanse:$B/02_instruments/select.log" .
