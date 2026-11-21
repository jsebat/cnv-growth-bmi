#!/bin/bash
# v10_v2draft_2026-11-20 (JS / Claude 2026-11-20): growth figure blocks for the BMI modifier paper v2.
# One run folder per job variant, each with its own frozen copies of the scripts and inputs. v9e is only read.
# Usage: bash expanse_setup_v10.sh <tag>   (tag is appended to the variant folder names; "" for the first submission)
set -euo pipefail
TAG="${1:-}"; ONLY="${2:-}"   # optional second argument: submit only this variant
BASE=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb
SRC=$BASE/v9e_chop_genotype_bands
ROOT=$BASE/v10_v2draft_2026-11-20
IN=$ROOT/inputs
mkdir -p "$IN"
# inputs: unchanged scripts and tables from v9e; 82 and the CHOP_OVERLAP=keep table were uploaded to $ROOT/upload
for f in 00_scoring.R 00_select_visit.R 17_main_with_ukbb.R 38_fig6_numbers.R growth_refs_lms.csv supp_bars_results.csv analysis_table_chop_2026-11-19.rds; do
  [ -f "$IN/$f" ] || cp -p "$SRC/$f" "$IN/$f"
done
cp -p "$ROOT/upload/82_growth_blocks_v2.R" "$IN/82_growth_blocks_v2$TAG.R"   # each submission round keeps its own frozen copy
[ -f "$IN/analysis_table_chop_overlapkeep_2026-11-20.rds" ] || cp -p "$ROOT/upload/analysis_table_chop_overlapkeep_2026-11-20.rds" "$IN/"
chmod -R g-w,o-rwx "$ROOT"

submit() {   # name, table, visit rule, drop cohort, make blocks
  if [ -n "$ONLY" ] && [ "$1" != "$ONLY" ]; then return; fi
  local name="$1$TAG" table="$2" rule="$3" drop="$4" blocks="$5"
  local D=$ROOT/$name
  if [ -e "$D" ]; then echo "exists, not touched: $D"; return; fi
  mkdir "$D"
  cp -p "$IN"/00_scoring.R "$IN"/00_select_visit.R "$IN"/17_main_with_ukbb.R "$IN"/38_fig6_numbers.R "$IN"/growth_refs_lms.csv "$IN"/supp_bars_results.csv "$IN/$table" "$D/"
  cp -p "$IN/82_growth_blocks_v2$TAG.R" "$D/82_growth_blocks_v2.R"
  cat > "$D/run.sbatch" <<EOF
#!/bin/bash
#SBATCH -J v10_$name
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nice=100
#SBATCH -N 1 -n 1 --cpus-per-task=4 --mem=32G -t 02:00:00
#SBATCH -o $D/v10_$name.%j.log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=\$HOME/R/dnmt3a-r
export ANALYSIS_TABLE=$D/$table VISIT_RULE=$rule VARIANT=$name DROP_COHORT=$drop MAKE_BLOCKS=$blocks
cd $D
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
md5sum *.R *.rds *.csv > inputs.md5
\$RS 82_growth_blocks_v2.R && echo STEP82_OK
EOF
  (cd "$D" && sbatch run.sbatch)
}
submit main                        analysis_table_chop_2026-11-19.rds             fill_gaps ""   1
submit sens_visit_latest           analysis_table_chop_2026-11-19.rds             latest    ""   0
submit sens_no_chop                analysis_table_chop_2026-11-19.rds             fill_gaps CHOP 0
submit sens_chop_overlap_keep      analysis_table_chop_overlapkeep_2026-11-20.rds fill_gaps ""   0
submit sens_no_chop_g2mh203_kept   analysis_table_chop_overlapkeep_2026-11-20.rds fill_gaps CHOP 0
squeue -u jsebat -o "%.10i %.30j %.8T %.10M %R" | head -20
