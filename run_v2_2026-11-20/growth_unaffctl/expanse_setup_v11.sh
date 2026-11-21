#!/bin/bash
# v11_unaffected_controls_2026-11-20 (JS / Claude 2026-11-20): the growth package rerun with the corrected control set
# (00_control_set.R: SSC designated siblings + SSC parents + Simons VIP non-carriers without autism; SSC probands out).
# One run folder per job variant, each with its own frozen copies of the scripts and inputs. v9e and v10 are only read.
# Usage: bash expanse_setup_v11.sh [tag] [only-this-variant]
set -euo pipefail
TAG="${1:-}"; ONLY="${2:-}"
BASE=/expanse/projects/sebat1/jsebat/growth_tables_20260920/main_with_ukbb
ROOT=$BASE/v11_unaffected_controls_2026-11-20
IN=$ROOT/inputs
mkdir -p "$IN"
for f in 00_scoring.R 00_select_visit.R 17_main_with_ukbb.R 38_fig6_numbers.R growth_refs_lms.csv supp_bars_results.csv analysis_table_chop_2026-11-19.rds; do
  [ -f "$IN/$f" ] || cp -p "$BASE/v10_v2draft_2026-11-20/inputs/$f" "$IN/$f"
done
[ -f "$IN/43_supp_main_effects_by_band.R" ] || cp -p "$BASE/v9e_chop_genotype_bands/43_supp_main_effects_by_band.R" "$IN/"
cp -p "$ROOT/upload/82_growth_blocks_v2.R" "$IN/82_growth_blocks_v2$TAG.R"   # each submission round keeps its own frozen copy
for f in analysis_table_chop_unaffctl_2026-11-20.rds analysis_table_chop_unaffctl_overlapkeep_2026-11-20.rds; do
  [ -f "$IN/$f" ] || cp -p "$ROOT/upload/$f" "$IN/"
done
chmod -R g-w,o-rwx "$ROOT"

submit() {   # name, table, visit rule, drop cohort, make blocks, run 43
  if [ -n "$ONLY" ] && [ "$1" != "$ONLY" ]; then return; fi
  local name="$1$TAG" table="$2" rule="$3" drop="$4" blocks="$5" s43="$6"
  local D=$ROOT/$name
  if [ -e "$D" ]; then echo "exists, not touched: $D"; return; fi
  mkdir "$D"
  cp -p "$IN"/00_scoring.R "$IN"/00_select_visit.R "$IN"/17_main_with_ukbb.R "$IN"/38_fig6_numbers.R "$IN"/43_supp_main_effects_by_band.R "$IN"/growth_refs_lms.csv "$IN"/supp_bars_results.csv "$IN/$table" "$D/"
  cp -p "$IN/82_growth_blocks_v2$TAG.R" "$D/82_growth_blocks_v2.R"
  cat > "$D/run.sbatch" <<EOF
#!/bin/bash
#SBATCH -J v11_$name
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nice=100
#SBATCH -N 1 -n 1 --cpus-per-task=4 --mem=32G -t 02:00:00
#SBATCH -o $D/v11_$name.%j.log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=\$HOME/R/dnmt3a-r
export ANALYSIS_TABLE=$D/$table VISIT_RULE=$rule VARIANT=$name DROP_COHORT=$drop MAKE_BLOCKS=$blocks
cd $D
RS=/home/jsebat/micromamba/envs/dnmt3a-r/bin/Rscript
md5sum *.R *.rds *.csv > inputs.md5
\$RS 82_growth_blocks_v2.R && echo STEP82_OK
if [ "$s43" = "1" ]; then   # Supplementary Fig. 6 on the job's own five-group table written by 82
  export ANALYSIS_TABLE=$D/analysis_table_job.rds DROP_COHORT=
  \$RS 43_supp_main_effects_by_band.R && echo STEP43_OK
fi
EOF
  (cd "$D" && sbatch run.sbatch)
}
U=analysis_table_chop_unaffctl_2026-11-20.rds; K=analysis_table_chop_unaffctl_overlapkeep_2026-11-20.rds
submit main                        $U fill_gaps ""   1 1
submit sens_visit_latest           $U latest    ""   0 0
submit sens_no_chop                $U fill_gaps CHOP 0 0
submit sens_chop_overlap_keep      $K fill_gaps ""   0 0
submit sens_no_chop_g2mh203_kept   $K fill_gaps CHOP 0 0
submit sens_controls_ssc_probands  analysis_table_chop_2026-11-19.rds fill_gaps "" 1 1
squeue -u jsebat -o "%.10i %.34j %.8T %.10M %R" | head -20
