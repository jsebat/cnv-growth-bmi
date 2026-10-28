# Figure 3: mirror and asymmetric CNV effects on height (oligogenic model)

Code and summary data for Figure 3 of Sacks et al. (2026), "Mirror effects of CNVs have an oligogenic basis. Asymmetric
effects are consistent with a rare single hit acting on the duplication axis." Jonathan Sebat (UCSD), analyses carried out
with Claude, July to October 2026. Copied from the lab Dropbox
(`Manuscripts/PGC_joint_submission/Molly/ContinuousVsSingleHit/`) on 2026-10-28.

The manuscript figure is `perp_offset_bimodality/figures/Figure_height_main_v2_me.{pdf,png}`; its legend and build notes are
`perp_offset_bimodality/Figure_height_main_v2_LEGEND.md`.

## Layout

The two folders keep their original names and relative positions, because the scripts find each other through `../`.

- `CNV_height_BMI_corrected_directional_models_v2/`: the dose-response model project. It holds the per-locus deletion and
  duplication effects on height and BMI from the Sacks et al. meta-analysis (`height/inputs/`, `bmi/inputs/`). It fits
  the mirror-gene model and the two directional models, continuous generalized dysfunction and sparse single hit
  (`code/run_corrected_analysis.py` and `code/run_axis_mixture_models.py`, with an R port in `code/*.R`). It also has the
  Macé et al. 2017 burden reference and the earlier report and figures.
- `perp_offset_bimodality/`: the Figure 3 analyses. These are the perpendicular offset from the mirror line and its
  bimodality (panels A and B), the size, gene-count and effect comparisons of mirror and asymmetric loci (C to H), and
  the simulations that score the three models (I to K). `README.md` there describes the bimodality test in full.
- `README_ContinuousVsSingleHit.md`: the original overview of both folders.

## Reproduce Figure 3

Run from `perp_offset_bimodality/` with Python 3 (numpy, pandas, scipy, matplotlib):

```bash
python ../CNV_height_BMI_corrected_directional_models_v2/code/run_axis_mixture_models.py   # model fits, if they have changed
python code/perp_offset_bimodality.py --nboot 5000 --starts 25                         # offsets, mirror-class calibration
python code/me_bimodality.py                                                           # measurement-error-aware mixture -> outputs/me/all_loci/
python code/score_models_on_criteria.py --me --out-dir outputs/me/all_loci             # simulated datasets per model
python code/make_height_figure_v2.py --bimodality-dir me/all_loci --criteria-dir me/all_loci --name Figure_height_main_v2_me
```

`make_supp_model_clouds.py` draws the supplementary model-cloud figure. `make_supp_table_sensitivity.py` writes the
sensitivity table that excludes loci with no nominal deletion or duplication effect (`tables/`).

## Not included

- `Sacks_supplementary_tables.xlsx`, the full manuscript supplement that `run_corrected_analysis.py` reads to build the
  per-locus input CSVs. The CSVs it produces are included, so every later step runs.
- Word copies of the reports (the PDFs are kept) and Python caches.

All inputs are locus-level summary statistics; there are no participant-level data.
