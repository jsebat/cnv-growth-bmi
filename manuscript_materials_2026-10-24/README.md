# Manuscript materials: body mass, growth and social-communication / cognitive phenotypes in CNV carriers

Snapshot of 2026-10-24, after plausibility QC of height / weight / head circumference and the SSC visit-age correction.
All numbers are provisional: the next round adds subjects and manually corrected records (see HANDOFF_2026-10-24.md).
Source scripts and live outputs are one level up (CNV_growth_trajectories); copies here are frozen.

- main_figure/            fig_main_v3.{png,pdf}, fig_main_v3_legend.txt, main_v3_bmi_slope_heterogeneity.csv
                          (Height, HC, SRS, SCQ x unaffected controls / 16p11.2 del / 16p11.2 dup / 22q11.2 del / 22q11.2 dup; script 25)
- supplementary_figures/  fig_supp_prospective            longitudinal: BMI at a visit -> change in height / HC Z (script 11)
                          fig_cognition_bmi_slopes        SRS/SCQ vs BMI with SSC ASD probands and a forest plot of slopes (script 24)
                          fig_supp_cnb_accuracy_total_bmi CNB accuracy total by CNV group (+ legend with heterogeneity test; script 26)
                          fig_supp_cnb_accuracy_bmi / fig_supp_cnb_speed_bmi   all CNB subtests (+ legend; script 26)
                          fig_supp_asd_by_cnv             height / HC / BMI by autism diagnosis within CNV group, SSC probands vs siblings (+ legend with PGS sensitivity; script 32)
- supplementary_tables/   supp_tables_bmi_effects.xlsx (sheets S1 growth, S2 cognition, S3 prospective, S4 sex and PGS strata, S5 PGS adjusted, S6 interaction power, S7 ASD by CNV, S8 ASD by CNV PGS adj, raw, notes) and .md;
                          supp_table_S1-S8 CSVs (scripts 28, 29, 30, 32); medication_sensitivity_results.md (methods, results, table for height, HC, SRS, SCQ, CNB; scripts 27 + medication_table.R)
- data_qc/                implausible_height_weight_corrections.csv (fill action / corrected_* to add records back),
                          implausible_height_weight.{csv,md} (flagged records), 00_qc_growth.R (rules)
- Analysis_report_2026-10-26.docx / .pdf   final report: summary of findings, methods, results, embedded figures and legends, table index (S1-S8, medication, QC), sample sizes, Dropbox and Expanse paths, PMID references (31_make_report.R, then data_qc/finish_report.ps1 for heading numbers, table layout and the PDF export via Word)
- HANDOFF_2026-10-26.md   current handoff (incl. the ASD-by-CNV supplementary figure and Tables S7-S8 added later that day) (pipeline, run order, caveats, Expanse resources, next steps); HANDOFF_2026-10-24.md kept for history
- reports/                HANDOFF_2026-10-24.md (pipeline, run order, caveats, next steps), unaffected_controls_report.md (full analysis report)
- figure_legends.md        all figure legends in one file (main figure, every supplementary figure, table notes); each figure folder also holds its own *_legend.txt
