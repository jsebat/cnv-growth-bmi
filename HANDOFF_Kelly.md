# Handoff for Kelly Lee (and her Claude): replicating the BMI modifier paper

Written 2026-11-20 by Jonathan Sebat's Claude session. Companion file: `HANDOFF_Kelly.md` in
`github.com/jsebat/cnv-dosage-architecture` (CNV main effects and 22q11.2 subregions for the other paper).

**Goal.** Rerun the non-genetic analyses of the paper on your cohort (CNV carriers and age-matched non-carrier controls)
and send back summary tables. Your results will be combined with ours by meta-analysis and added to Figures 1, 2, 3 and
the growth and autistic-trait parts of Figure 5. Figure 4 needs polygenic scores and is not part of this.

Nothing individual-level leaves your site. The script writes aggregate CSV files only.

## What the paper tests

For carriers of 16p11.2 BP4-BP5 and 22q11.2 A-D deletions and duplications:

1. Does the relationship between BMI and **height** differ between carriers and controls (CNV x BMI), by age band?
2. The same for **head circumference**.
3. In children with repeated visits: does BMI at one visit predict the **change** in height or head circumference by the next?
4. Do height, head circumference and BMI differ between carriers **with and without an autism diagnosis**, and does the
   proportion with a diagnosis change with BMI or head size?
5. If you have SRS or SCQ scores in carriers and controls: do **autistic traits** rise with BMI more in carriers?
6. Do any of these differ by **sex** (exploratory)?

Each is run per CNV, then pooled across the four CNVs with an indicator for CNV, with a test of whether the CNVs differ.

## Steps

1. Read `federated/README.md` (input table, exact commands, output files).
2. Check your R setup on synthetic data: `Rscript federated/tests/test_federated.R` (36 checks should pass).
3. Make one long table of visits: id, family id if you have one, sex, age at visit, height (cm), weight (kg), head
   circumference (cm), genotype group, and, if available, autism diagnosis (0/1), SRS total, SCQ total. Copy
   `federated/config_example.json` and set your column names and genotype labels.
4. Run `federated/site_export.R` with your config. It scores, picks one visit per person for the cross-sectional
   models, cleans repeated visits for the longitudinal models, fits everything and writes the output folder.
5. Send back the whole output folder, plus a few lines: how carriers and controls were recruited, how CNVs were
   genotyped, whether controls include anyone with a neurodevelopmental diagnosis, and how autism diagnosis was recorded.

## How trajectories are shared without individual data

The growth-by-age panels are drawn from two aggregates that your run exports for each CNV and BMI tail: mean Z in fixed
age bins (with n), and the coefficients of a curve fitted on a spline basis that is identical at every site.
`federated/meta_combine.R` pools the curves across sites and `federated/plot_federated.R` draws them. Small bins are
suppressed.

## Fixed choices (please keep them)

- **Scoring.** Children: WHO height and BMI, UK 1990 head circumference. Adults: Box-Cox with fixed exponents, relative
  to the non-carrier adults of your own cohort. A site with fewer than 100 non-carrier adults per sex gets no adult
  bands; that is expected for a paediatric cohort.
- **One visit per person** in cross-sectional models, chosen by the gap-filling rule (`00_select_visit.R`): people with
  visits in several age bands are placed in the band with the fewest people of their genotype, and one visit in that
  band is drawn at random with a fixed seed. Option `"latest"` gives the most recent visit instead.
- **Age bands** 2-9, 10-18, 19-39, 40-73 years; longitudinal windows 2-8, 8-13, 13-18 (height) and under 2, 2-5, 5-18
  (head circumference).
- **Minimum cell size** 5 by default (`min_cell` in the config; raise it if your site requires).
- **Controls** are non-carriers. On 2026-11-20 we found that the child controls of our own growth analysis had been
  autistic SSC probands and corrected it (`HANDOFF_2026-11-20_v2.md`). Please tell us exactly who your controls are.

## What we found, so you know what to expect

- At 2-9 years the BMI effect on growth differs by CNV: heavier 16p11.2 duplication carriers are shorter (CNV x BMI
  -0.41 SD per SD); head circumference tracks BMI more strongly in 16p11.2 deletion and 22q11.2 duplication carriers.
  The 22q11.2 duplication result rests on one hospital's records, so an independent estimate from you is valuable.
- Higher BMI precedes faster height gain before age 8, not after (pooled +0.089 SD per SD).
- No pooled all-CNV effect; the CNVs differ in direction.
- Sex differences are not robust in our data.

Full results: `run_v2_2026-11-20/` (one README per package). Do not tune anything toward these numbers; a clean
non-replication is as useful to us as a replication.

## If something does not fit

Change the config or your input table, not the models. If a model has to change for your data, write down what and why
and send that with the results. If your Claude finds a bug, open an issue or pull request here so every site gets the
fix. Questions go to Jon (jsebat).
