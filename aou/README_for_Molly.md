# All of Us adult anchor: height by BMI in 16p11.2 / 22q11.2 carriers

Molly, thank you for running this in the All of Us Researcher Workbench. It's the adult counterpart to the clinical
growth-trajectory figure (fig_main_growth): for adult carriers of the 16p11.2 BP4-BP5 and 22q11.2 A-D CNVs it gives mean
height Z by BMI tail and age band, and tests whether the BMI-height relation in carriers differs from non-carriers.

**Script:** `AoU_adult_anchor.R` (base R only; nothing to install).

**Steps:**
1. Copy `AoU_adult_anchor.R` into your Workbench workspace next to your Height/BMI follow-up scripts.
2. Edit the two SETTINGS lines: `ANCESTRY` (EUR, AAM or LATNAT) and `GENOMATRIX`, the path to that ancestry's genomatrix
   (the one used for the main-effect analyses, e.g. `genomatrix_CNV_EUR.20250715.csv`).
3. `source("AoU_adult_anchor.R")`, once per ancestry.
4. Check the two CSVs in `adult_anchor_output/`. Every group under 20 people is already suppressed (shown as "<20"),
   per the AoU dissemination policy, and only aggregate statistics are written.
5. Export them as usual and drop them in this folder (`.../CNV_growth_trajectories/aou/`).

**Columns it uses:** `person_id`, `sex_at_birth`, `Height_raw`, `Height_age`, `BMI_raw`, `BMI_age`,
`X16p11.2_proximal`, `X22q11.2_BP1_BP4_A_D` (1 = deletion, 2 = reference, 3 = duplication), as in
`Height_FollowUpAnalyses_01212025.AoU.EUR.R`. If any are named differently in your current genomatrix, edit the `cols`
line.

**Method (same as the UK Biobank run, `ukbb_anchor_reference.R`):**
- Height Z and BMI Z are residuals on age + age² within sex, scaled by the non-carrier residual SD.
- BMI tails are defined at the 70th and 30th centiles of BMI Z.
- Tests: `hZ ~ bmiZ * CNV + sex + age`, giving the CNV x BMI term, and the same with the high/low BMI tails.

UK Biobank results for comparison (EUR, 452,449 non-carriers):

| Locus | Carriers | CNV effect on height | CNV x BMI |
|---|---|---|---|
| 16p11.2 del | 117 | -1.10 SD | +0.04 (p = 0.41) |
| 16p11.2 dup | 149 | +0.24 SD | -0.27 (p = 0.009) |
| 22q11.2 dup | 270 | -0.26 SD | +0.04 (p = 0.36) |

Questions to Jon.
