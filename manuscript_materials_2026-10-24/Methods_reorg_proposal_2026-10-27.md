# Methods and Supplementary Information: proposed reorganisation and new Figure 4 methods

**Applied to the Google Doc on 2026-10-27** (batches A-C via the Docs connector; 83 Paperpile citation links and all 21 comment threads verified intact). Table numbering resolved: biobank Tables 1-9, GO:0008284 themes = Table 10, growth tables = 11-20 (merged workbook `Sacks Supplementary Tables v2.xlsx`).

Prepared 2026-10-27 from the Google Doc "Context-Dependent Effects of CNVs on Growth" (read in full, 21 comment threads) for
application in place once the Google Docs connector is enabled. Scope: reorder only, add the missing Figure 4 methods, fix
cross-references and formatting. No rewriting of existing text. Figure 5 (22q11.2 dissection) is left as it is.

Citation convention: `[n](https://paperpile.com/c/BnyKao6h/KEY)` reuses a reference already in the Paperpile library
(the key is taken from the Doc). `[PMID:nnnnnnnn, to add]` marks a reference that is not yet in the library and has to be
inserted through the Paperpile sidebar (search by PMID); I cannot add items to the library from here.

---

## 1. Proposed order of the main Methods

The order follows the figures. "Existing" = text already in the Doc, moved as a block; "new" = text supplied below;
"missing" = no text exists yet and a placeholder paragraph is inserted for the authors.

| # | Section | Figures | Status |
|---|---|---|---|
| 1 | Cohorts and ethics (one paragraph naming the five biobanks and the clinical cohorts; the EstBB ethics statement currently sitting in the Acknowledgements moves here) | Fig. 1 | new, short; cohort-level statements [to confirm] |
| 2 | Genotyping and CNV calling (CNV genotype assignment; UKBB, All of Us, EstBB, MyCode, MVP, G2MH, SSC, Simons Searchlight) | Fig. 1 | existing |
| 3 | SNP imputation and QC | | existing |
| 4 | Ancestry PCA | | existing |
| 5 | Polygenic scores (biobank PGS-CS scores) | Fig. 1, Ext. Data 5 | existing; the SBayesRC scores used for the clinical cohorts are described in the Supplementary Methods and cross-referenced here |
| 6 | Scaling of height and BMI (adults: residual on sex, age, age^2, Box-Cox within sex; children: childsds, WHO height, CDC BMI) | all | existing |
| 7 | CNV main effects, fixed-effects meta-analysis and variance explained | Fig. 1A-D | existing |
| 8 | Dose-response curves and model comparison (linear, quadratic, one- and two-component fits, LOO log-likelihood) | Fig. 1E-F, Ext. Data 1-2 | missing: placeholder |
| 9 | Classification of mirror and asymmetric CNVs (perpendicular offset from the mirror line, measurement-error-aware one- and two-component mixture, parametric bootstrap, size and gene-count comparisons, Theil-Sen fits) | Fig. 3A-H | missing: placeholder |
| 10 | Gene-level simulations of dose-response models (mirror, continuous generalized-dysfunction, sparse single-hit) | Fig. 3I-M | existing "Simulation framework" text, currently in Supplementary Materials; a two-sentence summary goes here and the full text moves to Supplementary Methods |
| 11 | Pathway analysis | Fig. 4 | new (section 2 below) |
| 12 | Dissection of the 22q11.2 A-D locus (additivity of subregion effects, observed versus expected) | Fig. 5, Ext. Data 7 | missing: placeholder; left for later as requested |
| 13 | Body mass as a modifier of growth and autistic traits (Participants; Anthropometric measures; Statistical analysis; Sensitivity analyses) | Fig. 6 | existing; formatting and cross-reference fixes (section 3 below) |
| 14 | Age-dependent effects and mediation analysis (cross-sectional trajectories, age-stratified mediation with the `mediation` package, ref. 76) | Ext. Data 4-6 | missing: placeholder |
| 15 | Sex-stratified analysis (slope of female versus male effects, F-test of the joint hypothesis slope = 1 and intercept = 0; per-CNV sex x CNV interaction with the >= 200-carrier rule; exclusion of sex-imbalanced cohorts) | sex figure (currently also labelled "Figure 5") | missing: placeholder |
| 16 | Statistics and reproducibility; Reporting summary | | existing (one paragraph) |
| 17 | Data availability | | existing, merged into one statement for all cohorts |
| 18 | Code availability | | existing, merged: https://github.com/mollysacks/Sacks2026 (biobank meta-analysis, Figures 1-3, 5 and the sex analysis) and https://github.com/jsebat/cnv-growth-bmi (part 1, Figure 4 pathway analysis; part 2, Figure 6 and its supplement) |

Items to decide before applying: (a) the figure numbers. The Doc has "Fig. 0" (study design), "Figure 1" (main effects),
"Fig. 2" cited in the text but the legend reads "Figure 3" (mirror/asymmetric), "Figure 4" (pathways), "Figure 5"
(22q11.2), "Figure 6" (BMI modifier) and a second "Figure 5" (sex). The reorganised Methods assume the order
1 main effects (with the design schematic folded in, as your comment proposes), 2 mirror versus asymmetric, 3 pathways,
4 22q11.2, 5 BMI modifier, 6 sex. Tell me if you prefer to keep the current numbers for now. (b) Supplementary Table
numbers. The biobank part already uses Supplementary Tables 1-4 and 9, and the growth part also uses 1-10; the growth
tables need renumbering (for example 10-19) or lettering, and every cross-reference in the Figure 6 text, legend and
supplement changes with them.

---

## 2. New Methods text: pathway analysis (Figure 4)

**Pathway analysis of asymmetric CNVs and single-gene height associations.** Protein-coding genes within the 34 recurrent
CNV intervals with reciprocal height effects (13 asymmetric and 21 mirror intervals; the nested 22q11.2 intervals A-B, A-D,
B-C, B-D and C-D were entered separately) were assigned the dose-response class of their interval (708 genes, GRCh38;
13 genes in nested 22q11.2 intervals belong to both classes). Gene sets were the GO Biological Process 2023 and KEGG 2021
Human libraries as distributed by Enrichr [PMID:27141961, to add] [PMID:36866529, to add] [PMID:33125081, to add],
restricted to sets of 5 to 500 genes (5,712 sets). Three tests were applied to every set. (1) Asymmetric versus mirror
CNVs: the set's genes within the CNV gene universe (215 genes mapped to the libraries: 124 in asymmetric intervals,
104 in mirror intervals) were cross-classified as asymmetric-class or not, giving an odds ratio (OR > 1, over-represented
in asymmetric CNVs; OR < 1, in mirror CNVs) with a two-sided Fisher exact *P*. Because the genes of one interval are
not independent, the statistic of record for this layer is a locus-block permutation: the asymmetric label was
reassigned among the 34 intervals, within tertiles of interval gene count and preserving the observed number of
asymmetric intervals per tertile, and the difference between the fraction of asymmetric-class genes and the fraction of
mirror-class genes belonging to the set was recomputed in each of 5,000 permutations (*P* = (number of permutations
at or above the observed difference + 1)/(5,001)). (2) Enrichment of the 207 genes associated with height in the exome
study of Kosmicki et al. [25](https://paperpile.com/c/BnyKao6h/Y1R6) (176 mapped to the libraries), one-sided
hypergeometric test against the library's gene universe. (3) Directional bias: among the Kosmicki genes in the set,
enrichment of height-decreasing associations relative to the fraction of height-decreasing genes among all mapped
Kosmicki genes (65%), one-sided hypergeometric test. The analysis set for all three tests was the 102 gene sets with at
least four genes in the CNV universe, the minimum at which the CNV test can reach *P* < 0.05 (the asymmetric test needs
at least six, and those 48 sets are a subset); Bonferroni correction (0.05/102) and Benjamini-Hochberg false-discovery
rates were computed within this set for each test. A set was called convergent when it was nominally significant in all
three tests; the same procedure run with mirror in place of asymmetric CNVs identifies no convergent set. Panels D and E
list the asymmetric-CNV genes of the convergent set by locus and its Kosmicki genes by effect size. Analyses were run in
R 4.6.1 with data.table; code, inputs and the full results table are at https://github.com/jsebat/cnv-growth-bmi
(`fig4_pathway_analysis/`).

Notes for the authors. (i) The CNV-class table used here is the "uncorrected 22q11.2 classification" of the July 2026
pipeline (A-B and A-D asymmetric; B-C, B-D and C-D mirror); please confirm it matches the final Figure 3 classes
(Figure 3 reports 12 asymmetric of 33 loci; the pathway table has 13 of 34 because A-B and A-D are both present).
(ii) The figure legend states the two-sided *P* for GO:0008284 as 0.047 and the one-sided pipeline value as 0.031; the
Methods above describe both. (iii) The three PMIDs are the Enrichr 2016 update (Kuleshov et al.), the Gene Ontology
knowledgebase 2023 paper and the KEGG 2021 update (Kanehisa et al.); check them when adding to Paperpile.

---

## 3. Fixes to the existing Figure 6 ("Effects of BMI on growth and autistic traits") Methods block

- Formatting: the five subsections are pasted as Heading 2 paragraphs (every line starts with `##`). They become body
  text with bold run-in headings (Participants, Anthropometric measures, Statistical analysis, Sensitivity analyses),
  matching the rest of the Methods.
- The superscript citations 5-12 in this block (SSC^5, Simons VIP^6, WHO^7,8, UK90^9, SRS^10, SCQ^11, mgcv^12, CNB^13,14)
  are numbers from the stand-alone Nature Genetics draft, not Paperpile citations. Simons VIP already exists in the
  library ([48](https://paperpile.com/c/BnyKao6h/hn9AD)). The others are new and need adding by PMID: SSC design
  (Fischbach & Lord 2010, PMID:20955926), WHO child growth standards (de Onis 2006, PMID:16817681) and WHO 2007
  reference (de Onis 2007, PMID:18026621), UK 1990 head-circumference reference (Freeman/Cole 1995, PMID:7639544;
  Cole 1998 LMS, PMID:9496720), SRS (Constantino 2003, PMID:12959421), SCQ (Berument 1999, PMID:10789276), mgcv
  (Wood 2011, J. R. Stat. Soc. B, no PMID; DOI 10.1111/j.1467-9868.2010.00749.x), Penn CNB (Gur 2010, PMID:19914906;
  Moore 2015, PMID:25180981). [to confirm] as in the earlier draft.
- Cross-references: "the 70th and 30th centile split in Fig. 1" becomes Fig. 6; "Supplementary Fig. X" (prospective)
  and "Supplementary Fig. Y" (autism within CNV) get numbers once the supplementary figure list is fixed (proposed:
  Supplementary Fig. 1 prospective growth, 2 CNB accuracy total, 3 CNB subtests, 4 autism diagnosis within CNV group;
  the "Supplementary Fig. 5" cited in the supplementary results is the same figure as 4 and is corrected).
- Participants: the sentence "[to confirm: CNV calling and inclusion criteria per cohort; ethics approvals ...]" can now
  point to the Genotyping and CNV calling section (G2MH, SSC, Simons Searchlight are described there); Cardiff is not
  and needs a sentence [to confirm].
- Data availability: merge with the biobank statement (section 17 above).

---

## 4. Proposed order of the Supplementary Information

| # | Section | Status |
|---|---|---|
| 1 | Supplementary Note 1: Model misspecification (cited in Methods, scaling) | missing from the Doc; placeholder |
| 2 | Supplementary Note 2: Statistical power (cited in the main text and the sex section) | missing from the Doc; placeholder |
| 3 | Supplementary Note 3: Comparison with Milind et al. | existing (the "Comparison of our results with Milind et al." text and its two numbered points); the trailing sketch lines ("A-B DEL suppresses ...", "Candidate genes for") stay as author notes |
| 4 | Supplementary Methods, part A: gene-level simulations (Simulation framework; Mirror-gene model; Continuous generalized-dysfunction model; Discrete single-hit model) | existing, moved from "Supplementary Materials" |
| 5 | Supplementary Methods, part B: body mass as a modifier (Ages of SSC relatives; Plausibility screen; Longitudinal models; Autism diagnosis within CNV group; Polygenic scores; Sex and polygenic-score sensitivity analyses; Medication; Cognitive measures) | new to the Doc: the main Methods already cite "Supplementary Methods" for these, but the text was never pasted; it is in `NatGenet_sections_2026-10-26.md` and is copied in verbatim |
| 6 | Supplementary Results, part A: simulations (Results for height; Results for BMI) | existing, moved |
| 7 | Supplementary Results, part B: body mass as a modifier (Longitudinal analysis; Head circumference by autism diagnosis; Cognitive performance; Sex and polygenic background; Psychotropic medication) | existing; three corrections: "seven were nominally significant" becomes "eleven were nominally significant (about two expected by chance, although the tests are correlated) and none survived correction for 47 tests (smallest P = 0.017)"; "Supplementary Fig. 5" becomes the autism-within-CNV figure number; the pre-screen CNB sentence ("had been nominally significant ... before the plausibility screen") is removed, as agreed for the supplement |
| 8 | Extended Data Figures 1-8 in numerical order (3 and 8 are cited but have no legend in the Doc; 8 = CNV length and gene count versus effect size) | existing; gaps flagged |
| 9 | Supplementary Figures 1-4 (growth part) | existing legends for 1 and 4; the two CNB figures (2, 3) need their legends pasted from `NatGenet_sections_2026-10-26.md` |
| 10 | Supplementary Tables index | new one-line index per table after the renumbering decision |
| 11 | The link to the simulation reports folder on Drive | existing, kept at the end of Supplementary Methods part A |

The "Supplementary Table SX. Biological themes of GO:0008284" table that now sits in the main text after the Figure 4
legend moves to the Supplementary Tables (it becomes the Figure 4 supplementary table) and the main text cites it by
number.

---

## 5. What stays untouched

Abstract, Main text results sections, Discussion, Acknowledgements, Figure legends 1-6, the Figure 5 (22q11.2) section and
its replication notes, and all reference entries. Open comment threads are left open.
