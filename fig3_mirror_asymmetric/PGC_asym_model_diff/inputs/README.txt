Copies of ../behavior_figure_revised/inputs/{height,bmi}_loci.csv with 22q11.21_A-D added
(the original package excluded it as "overlapping" while keeping other overlapping loci
such as 22q11.21_B-D, 16p11.2_BP1-BP3, 1q21.1_TAR_BP1-BP3, 16p13.11_BP1-BP3).

22q11.21_A-D row (added 2026-09-29):
- DEL/DUP beta and SE (and BMI carrier counts): MetaAnalysis/<trait>.v2/fixed_effects/2_CNV_main_effect.csv
- Length: chr22:19,049,819-21,112,437 (Supplementary Table T1) = 2.062618 Mb
- Gene_count = 109: union of the 22q11.21_A_B (58) and 22q11.21_B-D (43) gene lists in
  ../bmi/inputs/locus_gene_membership.csv (no overlap, 101 genes), plus 8 Ensembl genes
  overlapping the LCR22-B gap between them (chr22:20,250,000-20,400,000) that are not
  already in those lists: DGCR6L, PRODHLP, FAM230G, USP41P and 4 unnamed lncRNAs. Some of
  the unnamed lncRNAs may correspond to clone-named genes in the older lists, so the true
  count on the lists' annotation is 105-109.
- Height Dose_response_class = "unclassified" (not used by the mirror_only scripts).
