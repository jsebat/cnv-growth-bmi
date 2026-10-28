# Condensed legend — Figure_height_main_v2_clouds_me

Figure file: `perp_offset_bimodality/figures/Figure_height_main_v2_clouds_me.pdf`

---

**Figure 3. Asymmetric gene-dosage effects of CNVs on height are consistent with a rare
single hit acting on the duplication axis.** Analyses include 33 reciprocal CNV loci with
published height effects (21 mirror, blue; 12 asymmetric, red). (**A**) Perpendicular
offsets from the mirror line (DUP = −0.351 × DEL) follow a bimodal distribution. Curves
show measurement-error-aware one- and two-component mixture fits; rug marks indicate
published classes. (**B**) Parametric bootstrap test of bimodality (5,000 replicates;
observed likelihood ratio = 9.9, *p* = 0.016). (**C–F**) Asymmetric CNVs are longer
(*p* = 0.006), contain more genes (*p* = 0.019), have comparable deletion effects
(*p* = 0.651) but more negative duplication effects (*p* = 10⁻⁶). Bars indicate medians.
(**G, H**) After stratifying mirror CNVs from asymmetric CNVs, gene count correlates with
effect size negatively for deletion and positively for duplication (Spearman ρ = −0.61
and +0.51, respectively), but not for asymmetric loci. Lines show Theil–Sen fits,
including pooled estimates (black dotted). (**I, J**) Effects simulated under each
directional model, pooled over 150 datasets generated with the observed gene counts and
standard errors, with observed loci overlaid in black. The continuous model displaces
every locus in proportion to gene count (ΔAIC 2.6), filling the gap between classes; the
single-hit model displaces only loci carrying a hit (36% at the median gene count;
ΔAIC 0.0), generating two separated clouds. (**K–M**) Comparison of mirror-only,
continuous and single-hit models. (**K**) Both directional models reproduce the
association of duplication burden with short stature observed in previous studies
(Macé et al. 2017), unlike mirror-only. (**L**) Only the single-hit model predicts two
distinct clusters of mirror and asymmetric CNVs. (**M**) Simulated likelihood ratios show
that single-hit, but not continuous or mirror-only, typically reproduces the observed
bimodality (median *p* = 0.014, 0.308 and 0.423, respectively). Effects are in SD units
unless otherwise indicated.

---

Built by:
```bash
python code/make_height_figure_v2.py --with-clouds \
    --bimodality-dir me/all_loci --criteria-dir me/all_loci \
    --name Figure_height_main_v2_clouds_me
```
