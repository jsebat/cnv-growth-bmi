# Compose the main figures of the BMI modifier paper, v2 draft (JS 2026-11-20), from the PNG blocks written by
# 82_growth_blocks_v2.R (growth_unaffctl/), 81_fig_longitudinal_v2.R (longitudinal/), 88_fig_puberty_mr_v2.R (puberty_mr/),
# 84-87 (autism/). Drawing only: no statistics are computed here.
# Inputs : run_v2_2026-11-20/<package>/*.png ; Fig. 5 top = the sex-stratified biobank figure of the v1 draft (fig_sex_biobank_v1.png,
#          exported from the manuscript, panels A-D).
# Outputs: run_v2_2026-11-20/figures/Fig1..Fig5 (PNG) and copies in the manuscript Drive folder (argument 1), if given.
# null model: none; no tests.
import os, sys, shutil
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); RUN = os.path.join(HERE, "run_v2_2026-11-20"); OUT = os.path.join(RUN, "figures")
os.makedirs(OUT, exist_ok=True)
W = 4000; GAP = 60

def load(rel):
    im = Image.open(os.path.join(RUN, rel)).convert("RGB")
    if im.size[0] != W: im = im.resize((W, round(im.size[1] * W / im.size[0])), Image.LANCZOS)
    return im

def stack(name, parts):
    ims = [load(p) for p in parts]
    H = sum(i.size[1] for i in ims) + GAP * (len(ims) - 1)
    fig = Image.new("RGB", (W, H), "white"); y = 0
    for i in ims: fig.paste(i, (0, y)); y += i.size[1] + GAP
    f = os.path.join(OUT, name + ".png"); fig.save(f, dpi=(200, 200)); print(name, fig.size); return f

figs = {
    "Fig1_height": ["growth_unaffctl/block_height.png", "longitudinal/block_height_longitudinal.png", "puberty_mr/block_puberty_mr.png"],
    "Fig2_head_autistic_traits": ["growth_unaffctl/block_hc.png", "longitudinal/block_hc_longitudinal.png", "autism/block_srs_scq.png"],
    "Fig3_autism_diagnosis": ["autism/fig_asd_by_cnv_v2.png"],
    "Fig4_cnv_bmi_pgs": ["autism/fig_cnv_bmi_pgs_main.png"],
    "Fig5_sex": ["figures/fig_sex_biobank_v1.png", "autism/block_sex_autism.png"],
    "SuppFig_cnv_bmi_pgs_full": ["autism/fig_supp_cnv_bmi_pgs_full.png"],
    "SuppFig_main_effects_by_band": ["growth_unaffctl/fig_supp_main_effects_by_band.png"],
    "SuppFig_threeway_forest": ["autism/fig_supp_threeway_forest.png"],
    "Fig6_multimodal": ["multimodal/fig_multimodal_main.png"],
    "SuppFig_multimodal": ["multimodal/fig_supp_multimodal.png"],
}
made = []
for k, v in figs.items():
    if all(os.path.exists(os.path.join(RUN, p)) for p in v): made.append(stack(k, v))
    else: print("SKIPPED", k, [p for p in v if not os.path.exists(os.path.join(RUN, p))])
if len(sys.argv) > 1:
    for f in made: shutil.copy2(f, os.path.join(sys.argv[1], "v2_" + os.path.basename(f))); print("copied", os.path.basename(f))
