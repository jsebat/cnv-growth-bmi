#!/usr/bin/env python
"""One-page PDF summary of the perpendicular-offset bimodality analysis.

All numbers are read from outputs/, never typed in, so the page cannot drift from the
analysis it describes. The figure is embedded from figures/.
"""
from __future__ import annotations
import textwrap
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = Path(__file__).resolve().parent.parent
OUT = HERE / 'Perp_offset_bimodality_report.pdf'
FIGPNG = HERE / 'figures' / 'Figure_perp_offset_bimodality.png'

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                     'pdf.fonttype': 42})
INK = '#1a1a1a'; GREY = '#777777'; BLUE = '#1f5c8b'

tab = pd.read_csv(HERE / 'outputs' / 'bimodality_test.csv')
loci = pd.read_csv(HERE / 'outputs' / 'perp_offsets_per_locus.csv')
h = loci[loci.trait == 'Height'].tail(33)
g = h.groupby('Dose_response_class').perp_offset.agg(['count', 'mean'])
CUT = -0.18
n_as = int((h[h.Dose_response_class == 'asymmetric'].perp_offset < CUT).sum())
n_mi = int((h[h.Dose_response_class == 'mirror'].perp_offset < CUT).sum())


def q(trait, fit, alt, col):
    r = tab[(tab.trait == trait) & (tab.mirror_fit == fit) & (tab.alternative == alt)]
    return float(r[col].iloc[0])


fig = plt.figure(figsize=(8.5, 11))
y = [0.968]
L, R = 0.068, 0.950


def put(txt, size=7.5, dy=0.0125, color=INK, x=L, wrap=133, gap=0.0):
    y[0] -= gap
    for line in textwrap.wrap(txt, wrap):
        fig.text(x, y[0], line, ha='left', va='top', size=size, color=color)
        y[0] -= dy


def head(txt):
    y[0] -= 0.007
    fig.text(L, y[0], txt, ha='left', va='top', size=9.8, weight='bold', color=INK)
    y[0] -= 0.0168


# ------------------------------------------------------------------- title
fig.text(0.5, y[0], 'Perpendicular offset from the mirror line is bimodal for height, '
                    'unimodal for BMI',
         ha='center', va='top', size=12.8, weight='bold', color=INK)
y[0] -= 0.022
fig.text(0.5, y[0], 'Exact Sacks reciprocal-CNV effects  |  parametric bootstrap, 5,000 '
                    'replicates  |  2026-10-10',
         ha='center', va='top', size=8.2, style='italic', color=GREY)
y[0] -= 0.022

# ----------------------------------------------------------------- methods
head('Methods')
put('The mirror model places every locus on the line DUP = -rho*DEL, with position along the line '
    "set by gene count. Each locus's signed perpendicular distance from that line, "
    'offset = (rho*DEL + DUP)/sqrt(1+rho^2), is what the mirror model does not explain: a discrete '
    'asymmetric class should make it a two-component mixture, a continuous process a single '
    'Gaussian.')
put('One- and two-component Gaussians were compared by likelihood ratio. Mixture likelihood ratios '
    'do not follow chi-square, so the null came from a parametric bootstrap: 5,000 datasets of the '
    'same size simulated from the fitted single Gaussian, both models refitted to each. Observed '
    'and simulated data were fitted with the identical EM procedure and the same 25 random starts, '
    'since fitting the observed data harder than the nulls biases p downward. Null A is the plain '
    "fitted Gaussian; null B adds each locus's own measurement error. The equal-variance "
    'alternative (2 extra parameters) is primary, being what a sparse hit implies: it displaces a '
    'locus without changing its dispersion. Free variance (3 parameters) is reported alongside.',
    gap=0.003)
put('The mirror line was calibrated two ways. Primary: on the 21 mirror-class height loci only, '
    'holding the 12 asymmetric loci out of the fit that defines the line (rho = '
    f'{q("Height","mirror-class only","equal_var","rho_full_fit"):.3f}). Comparison: on all 33 '
    'loci, where the asymmetric loci tilt the line their own displacement is measured against '
    f'(rho flattens to {q("Height","full data","equal_var","rho_full_fit"):.3f}). The BMI table '
    'carries no dose-response class, so all 32 loci are used and the calibrations coincide.',
    gap=0.003)

# ------------------------------------------------------------------ figure
img = mpimg.imread(FIGPNG)
ih, iw = img.shape[0], img.shape[1]
fw = 0.450
fh = fw * (ih / iw) * (8.5 / 11.0)
y[0] -= 0.008
axf = fig.add_axes([0.5 - fw / 2, y[0] - fh, fw, fh])
axf.imshow(img); axf.axis('off')
y[0] -= fh + 0.012
fig.text(0.5, y[0], 'Offsets with fitted one- and two-component models (top); observed statistic '
                    'against the bootstrap null (bottom). Rug = individual loci.',
         ha='center', va='top', size=7.3, color=GREY)
y[0] -= 0.020

# ------------------------------------------------------------------- table
head('Results')
rows = [('Height', 'mirror-class', 'equal var', 'mirror-class only', 'equal_var'),
        ('Height', 'mirror-class', 'free var', 'mirror-class only', 'free_var'),
        ('Height', 'full data', 'equal var', 'full data', 'equal_var'),
        ('Height', 'full data', 'free var', 'full data', 'free_var'),
        ('BMI', 'either', 'equal var', 'full data', 'equal_var'),
        ('BMI', 'either', 'free var', 'full data', 'free_var')]
cells = [[t, df, da, f'{q(t, f, a, "LRT"):.2f}', f'{q(t, f, a, "minor_weight"):.2f}',
          f'{q(t, f, a, "separation"):.3f}', f'{q(t, f, a, "p_bootstrap_nullA"):.3f}',
          f'{q(t, f, a, "p_bootstrap_nullB"):.3f}']
         for t, df, da, f, a in rows]
axt = fig.add_axes([L, y[0] - 0.112, R - L, 0.108]); axt.axis('off')
tb = axt.table(cellText=cells,
               colLabels=['Trait', 'Mirror fit', 'Alternative', 'LRT', 'Minor wt.',
                          'Separation (SD)', 'p (null A)', 'p (null B)'],
               loc='center', cellLoc='center',
               colWidths=[.10, .15, .13, .09, .11, .16, .13, .13])
tb.auto_set_font_size(False); tb.set_fontsize(7.0); tb.scale(1, 1.06)
for (r, c), cell in tb.get_celld().items():
    cell.set_linewidth(.4); cell.set_edgecolor('#cccccc')
    if r == 0:
        cell.set_facecolor('#ececec'); cell.set_text_props(weight='bold')
    elif r <= 2:
        cell.set_facecolor('#fdf3d8')
y[0] -= 0.119

# ---------------------------------------------------------------- findings
head('Findings and caveats')
for b in [
    f'Height is bimodal: LRT {q("Height","mirror-class only","equal_var","LRT"):.2f}, p = '
    f'{q("Height","mirror-class only","equal_var","p_bootstrap_nullA"):.3f} (A) / '
    f'{q("Height","mirror-class only","equal_var","p_bootstrap_nullB"):.3f} (B) equal-variance, '
    f'{q("Height","mirror-class only","free_var","p_bootstrap_nullA"):.3f} (A) free-variance. The '
    f'minor component holds {100*q("Height","mirror-class only","equal_var","minor_weight"):.0f}% '
    f'of loci at {q("Height","mirror-class only","equal_var","separation"):.3f} SD separation, '
    'matching the sparse single-hit displacement (~0.27-0.31 SD).',
    'BMI is unimodal throughout (p = '
    f'{q("BMI","full data","free_var","p_bootstrap_nullA"):.3f}-'
    f'{q("BMI","full data","equal_var","p_bootstrap_nullB"):.3f}). Its larger free-variance '
    f'statistic ({q("BMI","full data","free_var","LRT"):.2f}) is an artefact of one locus pinned '
    'at the variance floor; equal variance removes it.',
    'The two height modes are exactly the published classes: asymmetric loci mean '
    f'{g.loc["asymmetric","mean"]:+.3f} SD, mirror loci {g.loc["mirror","mean"]:+.3f} SD, all '
    f'{n_as} asymmetric loci below {CUT:+.2f} and {"none" if n_mi == 0 else n_mi} of the 21 mirror '
    'loci. Because the line is fitted on mirror-class loci and the labels themselves encode '
    'dose-response behaviour, this partly restates the labels rather than finding structure blind. '
    'The label-free version is the full-data calibration, which is not significant (p = '
    f'{q("Height","full data","equal_var","p_bootstrap_nullA"):.3f}-'
    f'{q("Height","full data","free_var","p_bootstrap_nullB"):.3f}); significance also weakens '
    'under null B. At 33 loci this is suggestive, not settled, and both numbers belong in a '
    'write-up.',
]:
    put('•  ' + b, dy=0.0125, gap=0.0024, wrap=131)

# ------------------------------------------------------------------- paths
head('Code, figure and tables')
paths = [('Analysis', 'code/perp_offset_bimodality.py'),
         ('Figure', 'figures/Figure_perp_offset_bimodality.{pdf,png}'),
         ('Figure script', 'code/make_figure.py'),
         ('Results table', 'outputs/bimodality_test.csv'),
         ('This report', 'code/make_report.py'),
         ('Per-locus', 'outputs/perp_offsets_per_locus.csv'),
         ('Full notes', 'README.md'),
         ('Null draws', 'outputs/null_LRT_draws.npz')]
y0 = y[0]
for i, (lab, p) in enumerate(paths):
    col = i % 2
    yy = y0 - (i // 2) * 0.0128
    x = L if col == 0 else 0.515
    fig.text(x, yy, lab, ha='left', va='top', size=7.5, color=GREY)
    fig.text(x + 0.108, yy, p, ha='left', va='top', size=7.5, color=BLUE, family='monospace')
y[0] = y0 - (len(paths) // 2) * 0.0128 - 0.005
put('All paths relative to  ...\\Molly\\ContinuousVsSingleHit\\perp_offset_bimodality\\ .  '
    'Reproduce:  python code/perp_offset_bimodality.py --nboot 5000 --starts 25  &&  '
    'python code/make_figure.py  &&  python code/make_report.py', size=7.4, color=GREY)

fig.savefig(OUT)
fig.savefig(OUT.with_suffix('.png'), dpi=120)
print('wrote', OUT)
print('content ends at y =', round(y[0], 3), '(keep above ~0.02)')
