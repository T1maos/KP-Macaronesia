# -*- coding: utf-8 -*-
"""36-individual CADD: redraw Fig S5, splice new Fig5 e/f panels, build Table S8."""
import csv, re, zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from defusedxml import ElementTree as ET

# ---------- data ----------
rows = list(csv.DictReader(open('cadd_36.tsv', encoding='utf-8'), delimiter='\t'))
dunn = list(csv.DictReader(open('dunn36_results.tsv', encoding='utf-8'), delimiter='\t'))
POPS = ['Azores', 'Cabo Verde', 'Canary', 'Spain']
POPKEY = {'Azores': 'Azore', 'Cabo Verde': 'Verde', 'Canary': 'Canary', 'Spain': 'Spain'}
POS = {'Azore': 1, 'Verde': 2, 'Canary': 3, 'Spain': 4}
COLORS = ['#4a1486', '#d73027', '#f46d43', '#313695']

def series(col):
    return [[float(r[col]) for r in rows if r['Population'] == POPKEY[p]] for p in POPS]

def sig_pairs(metric):
    out = []
    for d in dunn:
        if d['Metric'] != metric or d['Signif'] == 'ns':
            continue
        a, b = [s.strip() for s in d['Comparison'].split(' - ')]
        out.append(((POS[a], POS[b]), d['Signif']))
    return out

def bracket_levels(pairs):
    """assign vertical levels, short spans low, no overlap per level"""
    pairs = sorted(pairs, key=lambda x: (x[0][1] - x[0][0], x[0][0]))
    levels = []
    for (a, b), s in pairs:
        for lv in range(len(levels) + 1):
            if lv == len(levels):
                levels.append([])
            if all(b < x[0] - 0.25 or a > x[1] + 0.25 for x, _ in levels[lv]):
                levels[lv].append(((a, b), s))
                break
    res = []
    for lv, lst in enumerate(levels):
        for item in lst:
            res.append((item[0], item[1], lv))
    return res

from scipy.stats import gaussian_kde

def violin_manual(ax, vals, positions, width=0.75):
    for x, v, c in zip(positions, vals, COLORS):
        v = np.asarray(v)
        kde = gaussian_kde(v)
        span = v.max() - v.min()
        grid = np.linspace(v.min() - 0.25 * span, v.max() + 0.25 * span, 200)
        dens = kde(grid)
        dens = dens / dens.max() * width / 2
        ax.fill_betweenx(grid, x - dens, x + dens, facecolor=c, alpha=0.45,
                         edgecolor='#666666', linewidth=0.8, zorder=2)

def panel(ax, vals, sigs, ylabel, ylim, y0, step, ticklabels=POPS, letter=None,
          ylab_size=11, tick_size=10):
    violin_manual(ax, vals, [1, 2, 3, 4])
    bp = ax.boxplot(vals, positions=[1, 2, 3, 4], widths=0.32, patch_artist=True, showfliers=False,
                    medianprops=dict(color='white', lw=1.4),
                    whiskerprops=dict(color='#333333', lw=0.9), capprops=dict(color='#333333', lw=0.9))
    for patch, c in zip(bp['boxes'], COLORS):
        patch.set_facecolor(c); patch.set_edgecolor('#333333'); patch.set_linewidth(0.9)
    rng = np.random.default_rng(7)
    for i, v in enumerate(vals):
        ax.scatter(1 + i + rng.uniform(-0.09, 0.09, len(v)), v, s=11, color='#4d4d4d',
                   alpha=0.65, zorder=3, linewidths=0)
    for (a, b), star, lv in bracket_levels(sigs):
        yy = y0 + lv * step
        ax.plot([a, a, b, b], [yy, yy + step * 0.25, yy + step * 0.25, yy], color='#333333', lw=0.9, clip_on=False)
        ax.text((a + b) / 2, yy + step * 0.3, star, ha='center', va='bottom', fontsize=9, color='#222222')
    ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(ticklabels, fontsize=tick_size)
    ax.set_ylabel(ylabel, fontsize=ylab_size)
    ax.set_ylim(*ylim)
    ax.yaxis.grid(True, color='#e0e0e0', lw=0.7); ax.set_axisbelow(True)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_color('#999999')
    ax.tick_params(colors='#555555', labelsize=tick_size)
    if letter:
        ax.text(-0.14, 1.03, letter, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

plt.rcParams['font.family'] = 'Arial'

# ---------- Fig S5 (300 dpi) ----------
v = series('Site_Count')
fig, ax = plt.subplots(figsize=(4.4, 3.9), dpi=300)
panel(ax, v, sig_pairs('Deleterious sites'), 'Deleterious sites', (18400, 23800), 21500, 520)
ax.set_xlabel('Population', fontsize=11)
fig.tight_layout(); fig.savefig('figS5_new.png', dpi=300); plt.close(fig)
print('figS5_new.png')

# ---------- Fig5 panels e/f (100 dpi, exact px) ----------
v = series('Total_Load')
fig, ax = plt.subplots(figsize=(5.75, 4.45), dpi=100)
panel(ax, v, sig_pairs('Total CADD load'), 'Total CADD', (333000, 352800), 347900, 800, letter='e')
ax.set_xlabel('Population', fontsize=11, fontweight='bold')
fig.tight_layout(); fig.savefig('panelE.png', dpi=100); plt.close(fig)
print('panelE.png', Image.open('panelE.png').size)

v = series('Ave_CADD')
fig, ax = plt.subplots(figsize=(5.86, 4.45), dpi=100)
panel(ax, v, sig_pairs('Average CADD score'), 'Average CADD', (15.9, 18.75), 17.98, 0.17, letter='f')
ax.set_xlabel('Population', fontsize=11, fontweight='bold')
fig.tight_layout(); fig.savefig('panelF.png', dpi=100); plt.close(fig)
print('panelF.png', Image.open('panelF.png').size)

# ---------- splice into composite ----------
comp = Image.open('cur_x/word/media/image6.tmp').convert('RGB')
from PIL import ImageDraw
d = ImageDraw.Draw(comp)
d.rectangle([0, 850, comp.width, comp.height], fill=(255, 255, 255))
comp.paste(Image.open('panelE.png'), (5, 855))
comp.paste(Image.open('panelF.png'), (584, 855))
comp.save('image6_new.png')
print('image6_new.png', comp.size)
