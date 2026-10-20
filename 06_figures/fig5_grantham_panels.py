# -*- coding: utf-8 -*-
"""Redraw Fig5 panels a-d (36 indiv, Dunn BH brackets), splice; rebuild Table S8 with 7 metrics."""
import csv, re, zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image, ImageDraw
from scipy.stats import gaussian_kde

# ---------- data ----------
rows = list(csv.DictReader(open('load39.txt', encoding='utf-8-sig'), delimiter='\t'))
EXC = {'Col_11', 'Col_13', 'Col_29'}
r36 = [r for r in rows if r['Sample'] not in EXC]
POPS = ['Azores', 'Cabo Verde', 'Canary', 'Spain']
POS = {'Azores': 1, 'Cabo Verde': 2, 'Canary': 3, 'Spain': 4}
COLORS = ['#4a1486', '#d73027', '#f46d43', '#313695']

dunn = list(csv.DictReader(open('dunn36_grantham.tsv', encoding='utf-8'), delimiter='\t'))

def series(col):
    return [[float(r[col]) for r in r36 if r['Population'].replace('Cape Verde', 'Cabo Verde') == p] for p in POPS]

def sig_pairs(metric):
    out = []
    for d in dunn:
        if d['Metric'] != metric or d['Signif'] == 'ns':
            continue
        a, b = [s.strip().replace('Cape Verde', 'Cabo Verde') for s in d['Comparison'].split(' - ')]
        out.append(((POS[a], POS[b]), d['Signif']))
    return out

def bracket_levels(pairs):
    pairs = sorted(pairs, key=lambda x: (x[0][1] - x[0][0], x[0][0]))
    levels = []
    for (a, b), s in pairs:
        for lv in range(len(levels) + 1):
            if lv == len(levels):
                levels.append([])
            if all(b < x[0] - 0.25 or a > x[1] + 0.25 for x, _ in levels[lv]):
                levels[lv].append(((a, b), s))
                break
    return [(it[0], it[1], lv) for lv, lst in enumerate(levels) for it in lst]

def violin_manual(ax, vals, width=0.75):
    for x, v, c in zip([1, 2, 3, 4], vals, COLORS):
        v = np.asarray(v)
        kde = gaussian_kde(v)
        span = v.max() - v.min()
        grid = np.linspace(v.min() - 0.25 * span, v.max() + 0.25 * span, 200)
        dens = kde(grid); dens = dens / dens.max() * width / 2
        ax.fill_betweenx(grid, x - dens, x + dens, facecolor=c, alpha=0.45,
                         edgecolor='#666666', linewidth=0.8, zorder=2)

def panel(ax, vals, sigs, ylabel, letter, ypad, ystep, legend_below_top=None, ytick_fmt=None):
    violin_manual(ax, vals)
    bp = ax.boxplot(vals, positions=[1, 2, 3, 4], widths=0.32, patch_artist=True, showfliers=False,
                    medianprops=dict(color='white', lw=1.4),
                    whiskerprops=dict(color='#333333', lw=0.9), capprops=dict(color='#333333', lw=0.9))
    for patch, c in zip(bp['boxes'], COLORS):
        patch.set_facecolor(c); patch.set_edgecolor('#333333'); patch.set_linewidth(0.9)
    rng = np.random.default_rng(7)
    for i, v in enumerate(vals):
        ax.scatter(1 + i + rng.uniform(-0.09, 0.09, len(v)), v, s=11, color='#4d4d4d',
                   alpha=0.65, zorder=3, linewidths=0)
    vmax = max(max(v) for v in vals)
    spans = max(max(v) - min(v) for v in vals)
    y0 = vmax + 0.25 * spans + ypad
    for (a, b), star, lv in bracket_levels(sigs):
        yy = y0 + lv * ystep
        ax.plot([a, a, b, b], [yy, yy + ystep * 0.25, yy + ystep * 0.25, yy], color='#333333', lw=0.9, clip_on=False)
        ax.text((a + b) / 2, yy + ystep * 0.3, star, ha='center', va='bottom', fontsize=9, color='#222222')
    nlv = max((lv for _, _, lv in bracket_levels(sigs)), default=-1) + 1
    vmin = min(min(v) for v in vals)
    ax.set_ylim(vmin - 0.28 * spans - ypad, y0 + nlv * ystep + ystep * 0.9)
    ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(POPS, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xlabel('Population', fontsize=11, fontweight='bold')
    ax.yaxis.grid(True, color='#e0e0e0', lw=0.7); ax.set_axisbelow(True)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_color('#999999')
    ax.tick_params(colors='#555555', labelsize=10)
    if ytick_fmt:
        ax.yaxis.set_major_formatter(plt.FuncFormatter(ytick_fmt))
    ax.text(-0.14, 1.03, letter, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
    if legend_below_top is not None:
        handles = [Patch(facecolor=c, edgecolor='#333333', label=p) for c, p in zip(COLORS, POPS)]
        leg = ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.99, legend_below_top),
                        frameon=False, fontsize=10, handlelength=1.4, handleheight=1.0, borderaxespad=0)
        leg.set_zorder(5)

plt.rcParams['font.family'] = 'Arial'

# panel regions in composite: (x, y, w, h) — 整行清除避免原图残影
REG = {'a': (0, 0, 578, 443), 'b': (578, 0, 601, 443), 'c': (0, 443, 578, 410), 'd': (578, 443, 601, 410)}
SPECS = {
    'a': ('mask', 'Masked load', 'Masked load', 60, 130, None, None),
    'b': ('realize', 'Realized load', 'Realized load', 35, 70, 0.80, None),
    'c': ('total', 'Total load', 'Total load', 90, 190, None, None),
    'd': ('dnds', 'dn/ds', 'N/S ratio', 0.0006, 0.0013, None, lambda v, p: f'{v:.3f}'),
}
comp = Image.open('image6_new.png').convert('RGB')
d = ImageDraw.Draw(comp)
for L, (col, metric, ylabel, ypad, ystep, legtop, yfmt) in SPECS.items():
    vals = series(col)
    x, y, w, h = REG[L]
    fig, ax = plt.subplots(figsize=(w / 100, h / 100), dpi=100)
    panel(ax, vals, sig_pairs(metric), ylabel, L, ypad, ystep, legend_below_top=legtop, ytick_fmt=yfmt)
    fig.tight_layout()
    fn = f'panel_{L}.png'
    fig.savefig(fn, dpi=100); plt.close(fig)
    d.rectangle([x, y, x + w, y + h], fill=(255, 255, 255))
    comp.paste(Image.open(fn), (x, y))
    print(L, 'done', Image.open(fn).size)
comp.save('image6_v2.png')
print('composite saved')
