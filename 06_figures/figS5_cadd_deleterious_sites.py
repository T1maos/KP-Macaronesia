# -*- coding: utf-8 -*-
"""Figure S5: per-individual deleterious site counts (CADD>10), Fig5 style, Dunn BH brackets"""
import zipfile, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------- 1. read data ----------
z = zipfile.ZipFile(r'D:\桌面\CADD_files\CADD_files\R\over10_cadd_ave_load.xlsx')
ss = [a or b for a, b in re.findall(r'<si><t[^>]*>([^<]*)</t></si>|<si><t[^>]*/>([^<]*)',
                                    z.read('xl/sharedStrings.xml').decode('utf-8'))]
x = z.read('xl/worksheets/sheet1.xml').decode('utf-8')
data = {}
for r in re.findall(r'<row[^>]*>(.*?)</row>', x, re.S)[1:]:
    cells = []
    for col, t, v in re.findall(r'<c r="([A-Z]+)\d+"(?: t="(\w+)")?>(?:<v>([^<]*)</v>)?', r):
        cells.append(ss[int(v)] if t == 's' else v)
    if len(cells) >= 6:
        data.setdefault(cells[1], []).append(float(cells[3]))  # Site_Count

pops = ['Azores', 'Cabo Verde', 'Canary', 'Spain']
key = {'Azores': 'Azore', 'Cabo Verde': 'Verde', 'Canary': 'Canary', 'Spain': 'Spain'}
vals = [data[key[p]] for p in pops]
for p, v in zip(pops, vals):
    print(p, len(v), f'mean={np.mean(v):.0f} range={min(v):.0f}-{max(v):.0f}')

# ---------- 2. plot ----------
plt.rcParams['font.family'] = 'Arial'
colors = ['#4a1486', '#d73027', '#f46d43', '#313695']
fig, ax = plt.subplots(figsize=(4.4, 3.9), dpi=300)

vp = ax.violinplot([v for v in vals], positions=[1, 2, 3, 4], widths=0.75,
                   showmeans=False, showmedians=False, showextrema=False)
for body, c in zip(vp['bodies'], colors):
    body.set_facecolor(c); body.set_alpha(0.45); body.set_edgecolor('#666666'); body.set_linewidth(0.8)

bp = ax.boxplot(vals, positions=[1, 2, 3, 4], widths=0.32, patch_artist=True,
                showfliers=False, medianprops=dict(color='white', lw=1.4),
                whiskerprops=dict(color='#333333', lw=0.9), capprops=dict(color='#333333', lw=0.9))
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c); patch.set_edgecolor('#333333'); patch.set_linewidth(0.9)

rng = np.random.default_rng(7)
for i, v in enumerate(vals):
    ax.scatter(1 + i + rng.uniform(-0.09, 0.09, len(v)), v,
               s=11, color='#4d4d4d', alpha=0.65, zorder=3, linewidths=0)

# Dunn BH brackets (Site amount): *** <0.001, ** <0.01
pairs = [((2, 3), '***'), ((1, 3), '***'), ((2, 4), '**'), ((1, 4), '**')]
y = 21550; step = 520
for (a, b), star in pairs:
    ax.plot([a, a, b, b], [y, y + 130, y + 130, y], color='#333333', lw=0.9, clip_on=False)
    ax.text((a + b) / 2, y + 160, star, ha='center', va='bottom', fontsize=9, color='#222222')
    y += step

ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(pops, fontsize=10)
ax.set_ylabel('Deleterious sites', fontsize=11)
ax.set_xlabel('Population', fontsize=11)
ax.set_ylim(18400, 23900)
ax.yaxis.grid(True, color='#e0e0e0', lw=0.7); ax.set_axisbelow(True)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_color('#999999')
ax.tick_params(colors='#555555', labelsize=10)
fig.tight_layout()
fig.savefig('figS5.png', dpi=300)
print('figS5.png saved')

# ---------- 3. insert into supplementary ----------
CX, CY = int(4.4 * 914400), int(3.9 * 914400)
RID = 'rId15'
img_para = (
    '<w:p><w:pPr><w:rPr><w:rFonts w:hint="eastAsia"/></w:rPr></w:pPr><w:r>'
    '<w:rPr><w:rFonts w:ascii="Times New Roman" w:eastAsia="宋体" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
    '<w:noProof/><w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr>'
    f'<w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0">'
    f'<wp:extent cx="{CX}" cy="{CY}"/><wp:effectExtent l="0" t="0" r="2540" b="0"/>'
    '<wp:docPr id="40" name="图片 40"/>'
    '<wp:cNvGraphicFramePr><a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/></wp:cNvGraphicFramePr>'
    '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:nvPicPr><pic:cNvPr id="40" name="图片 40"/><pic:cNvPicPr/></pic:nvPicPr>'
    f'<pic:blipFill><a:blip r:embed="{RID}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
    f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{CX}" cy="{CY}"/></a:xfrm>'
    '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>'
    '</a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')

legend = ('<w:p><w:pPr><w:rPr><w:rFonts w:hint="eastAsia"/></w:rPr></w:pPr>'
          '<w:r><w:rPr><w:rFonts w:hint="eastAsia"/></w:rPr>'
          '<w:t xml:space="preserve">Figure S5. Number of deleterious sites per individual (CADD score &gt; 10) '
          'in the three Macaronesian island populations and the mainland Spain population. '
          'Violin plots show per-population distributions with embedded boxplots (median and interquartile range) '
          'and individual data points; brackets denote significant pairwise differences '
          '(Dunn\u2019s post-hoc test, Benjamini\u2013Hochberg adjusted; ** P &lt; 0.01, *** P &lt; 0.001; Table S8).'
          '</w:t></w:r></w:p>')

p = 'supp2_doc.xml'
xml = open(p, encoding='utf-8').read()
m = re.search(r'<w:p [^>]*>(?:(?!</w:p>).)*?Figure S4\. Simulation of the best population.*?</w:p>', xml, re.S)
assert m, 'S4 legend not found'
xml = xml.replace(m.group(0), m.group(0) + img_para + legend, 1)
open(p, 'w', encoding='utf-8').write(xml)

# repack with new media + rel
zin = zipfile.ZipFile('supp2.docx', 'r')
zout = zipfile.ZipFile('supp3_new.docx', 'w', zipfile.ZIP_DEFLATED)
for item in zin.infolist():
    data = zin.read(item.filename)
    if item.filename == 'word/document.xml':
        data = open(p, 'rb').read()
    elif item.filename == 'word/_rels/document.xml.rels':
        rels = data.decode('utf-8')
        assert RID not in rels
        newrel = (f'<Relationship Id="{RID}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" '
                  'Target="media/image9.png"/>')
        rels = rels.replace('</Relationships>', newrel + '</Relationships>')
        data = rels.encode('utf-8')
    zout.writestr(item, data)
zout.writestr('word/media/image9.png', open('figS5.png', 'rb').read())
zout.close(); zin.close()

# verify
z2 = zipfile.ZipFile('supp3_new.docx')
from defusedxml import ElementTree as ET
x2 = z2.read('word/document.xml').decode('utf-8')
ET.fromstring(x2)
ET.fromstring(z2.read('word/_rels/document.xml.rels'))
print('良构 OK, image9.png' , 'word/media/image9.png' in z2.namelist(), ', Figure S5 legend:', 'Figure S5. Number of deleterious sites' in x2)
