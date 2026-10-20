# -*- coding: utf-8 -*-
"""Dunn test (BH) on 36-individual CADD metrics -> dunn36_results.tsv"""
import csv, re
import numpy as np
from scipy.stats import rankdata, norm
from dunn import dunn_test, bh_adjust, stars

rows = list(csv.DictReader(open('cadd_36.tsv', encoding='utf-8'), delimiter='\t'))
pops = ['Azore', 'Verde', 'Canary', 'Spain']
metrics = {'Total CADD load': 'Total_Load', 'Deleterious sites': 'Site_Count', 'Average CADD score': 'Ave_CADD'}

out = open('dunn36_results.tsv', 'w', encoding='utf-8')
out.write('Metric\tComparison\tZ\tP.unadj\tP.adj\tSignif\n')
for mname, col in metrics.items():
    groups = [np.array([float(r[col]) for r in rows if r['Population'] == p]) for p in pops]
    res = dunn_test(groups, pops)
    keys = list(res.keys())
    padj = bh_adjust([res[k][1] for k in keys])
    print(f'== {mname} ==')
    for k, pa in zip(keys, padj):
        z, p = res[k]
        comp = f'{k[0]} - {k[1]}'
        print(f'{comp:20s} z={z:+8.4f}  P.unadj={p:.4e}  P.adj={pa:.4e}  {stars(pa)}')
        out.write(f'{mname}\t{comp}\t{z:.4f}\t{p:.4e}\t{pa:.4e}\t{stars(pa)}\n')
out.close()
print('saved dunn36_results.tsv')
