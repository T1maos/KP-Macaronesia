#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cadd_load.py
============
CADD-based mutation load per individual (Lin et al., Kentish plover genomics).

Upstream processing:
  1. Chicken CADD scores for assembly GCF_016699485.2 (Lensing et al. 2026;
     https://data.4tu.nl/datasets/f2ff2a38-0766-48f0-99f1-65d875ba81d4).
  2. Chicken coordinates lifted galGal6 -> GRCg7b with UCSC liftOver and the
     official chain file galGal6ToGRCg7b_num.over.chain.gz:
         liftOver chicken_snps_galGal6.bed galGal6ToGRCg7b_num.over.chain.gz \
             chicken_snps_GRCg7b.bed unmapped.bed
  3. Kentish plover SNPs assigned the CADD score of their syntenic chicken
     position (synteny map from Satsuma v2.0, chicken GCA_000002315.5).
  4. Potentially deleterious sites retained at CADD > 10
     (Rentzsch et al. 2019).

This script computes, per individual:
  * number of retained deleterious sites (dosage > 0);
  * total CADD load  = sum(CADD * dosage), dosage in {0, 1, 2}
                       (Bertorelle et al. 2022; Bortoluzzi et al. 2024);
  * average CADD     = total load / retained sites.

Input:
  --sites : tab-separated  CHROM  POS  CADD   (retained sites, CADD > 10)
  --vcf   : SNP VCF with genotypes for the focal individuals

Usage:
  python cadd_load.py --sites cadd_gt10.tsv --vcf focal.vcf.gz --out cadd_load.tsv
"""
import argparse
import gzip
from collections import defaultdict

def open_any(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', required=True)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    cadd = {}
    with open_any(args.sites) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            chrom, pos, score = line.rstrip('\n').split('\t')[:3]
            cadd[(chrom, int(pos))] = float(score)

    samples = None
    total = defaultdict(float)
    nsites = defaultdict(int)

    with open_any(args.vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
                continue
            f = line.rstrip('\n').split('\t')
            key = (f[0], int(f[1]))
            if key not in cadd:
                continue
            score = cadd[key]
            for i, s in enumerate(samples):
                gt = f[9 + i].split(':')[0].replace('|', '/')
                if gt not in ('0/0', '0/1', '1/1'):
                    continue
                dosage = gt.count('1')
                if dosage > 0:
                    nsites[s] += 1
                    total[s] += score * dosage

    with open(args.out, 'w') as out:
        out.write('sample\tdeleterious_sites\ttotal_CADD_load\taverage_CADD\n')
        for s in sorted(nsites):
            out.write(f'{s}\t{nsites[s]}\t{total[s]:.3f}\t{total[s]/nsites[s]:.4f}\n')

if __name__ == '__main__':
    main()
