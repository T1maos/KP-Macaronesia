#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
count_load_per_individual.py
============================
Per-individual mutation-load counts (Lin et al., Kentish plover genomics).

  * Ancestral state: an allele is regarded as ancestral when it is fixed
    (homozygous) in any one of the three outgroup species (White-faced
    plover, Snowy plover, Kittlitz's plover); outgroups are present as
    samples in the input VCF and listed via --outgroups.
  * Input: SnpEff-annotated VCF restricted to canonical exons.
  * Grantham classification: benign < 50, deleterious > 50
    (Grantham 1974, Science 185: 862-864, Table 2), computed separately for
    homozygous and heterozygous derived sites.
  * Loss-of-function: stop-gain, stop-loss and frameshift variants.
  * Genetic load: homozygous sites counted twice; total load = realized +
    masked load.
  * N/S ratio per individual from derived allele counts.

Usage:
  python count_load_per_individual.py --vcf annotated.vcf.gz \
      --outgroups WfP-1 WfP-2 ... SP-1 ... KiP-1 ... \
      --out load_per_individual.tsv
"""
import argparse
import gzip
import re
from collections import defaultdict

# Grantham (1974) amino-acid distance matrix, Table 2.
# Order: S T P A G N D E Q H R K M I L V F Y C W
_AA = "STPAGNDEQHRKMILVFYCW"
_MATRIX = [
[0],
[57],
[77, 96],
[58, 61, 107],
[85, 101, 69, 60],
[121, 138, 152, 116, 46],
[126, 147, 126, 111, 94, 23],
[133, 107, 91, 107, 56, 42, 0],
[148, 121, 93, 110, 74, 68, 45, 27],
[149, 140, 128, 128, 98, 87, 81, 83, 24],
[155, 121, 103, 106, 112, 99, 101, 56, 29, 89],
[154, 109, 95, 101, 127, 94, 90, 56, 53, 82, 26],
[133, 115, 94, 76, 82, 142, 119, 98, 101, 87, 96, 110, 21],
[98, 93, 74, 68, 94, 135, 158, 142, 130, 127, 111, 101, 69, 22],
[99, 78, 61, 62, 87, 134, 158, 145, 131, 123, 102, 96, 33, 5],
[81, 75, 64, 63, 109, 158, 177, 152, 147, 130, 115, 100, 70, 28, 19],
[145, 147, 153, 134, 143, 131, 168, 155, 144, 147, 133, 121, 108, 81, 82, 60],
[148, 122, 114, 113, 116, 116, 163, 152, 143, 134, 110, 112, 95, 74, 59, 28, 22],
[194, 215, 169, 154, 184, 169, 181, 179, 163, 184, 174, 167, 182, 154, 145, 112, 60, 36],
[174, 196, 160, 135, 190, 177, 177, 154, 177, 155, 155, 142, 169, 142, 141, 115, 37, 37, 184],
[198, 240, 196, 155, 199, 226, 221, 194, 235, 212, 177, 178, 205, 172, 167, 143, 81, 77, 215, 67],
]
_GRANTHAM = {}
for _i, a in enumerate(_AA):
    for _j, b in enumerate(_AA[: _i + 1]):
        _GRANTHAM[frozenset((a, b))] = _MATRIX[_i][_j]

_AA3 = {
 'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
 'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
 'Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*','Sec':'U','Pyl':'O'}
_HGVS_P = re.compile(r'p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2}|\*=)')
_LOF_EFFECTS = ('stop_gained', 'stop_lost', 'frameshift_variant')

def grantham(hgvs_p):
    m = _HGVS_P.search(hgvs_p or '')
    if not m:
        return None
    a1, a2 = _AA3.get(m.group(1)), _AA3.get(m.group(2).replace('=*', 'Ter'))
    if a1 is None or a2 is None or a1 == a2 or '*' in (a1, a2):
        return None
    return _GRANTHAM.get(frozenset((a1, a2)))

def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--outgroups', nargs='+', required=True,
                    help='sample IDs of the three outgroup species')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    outgroups = set(args.outgroups)

    stats = defaultdict(lambda: defaultdict(int))
    samples = None
    og_idx, smp_idx = [], []

    with open_vcf(args.vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
                og_idx = [i for i, s in enumerate(samples) if s in outgroups]
                smp_idx = [i for i, s in enumerate(samples) if s not in outgroups]
                continue
            f = line.rstrip('\n').split('\t')
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1:      # biallelic SNPs only
                continue
            # ancestral state: fixed in any outgroup
            anc = None
            for i in og_idx:
                gt = f[9 + i].split(':')[0].replace('|', '/')
                if gt == '0/0':
                    anc = ref; break
                if gt == '1/1':
                    anc = alt; break
            if anc is None:
                continue
            der_code = '1' if anc == ref else '0'

            info = f[7]
            ann = next((x for x in info.split(';') if x.startswith('ANN=')), None)
            if ann is None:
                continue
            first = ann[4:].split(',')[0].split('|')
            effect = first[1] if len(first) > 1 else ''
            hgvs_p = first[10] if len(first) > 10 else ''

            is_syn = 'synonymous_variant' in effect
            is_nonsyn = ('missense_variant' in effect) or effect in _LOF_EFFECTS
            is_lof = effect in _LOF_EFFECTS
            gscore = grantham(hgvs_p) if 'missense_variant' in effect else None

            for i in smp_idx:
                gt = f[9 + i].split(':')[0].replace('|', '/')
                if gt not in ('0/0', '0/1', '1/1'):
                    continue
                n_der = gt.count(der_code)
                if n_der == 0:
                    continue
                s = stats[samples[i]]
                if is_syn:
                    s['syn_derived_alleles'] += n_der
                if is_nonsyn:
                    s['nonsyn_derived_alleles'] += n_der
                if is_lof:
                    s['lof_hom' if n_der == 2 else 'lof_het'] += 1
                if gscore is not None:
                    if gscore > 50:
                        s['del_hom' if n_der == 2 else 'del_het'] += 1
                    elif gscore < 50:
                        s['ben_hom' if n_der == 2 else 'ben_het'] += 1

    with open(args.out, 'w') as out:
        out.write('sample\tsyn_derived_alleles\tnonsyn_derived_alleles\tN_S_ratio\t'
                  'del_het\tdel_hom\tben_het\tben_hom\tlof_het\tlof_hom\t'
                  'realized_load_del\tmasked_load_del\ttotal_load_del\n')
        for s in sorted(stats):
            d = stats[s]
            ns = d['nonsyn_derived_alleles'] / d['syn_derived_alleles'] if d['syn_derived_alleles'] else float('nan')
            realized = 2 * d['del_hom']          # homozygous counted twice
            masked = d['del_het']
            out.write(f"{s}\t{d['syn_derived_alleles']}\t{d['nonsyn_derived_alleles']}\t{ns:.5f}\t"
                      f"{d['del_het']}\t{d['del_hom']}\t{d['ben_het']}\t{d['ben_hom']}\t"
                      f"{d['lof_het']}\t{d['lof_hom']}\t{realized}\t{masked}\t{realized+masked}\n")

if __name__ == '__main__':
    main()
