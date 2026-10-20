#!/usr/bin/env bash
# =============================================================================
# MSMC2 demographic inference (Figure S3).
# Standard msmc-tools workflow: per-sample consensus + callable mask,
# multihetsep merging, MSMC2 on 8 haplotypes per population.
# Scaling: mu = 2.03e-9 per site per generation; generation time 2.5 yr.
# =============================================================================
set -euo pipefail

# ---------- CONFIG ----------
REF=reference.fna
BAMDIR=bams
POPS=(Azores Verde Canary Spain)
SEX_CONTIGS=sex_contigs.list
OUT=msmc2_out; mkdir -p "$OUT"

# Four individuals per population (8 haplotypes); list the sample IDs here.
declare -A INDS=(
  [Azores]="Azore-3 Azore-4 Azore-5 Azore-6"
  [Verde]="Verde-1 Verde-2 Verde-3 Verde-4"
  [Canary]="Canary-1 Canary-2 Canary-3 Canary-4"
  [Spain]="Spain-1 Spain-2 Spain-3 Spain-4"
)

for pop in "${POPS[@]}"; do
  read -ra S4 <<< "${INDS[$pop]}"
  masks=(); vcfs=()
  for s in "${S4[@]}"; do
    bam="$BAMDIR/$s.dedup.bam"
    for chr in $(grep -vFf "$SEX_CONTIGS" <(cut -f1 "$REF".fai | sort -u)); do
      samtools mpileup -q 20 -Q 20 -C 50 -u -r "$chr" -f "$REF" "$bam" | \
        bcftools call -c -V indels | \
        bamCaller.py 4 "$OUT/$s.$chr.mask.bed" | gzip -c > "$OUT/$s.$chr.vcf.gz"
      bgzip -f "$OUT/$s.$chr.mask.bed"
      masks+=(--mask "$OUT/$s.$chr.mask.bed.gz")
      vcfs+=("$OUT/$s.$chr.vcf.gz")
    done
  done
  generate_multihetsep.py "${masks[@]}" "${vcfs[@]}" > "$OUT/$pop.multihetsep.txt"

  msmc2 -t 8 -I 0-1,2-3,4-5,6-7 -p 1*2+25*1+1*2+1*3 \
      -o "$OUT/$pop.msmc2" "$OUT/$pop.multihetsep.txt"
done
