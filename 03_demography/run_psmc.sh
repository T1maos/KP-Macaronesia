#!/usr/bin/env bash
# =============================================================================
# PSMC demographic inference (36 high-depth individuals:
# Azores 8, Cabo Verde 9, Canary 9, Spain 10)
# Parameters: -N30 -t5 -r5 -p "4+30*2+4+6+10", 100 bootstrap replicates.
# Sex-linked contigs and contigs < 10 kb are excluded.
# Scaling: mu = 2.03e-9 per site per generation; generation time 2.5 yr.
# =============================================================================
set -euo pipefail

# ---------- CONFIG ----------
REF=reference.fna
BAMLIST=psmc_bams.txt            # one deduplicated BAM per line (36 individuals)
SEX_CONTIGS=sex_contigs.list     # contigs mapping to chicken sex chromosomes
SHORT_CONTIGS=short_lt10kb.list  # contigs shorter than 10 kb
MU=2.03e-9
GEN=2.5
OUT=psmc_out; mkdir -p "$OUT"

cat "$SEX_CONTIGS" "$SHORT_CONTIGS" | sort -u > exclude_contigs.list

while read -r BAM; do
  s=$(basename "$BAM" .dedup.bam)
  # Diploid consensus restricted to autosomal contigs >= 10 kb;
  # depth limits (-d 4 -D 60) mirror the site-level coverage filter.
  bcftools mpileup -Ou -f "$REF" --regions-file ^exclude_contigs.list "$BAM" | \
    bcftools call -c - | vcfutils.pl vcf2fq -d 4 -D 60 > "$OUT/$s.fq"
  fq2psmcfa -q20 "$OUT/$s.fq" > "$OUT/$s.psmcfa"

  psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o "$OUT/$s.psmc" "$OUT/$s.psmcfa"

  # 100 bootstrap replicates
  splitfa "$OUT/$s.psmcfa" > "$OUT/$s.split.psmcfa"
  for i in $(seq 1 100); do
    psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -b \
      -o "$OUT/$s.boot.$i.psmc" "$OUT/$s.split.psmcfa"
  done

  # psmc_plot.pl -u $MU -g $GEN to scale and plot
done < "$BAMLIST"
