#!/usr/bin/env bash
# =============================================================================
# Stairway Plot 2 — folded SFS from ANGSD + easySFS (no MAF filter;
# sex-linked loci excluded). 85% of sites randomly selected for training.
# Scaling: mu = 8.11e-10 per site per year (x 2.5 yr = 2.03e-9 per generation).
# =============================================================================
set -euo pipefail

# ---------- CONFIG ----------
REF=reference.fna
AUTOSOME_BED=autosomes.bed          # sex-linked contigs excluded
POPS=(Azores Verde Canary Spain)
BAMDIR=bams
MU_GEN=2.03e-9                      # per site per generation
GEN=2.5                             # years per generation
OUT=stairway_out; mkdir -p "$OUT"

for pop in "${POPS[@]}"; do
  ls "$BAMDIR"/$pop/*.dedup.bam > "$OUT/$pop.bamlist"

  # 1) Site allele frequency likelihoods (ANGSD; no MAF filter)
  angsd -bam "$OUT/$pop.bamlist" -ref "$REF" -anc "$REF" \
        -rf "$AUTOSOME_BED" \
        -GL 1 -doSaf 1 -doCounts 1 -P 8 -out "$OUT/$pop"

  # 2) SFS estimation and folding
  realSFS "$OUT/$pop.saf.idx" > "$OUT/$pop.sfs"
  easySFS.py -i "$OUT/$pop.saf.idx" -p "$OUT/$pop.bamlist" -a -f \
      -o "$OUT/$pop.easySFS"

  # 3) Stairway Plot 2 blueprint
  #    nseq = 2 x (number of diploid individuals projected);
  #    L    = total callable sequence length after masking;
  #    mu   = $MU_GEN (per site per generation);
  #    pct_training = 0.85
  # java -cp stairway_plot_es Stairbuilder $pop.blueprint
  # bash $pop.blueprint.sh
done
