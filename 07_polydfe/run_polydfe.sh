#!/usr/bin/env bash
# =============================================================================
# polyDFE — DFE inference (input for the SLiM forward simulation, 05_slim/).
# Tataru & Bataillon (2020).
#
# The distribution of fitness effects provides the deleterious:neutral
# mutation ratio (r_d), the mean selection coefficient of deleterious
# mutations (s_d) and the shape parameter of the gamma DFE (b), which are
# used in 05_slim/plover_forward.slim.
# =============================================================================
set -euo pipefail

# ---------- CONFIG ----------
VCF=../04_load/snpeff_annotated.vcf.gz   # SnpEff-annotated VCF (canonical exons)
OUT=polydfe_out; mkdir -p "$OUT"

# 1) Annotated VCF -> polyDFE input
#    synonymous sites = neutral reference class, nonsynonymous = selected class
#    (parse_vcf.py is distributed with polyDFE)
python parse_vcf.py -i "$VCF" -o "$OUT/dfe_input.txt"

# 2) DFE inference — full model (model D), gamma DFE
polyDFE -d "$OUT/dfe_input.txt" -m D -e > "$OUT/dfe_modelD.txt"

# 3) Export the three values used by the SLiM script
grep -E "^(r_d|s_d|b|s_bar)" "$OUT/dfe_modelD.txt" | tee "$OUT/dfe_for_slim.txt"
echo "Copy r_d / s_d / b into the constants block of 05_slim/plover_forward.slim"
echo "(remember the x100 parameter rescaling: s values are multiplied by 100)."
