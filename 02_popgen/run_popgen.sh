#!/usr/bin/env bash
# =============================================================================
# Population genetics analyses (Kentish plover, Lin et al.)
# Input: joint.filtered.vcf.gz from 01_wgs_pipeline
# Software versions and citations: SOFTWARE.md (repository root).
# =============================================================================
set -euo pipefail

# ---------- CONFIG ----------
VCF=../01_wgs_pipeline/results/vcf/joint.filtered.vcf.gz
AUTOSOME_BED=autosomes.bed   # contigs not mapping to chicken sex chromosomes
                             # (Satsuma v2.0 vs chicken GCA_000002315.5)
POPS=(Azores Verde Canary Spain Tuzla Farasan XJ QH ZS TS LYG RD CZ TW FZ XM SW BH ZJ HN)
POPDIR=pops                  # one sample list per population: pops/<pop>.txt
OUT=results_popgen
mkdir -p "$OUT"

# ---------- 0. Second-pass VCFtools filters ----------
# biallelic SNPs, GQ > 30, mean depth 4-60X, <5% missing, MAF >= 0.1
# (the MAF filter is omitted for SFS analyses; see 03_demography)
vcftools --gzvcf "$VCF" \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 30 --min-meanDP 4 --max-meanDP 60 \
  --max-missing 0.95 --maf 0.1 \
  --bed "$AUTOSOME_BED" \
  --recode --stdout | bgzip -c > "$OUT/autosomal.biallelic.maf01.vcf.gz"
tabix "$OUT/autosomal.biallelic.maf01.vcf.gz"

# ---------- 1. Nucleotide diversity (pi): 10 kb windows, 5 kb step ----------
for pop in "${POPS[@]}"; do
  vcftools --gzvcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
    --keep "$POPDIR/$pop.txt" \
    --window-pi 10000 --window-pi-step 5000 \
    --out "$OUT/pi.$pop"
done

# ---------- 2. Whole-genome heterozygosity per individual ----------
# heterozygous SNPs / individual sequencing depth / reference genome length
bcftools query -f '[%SAMPLE\t%GT\n]' "$OUT/autosomal.biallelic.maf01.vcf.gz" \
  | awk '$2=="0/1"||$2=="0|1"{c[$1]++} END{for(s in c) print s, c[s]}' \
  > "$OUT/het_counts.tsv"

# ---------- 3. ROH and F_ROH ----------
plink --vcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
  --homozyg --homozyg-kb 100 --homozyg-snp 10 --homozyg-density 10 \
  --allow-extra-chr --out "$OUT/roh"
# F_ROH per individual = total ROH length (>=100 kb) / reference genome length

# ---------- 4. LD decay and LD pruning ----------
for pop in "${POPS[@]}"; do
  vcftools --gzvcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
    --keep "$POPDIR/$pop.txt" --recode --stdout | bgzip -c > "$OUT/tmp.$pop.vcf.gz"
  PopLDdecay -InVCF "$OUT/tmp.$pop.vcf.gz" -OutStat "$OUT/lddecay.$pop" -MaxDist 150
  rm "$OUT/tmp.$pop.vcf.gz"
done
plink --vcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
  --indep-pairwise 100kb 20 0.2 --allow-extra-chr --out "$OUT/ldprune"
plink --vcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
  --extract "$OUT/ldprune.prune.in" --make-bed --allow-extra-chr \
  --recode vcf --out "$OUT/unlinked"

# ---------- 5. Pairwise Fst (Weir & Cockerham) ----------
for i in "${!POPS[@]}"; do
  for j in $(seq $((i+1)) $((${#POPS[@]}-1))); do
    vcftools --gzvcf "$OUT/autosomal.biallelic.maf01.vcf.gz" \
      --weir-fst-pop "$POPDIR/${POPS[$i]}.txt" \
      --weir-fst-pop "$POPDIR/${POPS[$j]}.txt" \
      --out "$OUT/fst.${POPS[$i]}.${POPS[$j]}"
  done
done

# ---------- 6. Individual-level NJ tree (Vcf2Dis + PHYLIP neighbor) ----------
# White-faced plover samples serve as the outgroup.
Vcf2Dis -InPut "$OUT/unlinked.vcf" -OutPut "$OUT/unlinked.pdis"
# PHYLIP neighbor on the distance matrix (see PHYLIP documentation).

# ---------- 7. Population-level NJ tree (PAUP; 100 bootstrap replicates) ----------
# Input: pairwise Fst matrix from step 5.

# ---------- 8. PCA ----------
plink --bfile "$OUT/unlinked" --pca --allow-extra-chr --out "$OUT/pca"

# ---------- 9. ADMIXTURE (K = 1-7, cross-validation, 200 bootstraps) ----------
for K in 1 2 3 4 5 6 7; do
  admixture --cv -B200 "$OUT/unlinked.bed" $K -j8 | tee "$OUT/admixture.K$K.log"
done
# Optimal K = lowest cross-validation error.

# ---------- 10. fastSTRUCTURE (K = 1-7) ----------
for K in 1 2 3 4 5 6 7; do
  structure.py -K $K --input="$OUT/unlinked" --output="$OUT/faststructure" --cv 5
done
chooseK.py --input="$OUT/faststructure"
