# Software and versions

All tools were run with default settings unless noted. Parameter values below
are those stated in the Materials and Methods of the manuscript.

## WGS pipeline (`01_wgs_pipeline/`)

| Step | Tool | Version | Key parameters |
|---|---|---|---|
| Read QC | FastQC | 0.11.8 | both mates |
| Trimming | Trimmomatic | 0.39 | `LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36` |
| Mapping | bwa-mem | 0.7.12 | defaults |
| Duplicates | Picard MarkDuplicates | 1.140 | mark only |
| Variant calling | GATK | 4.1.4.1 | HaplotypeCaller `-ERC GVCF` per sample → CombineGVCFs → GenotypeGVCFs |
| Hard filters | GATK VariantFiltration + SelectVariants | 4.1.4.1 | `QD < 2.0`, `FS > 60.0`, `MQ < 40.0`, `MQRankSum < -12.5`, `ReadPosRankSum < -8.0`; `--exclude-filtered` |
| Second-pass filters | VCFtools | 0.1.17 | biallelic SNPs, `--minGQ 30`, `--min-meanDP 4 --max-meanDP 60`, `--max-missing 0.95`; `--maf 0.1` except for SFS analyses |
| Sex-chromosome contigs | Satsuma | v2.0 | synteny vs chicken (GCA_000002315.5) |

## Population genetics (`02_popgen/`)

| Analysis | Tool | Parameters |
|---|---|---|
| Nucleotide diversity (π) | VCFtools | 10 kb windows, 5 kb step, autosomes |
| ROH | PLINK | v1.9, `--homozyg --homozyg-kb 100 --homozyg-snp 10 --homozyg-density 10` |
| F<sub>ROH</sub> | — | total ROH length (≥100 kb) / genome length (McQuillan et al. 2008) |
| LD decay | PopLDdecay | pairwise r² within 150 kb, per population |
| LD pruning | PLINK | `--indep-pairwise 100kb 20 0.2` |
| Pairwise F<sub>ST</sub> | VCFtools | Weir & Cockerham |
| Individual NJ tree | Vcf2Dis + PHYLIP | neighbor; White-faced plover outgroup |
| Population NJ tree | PAUP | 100 bootstrap replicates |
| PCA | PLINK | `--pca` on LD-pruned SNPs |
| Ancestry | ADMIXTURE 1.3 | K = 1–7, cross-validation, 200 bootstrap replicates |
| Ancestry | fastSTRUCTURE 1.0 | K = 1–7, `chooseK.py` for optimal K |
| Kinship | KING | 2.2.5 |

## Demographic history (`03_demography/`)

| Analysis | Tool | Parameters |
|---|---|---|
| PSMC | PSMC | `-N30 -t5 -r5 -p "4+30*2+4+6+10"`, 100 bootstrap replicates; sex-linked contigs and contigs < 10 kb excluded; 36 high-depth individuals (Azores 8, Cabo Verde 9, Canary 9, Spain 10) |
| Stairway Plot 2 | ANGSD + realSFS + easySFS + Stairway Plot 2 | folded SFS, no MAF filter, sex-linked loci excluded; 85% of sites used for training |
| MSMC2 | msmc-tools + MSMC2 | 4 individuals (8 haploids) per population |
| fastsimcoal2 | fastsimcoal2 | 12 candidate models; 2D SFS via RealSFS; best model: successive colonization Spain → Canary → Cabo Verde → Azores with gene flow |

Mutation rate and generation time throughout: 8.11 × 10⁻¹⁰ per site per year
(Charadrius-specific estimate; Wang X et al. 2019) × 2.5 yr =
2.03 × 10⁻⁹ per site per generation.

## Mutation load (`04_load/`)

| Analysis | Tool | Parameters |
|---|---|---|
| Annotation | SnpEff | canonical exons |
| Ancestral state | custom (`count_load_per_individual.py`) | allele ancestral if homozygous in any of the three outgroups (White-faced, Snowy, Kittlitz's plover) |
| Grantham classification | custom | benign < 50, deleterious > 50 (Grantham 1974); homozygous sites counted twice |
| LoF | SnpEff classes | stop-gain, stop-loss, frameshift |
| CADD | UCSC liftOver + custom (`cadd_load.py`) | chicken CADD scores (GCF_016699485.2; Lensing et al. 2026), liftOver galGal6 → GRCg7b, retained sites CADD > 10 (Rentzsch et al. 2019); load = Σ CADD × dosage (Bertorelle et al. 2022; Bortoluzzi et al. 2024) |
| Statistics | R, dunn.test | Kruskal–Wallis + Dunn's test, Benjamini–Hochberg FDR |

## Forward simulations (`05_slim/`, `07_polydfe/`)

| Item | Tool | Parameters |
|---|---|---|
| DFE | polyDFE | full model; deleterious:neutral ratio, mean s, and gamma shape exported to SLiM |
| Simulation | SLiM | v3.6, Wright–Fisher; 20,000 genes × 1,300 bp (Huber et al. 2017); mutation rate 2.03 × 10⁻⁷ per site per generation; recombination among genes 1 × 10⁻³ (Xie et al. 2022); dominance h = 0 if s < −0.01, else h = 0.25 (Kyriazis et al. 2021); all parameters rescaled by ×100; demography from the best fastsimcoal2 model; ten burn-in runs; 8 individuals (16 haploids) sampled per generation for the N/S ratio; 100,000 years forward at constant Ne |
