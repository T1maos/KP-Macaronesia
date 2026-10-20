# Genomic characterization of Macaronesian island Kentish plover populations

Analysis scripts accompanying the manuscript by Lin et al. on genomic
characterization of Macaronesian island populations of the Kentish plover
(*Anarhynchus alexandrinus*): genetic diversity, demographic history,
deleterious mutation load, and forward simulations.

## Contents

| Directory | Contents |
|---|---|
| `01_wgs_pipeline/` | Snakemake workflow: raw reads → joint genotyped, hard-filtered VCF (FastQC, Trimmomatic, bwa-mem, Picard MarkDuplicates, GATK HaplotypeCaller / CombineGVCFs / GenotypeGVCFs, hard filters) |
| `02_popgen/` | Nucleotide diversity, individual heterozygosity, ROH/F<sub>ROH</sub> (PLINK), LD decay and pruning (PopLDdecay/PLINK), pairwise F<sub>ST</sub>, NJ trees (Vcf2Dis/PHYLIP/PAUP), PCA, ADMIXTURE, fastSTRUCTURE |
| `03_demography/` | PSMC (`run_psmc.sh`), Stairway Plot 2 (`stairway/`), MSMC2 (`msmc2/`), fastsimcoal2 model files (`fastsimcoal/`) |
| `04_load/` | Per-individual mutation-load counts (Grantham / LoF / N-S ratio), CADD-based load, Dunn's post-hoc test |
| `05_slim/` | SLiM v3.6 Wright–Fisher forward simulation (`plover_forward.slim`) |
| `06_figures/` | Plotting scripts for main and supplementary figures |
| `07_polydfe/` | polyDFE DFE inference (input for the SLiM simulation) |

## Requirements

Software versions and parameters are listed in [SOFTWARE.md](SOFTWARE.md).

## Usage

Run the directories in numeric order. The WGS pipeline is driven by
Snakemake; edit `01_wgs_pipeline/config.yaml` to point at the reference
genome and the sample list, then:

```bash
cd 01_wgs_pipeline
snakemake --cores 32
```

Downstream scripts (`02_popgen/run_popgen.sh`, `03_demography/`,
`04_load/`) consume the filtered VCF and its derivatives; adjust the paths
in the CONFIG block at the top of each script.

## Data

Raw resequencing data are deposited in the NCBI Sequence Read Archive
(BioProject accession to be added upon release). Sample information,
including ring numbers and collection years, is provided in Supplementary
Table S1 of the manuscript.

## Citation

If you use these scripts, please cite the manuscript (Lin et al., in
preparation).
