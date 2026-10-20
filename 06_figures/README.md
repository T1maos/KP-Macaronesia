# 06_figures — plotting scripts

Scripts producing the main and supplementary figures of the manuscript.

| File | Figure |
|---|---|
| `fig5_grantham_panels.py` | Fig. 5 — Grantham mutation-load panels (36 high-depth individuals; includes the N/S ratio panel) |
| `fig5_cadd_panels.py` | Fig. 5 — CADD panels |
| `figS5_cadd_deleterious_sites.py` | Fig. S5 — number of deleterious sites per individual (CADD) |
| `tableS7_dunn_36indiv.py` | Table S7 — Dunn's post-hoc test on the 36-individual load metrics |

Note: the per-individual summary values are embedded at the top of each
script (they are the outputs of `04_load/`); regenerate the input tables with
`04_load/count_load_per_individual.py` and `04_load/cadd_load.py` if you wish
to recompute them from the VCF.

Scripts for Figs. 1–3 and S1–S4 (map/PCA/ADMIXTURE, diversity panels,
PSMC/Stairway/fastsimcoal panels) are R scripts archived alongside this file.
