# fastsimcoal2 model files

`best_model.tpl` / `best_model.est` implement the best-supported demographic
model of the manuscript (Fig 3c / Fig S4): successive colonization
Spain → Canary → Cabo Verde → Azores, followed by gene flow between
neighbouring populations.

During model selection, twelve candidate models were compared; the alternative
topologies were generated from this template by permuting the colonization
order and the presence/absence of migration edges, keeping the search ranges
of `best_model.est` unchanged. Model likelihoods were compared with AIC after
100 independent optimization runs per model (`-n 100000 -N 40 -L 40`).

Input: 2D site frequency spectra built with RealSFS (ANGSD), following the
same filtering as the Stairway Plot analyses (`../stairway/run_stairway.sh`).
Mutation rate 2.03 × 10⁻⁹ per site per generation; generation time 2.5 yr.
