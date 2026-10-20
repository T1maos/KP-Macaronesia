// =============================================================================
// fastsimcoal2 template — best-supported demographic model (manuscript Fig 3c
// / Fig S4): successive colonization Spain -> Canary -> Cabo Verde -> Azores,
// with gene flow between neighbouring populations after the splits.
//
// Candidate models compared during model selection were generated from this
// template by permuting the colonization order and the migration edges.
// =============================================================================

//Number of population samples (demes)
4 samples to simulate
//Population effective sizes (number of genes)
NPOP_AZO
NPOP_CAV
NPOP_CAN
NPOP_ESP
//Sample sizes
8
9
9
10
//Growth rates: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0: ancestral epoch (Canary-CaboVerde-Azores pool vs Spain)
0	MIG_CAN_ESP	0	0
MIG_CAN_ESP	0	0	0
0	0	0	0
0	0	0	0
//Migration matrix 1: after the Canary split (CaboVerde-Azores pool vs Canary)
0	0	0	0
0	0	0	0
0	0	0	MIG_AZO_CAV
0	0	MIG_AZO_CAV	0
//Migration matrix 2: present epoch (all four demes)
0	0	MIG_AZO_CAV	0
0	0	0	0
MIG_AZO_CAV	0	0	MIG_CAV_CAN
0	0	MIG_CAV_CAN	0
//Historical event format: time source sink migrants new_size growth migr new_hist
4 historical event
TDIV_AZO 0 1 1 1 0 2 1
TDIV_CAV 1 2 1 1 0 1 1
TDIV_CAN 2 3 1 1 0 0 1
TEND    3 3 1 1 0 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut. rate + right margin
DNA 1000000 0 2.03E-9 0.33
