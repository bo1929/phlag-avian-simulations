# Avian gene tree simulations used to test Phlag
This repository contains the simulated gene trees, the input species trees, and the simulation scripts used to benchmark Phlag.

## Files
- `phlag-avian-simulations/estimated-genetrees`: Gene trees that were estimated from sequences simulated by msprime using IQTREE.
  * `avian-500Kb-{IDX}.gtrees`: gene trees for the entire avian tree, each `IDX` is independent, `concat-avian.gtrees` is the concatenated version of all these blocks.
  * `neoaves-500Kb-{IDX}.gtrees`: gene trees for neoaves, each `IDX` is independent, `concat-neoaves.gtrees` is the concatenated version of all these blocks.
  * `recombination_increase_10x-{BRANCH}-500Kb-{IDX}.gtrees`: 10x recombination rate increase on branch `BRANCH`, the rest is the same as the baseline simulations.
  * `recombination_suppression-{BRANCH}-500Kb-{IDX}.gtrees`: 1000x recombination rate decrease on branch `BRANCH`, the rest is the same as the baseline simulations.
  * `popsize_increase_10x-{BRANCH}-500Kb-{IDX}.gtrees`: 10x population size increase on branch `BRANCH`, the rest is the same as the baseline simulations.
  * `popsize_decrease_10x-{BRANCH}-500Kb-{IDX}.gtrees`: 10x population size decrease on branch `BRANCH`, the rest is the same as the baseline simulations.
  * For all scenarios, files with `concat-*` prefix in their name contain the same gene trees in `IDX` but in a concatenated form (in a random order).

- `phlag-avian-simulations/main-speciestrees`: Species trees with branch lengths in different units.
  * `63K.tre` and `63K_dated.tre`: The avian trees (main Stiller2024) with CU branch lengths and in time (million years), respectively.
  * `castlespro_stiller.rooted.tre` and `castlespro_stiller.tre` in substitution units, branch lengths estimated by CASTLES-pro.
  * `main-avian-numgen.nwk` and `main-neoaves-numgen.nwk` trees used in ARG simulations, branch lengths in number of generations.

- `phlag-avian-simulations/misc/`: estimated rates, population sizes, and number of generations in a tabular format. Internal node labels match the species tree. The generation time was kept fixed at 10 years across the tree.

## Description
We use msprime to simulate genealogies under the Hudson coalescent model, using the 363-taxon avian phylogeny from Stiller et al. as our null
demographic model. Model parameters are set using empirical branch length estimates in time, CU, and SU units, with effective sizes (2Ne) derived from
branch lengths and a generation time of 10. We simulate 6Mbp alignments under GTR with rate heterogeneity, then estimate gene trees for 1500 loci
(500bp each) using IQ-TREE.

To create replicates (i.e., `concat`), we mix null and alternative gene trees by substituting a random subsequence (2–25% of trees) with ones simulated under an alternative condition — keeping everything the same except one parameter on one branch. Changed parameters include population size (10x up or down, 40 branches) and recombination rate (10x increase or 1000x suppression, 25 branches each).
