# Nodiv

[![Build Status](https://github.com/mkborregaard/Nodiv.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mkborregaard/Nodiv.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mkborregaard/Nodiv.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mkborregaard/Nodiv.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

Node-based analysis of species distributions: a Julia implementation of the method in Borregaard et al. (2014, *Methods in Ecology and Evolution* 5: 1225–1235).
For every internal node of a phylogeny it compares the per-cell richness of the two descendant clades against a null model.
This gives a per-cell **SOS** (specific overrepresentation score) and per-node divergence scores.

The null model is matrix swapping (curveball).
It holds each species' range size and each cell's richness constant, so the divergence it measures is *spatial*: where the two clades occur, rather than how species-rich they are.

```julia
res = node_metrics(assemblage, tree)   # one pass; cache with JLD2 and explore
divergent_nodes(res)                    # nodes above the RMS-SOS cutoff, in tree order
res.sos["Node 123"]                     # per-cell SOS map of a node
sos_distances(res, divergent_nodes(res)) # how alike the nodes' SOS maps are
```

`node_metrics` analyses the nodes in parallel, so start Julia with several threads (`julia -t auto`).
The results carry Monte Carlo noise and differ slightly between runs.

## Divergence scores

`node_metrics` returns several per-node summaries of the SOS field:

- `rms`: **RMS-SOS**, the recommended score (see `sos_rms`).
  The root-mean-square SOS over the cells where the null model varies, in units of null standard deviations: ≈ 1 means no divergence, > 1 means divergence.
  It is an effect size, so it is stable across replication counts and grid resolutions, and it resolves the strongly divergent nodes that GND compresses.
- `sd`: the standard deviation of the SOS over the same cells, i.e. RMS-SOS without any uniform offset (see `sos_sd`).
- `ses`: the standardised effect size of the divergence (see `divergence_ses`).
  Useful within one analysis, but its magnitude inflates with clade size, so do **not** compare it across nodes or grains.
- `pval`: a null-calibrated Monte Carlo significance (see `divergence_pval`).
  Convenient as a threshold (`divergent_nodes(res; by=:pval)`), but, being a significance, larger clades clear a fixed cutoff more easily.
- `gnd`: the original GND of Borregaard et al. (2014), retained for continuity (see `gnd`).
  GND is a logit-mean of per-cell *p*-values rather than an effect size.
  Its no-divergence baseline is ≈ 0.5 (not 0), its maximum is `1 − O(1/nsims)`, and it saturates, so the most divergent nodes pile up near the ceiling.
  Prefer `rms`.
- `varying`: the share of the node's occupied cells where the null model varies (see below).

The same scores can be computed from a simulation matrix (`simulate_descendants`) with `sos`, `sos_rms`, `sos_sd`, `divergence_ses`, `divergence_pval` and `gnd`.
`sos_rms` and `sos_sd` also take a cached SOS vector, e.g. `sos_rms(res.sos[node])`.

### Cells where the null model cannot vary

In some occupied cells the null model cannot change the richness of the descendant clades.
The clearest case is a cell holding every species of the parent clade: as the null model keeps each cell's richness, both descendants are then fully present in every draw.
Such cells have no SOS (`NaN`), and RMS-SOS, the SD, the SES and the *P* value leave them out.

This can be unintuitive where the two clades co-occur fully.
A cell holding every species of both daughters is, biologically, a place where they do not diverge, yet it does not lower RMS-SOS: the score describes only the cells where divergence can be detected.
So a clade whose daughters co-occur almost everywhere, and segregate in a few cells, can get a high RMS-SOS from those few cells.
`varying` shows how much of the range the score rests on: a low value means most of the clade's occupied cells could not be tested.
`sqrt(varying) * rms` gives the RMS-SOS over all occupied cells, with the fixed ones counted as SOS = 0.

GND treats these cells as showing no divergence (a *P* value near 1).
`sos_distances` correlates SOS maps over the cells where both SOS values are defined, so it leaves them out too.

## Plotting

Nodiv has Plots recipes (through RecipesBase, so it does not depend on Plots itself):

- `plot_gnd(tree, res)`: the GND of the analysed nodes on the tree.
- `plot_node(assemblage, tree, node, res)`: the node's clade, its SOS map and its two descendant clades.

[NodivMakie](https://github.com/mkborregaard/NodivMakie.jl) has Makie versions of these, including an interactive explorer.

## Deprecated names

| Deprecated | Use instead |
|---|---|
| `calculate_SOS` | `sos` |
| `calculate_GND` | `gnd` |
| `calculate_GND_rms`, `gnd_rms` | `sos_rms` |
| `calculate_GND_spatial`, `gnd_spatial` | `sos_sd` |
| `calculate_GND_ses` | `divergence_ses` |
| `calculate_GND_pval` | `divergence_pval` |
| `node_based_analysis` | `node_analysis` |
| `sos_distances(assemblage, tree, nodes)` | `sos_distances(res, nodes)` |

`calculate_GND_rms` and `calculate_GND_spatial` count occupied cells where the null model cannot vary as SOS = 0; `sos_rms` and `sos_sd` leave them out.
