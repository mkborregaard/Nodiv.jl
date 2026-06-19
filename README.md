# Nodiv

[![Build Status](https://github.com/mkborregaard/Nodiv.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mkborregaard/Nodiv.jl/actions/workflows/CI.yml?query=branch%3Amain)

Node-based analysis of species distributions — a Julia implementation of the method in
Borregaard et al. (2014, *Methods in Ecology and Evolution* 5: 1225–1235). For every
internal node of a phylogeny it compares the per-cell richness of the two descendant
clades against a null model, giving a per-cell **SOS** (specific overrepresentation
score) and a single per-node divergence score.

The null model is matrix swapping (curveball): it holds each species' range size and
each cell's richness constant, so the divergence it measures is *spatial* — where the
two clades occur — rather than a difference in how species-rich they are.

```julia
res = node_metrics(assemblage, tree)   # one pass; cache with JLD2 and explore
divergent_nodes(res)                    # nodes above the RMS-SOS cutoff
res.sos["Node 123"]                     # per-cell SOS map for a node
```

## Divergence scores

`node_metrics` returns several per-node summaries of the SOS field:

- `rms` — **RMS-SOS**, the recommended score. The root-mean-square SOS over the focal
  clade's occupied cells, in units of null SDs: ≈ 1 means no divergence, > 1 means
  divergence. It is an effect size, so it is stable across replication counts and grid
  resolutions, and resolves the strongly divergent nodes that GND compresses.
- `spatial` — RMS-SOS with any uniform offset removed (spatial turnover only).
- `ses` — standardised effect size of the divergence. Useful within one analysis, but
  its magnitude inflates with clade size, so do **not** compare it across nodes or grains.
- `pval` — a null-calibrated Monte-Carlo significance. Convenient as a threshold
  (`divergent_nodes(res; by = :pval)`) but, being a significance, larger clades clear a
  fixed cutoff more easily.
- `gnd` — the original GND of Borregaard et al. (2014). Retained for continuity. Note
  that GND is a logit-mean of per-cell *p*-values rather than an effect size: its
  no-divergence baseline is ≈ 0.5 (not 0), its maximum is `1 − O(1/nsims)`, and it
  saturates, so the most divergent nodes pile up near the ceiling. Prefer `rms`.

Because divergence tends to scale with clade richness, `size_residual(res, tree)` gives
RMS-SOS residualised on log species count — "more divergent than expected for a clade
this size."
