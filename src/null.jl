# Build a sampling distribution of a descendant clade's per-site richness under the
# published null model: curveball (swap) randomization of the parent clade's
# presence/absence matrix, holding each species' range size and each site's richness
# constant, then recompute the descendant clade's richness. Fixing range size is what
# makes the divergence "spatial" - a uniform richness asymmetry between the two
# descendants is not flagged, only differences in WHERE they occur.
function simulate_descendants(clade, tree, descendant; nsims = 200)
    _simulate_descendants(clade, nodespecies(tree, descendant); nsims)
end

# Core sampler: the descendant is given as a precomputed species-name vector, so the
# hot loop touches no tree functions. That keeps it thread-safe (the parallel path in
# `node_metrics` does all tree access up front) and avoids recomputing the descendant
# species set on every one of the `nsims` draws.
function _simulate_descendants(clade, descsp::Vector{String}; nsims = 200)
    ret = zeros(nsims + 1, nsites(clade))                      # richness values from the draws
    ret[1, :] = richness(view(clade, species = descsp))        # empirical richness in row 1
    rmg = matrixrandomizer(clade)
    for i in 2:nsims + 1
        ret[i, :] .= richness(view(rand!(rmg), species = descsp))   # simulated richness
    end
    ret
end
