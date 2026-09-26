"""
    simulate_descendants(clade, tree, descendant; nsims=200) -> Matrix{Float64}

The richness of the `descendant` clade in each cell of `clade` (the parent clade's
assemblage, see [`get_clade`](@ref)), observed and under the null model. Row 1 holds the
observed richness, and each of the `nsims` further rows one draw of the null model, with
one column per cell. This is the `sims` matrix the divergence scores are computed from.

The null model is the curveball (swap) randomisation of the parent clade's
presence-absence matrix, which keeps each species' range size and each cell's richness.
So a uniform difference in richness between the two descendants is not divergence, only
a difference in where they occur.
"""
function simulate_descendants(clade, tree, descendant; nsims=200)
    return _simulate_descendants(clade, nodespecies(tree, descendant); nsims)
end

# Core sampler: the descendant is given as a precomputed species-name vector, so the
# hot loop touches no tree functions. That keeps it thread-safe (the parallel path in
# `node_metrics` does all tree access up front) and avoids recomputing the descendant
# species set on every one of the `nsims` draws.
function _simulate_descendants(clade, descsp::Vector{String}; nsims=200)
    ret = zeros(nsims + 1, nsites(clade))
    ret[1, :] = richness(view(clade; species=descsp))  # observed
    rmg = matrixrandomizer(clade)
    for i in 2:(nsims + 1)
        ret[i, :] .= richness(view(rand!(rmg); species=descsp))  # one null draw
    end
    return ret
end
