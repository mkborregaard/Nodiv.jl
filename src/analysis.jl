# A node can be analysed when it has two children with at least three species each in the
# assemblage (so neither child is a tip).
function _isanalysable(assemblage, tree, node)
    children = getchildren(tree, node)
    length(children) == 2 || return false
    return all(children) do child
        return !isleaf(tree, child) && nspecies(get_clade(assemblage, tree, child)) >= 3
    end
end

# The internal nodes in the order of a preorder traversal: the order they appear on the tree
function _internalnodes(tree)
    return [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
end

# Log progress every 100 analysed nodes; `done` is shared by the threads
function _progress!(done, total, label)
    n = Threads.atomic_add!(done, 1) + 1
    n % 100 == 0 && @info "$label: $n / $total analysable nodes done"
    return nothing
end

"""
    process_node(assemblage, tree, node; nsims=200) -> (sos_scores, gnd)

The per-cell [`sos`](@ref) and the [`gnd`](@ref) of one node, from `nsims` fresh draws
of the null model. A node that cannot be analysed gives all-`NaN` SOS and a `NaN` GND: a
node is analysed when it has two children with at least three species each in the
assemblage. To analyse the whole tree, use [`node_metrics`](@ref).
"""
function process_node(assemblage, tree, node; nsims=200)
    clade = get_clade(assemblage, tree, node)
    _isanalysable(assemblage, tree, node) || return (fill(NaN, nsites(clade)), NaN)
    occ = richness(clade) .> 0          # focal clade's occupied cells (deterministic)
    sims = simulate_descendants(clade, tree, first(getchildren(tree, node)); nsims)
    return sos(sims), gnd(sims, occ)
end

# Serial pre-pass shared by `node_analysis` and `node_metrics`: everything that touches
# the tree, which is not safe to share across threads. Returns the internal-node list
# and, per node, whether it is analysable plus the focal-clade and first-descendant
# species names, so the parallel heavy pass need only read the assemblage.
function _analysis_prepass(assemblage, tree)
    nodevec = _internalnodes(tree)
    N = length(nodevec)
    analysable = falses(N)
    parentsp = Vector{Vector{String}}(undef, N)   # focal clade species
    descsp = Vector{Vector{String}}(undef, N)   # first descendant's species
    for (i, node) in enumerate(nodevec)
        _isanalysable(assemblage, tree, node) || continue
        analysable[i] = true
        parentsp[i] = nodespecies(tree, node)
        descsp[i] = nodespecies(tree, first(getchildren(tree, node)))
    end
    return nodevec, analysable, parentsp, descsp
end

"""
    node_analysis(assemblage, tree; nsims=200) -> NodeAnalysis

The per-cell [`sos`](@ref) and the [`gnd`](@ref) of every internal node of `tree`, each
from `nsims` draws of the null model. [`node_metrics`](@ref) also gives the effect-size
scores from the same draws, and is usually the better choice.

The result is meant to be computed once, cached (e.g. with JLD2) and explored with
[`divergent_nodes`](@ref), [`sos_distances`](@ref), `plot_gnd` and `plot_node`. The
nodes are analysed in parallel: start Julia with several threads (`julia -t auto`).
"""
function node_analysis(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)
    gndv = fill(NaN, N)
    sosv = Vector{Vector{Float64}}(undef, N)
    done = Threads.Atomic{Int}(0)
    total = count(analysable)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage; species=parentsp[i])
        occ = richness(clade) .> 0
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        gndv[i] = gnd(sims, occ)
        sosv[i] = sos(sims)
        _progress!(done, total, "node_analysis")
    end
    gnd_scores = Dict{String,Float64}()
    sos_scores = Dict{String,Vector{Float64}}()
    for i in 1:N
        gnd_scores[nodevec[i]] = gndv[i]
        analysable[i] && (sos_scores[nodevec[i]] = sosv[i])
    end
    return NodeAnalysis(nodevec, gnd_scores, sos_scores)
end

"""
    node_metrics(assemblage, tree; nsims=200) -> NodeMetrics

The divergence of every internal node of `tree`: its per-cell [`sos`](@ref), the
original [`gnd`](@ref), and the effect-size scores [`sos_rms`](@ref) (`rms`, the
recommended score), [`sos_sd`](@ref) (`sd`), [`divergence_ses`](@ref) (`ses`) and
[`divergence_pval`](@ref) (`pval`), all from the same `nsims` draws of the null model.
`varying` is the share of each node's occupied cells where the null model varies, the
cells the SOS-based scores rest on; `sqrt(varying) * rms` is the RMS-SOS over all
occupied cells, with the others counted as 0.

The null model is the curveball (swap) randomisation of the parent clade's
presence-absence matrix, which keeps each species' range size and each cell's richness.
So the scores measure where the two descendant clades occur, not how rich they are.

Nodes that cannot be analysed (see [`process_node`](@ref)) have a `NaN` GND and no
other entries. The result is meant to be computed once, cached (e.g. with JLD2) and
explored with [`divergent_nodes`](@ref), [`sos_distances`](@ref), `plot_gnd` and
`plot_node`.

The nodes are analysed in parallel: start Julia with several threads (`julia -t auto`)
for a near-linear speed-up. The results carry Monte Carlo noise and differ between runs.
"""
function node_metrics(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)

    # parallel heavy pass: only assemblage reads + thread-local randomisers
    gndv = fill(NaN, N)
    rmsv = fill(NaN, N)
    sdv = fill(NaN, N)
    sesv = fill(NaN, N)
    pvalv = fill(NaN, N)
    varyv = fill(NaN, N)
    sosv = Vector{Vector{Float64}}(undef, N)
    done = Threads.Atomic{Int}(0)
    total = count(analysable)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage; species=parentsp[i])
        occ = richness(clade) .> 0                        # focal clade's occupied cells
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        me, sd = _moments(sims)
        sosv[i] = _sos(sims, me, sd)
        gndv[i] = gnd(sims, occ)
        rmsv[i] = sos_rms(sosv[i])
        sdv[i] = sos_sd(sosv[i])
        stat, nullstats = _msos_null(sims, me, sd, occ .& (sd .> 0))
        sesv[i] = _ses(stat, nullstats)
        pvalv[i] = _pval(stat, nullstats)
        varyv[i] = _varying_share(sosv[i], occ)
        _progress!(done, total, "node_metrics")
    end

    # assemble the result dicts (serial)
    gnd_scores = Dict{String,Float64}()
    rms = Dict{String,Float64}()
    sds = Dict{String,Float64}()
    ses = Dict{String,Float64}()
    pval = Dict{String,Float64}()
    varying = Dict{String,Float64}()
    sos_scores = Dict{String,Vector{Float64}}()
    for i in 1:N
        gnd_scores[nodevec[i]] = gndv[i]
        analysable[i] || continue
        rms[nodevec[i]] = rmsv[i]
        sds[nodevec[i]] = sdv[i]
        ses[nodevec[i]] = sesv[i]
        pval[nodevec[i]] = pvalv[i]
        varying[nodevec[i]] = varyv[i]
        sos_scores[nodevec[i]] = sosv[i]
    end
    return NodeMetrics(nodevec, gnd_scores, rms, sds, ses, pval, varying, sos_scores)
end

"""
    divergent_nodes(res::NodeMetrics; by=:rms, threshold) -> Vector{String}
    divergent_nodes(res::NodeAnalysis; threshold=0.8) -> Vector{String}
    divergent_nodes(scores::AbstractDict; threshold=0.8) -> Vector{String}

The divergent nodes of an analysis result, in the order they appear on the tree.

For a `NodeMetrics`, `by` picks the score: `:rms` (RMS-SOS, the default; divergent above
`threshold = 1.5`), `:pval` (below `threshold = 0.05`; being a significance, larger clades
pass it more easily) or `:gnd` (above `threshold = 0.8`). A `NodeAnalysis` only has the
GND. For a plain Dict of node => score, the nodes above `threshold` come most divergent
first, as there is no tree order. Nodes with a `NaN` score are never divergent.
"""
function divergent_nodes(res::NodeMetrics; by=:rms, threshold=_default_threshold(by))
    if by == :rms
        return _intreeorder(res, res.rms, >(threshold))
    elseif by == :pval
        return _intreeorder(res, res.pval, <(threshold))
    elseif by == :gnd
        return _intreeorder(res, res.gnd, >(threshold))
    else
        throw(ArgumentError("`by` must be :rms, :pval or :gnd; got $(repr(by))"))
    end
end
function divergent_nodes(res::NodeAnalysis; threshold=0.8)
    return _intreeorder(res, res.gnd, >(threshold))
end
function divergent_nodes(scores::AbstractDict; threshold=0.8)
    nodes = [node for (node, score) in scores if !isnan(score) && score > threshold]
    return sort!(nodes; by=node -> (-scores[node], node))
end

function _default_threshold(by)
    by == :rms && return 1.5
    by == :pval && return 0.05
    return 0.8
end

# The nodes of `res`, in tree order, whose score passes `isdivergent` (NaN scores excluded)
function _intreeorder(res, scores, isdivergent)
    passes(n) = haskey(scores, n) && !isnan(scores[n]) && isdivergent(scores[n])
    return filter(passes, res.nodes)
end
