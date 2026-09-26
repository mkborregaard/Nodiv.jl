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

# Calculate SOS and GND for a single node (NaN when the node can't be analysed).
function process_node(assemblage, tree, nodename; nsims=200)
    clade = get_clade(assemblage, tree, nodename)
    _isanalysable(assemblage, tree, nodename) || return (fill(NaN, nsites(clade)), NaN)
    occ = richness(clade) .> 0          # focal clade's occupied cells (deterministic)
    sims = simulate_descendants(clade, tree, first(getchildren(tree, nodename)); nsims)
    return calculate_SOS(sims), calculate_GND(sims, occ)
end

# Run the node-based analysis over every internal node of the tree.
# Recreates the main `Node_analysis` function of the nodiv R package
# (https://github.com/mkborregaard/nodiv).
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec = _internalnodes(tree)
    SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
    GNDs = Vector{Float64}(undef, length(nodevec))
    @progress for (i, node) in enumerate(nodevec)
        SOSs[:, i], GNDs[i] = process_node(assemblage, tree, node; nsims)
    end
    return SOSs, GNDs
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

# Compute the per-cell SOS pattern and the GND for every internal node. Returns a
# `NodeAnalysis` to hand to the explore functions - the "compute once, explore a lot"
# entry point; cache it (e.g. with JLD2) and reload it. See `node_metrics` for the
# effect-size scores (RMS/spatial/SES/pval). Multithreaded like `node_metrics`: launch
# Julia with `-t auto` for the speedup.
function node_analysis(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)
    gndv = fill(NaN, N)
    sosv = Vector{Vector{Float64}}(undef, N)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage; species=parentsp[i])
        occ = richness(clade) .> 0
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        gndv[i] = calculate_GND(sims, occ)
        sosv[i] = calculate_SOS(sims)
    end
    gnd = Dict{String,Float64}()
    sos = Dict{String,Vector{Float64}}()
    for i in 1:N
        gnd[nodevec[i]] = gndv[i]
        analysable[i] && (sos[nodevec[i]] = sosv[i])
    end
    return NodeAnalysis(nodevec, gnd, sos)
end

# One-pass analysis returning GND together with the effect-size alternatives, all
# from the same null draws (one randomisation per node). The published swap null
# fixes range size, so the `spatial`/GND scores measure spatial turnover rather than
# richness asymmetry. Cache the result with JLD2 as for `node_analysis`.
#
# Multithreaded: the per-node swap randomisation dominates the cost and is independent
# across nodes, so start Julia with `-t auto` (or set JULIA_NUM_THREADS) for a near-
# linear speedup. All tree access is done in a serial pre-pass; the parallel section
# only reads the assemblage and builds thread-local randomisers. Results carry the
# usual Monte-Carlo noise and are not bit-reproducible across runs (true single-thread
# too, as nothing is seeded).
function node_metrics(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)

    # parallel heavy pass: only assemblage reads + thread-local randomisers
    gndv = fill(NaN, N)
    rmsv = fill(NaN, N)
    spatv = fill(NaN, N)
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
        gndv[i] = calculate_GND(sims, occ)
        rmsv[i] = gnd_rms(sosv[i])
        spatv[i] = gnd_spatial(sosv[i])
        stat, nullstats = _msos_null(sims, me, sd, occ .& (sd .> 0))
        sesv[i] = _ses(stat, nullstats)
        pvalv[i] = _pval(stat, nullstats)
        varyv[i] = _varying_share(sosv[i], occ)
        n = Threads.atomic_add!(done, 1) + 1
        n % 100 == 0 && @info "node_metrics: $n / $total analysable nodes done"
    end

    # assemble the result dicts (serial)
    gnd = Dict{String,Float64}()
    rms = Dict{String,Float64}()
    spatial = Dict{String,Float64}()
    ses = Dict{String,Float64}()
    pval = Dict{String,Float64}()
    varying = Dict{String,Float64}()
    sos = Dict{String,Vector{Float64}}()
    for i in 1:N
        gnd[nodevec[i]] = gndv[i]
        analysable[i] || continue
        rms[nodevec[i]] = rmsv[i]
        spatial[nodevec[i]] = spatv[i]
        ses[nodevec[i]] = sesv[i]
        pval[nodevec[i]] = pvalv[i]
        varying[nodevec[i]] = varyv[i]
        sos[nodevec[i]] = sosv[i]
    end
    return NodeMetrics(nodevec, gnd, rms, spatial, ses, pval, varying, sos)
end

# Node names whose GND exceeds `threshold` (NaN GNDs excluded). For a `NodeAnalysis` the
# nodes come in tree order; a plain Dict of scores has no tree, so its nodes come most
# divergent first (ties by name).
function divergent_nodes(scores::AbstractDict; threshold=0.8)
    nodes = [node for (node, g) in scores if !isnan(g) && g > threshold]
    return sort!(nodes; by=node -> (-scores[node], node))
end
function divergent_nodes(res::NodeAnalysis; threshold=0.8)
    return _intreeorder(res, res.gnd, >(threshold))
end

# The nodes of `res`, in tree order, whose score passes `isdivergent` (NaN scores excluded)
function _intreeorder(res, scores, isdivergent)
    passes(n) = haskey(scores, n) && !isnan(scores[n]) && isdivergent(scores[n])
    return filter(passes, res.nodes)
end

function _default_threshold(by)
    by == :rms && return 1.5
    by == :pval && return 0.05
    return 0.8
end

# For a `NodeMetrics`, rank divergence by the size-robust RMS-SOS effect size by
# default (null = 1). Pass `by = :pval` for the null-calibrated Monte-Carlo
# significance (note: significance carries a clade-size/power bias), or `by = :gnd`
# for the original GND. The threshold default adapts to the chosen score.
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
