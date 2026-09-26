# Calculate SOS and GND for a single node (NaN when the node can't be analysed).
function process_node(assemblage, tree, nodename; nsims=200)
    clade = get_clade(assemblage, tree, nodename)
    children = getchildren(tree, nodename)

    if length(children) != 2 ||
        any(x -> isleaf(tree, x) || nspecies(get_clade(assemblage, tree, x)) < 3, children)
        return (fill(NaN, nsites(clade)), NaN)
    end

    occ = richness(clade) .> 0          # focal clade's occupied cells (deterministic)
    sims = simulate_descendants(clade, tree, children[1]; nsims)
    return calculate_SOS(sims), calculate_GND(sims, occ)
end

# Run the node-based analysis over every internal node of the tree.
# Recreates the main `Node_analysis` function of the nodiv R package
# (https://github.com/mkborregaard/nodiv).
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
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
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    N = length(nodevec)
    analysable = falses(N)
    parentsp = Vector{Vector{String}}(undef, N)   # focal clade species
    descsp = Vector{Vector{String}}(undef, N)   # first descendant's species
    for (i, node) in enumerate(nodevec)
        ch = getchildren(tree, node)
        if length(ch) == 2 &&
            all(x -> !isleaf(tree, x) && nspecies(get_clade(assemblage, tree, x)) >= 3, ch)
            analysable[i] = true
            parentsp[i] = nodespecies(tree, node)
            descsp[i] = nodespecies(tree, ch[1])
        end
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
    return NodeAnalysis(gnd, sos)
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
    sosv = Vector{Vector{Float64}}(undef, N)
    done = Threads.Atomic{Int}(0)
    total = count(analysable)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage; species=parentsp[i])
        occ = richness(clade) .> 0                        # focal clade's occupied cells
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        gndv[i] = calculate_GND(sims, occ)
        rmsv[i] = calculate_GND_rms(sims, occ)
        spatv[i] = calculate_GND_spatial(sims, occ)
        sesv[i] = calculate_GND_ses(sims, occ)
        pvalv[i] = calculate_GND_pval(sims, occ)
        sosv[i] = calculate_SOS(sims)
        n = Threads.atomic_add!(done, 1) + 1
        n % 100 == 0 && @info "node_metrics: $n / $total analysable nodes done"
    end

    # assemble the result dicts (serial)
    gnd = Dict{String,Float64}()
    rms = Dict{String,Float64}()
    spatial = Dict{String,Float64}()
    ses = Dict{String,Float64}()
    pval = Dict{String,Float64}()
    sos = Dict{String,Vector{Float64}}()
    for i in 1:N
        gnd[nodevec[i]] = gndv[i]
        analysable[i] || continue
        rms[nodevec[i]] = rmsv[i]
        spatial[nodevec[i]] = spatv[i]
        ses[nodevec[i]] = sesv[i]
        pval[nodevec[i]] = pvalv[i]
        sos[nodevec[i]] = sosv[i]
    end
    return NodeMetrics(gnd, rms, spatial, ses, pval, sos)
end

# Node names whose GND exceeds `threshold` (NaN GNDs excluded). Accepts a `NodeAnalysis`
# (from `node_analysis`) or a plain GND Dict.
function divergent_nodes(gnd::AbstractDict; threshold=0.8)
    return [node for (node, g) in gnd if !isnan(g) && g > threshold]
end
divergent_nodes(res::NodeAnalysis; threshold=0.8) = divergent_nodes(res.gnd; threshold)

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
        [n for (n, v) in res.rms if !isnan(v) && v > threshold]
    elseif by == :pval
        [n for (n, p) in res.pval if !isnan(p) && p < threshold]
    elseif by == :gnd
        divergent_nodes(res.gnd; threshold)
    else
        error("`by` must be :rms, :pval or :gnd")
    end
end
