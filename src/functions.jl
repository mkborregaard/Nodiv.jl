# Functions for node-based analysis of species distributions.
#
# Extracted from the SpatialEcology.jl docs:
# https://docs.ecojulia.org/SpatialEcology.jl/stable/examples/nodebased/
#
# Reimplements the method of Borregaard et al. (2014), Node-based analysis of
# species distributions, Methods in Ecology and Evolution 5: 1225-1235.

using SpatialEcology
using Phylo
using RecipesBase
using Random: rand!
using Statistics
using StatsBase: tiedrank
using ProgressLogging

# All tips/species descending from a node (or the node itself if it is a tip).
# Returns a concretely-typed Vector{String}: `getdescendants` is typed as
# Union{Nothing,String}, and an empty/tip result of that type misses the name
# matching method used by the downstream `view`.
function nodespecies(tree, node)
    isleaf(tree, node) && return String[getnodename(tree, node)]
    String[x for x in getdescendants(tree, node) if isleaf(tree, x)]
end

# Subset an Assemblage to the clade descending from a node.
get_clade(assemblage, tree, node) = view(assemblage, species = nodespecies(tree, node))

# Result of `node_analysis`: `gnd` maps every node name to its GND (NaN where the
# node cannot be analysed); `sos` maps the analysable nodes to their per-cell SOS
# vector. Pass it straight to the explore functions - `divergent_nodes`,
# `plot_gnd`, `plot_node`, `sos_distances` - to reuse cached results.
struct NodeAnalysis
    gnd::Dict{String, Float64}
    sos::Dict{String, Vector{Float64}}
end

# Plot a node in a 2x2 grid: the parent clade (top-left), the descendant SOS
# mapped over the assemblage's coordinates (top-right, RdYlBu), and the two child
# clades on the bottom row. Defined as a plot recipe so the package depends only
# on RecipesBase; the plotting backend (Plots) is supplied by the caller. Use as
# `plot_node(assemblage, tree, node)`.
@userplot Plot_Node

@recipe function f(pn::Plot_Node)
    assemblage, tree, node = pn.args[1:3]
    ch1, ch2 = getchildren(tree, node)[1:2]
    assm = get_clade(assemblage, tree, node)
    assmch1 = get_clade(assm, tree, ch1)
    assmch2 = get_clade(assm, tree, ch2)

    # SOS for the top-right panel, taken from the cached analysis result supplied as
    # the 4th argument - either a `NodeAnalysis`/`NodeMetrics` (looked up by node) or a
    # precomputed SOS vector. Pass the result of `node_metrics`/`node_analysis`; the
    # panel is never recomputed on the fly.
    length(pn.args) >= 4 ||
        error("plot_node needs the analysis result (or a precomputed SOS vector) as the " *
              "4th argument, e.g. plot_node(assemblage, tree, node, res)")
    cached = pn.args[4]
    sos = cached isa Union{NodeAnalysis, NodeMetrics} ? cached.sos[node] : cached

    layout := (2, 2)
    size --> (900, 800)

    @series begin              # top-left: parent clade
        subplot := 1
        title := string(node)
        assm
    end
    @series begin              # top-right: SOS in environmental space
        subplot := 2
        title := "SOS"
        fillcolor := :RdYlBu
        clim := (-8, 8)
        sos, assm
    end
    @series begin              # bottom-left: child 1
        subplot := 3
        title := string(ch1)
        assmch1
    end
    @series begin              # bottom-right: child 2
        subplot := 4
        title := string(ch2)
        assmch2
    end
end

# Plot per-node values (e.g. GND) on the tree. `gndvals` is a Dict of node name
# => value; a coloured marker is drawn at every node present in it (NaN values
# skipped) and nothing elsewhere. Pass a NodeAnalysis (or its `gnd` Dict) to show
# all analysable nodes, or a filtered Dict (e.g. only divergent nodes) for a subset.
# Use as `plot_gnd(tree, gnd)`.
@userplot Plot_Gnd

@recipe function f(pg::Plot_Gnd)
    tree, gndvals = pg.args
    gndvals isa Union{NodeAnalysis, NodeMetrics} && (gndvals = gndvals.gnd)
    shown = Dict(k => v for (k, v) in gndvals if !isnan(v))

    # The tree recipe draws a marker at every node and renders NaN-marker_z nodes
    # solid, so hide the non-shown nodes with size 0. markersize must follow the
    # recipe's node order, which is Phylo's layout helper _findxy.
    layoutnodes = Phylo._findxy(tree)[3]

    treetype --> :fan
    showtips --> false
    markerstrokewidth --> 0
    color --> :YlOrRd
    clim --> (0, 1)
    size --> (1000, 1000)
    marker_z := shown
    markersize := [haskey(shown, node) ? 6 : 0 for node in layoutnodes]
    tree
end

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

# Standardized effect size (SOS metric) per grid cell from a simulation matrix.
function calculate_SOS(sims)
    sd = std.(eachcol(sims))
    me = mean.(eachcol(sims))
    (sims[1, :] .- me) ./ sd
end

# Default occupancy mask when the parent's occupancy is not supplied: cells whose
# descendant-richness column varies across the draws. The analysis entry points pass
# the deterministic, nsims-independent "parent clade present" mask instead.
_occupied(sims) = [!all(==(first(c)), c) for c in eachcol(sims)]

# Geographic node divergence (GND metric) from a simulation matrix, averaged over the
# OCCUPIED sites of the focal (parent) clade (Borregaard et al. 2014). Cells where the
# parent is absent carry no divergence signal and are excluded; a constant occupied
# column (descendant richness fixed across draws) contributes P ~ 1, i.e. no
# divergence. `occupied` is a per-cell boolean mask (default: the non-constant columns).
function calculate_GND(sims, occupied = _occupied(sims))
  # two internal convenience functions
  logit(p) = log(p/(1-p))
  invlogit(p) = exp(p)/(1+exp(p))

  n = size(sims, 1)
  idx = findall(occupied)
  isempty(idx) && return NaN
  r = [tiedrank(view(sims, :, j))[1]/(n + 1) for j in idx]
  # two-sided P (eqn 3); the -1/n keeps P off the 0/1 boundary so logit stays finite
  p = 1 .- 2 .* abs.(r .- 0.5) .- 1/n
  α = mean(logit.(p))
  1-invlogit(α)
end

# ---- effect-size alternatives to GND -----------------------------------------
# GND (eqn 4) combines per-cell two-sided P values in logit space. That makes it a
# p-value summary, not an effect size: its scale rides on `nsims` (P is bounded by
# ~1/nsims, so the maximum GND is 1 - O(1/nsims)), its no-divergence baseline is ~0.5
# rather than 0, and it saturates - the most divergent nodes pile against the ceiling.
# The functions below summarise the SOS field directly instead, giving a
# replication-stable magnitude in units of null SDs.

# Total spatial-divergence intensity: the root-mean-square SOS over the focal clade's
# occupied cells. Each SOS is ~standardised, so this is ~1 under the null and >1 under
# divergence. Replication-stable and not dominated by boundary cells. Constant occupied
# cells (sd = 0) contribute 0. Under the swap null, which fixes range size, a uniform
# richness/occupancy asymmetry between the clades is not flagged (use
# `calculate_GND_spatial` if you ever need to strip an offset explicitly).
function calculate_GND_rms(sims, occupied = _occupied(sims))
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    idx = findall(occupied)
    isempty(idx) && return NaN
    sqrt(mean(j -> isfinite(sos[j]) ? sos[j]^2 : 0.0, idx))
end

# Spatial-only intensity: as `calculate_GND_rms` but with the uniform offset
# (mean SOS = the clades' richness/occupancy asymmetry) removed, so it reflects only
# how over/under-representation varies ACROSS cells. Equals the SD of the SOS field.
function calculate_GND_spatial(sims, occupied = _occupied(sims))
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    idx = findall(occupied)
    isempty(idx) && return NaN
    std([isfinite(sos[j]) ? sos[j] : 0.0 for j in idx])
end

# Mean-square-SOS statistic of a richness row, standardised against the per-cell null
# moments, over occupied cells that actually vary. Internal helper for the SES and the
# Monte-Carlo P value, which share the null distribution of this statistic.
function _msos_null(sims, occupied)
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    keep = occupied .& (sd .> 0)
    any(keep) || return (NaN, Float64[])
    mek = me[keep]; sdk = sd[keep]
    msos(row) = (z = (view(row, keep) .- mek) ./ sdk; mean(abs2, z))
    msos(view(sims, 1, :)), [msos(view(sims, i, :)) for i in 2:size(sims, 1)]
end

# Standardised effect size of the divergence: the mean-square-SOS statistic of the
# empirical row, expressed in SDs of its null distribution. 0 = no divergence;
# replication-stable. NB the SES magnitude inflates with clade size (more cells ->
# tighter null), so prefer RMS-SOS for cross-node / cross-grain comparison.
function calculate_GND_ses(sims, occupied = _occupied(sims))
    Temp, Tnull = _msos_null(sims, occupied)
    isempty(Tnull) && return NaN
    s = std(Tnull)
    s == 0 ? NaN : (Temp - mean(Tnull)) / s
end

# Null-calibrated significance: the one-sided Monte-Carlo P value of the mean-square-
# SOS statistic (add-one estimator, bounded by 1/(nsims+1)). Small = divergent. NB this
# is a significance, so it carries a clade-size/power bias - larger clades clear a fixed
# P threshold more easily; `calculate_GND_rms` does not.
function calculate_GND_pval(sims, occupied = _occupied(sims))
    Temp, Tnull = _msos_null(sims, occupied)
    isempty(Tnull) && return NaN
    (1 + count(>=(Temp), Tnull)) / (length(Tnull) + 1)
end

# Same summaries from an already-computed per-cell SOS vector (e.g. a cached
# `NodeAnalysis.sos[node]`), so they can be derived without re-running the null.
# NaN/Inf cells (clade absent / constant) are dropped. `gnd_rms` = total intensity,
# `gnd_spatial` = spatial-only. The SES/P value need the null draws, so have no SOS-only form.
gnd_rms(sos::AbstractVector)     = (f = filter(isfinite, sos); isempty(f) ? NaN : sqrt(mean(abs2, f)))
gnd_spatial(sos::AbstractVector) = (f = filter(isfinite, sos); isempty(f) ? NaN : std(f))

# Calculate SOS and GND for a single node (NaN when the node can't be analysed).
function process_node(assemblage, tree, nodename; nsims = 200)
    clade = get_clade(assemblage, tree, nodename)
    children = getchildren(tree, nodename)

    if length(children) != 2 || any(x -> isleaf(tree, x) || nspecies(get_clade(assemblage, tree, x)) < 3, children)
        return (fill(NaN, nsites(clade)), NaN)
    end

    occ = richness(clade) .> 0          # focal clade's occupied cells (deterministic)
    sims = simulate_descendants(clade, tree, children[1]; nsims)
    calculate_SOS(sims), calculate_GND(sims, occ)
end

# Run the node-based analysis over every internal node of the tree.
# Recreates the main `Node_analysis` function of the nodiv R package
# (https://github.com/mkborregaard/nodiv).
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 200)
   nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)] #shuffle!(collect(nodenamefilter(!isleaf, tree)))
   SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
   GNDs = Vector{Float64}(undef, length(nodevec))
   @progress for (i, node) in enumerate(nodevec)
       SOSs[:,i], GNDs[i] = process_node(assemblage, tree, node; nsims)
   end
   SOSs, GNDs
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
    descsp   = Vector{Vector{String}}(undef, N)   # first descendant's species
    for (i, node) in enumerate(nodevec)
        ch = getchildren(tree, node)
        if length(ch) == 2 && all(x -> !isleaf(tree, x) && nspecies(get_clade(assemblage, tree, x)) >= 3, ch)
            analysable[i] = true
            parentsp[i] = nodespecies(tree, node)
            descsp[i]   = nodespecies(tree, ch[1])
        end
    end
    nodevec, analysable, parentsp, descsp
end

# Compute the per-cell SOS pattern and the GND for every internal node. Returns a
# `NodeAnalysis` to hand to the explore functions - the "compute once, explore a lot"
# entry point; cache it (e.g. with JLD2) and reload it. See `node_metrics` for the
# effect-size scores (RMS/spatial/SES/pval). Multithreaded like `node_metrics`: launch
# Julia with `-t auto` for the speedup.
function node_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)
    gndv = fill(NaN, N)
    sosv = Vector{Vector{Float64}}(undef, N)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage, species = parentsp[i])
        occ = richness(clade) .> 0
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        gndv[i] = calculate_GND(sims, occ)
        sosv[i] = calculate_SOS(sims)
    end
    gnd = Dict{String, Float64}(); sos = Dict{String, Vector{Float64}}()
    for i in 1:N
        gnd[nodevec[i]] = gndv[i]
        analysable[i] && (sos[nodevec[i]] = sosv[i])
    end
    NodeAnalysis(gnd, sos)
end

# Like `NodeAnalysis`, but also carries the effect-size scores: `rms` (total
# intensity), `spatial` (spatial-only), `ses` (standardised effect size) and `pval`
# (null-calibrated Monte-Carlo significance). Kept separate from `NodeAnalysis` so
# existing cached results still load unchanged.
struct NodeMetrics
    gnd::Dict{String, Float64}
    rms::Dict{String, Float64}
    spatial::Dict{String, Float64}
    ses::Dict{String, Float64}
    pval::Dict{String, Float64}
    sos::Dict{String, Vector{Float64}}
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
function node_metrics(assemblage::Assemblage, tree::AbstractTree; nsims = 200)
    nodevec, analysable, parentsp, descsp = _analysis_prepass(assemblage, tree)
    N = length(nodevec)

    # parallel heavy pass: only assemblage reads + thread-local randomisers
    gndv = fill(NaN, N); rmsv = fill(NaN, N); spatv = fill(NaN, N)
    sesv = fill(NaN, N); pvalv = fill(NaN, N)
    sosv = Vector{Vector{Float64}}(undef, N)
    done = Threads.Atomic{Int}(0); total = count(analysable)
    Threads.@threads :dynamic for i in 1:N
        analysable[i] || continue
        clade = view(assemblage, species = parentsp[i])
        occ = richness(clade) .> 0                        # focal clade's occupied cells
        sims = _simulate_descendants(clade, descsp[i]; nsims)
        gndv[i]  = calculate_GND(sims, occ);     rmsv[i]  = calculate_GND_rms(sims, occ)
        spatv[i] = calculate_GND_spatial(sims, occ); sesv[i] = calculate_GND_ses(sims, occ)
        pvalv[i] = calculate_GND_pval(sims, occ);    sosv[i] = calculate_SOS(sims)
        n = Threads.atomic_add!(done, 1) + 1
        n % 100 == 0 && @info "node_metrics: $n / $total analysable nodes done"
    end

    # assemble the result dicts (serial)
    gnd = Dict{String, Float64}(); rms = Dict{String, Float64}()
    spatial = Dict{String, Float64}(); ses = Dict{String, Float64}()
    pval = Dict{String, Float64}(); sos = Dict{String, Vector{Float64}}()
    for i in 1:N
        gnd[nodevec[i]] = gndv[i]
        analysable[i] || continue
        rms[nodevec[i]] = rmsv[i]; spatial[nodevec[i]] = spatv[i]; ses[nodevec[i]] = sesv[i]
        pval[nodevec[i]] = pvalv[i]; sos[nodevec[i]] = sosv[i]
    end
    NodeMetrics(gnd, rms, spatial, ses, pval, sos)
end

# Prune `tree` in place to the tips it shares with all the given assemblage(s),
# so clade subsetting never references a species absent from the data. Returns
# the tree.
function prune_to_shared!(tree, assemblages...)
    shared = intersect(getleafnames(tree), speciesnames.(assemblages)...)
    keeptips!(tree, shared)
    tree
end

# Node names whose GND exceeds `threshold` (NaN GNDs excluded). Accepts a `NodeAnalysis`
# (from `node_analysis`) or a plain GND Dict.
divergent_nodes(gnd::AbstractDict; threshold = 0.8) =
    [node for (node, g) in gnd if !isnan(g) && g > threshold]
divergent_nodes(res::NodeAnalysis; threshold = 0.8) = divergent_nodes(res.gnd; threshold)

# For a `NodeMetrics`, rank divergence by the size-robust RMS-SOS effect size by
# default (null = 1). Pass `by = :pval` for the null-calibrated Monte-Carlo
# significance (note: significance carries a clade-size/power bias), or `by = :gnd`
# for the original GND. The threshold default adapts to the chosen score.
function divergent_nodes(res::NodeMetrics; by = :rms,
                         threshold = by == :rms ? 1.5 : by == :pval ? 0.05 : 0.8)
    if by == :rms
        [n for (n, v) in res.rms  if !isnan(v) && v > threshold]
    elseif by == :pval
        [n for (n, p) in res.pval if !isnan(p) && p < threshold]
    elseif by == :gnd
        divergent_nodes(res.gnd; threshold)
    else
        error("`by` must be :rms, :pval or :gnd")
    end
end

# Size-corrected divergence: RMS-SOS residualised on log clade richness (the clade-size
# signal is mostly a species-count effect). Positive residual = more divergent than a
# clade of that size typically is. Returns a Dict node => residual over analysable nodes.
function size_residual(res::NodeMetrics, tree)
    nodes = [n for (n, v) in res.rms if !isnan(v)]
    y = [res.rms[n] for n in nodes]
    x = [log(length(nodespecies(tree, n))) for n in nodes]
    X = hcat(ones(length(x)), x)
    resid = y .- X * (X \ y)
    Dict(nodes[i] => resid[i] for i in eachindex(nodes))
end

# Pairwise distance matrix (1 - |Pearson r|) between per-cell SOS patterns,
# correlated over cells where both patterns are defined (SOS is NaN where a clade
# is absent). Feed to an ordination (e.g. MDS) to see which nodes have similar SOS
# maps. Either compute the SOS on the fly from an assemblage + nodes, or pass
# precomputed SOS vectors (e.g. cached from `node_analysis`).
function _sos_distances(sosmat)
    n = size(sosmat, 2)
    D = zeros(n, n)
    for i in 1:n, j in i+1:n
        a, b = view(sosmat, :, i), view(sosmat, :, j)
        ok = .!(isnan.(a) .| isnan.(b))
        r = count(ok) > 2 ? cor(a[ok], b[ok]) : 0.0
        D[i, j] = D[j, i] = 1 - abs(isnan(r) ? 0.0 : r)
    end
    D
end
sos_distances(sosvectors::AbstractVector) = _sos_distances(reduce(hcat, sosvectors))
sos_distances(res::NodeAnalysis, nodes) = sos_distances([res.sos[n] for n in nodes])
sos_distances(res::NodeMetrics, nodes) = sos_distances([res.sos[n] for n in nodes])
sos_distances(assemblage, tree, nodes; nsims = 200) =
    _sos_distances(reduce(hcat, process_node(assemblage, tree, node; nsims)[1] for node in nodes))
