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

    # SOS for the top-right panel. A 4th argument supplies a cached SOS - either a
    # `NodeAnalysis` (looked up by node) or a precomputed SOS vector; otherwise it
    # is recomputed here (a fresh randomization, so the panel varies between calls).
    sos = if length(pn.args) >= 4
        cached = pn.args[4]
        cached isa Union{NodeAnalysis, NodeMetrics} ? cached.sos[node] : cached
    else
        calculate_SOS(simulate_descendants(assm, tree, ch1; nsims = 1000))
    end

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
function simulate_descendants(clade, tree, descendant; nsims = 99)
    ret = zeros(nsims + 1, nsites(clade))  # a matrix to hold the richness values from the simulations
    ret[1, :] = richness(get_clade(clade, tree, descendant))   # empirical richness in the first row
    rmg = matrixrandomizer(clade)
    for i in 2:nsims + 1
        # simulated richness in the remaining rows
        ret[i, :] .= richness(get_clade(rand!(rmg), tree, descendant))
    end
    ret
end

# Standardized effect size (SOS metric) per grid cell from a simulation matrix.
function calculate_SOS(sims)
    sd = std.(eachcol(sims))
    me = mean.(eachcol(sims))
    (sims[1, :] .- me) ./ sd
end

# Geographic node divergence (GND metric) from a simulation matrix. Summarised over
# OCCUPIED sites only (Borregaard et al. 2014): a cell where the clade is absent has
# identical richness in every simulation (a constant column) and carries no
# divergence signal. Including such cells dilutes GND toward 0 by a factor that
# scales with the fraction of empty cells, so it deflates GND on finer grids - drop
# the constant columns before averaging.
function calculate_GND(sims)
  # two internal convenience functions
  logit(p) = log(p/(1-p))
  invlogit(p) = exp(p)/(1+exp(p))

  n = size(sims, 1)
  occupied = Iterators.filter(x -> !all(==(first(x)), x), eachcol(sims))
  r = [tiedrank(x)[1]/(n + 1) for x in occupied]
  isempty(r) && return NaN
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

# Total spatial-divergence intensity: the root-mean-square SOS over occupied
# (non-constant) cells. Each SOS is ~standardised, so this is ~1 under the null and
# >1 under divergence. Replication-stable and not dominated by boundary cells. NB:
# under a null that does not fix per-clade range size (e.g. :tipshuffle), this also
# responds to a uniform richness/occupancy asymmetry between the clades; use
# `calculate_GND_spatial` to strip that out, or the :swap null which fixes range size.
function calculate_GND_rms(sims)
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    keep = isfinite.(sos)                 # drops constant columns (sd = 0 -> NaN/Inf)
    any(keep) ? sqrt(mean(abs2, sos[keep])) : NaN
end

# Spatial-only intensity: as `calculate_GND_rms` but with the uniform offset
# (mean SOS = the clades' richness/occupancy asymmetry) removed, so it reflects only
# how over/under-representation varies ACROSS cells. Equals the SD of the SOS field.
function calculate_GND_spatial(sims)
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    keep = isfinite.(sos)
    any(keep) ? std(sos[keep]) : NaN
end

# Standardised effect size of the divergence: the mean-square-SOS statistic of the
# empirical row, expressed in SDs of its null distribution (the simulated rows scored
# against the same per-cell moments). 0 = no divergence; replication-stable scale.
function calculate_GND_ses(sims)
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    keep = sd .> 0
    any(keep) || return NaN
    mek = me[keep]; sdk = sd[keep]
    msos(row) = (z = (row[keep] .- mek) ./ sdk; mean(abs2, z))
    Temp = msos(view(sims, 1, :))
    Tnull = [msos(view(sims, i, :)) for i in 2:size(sims, 1)]
    s = std(Tnull)
    s == 0 ? NaN : (Temp - mean(Tnull)) / s
end

# Same summaries from an already-computed per-cell SOS vector (e.g. a cached
# `NodeAnalysis.sos[node]`), so they can be derived without re-running the null.
# NaN/Inf cells (clade absent / constant) are dropped. `gnd_rms` = total intensity,
# `gnd_spatial` = spatial-only. The SES needs the null draws, so it has no SOS-only form.
gnd_rms(sos::AbstractVector)     = (f = filter(isfinite, sos); isempty(f) ? NaN : sqrt(mean(abs2, f)))
gnd_spatial(sos::AbstractVector) = (f = filter(isfinite, sos); isempty(f) ? NaN : std(f))

# Calculate SOS and GND for a single node (NaN when the node can't be analysed).
function process_node(assemblage, tree, nodename; nsims = 100)
    clade = get_clade(assemblage, tree, nodename)
    children = getchildren(tree, nodename)

    if length(children) != 2 || any(x -> isleaf(tree, x) || nspecies(get_clade(assemblage, tree, x)) < 4, children)
        return (fill(NaN, nsites(clade)), NaN)
    end

    sims = simulate_descendants(clade, tree, children[1]; nsims)
    calculate_SOS(sims), calculate_GND(sims)
end

# Run the node-based analysis over every internal node of the tree.
# Recreates the main `Node_analysis` function of the nodiv R package
# (https://github.com/mkborregaard/nodiv).
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 100)
   nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)] #shuffle!(collect(nodenamefilter(!isleaf, tree)))
   SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
   GNDs = Vector{Float64}(undef, length(nodevec))
   @progress for (i, node) in enumerate(nodevec)
       SOSs[:,i], GNDs[i] = process_node(assemblage, tree, node; nsims)
   end
   SOSs, GNDs
end

# Compute both the per-cell SOS pattern and the GND for every internal node in a
# single pass (the SOS is calculated anyway when getting GND, so this avoids
# recomputing it later). Returns a `NodeAnalysis` to hand to the explore functions.
# The "compute once, explore a lot" entry point - cache the result (e.g. with JLD2)
# and reload it. See `node_metrics` for the effect-size scores (RMS/spatial/SES).
function node_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 100)
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    gnd = Dict{String, Float64}()
    sos = Dict{String, Vector{Float64}}()
    @progress for node in nodevec
        s, g = process_node(assemblage, tree, node; nsims)
        gnd[node] = g
        isnan(g) || (sos[node] = s)
    end
    NodeAnalysis(gnd, sos)
end

# Like `NodeAnalysis`, but also carries the effect-size scores: `rms` (total
# intensity), `spatial` (spatial-only), and `ses` (standardised effect size). Kept
# separate from `NodeAnalysis` so existing cached results still load unchanged.
struct NodeMetrics
    gnd::Dict{String, Float64}
    rms::Dict{String, Float64}
    spatial::Dict{String, Float64}
    ses::Dict{String, Float64}
    sos::Dict{String, Vector{Float64}}
end

# One-pass analysis returning GND together with the effect-size alternatives, all
# from the same null draws (one randomisation per node). The published swap null
# fixes range size, so the `spatial`/GND scores measure spatial turnover rather than
# richness asymmetry. Cache the result with JLD2 as for `node_analysis`.
function node_metrics(assemblage::Assemblage, tree::AbstractTree; nsims = 100)
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    gnd = Dict{String, Float64}(); rms = Dict{String, Float64}()
    spatial = Dict{String, Float64}(); ses = Dict{String, Float64}()
    sos = Dict{String, Vector{Float64}}()
    @progress for node in nodevec
        clade = get_clade(assemblage, tree, node)
        children = getchildren(tree, node)
        if length(children) != 2 || any(x -> isleaf(tree, x) || nspecies(get_clade(assemblage, tree, x)) < 4, children)
            gnd[node] = NaN
            continue
        end
        sims = simulate_descendants(clade, tree, children[1]; nsims)
        gnd[node]     = calculate_GND(sims)
        rms[node]     = calculate_GND_rms(sims)
        spatial[node] = calculate_GND_spatial(sims)
        ses[node]     = calculate_GND_ses(sims)
        sos[node]     = calculate_SOS(sims)
    end
    NodeMetrics(gnd, rms, spatial, ses, sos)
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
divergent_nodes(res::NodeMetrics; threshold = 0.8) = divergent_nodes(res.gnd; threshold)

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
sos_distances(assemblage, tree, nodes; nsims = 100) =
    _sos_distances(reduce(hcat, process_node(assemblage, tree, node; nsims)[1] for node in nodes))
