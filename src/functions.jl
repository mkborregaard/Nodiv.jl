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
using StatsBase: tiedrank, sample
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

# Plot a node in a 2x2 grid: the parent clade (top-left), the descendant SOS
# mapped over the assemblage's coordinates (top-right, RdYlBu), and the two child
# clades on the bottom row. Defined as a plot recipe so the package depends only
# on RecipesBase; the plotting backend (Plots) is supplied by the caller. Use as
# `plot_node(assemblage, tree, node)`.
@userplot Plot_Node

@recipe function f(pn::Plot_Node)
    assemblage, tree, node = pn.args
    ch1, ch2 = getchildren(tree, node)[1:2]
    assm = get_clade(assemblage, tree, node)
    assmch1 = get_clade(assm, tree, ch1)
    assmch2 = get_clade(assm, tree, ch2)

    # SOS of the first descendant over the parent clade's cells. NB this runs a
    # randomization on every plot, so the SOS panel varies between calls.
    sos = calculate_SOS(simulate_descendants(assm, tree, ch1; method = :tipshuffle, nsims = 1000))

    layout := (2, 2)
    size --> (900, 800)

    @series begin              # top-left: parent clade
        subplot := 1
        title := "parent"
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
        title := "child 1"
        assmch1
    end
    @series begin              # bottom-right: child 2
        subplot := 4
        title := "child 2"
        assmch2
    end
end

# Plot per-node values (e.g. GND) on the tree. `gndvals` is a Dict of node name
# => value; a coloured marker is drawn at every node present in it (NaN values
# skipped) and nothing elsewhere. Pass the full `node_gnd` Dict to show all
# analysable nodes, or a filtered Dict (e.g. only divergent nodes) for a subset.
# Use as `plot_gnd(tree, gnd)`.
@userplot Plot_Gnd

@recipe function f(pg::Plot_Gnd)
    tree, gndvals = pg.args
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

# Build a sampling distribution of a descendant clade's per-site richness.
#
# :swap       - curveball randomization of the clade's presence/absence matrix
#               (keeps species ranges and site richness constant), then recompute
#               the descendant's richness. The published null model.
# :tipshuffle - keep the presence/absence matrix intact and randomize only which
#               species belong to the focal descendant, holding the descendant's
#               species count fixed. Much faster, as it avoids matrix swapping.
function simulate_descendants(clade, tree, descendant; method = :swap, nsims = 99)
    ret = zeros(nsims + 1, nsites(clade))  # a matrix to hold the richness values from the simulations
    # the empirical richness in the first row
    ret[1, :] = richness(get_clade(clade, tree, descendant))
    if method == :swap
        rmg = matrixrandomizer(clade)
        for i in 2:nsims + 1
            # and simulated richness in the rest of the nsims rows
            ret[i, :] .= richness(get_clade(rand!(rmg), tree, descendant))
        end
    elseif method == :tipshuffle
        # Materialize the clade once so the per-sim views are single-level (a
        # nested view of a view would not hit SpatialEcology's colsum method).
        cl = Assemblage(clade)
        nclade = nspecies(cl)
        ndesc = nspecies(get_clade(cl, tree, descendant))
        for i in 2:nsims + 1
            # draw a random set of `ndesc` species and take their per-site richness
            ret[i, :] .= richness(view(cl, species = sample(1:nclade, ndesc; replace = false)))
        end
    else
        error("unrecognized method")
    end
    ret
end

# Standardized effect size (SOS metric) per grid cell from a simulation matrix.
function calculate_SOS(sims)
    sd = std.(eachcol(sims))
    me = mean.(eachcol(sims))
    (sims[1, :] .- me) ./ sd
end

# Geographic node divergence (GND metric) from a simulation matrix.
function calculate_GND(sims)
  # two internal convenience functions
  logit(p) = log(p/(1-p))
  invlogit(p) = exp(p)/(1+exp(p))

  n = size(sims, 1)
  r = [tiedrank(x)[1]/(n + 1) for x in eachcol(sims)]
  p = 1 .- 2 .* abs.(r .- 0.5) .- 1/n
  α = mean(logit.(p))
  1-invlogit(α)
end

# Calculate SOS and GND for a single node (NaN when the node can't be analysed).
function process_node(assemblage, tree, nodename; nsims = 100, method = :swap)
    clade = get_clade(assemblage, tree, nodename)
    children = getchildren(tree, nodename)

    if length(children) != 2 || any(x -> isleaf(tree, x) || nspecies(get_clade(assemblage, tree, x)) < 4, children)
        return (fill(NaN, nsites(clade)), NaN)
    end

    sims = simulate_descendants(clade, tree, children[1]; nsims, method)
    calculate_SOS(sims), calculate_GND(sims)
end

# Run the node-based analysis over every internal node of the tree.
# Recreates the main `Node_analysis` function of the nodiv R package
# (https://github.com/mkborregaard/nodiv).
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 100, method = :swap)
   nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)] #shuffle!(collect(nodenamefilter(!isleaf, tree)))
   SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
   GNDs = Vector{Float64}(undef, length(nodevec))
   @progress for (i, node) in enumerate(nodevec)
       SOSs[:,i], GNDs[i] = process_node(assemblage, tree, node; nsims, method)
   end
   SOSs, GNDs
end

# Calculate the GND value for every internal node, returned as a Dict keyed by
# node name (GND = NaN where a node cannot be analysed). Lighter than
# `node_based_analysis`, which also builds the per-cell SOS maps - use this for a
# fast divergence scan. Defaults to the :tipshuffle null. The Dict can be passed
# straight to Phylo's tree plot recipe as `marker_z`, which looks values up by
# node name.
function node_gnd(assemblage::Assemblage, tree::AbstractTree; nsims = 100, method = :tipshuffle)
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    gnd = Dict{eltype(nodevec), Float64}()
    @progress for node in nodevec
        gnd[node] = process_node(assemblage, tree, node; nsims, method)[2]
    end
    gnd
end

# Compute both the per-cell SOS pattern and the GND for every internal node in a
# single pass (the SOS is calculated anyway when getting GND, so this avoids
# recomputing it later). Returns a NamedTuple `(; gnd, sos)`: `gnd` maps every
# node name to its GND (NaN where it cannot be analysed); `sos` maps the
# analysable nodes to their per-cell SOS vector. Defaults to the :tipshuffle null.
function node_analysis(assemblage::Assemblage, tree::AbstractTree; nsims = 100, method = :tipshuffle)
    nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    gnd = Dict{eltype(nodevec), Float64}()
    sos = Dict{eltype(nodevec), Vector{Float64}}()
    @progress for node in nodevec
        s, g = process_node(assemblage, tree, node; nsims, method)
        gnd[node] = g
        isnan(g) || (sos[node] = s)
    end
    (; gnd, sos)
end

# Prune `tree` in place to the tips it shares with all the given assemblage(s),
# so clade subsetting never references a species absent from the data. Returns
# the tree.
function prune_to_shared!(tree, assemblages...)
    shared = intersect(getleafnames(tree), speciesnames.(assemblages)...)
    keeptips!(tree, shared)
    tree
end

# Node names whose GND exceeds `threshold` (NaN GNDs excluded). Companion to
# `node_gnd`; pass its Dict.
divergent_nodes(gnd::AbstractDict; threshold = 0.8) =
    [node for (node, g) in gnd if !isnan(g) && g > threshold]

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
sos_distances(sosvectors) = _sos_distances(reduce(hcat, sosvectors))
sos_distances(assemblage, tree, nodes; nsims = 100, method = :tipshuffle) =
    _sos_distances(reduce(hcat, process_node(assemblage, tree, node; nsims, method)[1] for node in nodes))
