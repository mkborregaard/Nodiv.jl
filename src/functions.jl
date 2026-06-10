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

    # SOS of the first descendant over the parent clade's cells (fast :tipshuffle null)
    sos = calculate_SOS(simulate_descendants(assm, tree, ch1; method = :swap))

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
