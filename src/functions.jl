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

# All tips/species descending from a node.
nodespecies(tree, node) = filter(x -> isleaf(tree, x), getdescendants(tree, node))

# Subset an Assemblage to the clade descending from a node.
get_clade(assemblage, tree, node) = view(assemblage, species = nodespecies(tree, node))

# Plot a parent clade alongside its two descendant clades.
# Defined as a plot recipe so the package depends only on RecipesBase; the
# actual plotting backend (Plots) is supplied by the caller. Use as
# `plot_node(assemblage, tree, node)`.
@userplot Plot_Node

@recipe function f(pn::Plot_Node)
    assemblage, tree, node = pn.args
    ch1, ch2 = getchildren(tree, node)[1:2]
    assm = get_clade(assemblage, tree, node)
    assmch1 = get_clade(assm, tree, ch1)
    assmch2 = get_clade(assm, tree, ch2)

    layout := (1, 3)
    size --> (1000, 350)

    @series begin
        subplot := 1
        title := "parent"
        assm
    end
    @series begin
        subplot := 2
        title := "child 1"
        assmch1
    end
    @series begin
        subplot := 3
        title := "child 2"
        assmch2
    end
end

# Build a sampling distribution of a descendant clade's richness via curveball randomization.
function simulate_descendants(clade, tree, descendant; nsims = 99)
    rmg = matrixrandomizer(clade)
    ret = zeros(nsims + 1, nsites(clade))  # a matrix to hold the richness values from the simulations
    # the empirical richness in the first row
    ret[1, :] = richness(get_clade(clade, tree, descendant))
    for i in 2:nsims + 1
        # and simulated richness in the rest of the nsims rows
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
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree)
   nodevec = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)] #shuffle!(collect(nodenamefilter(!isleaf, tree)))
   SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
   GNDs = Vector{Float64}(undef, length(nodevec))
   @progress for (i, node) in enumerate(nodevec)
       SOSs[:,i], GNDs[i] = process_node(assemblage, tree, node)
   end
   SOSs, GNDs
end
