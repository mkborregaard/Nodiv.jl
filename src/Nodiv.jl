# Node-based analysis of species distributions.
#
# Extracted from the SpatialEcology.jl docs:
# https://docs.ecojulia.org/SpatialEcology.jl/stable/examples/nodebased/
#
# Reimplements the method of Borregaard et al. (2014), Node-based analysis of
# species distributions, Methods in Ecology and Evolution 5: 1225-1235.
module Nodiv

using Phylo: Phylo
using Phylo: AbstractTree, getchildren, getdescendants, getleafnames, getnodename, isleaf,
    keeptips!, preorder, traversal
using ProgressLogging: @progress
using Random: rand!
using RecipesBase: @recipe, @series, @userplot
using SpatialEcology: Assemblage, matrixrandomizer, nsites, nspecies, richness, speciesnames
using Statistics: cor, mean, std
using StatsBase: corspearman, tiedrank

export NodeAnalysis, NodeMetrics
export get_clade, nodespecies, prune_to_shared!
export simulate_descendants
export calculate_GND, calculate_SOS
export calculate_GND_pval, calculate_GND_rms, calculate_GND_ses, calculate_GND_spatial
export gnd_rms, gnd_spatial
export node_analysis, node_based_analysis, node_metrics, process_node
export divergent_nodes, sos_distances
export plot_gnd, plot_node

include("types.jl")
include("clades.jl")
include("null.jl")
include("metrics.jl")
include("analysis.jl")
include("distances.jl")
include("recipes.jl")

end
