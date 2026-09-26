"""
    Nodiv

Node-based analysis of species distributions: the method of Borregaard et al. (2014),
Node-based analysis of species distributions, *Methods in Ecology and Evolution* 5:
1225-1235. For every internal node of a phylogeny, the richness of its two descendant
clades in each cell is compared with a null model that keeps range sizes and cell
richness. Start with [`node_metrics`](@ref).
"""
module Nodiv

using Clustering: Hclust, cutree, hclust
using LinearAlgebra: Hermitian, eigen!
using Phylo: Phylo
using Phylo: AbstractTree, getchildren, getdescendants, getleafnames, getnodename
using Phylo: isleaf, keeptips!, preorder, traversal
using Random: AbstractRNG, Xoshiro, default_rng, rand!
using RecipesBase: @recipe, @series, @userplot
using SparseArrays: nonzeros, nzrange, rowvals, sparse
using SpatialEcology: Assemblage, matrixrandomizer, nsites, nspecies, occurrences, richness
using SpatialEcology: speciesnames
using Statistics: cor, mean, std
using StatsBase: corspearman, tiedrank

export AbstractNodeResult, NodeAnalysis, NodeMetrics
export clade_richness, get_clade, nodespecies, prune_to_shared!
export simulate_descendants
export divergence_pval, divergence_ses, gnd, sos, sos_rms, sos_sd
export node_analysis, node_metrics, process_node
export default_score, divergent_nodes, most_divergent, node_scores
export sos_distances, sos_ordination, SOSOrdination
export sos_clusters, sos_cluster_sizes, SOSClusters
export sos_similarity_communities, SOSCommunities
export plot_gnd, plot_node

# Deprecated
export calculate_GND, calculate_GND_pval, calculate_GND_rms, calculate_GND_ses
export calculate_GND_spatial, calculate_SOS, gnd_rms, gnd_spatial, node_based_analysis

include("types.jl")
include("clades.jl")
include("null.jl")
include("metrics.jl")
include("analysis.jl")
include("distances.jl")
include("ordination.jl")
include("grouping.jl")
include("recipes.jl")
include("deprecated.jl")

end
