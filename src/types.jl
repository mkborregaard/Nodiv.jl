# Supertype of the analysis results, so the explore functions accept either.
abstract type AbstractNodeResult end

# Result of `node_analysis`: `nodes` lists every internal node in the order of a preorder
# traversal (the order they appear on the tree); `gnd` maps every node name to its GND
# (NaN where the node cannot be analysed); `sos` maps the analysable nodes to their
# per-cell SOS vector. Pass it straight to the explore functions - `divergent_nodes`,
# `plot_gnd`, `plot_node`, `sos_distances` - to reuse cached results.
struct NodeAnalysis <: AbstractNodeResult
    nodes::Vector{String}
    gnd::Dict{String,Float64}
    sos::Dict{String,Vector{Float64}}
end

# Like `NodeAnalysis`, but also carries the effect-size scores: `rms` (total
# intensity), `spatial` (spatial-only), `ses` (standardised effect size), `pval`
# (null-calibrated Monte-Carlo significance), and `varying`, the share of the focal
# clade's occupied cells where the null model varies (the cells the scores are based on).
struct NodeMetrics <: AbstractNodeResult
    nodes::Vector{String}
    gnd::Dict{String,Float64}
    rms::Dict{String,Float64}
    spatial::Dict{String,Float64}
    ses::Dict{String,Float64}
    pval::Dict{String,Float64}
    varying::Dict{String,Float64}
    sos::Dict{String,Vector{Float64}}
end
