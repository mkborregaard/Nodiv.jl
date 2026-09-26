"""
    AbstractNodeResult

Supertype of the results of a node-based analysis, [`NodeAnalysis`](@ref) and
[`NodeMetrics`](@ref). Both have the fields `nodes` (every internal node, in the order
they appear on the tree), `gnd` and `sos`.
"""
abstract type AbstractNodeResult end

"""
    NodeAnalysis

The result of [`node_analysis`](@ref).

# Fields
- `nodes::Vector{String}`: every internal node, in the order they appear on the tree (a
  preorder traversal)
- `gnd::Dict{String,Float64}`: the GND of every internal node, `NaN` where the node cannot
  be analysed
- `sos::Dict{String,Vector{Float64}}`: the per-cell SOS of every analysed node
"""
struct NodeAnalysis <: AbstractNodeResult
    nodes::Vector{String}
    gnd::Dict{String,Float64}
    sos::Dict{String,Vector{Float64}}
end

"""
    NodeMetrics

The result of [`node_metrics`](@ref): the fields of a [`NodeAnalysis`](@ref) and the
effect-size scores of every analysed node.

# Fields
- `nodes::Vector{String}`: every internal node, in the order they appear on the tree
- `gnd::Dict{String,Float64}`: the GND of every internal node, `NaN` where the node cannot
  be analysed
- `rms::Dict{String,Float64}`: RMS-SOS, see [`sos_rms`](@ref)
- `spatial::Dict{String,Float64}`: the SD of the SOS, see [`sos_sd`](@ref)
- `ses::Dict{String,Float64}`: see [`divergence_ses`](@ref)
- `pval::Dict{String,Float64}`: see [`divergence_pval`](@ref)
- `varying::Dict{String,Float64}`: the share of the node's occupied cells where the null
  model varies, the cells the SOS-based scores rest on
- `sos::Dict{String,Vector{Float64}}`: the per-cell SOS
"""
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

function Base.show(io::IO, res::AbstractNodeResult)
    ncells = isempty(res.sos) ? 0 : length(first(values(res.sos)))
    print(io, nameof(typeof(res)), "(", length(res.nodes), " internal nodes, ")
    print(io, length(res.sos), " analysed, ", ncells, " cells)")
    return nothing
end
