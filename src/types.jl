# Result of `node_analysis`: `gnd` maps every node name to its GND (NaN where the
# node cannot be analysed); `sos` maps the analysable nodes to their per-cell SOS
# vector. Pass it straight to the explore functions - `divergent_nodes`,
# `plot_gnd`, `plot_node`, `sos_distances` - to reuse cached results.
struct NodeAnalysis
    gnd::Dict{String, Float64}
    sos::Dict{String, Vector{Float64}}
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
