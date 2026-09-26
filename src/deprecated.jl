# The result types before they carried the node order (and, for `NodeMetrics`, the
# share of varying cells). Without a tree the nodes are taken in alphabetical order.
function NodeAnalysis(gnd::AbstractDict, sos::AbstractDict)
    Base.depwarn(
        "`NodeAnalysis(gnd, sos)` is deprecated; use `NodeAnalysis(nodes, gnd, sos)` " *
        "with `nodes` the internal nodes in tree order.",
        :NodeAnalysis,
    )
    return NodeAnalysis(sort!(collect(keys(gnd))), gnd, sos)
end

function NodeMetrics(
    gnd::AbstractDict,
    rms::AbstractDict,
    spatial::AbstractDict,
    ses::AbstractDict,
    pval::AbstractDict,
    sos::AbstractDict,
)
    Base.depwarn(
        "`NodeMetrics(gnd, rms, spatial, ses, pval, sos)` is deprecated; use " *
        "`NodeMetrics(nodes, gnd, rms, spatial, ses, pval, varying, sos)`. " *
        "The share of varying cells is unknown here and set to NaN.",
        :NodeMetrics,
    )
    varying = Dict{String,Float64}(n => NaN for n in keys(sos))
    return NodeMetrics(
        sort!(collect(keys(gnd))), gnd, rms, spatial, ses, pval, varying, sos
    )
end
