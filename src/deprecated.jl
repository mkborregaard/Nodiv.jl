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
    sd::AbstractDict,
    ses::AbstractDict,
    pval::AbstractDict,
    sos::AbstractDict,
)
    Base.depwarn(
        "`NodeMetrics(gnd, rms, sd, ses, pval, sos)` is deprecated; use " *
        "`NodeMetrics(nodes, gnd, rms, sd, ses, pval, varying, sos)`. " *
        "The share of varying cells is unknown here and set to NaN.",
        :NodeMetrics,
    )
    varying = Dict{String,Float64}(n => NaN for n in keys(sos))
    return NodeMetrics(sort!(collect(keys(gnd))), gnd, rms, sd, ses, pval, varying, sos)
end

# Renamed metrics
@deprecate calculate_SOS(sims) sos(sims) false
@deprecate calculate_GND(sims) gnd(sims) false
@deprecate calculate_GND(sims, occupied) gnd(sims, occupied) false
@deprecate gnd_rms(sos_scores::AbstractVector) sos_rms(sos_scores) false
@deprecate gnd_spatial(sos_scores::AbstractVector) sos_sd(sos_scores) false

# The RMS-SOS and SD-SOS over all occupied cells, counting those where the null model
# cannot vary as SOS = 0. `sos_rms`/`sos_sd` leave those cells out instead.
function calculate_GND_rms(sims, occupied=_occupied(sims))
    Base.depwarn(
        "`calculate_GND_rms` is deprecated; use `sos_rms`, which leaves out the occupied " *
        "cells where the null model cannot vary instead of counting them as 0. The old " *
        "value is `sqrt(varying) * sos_rms(sims)`, with `varying` the share of occupied " *
        "cells where the null varies.",
        :calculate_GND_rms,
    )
    sos_scores = sos(sims)
    idx = findall(occupied)
    isempty(idx) && return NaN
    return sqrt(mean(j -> isfinite(sos_scores[j]) ? sos_scores[j]^2 : 0.0, idx))
end

function calculate_GND_spatial(sims, occupied=_occupied(sims))
    Base.depwarn(
        "`calculate_GND_spatial` is deprecated; use `sos_sd`, which leaves out the " *
        "occupied cells where the null model cannot vary instead of counting them as 0.",
        :calculate_GND_spatial,
    )
    sos_scores = sos(sims)
    idx = findall(occupied)
    isempty(idx) && return NaN
    return std([isfinite(sos_scores[j]) ? sos_scores[j] : 0.0 for j in idx])
end

function calculate_GND_ses(sims, occupied=_occupied(sims))
    Base.depwarn(
        "`calculate_GND_ses` is deprecated; use `divergence_ses`.", :calculate_GND_ses
    )
    return _ses(_msos_null(sims, occupied)...)
end

function calculate_GND_pval(sims, occupied=_occupied(sims))
    Base.depwarn(
        "`calculate_GND_pval` is deprecated; use `divergence_pval`.", :calculate_GND_pval
    )
    return _pval(_msos_null(sims, occupied)...)
end

# The serial analysis that recreates `Node_analysis` of the nodiv R package
function node_based_analysis(assemblage::Assemblage, tree::AbstractTree; nsims=200)
    Base.depwarn(
        "`node_based_analysis` is deprecated; use `node_analysis`, which returns a " *
        "`NodeAnalysis` and runs multithreaded.",
        :node_based_analysis,
    )
    nodevec = _internalnodes(tree)
    SOSs = Matrix{Float64}(undef, nsites(assemblage), length(nodevec))
    GNDs = Vector{Float64}(undef, length(nodevec))
    for (i, node) in enumerate(nodevec)
        SOSs[:, i], GNDs[i] = process_node(assemblage, tree, node; nsims)
    end
    return SOSs, GNDs
end

# SOS distances that rerun the null model for every node
function sos_distances(assemblage, tree, nodes; nsims=200, overlapweight::Bool=false, kw...)
    Base.depwarn(
        "`sos_distances(assemblage, tree, nodes)` is deprecated, as it reruns the null " *
        "model; use `sos_distances(res, nodes)` on the result of `node_metrics` or " *
        "`node_analysis`.",
        :sos_distances,
    )
    sosmat = reduce(hcat, process_node(assemblage, tree, node; nsims)[1] for node in nodes)
    occupied = overlapweight ? _occupied_matrix(assemblage, tree, nodes) : nothing
    return _sos_distances(sosmat; overlapweight, occupied, kw...)
end
