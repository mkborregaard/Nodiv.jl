# Pairwise distance matrix between per-cell SOS patterns, for grouping nodes by SOS-map
# similarity (see docs/sos_pattern_grouping_design.md). D(k,l) = 1 - |r|, with r the
# correlation of the two SOS vectors over the cells where BOTH are finite. A finite SOS means
# the parent clade is present and the null model varies there, the same cells the RMS-SOS
# is based on. Occupied cells where the null cannot vary (e.g. every species of the parent
# clade present) have no SOS and are left out; they can be most of a small, sympatric
# clade's range, so `minoverlap` counts finite cells, not occupied ones. |r| (not r) folds
# the arbitrary per-node daughter
# labelling: a mirror-image SOS map is the same divergence geography with the labels swapped.
#
# `minoverlap` guards the correlation's own sample size: pairs sharing fewer than `minoverlap`
# finite cells get D = 1 rather than a correlation fit on a handful of cells. This is what pins
# disjoint pairs (no shared occupied cells) at the maximum distance. Set it per space - the
# geographic scan has ~18k cells (a floor of 5-10 is reasonable), the environmental scan only
# tens of bins, so the floor must be chosen separately for each.
#
# `method` selects Pearson (default) or Spearman; Spearman is the robust choice if the marginal
# SOS distributions are heavy-tailed, as they tend to be at strongly divergent nodes.
#
# `overlapweight` is the RETAINED ALTERNATIVE (off by default): it multiplies |r| by a Sorensen
# index O of the two nodes' occupied cells (D = 1 - O*|r|), so high correlation over a small
# shared range does not read as full similarity. Occupancy enters ONLY here, never the
# correlation mask, and must be a genuine `richness(clade) .> 0` mask (supply `assemblage`/`tree`).
function _sos_distances(
    sosmat;
    minoverlap::Integer=3,
    method::Symbol=:pearson,
    overlapweight::Bool=false,
    occupied=nothing,
)
    corfun = if method === :pearson
        cor
    elseif method === :spearman
        corspearman
    else
        throw(ArgumentError("`method` must be :pearson or :spearman; got $(repr(method))"))
    end
    if overlapweight && occupied === nothing
        msg =
            "overlapweight = true needs an `occupied` mask: pass `assemblage` and " *
            "`tree`, or an explicit `occupied` matrix"
        throw(ArgumentError(msg))
    end
    n = size(sosmat, 2)
    D = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        a, b = view(sosmat, :, i), view(sosmat, :, j)
        ok = .!(isnan.(a) .| isnan.(b))          # shared finite == shared occupied cells
        if count(ok) < minoverlap
            D[i, j] = D[j, i] = 1.0
            continue
        end
        r = corfun(a[ok], b[ok])
        sim = isnan(r) ? 0.0 : abs(r)            # constant column -> r is NaN -> D = 1
        if overlapweight
            oi, oj = view(occupied, :, i), view(occupied, :, j)
            denom = count(oi) + count(oj)
            sim *= denom == 0 ? 0.0 : 2count(oi .& oj) / denom
        end
        D[i, j] = D[j, i] = 1 - sim
    end
    return D
end

# Occupied-cell mask for `nodes`, columns aligned with their SOS vectors. Deterministic and
# cheap (no null model), so it is safe to derive on the fly for the overlap weighting.
function _occupied_matrix(assemblage, tree, nodes)
    return reduce(hcat, (richness(get_clade(assemblage, tree, n)) .> 0) for n in nodes)
end

function sos_distances(sosvectors::AbstractVector; kw...)
    return _sos_distances(reduce(hcat, sosvectors); kw...)
end

function sos_distances(
    res::AbstractNodeResult,
    nodes;
    overlapweight::Bool=false,
    assemblage=nothing,
    tree=nothing,
    kw...,
)
    occupied = nothing
    if overlapweight
        if assemblage === nothing || tree === nothing
            msg =
                "overlapweight = true on a `res` needs `assemblage` and `tree` to " *
                "derive occupancy"
            throw(ArgumentError(msg))
        end
        occupied = _occupied_matrix(assemblage, tree, nodes)
    end
    return _sos_distances(
        reduce(hcat, res.sos[n] for n in nodes); overlapweight, occupied, kw...
    )
end

function sos_distances(assemblage, tree, nodes; nsims=200, overlapweight::Bool=false, kw...)
    sosmat = reduce(hcat, process_node(assemblage, tree, node; nsims)[1] for node in nodes)
    occupied = overlapweight ? _occupied_matrix(assemblage, tree, nodes) : nothing
    return _sos_distances(sosmat; overlapweight, occupied, kw...)
end
