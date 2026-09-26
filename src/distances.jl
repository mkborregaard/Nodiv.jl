# The distances between the columns of `sosmat`; see `sos_distances`
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
    richness_of = clade_richness(assemblage, tree)
    return reduce(hcat, (richness_of(n) .> 0) for n in nodes)
end

"""
    sos_distances(res, nodes; minoverlap=3, method=:pearson, overlapweight=false,
                  assemblage=nothing, tree=nothing) -> Matrix{Float64}
    sos_distances(sosvectors; minoverlap=3, method=:pearson, overlapweight=false,
                  occupied=nothing) -> Matrix{Float64}

Pairwise distances between the SOS maps of `nodes`, for grouping nodes by how alike
their divergence geography is. `res` is the result of [`node_metrics`](@ref) or
[`node_analysis`](@ref), or a Dict of node name => SOS vector; or pass the SOS vectors
themselves. The result is symmetric with a zero diagonal, one row and column per node. It
is what [`sos_ordination`](@ref), [`sos_clusters`](@ref) and
[`sos_similarity_communities`](@ref) group the nodes by.

The distance is `1 - |r|`, with `r` the correlation of the two SOS maps over the cells
where both are defined: where both parent clades are present and the null model varies.
Using `|r|` makes a mirror image the same pattern, as which daughter is "first" at a node
is arbitrary. A constant map has no correlation, so distance 1.

# Keywords
- `minoverlap::Integer=3`: pairs sharing fewer cells with both SOS defined get distance 1
  rather than a correlation from a handful of cells; this also puts pairs with no shared
  cells at the maximum. Choose it per space: a floor of about 5-10 cells suits a
  geographic grid of thousands of cells, a lower one an environmental space of tens of
  bins. Occupied cells where the null model cannot vary have no SOS and do not count;
  they can be most of a small, sympatric clade's range.
- `method::Symbol=:pearson`: or `:spearman`, the robust choice if the SOS values are
  heavy-tailed, as they tend to be at strongly divergent nodes.
- `overlapweight::Bool=false`: the retained alternative, `D = 1 - O * |r|`, with `O` the
  Sørensen index of the two nodes' occupied cells, so a high correlation over a small
  shared range does not read as full similarity. Occupancy is derived from `assemblage`
  and `tree` (or given as an `occupied` matrix, one column per node) and enters only here,
  never the choice of cells for the correlation.
"""
function sos_distances(sosvectors::AbstractVector; kw...)
    return _sos_distances(reduce(hcat, sosvectors); kw...)
end

sos_distances(res::AbstractNodeResult, nodes; kw...) = sos_distances(res.sos, nodes; kw...)

function sos_distances(
    sos::AbstractDict,
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
                "overlapweight = true on a result needs `assemblage` and `tree` to " *
                "derive occupancy"
            throw(ArgumentError(msg))
        end
        occupied = _occupied_matrix(assemblage, tree, nodes)
    end
    return _sos_distances(
        reduce(hcat, sos[n] for n in nodes); overlapweight, occupied, kw...
    )
end
