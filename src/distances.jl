# The distance between two SOS maps of similarity `s` (an |r| or |cos| in [0, 1]), and back.
# sqrt(1 - s^2) is the sine of the angle between the two (standardised) maps. It ignores the
# sign, like |r|, and, unlike 1 - |r|, is a Euclidean distance when every pair's similarity
# is taken over the same cells: it is 1/sqrt(2) times the distance between the maps' outer
# products x*x'/|x|^2, which are the same for a map and its mirror image.
_sos_distance(s::Real) = sqrt(max(1 - s^2, 0.0))

"""
    sos_similarity(d)

The similarity `|r|` (or `|cos|`, for `method=:cosine`) of two SOS maps at the distance `d`
of [`sos_distances`](@ref): `sqrt(1 - d^2)`, the inverse of `d = sqrt(1 - |r|^2)`. Broadcast
it over a distance matrix, `sos_similarity.(D)`, for the similarities.
"""
sos_similarity(d::Real) = sqrt(max(1 - d^2, 0.0))

# The distances between the columns of `sosmat`; see `sos_distances`
function _sos_distances(
    sosmat;
    minoverlap::Integer=3,
    method::Symbol=:pearson,
    overlapweight::Bool=false,
    occupied=nothing,
)
    if !(method in (:pearson, :spearman, :cosine))
        msg = "`method` must be :pearson, :spearman or :cosine; got $(repr(method))"
        throw(ArgumentError(msg))
    end
    if overlapweight && occupied === nothing
        msg =
            "overlapweight = true needs an `occupied` mask: pass `assemblage` and " *
            "`tree`, or an explicit `occupied` matrix"
        throw(ArgumentError(msg))
    end
    n = size(sosmat, 2)
    S = if method === :cosine
        _abs_cosines(sosmat)
    else
        _abs_correlations(sosmat, method, minoverlap)
    end
    D = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        sim = S[i, j]
        if overlapweight
            oi, oj = view(occupied, :, i), view(occupied, :, j)
            denom = count(oi) + count(oj)
            sim *= denom == 0 ? 0.0 : 2count(oi .& oj) / denom
        end
        D[i, j] = D[j, i] = _sos_distance(sim)
    end
    return D
end

# |r| between the columns of `sosmat`, each pair over the cells where both are defined; 0
# for a pair sharing fewer than `minoverlap` such cells, or with a constant column there
function _abs_correlations(sosmat, method, minoverlap)
    corfun = method === :pearson ? cor : corspearman
    n = size(sosmat, 2)
    S = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        a, b = view(sosmat, :, i), view(sosmat, :, j)
        ok = .!(isnan.(a) .| isnan.(b))          # shared finite == shared occupied cells
        count(ok) < minoverlap && continue
        r = corfun(a[ok], b[ok])
        S[i, j] = S[j, i] = isnan(r) ? 0.0 : abs(r)  # constant column -> r is NaN -> 0
    end
    return S
end

# |cos| between the columns of `sosmat`, with the cells where an SOS is not defined taken as
# 0, i.e. as showing no divergence. Cells that are 0 in both columns add nothing, so each
# pair is in effect compared over the union of its cells. A column with no SOS anywhere has
# similarity 0 to all others.
function _abs_cosines(sosmat)
    Z = map(x -> isnan(x) ? 0.0 : Float64(x), sosmat)
    norms = vec(sqrt.(sum(abs2, Z; dims=1)))
    G = Z' * Z
    n = size(Z, 2)
    S = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        denom = norms[i] * norms[j]
        S[i, j] = S[j, i] = denom == 0 ? 0.0 : min(abs(G[i, j]) / denom, 1.0)
    end
    return S
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

The distance is `sqrt(1 - r^2)`, with `r` the correlation of the two SOS maps over the
cells where both are defined: where both parent clades are present and the null model
varies. It depends on `r` only through `|r|`, which makes a mirror image the same pattern,
as which daughter is "first" at a node is arbitrary. It runs from 0 (`|r| = 1`) to 1
(`r = 0`), and [`sos_similarity`](@ref) turns it back into `|r|`. A constant map has no
correlation, so distance 1.

`sqrt(1 - r^2)` rather than `1 - |r|` because it is the sine of the angle between the two
maps, a proper distance: `1 - |r|` breaks the triangle inequality (two maps can each
correlate at 0.7 with a third and not at all with each other), which an ordination of it
shows by spreading nodes that are similar apart. Both order the pairs the same way, so
complete-linkage clusters cut at a given `|r|` are the same with either.

# Keywords
- `minoverlap::Integer=3`: pairs sharing fewer cells with both SOS defined get distance 1
  rather than a correlation from a handful of cells; this also puts pairs with no shared
  cells at the maximum. Choose it per space: a floor of about 5-10 cells suits a
  geographic grid of thousands of cells, a lower one an environmental space of tens of
  bins. Occupied cells where the null model cannot vary have no SOS and do not count;
  they can be most of a small, sympatric clade's range. Not used by `method=:cosine`.
- `method::Symbol=:pearson`: or `:spearman`, the robust choice if the SOS values are
  heavy-tailed, as they tend to be at strongly divergent nodes; or `:cosine`, the
  alternative below.
- `overlapweight::Bool=false`: the retained alternative, `D = sqrt(1 - (O * |r|)^2)`, with
  `O` the Sørensen index of the two nodes' occupied cells, so a high correlation over a
  small shared range does not read as full similarity. Occupancy is derived from
  `assemblage` and `tree` (or given as an `occupied` matrix, one column per node) and
  enters only here, never the choice of cells for the correlation.

# `method = :cosine`
The similarity is then `|cos|`, the uncentred correlation of the two maps over all cells,
with a cell where a map has no SOS counted as 0: no divergence, which is what an SOS of 0
means. Cells where neither map has an SOS add nothing, so this compares the two maps over
the union of their cells rather than the cells they share. Two clades that agree closely
where they overlap but each have much of their divergence outside the overlap are then
not similar; with a correlation they would be. It also needs no `minoverlap` floor, as a
small overlap already gives a small `|cos|`. As every map is then a vector over the same
cells, the distances are exactly Euclidean. Note that "no SOS" includes the occupied cells
where the null model cannot vary, not only the cells outside a clade's range.
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
