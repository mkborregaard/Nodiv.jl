"""
    SOSOrdination

The result of [`sos_ordination`](@ref).

# Fields
- `nodes::Vector{String}`: the node names, in the order of the columns of `coords`
- `coords::Matrix{Float64}`: the MDS coordinates, one row per axis and one column per node
- `eigenvalues::Vector{Float64}`: the eigenvalues of the axes, largest first
- `distances::Matrix{Float64}`: the SOS distance matrix the ordination was fitted to
"""
struct SOSOrdination
    nodes::Vector{String}
    coords::Matrix{Float64}
    eigenvalues::Vector{Float64}
    distances::Matrix{Float64}
end

function Base.show(io::IO, o::SOSOrdination)
    print(io, "SOSOrdination(", length(o.nodes), " nodes, ", size(o.coords, 1), " axes)")
    return nothing
end

"""
    sos_ordination(res, nodes; maxoutdim=2, kw...) -> SOSOrdination
    sos_ordination(D::AbstractMatrix, nodes; maxoutdim=2) -> SOSOrdination

Classical (metric) multidimensional scaling of `nodes` by the similarity of their SOS maps.
The distances are [`sos_distances`](@ref)`(res, nodes; kw...)`, with `res` a
[`NodeMetrics`](@ref)/[`NodeAnalysis`](@ref) or a Dict of node name => SOS vector; nothing
is recomputed. Or pass the distance matrix `D` of `nodes` itself, e.g. one computed once
with `sos_distances` and shared with other analyses of the same nodes.

`maxoutdim` is the number of axes. The `eigenvalues` of an ordination with more axes (say
10) tell whether the first two show the structure. Where the nodes are mostly unrelated in
SOS pattern, the distances are all near 1 and the points form a ring. That is the finding,
not a failure of the method.

The ordination is that of `MultivariateStats.fit(MDS, D; distances=true, maxoutdim)`:
axes with a non-positive eigenvalue are left out of `eigenvalues` and get zero
coordinates.
"""
function sos_ordination(res, nodes; maxoutdim=2, kw...)
    nodes = _ordination_nodes(nodes)
    return sos_ordination(sos_distances(res, nodes; kw...), nodes; maxoutdim)
end

function sos_ordination(D::AbstractMatrix{<:Real}, nodes; maxoutdim=2)
    nodes = _ordination_nodes(nodes)
    n = length(nodes)
    if size(D) != (n, n)
        msg =
            "The distance matrix must be $n x $n, one row and column per node; got " *
            join(size(D), " x ")
        throw(ArgumentError(msg))
    end
    D = Matrix{Float64}(D)
    eigenvalues, coords = _classical_mds(D, maxoutdim)
    return SOSOrdination(nodes, coords, eigenvalues, D)
end

function _ordination_nodes(nodes)
    nodes = String[n for n in nodes]
    if length(nodes) < 3
        throw(ArgumentError("An ordination needs at least 3 nodes; got $(length(nodes))"))
    end
    return nodes
end

# Classical MDS of the distance matrix `D`, computed as MultivariateStats does, so the
# coordinates are the same to the last bit: the double-centred Gram matrix, its largest
# eigenvalues (all negated if the largest in magnitude is negative), and the positive ones
# kept. Returns the eigenvalues kept and the coordinates, `maxoutdim` x n.
function _classical_mds(D, maxoutdim)
    n = size(D, 1)
    u = [sum(abs2, view(D, :, j)) / n for j in 1:n]
    s = foldl(+, u; init=0.0) / n
    G = Matrix{Float64}(undef, n, n)
    for j in 1:n, i in j:n
        G[i, j] = G[j, i] = (u[i] + u[j] - abs2(D[i, j]) - s) / 2
    end
    E = eigen!(Hermitian(G))
    m = min(maxoutdim, n)
    mineig, maxeig = extrema(E.values)
    if mineig < 0 && abs(mineig) > abs(maxeig)
        order = sortperm(E.values)
        λ = -E.values[order[1:m]]
    else
        order = sortperm(E.values; rev=true)
        λ = E.values[order[1:m]]
    end
    kept = findfirst(<=(0), λ)
    if kept !== nothing
        m = kept - 1
        λ = λ[1:m]
    end
    coords = zeros(maxoutdim, n)
    coords[1:m, :] = sqrt.(λ) .* permutedims(E.vectors[:, order[1:m]])
    return λ, coords
end
