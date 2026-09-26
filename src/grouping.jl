# Grouping nodes by the similarity of their SOS maps, from the distances of `sos_distances`:
# clusters from hierarchical clustering, and communities in a thresholded similarity graph.

"""
    SOSClusters

The result of [`sos_clusters`](@ref). Shown, it lists the clusters of more than one node.

# Fields
- `nodes::Vector{String}`: the nodes clustered, in the order of the rows of `distances`
- `distances::Matrix{Float64}`: their SOS distance matrix
- `hclust`: the hierarchical clustering (a `Clustering.Hclust`); `nodes[hclust.order]` is
  the order of the dendrogram's leaves
- `simcut::Float64`: the similarity `|r|` the dendrogram is cut at
- `groups::Dict{String,Int}`: the cluster of every node, numbered as by `Clustering.cutree`
- `labels::Dict{Int,Int}`: the clusters of more than one node, from their `groups` number to
  a label `1:m` in the order of those numbers. The co-patterned groups, as opposed to the
  idiosyncratic nodes left on their own, are referred to by these labels.
"""
struct SOSClusters
    nodes::Vector{String}
    distances::Matrix{Float64}
    hclust::Hclust
    simcut::Float64
    groups::Dict{String,Int}
    labels::Dict{Int,Int}
end

"""
    sos_clusters(D, nodes; simcut, linkage=:complete) -> SOSClusters

Hierarchical clustering of `nodes` by the SOS distances `D` (from
[`sos_distances`](@ref)), cut where the similarity `|r|` falls below `simcut`, i.e. at the
distance `1 - simcut`.

Complete linkage, the default, is the conservative choice: it puts nodes in a cluster only
if all pairs in it are similar, and does not chain marginal pairs. Genuine groups, if any,
show up as clusters of more than one node against a background of nodes on their own,
which [`sos_cluster_sizes`](@ref) counts. Few such clusters is a result, not a failure: the
divergent nodes are then largely idiosyncratic in their SOS pattern.
"""
function sos_clusters(D::AbstractMatrix{<:Real}, nodes; simcut, linkage=:complete)
    nodes = String[n for n in nodes]
    D = Matrix{Float64}(D)
    hc = hclust(D; linkage)
    cut = cutree(hc; h=1 - simcut)
    groups = Dict(node => cut[i] for (i, node) in enumerate(nodes))
    ids = sort(collect(keys(sos_cluster_sizes(groups).members)))
    labels = Dict(c => i for (i, c) in enumerate(ids))
    return SOSClusters(nodes, D, hc, simcut, groups, labels)
end

"""
    sos_cluster_sizes(clusters) -> (; nclusters, nsingletons, members)

How many nodes fall in each cluster of an [`SOSClusters`](@ref) (or of a Dict of node =>
cluster). A near-flat table of singletons is the "largely idiosyncratic" result; a few
clusters of several nodes are the co-patterned exceptions.

Returns the number of clusters, the number of them with a single node, and `members`, a
Dict from each cluster of more than one node to its sorted node names.
"""
sos_cluster_sizes(clusters::SOSClusters) = sos_cluster_sizes(clusters.groups)

function sos_cluster_sizes(groups::AbstractDict)
    counts = Dict{Int,Int}()
    for c in values(groups)
        counts[c] = get(counts, c, 0) + 1
    end
    nontrivial = [c for (c, k) in counts if k > 1]
    members = Dict(c => sort([n for (n, g) in groups if g == c]) for c in nontrivial)
    return (; nclusters=length(counts), nsingletons=count(==(1), values(counts)), members)
end

function Base.show(io::IO, c::SOSClusters)
    sizes = sos_cluster_sizes(c)
    print(io, "SOSClusters(", length(c.nodes), " nodes, ", sizes.nclusters, " clusters, ")
    print(io, length(c.labels), " of more than one node)")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", c::SOSClusters)
    sizes = sos_cluster_sizes(c)
    print(io, "SOSClusters of ", length(c.nodes), " nodes, cut at |r| >= ", c.simcut, ": ")
    print(io, sizes.nclusters, " clusters, ", sizes.nsingletons, " singletons")
    for (id, label) in sort(collect(c.labels); by=last)
        print(io, "\n  cluster ", label, ": ", join(sizes.members[id], ", "))
    end
    return nothing
end

"""
    SOSCommunities

The result of [`sos_similarity_communities`](@ref). Shown, it lists the communities.

# Fields
- `nodes::Vector{String}`: the nodes of the graph
- `simthresh::Float64`: the similarity `|r|` at or above which two nodes are linked
- `communities::Vector{Vector{String}}`: the communities of more than one node
- `modularity::Float64`: the modularity Q of the division of the graph into communities
"""
struct SOSCommunities
    nodes::Vector{String}
    simthresh::Float64
    communities::Vector{Vector{String}}
    modularity::Float64
end

"""
    sos_similarity_communities(D, nodes; simthresh) -> SOSCommunities

Communities of `nodes` in their thresholded similarity graph, a check on
[`sos_clusters`](@ref). Two nodes are linked when the similarity `|r|` of their SOS maps
(`1 - D` for the distances `D` from [`sos_distances`](@ref)) is at least `simthresh`. The
overlap floor is already in `D`: pairs sharing too few cells, or none, are at distance 1, so
never linked.

The communities are those that maximise the graph's modularity Q, found by greedy
agglomeration (Clauset, Newman & Moore 2004): groups more densely linked inside than a
graph with the same degrees would be by chance. Unlike the connected components, which
chain any path of links into one group, this does not merge groups joined by a few links.
"""
function sos_similarity_communities(D::AbstractMatrix{<:Real}, nodes; simthresh)
    nodes = String[n for n in nodes]
    n = length(nodes)
    linked = falses(n, n)
    for i in 1:n, j in (i + 1):n
        linked[i, j] = linked[j, i] = (1 - D[i, j]) >= simthresh
    end
    m = count(linked) ÷ 2
    m == 0 && return SOSCommunities(nodes, simthresh, Vector{String}[], 0.0)
    degree = vec(sum(linked; dims=2))
    comm = collect(1:n)
    while true
        best, bestpair = 0.0, nothing
        labels = unique(comm)
        members = Dict(c => findall(==(c), comm) for c in labels)
        for a in labels, b in labels
            a < b || continue
            links = count(view(linked, members[a], members[b]))
            links == 0 && continue
            dQ = links / m - sum(degree[members[a]]) * sum(degree[members[b]]) / (2m^2)
            if dQ > best
                best, bestpair = dQ, (a, b)
            end
        end
        bestpair === nothing && break
        comm[comm .== bestpair[2]] .= bestpair[1]
    end
    groups = filter(c -> length(c) > 1, [findall(==(c), comm) for c in unique(comm)])
    communities = [[nodes[i] for i in c] for c in groups]
    return SOSCommunities(nodes, simthresh, communities, _modularity(linked, degree, comm))
end

# Newman's modularity of the division `comm` of the graph with adjacency matrix `linked`
function _modularity(linked, degree, comm)
    twom = sum(degree)
    return sum(unique(comm)) do c
        inside = comm .== c
        return count(view(linked, inside, inside)) / twom - (sum(degree[inside]) / twom)^2
    end
end

function Base.show(io::IO, c::SOSCommunities)
    print(io, "SOSCommunities(", length(c.nodes), " nodes, ", length(c.communities))
    print(io, " communities, Q = ", round(c.modularity; digits=3), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", c::SOSCommunities)
    print(io, "SOSCommunities of ", length(c.nodes), " nodes, linked at |r| >= ")
    print(io, c.simthresh, ": modularity Q = ", round(c.modularity; digits=3))
    for (i, members) in enumerate(c.communities)
        print(io, "\n  community ", i, ": ", join(members, ", "))
    end
    return nothing
end
