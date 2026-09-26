# Two tight groups of patterns (a1-a3 and b1-b2) and two patterns on their own
function grouped_sos()
    rng = Xoshiro(5)
    a, b = randn(rng, 200), randn(rng, 200)
    noisy(x, s) = x .+ s .* randn(rng, 200)
    sosvectors = [
        noisy(a, 0.1),
        noisy(a, 0.1),
        -noisy(a, 0.1),  # A mirror image is the same pattern
        noisy(b, 0.1),
        noisy(b, 0.1),
        randn(rng, 200),
        randn(rng, 200),
    ]
    nodes = ["a1", "a2", "a3", "b1", "b2", "c", "d"]
    return sos_distances(sosvectors), nodes
end

@testset "sos_ordination" begin
    # Points in the plane: their distances are reproduced exactly by two axes
    rng = Xoshiro(2)
    points = randn(rng, 2, 6)
    D = [sqrt(sum(abs2, points[:, i] - points[:, j])) for i in 1:6, j in 1:6]
    nodes = ["n$i" for i in 1:6]
    o = sos_ordination(D, nodes)
    @test o isa SOSOrdination
    @test o.nodes == nodes
    @test size(o.coords) == (2, 6)
    @test o.distances == D
    fitted = [sqrt(sum(abs2, o.coords[:, i] - o.coords[:, j])) for i in 1:6, j in 1:6]
    @test fitted ≈ D
    @test issorted(o.eigenvalues; rev=true)

    # More axes than the data have: only the positive eigenvalues are kept (two, and any
    # more at round-off level), and the other axes are zero
    o5 = sos_ordination(D, nodes; maxoutdim=5)
    @test size(o5.coords) == (5, 6)
    k = length(o5.eigenvalues)
    @test k >= 2 && all(>(0), o5.eigenvalues)
    @test all(<(1e-10), o5.eigenvalues[3:end])
    @test all(iszero, o5.coords[(k + 1):5, :])
    @test o5.coords[1:2, :] ≈ o.coords

    @test_throws ArgumentError sos_ordination(D, nodes[1:5])
    @test_throws ArgumentError sos_ordination(D[1:2, 1:2], nodes[1:2])

    # From an analysis result or a Dict of SOS vectors
    assemblage, tree = toy_data()
    res = node_metrics(assemblage, tree; nsims=99)
    o = sos_ordination(res, TOY_ANALYSABLE; minoverlap=2)
    @test o.distances == sos_distances(res, TOY_ANALYSABLE; minoverlap=2)
    @test sos_ordination(res.sos, TOY_ANALYSABLE; minoverlap=2).coords == o.coords
    @test sprint(show, o) == "SOSOrdination(3 nodes, 2 axes)"
end

@testset "sos_clusters" begin
    D, nodes = grouped_sos()
    c = sos_clusters(D, nodes; simcut=0.7)
    @test c isa SOSClusters
    @test c.nodes == nodes
    @test c.distances == D
    @test c.simcut == 0.7
    g = c.groups
    @test g["a1"] == g["a2"] == g["a3"]
    @test g["b1"] == g["b2"]
    @test length(unique(values(g))) == 4
    # Labels 1:m for the clusters of more than one node, in the order of their numbers
    @test sort(collect(keys(c.labels))) == sort(unique([g["a1"], g["b1"]]))
    @test sort(collect(values(c.labels))) == [1, 2]
    @test issorted(sort(collect(c.labels)); by=last)

    sizes = sos_cluster_sizes(c)
    @test sizes.nclusters == 4
    @test sizes.nsingletons == 2
    @test sizes.members[g["a1"]] == ["a1", "a2", "a3"]
    @test sos_cluster_sizes(c.groups) == sizes

    shown = sprint(show, MIME("text/plain"), c)
    @test startswith(shown, "SOSClusters of 7 nodes, cut at |r| >= 0.7: 4 clusters, 2 ")
    @test occursin("cluster $(c.labels[g["a1"]]): a1, a2, a3", shown)
    @test sprint(show, c) == "SOSClusters(7 nodes, 4 clusters, 2 of more than one node)"

    # Nothing similar enough: every node on its own
    alone = sos_clusters(D, nodes; simcut=0.9999)
    @test isempty(alone.labels)
    @test sos_cluster_sizes(alone).nsingletons == 7
end

@testset "sos_similarity_communities" begin
    D, nodes = grouped_sos()
    c = sos_similarity_communities(D, nodes; simthresh=0.7)
    @test c isa SOSCommunities
    @test Set(Set.(c.communities)) == Set([Set(["a1", "a2", "a3"]), Set(["b1", "b2"])])
    # Two cliques, 3 and 1 links: Q = sum over them of links / m - (degree sum / 2m)^2
    @test c.modularity ≈ (3 / 4 - (6 / 8)^2) + (1 / 4 - (2 / 8)^2)
    shown = sprint(show, MIME("text/plain"), c)
    @test startswith(shown, "SOSCommunities of 7 nodes, linked at |r| >= 0.7")
    @test occursin("community 1: a1, a2, a3", shown)

    # No links: no communities
    none = sos_similarity_communities(D, nodes; simthresh=0.9999)
    @test isempty(none.communities)
    @test none.modularity == 0

    # A chain of two cliques joined by one link is split, where connected components
    # would make it one group
    n = 8
    linkedpairs = [(i, j) for i in 1:4 for j in (i + 1):4]
    append!(linkedpairs, [(i, j) for i in 5:8 for j in (i + 1):8])
    push!(linkedpairs, (4, 5))
    Dchain = ones(n, n)
    for i in 1:n
        Dchain[i, i] = 0
    end
    for (i, j) in linkedpairs
        Dchain[i, j] = Dchain[j, i] = 0.1
    end
    chain = sos_similarity_communities(Dchain, string.(1:n); simthresh=0.5)
    @test Set(Set.(chain.communities)) == Set([Set(string.(1:4)), Set(string.(5:8))])
end
