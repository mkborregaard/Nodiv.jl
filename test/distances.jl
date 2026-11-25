@testset "sos_distances" begin
    rng = Xoshiro(3)
    a = randn(rng, 50)
    b = randn(rng, 50)

    D = sos_distances([a, b, -a])
    @test D == D'
    @test all(iszero, D[i, i] for i in 1:3)
    @test all(0 .<= D .<= 1)
    @test D[1, 2] ≈ sqrt(1 - cor(a, b)^2)
    @test sos_similarity.(D) ≈
        [1 abs(cor(a, b)) 1; abs(cor(a, b)) 1 abs(cor(a, b)); 1 abs(cor(a, b)) 1]
    @test D[1, 3] ≈ 0 atol = 1e-7  # a mirror image is the same pattern

    # NaN cells are left out pair by pair
    gappy = copy(a)
    gappy[1:10] .= NaN
    @test sos_distances([gappy, b])[1, 2] ≈ sqrt(1 - cor(a[11:end], b[11:end])^2)

    # fewer than `minoverlap` shared cells -> distance 1
    few = fill(NaN, 50)
    few[1:4] = a[1:4]
    @test sos_distances([few, a]; minoverlap=5)[1, 2] == 1
    @test sos_distances([few, a]; minoverlap=3)[1, 2] ≈ 0 atol = 1e-7

    # a constant pattern has no correlation -> distance 1
    @test sos_distances([fill(2.0, 50), a])[1, 2] == 1

    # Spearman sees a monotone transform as the same pattern
    @test sos_distances([a, exp.(a)]; method=:spearman)[1, 2] ≈ 0 atol = 1e-7
    @test sos_distances([a, exp.(a)])[1, 2] > 0.01
    @test_throws ArgumentError sos_distances([a, b]; method=:kendall)

    # overlap weighting by the Sorensen index of the occupied cells
    occupied = [trues(50) vcat(trues(25), falses(25))]
    weighted = sos_distances([a, b]; overlapweight=true, occupied)
    @test weighted[1, 2] ≈ sqrt(1 - (2 / 3 * cor(a, b))^2)
    @test_throws ArgumentError sos_distances([a, b]; overlapweight=true)

    # sqrt(1 - r^2) is a proper distance: the triangle inequality holds, where 1 - |r| breaks
    # it (two maps correlated at 0.7 with a third and not with each other)
    c = a + randn(rng, 50)
    trio = sos_distances([a, b, c, a + b, b - c])
    @test all(
        trio[i, j] <= trio[i, k] + trio[k, j] + 1e-12 for i in 1:5, j in 1:5, k in 1:5
    )
    @test sos_similarity(sos_distances([a, b])[1, 2]) ≈ abs(cor(a, b))
    @test sos_similarity(0) == 1 && sos_similarity(1) == 0
end

@testset "sos_distances with method = :cosine" begin
    rng = Xoshiro(4)
    a = randn(rng, 60)
    b = randn(rng, 60)
    cosine(x, y) = abs(sum(x .* y)) / sqrt(sum(abs2, x) * sum(abs2, y))

    D = sos_distances([a, b, -a]; method=:cosine)
    @test D[1, 2] ≈ sqrt(1 - cosine(a, b)^2)
    @test D[1, 3] ≈ 0 atol = 1e-7  # a mirror image is the same pattern

    # Cells without an SOS count as 0, so the comparison is over the pair's union of cells:
    # the same values on half of the range give |cos| = sqrt(1/2), not 1
    half = copy(a)
    half[31:end] .= NaN
    @test sos_similarity(sos_distances([half, a]; method=:cosine)[1, 2]) ≈
        sqrt(sum(abs2, a[1:30]) / sum(abs2, a))
    @test sos_distances([half, a])[1, 2] ≈ 0 atol = 1e-7  # a correlation would see 1

    # Disjoint ranges are unrelated; no SOS anywhere is unrelated to everything; no floor
    other = fill(NaN, 60)
    other[31:end] = b[31:end]
    @test sos_distances([half, other]; method=:cosine)[1, 2] == 1
    @test sos_distances([fill(NaN, 60), a]; method=:cosine)[1, 2] == 1
    few = fill(NaN, 60)
    few[1:2] = a[1:2]
    @test sos_distances([few, a]; method=:cosine, minoverlap=50)[1, 2] < 1

    # Every map is a vector over the same cells, so the distances are Euclidean: classical
    # MDS has no negative eigenvalues, even with gappy maps
    maps = [randn(rng, 60) for _ in 1:12]
    for (k, m) in enumerate(maps)
        m[rand(rng, 1:60, 5k)] .= NaN
    end
    Dc = sos_distances(maps; method=:cosine)
    n = size(Dc, 1)
    J = I - fill(1 / n, n, n)
    @test minimum(eigvals(Symmetric(-J * (Dc .^ 2) * J / 2))) > -1e-10
end

@testset "sos_distances on an analysis result" begin
    assemblage, tree = toy_data()
    res = node_metrics(assemblage, tree; nsims=99)
    nodes = TOY_ANALYSABLE
    @test sos_distances(res, nodes) == sos_distances([res.sos[n] for n in nodes])
    @test sos_distances(res, nodes; minoverlap=100) == [0 1 1; 1 0 1; 1 1 0]
    weighted = sos_distances(res, nodes; overlapweight=true, assemblage, tree)
    @test size(weighted) == (3, 3)
    @test all(weighted .>= sos_distances(res, nodes) .- 1e-12)
    @test_throws ArgumentError sos_distances(res, nodes; overlapweight=true)
end
