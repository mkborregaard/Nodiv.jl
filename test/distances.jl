@testset "sos_distances" begin
    rng = Xoshiro(3)
    a = randn(rng, 50)
    b = randn(rng, 50)

    D = sos_distances([a, b, -a])
    @test D == D'
    @test all(iszero, D[i, i] for i in 1:3)
    @test all(0 .<= D .<= 1)
    @test D[1, 2] ≈ 1 - abs(cor(a, b))
    @test D[1, 3] ≈ 0 atol = 1e-12  # a mirror image is the same pattern

    # NaN cells are left out pair by pair
    gappy = copy(a)
    gappy[1:10] .= NaN
    @test sos_distances([gappy, b])[1, 2] ≈ 1 - abs(cor(a[11:end], b[11:end]))

    # fewer than `minoverlap` shared cells -> distance 1
    few = fill(NaN, 50)
    few[1:4] = a[1:4]
    @test sos_distances([few, a]; minoverlap=5)[1, 2] == 1
    @test sos_distances([few, a]; minoverlap=3)[1, 2] ≈ 0 atol = 1e-12

    # a constant pattern has no correlation -> distance 1
    @test sos_distances([fill(2.0, 50), a])[1, 2] == 1

    # Spearman sees a monotone transform as the same pattern
    @test sos_distances([a, exp.(a)]; method=:spearman)[1, 2] ≈ 0 atol = 1e-12
    @test sos_distances([a, exp.(a)])[1, 2] > 0.01
    @test_throws ArgumentError sos_distances([a, b]; method=:kendall)

    # overlap weighting by the Sorensen index of the occupied cells
    occupied = [trues(50) vcat(trues(25), falses(25))]
    weighted = sos_distances([a, b]; overlapweight=true, occupied)
    @test weighted[1, 2] ≈ 1 - 2 / 3 * abs(cor(a, b))
    @test_throws ArgumentError sos_distances([a, b]; overlapweight=true)
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
