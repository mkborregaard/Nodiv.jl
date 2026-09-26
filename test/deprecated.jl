@testset "deprecations" begin
    gnd = Dict("b" => 0.9, "a" => NaN)
    sos = Dict("b" => [1.0, -1.0])
    ana = @test_deprecated NodeAnalysis(gnd, sos)
    @test ana.nodes == ["a", "b"]
    @test divergent_nodes(ana) == ["b"]

    one(x) = Dict("b" => x)
    met = @test_deprecated NodeMetrics(gnd, one(2.0), one(1.0), one(3.0), one(0.01), sos)
    @test met.nodes == ["a", "b"]
    @test isnan(met.varying["b"])
    @test divergent_nodes(met) == ["b"]
end

@testset "renamed metrics" begin
    Random.seed!(3)
    sims = Float64.(rand(0:4, 41, 10))
    sims[1, 1:4] .= 9.0
    occupied = trues(10)
    @test isequal(@test_deprecated(calculate_SOS(sims)), sos(sims))
    @test @test_deprecated(calculate_GND(sims, occupied)) == gnd(sims, occupied)
    @test @test_deprecated(gnd_rms(sos(sims))) == sos_rms(sims)
    @test @test_deprecated(gnd_spatial(sos(sims))) == sos_sd(sims)
    @test @test_deprecated(calculate_GND_ses(sims)) ≈ divergence_ses(sims)
    @test @test_deprecated(calculate_GND_pval(sims)) ≈ divergence_pval(sims)
    # the old RMS over all occupied cells, a fixed cell counted as SOS = 0
    fixed = hcat(sims, fill(2.0, 41))
    @test @test_deprecated(calculate_GND_rms(fixed, trues(11))) ≈
        sqrt(10 / 11) * sos_rms(fixed)
    @test isfinite(@test_deprecated(calculate_GND_spatial(fixed, trues(11))))
end

@testset "deprecated analysis functions" begin
    assemblage, tree = toy_data()
    SOSs, GNDs = @test_deprecated node_based_analysis(assemblage, tree; nsims=20)
    @test size(SOSs) == (nsites(assemblage), 11)
    @test count(isfinite, GNDs) == length(TOY_ANALYSABLE)
    D = @test_deprecated sos_distances(assemblage, tree, TOY_ANALYSABLE; nsims=20)
    @test size(D) == (3, 3)
end
