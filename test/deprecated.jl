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
