# The Plots recipes, applied with RecipesBase alone (no plotting backend needed)
@testset "plot recipes" begin
    assemblage, tree = toy_data()
    res = node_metrics(assemblage, tree; nsims=99)

    series = RecipesBase.apply_recipe(Dict{Symbol,Any}(), Nodiv.Plot_Gnd((tree, res)))
    attributes = only(series).plotattributes
    @test only(only(series).args) === tree
    @test attributes[:treetype] == :fan
    @test attributes[:marker_z] == Dict(n => g for (n, g) in res.gnd if !isnan(g))
    # markers only at the analysed nodes
    @test count(>(0), attributes[:markersize]) == length(TOY_ANALYSABLE)

    series = RecipesBase.apply_recipe(
        Dict{Symbol,Any}(), Nodiv.Plot_Node((assemblage, tree, "root", res))
    )
    @test length(series) == 4
    @test series[1].plotattributes[:title] == "root"
    @test series[2].plotattributes[:title] == "SOS"
    @test isequal(series[2].args[1], res.sos["root"])
    @test Set(s.plotattributes[:title] for s in series[3:4]) == Set(["X", "Y"])
    # a precomputed SOS vector works in place of the result
    series = RecipesBase.apply_recipe(
        Dict{Symbol,Any}(), Nodiv.Plot_Node((assemblage, tree, "X", res.sos["X"]))
    )
    @test isequal(series[2].args[1], res.sos["X"])
    @test_throws ArgumentError RecipesBase.apply_recipe(
        Dict{Symbol,Any}(), Nodiv.Plot_Node((assemblage, tree, "root"))
    )
end
